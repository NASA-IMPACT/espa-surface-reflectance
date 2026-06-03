//! Lookup table interpolation routines ported from `lut_subr.c`.

use crate::constants::{
    AOT550NM, LOG_AOT550NM, NAOT_VALS, NAOT_X_NSOLAR, NAOT_X_NSUNANGLE, NPRES_VALS,
    NSOLAR_ZEN_VALS, NSUNANGLE_VALS, NVIEW_ZEN_VALS, TPRES, XTS_MIN, XTS_STEP, XTV_MIN,
    XTV_STEP,
};
use crate::geometry::scattering_angle;

/// Pre-computed atmospheric lookup tables.
/// All float data is stored as f32 to match C's `float` precision.
pub struct LookupTables {
    /// Intrinsic atmospheric reflectance: [NSR_BANDS][NPRES_VALS][NAOT_VALS][NSOLAR_VALS]
    pub rolutt: Vec<f32>,
    /// Atmospheric transmission: [NSR_BANDS][NPRES_VALS][NAOT_VALS][NSUNANGLE_VALS]
    pub transt: Vec<f32>,
    /// Spherical albedo: [NSR_BANDS][NPRES_VALS][NAOT_VALS]
    pub sphalbt: Vec<f32>,
    /// Aerosol extinction normalization: [NSR_BANDS][NPRES_VALS][NAOT_VALS]
    pub normext: Vec<f32>,
    // Angle tables: [NVIEW_ZEN_VALS][NSOLAR_ZEN_VALS]
    pub tsmax: Vec<f32>,
    pub tsmin: Vec<f32>,
    pub nbfic: Vec<f32>,
    pub nbfi: Vec<i32>,
    pub ttv: Vec<f32>,
    pub tts: [f32; NSOLAR_ZEN_VALS],
    /// Cumulative azimuth angle offset for each solar zenith index.
    /// Used to index into the rolutt NSOLAR_VALS block.
    /// Loaded from the INDTS SDS in the angle HDF file.
    pub indts: Vec<i32>,
    /// Number of SR bands (8 for Landsat, 11 for Sentinel)
    pub nsr_bands: usize,
}

/// Indices into the pressure and AOT lookup tables.
#[derive(Debug, Clone, Copy)]
pub struct LutIndices {
    pub ip1: usize,
    pub ip2: usize,
    pub iaot1: usize,
    pub iaot2: usize,
}

impl LookupTables {
    /// Find pressure and AOT table indices for the given values.
    ///
    /// Searches TPRES[7] to find ip1 where pressure < TPRES[ip1], ip2 = ip1+1.
    /// Searches AOT550NM[22] to find iaot1 where raot550nm > AOT550NM[iaot1], iaot2 = iaot1+1.
    /// Both are clamped to valid ranges.
    pub fn find_indices(&self, pressure: f64, raot550nm: f64) -> LutIndices {
        let pres = pressure as f32;
        let raot = raot550nm as f32;
        let mut ip1 = 0usize;
        for ip in 0..NPRES_VALS - 1 {
            if pres < TPRES[ip] as f32 {
                ip1 = ip;
            }
        }
        let ip2 = ip1 + 1;

        let mut iaot1 = 0usize;
        for i in 0..NAOT_VALS {
            if raot > AOT550NM[i] as f32 {
                iaot1 = i;
            }
        }
        if iaot1 >= NAOT_VALS - 1 {
            iaot1 = NAOT_VALS - 2;
        }
        let iaot2 = iaot1 + 1;

        LutIndices { ip1, ip2, iaot1, iaot2 }
    }

    /// Compute spherical albedo via bilinear interpolation across AOT and pressure.
    /// All arithmetic in f32 matching C's float.
    pub fn interp_spherical_albedo(
        &self,
        indices: &LutIndices,
        iband: usize,
        pressure: f64,
        raot550nm: f64,
    ) -> f64 {
        let LutIndices { ip1, ip2, iaot1, iaot2 } = *indices;
        let raot = raot550nm as f32;
        let pres = pressure as f32;

        let aot1 = AOT550NM[iaot1] as f32;
        let aot2 = AOT550NM[iaot2] as f32;
        let deltaaot: f32 = (raot - aot1) / (aot2 - aot1);

        let iband_base = iband * NPRES_VALS * NAOT_VALS;
        let ip1_base = ip1 * NAOT_VALS;
        let ip2_base = ip2 * NAOT_VALS;

        let v_ip1_iaot1 = self.sphalbt[iband_base + ip1_base + iaot1];
        let v_ip1_iaot2 = self.sphalbt[iband_base + ip1_base + iaot2];
        let satm1: f32 = v_ip1_iaot1 + (v_ip1_iaot2 - v_ip1_iaot1) * deltaaot;

        let v_ip2_iaot1 = self.sphalbt[iband_base + ip2_base + iaot1];
        let v_ip2_iaot2 = self.sphalbt[iband_base + ip2_base + iaot2];
        let satm2: f32 = v_ip2_iaot1 + (v_ip2_iaot2 - v_ip2_iaot1) * deltaaot;

        let tpres1 = TPRES[ip1] as f32;
        let tpres2 = TPRES[ip2] as f32;
        let dpres: f32 = (pres - tpres1) / (tpres2 - tpres1);
        (satm1 + (satm2 - satm1) * dpres) as f64
    }

    /// Compute atmospheric transmission via 3D interpolation (zenith, AOT, pressure).
    /// All arithmetic in f32 matching C's float.
    pub fn interp_transmission(
        &self,
        indices: &LutIndices,
        iband: usize,
        pressure: f64,
        raot550nm: f64,
        xts: f64,
    ) -> f64 {
        let LutIndices { ip1, ip2, iaot1, iaot2 } = *indices;
        let xts_f = xts as f32;
        let raot = raot550nm as f32;
        let pres = pressure as f32;

        let its = if xts_f <= XTS_MIN as f32 {
            0
        } else {
            ((xts_f - XTS_MIN as f32) / XTS_STEP as f32) as usize
        };
        let its = its.min(NSUNANGLE_VALS - 2);

        let xmts: f32 = (xts_f - self.tts[its]) * 0.25f32;

        let aot1 = AOT550NM[iaot1] as f32;
        let aot2 = AOT550NM[iaot2] as f32;
        let deltaaot: f32 = (raot - aot1) / (aot2 - aot1);

        let iband_base = iband * NPRES_VALS * NAOT_X_NSUNANGLE;
        let ip1_base = ip1 * NAOT_X_NSUNANGLE;
        let ip2_base = ip2 * NAOT_X_NSUNANGLE;
        let iaot1_base = iaot1 * NSUNANGLE_VALS;
        let iaot2_base = iaot2 * NSUNANGLE_VALS;

        let base = iband_base + ip1_base + iaot1_base;
        let xtranst = self.transt[base + its];
        let xtiaot1_ip1: f32 = xtranst + (self.transt[base + its + 1] - xtranst) * xmts;

        let base = iband_base + ip1_base + iaot2_base;
        let xtranst = self.transt[base + its];
        let xtiaot2_ip1: f32 = xtranst + (self.transt[base + its + 1] - xtranst) * xmts;

        let xtts1: f32 = xtiaot1_ip1 + (xtiaot2_ip1 - xtiaot1_ip1) * deltaaot;

        let base = iband_base + ip2_base + iaot1_base;
        let xtranst = self.transt[base + its];
        let xtiaot1_ip2: f32 = xtranst + (self.transt[base + its + 1] - xtranst) * xmts;

        let base = iband_base + ip2_base + iaot2_base;
        let xtranst = self.transt[base + its];
        let xtiaot2_ip2: f32 = xtranst + (self.transt[base + its + 1] - xtranst) * xmts;

        let xtts2: f32 = xtiaot1_ip2 + (xtiaot2_ip2 - xtiaot1_ip2) * deltaaot;

        let tpres1 = TPRES[ip1] as f32;
        let tpres2 = TPRES[ip2] as f32;
        let dpres: f32 = (pres - tpres1) / (tpres2 - tpres1);
        (xtts1 + (xtts2 - xtts1) * dpres) as f64
    }

    /// Private helper: interpolate intrinsic reflectance at scattering angle.
    /// All arithmetic in f32 matching C's float.
    fn interp_refl_at_scat_angle(
        &self,
        iband: usize,
        ip: usize,
        iaot: usize,
        scaa: f32,
        its: usize,
        itv: usize,
        t: f32,
        u: f32,
    ) -> f32 {
        let rolutt_base = iband * NPRES_VALS * NAOT_X_NSOLAR
            + ip * NAOT_X_NSOLAR
            + iaot * crate::constants::NSOLAR_VALS;

        let mut ro = [0.0f32; 4];

        for i in 0..4usize {
            let is = its + (i % 2);
            let iv = if i < 2 { itv } else { itv + 1 };

            let angle_idx = iv * NSOLAR_ZEN_VALS + is;
            let xtsmax_i = self.tsmax[angle_idx];
            let xtsmin_i = self.tsmin[angle_idx];
            let nbfi_i = self.nbfi[angle_idx];
            let nbfic_i = self.nbfic[angle_idx];

            let j = self.indts[is] as usize + (nbfic_i - nbfi_i as f32) as usize;

            if is != 0 && iv != 0 {
                let mut isca = ((xtsmax_i - scaa) * 0.25f32 + 1.0f32) as usize;
                if isca == 0 {
                    isca = 1;
                }

                let (sca1, sca2, isca_used) = if isca + 1 < nbfi_i as usize {
                    let s1 = xtsmax_i - (isca as f32 - 1.0f32) * 4.0f32;
                    let s2 = s1 - 4.0f32;
                    (s1, s2, isca)
                } else {
                    let isca_c = (nbfi_i as usize).saturating_sub(1);
                    let s1 = xtsmax_i - (isca_c as f32 - 1.0f32) * 4.0f32;
                    let s2 = xtsmin_i;
                    (s1, s2, isca_c)
                };

                let roinf_idx = rolutt_base + j + isca_used.saturating_sub(1);
                let rosup_idx = rolutt_base + j + isca_used;
                let roinf: f32 = if roinf_idx < self.rolutt.len() {
                    self.rolutt[roinf_idx]
                } else {
                    0.0
                };
                let rosup: f32 = if rosup_idx < self.rolutt.len() {
                    self.rolutt[rosup_idx]
                } else {
                    0.0
                };

                ro[i] = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
            } else {
                let idx = rolutt_base + j;
                ro[i] = if idx < self.rolutt.len() { self.rolutt[idx] } else { 0.0 };
            }
        }

        ro[3] + u * (ro[1] - ro[3]) + t * (ro[2] - ro[3]) + u * t * (ro[0] - ro[1] - ro[2] + ro[3])
    }

    /// Compute intrinsic atmospheric reflectance via complex 5D interpolation.
    /// All arithmetic in f32 matching C's float.
    #[allow(clippy::too_many_arguments)]
    pub fn interp_atmospheric_reflectance(
        &self,
        indices: &LutIndices,
        iband: usize,
        pressure: f64,
        raot550nm: f64,
        xts: f64,
        xtv: f64,
        xmus: f64,
        xmuv: f64,
        cosxfi: f64,
    ) -> f64 {
        let LutIndices { ip1, ip2, iaot1, iaot2 } = *indices;
        let pres = pressure as f32;

        // Compute scattering angle (in f64, truncated to f32 like C)
        let scaa = scattering_angle(xmus, xmuv, cosxfi) as f32;

        let xtv_f = xtv as f32;
        let xts_f = xts as f32;

        // View zenith index
        let itv_f = (xtv_f - XTV_MIN as f32) / XTV_STEP as f32 + 1.0f32;
        let itv = (itv_f as usize).min(NVIEW_ZEN_VALS - 2);

        // Solar zenith index
        let its = if xts_f <= XTS_MIN as f32 {
            0
        } else {
            let its_f = (xts_f - XTS_MIN as f32) / XTS_STEP as f32;
            (its_f as usize).min(NSOLAR_ZEN_VALS - 2)
        };

        // Interpolation parameters t (solar) and u (view)
        let itv_its_indx = itv * NSOLAR_ZEN_VALS + its;
        let itv1_its_indx = itv_its_indx + NSOLAR_ZEN_VALS;

        let tts_denom = self.tts[its + 1] - self.tts[its];
        let t: f32 = if tts_denom.abs() > f32::EPSILON {
            (self.tts[its + 1] - xts_f) / tts_denom
        } else {
            0.0
        };

        let ttv_denom = self.ttv[itv1_its_indx] - self.ttv[itv_its_indx];
        let u: f32 = if ttv_denom.abs() > f32::EPSILON {
            (self.ttv[itv1_its_indx] - xtv_f) / ttv_denom
        } else {
            0.0
        };

        // AOT interpolation in LOG space (C uses float for log)
        let raot_f = raot550nm as f32;
        let log_raot = raot_f.ln();
        let log_aot1 = LOG_AOT550NM[iaot1] as f32;
        let log_aot2 = LOG_AOT550NM[iaot2] as f32;
        let deltaaot: f32 = (log_raot - log_aot1) / (log_aot2 - log_aot1);

        let roiaot1 = self.interp_refl_at_scat_angle(iband, ip1, iaot1, scaa, its, itv, t, u);
        let roiaot2 = self.interp_refl_at_scat_angle(iband, ip1, iaot2, scaa, its, itv, t, u);
        let rop1: f32 = roiaot1 + (roiaot2 - roiaot1) * deltaaot;

        let roiaot1 = self.interp_refl_at_scat_angle(iband, ip2, iaot1, scaa, its, itv, t, u);
        let roiaot2 = self.interp_refl_at_scat_angle(iband, ip2, iaot2, scaa, its, itv, t, u);
        let rop2: f32 = roiaot1 + (roiaot2 - roiaot1) * deltaaot;

        let tpres1 = TPRES[ip1] as f32;
        let tpres2 = TPRES[ip2] as f32;
        let dpres: f32 = (pres - tpres1) / (tpres2 - tpres1);
        (rop1 + (rop2 - rop1) * dpres) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{NAOT_X_NSOLAR, NAOT_VALS, NPRES_VALS, NSOLAR_ZEN_VALS, NVIEW_ZEN_VALS};

    fn make_test_lut() -> LookupTables {
        let nsr = 1;
        let mut sphalbt = vec![0.0f32; nsr * NPRES_VALS * NAOT_VALS];
        for ip in 0..NPRES_VALS {
            for ia in 0..NAOT_VALS {
                sphalbt[ip * NAOT_VALS + ia] = AOT550NM[ia] as f32 * 0.1;
            }
        }
        let mut transt = vec![0.0f32; nsr * NPRES_VALS * NAOT_X_NSUNANGLE];
        for ip in 0..NPRES_VALS {
            for ia in 0..NAOT_VALS {
                for is_ in 0..NSUNANGLE_VALS {
                    let idx = ip * NAOT_X_NSUNANGLE + ia * NSUNANGLE_VALS + is_;
                    transt[idx] = 1.0 - AOT550NM[ia] as f32 * 0.1;
                }
            }
        }
        let normext = vec![1.0f32; nsr * NPRES_VALS * NAOT_VALS];
        let rolutt = vec![0.01f32; nsr * NPRES_VALS * NAOT_X_NSOLAR];
        let tsmax = vec![180.0f32; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let tsmin = vec![0.0f32; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let nbfic = vec![45.0f32; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let nbfi = vec![45; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let ttv = vec![3.0f32; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let mut tts = [0.0f32; NSOLAR_ZEN_VALS];
        for i in 0..NSOLAR_ZEN_VALS {
            tts[i] = i as f32 * 4.0;
        }
        LookupTables {
            rolutt,
            transt,
            sphalbt,
            normext,
            tsmax,
            tsmin,
            nbfic,
            nbfi,
            ttv,
            tts,
            indts: vec![0; NSUNANGLE_VALS],
            nsr_bands: nsr,
        }
    }

    #[test]
    fn test_find_indices_sea_level() {
        let lut = make_test_lut();
        let idx = lut.find_indices(1013.0, 0.15);
        assert_eq!(idx.ip1, 0);
        assert_eq!(idx.ip2, 1);
        assert_eq!(idx.iaot1, 2);
        assert_eq!(idx.iaot2, 3);
    }

    #[test]
    fn test_interp_spherical_albedo_monotonic() {
        let lut = make_test_lut();
        let idx1 = lut.find_indices(1013.0, 0.1);
        let sa1 = lut.interp_spherical_albedo(&idx1, 0, 1013.0, 0.1);
        let idx2 = lut.find_indices(1013.0, 0.5);
        let sa2 = lut.interp_spherical_albedo(&idx2, 0, 1013.0, 0.5);
        assert!(sa2 > sa1, "Higher AOT = higher albedo: {sa1} vs {sa2}");
    }

    #[test]
    fn test_interp_transmission_range() {
        let lut = make_test_lut();
        let idx = lut.find_indices(1013.0, 0.1);
        let trans = lut.interp_transmission(&idx, 0, 1013.0, 0.1, 30.0);
        assert!(trans > 0.0 && trans <= 1.0, "trans={trans}");
    }

    #[test]
    fn test_interp_atmospheric_reflectance_positive() {
        let lut = make_test_lut();
        let idx = lut.find_indices(1013.0, 0.1);
        let roatm = lut.interp_atmospheric_reflectance(
            &idx, 0, 1013.0, 0.1, 30.0, 5.0, 0.866, 0.996, 0.5,
        );
        assert!(roatm >= 0.0, "roatm={roatm}");
    }
}
