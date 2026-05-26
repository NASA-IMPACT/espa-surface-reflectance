//! Lookup table interpolation routines ported from `lut_subr.c`.

use crate::constants::{
    AOT550NM, LOG_AOT550NM, NAOT_VALS, NAOT_X_NSOLAR, NAOT_X_NSUNANGLE, NPRES_VALS,
    NSOLAR_ZEN_VALS, NSUNANGLE_VALS, NVIEW_ZEN_VALS, TPRES, XTS_MIN, XTS_STEP, XTV_MIN,
    XTV_STEP,
};
use crate::geometry::scattering_angle;

/// Pre-computed atmospheric lookup tables.
pub struct LookupTables {
    /// Intrinsic atmospheric reflectance: [NSR_BANDS][NPRES_VALS][NAOT_VALS][NSOLAR_VALS]
    pub rolutt: Vec<f64>,
    /// Atmospheric transmission: [NSR_BANDS][NPRES_VALS][NAOT_VALS][NSUNANGLE_VALS]
    pub transt: Vec<f64>,
    /// Spherical albedo: [NSR_BANDS][NPRES_VALS][NAOT_VALS]
    pub sphalbt: Vec<f64>,
    /// Aerosol extinction normalization: [NSR_BANDS][NPRES_VALS][NAOT_VALS]
    pub normext: Vec<f64>,
    // Angle tables: [NVIEW_ZEN_VALS][NSOLAR_ZEN_VALS]
    pub tsmax: Vec<f64>,
    pub tsmin: Vec<f64>,
    pub nbfic: Vec<f64>,
    pub nbfi: Vec<i32>,
    pub ttv: Vec<f64>,
    pub tts: [f64; NSOLAR_ZEN_VALS],
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
        // Find pressure index: ip1 is first index where pressure < TPRES[ip1]
        let mut ip1 = 0usize;
        for i in 0..NPRES_VALS {
            if pressure < TPRES[i] {
                ip1 = i;
                break;
            }
            ip1 = i;
        }
        // Clamp ip1 so ip2 = ip1+1 is valid
        if ip1 >= NPRES_VALS - 1 {
            ip1 = NPRES_VALS - 2;
        }
        let ip2 = ip1 + 1;

        // Find AOT index: iaot1 is last index where raot550nm > AOT550NM[iaot1]
        let mut iaot1 = 0usize;
        for i in 0..NAOT_VALS {
            if raot550nm > AOT550NM[i] {
                iaot1 = i;
            }
        }
        // Clamp iaot1 so iaot2 = iaot1+1 is valid
        if iaot1 >= NAOT_VALS - 1 {
            iaot1 = NAOT_VALS - 2;
        }
        let iaot2 = iaot1 + 1;

        LutIndices { ip1, ip2, iaot1, iaot2 }
    }

    /// Compute spherical albedo via bilinear interpolation across AOT and pressure.
    ///
    /// Ports `compsalb` from `lut_subr.c`.
    pub fn interp_spherical_albedo(
        &self,
        indices: &LutIndices,
        iband: usize,
        pressure: f64,
        raot550nm: f64,
    ) -> f64 {
        let LutIndices { ip1, ip2, iaot1, iaot2 } = *indices;

        let deltaaot =
            (raot550nm - AOT550NM[iaot1]) / (AOT550NM[iaot2] - AOT550NM[iaot1]);

        let iband_base = iband * NPRES_VALS * NAOT_VALS;
        let ip1_base = ip1 * NAOT_VALS;
        let ip2_base = ip2 * NAOT_VALS;

        // Interpolate along AOT at pressure ip1
        let v_ip1_iaot1 = self.sphalbt[iband_base + ip1_base + iaot1];
        let v_ip1_iaot2 = self.sphalbt[iband_base + ip1_base + iaot2];
        let satm1 = v_ip1_iaot1 + (v_ip1_iaot2 - v_ip1_iaot1) * deltaaot;

        // Interpolate along AOT at pressure ip2
        let v_ip2_iaot1 = self.sphalbt[iband_base + ip2_base + iaot1];
        let v_ip2_iaot2 = self.sphalbt[iband_base + ip2_base + iaot2];
        let satm2 = v_ip2_iaot1 + (v_ip2_iaot2 - v_ip2_iaot1) * deltaaot;

        // Interpolate along pressure
        let dpres = (pressure - TPRES[ip1]) / (TPRES[ip2] - TPRES[ip1]);
        satm1 + (satm2 - satm1) * dpres
    }

    /// Compute atmospheric transmission via 3D interpolation (zenith, AOT, pressure).
    ///
    /// Ports `comptrans` from `lut_subr.c`.
    pub fn interp_transmission(
        &self,
        indices: &LutIndices,
        iband: usize,
        pressure: f64,
        raot550nm: f64,
        xts: f64,
    ) -> f64 {
        let LutIndices { ip1, ip2, iaot1, iaot2 } = *indices;

        // Find sun angle index
        let its = if xts <= XTS_MIN {
            0
        } else {
            ((xts - XTS_MIN) / XTS_STEP) as usize
        };
        let its = its.min(NSUNANGLE_VALS - 2);

        // Fractional angle
        let xmts = (xts - self.tts[its]) * 0.25;

        let deltaaot =
            (raot550nm - AOT550NM[iaot1]) / (AOT550NM[iaot2] - AOT550NM[iaot1]);

        let iband_base = iband * NPRES_VALS * NAOT_X_NSUNANGLE;
        let ip1_base = ip1 * NAOT_X_NSUNANGLE;
        let ip2_base = ip2 * NAOT_X_NSUNANGLE;
        let iaot1_base = iaot1 * NSUNANGLE_VALS;
        let iaot2_base = iaot2 * NSUNANGLE_VALS;

        // ip1, iaot1: interpolate along angle
        let base = iband_base + ip1_base + iaot1_base;
        let xtranst = self.transt[base + its];
        let xtiaot1_ip1 = xtranst + (self.transt[base + its + 1] - xtranst) * xmts;

        // ip1, iaot2: interpolate along angle
        let base = iband_base + ip1_base + iaot2_base;
        let xtranst = self.transt[base + its];
        let xtiaot2_ip1 = xtranst + (self.transt[base + its + 1] - xtranst) * xmts;

        // Interpolate along AOT for ip1
        let xtts1 = xtiaot1_ip1 + (xtiaot2_ip1 - xtiaot1_ip1) * deltaaot;

        // ip2, iaot1: interpolate along angle
        let base = iband_base + ip2_base + iaot1_base;
        let xtranst = self.transt[base + its];
        let xtiaot1_ip2 = xtranst + (self.transt[base + its + 1] - xtranst) * xmts;

        // ip2, iaot2: interpolate along angle
        let base = iband_base + ip2_base + iaot2_base;
        let xtranst = self.transt[base + its];
        let xtiaot2_ip2 = xtranst + (self.transt[base + its + 1] - xtranst) * xmts;

        // Interpolate along AOT for ip2
        let xtts2 = xtiaot1_ip2 + (xtiaot2_ip2 - xtiaot1_ip2) * deltaaot;

        // Interpolate along pressure
        let dpres = (pressure - TPRES[ip1]) / (TPRES[ip2] - TPRES[ip1]);
        xtts1 + (xtts2 - xtts1) * dpres
    }

    /// Private helper: interpolate intrinsic reflectance at scattering angle.
    ///
    /// Ports `interp_refl_using_scat_angle` from `lut_subr.c`.
    ///
    /// # Arguments
    /// * `iband` - Band index
    /// * `ip` - Pressure index
    /// * `iaot` - AOT index
    /// * `scaa` - Scattering angle (degrees)
    /// * `its` - Solar zenith table index
    /// * `itv` - View zenith table index
    /// * `t` - Solar angle interpolation parameter
    /// * `u` - View angle interpolation parameter
    fn interp_refl_at_scat_angle(
        &self,
        iband: usize,
        ip: usize,
        iaot: usize,
        scaa: f64,
        its: usize,
        itv: usize,
        t: f64,
        u: f64,
    ) -> f64 {
        let rolutt_base = iband * NPRES_VALS * NAOT_X_NSOLAR
            + ip * NAOT_X_NSOLAR
            + iaot * crate::constants::NSOLAR_VALS;

        let mut ro = [0.0f64; 4];

        // 4 corners: (itv, its), (itv, its+1), (itv+1, its), (itv+1, its+1)
        // i=0: (itv,   its  )  is = its,   iv = itv
        // i=1: (itv,   its+1)  is = its+1, iv = itv
        // i=2: (itv+1, its  )  is = its,   iv = itv+1
        // i=3: (itv+1, its+1)  is = its+1, iv = itv+1
        for i in 0..4usize {
            let is = its + (i % 2); // its or its+1
            let iv = if i < 2 { itv } else { itv + 1 }; // itv or itv+1

            let angle_idx = iv * NSOLAR_ZEN_VALS + is;
            let xtsmax_i = self.tsmax[angle_idx];
            let xtsmin_i = self.tsmin[angle_idx];
            let nbfi_i = self.nbfi[angle_idx];
            let nbfic_i = self.nbfic[angle_idx];

            // offset within the NSOLAR_VALS block for this (iv, is) combination
            // j corresponds to the cumulative azimuth offset
            let j = (nbfic_i - nbfi_i as f64) as usize;

            if is != 0 && iv != 0 {
                let mut isca = ((xtsmax_i - scaa) * 0.25 + 1.0) as usize;
                if isca == 0 {
                    isca = 1;
                }

                let (sca1, sca2, isca_used) = if isca + 1 < nbfi_i as usize {
                    let s1 = xtsmax_i - (isca as f64 - 1.0) * 4.0;
                    let s2 = s1 - 4.0;
                    (s1, s2, isca)
                } else {
                    let isca_c = (nbfi_i as usize).saturating_sub(1);
                    let s1 = xtsmax_i - (isca_c as f64 - 1.0) * 4.0;
                    let s2 = xtsmin_i;
                    (s1, s2, isca_c)
                };

                let roinf_idx = rolutt_base + j + isca_used.saturating_sub(1);
                let rosup_idx = rolutt_base + j + isca_used;
                let roinf = if roinf_idx < self.rolutt.len() {
                    self.rolutt[roinf_idx]
                } else {
                    0.0
                };
                let rosup = if rosup_idx < self.rolutt.len() {
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

        // Bilinear interpolation of the 4 corner values
        ro[3] + u * (ro[1] - ro[3]) + t * (ro[2] - ro[3]) + u * t * (ro[0] - ro[1] - ro[2] + ro[3])
    }

    /// Compute intrinsic atmospheric reflectance via complex 5D interpolation.
    ///
    /// Ports `comproatm` from `lut_subr.c`.
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

        // Compute scattering angle
        let scaa = scattering_angle(xmus, xmuv, cosxfi);

        // View zenith index
        let itv_f = (xtv - XTV_MIN) / XTV_STEP + 1.0;
        let itv = (itv_f as usize).min(NVIEW_ZEN_VALS - 2);

        // Solar zenith index
        let its_f = (xts - XTS_MIN) / XTS_STEP;
        let its = if xts <= XTS_MIN {
            0
        } else {
            (its_f as usize).min(NSOLAR_ZEN_VALS - 2)
        };

        // Interpolation parameters t (solar) and u (view)
        let itv_its_indx = itv * NSOLAR_ZEN_VALS + its;
        let itv1_its_indx = itv_its_indx + NSOLAR_ZEN_VALS;

        let tts_denom = self.tts[its + 1] - self.tts[its];
        let t = if tts_denom.abs() > f64::EPSILON {
            (self.tts[its + 1] - xts) / tts_denom
        } else {
            0.0
        };

        let ttv_denom = self.ttv[itv1_its_indx] - self.ttv[itv_its_indx];
        let u = if ttv_denom.abs() > f64::EPSILON {
            (self.ttv[itv1_its_indx] - xtv) / ttv_denom
        } else {
            0.0
        };

        // AOT interpolation in LOG space
        let deltaaot = (raot550nm.ln() - LOG_AOT550NM[iaot1])
            / (LOG_AOT550NM[iaot2] - LOG_AOT550NM[iaot1]);

        // ip1: interpolate at iaot1 and iaot2, then in log-AOT space
        let roiaot1 = self.interp_refl_at_scat_angle(iband, ip1, iaot1, scaa, its, itv, t, u);
        let roiaot2 = self.interp_refl_at_scat_angle(iband, ip1, iaot2, scaa, its, itv, t, u);
        let rop1 = roiaot1 + (roiaot2 - roiaot1) * deltaaot;

        // ip2: interpolate at iaot1 and iaot2, then in log-AOT space
        let roiaot1 = self.interp_refl_at_scat_angle(iband, ip2, iaot1, scaa, its, itv, t, u);
        let roiaot2 = self.interp_refl_at_scat_angle(iband, ip2, iaot2, scaa, its, itv, t, u);
        let rop2 = roiaot1 + (roiaot2 - roiaot1) * deltaaot;

        // Interpolate along pressure
        let dpres = (pressure - TPRES[ip1]) / (TPRES[ip2] - TPRES[ip1]);
        rop1 + (rop2 - rop1) * dpres
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{NAOT_X_NSOLAR, NAOT_VALS, NPRES_VALS, NSOLAR_ZEN_VALS, NVIEW_ZEN_VALS};

    fn make_test_lut() -> LookupTables {
        let nsr = 1;
        let mut sphalbt = vec![0.0; nsr * NPRES_VALS * NAOT_VALS];
        for ip in 0..NPRES_VALS {
            for ia in 0..NAOT_VALS {
                sphalbt[ip * NAOT_VALS + ia] = AOT550NM[ia] * 0.1;
            }
        }
        let mut transt = vec![0.0; nsr * NPRES_VALS * NAOT_X_NSUNANGLE];
        for ip in 0..NPRES_VALS {
            for ia in 0..NAOT_VALS {
                for is_ in 0..NSUNANGLE_VALS {
                    let idx = ip * NAOT_X_NSUNANGLE + ia * NSUNANGLE_VALS + is_;
                    transt[idx] = 1.0 - AOT550NM[ia] * 0.1;
                }
            }
        }
        let normext = vec![1.0; nsr * NPRES_VALS * NAOT_VALS];
        let rolutt = vec![0.01; nsr * NPRES_VALS * NAOT_X_NSOLAR];
        let tsmax = vec![180.0; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let tsmin = vec![0.0; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let nbfic = vec![45.0; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let nbfi = vec![45; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let ttv = vec![3.0; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let mut tts = [0.0f64; NSOLAR_ZEN_VALS];
        for i in 0..NSOLAR_ZEN_VALS {
            tts[i] = i as f64 * 4.0;
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
