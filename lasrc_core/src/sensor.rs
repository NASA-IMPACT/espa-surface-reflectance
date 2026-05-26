//! Sensor trait and implementations for Landsat 8/9 and Sentinel-2A/B/C.

use crate::constants::{
    LAERO_WINDOW, LHALF_AERO_WINDOW, LFIX_AERO_WINDOW, LHALF_FIX_AERO_WINDOW, LMIN_CLEAR_PIX,
    SAERO_WINDOW, SFIX_AERO_WINDOW, SHALF_FIX_AERO_WINDOW, SMIN_CLEAR_PIX,
    TAURAY_LANDSAT, LAMBDA_LANDSAT, LAMBDA_SENTINEL,
};

/// Configuration for a single spectral band.
#[derive(Debug, Clone)]
pub struct BandConfig {
    pub name: &'static str,
    pub wavelength_um: f64,
    pub native_resolution_m: f64,
}

/// Logical band index mapping for standard spectral positions.
#[derive(Debug, Clone)]
pub struct BandIndices {
    pub coastal: usize,
    pub blue: usize,
    pub green: usize,
    pub red: usize,
    pub nir: usize,
    pub swir1: usize,
    pub swir2: usize,
}

/// Trait representing a sensor's physical and algorithmic configuration.
pub trait Sensor: Send + Sync {
    fn name(&self) -> &str;
    fn reflectance_bands(&self) -> &[BandConfig];
    fn thermal_bands(&self) -> &[BandConfig];
    fn output_resolution_m(&self) -> f64;
    fn aerosol_window(&self) -> usize;
    fn half_aerosol_window(&self) -> usize;
    fn fix_aerosol_window(&self) -> usize;
    fn half_fix_aerosol_window(&self) -> usize;
    fn min_clear_pix(&self) -> usize;
    fn band_indices(&self) -> &BandIndices;
    fn tauray(&self) -> &[f64];
    fn lambda(&self) -> &[f64];
    fn num_refl_bands(&self) -> usize {
        self.reflectance_bands().len()
    }
}

// ── Landsat 8/9 constants ──────────────────────────────────────────────────

const LANDSAT_REFL_BANDS: [BandConfig; 8] = [
    BandConfig { name: "band1", wavelength_um: 0.443, native_resolution_m: 30.0 },
    BandConfig { name: "band2", wavelength_um: 0.480, native_resolution_m: 30.0 },
    BandConfig { name: "band3", wavelength_um: 0.585, native_resolution_m: 30.0 },
    BandConfig { name: "band4", wavelength_um: 0.655, native_resolution_m: 30.0 },
    BandConfig { name: "band5", wavelength_um: 0.865, native_resolution_m: 30.0 },
    BandConfig { name: "band6", wavelength_um: 1.610, native_resolution_m: 30.0 },
    BandConfig { name: "band7", wavelength_um: 2.200, native_resolution_m: 30.0 },
    BandConfig { name: "band9", wavelength_um: 1.370, native_resolution_m: 30.0 },
];

const LANDSAT_THM_BANDS: [BandConfig; 2] = [
    BandConfig { name: "band10", wavelength_um: 10.9, native_resolution_m: 100.0 },
    BandConfig { name: "band11", wavelength_um: 12.0, native_resolution_m: 100.0 },
];

const LANDSAT_BAND_INDICES: BandIndices = BandIndices {
    coastal: 0,
    blue: 1,
    green: 2,
    red: 3,
    nir: 4,
    swir1: 5,
    swir2: 6,
};

// ── Sentinel-2 constants ───────────────────────────────────────────────────

const SENTINEL_REFL_BANDS: [BandConfig; 11] = [
    BandConfig { name: "B01", wavelength_um: 0.443, native_resolution_m: 60.0 },
    BandConfig { name: "B02", wavelength_um: 0.490, native_resolution_m: 10.0 },
    BandConfig { name: "B03", wavelength_um: 0.560, native_resolution_m: 10.0 },
    BandConfig { name: "B04", wavelength_um: 0.665, native_resolution_m: 10.0 },
    BandConfig { name: "B05", wavelength_um: 0.705, native_resolution_m: 20.0 },
    BandConfig { name: "B06", wavelength_um: 0.740, native_resolution_m: 20.0 },
    BandConfig { name: "B07", wavelength_um: 0.783, native_resolution_m: 20.0 },
    BandConfig { name: "B08", wavelength_um: 0.842, native_resolution_m: 10.0 },
    BandConfig { name: "B8A", wavelength_um: 0.865, native_resolution_m: 20.0 },
    BandConfig { name: "B11", wavelength_um: 1.610, native_resolution_m: 20.0 },
    BandConfig { name: "B12", wavelength_um: 2.190, native_resolution_m: 20.0 },
];

const SENTINEL_BAND_INDICES: BandIndices = BandIndices {
    coastal: 0,
    blue: 1,
    green: 2,
    red: 3,
    nir: 8,   // B8A
    swir1: 9,
    swir2: 10,
};

// ── Concrete sensor types ──────────────────────────────────────────────────

/// Landsat 8 sensor.
pub struct Landsat8;

/// Landsat 9 sensor (same band configuration as Landsat 8).
pub struct Landsat9;

/// Sentinel-2A sensor.
pub struct Sentinel2A;

/// Sentinel-2B sensor.
pub struct Sentinel2B;

/// Sentinel-2C sensor.
pub struct Sentinel2C;

// ── Landsat 8 impl ─────────────────────────────────────────────────────────

impl Sensor for Landsat8 {
    fn name(&self) -> &str { "LANDSAT_8" }
    fn reflectance_bands(&self) -> &[BandConfig] { &LANDSAT_REFL_BANDS }
    fn thermal_bands(&self) -> &[BandConfig] { &LANDSAT_THM_BANDS }
    fn output_resolution_m(&self) -> f64 { 30.0 }
    fn aerosol_window(&self) -> usize { LAERO_WINDOW }
    fn half_aerosol_window(&self) -> usize { LHALF_AERO_WINDOW }
    fn fix_aerosol_window(&self) -> usize { LFIX_AERO_WINDOW }
    fn half_fix_aerosol_window(&self) -> usize { LHALF_FIX_AERO_WINDOW }
    fn min_clear_pix(&self) -> usize { LMIN_CLEAR_PIX }
    fn band_indices(&self) -> &BandIndices { &LANDSAT_BAND_INDICES }
    fn tauray(&self) -> &[f64] { &TAURAY_LANDSAT }
    fn lambda(&self) -> &[f64] { &LAMBDA_LANDSAT }
}

// ── Landsat 9 impl ─────────────────────────────────────────────────────────

impl Sensor for Landsat9 {
    fn name(&self) -> &str { "LANDSAT_9" }
    fn reflectance_bands(&self) -> &[BandConfig] { &LANDSAT_REFL_BANDS }
    fn thermal_bands(&self) -> &[BandConfig] { &LANDSAT_THM_BANDS }
    fn output_resolution_m(&self) -> f64 { 30.0 }
    fn aerosol_window(&self) -> usize { LAERO_WINDOW }
    fn half_aerosol_window(&self) -> usize { LHALF_AERO_WINDOW }
    fn fix_aerosol_window(&self) -> usize { LFIX_AERO_WINDOW }
    fn half_fix_aerosol_window(&self) -> usize { LHALF_FIX_AERO_WINDOW }
    fn min_clear_pix(&self) -> usize { LMIN_CLEAR_PIX }
    fn band_indices(&self) -> &BandIndices { &LANDSAT_BAND_INDICES }
    fn tauray(&self) -> &[f64] { &TAURAY_LANDSAT }
    fn lambda(&self) -> &[f64] { &LAMBDA_LANDSAT }
}

// ── Sentinel-2 shared impl helper (macro) ─────────────────────────────────

macro_rules! impl_sentinel_sensor {
    ($type:ty, $name:expr) => {
        impl Sensor for $type {
            fn name(&self) -> &str { $name }
            fn reflectance_bands(&self) -> &[BandConfig] { &SENTINEL_REFL_BANDS }
            fn thermal_bands(&self) -> &[BandConfig] { &[] }
            fn output_resolution_m(&self) -> f64 { 10.0 }
            fn aerosol_window(&self) -> usize { SAERO_WINDOW }
            fn half_aerosol_window(&self) -> usize { SAERO_WINDOW / 2 }
            fn fix_aerosol_window(&self) -> usize { SFIX_AERO_WINDOW }
            fn half_fix_aerosol_window(&self) -> usize { SHALF_FIX_AERO_WINDOW }
            fn min_clear_pix(&self) -> usize { SMIN_CLEAR_PIX }
            fn band_indices(&self) -> &BandIndices { &SENTINEL_BAND_INDICES }
            fn tauray(&self) -> &[f64] {
                todo!("Sentinel Rayleigh values loaded from LUT")
            }
            fn lambda(&self) -> &[f64] { &LAMBDA_SENTINEL }
        }
    };
}

impl_sentinel_sensor!(Sentinel2A, "SENTINEL_2A");
impl_sentinel_sensor!(Sentinel2B, "SENTINEL_2B");
impl_sentinel_sensor!(Sentinel2C, "SENTINEL_2C");

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    #[allow(unused_imports)]
    use crate::constants::*;

    #[test]
    fn test_landsat8_band_count() {
        let sensor = Landsat8;
        assert_eq!(sensor.reflectance_bands().len(), 8);
        assert_eq!(sensor.thermal_bands().len(), 2);
    }

    #[test]
    fn test_landsat8_resolution() {
        let sensor = Landsat8;
        assert_eq!(sensor.output_resolution_m(), 30.0);
    }

    #[test]
    fn test_landsat8_aerosol_window() {
        let sensor = Landsat8;
        assert_eq!(sensor.aerosol_window(), LAERO_WINDOW);
        assert_eq!(sensor.half_aerosol_window(), LHALF_AERO_WINDOW);
    }

    #[test]
    fn test_landsat8_band_indices() {
        let sensor = Landsat8;
        let idx = sensor.band_indices();
        assert_eq!(idx.coastal, 0);
        assert_eq!(idx.blue, 1);
        assert_eq!(idx.red, 3);
        assert_eq!(idx.nir, 4);
        assert_eq!(idx.swir2, 6);
    }

    #[test]
    fn test_landsat8_tauray() {
        let sensor = Landsat8;
        assert_eq!(sensor.tauray().len(), 8);
        assert!((sensor.tauray()[0] - 0.23638).abs() < 1e-5);
    }

    #[test]
    fn test_landsat9_same_config() {
        let l8 = Landsat8;
        let l9 = Landsat9;
        assert_eq!(l9.reflectance_bands().len(), l8.reflectance_bands().len());
        assert_eq!(l9.output_resolution_m(), l8.output_resolution_m());
        assert_eq!(l9.name(), "LANDSAT_9");
    }

    #[test]
    fn test_sentinel2a_band_count() {
        let sensor = Sentinel2A;
        assert_eq!(sensor.reflectance_bands().len(), 11);
        assert_eq!(sensor.thermal_bands().len(), 0);
    }

    #[test]
    fn test_sentinel2a_resolution() {
        let sensor = Sentinel2A;
        assert_eq!(sensor.output_resolution_m(), 10.0);
    }

    #[test]
    fn test_sentinel2a_aerosol_window() {
        let sensor = Sentinel2A;
        assert_eq!(sensor.aerosol_window(), SAERO_WINDOW);
        assert_eq!(sensor.half_aerosol_window(), 3);
    }

    #[test]
    fn test_sentinel2a_band_indices() {
        let sensor = Sentinel2A;
        let idx = sensor.band_indices();
        assert_eq!(idx.coastal, 0);
        assert_eq!(idx.blue, 1);
        assert_eq!(idx.red, 3);
        assert_eq!(idx.nir, 8);
        assert_eq!(idx.swir1, 9);
        assert_eq!(idx.swir2, 10);
    }

    #[test]
    fn test_sentinel2a_lambda() {
        let sensor = Sentinel2A;
        assert_eq!(sensor.lambda().len(), 11);
        assert!((sensor.lambda()[0] - 0.443).abs() < 1e-5);
    }

    #[test]
    fn test_sentinel_names() {
        assert_eq!(Sentinel2A.name(), "SENTINEL_2A");
        assert_eq!(Sentinel2B.name(), "SENTINEL_2B");
        assert_eq!(Sentinel2C.name(), "SENTINEL_2C");
    }

    #[test]
    fn test_num_refl_bands_default_impl() {
        assert_eq!(Landsat8.num_refl_bands(), 8);
        assert_eq!(Sentinel2A.num_refl_bands(), 11);
    }
}
