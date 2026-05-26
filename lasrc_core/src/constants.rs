//! Physical constants, LUT dimensions, and thresholds from the C LaSRC code.

// Version
pub const SR_VERSION: &str = "3.5.1.0 (Collection 2)";

// Angle conversions
pub const DEG2RAD: f64 = 0.017453293;
pub const RAD2DEG: f64 = 57.29577951;

// Atmospheric pressure
pub const ATMOS_PRES_0: f64 = 1013.0;
pub const ONE_DIV_ATMOS_PRES_0: f64 = 0.000987166;
pub const ONE_DIV_8500: f64 = 0.000117647;

// Default aerosol values
pub const DEFAULT_AERO: f64 = 0.05;
pub const DEFAULT_EPS: f64 = 1.5;

// Angstrom coefficient categories
pub const LOW_EPS: f64 = 1.0;
pub const WATER_EPS: f64 = 1.5;
pub const MOD_EPS: f64 = 1.75;
pub const HIGH_EPS: f64 = 2.5;

// Aerosol thresholds
pub const LOW_AERO_THRESH: f64 = 0.015;
pub const AVG_AERO_THRESH: f64 = 0.03;

// Aerosol window sizes
pub const LAERO_WINDOW: usize = 3;
pub const LHALF_AERO_WINDOW: usize = 1;
pub const SAERO_WINDOW: usize = 6;

pub const LFIX_AERO_WINDOW: usize = 15;
pub const LHALF_FIX_AERO_WINDOW: usize = 7;
pub const LMIN_CLEAR_PIX: usize = 4;

pub const SFIX_AERO_WINDOW: usize = 45;
pub const SHALF_FIX_AERO_WINDOW: usize = 22;
pub const SMIN_CLEAR_PIX: usize = 8;

// LUT dimensions
pub const NPRES_VALS: usize = 7;
pub const NAOT_VALS: usize = 22;
pub const NSOLAR_VALS: usize = 8000;
pub const NSUNANGLE_VALS: usize = 22;
pub const NVIEW_ZEN_VALS: usize = 20;
pub const NSOLAR_ZEN_VALS: usize = 22;
pub const NCOEF: usize = 4;

// Derived LUT dimensions
pub const NAOT_X_NSOLAR: usize = NAOT_VALS * NSOLAR_VALS;
pub const NAOT_X_NSUNANGLE: usize = NAOT_VALS * NSUNANGLE_VALS;

// CMG/DEM/Ratio grid dimensions (0.05 degree resolution)
pub const CMG_NBLAT: usize = 3600;
pub const CMG_NBLON: usize = 7200;
pub const DEM_NBLAT: usize = 3600;
pub const DEM_NBLON: usize = 7200;
pub const RATIO_NBLAT: usize = 3600;
pub const RATIO_NBLON: usize = 7200;

// Band counts
pub const NREFLL_BANDS: usize = 7;
pub const NSRL_BANDS: usize = 8;
pub const NBANDL_THM: usize = 2;

// Reflectance valid range
pub const MIN_VALID_REFL: f64 = -0.2;
pub const MAX_VALID_REFL: f64 = 1.60;
pub const MIN_VALID_TH: f64 = 150.0;
pub const MAX_VALID_TH: f64 = 350.0;

// Output scale factors
pub const SCALE_FACTOR_REFL: f64 = 0.0000275;
pub const OFFSET_REFL: f64 = -0.20;
pub const MULT_FACTOR_REFL: f64 = 36363.636363636;
pub const BAND_OFFSET_REFL: f64 = 0.20;
pub const SCALE_FACTOR_TH: f64 = 0.00341802;
pub const OFFSET_TH: f64 = 149.0;
pub const MULT_FACTOR_TH: f64 = 292.668;
pub const BAND_OFFSET_TH: f64 = 149.0;
pub const SCALE_FACTOR_AERO: f64 = 0.001;
pub const MULT_FACTOR_AERO: f64 = 1000.0;
pub const AERO_FILL: i16 = -9999;

// AOT values at 550nm (22 values)
pub const AOT550NM: [f64; NAOT_VALS] = [
    0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60,
    0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.30,
    2.60, 3.00, 3.50, 4.00, 4.50, 5.00,
];

// Log of AOT values (for atmospheric reflectance interpolation)
pub const LOG_AOT550NM: [f64; NAOT_VALS] = [
    -4.605170186, -2.995732274, -2.302585093, -1.897119985,
    -1.609437912, -1.203972804, -0.916290732, -0.510825624,
    -0.223143551, 0.000000000, 0.182321557, 0.336472237,
    0.470003629, 0.587786665, 0.693157181, 0.832909123,
    0.955511445, 1.098612289, 1.252762969, 1.386294361,
    1.504077397, 1.609437912,
];

// Pressure table (millibars)
pub const TPRES: [f64; NPRES_VALS] = [1050.0, 1013.0, 900.0, 800.0, 700.0, 600.0, 500.0];

// Rayleigh optical depth per Landsat band
pub const TAURAY_LANDSAT: [f64; NSRL_BANDS] = [
    0.23638, 0.16933, 0.09070, 0.04827, 0.01563, 0.00129, 0.00037, 0.07984,
];

// Band wavelengths (micrometers)
pub const LAMBDA_LANDSAT: [f64; NREFLL_BANDS] = [0.443, 0.480, 0.585, 0.655, 0.865, 1.61, 2.2];
pub const LAMBDA_SENTINEL: [f64; 11] = [
    0.443, 0.490, 0.560, 0.665, 0.705, 0.740, 0.783, 0.842, 0.865, 1.61, 2.19,
];

// LUT angle parameters
pub const XTS_MIN: f64 = 0.0;
pub const XTS_STEP: f64 = 4.0;
pub const XTV_MIN: f64 = 2.84090;
pub const XTV_STEP: f64 = 3.68017;

// QA flag bit positions
pub const IPFLAG_FILL: u8 = 0;
pub const IPFLAG_CLEAR: u8 = 1;
pub const IPFLAG_WATER: u8 = 2;
pub const IPFLAG_FAILED: u8 = 3;
pub const IPFLAG_FIXED: u8 = 4;
pub const IPFLAG_INTERP_WINDOW: u8 = 5;
pub const AERO1_QA: u8 = 6;
pub const AERO2_QA: u8 = 7;

// Sentinel-specific aerosol window constants
pub const HALF_EXPAND_WIN: usize = 12;
pub const HALF_FAILED_WIN: usize = 30;
pub const MIN_VALID_WINDOW_PIX: usize = 20;

// Fill values
pub const INPUT_FILL: u16 = 0;
pub const ANGLE_FILL: f64 = -999.0;
