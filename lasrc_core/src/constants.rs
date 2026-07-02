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

// Landsat output-band index for the cirrus band (B9). Not atmospherically
// corrected; the TOA reflectance is copied through instead, matching C's
// SRL_BAND9 handling in lasrc.c.
pub const SRL_BAND9: usize = 7;

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

// Fill/nodata value written into output SR and BT bands for fill pixels.
pub const SR_FILL_VALUE: u16 = 0;

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

// Surface reflectance threshold arrays for subaeroret_new (per band)
// C: landsat_tth[NSRL_BANDS] — bands 1,2,3,4,5,6,7,9
pub const LANDSAT_TTH: [f64; NSRL_BANDS] = [
    1.0e-3, 1.0e-3, 0.0, 1.0e-3, 0.0, 0.0, 1.0e-4, 0.0,
];
pub const LANDSAT_TTH_WATER: [f64; NSRL_BANDS] = [
    1.0e-3, 1.0e-3, 0.0, 1.0e-3, 1.0e-3, 0.0, 1.0e-4, 0.0,
];

// Band wavelengths (micrometers)
pub const LAMBDA_LANDSAT: [f64; NREFLL_BANDS] = [0.443, 0.480, 0.585, 0.655, 0.865, 1.61, 2.2];
pub const LAMBDA_SENTINEL: [f64; 11] = [
    0.443, 0.490, 0.560, 0.665, 0.705, 0.740, 0.783, 0.842, 0.865, 1.61, 2.19,
];

// Gas transmission coefficients per Landsat band (from 6S model, gascoef-ldcm.ASC)
// Array order: bands 1,2,3,4,5,6,7,9
pub const OZTRANSA_LANDSAT: [f64; NSRL_BANDS] = [
    -0.00255649, -0.0177861, -0.0969872, -0.0611428,
    0.0001, 0.0001, 0.0001, -0.0834061,
];
pub const WVTRANSA_LANDSAT: [f64; NSRL_BANDS] = [
    2.29849e-27, 2.29849e-27, 0.00194772, 0.00404159,
    0.000729136, 0.00067324, 0.0177533, 0.00279738,
];
pub const WVTRANSB_LANDSAT: [f64; NSRL_BANDS] = [
    0.999742, 0.999742, 0.775024, 0.774482,
    0.893085, 0.939669, 0.65094, 0.759952,
];
pub const OGTRANSA1_LANDSAT: [f64; NSRL_BANDS] = [
    4.91586e-20, 4.91586e-20, 4.91586e-20, 1.04801e-05,
    1.35216e-05, 0.0205425, 0.0256526, 0.000214329,
];
pub const OGTRANSB0_LANDSAT: [f64; NSRL_BANDS] = [
    0.000197019, 0.000197019, 0.000197019, 0.640215,
    -0.195998, 0.326577, 0.243961, 0.396322,
];
pub const OGTRANSB1_LANDSAT: [f64; NSRL_BANDS] = [
    9.57011e-16, 9.57011e-16, 9.57011e-16, -0.348785,
    0.275239, 0.0117192, 0.0616101, 0.04728,
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
pub const ANGLE_FILL: f64 = -999.0;

// Sentinel-2 Rayleigh optical depth per band (13 bands, from tauray-msi.ASC)
// Order: B01, B02, B03, B04, B05, B06, B07, B08, B8A, B09, B10, B11, B12
pub const TAURAY_SENTINEL: [f64; 13] = [
    0.23432, 0.15106, 0.09102, 0.04535, 0.03584, 0.02924, 0.02338, 0.01847,
    0.01560, 0.01092, 0.00243, 0.00128, 0.00037,
];

// Sentinel-2 band wavelengths (micrometers), all 13 bands
pub const LAMBDA_SENTINEL_ALL: [f64; 13] = [
    0.443, 0.490, 0.560, 0.665, 0.705, 0.740, 0.783, 0.842,
    0.865, 0.945, 1.375, 1.610, 2.190,
];

// Sentinel-2 band count (all 13 bands including B09/B10)
pub const NSRS_BANDS: usize = 13;

// Gas transmission coefficients per Sentinel band (from gascoef-msi.ASC)
// Order: B01, B02, B03, B04, B05, B06, B07, B08, B8A, B09, B10, B11, B12
pub const OZTRANSA_SENTINEL: [f64; 13] = [
    -0.00264691, -0.0272572, -0.0986512, -0.0500348, -0.0204295,
    -0.0108641, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
];
pub const WVTRANSA_SENTINEL: [f64; 13] = [
    2.29849e-27, 2.29849e-27, 0.000777307, 0.00361051, 0.0141249,
    0.0137067, 0.00410217, 0.0285871, 0.000390755, 0.00001, 0.01,
    0.000640155, 0.018006,
];
pub const WVTRANSB_SENTINEL: [f64; 13] = [
    0.999742, 0.999742, 0.891099, 0.754895, 0.75596, 0.763497, 0.74117,
    0.578722, 0.900899, 0.45818, 1.0, 0.943712, 0.647517,
];
pub const OGTRANSA1_SENTINEL: [f64; 13] = [
    4.91586e-20, 4.91586e-20, 4.91586e-20, 4.91586e-20, 5.3367e-06,
    4.91586e-20, 9.03583e-05, 1.64109e-09, 1.90458e-05, 4.91586e-20,
    7.62429e-06, 0.0212751, 0.0243065,
];
pub const OGTRANSB0_SENTINEL: [f64; 13] = [
    0.000197019, 0.000197019, 0.000197019, 0.000197019, -0.980313,
    0.000197019, 0.0265393, 1.0e-10, 0.0322844, 0.000197019, 0.000197019,
    0.000197019, 0.000197019,
];
pub const OGTRANSB1_SENTINEL: [f64; 13] = [
    9.57011e-16, 9.57011e-16, 9.57011e-16, 9.57011e-16, 1.33639,
    9.57011e-16, 0.0532256, 1.0e-10, -0.0219907, 9.57011e-16, -0.216849,
    0.0116062, 0.0604312,
];

// Sentinel-specific QA flag (shares bit 5 with IPFLAG_INTERP_WINDOW, different sensor)
pub const IPFLAG_FAILED_TMP: u8 = 5;

// Sentinel-2 band indices (for the 13-band array)
pub const DNS_BAND1: usize = 0;   // B01 - Coastal aerosol
pub const DNS_BAND2: usize = 1;   // B02 - Blue
pub const DNS_BAND3: usize = 2;   // B03 - Green
pub const DNS_BAND4: usize = 3;   // B04 - Red (aerosol reference)
pub const DNS_BAND5: usize = 4;   // B05
pub const DNS_BAND6: usize = 5;   // B06
pub const DNS_BAND7: usize = 6;   // B07
pub const DNS_BAND8: usize = 7;   // B08 - NIR
pub const DNS_BAND8A: usize = 8;  // B8A - NIR narrow (NDWI/NDVI)
pub const DNS_BAND9: usize = 9;   // B09 - Water vapor (tgo=1 bypass)
pub const DNS_BAND10: usize = 10; // B10 - Cirrus (copy TOA)
pub const DNS_BAND11: usize = 11; // B11 - SWIR1
pub const DNS_BAND12: usize = 12; // B12 - SWIR2 (NDWI, aerosol)

// Sentinel-2 surface reflectance threshold arrays (13 bands, PROC_ALL_BANDS)
// Used in subaeroret_new convergence loop to stop when roslamb < tth[ib]
// Land thresholds: B01=1e-3, B02=1e-3, B03=0, B04=1e-3, B05-B12=0 except B12=1e-4
pub const SENTINEL_TTH: [f64; 13] = [
    1.0e-03, 1.0e-03, 0.0, 1.0e-03,  // B01, B02, B03, B04
    0.0, 0.0, 0.0, 0.0,              // B05, B06, B07, B08
    0.0, 0.0, 0.0, 0.0,              // B8A, B09, B10, B11
    1.0e-04,                          // B12
];
// Water thresholds: B01=1e-3, B04=1e-3, B8A=1e-3, B12=1e-4
pub const SENTINEL_TTH_WATER: [f64; 13] = [
    1.0e-03, 0.0, 0.0, 1.0e-03,      // B01, B02, B03, B04
    0.0, 0.0, 0.0, 0.0,              // B05, B06, B07, B08
    1.0e-03, 0.0, 0.0, 0.0,          // B8A, B09, B10, B11
    1.0e-04,                          // B12
];

/// Landsat Collection 2 QA_PIXEL: bit 0 = designated fill.
/// Matches C code's `level1_qa_is_fill()`.
#[inline]
pub fn is_fill_pixel(qa: u16) -> bool {
    (qa & 1) == 1
}
