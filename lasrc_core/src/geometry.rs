//! Geometry utilities ported from the C LaSRC code.
//!
//! Includes UTM-to-degrees coordinate conversion (from `utmtodeg.c`) and
//! scattering angle computation.

use crate::constants::RAD2DEG;

/// Image space definition used for coordinate transformations.
pub struct SpaceDef {
    pub ul_corner_x: f64,
    pub ul_corner_y: f64,
    pub pixel_size: [f64; 2],
    /// UTM zone number; negative indicates southern hemisphere.
    pub zone: i32,
}

/// Convert a UTM image (line, sample) position to WGS84 latitude/longitude.
///
/// Returns `(lat_deg, lon_deg)`.
///
/// Algorithm ported directly from `utmtodeg.c` using WGS84 ellipsoid parameters.
pub fn utm_to_deg(space_def: &SpaceDef, line: i32, samp: i32) -> (f64, f64) {
    utm_to_deg_f(space_def, line as f64, samp as f64)
}

/// Convert fractional UTM image (line, sample) to WGS84 latitude/longitude using
/// the GCTP inverse transverse Mercator, as C reaches it through ESPA
/// `from_space` -> `inv_trans`.
///
/// The Landsat aerosol window lookup uses this, while the scene center and the
/// Sentinel window lookup use [`utm_to_deg`] (C's `utmtodeg`, since ESPA LaSRC is
/// not built with USE_GCTP). The two differ by ~1e-7 degrees, which is enough to
/// move the CMG interpolation weights.
///
/// Ported from `tm.c` (`inverse_transform`) and `gctp_utility.c` in ESPA's GCTP3,
/// with WGS84 axes from the `sphdz.c` spheroid table (index 12).
pub fn utm_to_deg_gctp(space_def: &SpaceDef, line: f64, samp: f64) -> (f64, f64) {
    const R_MAJOR: f64 = 6_378_137.0;
    const R_MINOR: f64 = 6_356_752.314_245;
    const SCALE_FACTOR: f64 = 0.9996;
    const EPSLN: f64 = 1.0e-10;
    const MAX_ITER: usize = 6;

    let zone = space_def.zone.abs();
    let lon_center = ((6 * zone) as f64 - 183.0).to_radians();
    let false_easting = 500_000.0;
    let false_northing = if space_def.zone < 0 { 10_000_000.0 } else { 0.0 };

    // Eccentricity constants (gctp_calc_e0..e3, gctp_calc_dist_from_equator)
    let temp = R_MINOR / R_MAJOR;
    let es = 1.0 - temp * temp;
    let e0 = 1.0 - 0.25 * es * (1.0 + es / 16.0 * (3.0 + 1.25 * es));
    let e1 = 0.375 * es * (1.0 + 0.25 * es * (1.0 + 0.46875 * es));
    let e2 = 0.05859375 * es * es * (1.0 + 0.75 * es);
    let e3 = es * es * es * (35.0 / 3072.0);
    // lat_origin is 0.0 for UTM, so ml0 = 0.0 (kept explicit to mirror the C)
    let lat_origin = 0.0f64;
    let ml0 = R_MAJOR
        * (e0 * lat_origin - e1 * (2.0 * lat_origin).sin() + e2 * (4.0 * lat_origin).sin()
            - e3 * (6.0 * lat_origin).sin());
    let esp = es / (1.0 - es);

    // from_space: projection coordinates for the (fractional) line/sample
    let x = space_def.ul_corner_x + samp * space_def.pixel_size[0] - false_easting;
    let y = space_def.ul_corner_y - line * space_def.pixel_size[1] - false_northing;

    let con = (ml0 + y / SCALE_FACTOR) / R_MAJOR;
    let mut phi = con;
    for i in 0.. {
        let delta_phi = ((con + e1 * (2.0 * phi).sin() - e2 * (4.0 * phi).sin()
            + e3 * (6.0 * phi).sin())
            / e0)
            - phi;
        phi += delta_phi;
        if delta_phi.abs() <= EPSLN || i >= MAX_ITER {
            break;
        }
    }

    let (sin_phi, cos_phi) = phi.sin_cos();
    let tan_phi = phi.tan();
    let c = esp * cos_phi * cos_phi;
    let cs = c * c;
    let t = tan_phi * tan_phi;
    let ts = t * t;
    let con = 1.0 - es * sin_phi * sin_phi;
    let n = R_MAJOR / con.sqrt();
    let r = n * (1.0 - es) / con;
    let d = x / (n * SCALE_FACTOR);
    let ds = d * d;

    let lat = phi
        - (n * tan_phi * ds / r)
            * (0.5
                - ds / 24.0
                    * (5.0 + 3.0 * t + 10.0 * c - 4.0 * cs - 9.0 * esp
                        - ds / 30.0
                            * (61.0 + 90.0 * t + 298.0 * c + 45.0 * ts - 252.0 * esp - 3.0 * cs)));
    let lon = lon_center
        + (d * (1.0
            - ds / 6.0
                * (1.0 + 2.0 * t + c
                    - ds / 20.0 * (5.0 - 2.0 * c + 28.0 * t - 3.0 * cs + 8.0 * esp + 24.0 * ts)))
            / cos_phi);
    // C: adjust_lon
    let lon = if lon.abs() < std::f64::consts::PI {
        lon
    } else {
        lon - lon.signum() * 2.0 * std::f64::consts::PI
    };

    (lat.to_degrees(), lon.to_degrees())
}

/// Like [`utm_to_deg`] but for fractional image coordinates.
pub fn utm_to_deg_f(space_def: &SpaceDef, line: f64, samp: f64) -> (f64, f64) {
    // WGS84 ellipsoid constants
    let sa: f64 = 6378137.0;
    let inv_flattening: f64 = 298.257223563;
    let false_easting: f64 = 500_000.0;
    let false_northing: f64 = 10_000_000.0;
    let scale_fact: f64 = 0.9996;

    // Derived ellipsoid quantities
    let sb = sa - (sa / inv_flattening);
    let e2 = ((sa * sa - sb * sb).sqrt()) / sb;
    let e2sq = e2 * e2;
    let c = (sa * sa) / sb;

    // Projection coordinates for the given line/sample
    let mut x = space_def.ul_corner_x + (samp * space_def.pixel_size[0]);
    let mut y = space_def.ul_corner_y - (line * space_def.pixel_size[1]);

    x -= false_easting;
    if space_def.zone < 0 {
        y -= false_northing;
    }

    let zone = space_def.zone.abs();
    let central_meridian = (zone as f64) * 6.0 - 183.0;

    // Initial latitude estimate from northing
    let mut lat = y / (6_366_197.724 * scale_fact);
    let cos_lat = lat.cos();
    let sqr_cos_lat = cos_lat * cos_lat;

    // Intermediate variables
    let v = (c / (1.0 + e2sq * sqr_cos_lat).sqrt()) * scale_fact;
    let a = x / v;
    let a1 = (2.0 * lat).sin();
    let a2 = a1 * sqr_cos_lat;
    let j2 = lat + a1 / 2.0;
    let j4 = (3.0 * j2 + a2) / 4.0;
    let j6 = (5.0 * j4 + a2 * sqr_cos_lat) / 3.0;
    let alpha = 0.75 * e2sq;
    let beta = (5.0 / 3.0) * alpha * alpha;
    let gama = (35.0 / 27.0) * alpha * alpha * alpha;
    let bm = scale_fact * c * (lat - alpha * j2 + beta * j4 - gama * j6);
    let b = (y - bm) / v;
    let epsi = e2sq * a * a / 2.0 * sqr_cos_lat;
    let eps = a * (1.0 - epsi / 3.0);
    let nab = b * (1.0 - epsi) + lat;
    let senoheps = (eps.exp() - (-eps).exp()) / 2.0;
    let delta = (senoheps / nab.cos()).atan();
    let ta0 = (delta.cos() * nab.tan()).atan();

    let lon = delta * RAD2DEG + central_meridian;
    lat = (lat
        + (1.0 + e2sq * sqr_cos_lat - 1.5 * e2sq * lat.sin() * cos_lat * (ta0 - lat))
            * (ta0 - lat))
        * RAD2DEG;

    (lat, lon)
}

/// Compute the scattering angle in degrees.
///
/// - `xmus`: cosine of solar zenith angle
/// - `xmuv`: cosine of view zenith angle
/// - `cosxfi`: cosine of relative azimuth angle
pub fn scattering_angle(xmus: f64, xmuv: f64, cosxfi: f64) -> f64 {
    // C: float cscaa = -xmus * xmuv - cosxfi * sqrt(1.0 - xmus*xmus) * sqrt(1.0 - xmuv*xmuv);
    // xmus, xmuv, cosxfi are float; sqrt is double; cscaa is float (truncated)
    // scaa = acos(cscaa) * RAD2DEG; — acos is double, scaa is float (truncated)
    let xmus_f = xmus as f32;
    let xmuv_f = xmuv as f32;
    let cosxfi_f = cosxfi as f32;
    // -xmus * xmuv and xmus * xmus are float products; the rest is double.
    let cscaa: f32 = ((-xmus_f * xmuv_f) as f64
        - cosxfi_f as f64
            * (1.0 - (xmus_f * xmus_f) as f64).sqrt()
            * (1.0 - (xmuv_f * xmuv_f) as f64).sqrt()) as f32;
    let cscaa = cscaa.clamp(-1.0, 1.0);
    let scaa: f32 = ((cscaa as f64).acos() * RAD2DEG) as f32;
    scaa as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_utm_to_deg_gctp_matches_c() {
        // Reference from ESPA GCTP3 tm.c inverse_transform compiled standalone,
        // for LC09_L1TP_035024_20250104 (zone 13N, UL 516585, 5846415, 30 m)
        // at the aerosol window center line 787, samp 1480 (C uses i-0.5, j+0.5).
        let space_def = SpaceDef {
            ul_corner_x: 516_585.0,
            ul_corner_y: 5_846_415.0,
            pixel_size: [30.0, 30.0],
            zone: 13,
        };
        let (lat, lon) = utm_to_deg_gctp(&space_def, 786.5, 1480.5);
        assert!((lat - 52.552_032_861_723).abs() < 1e-11, "lat={lat}");
        assert!((lon - -104.100_323_871_297).abs() < 1e-11, "lon={lon}");
    }

    #[test]
    fn test_utm_to_deg_zone_10n() {
        let space_def = SpaceDef {
            ul_corner_x: 545_000.0,
            ul_corner_y: 4_185_000.0,
            pixel_size: [30.0, 30.0],
            zone: 10,
        };
        let (lat, lon) = utm_to_deg(&space_def, 0, 0);
        assert!((lat - 37.79).abs() < 0.1, "lat={lat}");
        assert!((lon - (-122.40)).abs() < 0.1, "lon={lon}");
    }

    #[test]
    fn test_utm_to_deg_southern_hemisphere() {
        let space_def = SpaceDef {
            ul_corner_x: 500_000.0,
            ul_corner_y: 6_200_000.0,
            pixel_size: [30.0, 30.0],
            zone: -23,
        };
        let (lat, _lon) = utm_to_deg(&space_def, 0, 0);
        assert!(lat < 0.0, "Expected southern hemisphere, got lat={lat}");
    }

    #[test]
    fn test_utm_to_deg_pixel_offset() {
        let space_def = SpaceDef {
            ul_corner_x: 500_000.0,
            ul_corner_y: 4_000_000.0,
            pixel_size: [30.0, 30.0],
            zone: 11,
        };
        let (lat0, lon0) = utm_to_deg(&space_def, 0, 0);
        let (lat1, lon1) = utm_to_deg(&space_def, 100, 100);
        assert!(lat1 < lat0, "lat should decrease: {lat1} vs {lat0}");
        assert!(lon1 > lon0, "lon should increase: {lon1} vs {lon0}");
    }

    #[test]
    fn test_scattering_angle_backscatter() {
        let xmus = (0.5_f64).cos();
        let xmuv = (0.0_f64).cos();
        let cosxfi = 1.0;
        let sca = scattering_angle(xmus, xmuv, cosxfi);
        assert!(sca > 90.0, "Expected backscatter, got sca={sca}");
    }

    #[test]
    fn test_scattering_angle_range() {
        let sca = scattering_angle(0.7, 0.95, 0.5);
        assert!(sca >= 0.0 && sca <= 180.0, "sca={sca}");
    }
}
