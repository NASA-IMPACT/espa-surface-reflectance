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
    let mut x = space_def.ul_corner_x + (samp as f64 * space_def.pixel_size[0]);
    let mut y = space_def.ul_corner_y - (line as f64 * space_def.pixel_size[1]);

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
    let cscaa: f32 = (-(xmus_f as f64) * xmuv_f as f64
        - cosxfi_f as f64
            * (1.0 - xmus_f as f64 * xmus_f as f64).sqrt()
            * (1.0 - xmuv_f as f64 * xmuv_f as f64).sqrt()) as f32;
    let cscaa = cscaa.clamp(-1.0, 1.0);
    let scaa: f32 = ((cscaa as f64).acos() * RAD2DEG) as f32;
    scaa as f64
}

#[cfg(test)]
mod tests {
    use super::*;

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
