use numpy::{IntoPyArray, PyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use lasrc_core::constants::NSOLAR_ZEN_VALS;
use lasrc_core::correction::{AuxiliaryData, compute_surface_reflectance};
use lasrc_core::geometry::SpaceDef;
use lasrc_core::lut::LookupTables;
use lasrc_core::sensor::{Landsat8, Landsat9, Sentinel2A, Sentinel2B, Sentinel2C, Sensor};

#[pyclass]
struct PyLookupTables {
    inner: LookupTables,
}

#[pymethods]
impl PyLookupTables {
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
        rolutt: Vec<f64>,
        transt: Vec<f64>,
        sphalbt: Vec<f64>,
        normext: Vec<f64>,
        tsmax: Vec<f64>,
        tsmin: Vec<f64>,
        nbfic: Vec<f64>,
        nbfi: Vec<i32>,
        ttv: Vec<f64>,
        tts: Vec<f64>,
        nsr_bands: usize,
    ) -> PyResult<Self> {
        if tts.len() != NSOLAR_ZEN_VALS {
            return Err(PyValueError::new_err(format!(
                "tts must have {} elements, got {}",
                NSOLAR_ZEN_VALS,
                tts.len()
            )));
        }
        let mut tts_arr = [0.0f64; NSOLAR_ZEN_VALS];
        tts_arr.copy_from_slice(&tts);
        Ok(Self {
            inner: LookupTables {
                rolutt,
                transt,
                sphalbt,
                normext,
                tsmax,
                tsmin,
                nbfic,
                nbfi,
                ttv,
                tts: tts_arr,
                nsr_bands,
            },
        })
    }
}

#[pyclass]
struct PyAuxiliaryData {
    inner: AuxiliaryData,
}

#[pymethods]
impl PyAuxiliaryData {
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
        dem: Vec<i16>,
        wv: Vec<i16>,
        oz: Vec<i16>,
        ratiob1: Vec<i16>,
        ratiob2: Vec<i16>,
        ratiob7: Vec<i16>,
        intratiob1: Vec<i16>,
        intratiob2: Vec<i16>,
        intratiob7: Vec<i16>,
        slpratiob1: Vec<i16>,
        slpratiob2: Vec<i16>,
        slpratiob7: Vec<i16>,
        andwi: Vec<i16>,
        sndwi: Vec<i16>,
        wv_scale: f64,
        oz_scale: f64,
        wv_default: f64,
        oz_default: f64,
    ) -> Self {
        Self {
            inner: AuxiliaryData {
                dem,
                wv,
                oz,
                ratiob1,
                ratiob2,
                ratiob7,
                intratiob1,
                intratiob2,
                intratiob7,
                slpratiob1,
                slpratiob2,
                slpratiob7,
                andwi,
                sndwi,
                wv_scale,
                oz_scale,
                wv_default,
                oz_default,
            },
        }
    }
}

/// Compute surface reflectance for a full scene.
///
/// Parameters
/// ----------
/// sensor_name : str
///     One of "LANDSAT_8", "LANDSAT_9", "SENTINEL_2A", "SENTINEL_2B", "SENTINEL_2C"
/// toa_bands : list of 2-D float32 arrays
///     TOA reflectance per reflectance band.
/// bt_bands : list of 2-D float32 arrays
///     Brightness temperature bands (Landsat only; pass [] for Sentinel).
/// solar_zenith, solar_azimuth, view_zenith, view_azimuth : 2-D float32 arrays
///     Per-pixel angle grids (degrees).
/// qa_band : 2-D uint16 array
///     Input QA / fill mask (0 = fill).
/// lut : PyLookupTables
/// aux : PyAuxiliaryData
/// ul_corner_x, ul_corner_y : float
///     Upper-left corner coordinates in UTM meters.
/// pixel_size_x, pixel_size_y : float
///     Pixel size in meters (e.g. 30.0 for Landsat).
/// utm_zone : int
///     UTM zone number; negative for southern hemisphere.
/// use_orig_aero : bool
///
/// Returns
/// -------
/// dict with keys:
///   "sr_bands"  – list of 1-D int16 arrays (flattened), one per reflectance band
///   "bt_bands"  – list of 1-D uint16 arrays (flattened), one per thermal band
///   "aerosol"   – 1-D int16 array (flattened)
///   "qa"        – 1-D uint8 array (flattened)
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn process_surface_reflectance<'py>(
    py: Python<'py>,
    sensor_name: &str,
    toa_bands: Vec<PyReadonlyArray2<'py, f32>>,
    bt_bands: Vec<PyReadonlyArray2<'py, f32>>,
    solar_zenith: PyReadonlyArray2<'py, f32>,
    solar_azimuth: PyReadonlyArray2<'py, f32>,
    view_zenith: PyReadonlyArray2<'py, f32>,
    view_azimuth: PyReadonlyArray2<'py, f32>,
    qa_band: PyReadonlyArray2<'py, u16>,
    lut: &PyLookupTables,
    aux: &PyAuxiliaryData,
    ul_corner_x: f64,
    ul_corner_y: f64,
    pixel_size_x: f64,
    pixel_size_y: f64,
    utm_zone: i32,
    use_orig_aero: bool,
) -> PyResult<PyObject> {
    let sensor: Box<dyn Sensor> = match sensor_name {
        "LANDSAT_8" => Box::new(Landsat8),
        "LANDSAT_9" => Box::new(Landsat9),
        "SENTINEL_2A" => Box::new(Sentinel2A),
        "SENTINEL_2B" => Box::new(Sentinel2B),
        "SENTINEL_2C" => Box::new(Sentinel2C),
        _ => {
            return Err(PyValueError::new_err(format!(
                "Unknown sensor: {sensor_name}. Expected one of LANDSAT_8, LANDSAT_9, \
                 SENTINEL_2A, SENTINEL_2B, SENTINEL_2C"
            )))
        }
    };

    let toa_views: Vec<_> = toa_bands.iter().map(|a| a.as_array()).collect();
    let bt_views: Vec<_> = bt_bands.iter().map(|a| a.as_array()).collect();

    let space_def = SpaceDef {
        ul_corner_x,
        ul_corner_y,
        pixel_size: [pixel_size_x, pixel_size_y],
        zone: utm_zone,
    };

    let result = compute_surface_reflectance(
        sensor.as_ref(),
        &toa_views,
        &bt_views,
        &solar_zenith.as_array(),
        &solar_azimuth.as_array(),
        &view_zenith.as_array(),
        &view_azimuth.as_array(),
        &qa_band.as_array(),
        &lut.inner,
        &aux.inner,
        &space_def,
        use_orig_aero,
    );

    // Build output dict.
    // Array2<T> implements IntoPyArray directly, so no intermediate Vec needed.
    let dict = pyo3::types::PyDict::new(py);

    // SR bands: Vec<Array2<i16>> -> list of 1-D PyArray<i16>
    let sr_list: Vec<Bound<'py, PyArray1<i16>>> = result
        .sr_bands
        .into_iter()
        .map(|a| {
            let (v, _offset) = a.into_raw_vec_and_offset();
            v.into_pyarray(py)
        })
        .collect();
    dict.set_item("sr_bands", sr_list)?;

    // BT bands: Vec<Array2<u16>> -> list of 1-D PyArray<u16>
    let bt_list: Vec<Bound<'py, PyArray1<u16>>> = result
        .bt_bands
        .into_iter()
        .map(|a| {
            let (v, _offset) = a.into_raw_vec_and_offset();
            v.into_pyarray(py)
        })
        .collect();
    dict.set_item("bt_bands", bt_list)?;

    // Aerosol: Array2<i16> -> 1-D PyArray<i16>
    let (aerosol_vec, _) = result.aerosol.into_raw_vec_and_offset();
    dict.set_item("aerosol", aerosol_vec.into_pyarray(py))?;

    // QA: Array2<u8> -> 1-D PyArray<u8>
    let (qa_vec, _) = result.qa.into_raw_vec_and_offset();
    dict.set_item("qa", qa_vec.into_pyarray(py))?;

    Ok(dict.into())
}

#[pymodule]
fn lasrc(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", "0.1.0")?;
    m.add_class::<PyLookupTables>()?;
    m.add_class::<PyAuxiliaryData>()?;
    m.add_function(wrap_pyfunction!(process_surface_reflectance, m)?)?;
    Ok(())
}
