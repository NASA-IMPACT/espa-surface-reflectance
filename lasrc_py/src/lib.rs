use pyo3::prelude::*;

#[pymodule]
fn lasrc(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", "0.1.0")?;
    Ok(())
}
