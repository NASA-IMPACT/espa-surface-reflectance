"""Command-line interface for LaSRC surface reflectance processing."""

import click


@click.command()
@click.option("--input", "input_path", required=True, type=click.Path(exists=True),
              help="Path to input scene directory")
@click.option("--output", "output_path", required=True, type=click.Path(),
              help="Path for output file or directory")
@click.option("--angle-hdf", required=True, type=click.Path(exists=True),
              help="Path to angle LUT HDF file (e.g. ANGLE_NEW.hdf)")
@click.option("--intref-hdf", required=True, type=click.Path(exists=True),
              help="Path to intrinsic reflectance LUT (e.g. RES_LUT_V3.0-LANDSAT.hdf)")
@click.option("--transm-hdf", required=True, type=click.Path(exists=True),
              help="Path to transmission LUT (e.g. TRANS_LUT_V3.0-LANDSAT.hdf)")
@click.option("--sphera-hdf", required=True, type=click.Path(exists=True),
              help="Path to spherical albedo LUT (e.g. AERO_LUT_V3.0-LANDSAT.hdf)")
@click.option("--wv-oz-hdf", required=True, type=click.Path(exists=True),
              help="Path to water vapor / ozone HDF file (date-specific)")
@click.option("--dem-hdf", required=True, type=click.Path(exists=True),
              help="Path to CMG DEM HDF file (e.g. CMGDEM.hdf)")
@click.option("--ratio-hdf", required=True, type=click.Path(exists=True),
              help="Path to band ratio / NDWI HDF file (e.g. ratiomapndwiexp.hdf)")
@click.option("--sensor", default=None,
              type=click.Choice(["LANDSAT_8", "LANDSAT_9", "SENTINEL_2A", "SENTINEL_2B", "SENTINEL_2C"]),
              help="Sensor name (auto-detected if omitted)")
@click.option("--output-format", default="cog", type=click.Choice(["cog", "espa"]),
              help="Output format (default: cog)")
@click.option("--aux-source", default="VIIRS", type=click.Choice(["VIIRS", "MODIS"]),
              help="Auxiliary data source (default: VIIRS)")
def main(input_path, output_path, angle_hdf, intref_hdf, transm_hdf, sphera_hdf,
         wv_oz_hdf, dem_hdf, ratio_hdf, sensor, output_format, aux_source):
    """LaSRC: Compute surface reflectance from satellite imagery."""
    from lasrc.pipeline import AuxFilePaths, process_scene

    if sensor is None:
        sensor = _detect_sensor(input_path)

    aux_files = AuxFilePaths(
        angle_hdf=angle_hdf,
        intref_hdf=intref_hdf,
        transm_hdf=transm_hdf,
        sphera_hdf=sphera_hdf,
        wv_oz_hdf=wv_oz_hdf,
        dem_hdf=dem_hdf,
        ratio_hdf=ratio_hdf,
        aux_source=aux_source,
    )

    click.echo(f"Processing {input_path} with sensor {sensor}")
    process_scene(
        input_path=input_path,
        aux_files=aux_files,
        output_path=output_path,
        sensor_name=sensor,
        output_format=output_format,
    )
    click.echo(f"Output written to {output_path}")


def _detect_sensor(input_path: str) -> str:
    """Auto-detect sensor from input scene metadata."""
    from pathlib import Path
    p = Path(input_path)
    name = p.name.upper()
    if "LC08" in name or "LO08" in name:
        return "LANDSAT_8"
    elif "LC09" in name or "LO09" in name:
        return "LANDSAT_9"
    elif "S2A" in name:
        return "SENTINEL_2A"
    elif "S2B" in name:
        return "SENTINEL_2B"
    elif "S2C" in name:
        return "SENTINEL_2C"
    else:
        raise click.ClickException(f"Cannot detect sensor from {name}. Use --sensor flag.")


if __name__ == "__main__":
    main()
