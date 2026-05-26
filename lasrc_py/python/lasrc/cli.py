"""Command-line interface for LaSRC surface reflectance processing."""

import click


@click.command()
@click.option("--input", "input_path", required=True, type=click.Path(exists=True),
              help="Path to input scene directory")
@click.option("--aux-dir", required=True, type=click.Path(exists=True),
              help="Path to auxiliary data directory")
@click.option("--output", "output_path", required=True, type=click.Path(),
              help="Path for output file or directory")
@click.option("--sensor", default=None,
              type=click.Choice(["LANDSAT_8", "LANDSAT_9", "SENTINEL_2A", "SENTINEL_2B", "SENTINEL_2C"]),
              help="Sensor name (auto-detected if omitted)")
@click.option("--output-format", default="cog", type=click.Choice(["cog", "espa"]),
              help="Output format (default: cog)")
@click.option("--use-orig-aero", is_flag=True, default=False,
              help="Use original aerosol algorithm (slower)")
@click.option("--aux-source", default="VIIRS", type=click.Choice(["VIIRS", "MODIS"]),
              help="Auxiliary data source (default: VIIRS)")
def main(input_path, aux_dir, output_path, sensor, output_format, use_orig_aero, aux_source):
    """LaSRC: Compute surface reflectance from satellite imagery."""
    from lasrc.pipeline import process_scene

    if sensor is None:
        sensor = _detect_sensor(input_path)

    click.echo(f"Processing {input_path} with sensor {sensor}")
    process_scene(
        input_path=input_path,
        aux_dir=aux_dir,
        output_path=output_path,
        sensor_name=sensor,
        output_format=output_format,
        use_orig_aero=use_orig_aero,
        aux_source=aux_source,
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
