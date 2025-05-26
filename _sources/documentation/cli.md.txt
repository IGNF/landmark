# Command Line Interface (CLI)

The command-line interface (CLI) for **LANDMARK** allows users to extract thalweg and ridge networks directly from Digital Elevation Models (DEMs) via the terminal. This document details how to use the CLI, including all commands, options, requirements, and practical examples.

## Usage

After installation, run the CLI with:

```bash
landmark -h
```

Use the `-h` option to display help and see all available arguments and their descriptions.

To process a DEM, use:

```bash
landmark <MNT_input_path> [options]
```

### Positional Arguments

- **MNT_input_path**:  
  Path to the input GeoTIFF file containing the DEM to be processed.

### Optional Arguments

You may customize the analysis with the following optional arguments:

- `--output_dir` (str):  
  Directory where the output files (GeoPackages) will be saved.  
  *(default: "../../outputs")*
- `--a_spread` (float):  
  Threshold for ridge dispersion area (in m²).  
  *(default: 1e5)*
- `--a_out` (float):  
  Threshold for drainage contributing area (in m²).  
  *(default: 1e5)*
- `--hso_th` (int):  
  Horton stream order threshold for line export.  
  *(default: 5)*
- `--curvature_slope` (flag):  
  Enable computation of slope and curvature.  
  *(default: False)*
- `--n_pts_calc_slope` (int):  
  Number of points to use for local slope estimation.  
  *(default: 5, min: 3)*
- `--no_data_values` (list of float):  
  List of values to be treated as NoData in the DEM (space-separated).  
  *(default: -9999 0)*

### Example

Process a DEM with default parameters:

```bash
landmark my_dem.tif
```

Change the output directory and use a custom drainage area threshold:

```bash
landmark my_dem.tif --output_dir results --a_out 50000
```

Enable slope/curvature computation and specify NoData values:

```bash
landmark my_dem.tif --curvature_slope --no_data_values -9999 -32767
```

### Output

Two GeoPackage files will be created in the specified output directory:
- `slopelines_se_HSO_<input_name>.gpkg` – Extracted thalweg (slopeline) network.
- `ridgelines_se_HSO_<input_name>.gpkg` – Extracted ridge network.

These vector files can be directly opened in QGIS or any GIS software for visualization and analysis.

---

## Input DEM requirements

The input Digital Elevation Model must:
- Be a **GeoTIFF** raster file.
- Use a **projected coordinate system with metric units** (e.g., UTM, Lambert 93).  
  **Do not use a geographic CRS in degrees** (e.g., WGS84/EPSG:4326).  
  The algorithm relies on metric distances for all calculations (slope, curvature, contributing area, network thresholds). Using degrees will result in incorrect or unusable outputs.

---

## Error Handling

If the input file does not exist or is invalid, the CLI will raise an error:

```bash
FileNotFoundError: Input file not found: my_dem.tif
```

If an invalid value is provided for a parameter, an error will be raised:

```bash
ValueError: a_spread must be non-negative.
```

Or, for too few points in slope estimation:

```bash
ValueError: n_pts_calc_slope must be at least 3.
```

---

## Development Notes

The CLI is implemented in the `cli_landmark_geotiff.py` module (or as the package entry point).  
It validates user inputs, parses command-line arguments, and calls the `landmark_processing` function from `main_landmark_geotiff.py`.  
For detailed algorithmic explanations, please refer to the code documentation and the reference article by Moretti & Orlandini (2023).
