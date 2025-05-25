# LANDMARK User Manual

## Introduction

This document provides a detailed user manual for **LANDMARK**, a Python package for the extraction of meaningful ridge and thalweg (valley line) networks from high-resolution Digital Elevation Models (DEMs). LANDMARK enables researchers and practitioners to analyze and partition terrain based on the intrinsic structure of the land surface, supporting applications in hydrology, geomorphology, and environmental modeling.

---

## What Does This Program Do?

LANDMARK identifies and extracts the main ridge and thalweg networks from an input DEM without modifying the original elevation values. It leverages advanced slopeline analysis to provide a physically meaningful representation of convergent (thalweg) and divergent (ridge) landscape features. These outputs are suitable for hydrological modeling, basin partitioning, and geomorphological interpretation.

---

## DEM Input Requirements

- The input must be a **GeoTIFF** file containing elevation values.
- The DEM must use a **projected coordinate system with metric units** (such as UTM or Lambert 93).  
  **Do not use a DEM in geographic (degree) coordinates (e.g., WGS84/EPSG:4326).**  
  The algorithm relies on metric distances for all calculations (e.g., slope, curvature, contributing area).
- The method is best suited to landscapes with **significant relief**. Very flat areas (such as river deltas or wide floodplains) may yield less relevant results, as the hydrological flow model cannot be computed robustly in such settings.

---

## Program Status

This program is currently provided as a **demonstrator**:  
While it is fully functional and produces scientifically meaningful results, it has not yet been optimized for speed or memory usage. Users can expect moderate to high resource consumption, especially with large DEMs. Future optimizations will be considered if there is interest from the user community.

---

## Main Processing Function

The primary entry point for LANDMARK is the `landmark_processing` function, defined in `main_landmark_geotiff.py`. This function runs the full processing pipeline on a single DEM.

### Function Signature

```python
def landmark_processing(
    mnt_path: str,
    output_dir: str = "../../outputs",
    a_spread: float = 1e5,
    a_out: float = 1e5,
    hso_th: int = 5,
    curvature_slope: bool = False,
    n_pts_calc_slope: int = 5,
    no_data_values: list = [-9999, 0]
) -> None:
    ...
```

### Parameters

- **mnt_path** (`str`):  
  Path to the input DEM (GeoTIFF format).
- **output_dir** (`str`, optional):  
  Directory where the output files (GeoPackages) will be saved. Default is `"../../outputs"`.
- **a_spread** (`float`, optional):  
  Threshold for the minimum dispersal area used to identify ridges (in square meters). Default is `1e5`.
- **a_out** (`float`, optional):  
  Threshold for the minimum drainage contributing area (in square meters). Default is `1e5`.
- **hso_th** (`int`, optional):  
  Horton stream order threshold for network extraction. Only lines of this order or higher are exported. Default is `5`.
- **curvature_slope** (`bool`, optional):  
  If `True`, computes local slope and curvature metrics. Default is `False`.
- **n_pts_calc_slope** (`int`, optional):  
  Number of points used to estimate slope. Must be at least `3`. Default is `5`.
- **no_data_values** (`list` of `int`/`float`, optional):  
  List of values in the DEM to treat as NoData (invalid). Default is `[-9999, 0]`.

### Output

Two files are produced in the output directory, both in GeoPackage format (`.gpkg`):

1. `slopelines_se_HSO_<input_name>.gpkg` – Contains the extracted thalweg (slopeline) network.
2. `ridgelines_se_HSO_<input_name>.gpkg` – Contains the extracted ridge network.

For details about the structure and attributes of these output files, see the next section.

---

## Output File Structure and Attributes

The structure and detailed explanation of the output files—including all attributes and practical tips for GIS visualization—are provided in a dedicated page:

[See: Outputs and Visualization Guide (outputs.md)](outputs.md)



