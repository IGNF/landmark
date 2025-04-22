# DEM Generalization with LIC

[![License: MIT](./docs/_static/License-MIT.svg)](./LICENSE.txt)
[![Python](./docs/_static/Python-3.10.svg)](https://www.python.org/)


This Python package provides a python implementation of the method proposed by Moretti and Orlandini (2023) for extracting physically meaningful thalweg and ridge networks directly from high-resolution, unaltered digital elevation models (DEMs). The extracted networks form the basis for a new approach to terrain partitioning that preserves essential topographic features without relying on depression filling or grid coarsening. This work is a direct translation of the original Fortran implementation described in the published article :



```{note}
Moretti, G., & Orlandini, S. (2023). Thalweg and ridge network extraction from unaltered topographic data as a basis for terrain partitioning. Journal of Geophysical Research: Earth Surface, 128, e2022JF006943. 
```
[https://doi.org/10.1029/2022JF006943](https://doi.org/10.1029/2022JF006943)


## Documentation

Full documentation is available at:

[https://esaint-denis.github.io/landmark/](https://esaint-denis.github.io/landmark/)

## What does it do?

`landmark` processes high-resolution DEMs to identify and extract interconnected thalweg and ridge networks without modifying the original elevation data. These networks are derived from a slopeline-based analysis and provide a physically consistent representation of convergent (thalweg) and divergent (ridge) terrain structures.

Key functionalities include:
- Detection of **ridge points** as local high-divide locations uncrossed by slopelines.
- Identification of **endorheic** and **exorheic** basins, and their connections through **spilling saddles**.
- Construction of a fully connected **thalweg network**, including within nested depressions.
- Derivation of a **ridge network** structured around dispersal areas.
- Support for **terrain partitioning** using extracted networks as breaklines, enabling unstructured mesh generation for hydrological or geomorphological modeling.

This tool is particularly useful in mountainous or lowland terrains where preserving depressions is essential for a meaningful representation of surface processes.


![Exemple of Thalweg and Ridge Network extraction](docs/images/extraction_zone_O1.png)

## Installation

To install the package, clone the repository and install it in editable mode.

### HTTP method

```bash
git clone https://github.com/ESaint-Denis/landmark.git
cd landmark
pip install -e .
```

### SSH method

```bash
git clone git@github.com:ESaint-Denis/landmark.git
cd landmark
pip install -e .
```


## Python Quick Start

You can use **Landmark** directly in Python. The main script is "main_landmark_geotiff.py" or "main_landmark_geotiff_batch.py" for process many DEM geotiff files.
```

## License

This project is licensed under the MIT License. See the LICENSE.TXT file for details.