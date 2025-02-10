# -*- coding: utf-8 -*-
"""
Script to read and transform a digital terrain model file in .dat format into a numpy array and geotiff file
Function to do the reverse : geotiff to .dat file, with the header.dat associate
"""

#External import
import numpy as np
import re
from typing import Dict
import rasterio
from rasterio.transform import from_origin
from rasterio.crs import CRS

def parse_header(file_path: str) -> Dict[str, float]:
    """Parses the header.dat file to extract useful metadata.

    Parameters
    ----------
    file_path : str
        Path to the header.dat file.

    Returns
    -------
    Dict[str, float]
        A dictionary containing the extracted metadata, including grid spacing,
        DEM size, coordinates of the lower-left corner, and the NoData value.
    """
    metadata = {}

    with open(file_path, 'r') as file:
        for line in file:
            # Match and extract specific metadata lines
            if "Grid spacing along the x-direction" in line:
                metadata["grid_spacing_x"] = float(re.search(r"[-+]?[0-9]*\.?[0-9]+", line).group())
            elif "Grid spacing along the y-direction" in line:
                metadata["grid_spacing_y"] = float(re.search(r"[-+]?[0-9]*\.?[0-9]+", line).group())
            elif "Number of row" in line:
                metadata["num_rows"] = int(re.search(r"[-+]?[0-9]+", line).group())
            elif "Number of column" in line:
                metadata["num_cols"] = int(re.search(r"[-+]?[0-9]+", line).group())
            elif "X low left corner coordinate" in line:
                metadata["x_low_left"] = float(re.search(r"[-+]?[0-9]*\.?[0-9]+", line).group())
            elif "Y low left corner coordinate" in line:
                metadata["y_low_left"] = float(re.search(r"[-+]?[0-9]*\.?[0-9]+", line).group())
            elif "NoData" in line:
                metadata["nodata"] = float(re.search(r"[-+]?[0-9]*\.?[0-9]+", line).group())

    return metadata


def read_dat_dem(file_path: str, num_rows: int, num_cols: int, nodata: float) -> np.ndarray:
    """Reads the dtm.dat file and converts it into a NumPy array.

    Parameters
    ----------
    file_path : str
        Path to the dtm.dat file.
    num_rows : int
        Number of rows in the DEM grid (from header).
    num_cols : int
        Number of columns in the DEM grid (from header).
    nodata : float
        NoData value to handle missing data in the DEM.

    Returns
    -------
    np.ndarray
        A 2D NumPy array representing the DEM, with NoData values set to np.nan.
    """
    dem = np.loadtxt(file_path, dtype=float).reshape((num_rows, num_cols))
    dem[dem == nodata] = np.nan  # Replace NoData values with NaN
    return np.flipud(dem)


def export_to_geotiff(output_path: str, dem: np.ndarray, metadata: Dict[str, float], epsg_code: int = 6707):
    """Exports the DEM array to a GeoTIFF file.

    Parameters
    ----------
    output_path : str
        Path to save the GeoTIFF file.
    dem : np.ndarray
        The DEM array to export.
    metadata : Dict[str, float]
        Metadata dictionary containing grid spacing and corner coordinates.
    epsg_code : int, optional
        EPSG code of the projection (default is 6707).
    """
    transform = from_origin(
        metadata["x_low_left"],
        metadata["y_low_left"],
        metadata["grid_spacing_x"],
        -metadata["grid_spacing_y"]
    )

    with rasterio.open(
        output_path,
        "w",
        driver="GTiff",
        height=dem.shape[0],
        width=dem.shape[1],
        count=1,
        dtype=dem.dtype,
        crs=CRS.from_epsg(epsg_code),
        transform=transform,
        nodata=np.nan
    ) as dst:
        dst.write(dem, 1)



def geotiff_to_dat(geotiff_path: str, output_dat_path: str, output_header_path: str):
    """
    Converts a GeoTIFF file to .dat format and creates a corresponding header.dat file.

    Parameters
    ----------
    geotiff_path : str
        Path to the input GeoTIFF file.
    output_dat_path : str
        Path to save the .dat file.
    output_header_path : str
        Path to save the header.dat file.
    """
    with rasterio.open(geotiff_path) as src:
        dem = src.read(1)  # Read the first band
        nodata = src.nodata if src.nodata is not None else -9999.0
        dem[dem == nodata] = -9999.00  # Replace NoData values with -9999.00

        transform = src.transform
        grid_spacing_x = transform[0]
        grid_spacing_y = -transform[4]

        num_rows, num_cols = dem.shape
        x_low_left = transform[2]
        y_low_left = transform[5] - (grid_spacing_y * (num_rows ))


        # Write .dat file
        with open(output_dat_path, "w") as f:
            for row in dem:
                f.write(" ".join(map(lambda v: f"{v:.2f}" if not np.isnan(v) else f"{int(nodata)}", row)) + "\n")

        # Write header.dat file
        with open(output_header_path, "w") as f:
            f.write("""
------------------------------------------------------------------------------
! HEADER
------------------------------------------------------------------------------
Grid spacing along the x-direction = {grid_spacing_x:.2f}
Grid spacing along the y-direction = {grid_spacing_y:.2f}
Number of row, DEM size along the y-direction = {num_rows}
Number of column, DEM size along the x-direction = {num_cols}
X low left corner coordinate = {x_low_left:.8f}
Y low left corner coordinate = {y_low_left:.8f}
NoData = -9999.00
------------------------------------------------------------------------------
""".format(
                grid_spacing_x=grid_spacing_x,
                grid_spacing_y=grid_spacing_y,
                num_rows=num_rows,
                num_cols=num_cols,
                x_low_left=x_low_left,
                y_low_left=y_low_left,
            ))




if __name__ == "__main__":
    
    """.dat to geotiff"""
    # #file path
    # header_path  = "../../out_scripts/cordevole_extrait_coord_header.dat"
    # dat_dem_path = "../../out_scripts/cordevole_extrait_coord.dat"
    # tiff_dem_path = "../../out_scripts_test_temp/cordevole_extrait_cooord_from_dat.tif"
    
    
    # #header info extraction
    # header_info = parse_header(header_path)
    
    # #reading .dat dem and transform in numpy array
    # dem = read_dat_dem(dat_dem_path, header_info["num_rows"], header_info["num_cols"], header_info["nodata"])
    
    # #Export dem in geotiff format
    # export_to_geotiff(tiff_dem_path, dem, header_info, epsg_code=6707)
    
    
    """geotiff to .dat"""
    #File path
    tiff_dem_path = "../../QGIS/out/cordevole_extrait_minimum2_6.tif"
    dat_dem_path = "../../out_scripts/cordevole_extrait_minimum2_6.dat"
    dat_header_path = "../../out_scripts/cordevole_extrait_minimum2_6_header.dat"
    
    #Convert geotiff to .dat format
    geotiff_to_dat(tiff_dem_path, dat_dem_path, dat_header_path)


