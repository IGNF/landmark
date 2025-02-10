# -*- coding: utf-8 -*-
"""
Geotiff loader for Landmark algorithm
"""

#Extern import
import numpy as np
import rasterio
from typing import Dict
from tqdm import tqdm
import geopandas as gpd
from shapely.geometry import LineString

#intern imports
from data_structures import DrainagePoint


class LoadData:
    """
    Class to load a DEM from a GeoTIFF file using rasterio.
    Stores the DEM data and its metadata for further processing.
    """
    def __init__(self):
        self.dem = None  # DEM as a NumPy array
        self.transform = None  # Affine transform for georeferencing
        self.crs = None  # Coordinate reference system
        self.nodata = None  # NoData value
        self.delta_x = None # Grid spacing along the x-direction
        self.delta_y = None # Grid spacing along the y-direction
        self.N = None # Number of row, DEM size along the y-direction
        self.M = None # Number of column, DEM size along the x-direction
        self.xllcorner = None # X low left corner coordinate
        self.yllcorner = None # Y low left corner coordinate
        self.nodata = None # NoData value
        
        # Data structures for the DEM
        self.dr_pt = []  # List of drainage points (DrainagePoint objects)
        # Matrix for direct access to points (size: (N*2-1, M*2-1))
        # Points will be placed at indices (i*2, j*2) to allow for interpolation/neighborhood processing.
        self.mat_id = None



    def read_geotiff(self, geotiff_file: str):
        """
        Reads a DEM from a GeoTIFF file and stores its metadata.

        Parameters
        ----------
        geotiff_file : str
            Path to the GeoTIFF file.
        """
        try:
            with rasterio.open(geotiff_file) as src:
                self.dem = src.read(1)  # Read the first band
                self.transform = src.transform  # Store georeferencing transform
                self.crs = src.crs  # Store coordinate reference system
                self.nodata = -9999.00
                
                # Replace NoData values with NaN for easier processing
                self.dem = np.where(self.dem == self.nodata, np.nan, self.dem)
                
                self.delta_x = self.transform[0]
                self.delta_y = -self.transform[4]
                self.N, self.M = self.dem.shape
                self.xllcorner = self.transform[2]
                self.yllcorner = self.transform[5] - self.delta_y * self.N
                self.mat_id = np.full((self.N * 2 - 1, self.M * 2 - 1), None, dtype=object)

                
                for id_i in tqdm(range(self.N), desc="Loading DEM", unit="row"):
                    for id_j in range(self.M):
                        elevation = self.dem[id_i, id_j]
                        # Create a DrainagePoint (id_pnt starts at 1)
                        dp = DrainagePoint(i=id_i, j=id_j, Z=float(elevation), id_pnt=len(self.dr_pt) + 1)
                        self.dr_pt.append(dp)
                        # Place the point in the matrix.
                        self.mat_id[id_i * 2, id_j * 2] = dp


                # Sorts the list of drainage points by elevation (ascending order).
                # This sort is essential for further processing (e.g., flow calculations).
                self.dr_pt.sort(key=lambda dp: dp.Z)

                
        except Exception as e:
            raise RuntimeError(f"Error reading GeoTIFF file '{geotiff_file}': {e}")
            
 

    def get_metadata(self) -> Dict[str, any]:
        """
        Returns metadata of the loaded DEM.

        Returns
        -------
        Dict[str, any]
            Metadata dictionary containing transform, CRS, and NoData value.
        """
        return {
            "transform": self.transform,
            "crs": self.crs,
            "nodata": self.nodata
        }
    
    def export_slopelines_to_shapefile(self, output_shapefile: str):
        """
        Exports the computed drainage networks (slopelines) to a Shapefile.
    
        Parameters
        ----------
        output_shapefile : str
            Path to save the Shapefile.
        """
        lines = []
        ids = []
    
        for net in self.dr_net:
            coords = []
            for pnt_id in net.id_pnts:
                dp = self.dr_pt_by_id.get(pnt_id)
                if dp:
                    x, y = self.transform * (dp.j, dp.i)  # Convert to projected coordinates
                    coords.append((x, y))
            if len(coords) > 1:
                lines.append(LineString(coords))
                ids.append(net.id_ch)
    
        gdf = gpd.GeoDataFrame({"id_ch": ids, "geometry": lines}, crs=self.crs)
        gdf.to_file(output_shapefile)


if __name__ == "__main__":
    dtm_path = "../../QGIS/out/cordevole_extrait_minimum2_6.tif"
    
    
    loader = LoadData()
    loader.read_geotiff(dtm_path)
    metadata = loader.get_metadata()
    # print(metadata)
