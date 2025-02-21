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
from shapely.geometry import LineString, Point

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
                        dp = DrainagePoint(
                            i=id_i,
                            j=id_j,
                            Z=float(elevation),
                            id_pnt=len(self.dr_pt) + 1
                            )
                        self.dr_pt.append(dp)
                        # Place the point in the matrix.
                        self.mat_id[id_i * 2, id_j * 2] = dp


                # Build the list 'qoi' by sorting DrainagePoints in ascending Z
                # qoi will store the 'id_pnt' for each DrainagePoint
                self.qoi = [dp.id_pnt.value for dp in sorted(self.dr_pt, key=lambda dp: dp.Z)]

                
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
        Exporte le réseau de drainage (slopelines) sous forme de Shapefile.
    
        Parameters
        ----------
        output_shapefile : str
            Chemin pour enregistrer le fichier Shapefile.
        """
        lines = []
        attributes = []
        
        for net in self.l_dr_net:
            coords = []
            for pnt_id in net.id_pnts.value:
                dp = self.dr_pt[pnt_id - 1]
                if dp:
                    x, y = self.transform * (dp.j, dp.i)
                    coords.append((x + self.delta_x / 2, y - self.delta_y / 2))
            
            if len(coords) > 1:
                lines.append(LineString(coords))
    
                # Nettoyage et conversion des valeurs
                attributes.append({
                    "id_ch": int(net.id_ch.value) if net.id_ch.value is not None else 0,
                    "nel": int(net.nel),
                    "id_pnts": str([int(p) if p is not None else 0 for p in net.id_pnts.value])[:254],
                    "id_start_pt": int(net.id_start_pt.value) if net.id_start_pt.value is not None else 0,
                    "id_end_pt": int(net.id_end_pt.value) if net.id_end_pt.value is not None else 0,
                    "length": float(net.length),
                    "id_ch_out": int(net.id_ch_out.value) if net.id_ch_out.value is not None else 0,
                    "n_jun": int(net.n_jun),
                    "id_in": [int(i.value) if i.value is not None else 0 for i in net.id_in.value] if net.id_in is not None else [],
                    "n_path": int(net.n_path),
                    "id_path": str([int(i) if i is not None else 0 for i in net.id_path] if net.id_path is not None else [])[:254],
                    "id_endo": int(net.id_endo.value) if net.id_endo.value is not None else -1,
                    "sso": int(net.sso) if net.sso is not None else 0,
                    "hso": int(net.hso) if net.hso is not None else 0
                })
    
        # Création du GeoDataFrame et exportation en fichier Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=lines, crs=self.crs)
        gdf.to_file(output_shapefile)
        
        
    def export_drainage_point(self, output_shapefile: str):
        """
        Exporte les points de drainage sous forme de Shapefile.
    
        Parameters
        ----------
        output_shapefile : str
            Chemin pour enregistrer le fichier Shapefile.
        """
        points = []
        attributes = []
        
        for dp in self.dr_pt:
            x, y = self.transform * (dp.j, dp.i)
            points.append(Point(x + self.delta_x / 2, y - self.delta_y / 2))
            
            # Nettoyage des valeurs avant d'exporter
            attributes.append({
                "id_pnt": int(dp.id_pnt.value) if dp.id_pnt.value is not None else 0,
                "i": int(dp.i),
                "j": int(dp.j),
                "Z": float(dp.Z),
                "fldir": int(dp.fldir.value) if dp.fldir.value is not None else 0,
                "fldir_ss": int(dp.fldir_ss.value) if dp.fldir_ss.value is not None else 0,
                "A_in": int(dp.A_in),
                "upl": float(dp.upl),
                "dpl": float(dp.dpl),
                "sumdev": float(dp.sumdev),
                "id_endo": int(dp.id_endo.value) if dp.id_endo.value is not None else 0,
                "ninf": int(dp.ninf) if dp.ninf is not None else 0,
                "inflow": str(dp.inflow.value) if dp.inflow.value is not None else "NULL",
                "Linflow": str(dp.Linflow.value) if dp.Linflow.value is not None else "NULL",
                "Sinflow": str(dp.Sinflow.value) if dp.Sinflow.value is not None else "NULL",
                "id_ch": int(dp.id_ch.value) if dp.id_ch.value is not None else 0
            })
        
        # Créer un GeoDataFrame et écrire en fichier Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=points, crs=self.crs)
        gdf.to_file(output_shapefile)        
            
            
    def export_ridge_point(self, output_shapefile: str):
        """
        Exporte les points de crête sous forme de Shapefile.
    
        Parameters
        ----------
        output_shapefile : str
            Chemin pour enregistrer le fichier Shapefile.
        """
        points = []
        attributes = []
        
        for rp in self.rd_pt:
            x, y = self.transform * (rp.j / 2, rp.i / 2)
            points.append(Point(x + self.delta_x / 2, y - self.delta_y / 2))
    
            # Nettoyage et conversion des valeurs
            attributes.append({
                "id_pnt": int(rp.id_pnt) if rp.id_pnt is not None else 0,
                "i": int(rp.i),
                "j": int(rp.j),
                "Z": float(rp.Z),
                "md": min(float(rp.md), 99999999) if rp.md is not None else 0.0,  # Limite de stockage
                "id_drpt1": int(rp.id_drpt1.value) if rp.id_drpt1 and rp.id_drpt1.value is not None else 0,
                "id_drpt2": int(rp.id_drpt2.value) if rp.id_drpt2 and rp.id_drpt2.value is not None else 0,
                "A_in": int(rp.A_in) if rp.A_in is not None else 0,
                "A_in_min": int(rp.A_in_min) if rp.A_in_min is not None else 0,
                "nen": int(rp.nen) if rp.nen is not None else 0,  
                "n_jun": int(rp.n_jun) if rp.n_jun is not None else 0,
                "id_neigh": str(rp.id_neigh.value) if rp.id_neigh and rp.id_neigh.value is not None else 0,
                "id_sdl": int(rp.id_sdl.value) if rp.id_sdl and rp.id_sdl.value is not None else 0,
                "nrdl": int(rp.nrdl) if rp.nrdl is not None else 0,
                "id_rdl": str(rp.id_rdl.value) if rp.id_rdl and rp.id_rdl.value is not None else 0,
                "junc": int(rp.junc) if rp.junc is not None else 0,
                "n_ptsa": int(rp.n_ptsa) if rp.n_ptsa is not None else 0
            })
    
        # Créer un GeoDataFrame et écrire en fichier Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=points, crs=self.crs)
        gdf.to_file(output_shapefile)
    

if __name__ == "__main__":
    dtm_path = "../../QGIS/out/cordevole_extrait_minimum2_6.tif"
    
    
    loader = LoadData()
    loader.read_geotiff(dtm_path)
    metadata = loader.get_metadata()
    # print(metadata)
