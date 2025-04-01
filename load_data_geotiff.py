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
        a_unit = self.delta_x * self.delta_y
        
        for net in tqdm(self.dr_net):
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
                    "id_ch_main": int(net.id_path[-1]),
                    "A_out": int(self.dr_pt[net.id_pnts.value[-2]-1].A_in*a_unit),
                    "n_jun": int(net.n_jun),
                    "id_in": [int(i.value) if i.value is not None else 0 for i in net.id_in.value] if net.id_in is not None else [],
                    "n_path": int(net.n_path),
                    "id_path": str([int(i) if i is not None else 0 for i in net.id_path] if net.id_path is not None else [])[:254],
                    "id_endo": int(net.id_endo.value) if net.id_endo.value is not None else 0,
                    "sso": int(net.sso) if net.sso is not None else 0,
                    "hso": int(net.hso) if net.hso is not None else 0
                })
    
        # Création du GeoDataFrame et exportation en fichier Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=lines, crs=self.crs)
        gdf.to_file(output_shapefile)
        
    
    def export_slopelines_single_element_to_shapefile(self, output_shapefile: str):
        """
        Exporte le réseau de drainage (slopelines), tronçon par tronçon, sous forme de Shapefile.
    
        Parameters
        ----------
        output_shapefile : str
            Chemin pour enregistrer le fichier Shapefile.
        """
        lines = []
        attributes = []
        a_unit = self.delta_x * self.delta_y
    
        for net in tqdm(self.dr_net):
            points_ids = net.id_pnts.value
            for idx in range(len(points_ids) - 1):
                dp1 = self.dr_pt[points_ids[idx] - 1]
                dp2 = self.dr_pt[points_ids[idx + 1] - 1]
    
                if dp1 and dp2:
                    a_out = dp1.A_in * a_unit
                    if a_out >= self.a_out_threshold:
                        x1, y1 = self.transform * (dp1.j, dp1.i)
                        x2, y2 = self.transform * (dp2.j, dp2.i)
                        line = LineString([
                            (x1 + self.delta_x / 2, y1 - self.delta_y / 2),
                            (x2 + self.delta_x / 2, y2 - self.delta_y / 2)
                        ])
        
                        lines.append(line)
        
                        # Nettoyage et conversion des valeurs
                        attributes.append({
                            "id_ch": int(net.id_ch.value) if net.id_ch.value is not None else 0,
                            "start_pnt": int(dp1.id_pnt.value) if dp1.id_pnt.value is not None else 0,
                            "end_pnt": int(dp2.id_pnt.value) if dp2.id_pnt.value is not None else 0,
                            "id_ch_out": int(net.id_ch_out.value) if net.id_ch_out.value is not None else 0,
                            "id_ch_main": int(net.id_path[-1]),
                            "A_out": int(dp1.A_in * a_unit),
                            "length": float(line.length),
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
        
        for dp in tqdm(self.dr_pt):
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
        
        for rp in tqdm(self.rd_pt):
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
                "id_neigh": str(rp.id_neigh) if rp.id_neigh is not None else 0,
                "id_sdl": int(rp.id_sdl) if rp.id_sdl is not None else 0,
                "nrdl": int(rp.nrdl) if rp.nrdl is not None else 0,
                "id_rdl": str(rp.id_rdl) if rp.id_rdl is not None else 0,
                "junc": int(rp.junc) if rp.junc is not None else 0,
                "n_ptsa": float(rp.n_ptsa*0.5*self.delta_x*self.delta_y ) if rp.n_ptsa is not None else 0
            })
    
        # Créer un GeoDataFrame et écrire en fichier Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=points, crs=self.crs)
        gdf.to_file(output_shapefile)
   
    
    def export_saddle_points(self, output_shapefile: str):
        """
        Exporte les points de selle (`saddle points`) sous forme de Shapefile.
    
        Parameters
        ----------
        output_shapefile : str
            Chemin pour enregistrer le fichier Shapefile.
        """
        points = []
        attributes = []
        
        for sp in tqdm(self.sdl_pt):  # Parcours tous les Saddle Points
            if sp.id_cis_endo.value != None :
                x, y = self.transform * (sp.j / 2, sp.i / 2)
                points.append(Point(x + self.delta_x / 2, y - self.delta_y / 2))
        
                # Nettoyage et conversion des valeurs
                attributes.append({
                    "id_pnt": int(sp.id_pnt) if sp.id_pnt is not None else 0,
                    "id_rdpt": int(sp.id_rdpt) if sp.id_rdpt is not None else 0,
                    "id_rdpt2": int(sp.id_rdpt2) if sp.id_rdpt2 is not None else 0,
                    "Z": float(sp.Z),
                    "id_cis_endo": int(sp.id_cis_endo.value) if sp.id_cis_endo and sp.id_cis_endo.value is not None else 0,
                    "id_trans_out": int(sp.id_trans_out.value) if sp.id_trans_out and sp.id_trans_out.value is not None else 0,
                    "A_endo": float(sp.A_endo) if sp.A_endo is not None else 0.0
                })
    
        # Création du GeoDataFrame et exportation en fichier Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=points, crs=self.crs)
        gdf.to_file(output_shapefile) 


    def export_endo_points(self, output_shapefile: str):
        """
        Exporte les points endohériques sous forme de Shapefile.
    
        Parameters
        ----------
        output_shapefile : str
            Chemin pour enregistrer le fichier Shapefile.
        """
        points = []
        attributes = []
        
        for ep in tqdm(self.l_endo_pt):  # Parcours tous les Endo Points
            dp = self.dr_pt[ep.id_pnt.value-1]
    
            x, y = self.transform * (dp.j, dp.i)
            points.append(Point(x + self.delta_x / 2, y - self.delta_y / 2))
            
    
            # Nettoyage et conversion des valeurs
            attributes.append({
                "id_eo": int(ep.id_eo.value) if ep.id_eo.value is not None else 0,
                "id_pnt": int(ep.id_pnt.value) if ep.id_pnt.value is not None else 0,
                "bas_type": int(ep.bas_type) if ep.bas_type is not None else 0,
                "nsaddle": int(ep.nsaddle) if ep.nsaddle is not None else 0,
                "beyo_sad": str(ep.beyo_sad.value) if ep.beyo_sad.value is not None else "NULL",
                "idms": str(ep.idms.value) if ep.idms.value is not None else "NULL",
            })
    
        # Création du GeoDataFrame et exportation en fichier Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=points, crs=self.crs)
        gdf.to_file(output_shapefile) 
        
        
    def export_ridgelines_to_shapefile(self, output_shapefile: str):
        """
        Exporte les réseaux de crêtes (ridgelines) sous forme de Shapefile.
    
        Parameters
        ----------
        output_shapefile : str
            Chemin pour enregistrer le fichier Shapefile.
        """
        lines = []
        attributes = []
    
        for net in tqdm(self.rd_net):
            coords = []
            for pnt_id in net.id_pnts:
                rp = self.rd_pt[pnt_id - 1]  # RidgePoint
                if rp:
                    x, y = self.transform * (rp.j / 2, rp.i / 2)
                    coords.append((x + self.delta_x / 2, y - self.delta_y / 2))
    
            if len(coords) > 1:
                lines.append(LineString(coords))
    
                attributes.append({
                    "id_rdl": int(net.id_rdl),
                    "nel": int(net.nel),
                    "id_pnts": str(net.id_pnts)[:254],
                    "length": float(net.length),
                    "Zmean": float(net.Zmean),
                    "nrdpt_down": int(net.nrdpt_down),
                    "n_down": int(net.n_down),
                    "Zmean_down": float(net.Zmean_down),
                    "Z_diff": float(net.Zmean-net.Zmean_down) if net.Zmean_down > 0 else 0,
                    "A_in": float((self.rd_pt[net.id_pnts[0]-1].A_in+1)*self.delta_x*self.delta_y),
                    "A_in_min": float((self.rd_pt[net.id_pnts[0]-1].A_in_min+1)*self.delta_x*self.delta_y),
                    "jun_el": int(net.jun_el)
                })
    
        # Création du GeoDataFrame et exportation en fichier Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=lines, crs=self.crs)
        gdf.to_file(output_shapefile)


    def export_ridgelines_single_element_to_shapefile(self, output_shapefile: str):
        """
        Exporte les réseaux de crêtes (ridgelines), tronçon par tronçon, sous forme de Shapefile.
    
        Parameters
        ----------
        output_shapefile : str
            Chemin pour enregistrer le fichier Shapefile.
        """
        lines = []
        attributes = []
    
        for net in tqdm(self.rd_net):
            points_ids = net.id_pnts
            for idx in range(len(points_ids) - 1):
                rp1 = self.rd_pt[points_ids[idx] - 1]
                rp2 = self.rd_pt[points_ids[idx + 1] - 1]
    
                if rp1 and rp2:
                    A_last = (1 / (rp2.n_jun)) * self.delta_x * self.delta_y
                    a_spread = (rp1.n_ptsa*0.5*self.delta_x*self.delta_y ) + A_last 
                    
                    if a_spread >= self.a_spread_threshold:
                        x1, y1 = self.transform * (rp1.j / 2, rp1.i / 2)
                        x2, y2 = self.transform * (rp2.j / 2, rp2.i / 2)
                        line = LineString([
                            (x1 + self.delta_x / 2, y1 - self.delta_y / 2),
                            (x2 + self.delta_x / 2, y2 - self.delta_y / 2)
                        ])
        
                        lines.append(line)
        
                        attributes.append({
                            "id_rdl": int(net.id_rdl),
                            "start_pnt": int(rp1.id_pnt),
                            "end_pnt": int(rp2.id_pnt),
                            "length": float(line.length),
                            "jun_el": int(net.jun_el),
                            "A_spread": float(a_spread ) 
                        })
    
        # Création du GeoDataFrame et exportation en fichier Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=lines, crs=self.crs)
        gdf.to_file(output_shapefile)



if __name__ == "__main__":
    dtm_path = "../../QGIS/out/cordevole_extrait_minimum2_6.tif"
    
    
    loader = LoadData()
    loader.read_geotiff(dtm_path)
    metadata = loader.get_metadata()
    # print(metadata)
