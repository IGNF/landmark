# -*- coding: utf-8 -*-
"""
Transcript in python of the original landmark.f90

! The program LENDMARK, starting from the DTM,
! extracts drainage network by using D8-LTD method.
! Calculates mutual distance for each ridgeline point based on downslope path
! lenth.
! Defines connections between endorheic basins, even nested, and esorheic basins.
! For each endorheic basin defines the outflow path from the lowest saddle spilling
! to an esorheic basin.
! Constructs a hierarchically ordered ridglines network.
"""

from datetime import datetime
import sys
from pathlib import Path
import os
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

sys.setrecursionlimit(5000)

#Internal import
from landmark.geomorph_tools.load_data_geotiff import LoadData
from landmark.geomorph_tools.D8_LTD import SlopelineMixin
from landmark.geomorph_tools.slopeline import calculate_slopelines
from landmark.geomorph_tools.dpl import dpl
from landmark.geomorph_tools.mutual_dist import mutual_dist
from landmark.geomorph_tools.endo_del import endo_del
from landmark.geomorph_tools.saddle_spill import saddle_spill
from landmark.geomorph_tools.ridge_point import find_ridge_neighbors
from landmark.geomorph_tools.ridge_hier import ridge_hier


class HydroModel(LoadData, SlopelineMixin):
    pass


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
    """
    Full LANDMARK processing pipeline applied to a single GeoTIFF DEM.

    Parameters
    ----------
    mnt_path : str
        Path to the input DEM (GeoTIFF format).
    output_dir : str, optional
        Directory where output GeoPackages will be saved. Default is "../../outputs".
    a_spread : float, optional
        Threshold for ridge dispersion area (in m^2). Default is 1e5.
    a_out : float, optional
        Threshold for drainage contributing area (in m^2). Default is 1e5.
    hso_th : int, optional
        Horton stream order threshold for exporting lines. Default is 5.

    curvature_slope : bool, optional
        Whether to compute slope and curvature. Default is False.
    n_pts_calc_slope : int, optional
        Number of points for slope estimation. Default is 5.
    no_data_values : list, optional
        List of NoData values to consider as invalid. Default is [-9999, 0].
    """

    # Input validation
    if not os.path.isfile(mnt_path):
        raise FileNotFoundError(f"DEM file not found: {mnt_path}")
    if not isinstance(output_dir, str) or not output_dir:
        raise ValueError("Output directory must be a non-empty string.")
    if not isinstance(a_spread, (int, float)) or a_spread < 0:
        raise ValueError("a_spread must be a non-negative number.")
    if not isinstance(a_out, (int, float)) or a_out < 0:
        raise ValueError("a_out must be a non-negative number.")
    if not isinstance(hso_th, int) or hso_th < 0:
        raise ValueError("hso_th must be a non-negative integer.")
    if not isinstance(curvature_slope, bool):
        raise ValueError("curvature_slope must be a boolean value.")
    if not isinstance(n_pts_calc_slope, int) or n_pts_calc_slope < 3:
        raise ValueError("n_pts_calc_slope must be an integer >= 3.")
    if not isinstance(no_data_values, list) or not all(isinstance(v, (int, float)) for v in no_data_values):
        raise ValueError("no_data_values must be a list of numeric values.")

    file_name = os.path.splitext(os.path.basename(mnt_path))[0]
    os.makedirs(output_dir, exist_ok=True)

    start_time = datetime.now()
    print(f"[START] {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    model = HydroModel()
    model.main_channel_choice = 2
    model.type_of_landscape = 0
    model.a_spread_threshold = a_spread
    model.a_out_threshold = a_out
    model.hso_th = hso_th
    model.curvature_slope = curvature_slope
    model.n_pts_calc_slope = n_pts_calc_slope
    model.noData = no_data_values

    print("Loading DEM...")
    model.read_geotiff(mnt_path)

    print("Calculating slopelines...")
    calculate_slopelines(model)

    print("Calculating downslope paths...")
    dpl(model)

    print("Computing mutual distances...")
    mutual_dist(model)

    print("Delineating endorheic basins...")
    endo_del(model)

    print("Connecting basins via saddles...")
    saddle_spill(model)

    print("Building ridge point neighborhoods...")
    find_ridge_neighbors(model)

    print("Constructing hierarchical ridge network...")
    ridge_hier(model)

    print("Exporting ridgeline and slopeline segments...")
    model.export_slopelines_single_element(os.path.join(output_dir, f"slopelines_se_HSO_{file_name}.gpkg"))
    model.export_ridgelines_single_element(os.path.join(output_dir, f"ridgelines_se_HSO_{file_name}.gpkg"))

    end_time = datetime.now()
    print(f"[END] {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[INFO] Total elapsed time: {end_time - start_time}")


if __name__ == "__main__":

    start_time = datetime.now()
    print(f"[START] Processing started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    
    #Input DTM in geoTiff format
    dtm_path = ""
 
    
    file_name = dtm_path.split("/")[-1]
    
    #Output Geopackage path
    slopelines_HSO_path = f"../../outputs/slopelines_HSO_{file_name[:-4]}"
    slopelines_se_HSO_path = f"../../outputs/slopelines_se_HSO_{file_name[:-4]}"

    drainage_points_HSO_path = f"../../outputs/drain_points_HSO_{file_name[:-4]}"
    
    ridge_points_HSO_path = f"../../outputs/ridge_points_HSO{file_name[:-4]}"
    
    ridgelines_HSO_path = f"../../outputs/ridgelines_HSO_{file_name[:-4]}"
    ridgelines_se_HSO_path = f"../../outputs/ridgelines_se_HSO_{file_name[:-4]}"
    
    
    saddle_points_HSO_path = f"../../outputs/saddle_points_HSO_filtre_{file_name[:-4]}"
    
    endo_points_HSO_path  = f"../../outputs/endo_points_HSO{file_name[:-4]}"

    #NoData values in the DEM
    noData = [-9999, 0]
    
    #Landmarks option (for now, only the choices indicated are coded)
    main_channel_choice =2 #area:0, length:1, hso:2
    type_of_landscape = 0 #Natural drainage basin:0, Flood Plane:1
    
    #Export threshold value:
    A_spread = 1e5 #Ridge dispersion area
    A_out = 1e5 #Thalweg drainage area
    HSO_th = 5 # Horton stream order theshold
    
    #Calculate cuvature and slope
    calculate_curvature_slope = False
    
    #Calculate slope
    n_pts_calc_slope = 5 # Number of points to use to calculate the slope


    model_geotiff = HydroModel()
    model_geotiff.main_channel_choice = main_channel_choice
    model_geotiff.type_of_landscape = type_of_landscape
    model_geotiff.a_spread_threshold = A_spread
    model_geotiff.a_out_threshold = A_out
    model_geotiff.hso_th = HSO_th
    model_geotiff.curvature_slope = calculate_curvature_slope
    model_geotiff.n_pts_calc_slope = n_pts_calc_slope
    model_geotiff.noData = noData
    
    
    print("Loading data...")
    model_geotiff.read_geotiff(dtm_path)
    
    print("Calculating slopelines...")
    calculate_slopelines(model_geotiff)
    
        
    print("Calculating the length of the path between each DTM cell and the outflow point even if the basin is endorheic ")
    dpl(model_geotiff)
    
    print ("calculates the mutual distance between the two neighbor drainage points")
    mutual_dist(model_geotiff)
    
    
    print("Delineating endorheic basins")
    endo_del(model_geotiff)
        
    
    print("Connect basin by sadlle spill")
    saddle_spill(model_geotiff)
    
     
    print("Define the relationship between ridge points")
    find_ridge_neighbors(model_geotiff)

    print("Ridge hierarchization")
    ridge_hier(model_geotiff)
      
    
    #---------------Geopackage exports-------------------------------------
    print("Export slopelines single element HSO")
    model_geotiff.export_slopelines_single_element(slopelines_se_HSO_path)
    
    print("Export ridgelines single element HSO")
    model_geotiff.export_ridgelines_single_element(ridgelines_se_HSO_path)

    #--------------Export for debug only-----------------------------------
    # print("Export drainage points HSO")
    # model_geotiff.export_drainage_point(drainage_points_HSO_path)
    
    # print("\nExport endo points")
    # model_geotiff.export_endo_points(endo_points_HSO_path)

    # print("Export slopelines HSO")
    # model_geotiff.export_slopelines(slopelines_HSO_path)
    
    # print("Export saddle points")
    # model_geotiff.export_saddle_points(saddle_points_HSO_path)

    # print("\nExport ridges points HSO")
    # model_geotiff.export_ridge_point(ridge_points_HSO_path)
    
    # print("Export ridgelines HSO")
    # model_geotiff.export_ridgelines(ridgelines_HSO_path)
    
    end_time = datetime.now()
    duration = end_time - start_time
    print(f"[END] Processing finished at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[INFO] Total elapsed time: {duration}")


    