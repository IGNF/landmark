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

import sys
sys.setrecursionlimit(5000)

#Internal import
from load_data_geotiff import LoadData
from D8_LTD import SlopelineMixin
from slopeline import calculate_slopelines
from dpl import dpl
from mutual_dist import mutual_dist
from endo_del import endo_del
from saddle_spill import saddle_spill
from ridge_point import find_ridge_neighbors
from ridge_hier import ridge_hier


class HydroModel(LoadData, SlopelineMixin):
    pass


if __name__ == "__main__":

    dtm_path = "../../QGIS/out/cordevole_extrait_minimum2_6.tif"
    # dtm_path = "../../QGIS/out/cordevole_extrait_coord.tif"
    # dtm_path = "../../QGIS/out/cordevole_extrait/cordevole_extrait_extrait.tif"
    # dtm_path = "../../QGIS/out/cordevole_debug_mini.tif"
    # dtm_path = "../../QGIS/out/cordevole_debug_riquiqui.tif"
    # dtm_path = "../../QGIS/out/cordevole_debug_mini_mini.tif"
    # dtm_path = "../../QGIS/out/cordevole_debug.tif"




    
    file_name = dtm_path.split("/")[-1]
    
    #Geopackage path
    slopelines_HSO_path = f"../../out_scripts_test_temp/slopelines_HSO_{file_name[:-4]}_test_nocython"
    slopelines_se_HSO_path = f"../../out_scripts_test_temp/slopelines_se_HSO_{file_name[:-4]}_test_nocython"

    drainage_points_HSO_path = f"../../out_scripts_test_temp/drain_points_HSO_{file_name[:-4]}_test"
    
    ridge_points_HSO_path = f"../../out_scripts_test_temp/ridge_points_HSO{file_name[:-4]}_test"
    
    ridgelines_HSO_path = f"../../out_scripts_test_temp/ridgelines_HSO_{file_name[:-4]}_test_nocython"
    ridgelines_se_HSO_path = f"../../out_scripts_test_temp/ridgelines_se_HSO_{file_name[:-4]}_test_nocython"
    
    
    saddle_points_HSO_path = f"../../out_scripts_test_temp/saddle_points_HSO_filtre_{file_name[:-4]}_test"
    
    endo_points_HSO_path  = f"../../out_scripts_test_temp/endo_points_HSO{file_name[:-4]}_test"

    
    #Landmarks option (for now, only the choices indicated are coded)
    main_channel_choice =2 #area:0, length:1, hso:2
    type_of_landscape = 0 #Natural drainage basin:0, Flood Plane:1
    
    #Export threshold value:
    A_spread = 0 #Ridge dispersion area
    A_out = 0 #Thalweg drainage area
    HSO_th = 0 # Horton stream order theshold
    

        
    model_geotiff = HydroModel()
    model_geotiff.main_channel_choice = main_channel_choice
    model_geotiff.type_of_landscape = type_of_landscape
    model_geotiff.a_spread_threshold = A_spread
    model_geotiff.a_out_threshold = A_out
    model_geotiff.hso_th = HSO_th
    
    
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
    
    
    
    print("Export drainage points HSO")
    model_geotiff.export_drainage_point(drainage_points_HSO_path)
    
    print("\nExport endo points")
    model_geotiff.export_endo_points(endo_points_HSO_path)


    print("Export slopelines HSO")
    model_geotiff.export_slopelines(slopelines_HSO_path)
    
    print("Export slopelines single element HSO")
    model_geotiff.export_slopelines_single_element(slopelines_se_HSO_path)

    
    print("Export saddle points")
    model_geotiff.export_saddle_points(saddle_points_HSO_path)

    
    print("\nExport ridges points HSO")
    model_geotiff.export_ridge_point(ridge_points_HSO_path)
    
    print("Export ridgelines HSO")
    model_geotiff.export_ridgelines(ridgelines_HSO_path)
    
    print("Export ridgelines single element HSO")
    model_geotiff.export_ridgelines_single_element(ridgelines_se_HSO_path)

    