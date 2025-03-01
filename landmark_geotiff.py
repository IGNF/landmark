# -*- coding: utf-8 -*-
"""
Transcript in python of the original landmark.f90

! The program LENDMARK, starting from the DTM,
! extracts drainage network by using D8-TD method.
! Calculates mutual distance for each ridgeline point based on downslope path
! lenth.
! Defines connections between endorheic basins, even nested, and esorheic basins.
! For each endorheic basin defines the outflow path from the lowest saddle spilling
! to an esorheic basin.
! Constructs a hierarchically ordered ridglines network.
"""

#Internal import
from load_data_geotiff import LoadData
from D8_LTD import SlopelineMixin
from slopeline_recode import calculate_slopelines
from dpl_recode import dpl
from mutual_dist_recode import mutual_dist
from endo_del_recode import endo_del
from saddle_spill_recode import saddle_spill
from ridge_point import find_ridge_neighbors


class HydroModel(LoadData, SlopelineMixin):
    pass


if __name__ == "__main__":

    dtm_path = "../../QGIS/out/cordevole_extrait_minimum2_6.tif"
    
    file_name = dtm_path.split("/")[-1]
    
    #Shapefiles path
    slopelines_shapefile_path = f"../../out_scripts_test_temp/slopelines_{file_name[:-4]}_test"
    slopelines_HSO_shapefile_path = f"../../out_scripts_test_temp/slopelines_HSO2_{file_name[:-4]}_test"

    drainage_points_shapefile_path = f"../../out_scripts_test_temp/drain_points_{file_name[:-4]}_test"
    drainage_points_HSO_shapefile_path = f"../../out_scripts_test_temp/drain_points_HSO_{file_name[:-4]}_test"

    
    ridge_points_shapefile_path = f"../../out_scripts_test_temp/ridge_points_{file_name[:-4]}_test"
    ridge_points_shapefile_path_avant_relation = f"../../out_scripts_test_temp/ridge_points__avant_relation{file_name[:-4]}_test"
    saddle_points_shapefile_path = f"../../out_scripts_test_temp/saddle_points_filtre_{file_name[:-4]}_test"
    
    #Landmarks option
    main_channel_choice =2 #area:0, length:1, hso:2


        
    model_geotiff = HydroModel()
    model_geotiff.main_channel_choice = main_channel_choice
    
    
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
    
    print("Export drainage points in shapefile")
    # model_geotiff.export_drainage_point(drainage_points_shapefile_path)

    
    print("Export saddle points in shapefile")
    # model_geotiff.export_saddle_points(saddle_points_shapefile_path)

    
    print("Export ridges points in shapefile")
    # model_geotiff.export_ridge_point(ridge_points_shapefile_path_avant_relation)

    print("Export slopelines to shapefile")
    # model_geotiff.export_slopelines_to_shapefile(slopelines_shapefile_path)

    
    print("Define the relationship between ridge points")
    find_ridge_neighbors(model_geotiff)
    
    print("Export drainage points HSO in shapefile")
    # model_geotiff.export_drainage_point(drainage_points_HSO_shapefile_path)

    print("Export slopelines HSO to shapefile")
    # model_geotiff.export_slopelines_to_shapefile(slopelines_HSO_shapefile_path)
    
    print("Export ridges points in shapefile")
    # model_geotiff.export_ridge_point(ridge_points_shapefile_path)

    

        
    
        