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
# from endo_del import endo_del


class HydroModel(LoadData, SlopelineMixin):
    pass


if __name__ == "__main__":

    dtm_path = "../../QGIS/out/cordevole_extrait_minimum2_6.tif"
    
    file_name = dtm_path.split("/")[-1]
    
    #Shapefiles path
    slopelines_shapefile_path = f"../../out_scripts_test_temp/slopelines_{file_name[:-4]}_test"
    drainage_points_shapefile_path = f"../../out_scripts_test_temp/drain_points_{file_name[:-4]}_test"
    ridge_points_shapefile_path = f"../../out_scripts_test_temp/ridge_points_{file_name[:-4]}_test"



        
    model_geotiff = HydroModel()    
    
    print("Loading data...")
    model_geotiff.read_geotiff(dtm_path)
    
    # print("First 5 drainage points:")
    # for point in model_geotiff.dr_pt[:5]:
    #     print("\n",point)
    
    
    print("Calculating slopelines...")
    calculate_slopelines(model_geotiff)
    
        
    print("Calculating the length of the path between each DTM cell and the outflow point even if the basin is endorheic ")
    dpl(model_geotiff)
    
    print("Export drainage points in shapefile")
    # model_geotiff.export_drainage_point(drainage_points_shapefile_path)

    
    print("Export slopelines to shapefile")
    # model_geotiff.export_slopelines_to_shapefile(slopelines_shapefile_path)

    
    # for pt in model.dr_pt[100:110]:
    #     print(f"Downslope path length of point {pt.id_pnt} is {pt.dpl} meters long")        
    
    
    print ("calculates the mutual distance between the two neighbor drainage points")
    mutual_dist(model_geotiff)
    
    print("Export ridges points in shapefile")
    model_geotiff.export_ridge_point(ridge_points_shapefile_path)
    
    
    # for pt in model.rd_pt[100:110]:
    #     print(pt)
        
    
    # print("Delineating endorheic basins")
    # endo_del(model)
    
    # print(model.endo_pt[10])
    
    
        