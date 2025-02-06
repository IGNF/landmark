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
from load_data import LoadData
from D8_LTD import SlopelineMixin
from slopeline import calculate_slopelines

import cProfile
import pstats

class HydroModel(LoadData, SlopelineMixin):
    pass


if __name__ == "__main__":

    # Paths to input files
    # header_path  = "../../cordevole/cordevole_extrait_coord_header.dat"
    # dtm_path = "../../cordevole/cordevole_extrait_coord.dat"
    
    header_path  = "../../out_scripts/cordevole_extrait_minimum2_6_header.dat"
    dtm_path = "../../out_scripts/cordevole_extrait_minimum2_6.dat"

    rs_path = "rs.dat"           # Optionnel : chemin vers rs.dat
    
    model = HydroModel()    
    
    print("Loading data...")
    model.process(header_path, dtm_path, rs_file=rs_path)
    
    # print("First 5 drainage points:")
    # for point in model.dr_pt[:5]:
    #     print(point)
    
    print("Calculating slopelines...")
    # calculate_slopelines(model)
    cProfile.run("calculate_slopelines(model)", "profile_results")
    
    
    # Affichage trié des résultats
    stats = pstats.Stats("profile_results")
    stats.strip_dirs().sort_stats("cumulative").print_stats(20)  # Afficher les 20 fonctions les plus lentes
    
    # print("\nDrainage Networks:")
    # for net in model.dr_net:
    #     print(f"Channel ID {net.id_ch}: {len(net.id_pnts)} points, length = {net.length}")
    
    # print("\nEndorheic Points:")
    # for endo in model.endo_pt:
    #     print(f"EndoPoint ID {endo.id_eo} at drainage point {endo.id_pnt}, basin type {endo.bas_type}")