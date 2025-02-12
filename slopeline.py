# -*- coding: utf-8 -*-

"""
Module for constructing drainage networks (slopelines) based on the D8-LTD algorithm.
This module uses the previously loaded DEM data (via the model) and calls the D8_LTD algorithm.
"""

import time
from math import sqrt
import numpy as np
from tqdm import tqdm

# Import the data structures.
from data_structures import DrainageNetwork, EndoPoint



def calculate_slopelines(model):
    """
    Constructs drainage networks (slopelines) by processing each drainage point.
    For each point (processed in descending order of elevation), the function:
      - If the point has upstream inflows (ninf > 0), updates the corresponding drainage network.
      - If the point has no inflows, it is treated as a channel head and a new network is created.
      - Calls d8_ltd() for the point to compute its outflow direction and cumulative deviation.
      - Updates the outflow point's inflow information or, if no valid outflow is found,
        classifies the point as belonging to an endorheic basin.
      
    Parameters
    ----------
    model : object
        An instance of your data model class which must have at least:
           - dr_pt: list of DrainagePoint objects (sorted in ascending order of elevation)
           - mat_id: 2D numpy array with DrainagePoint objects at positions (i*2, j*2)
           - N, M: dimensions of the DEM (number of rows, number of columns)
           - delta_x, delta_y: grid spacings
           - nodata: no-data value for the DEM
    """
    start_time = time.process_time()
    
    # Initialize drainage network and endorheic lists.
    model.dr_net = []      # List of DrainageNetwork objects.
    model.endo_pt = []     # List of EndoPoint objects.
    model.endorheic_count = 0  # Counter for endorheic points.
    
    
    # Process drainage points in descending order (assuming model.dr_pt is sorted in ascending order)
    dr_pt_dict = {p.id_pnt: p for p in model.dr_pt}  
    dr_net_dict = {net.id_ch: net for net in model.dr_net}  
    
    for dp in tqdm(model.dr_pt[::-1], desc="Building drainage networks", unit="point"):
    
        if dp.ninf > 0:
            if dp.Linflow:
                l_max = max(dp.Linflow)
                i_max = np.argmax(dp.Linflow)  # Remplace index(l_max) par argmax()
                dp.upl = l_max
                dp.sumdev = dp.Sinflow[i_max] 
                id_up_max = dp.inflow[i_max] 
    
                up_dp = dr_pt_dict.get(id_up_max)
                if up_dp:
                    dp.id_ch = up_dp.id_ch
                    id_main = dp.id_ch
                else:
                    dp.id_ch = None
                    id_main = None
    
            for inflow_id in dp.inflow:
                inflow_dp = dr_pt_dict.get(inflow_id)
                if inflow_dp is None:
                    continue
    
                curr_ch = inflow_dp.id_ch
                net = dr_net_dict.get(curr_ch)
    
                if net is None:
                    # Créer un nouveau réseau si inexistant
                    net = DrainageNetwork(
                        id_ch=curr_ch,
                        id_pnts=[inflow_dp.id_pnt],
                        id_start_pt=inflow_dp.id_pnt,
                        id_end_pt=inflow_dp.id_pnt,
                        length=0.0,
                        id_ch_out=curr_ch,
                        n_jun=0,
                        id_in=[]
                    )
                    # Ajouter le réseau au dictionnaire et à la liste
                    dr_net_dict[curr_ch] = net
                    model.dr_net.append(net)
    
                net.id_pnts.append(dp.id_pnt)
                net.length = dp.upl if dp.Linflow else net.length
                net.id_ch_out = dp.id_ch
                net.id_end_pt = dp.id_pnt
    
                # Vérifier si le canal courant est un tributaire
                if id_main is not None and curr_ch != id_main:
                    main_net = dr_net_dict.get(id_main)
                    if main_net and curr_ch not in main_net.id_in:
                        main_net.id_in.append(curr_ch)
    
        else:
            # GESTION DES POINTS SANS AFFLUENTS
            new_channel_id = len(model.dr_net) + 1
            dp.id_ch = new_channel_id
            net = DrainageNetwork(
                id_ch=new_channel_id,
                nel = 1,
                id_pnts=[dp.id_pnt],
                id_start_pt=dp.id_pnt,
                id_end_pt=dp.id_pnt,
                length=0.0,
                id_ch_out=new_channel_id,
                n_jun=0,
                id_in=[]
            )
            dr_net_dict[new_channel_id] = net
            model.dr_net.append(net)

        
        # --- Compute drainage direction using D8_LTD ---
        # Call the D8_LTD function on the current drainage point.
        i_out, j_out, ndfl, sumdev = model.d8_ltd(dp)
        
        if i_out is not None and j_out is not None and i_out > 0 and j_out > 0:
            # Retrieve the drainage point corresponding to the outflow direction.
            out_dp = model.mat_id[i_out * 2, j_out * 2]
            if out_dp is not None:
                dp.fldir = out_dp.id_pnt
                dp.fldir_ss = None
                # Update outflow point: increment inflow count and record the current point.
                out_dp.ninf += 1
                out_dp.inflow.append(dp.id_pnt)
                # Compute the Euclidean distance between dp and out_dp.
                distance = sqrt(((dp.i - out_dp.i) * model.delta_x)**2 +
                                ((dp.j - out_dp.j) * model.delta_y)**2)
                out_dp.Linflow.append(dp.upl + distance)
                out_dp.Sinflow.append(sumdev)
                # out_dp.sumdev = sumdev
                
                if dp.id_pnt == 3474:
                    print("=" * 40)
                    print(f"DEBUG - Sinflow values for id_pnt = {dp.id_pnt}")
                    print("out_dp.Sinflow content:")
                    for idx, val in enumerate(out_dp.Sinflow):
                        print(f"Sinflow[{idx + 1}] = {val:.6f}")  # Pour correspondre à l'indexation Fortran (1-based)
                    print("=" * 40)

                
        else:
            # If no valid outflow is found, classify dp as a low point / endorheic.
            model.endorheic_count += 1
            dp.id_endo = model.endorheic_count
            # Create an EndoPoint for dp.
            endo = EndoPoint(
                id_eo=model.endorheic_count,
                id_pnt=dp.id_pnt,
                bas_type=ndfl,
                nsaddle=0
            )
            model.endo_pt.append(endo)
        
    model.dr_pt_by_id = {dp.id_pnt: dp for dp in model.dr_pt}
    model.dr_net_by_id = {net.id_ch: net for net in model.dr_net}
    
    
    n_drnet = len(model.dr_net)
    print ('n_drnet, size(dr_net) = ', n_drnet)

    
    finish_time = time.process_time()
    elapsed_time = finish_time - start_time
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = elapsed_time % 60
    print(f"Elapsed time: {hours}h {minutes}m {seconds:.2f}s")



