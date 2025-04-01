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
from data_structures import DrainageNetwork, EndoPoint, IDPointer, ListPointer



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
    dr_net_dict = {}  
    
    for dp in tqdm(model.dr_pt[::-1], desc="Building drainage networks", unit="point"):
    
        if dp.ninf > 0:
            if dp.Linflow:
                i_max = np.argmax(dp.Linflow.value)  
                l_max = dp.Linflow.value[i_max]
                dp.upl = l_max
                dp.sumdev = dp.Sinflow.value[i_max] 
                id_up_max = dp.inflow.value[i_max] 
    
                up_dp = dr_pt_dict.get(id_up_max)
                if up_dp:
                    dp.id_ch = up_dp.id_ch
                    id_main = dp.id_ch
                else:
                    dp.id_ch = IDPointer(None)  # Ensure it's an IDPointer instance
                    id_main = None
    
            for inflow_id in dp.inflow.value:  # Extract values from ListPointer
                inflow_dp = dr_pt_dict.get(inflow_id)
                if inflow_dp is None:
                    continue
            
                curr_ch = inflow_dp.id_ch  # Keep reference to IDPointer
                net = dr_net_dict.get(curr_ch.value)
            
                if net is None:
                    # Create a new network if it doesn't exist
                    net = DrainageNetwork(
                        id_ch=curr_ch.value,  # Channel ID
                        nel=0,  # Number of points in the network
                    )
                    net.id_pnts = ListPointer([inflow_dp.id_pnt])  # List of drainage points
                    net.id_start_pt = IDPointer(inflow_dp.id_pnt)  # Start point
                    net.id_end_pt = IDPointer(inflow_dp.id_pnt)  # End point
                    net.length = 0.0  # Channel length
                    net.id_ch_out = curr_ch  # Outflow channel ID
                    net.n_jun = 0  # Number of tributaries
                    net.id_in = ListPointer([])  # List of tributary IDs
                    net.n_path = 0  # Number of downstream paths
                    net.id_path = ListPointer([])  # List of downstream paths
                    net.id_endo = IDPointer(None)  # Endorheic basin ID
                    net.sso = None  # Strahler Stream Order
                    net.hso = None  # Horton Stream Order
            
                    # Add to the dictionary and list
                    dr_net_dict[curr_ch.value] = net
                    model.dr_net.append(net)
            
                # Expand the list of points dynamically
                net.id_pnts.value.append(dp.id_pnt)  # Add drainage point
                net.nel += 1  # Update the number of elements
                if inflow_id in dp.inflow.value:
                    index = dp.inflow.value.index(inflow_id)  # Trouver la position de inflow_id
                    net.length = dp.Linflow.value[index]  # Prendre la même position dans Linflow

                net.id_ch_out = dp.id_ch
                net.id_end_pt = IDPointer(dp.id_pnt)
            
                # Handle tributary channels
                if id_main is not None and curr_ch.value != id_main.value:
                    main_net = dr_net_dict.get(id_main.value)
                    if main_net:
                        if curr_ch.value not in main_net.id_in.value:
                            main_net.id_in.value.append(curr_ch.value)  # Append tributary ID
            
                            # Increase the number of junctions
                            main_net.n_jun += 1
    
        else:
            # Handling drainage points with no inflows (channel head)
            new_channel_id = len(model.dr_net) + 1
            dp.id_ch = IDPointer(new_channel_id)
            net = DrainageNetwork(
                id_ch=new_channel_id,
                nel=1)
            net.id_pnts=ListPointer([dp.id_pnt])
            net.id_start_pt=IDPointer(dp.id_pnt)
            net.id_end_pt=IDPointer(dp.id_pnt)
            net.length=0.0
            net.id_ch_out=IDPointer(new_channel_id)
            net.n_jun=0
            net.id_in=ListPointer([])
            net.n_path=0
            net.id_path=ListPointer([])
            net.id_endo=IDPointer(None)
            net.sso=None
            net.hso=None
            
            dr_net_dict[new_channel_id] = net
            model.dr_net.append(net)
        
        # --- Compute drainage direction using D8_LTD ---
        i_out, j_out, ndfl, sumdev = model.d8_ltd(dp)
        
        if i_out is not None and j_out is not None:
            # Retrieve the drainage point corresponding to the outflow direction.
            out_dp = model.mat_id[i_out * 2, j_out * 2]
            if out_dp is not None:
                dp.fldir = IDPointer(out_dp.id_pnt)
                dp.fldir_ss = IDPointer(None)
                
                # Update outflow point: increment inflow count and record the current point.
                out_dp.ninf += 1
                out_dp.inflow.value.append(dp.id_pnt)
                
                # Compute the Euclidean distance between dp and out_dp.
                distance = sqrt(((dp.i - out_dp.i) * model.delta_x) ** 2 +
                                ((dp.j - out_dp.j) * model.delta_y) ** 2)
                out_dp.Linflow.value.append(dp.upl + distance)
                out_dp.Sinflow.value.append(sumdev)
        else:
            # If no valid outflow is found, classify dp as a low point / endorheic.
            model.endorheic_count += 1
            dp.id_endo = IDPointer(model.endorheic_count)
            
            # Create an EndoPoint for dp.
            endo = EndoPoint()
            endo.id_eo=IDPointer(model.endorheic_count)
            endo.id_pnt=IDPointer(dp.id_pnt)
            endo.bas_type=ndfl
            endo.nsaddle=0
            
            model.endo_pt.append(endo)
        
    # Update lookup dictionaries
    model.dr_pt_by_id = {dp.id_pnt: dp for dp in model.dr_pt}
    model.dr_net_by_id = {net.id_ch.value: net for net in model.dr_net}
    
    n_drnet = len(model.dr_net)
    print('n_drnet, size(dr_net) = ', n_drnet)
    
    finish_time = time.process_time()
    elapsed_time = finish_time - start_time
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = elapsed_time % 60
    print(f"Elapsed time: {hours}h {minutes}m {seconds:.2f}s")



