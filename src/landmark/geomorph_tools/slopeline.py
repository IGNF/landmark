# -*- coding: utf-8 -*-

"""
Module for constructing drainage networks (slopelines) based on the D8-LTD algorithm.
This module uses the previously loaded DEM data (via the model) and calls the D8_LTD algorithm.
"""

from math import sqrt
import numpy as np
from tqdm import tqdm

# Import the data structures.
from landmark.geomorph_tools.data_structures import DrainageNetwork, EndoPoint


def calculate_slopelines(model):
    """Constructs drainage networks (slopelines) by processing each drainage point.
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
           - dr_pt: list of DrainagePoint objects 
           - mat_id: 2D numpy array with DrainagePoint objects at positions (i*2, j*2)
           - N, M: dimensions of the DEM (number of rows, number of columns)
           - delta_x, delta_y: grid spacings
           - nodata: no-data value for the DEM
    """
    
    # Initialize drainage network and endorheic lists.
    model.dr_net = [None] * len(model.dr_pt)      # List of DrainageNetwork objects.
    model.l_endo_pt = []     # List of EndoPoint objects.
    model.endorheic_count = 0  # Counter for endorheic points.
    inet = 0
    
    qoi = [dp.id_pnt for dp in sorted(model.dr_pt, key=lambda dp: dp.Z)]

    for id_dr in tqdm(qoi[::-1], desc="Building drainage networks", unit="point"):
        dp = model.dr_pt[id_dr-1]
    
        if dp.ninf > 0:
            # channel points
            # all upslope points are already considered
            i_max = np.argmax(dp.Linflow)  
            l_max = dp.Linflow[i_max]
            dp.upl = l_max
            dp.sumdev = dp.Sinflow[i_max] 
            id_up_max = dp.inflow[i_max] #Id of the point with maximum upslope length
            dp.id_ch = model.dr_pt[id_up_max-1].id_ch 
            id_main = dp.id_ch

            
            for length_inflow, id_inflow in zip(dp.Linflow,dp.inflow):  
                curr_ch = model.dr_pt[id_inflow-1].id_ch
                model.dr_net[curr_ch-1].id_pnts.append(dp.id_pnt)
                model.dr_net[curr_ch-1].nel+=1
                model.dr_net[curr_ch-1].length = length_inflow
                model.dr_net[curr_ch-1].id_ch_out = dp.id_ch
                model.dr_net[curr_ch-1].id_end_pt = dp.id_pnt
                

            
                # Handle tributary channels
                if curr_ch != id_main:
                    #Add triburaties index to main channel
                    if curr_ch not in model.dr_net[id_main-1].id_in:
                        model.dr_net[id_main-1].n_jun += 1
                        model.dr_net[id_main-1].id_in.append(curr_ch)
                        
    
        else:
            # Handling drainage points with no inflows (channel head)
            inet += 1
            new_channel_id = inet
            dp.id_ch = new_channel_id
            net = DrainageNetwork(
                id_ch=new_channel_id,
                nel=1)
            net.id_pnts=[dp.id_pnt]
            net.id_start_pt=dp.id_pnt
            net.id_end_pt=dp.id_pnt
            net.id_ch_out=new_channel_id
            
            model.dr_net[inet-1] = net
        
        i_out, j_out, ndfl, sumdev = model.d8_ltd(dp)
        
        if i_out is not None and j_out is not None:
            # Retrieve the drainage point corresponding to the outflow direction.
            out_dp = model.mat_id[i_out * 2, j_out * 2]
            if out_dp is not None:
                dp.fldir = out_dp.id_pnt
                dp.fldir_ss = None
                
                # Update outflow point: increment inflow count and record the current point.
                out_dp.ninf += 1
                out_dp.inflow.append(dp.id_pnt)
                
                # Compute the Euclidean distance between dp and out_dp.
                distance = sqrt(((dp.i - out_dp.i) * model.delta_x) ** 2 +
                                ((dp.j - out_dp.j) * model.delta_y) ** 2)
                out_dp.Linflow.append(dp.upl + distance)
                out_dp.Sinflow.append(sumdev)
        else:
            # If no valid outflow is found, classify dp as a low point / endorheic.
            model.endorheic_count += 1
            dp.id_endo = model.endorheic_count
            
            # Create an EndoPoint for dp.
            endo = EndoPoint()
            endo.id_eo=model.endorheic_count
            endo.id_pnt=dp.id_pnt
            endo.bas_type=ndfl
            endo.nsaddle=0
            
            model.l_endo_pt.append(endo)
            
            
    model.dr_net = model.dr_net[:inet]
        
