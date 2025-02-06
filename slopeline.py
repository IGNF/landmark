# -*- coding: utf-8 -*-

"""
Module for constructing drainage networks (slopelines) based on the D8-LTD algorithm.
This module uses the previously loaded DEM data (via the loader) and calls the D8_LTD algorithm.
"""

import time
from math import sqrt
import numpy as np
from tqdm import tqdm

# Import the data structures.
from data_structures import DrainageNetwork, EndoPoint



def calculate_slopelines(loader):
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
    loader : object
        An instance of your data loader class which must have at least:
           - dr_pt: list of DrainagePoint objects (sorted in ascending order of elevation)
           - mat_id: 2D numpy array with DrainagePoint objects at positions (i*2, j*2)
           - N, M: dimensions of the DEM (number of rows, number of columns)
           - delta_x, delta_y: grid spacings
           - nodata: no-data value for the DEM
    """
    start_time = time.process_time()
    
    # Initialize drainage network and endorheic lists.
    loader.dr_net = []      # List of DrainageNetwork objects.
    loader.endo_pt = []     # List of EndoPoint objects.
    loader.endorheic_count = 0  # Counter for endorheic points.
    
    dr_pt_dict = {p.id_pnt: p for p in loader.dr_pt}

    
    # Process drainage points in descending order (assuming loader.dr_pt is sorted in ascending order)
    for dp in tqdm(loader.dr_pt[::-1], desc="Building drainage networks", unit="point"):
        
        # --- Drainage network update ---
        if dp.ninf > 0:
            # The point has upstream inflows: treat it as a channel point.
            # (Assuming that dp.Linflow and dp.Sinflow have been filled previously during loading.)
            if dp.Linflow:
                l_max = max(dp.Linflow)
                i_max = dp.Linflow.index(l_max)
                dp.upl = l_max
                # Update cumulative deviation from upstream.
                dp.sumdev = dp.Sinflow[i_max] if dp.Sinflow else 0.0
                # Retrieve the id of the inflow that has the maximum upstream length.
                id_up_max = dp.inflow[i_max] if dp.inflow else None
                if id_up_max is not None:
                    # Find the upstream drainage point.
                    up_dp = dr_pt_dict.get(id_up_max)

                    if up_dp is not None:
                        dp.id_ch = up_dp.id_ch
                        id_main = dp.id_ch
                    else:
                        dp.id_ch = None
                        id_main = None
                else:
                    id_main = None
            
            # Update the drainage network for each inflow.
            for inflow_id in dp.inflow:
                inflow_dp = next((p for p in loader.dr_pt if p.id_pnt == inflow_id), None)
                if inflow_dp is None:
                    continue
                curr_ch = inflow_dp.id_ch
                # Check if a network with channel id curr_ch already exists.
                existing_nets = [net for net in loader.dr_net if net.id_ch == curr_ch]
                if not existing_nets:
                    # Create a new network for this channel.
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
                    loader.dr_net.append(net)
                else:
                    net = existing_nets[0]
                # Append current point id to the network.
                net.id_pnts.append(dp.id_pnt)
                # For simplicity, assign the network's length as dp.upl.
                if dp.Linflow:
                    net.length = dp.upl
                net.id_ch_out = dp.id_ch
                net.id_end_pt = dp.id_pnt
                # If this channel is not the main one, add it as a tributary to the main channel.
                if id_main is not None and curr_ch != id_main:
                    main_net = next((net for net in loader.dr_net if net.id_ch == id_main), None)
                    if main_net is not None and curr_ch not in main_net.id_in:
                        main_net.id_in.append(curr_ch)
        else:
            # The point has no upstream inflows: it is a channel head.
            new_channel_id = len(loader.dr_net) + 1
            dp.id_ch = new_channel_id
            net = DrainageNetwork(
                id_ch=new_channel_id,
                id_pnts=[dp.id_pnt],
                id_start_pt=dp.id_pnt,
                id_end_pt=dp.id_pnt,
                length=0.0,
                id_ch_out=new_channel_id,
                n_jun=0,
                id_in=[]
            )
            loader.dr_net.append(net)
        
        # --- Compute drainage direction using D8_LTD ---
        # Call the D8_LTD function on the current drainage point.
        i_out, j_out, ndfl, sumdev = loader.d8_ltd(dp)
        
        if i_out is not None and j_out is not None and i_out > 0 and j_out > 0:
            # Retrieve the drainage point corresponding to the outflow direction.
            out_dp = loader.mat_id[i_out * 2, j_out * 2]
            if out_dp is not None:
                dp.fdir = out_dp.id_pnt
                # Update outflow point: increment inflow count and record the current point.
                out_dp.ninf += 1
                out_dp.inflow.append(dp.id_pnt)
                # Compute the Euclidean distance between dp and out_dp.
                distance = sqrt(((dp.i - out_dp.i) * loader.delta_x)**2 +
                                ((dp.j - out_dp.j) * loader.delta_y)**2)
                out_dp.Linflow.append(dp.upl + distance)
                out_dp.Sinflow.append(sumdev)
        else:
            # If no valid outflow is found, classify dp as a low point / endorheic.
            loader.endorheic_count += 1
            dp.id_endo = loader.endorheic_count
            # Create an EndoPoint for dp.
            endo = EndoPoint(
                id_eo=loader.endorheic_count,
                id_pnt=dp.id_pnt,
                bas_type=ndfl,
                nsaddle=0
            )
            loader.endo_pt.append(endo)
        
    
    finish_time = time.process_time()
    elapsed_time = finish_time - start_time
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = elapsed_time % 60
    print(f"Elapsed time: {hours}h {minutes}m {seconds:.2f}s")




# def facet(e0, e1, delta_x, delta_y):
#     """
#     Compute slope and angle of steepest descent between three elevation points.

#     Parameters
#     ----------
#     e0 : float
#         Elevation of the central cell.
#     e1 : float
#         Elevation of the neighboring cell.
#     delta_x : float
#         Grid spacing in the x-direction.
#     delta_y : float
#         Grid spacing in the y-direction.

#     Returns
#     -------
#     r_val : float
#         Angle of steepest descent.
#     """
#     dzdx = (e1 - e0) / delta_x
#     dzdy = (e1 - e0) / delta_y
#     r_val = np.arctan2(dzdy, dzdx)
#     return r_val

# def D8_flow_direction(dem, delta_x, delta_y):
#     """
#     Calculate the D8 flow direction and deviation values for a given DEM.

#     Parameters
#     ----------
#     dem : numpy.ndarray
#         2D array representing the digital elevation model (DEM).
#     delta_x : float
#         Grid spacing in the x-direction.
#     delta_y : float
#         Grid spacing in the y-direction.

#     Returns
#     -------
#     flow_dir : numpy.ndarray
#         2D array where each cell contains the flow direction (D8 coding).
#     ndfl : numpy.ndarray
#         2D array indicating valid flow cells (1) or boundaries (0).
#     dev : numpy.ndarray
#         2D array with deviation values computed between connected points.
#     """
#     d8_offsets = [
#         (0, 0),  # Central (0, "no flow")
#         (-1, 0),  # North (64)
#         (-1, 1),  # North-East (128)
#         (0, 1),  # East (1)
#         (1, 1),  # South-East (2)
#         (1, 0),  # South (4)
#         (1, -1),  # South-West (8)
#         (0, -1),  # West (16)
#         (-1, -1),  # North-West (32)
#     ]
#     d8_codes = [0, 64, 128, 1, 2, 4, 8, 16, 32]
#     d8_codes_inv = [0, 4, 8, 16, 32, 64, 128, 1, 2]  # Directions reversed to match movement

#     dem_padded = np.pad(dem, pad_width=1, mode="constant", constant_values=np.nan)
#     flow_dir = np.zeros_like(dem, dtype=np.uint8)
#     dev = np.full_like(dem, np.nan, dtype=float)
#     min_altitude = np.full_like(dem, np.inf, dtype=dem.dtype)

#     for (dy, dx), code in tqdm(zip(d8_offsets, d8_codes_inv)):
#         shifted = np.roll(dem_padded, (dy, dx), axis=(0, 1))[1:-1, 1:-1]
#         mask = shifted < min_altitude
#         flow_dir[mask] = code
#         min_altitude[mask] = shifted[mask]

#     ndfl = np.ones_like(dem, dtype=np.uint8)
#     ndfl[0, :] = 0
#     ndfl[-1, :] = 0
#     ndfl[:, 0] = 0
#     ndfl[:, -1] = 0

#     for i in range(1, dem.shape[0] - 1):
#         for j in range(1, dem.shape[1] - 1):
#             if flow_dir[i, j] != 0:
#                 dy, dx = d8_offsets[d8_codes.index(flow_dir[i, j])]
#                 i_out, j_out = i + dy, j + dx
#                 e0, e1 = dem[i, j], dem[i_out, j_out]
#                 r_val = facet(e0, e1, delta_x, delta_y)
#                 if abs(dx) + abs(dy) == 1:
#                     dev[i, j] = delta_x * np.sin(r_val)
#                 else:
#                     dev[i, j] = delta_x * np.sqrt(2) * np.sin(np.pi / 4 - r_val)

#     return flow_dir, ndfl, dev



# def calculate_slopelines(loader):
#     """
#     Constructs drainage networks (slopelines) by processing each drainage point.
#     For each point (processed in descending order of elevation), the function:
#       - If the point has upstream inflows (ninf > 0), updates the corresponding drainage network.
#       - If the point has no inflows, it is treated as a channel head and a new network is created.
#       - Calls d8_ltd() for the point to compute its outflow direction and cumulative deviation.
#       - Updates the outflow point's inflow information or, if no valid outflow is found,
#         classifies the point as belonging to an endorheic basin.
      
#     Parameters
#     ----------
#     loader : object
#         An instance of your data loader class which must have at least:
#            - dr_pt: list of DrainagePoint objects (sorted in ascending order of elevation)
#            - mat_id: 2D numpy array with DrainagePoint objects at positions (i*2, j*2)
#            - N, M: dimensions of the DEM (number of rows, number of columns)
#            - delta_x, delta_y: grid spacings
#            - nodata: no-data value for the DEM
#     """
#     start_time = time.process_time()
    
#     # Initialize drainage network and endorheic lists.
#     loader.dr_net = []      # List of DrainageNetwork objects.
#     loader.endo_pt = []     # List of EndoPoint objects.
#     loader.endorheic_count = 0  # Counter for endorheic points.
    
#     #Initialize arrays
#     dem  = loader.dem
#     sumdev_array = np.zeros_like(dem)
#     flow_dir_array, ndfl_array, dev_array = D8_flow_direction(dem, loader.delta_x, loader.delta_y ) 
    
    
#     d8_offsets = [
#         (0, 0),  # Central (0, "no flow")
#         (-1, 0),  # North (64)
#         (-1, 1),  # North-East (128)
#         (0, 1),  # East (1)
#         (1, 1),  # South-East (2)
#         (1, 0),  # South (4)
#         (1, -1),  # South-West (8)
#         (0, -1),  # West (16)
#         (-1, -1),  # North-West (32)
#     ]
#     d8_codes = [0, 64, 128, 1, 2, 4, 8, 16, 32]

    
    
    
#     # Process drainage points in descending order (assuming loader.dr_pt is sorted in ascending order)
#     for dp in tqdm(reversed(loader.dr_pt), desc="Building drainage networks", unit="point"):
        
#         # --- Drainage network update ---
#         if dp.ninf > 0:
#             # The point has upstream inflows: treat it as a channel point.
#             # (Assuming that dp.Linflow and dp.Sinflow have been filled previously during loading.)
#             if dp.Linflow:
#                 l_max = max(dp.Linflow)
#                 i_max = dp.Linflow.index(l_max)
#                 dp.upl = l_max
#                 # Update cumulative deviation from upstream.
#                 dp.sumdev = dp.Sinflow[i_max] if dp.Sinflow else 0.0
#                 # Retrieve the id of the inflow that has the maximum upstream length.
#                 id_up_max = dp.inflow[i_max] if dp.inflow else None
#                 if id_up_max is not None:
#                     # Find the upstream drainage point.
#                     up_dp = next((p for p in loader.dr_pt if p.id_pnt == id_up_max), None)
#                     if up_dp is not None:
#                         dp.id_ch = up_dp.id_ch
#                         id_main = dp.id_ch
#                     else:
#                         dp.id_ch = None
#                         id_main = None
#                 else:
#                     id_main = None
            
#             # Update the drainage network for each inflow.
#             for inflow_id in dp.inflow:
#                 inflow_dp = next((p for p in loader.dr_pt if p.id_pnt == inflow_id), None)
#                 if inflow_dp is None:
#                     continue
#                 curr_ch = inflow_dp.id_ch
#                 # Check if a network with channel id curr_ch already exists.
#                 existing_nets = [net for net in loader.dr_net if net.id_ch == curr_ch]
#                 if not existing_nets:
#                     # Create a new network for this channel.
#                     net = DrainageNetwork(
#                         id_ch=curr_ch,
#                         id_pnts=[inflow_dp.id_pnt],
#                         id_start_pt=inflow_dp.id_pnt,
#                         id_end_pt=inflow_dp.id_pnt,
#                         length=0.0,
#                         id_ch_out=curr_ch,
#                         n_jun=0,
#                         id_in=[]
#                     )
#                     loader.dr_net.append(net)
#                 else:
#                     net = existing_nets[0]
#                 # Append current point id to the network.
#                 net.id_pnts.append(dp.id_pnt)
#                 # For simplicity, assign the network's length as dp.upl.
#                 if dp.Linflow:
#                     net.length = dp.upl
#                 net.id_ch_out = dp.id_ch
#                 net.id_end_pt = dp.id_pnt
#                 # If this channel is not the main one, add it as a tributary to the main channel.
#                 if id_main is not None and curr_ch != id_main:
#                     main_net = next((net for net in loader.dr_net if net.id_ch == id_main), None)
#                     if main_net is not None and curr_ch not in main_net.id_in:
#                         main_net.id_in.append(curr_ch)
#         else:
#             # The point has no upstream inflows: it is a channel head.
#             new_channel_id = len(loader.dr_net) + 1
#             dp.id_ch = new_channel_id
#             net = DrainageNetwork(
#                 id_ch=new_channel_id,
#                 id_pnts=[dp.id_pnt],
#                 id_start_pt=dp.id_pnt,
#                 id_end_pt=dp.id_pnt,
#                 length=0.0,
#                 id_ch_out=new_channel_id,
#                 n_jun=0,
#                 id_in=[]
#             )
#             loader.dr_net.append(net)
        
#         # --- Compute drainage direction using D8_LTD ---
#         dy, dx = d8_offsets[d8_codes.index(flow_dir_array[dp.i, dp.j])]
#         i_out, j_out = dp.i + dy, dp.j + dx
#         ndfl = ndfl_array[i_out, j_out]
#         dp.sumdev = sumdev_array[i_out, j_out] + dev_array[i_out, j_out]
#         sumdev_array[i_out, j_out] = dp.sumdev
        
        
#         if i_out is not None and j_out is not None and i_out > 0 and j_out > 0:
#             # Retrieve the drainage point corresponding to the outflow direction.
#             out_dp = loader.mat_id[i_out * 2, j_out * 2]
#             if out_dp is not None:
#                 dp.fdir = out_dp.id_pnt
#                 # Update outflow point: increment inflow count and record the current point.
#                 out_dp.ninf += 1
#                 out_dp.inflow.append(dp.id_pnt)
#                 # Compute the Euclidean distance between dp and out_dp.
#                 distance = sqrt(((dp.i - out_dp.i) * loader.delta_x)**2 +
#                                 ((dp.j - out_dp.j) * loader.delta_y)**2)
#                 out_dp.Linflow.append(dp.upl + distance)
#                 out_dp.Sinflow.append(dp.sumdev)
#         else:
#             # If no valid outflow is found, classify dp as a low point / endorheic.
#             loader.endorheic_count += 1
#             dp.id_endo = loader.endorheic_count
#             # Create an EndoPoint for dp.
#             endo = EndoPoint(
#                 id_eo=loader.endorheic_count,
#                 id_pnt=dp.id_pnt,
#                 bas_type=ndfl,
#                 nsaddle=0
#             )
#             loader.endo_pt.append(endo)
        
    
#     finish_time = time.process_time()
#     elapsed_time = finish_time - start_time
#     hours = int(elapsed_time // 3600)
#     minutes = int((elapsed_time % 3600) // 60)
#     seconds = elapsed_time % 60
#     print(f"Elapsed time: {hours}h {minutes}m {seconds:.2f}s")
