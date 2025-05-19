# -*- coding: utf-8 -*-
"""
 THALweg NETwork Horton Stream Order LENGTH method construction

 The subroutine thal_net_hso_length define the thalweg network. If necessary
 the dowslope path that spill out from saddles of endorheic basins is traced.
 At each junction the main channel of the thalweg newtwork is defined as 
 the channel with the maximum upslope length.
"""

from tqdm import tqdm
import numpy as np


from landmark.geomorph_tools.data_structures import DrainageNetwork, DrainagePointInflow


def thal_net_hso_length(model, river_mask):
    """
    Builds the drainage network and calculates Horton Stream Order (HSO)
    based on a river mask and upslope/downstream relationships.

    The function processes drainage points in decreasing elevation order,
    assigns them to drainage networks, computes channel length,
    and updates subsurface flow direction (fldir_ss).

    Parameters
    ----------
    model : object
        Hydrological model containing:
        - dr_pt : list of DrainagePoint
        - dr_net : list to store DrainageNetwork objects
        - mat_id : 2D np.ndarray for visualization
        - delta_x : float, grid spacing
        - etc.

    river_mask : np.ndarray (bool)
        Boolean mask indicating which cells are already considered rivers.
        Used to reset `mat_id` in those cells before reconstruction.
    """

    # ▸ Step 1: Sort drainage points by elevation (lowest to highest)
    qoi = [dp.id_pnt for dp in sorted(model.dr_pt, key=lambda dp: dp.Z)]

    # ▸ Step 2: Reinitialize pointers and structures
    print("Reset id_ch, inflow, Linflow")
    for dp in tqdm(model.dr_pt):
        dp.reset_flow_data()

    # ▸ Step 3: Prepare clean containers
    del(model.dr_net)
    dr_net = [DrainageNetwork() for _ in tqdm(range(len(model.dr_pt)), desc="Drainage network reset")]
    dr_pt_in = [DrainagePointInflow() for _ in tqdm(range(len(model.dr_pt)), desc="Temporary drainage points")]
    model.mat_id = np.where(river_mask, None, model.mat_id)

    inet = 0  # Counter for active drainage networks

    # ▸ Step 4: Main loop — from highest Z to lowest
    for id_dr in tqdm(qoi[::-1], desc="Building drainage networks HSO length", unit="point"):
        dp = model.dr_pt[id_dr - 1]
        i_curr, j_curr = dp.i, dp.j

        # Only consider points with no inflow and valid basin
        if dp.ninf == 0 and dp.id_endo >= 0:
            if dp.fldir != None or dp.fldir_ss != None:
                if dp.fldir_ss != None:
                    curr_fldir = dp.fldir_ss
                    dp.fldir = dp.fldir_ss
                    
                else:
                    curr_fldir = dp.fldir

                max_Z = dp.Z

                if curr_fldir is not None:
                    dp_fldir = model.dr_pt[curr_fldir - 1]

                    if dp.A_in > 0:
                        # ▸ Case: point already belongs to an existing drainage network
                        curr_ch = dp.id_ch
                        dr_net[curr_ch-1].nel += 1
                        dr_net[curr_ch-1].sso = dr_net[curr_ch-1].hso
                        dr_net[curr_ch-1].id_pnts.append(curr_fldir)
                        dr_net[curr_ch-1].id_end_pt = curr_fldir #id of the end point
                        dr_net[curr_ch-1].length += model.delta_x*((i_curr-dp_fldir.i)**2 + (j_curr-dp_fldir.j)**2)**0.5
                        i_mat = i_curr*2 + (dp_fldir.i - i_curr)
                        j_mat = j_curr*2 + (dp_fldir.j - j_curr)
                        model.mat_id[i_mat, j_mat] = dr_net[curr_ch-1].id_ch

                        # ▸ If downstream point already processed
                        if dp_fldir.A_in > 0:
                            dr_pt_in[curr_fldir - 1].ninf += 1
                            dr_pt_in[curr_fldir - 1].inflow.append(dp.id_ch)

                            # ▸ Compare stream orders for hierarchical assignment
                            if dr_net[curr_ch-1].sso == dr_net[dp_fldir.id_ch-1].sso:
                                # upgrading order of the channel with greater length
                                if dr_net[curr_ch-1].length > dp_fldir.upl:
                                    #feature point update as the main channel is the current
                                    dr_net[dp_fldir.id_ch-1].hso = dr_net[dp_fldir.id_ch-1].sso
                                    dr_net[curr_ch-1].hso += 1
                                    dp_fldir.upl = dr_net[curr_ch-1].length
                                    dp_fldir.id_ch = dr_net[curr_ch-1].id_ch
                                    dp_fldir.A_in += dp.A_in+1
                                    dp_fldir.ninf -= 1
                                    dr_net[curr_ch-1].id_ch_out = dp.id_ch
                                    if dp_fldir.Z >= max_Z:
                                        dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
                                else:
                                    dr_net[curr_ch-1].hso = dr_net[curr_ch-1].sso
                                    dr_net[dp_fldir.id_ch-1].hso += 1
                                    dp_fldir.A_in += dp.A_in + 1
                                    dp_fldir.ninf -= 1
                                    if dp_fldir.Z >= max_Z:
                                        dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
                            
                            if dr_net[curr_ch-1].sso > dr_net[dp_fldir.id_ch-1].sso:                            
                                # upgrading order of the channel with greater length                  
                                # feature point update as the main channel is the current
                                dr_net[dp_fldir.id_ch-1].hso = dr_net[dp_fldir.id_ch-1].sso
                                dp_fldir.upl = dr_net[curr_ch-1].length
                                dp_fldir.id_ch = dr_net[curr_ch-1].id_ch
                                dp_fldir.A_in += dp.A_in+1
                                dp_fldir.ninf -= 1
                                dr_net[curr_ch-1].id_ch_out = dp.id_ch
                                if dp_fldir.Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
                                    
                            if dr_net[curr_ch-1].sso < dr_net[dp_fldir.id_ch-1].sso:                            
                                dp_fldir.A_in += dp.A_in+1
                                dp_fldir.ninf -= 1
                                if dp_fldir.Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
                        
                        else:
                            dp_fldir.upl = dr_net[curr_ch-1].length
                            dp_fldir.id_ch = dr_net[curr_ch-1].id_ch
                            dp_fldir.A_in += dp.A_in+1
                            dp_fldir.ninf -= 1
                            dr_net[curr_ch-1].id_ch_out = dp.id_ch
                            dr_pt_in[curr_fldir-1].ninf = 1
                            dr_pt_in[curr_fldir-1].inflow.append(dp.id_ch)
                            if dp_fldir.Z >= max_Z:
                                dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
                    

                    else:
                        # ▸ New head channel
                        inet += 1
                        dr_net_inet = dr_net[inet-1]
                        if dr_net_inet.id_ch == 0 or dr_net_inet.id_ch == None:
                            dr_net_inet.id_ch = inet
                            dp.id_ch = dr_net_inet.id_ch
                            dr_net_inet.nel = 2
                            dr_net_inet.id_pnts.append(dp.id_pnt)
                            dr_net_inet.id_start_pt = dp.id_pnt #id of the start point
                            dr_net_inet.id_ch_out = dr_net_inet.id_ch #first assignement
                            dr_net_inet.id_pnts.append(curr_fldir)
                            dr_net_inet.id_end_pt = curr_fldir
                            dr_net_inet.length = model.delta_x*((i_curr-dp_fldir.i)**2 + (j_curr-dp_fldir.j)**2)**0.5
                            dr_net_inet.sso = 1
                            dr_net_inet.hso = 1
                            i_mat = i_curr*2 + (dp_fldir.i - i_curr)
                            j_mat = j_curr*2 + (dp_fldir.j - j_curr)
                            model.mat_id[i_mat, j_mat] = dr_net[inet-1].id_ch
                            dr_pt_in[id_dr-1].ninf = 1
                            dr_pt_in[id_dr-1].inflow.append(dr_net_inet.id_ch)
                            if dp_fldir.A_in > 0:
                                dr_pt_in[curr_fldir-1].ninf += 1
                                dr_pt_in[curr_fldir-1].inflow.append(dr_net_inet.id_ch)
                                if dr_net_inet.sso == dr_net[dp_fldir.id_ch-1].sso:
                                    #upgrading order of the channel with greater length
                                    if dr_net_inet.length > dp_fldir.upl:
                                        #feature point update as the main channel is the current
                                        dr_net[dp_fldir.id_ch-1].hso = dr_net[dp_fldir.id_ch-1].sso
                                        dr_net_inet.hso += 1
                                        dp_fldir.upl = dr_net_inet.length
                                        dp_fldir.id_ch = dr_net_inet.id_ch
                                        dp_fldir.A_in += dp.A_in+1
                                        dp_fldir.ninf -= 1
                                        if dp_fldir.Z >= max_Z:
                                            dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
                                    else:
                                        dr_net_inet.hso = dr_net_inet.sso
                                        if dr_net[dp_fldir.id_ch-1].hso == dr_net[dp_fldir.id_ch-1].sso:
                                            dr_net[dp_fldir.id_ch-1].hso += 1
                                        dp_fldir.A_in += dp.A_in+1
                                        dp_fldir.ninf -= 1
                                        if dp_fldir.Z >= max_Z:
                                            dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
                                if dr_net_inet.sso > dr_net[dp_fldir.id_ch-1].sso:
                                    # upgrading order of the channel with greater length                  
                                    # feature point update as the main channel is the current
                                    dp_fldir.upl = dr_net_inet.length
                                    dp_fldir.id_ch = dr_net_inet.id_ch
                                    dp_fldir.A_in += dp.A_in+1
                                    dp_fldir.ninf -= 1
                                    if dp_fldir.Z >= max_Z:
                                        dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
                                        
                                if dr_net_inet.sso < dr_net[dp_fldir.id_ch-1].sso:
                                    dp_fldir.A_in += dp.A_in+1
                                    dp_fldir.ninf -= 1
                                    if dp_fldir.Z >= max_Z:
                                        dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
                            else:
                                #downslope point pointed by curr_fldir never processed
                                dp_fldir.upl = dr_net_inet.length
                                dp_fldir.id_ch = dr_net_inet.id_ch
                                dp_fldir.A_in += dp.A_in+1
                                dp_fldir.ninf -= 1
                                dr_pt_in[curr_fldir-1].ninf = 1
                                dr_pt_in[curr_fldir-1].inflow.append(dr_net_inet.id_ch)
                                if dp_fldir.Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)

                dp.id_endo = -1

    # ▸ Final assignment
    model.dr_net = dr_net[:inet]
    model.dr_pt_in = dr_pt_in

    # ▸ Update downstream channel references
    for curr_net in tqdm(model.dr_net, desc=" Updating downstream channel references", unit="channel"):
        curr_net.id_ch_out = model.dr_pt[curr_net.id_end_pt - 1].id_ch


                                
def dwnslp_hso(model, id_dr, dr_net, dr_pt_in, max_Z):
    """Recursive subroutine to trace downslope flow paths from points that spill from endorheic basins,
    while maintaining Horton Stream Order (HSO) consistency.

    This function is called when a spill occurs across a saddle. It updates channel hierarchy,
    flow lengths, inflow links, and recursively propagates downstream if needed.

    Parameters
    ----------
    model : object
        Hydrological model with DrainagePoints, DrainageNetworks, and grid info.
    id_dr : int
        DrainagePoint ID (1-based) from which to continue the spill trace.
    dr_net : list of DrainageNetwork
        Current drainage networks being updated.
    dr_pt_in : list of DrainagePointInflow
        Temporary inflow record for each drainage point.
    max_Z : float
        Maximum elevation threshold above which recursive tracing is allowed.

    Returns
    -------
    tuple
        (dr_net, dr_pt_in) updated after processing the spill path.
    """

    dp = model.dr_pt[id_dr - 1]
    i_curr, j_curr = dp.i, dp.j

    if dp.ninf == 0 and dp.id_endo >= 0:
        # Select the correct downstream direction
            if dp.fldir_ss != None:
                curr_fldir = dp.fldir_ss
                dp.fldir = dp.fldir_ss
            else:
                curr_fldir = dp.fldir

            if curr_fldir is not None:
                dp_fldir = model.dr_pt[curr_fldir - 1]

                # ▸ If upslope area is already known: part of existing network
                if dp.A_in > 0:
                    curr_ch = dp.id_ch
                    dr_net[curr_ch-1].nel += 1
                    dr_net[curr_ch-1].sso = dr_net[curr_ch-1].hso
                    dr_net[curr_ch-1].id_pnts.append(curr_fldir)
                    dr_net[curr_ch-1].id_end_pt = curr_fldir #id of the end point
                    dr_net[curr_ch-1].length += model.delta_x*((i_curr-dp_fldir.i)**2 + (j_curr-dp_fldir.j)**2)**0.5
                    i_mat = i_curr*2 + (dp_fldir.i - i_curr)
                    j_mat = j_curr*2 + (dp_fldir.j - j_curr)
                    model.mat_id[i_mat, j_mat] = dr_net[curr_ch-1].id_ch
                    
                    if dp_fldir.A_in > 0:
                        # This means that it has already been processed in dr_net_ss and the upl value has been updated
                        dr_pt_in[curr_fldir-1].ninf += 1
                        dr_pt_in[curr_fldir-1].inflow.append(dr_net[curr_ch-1].id_ch) 
                        if dr_net[curr_ch-1].sso == dr_net[dp_fldir.id_ch-1].sso:
                            # upgrading order of the channel with greater length
                            if dr_net[curr_ch-1].length > dp_fldir.upl:
                                #feature point update as the main channel is the current
                                dr_net[dp_fldir.id_ch-1].hso = dr_net[dp_fldir.id_ch-1].sso
                                dr_net[curr_ch-1].hso += 1
                                dp_fldir.upl = dr_net[curr_ch-1].length
                                dp_fldir.id_ch = dr_net[curr_ch-1].id_ch
                                dp_fldir.A_in += dp.A_in+1
                                dp_fldir.ninf -= 1
                                dr_net[curr_ch-1].id_ch_out = dp_fldir.id_ch
                                if dp_fldir.Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
                            else:
                                dr_net[curr_ch-1].hso = dr_net[curr_ch-1].sso
                                dr_net[dp_fldir.id_ch-1].hso += 1
                                dp_fldir.A_in += dp.A_in + 1
                                dp_fldir.ninf -= 1
                                if dp_fldir.Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
                        
                        if dr_net[curr_ch-1].sso > dr_net[dp_fldir.id_ch-1].sso:                            
                            # upgrading order of the channel with greater length                  
                            # feature point update as the main channel is the current
                            dr_net[dp_fldir.id_ch-1].hso = dr_net[dp_fldir.id_ch-1].sso
                            dp_fldir.upl = dr_net[curr_ch-1].length
                            dp_fldir.id_ch = dr_net[curr_ch-1].id_ch
                            dp_fldir.A_in += dp.A_in+1
                            dp_fldir.ninf -= 1
                            dr_net[curr_ch-1].id_ch_out = dp.id_ch
                            if dp_fldir.Z >= max_Z:
                                dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
                                
                        if dr_net[curr_ch-1].sso < dr_net[dp_fldir.id_ch-1].sso:                            
                            dp_fldir.A_in += dp.A_in+1
                            dp_fldir.ninf -= 1
                            if dp_fldir.Z >= max_Z:
                                dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
                    
                    else:
                        dp_fldir.upl = dr_net[curr_ch-1].length
                        dp_fldir.id_ch = dr_net[curr_ch-1].id_ch
                        dp_fldir.A_in += dp.A_in+1
                        dp_fldir.ninf -= 1
                        dr_net[curr_ch-1].id_ch_out = dp_fldir.id_ch
                        dr_pt_in[curr_fldir-1].ninf = 1
                        dr_pt_in[curr_fldir-1].inflow.append(dr_net[curr_ch-1].id_ch)
                        if dp_fldir.Z >= max_Z:
                            dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir, dr_net, dr_pt_in, max_Z)
            dp.id_endo = -1

    return dr_net, dr_pt_in
                                
