# -*- coding: utf-8 -*-
"""
Transcript in python of the original dpl.f90

! The subroutine Downslope Path Legth calculates the legth of the path 
! between each DTM cell and the outflow point even if the basin is 
! endorheic 
"""

import numpy as np
from tqdm import tqdm
import time


def dpl(model):
    """
    Calculates the downslope path length (dpl) for each drainage point and updates the
    drainage network (dr_net) information. This routine works over the endorheic points
    (model.endo_pt) and then, for each channel, computes the downslope path length for each
    drainage point, and then recursively updates tributary channels.
    
    Parameters
    ----------
    model : object
        The hydrological model object that contains at least:
           - dr_pt: list of DrainagePoint objects.
           - dr_pt_by_id: dict mapping id_pnt -> DrainagePoint.
           - dr_net: list of DrainageNetwork objects.
           - dr_net_by_id: dict mapping channel id (id_ch) -> DrainageNetwork.
           - mat_id: 2D numpy array (of size (N*2-1, M*2-1)) containing DrainagePoint objects.
           - endo_pt: list of EndoPoint objects.
           - N, M, delta_x, delta_y, nodata: global parameters.
    """
    start_time = time.process_time()
    
    # Process each endorheic point.
    for endo in tqdm(model.endo_pt, desc="Processing endorheic channels", unit="channel"):
        # Determine the main channel for this endorheic basin:
        # main is obtained from the DrainagePoint corresponding to the endo_pt.
        # (En Fortran : main => dr_pt(endo_pt(cnt_endo)%id_pnt)%id_ch)
        main_channel = model.dr_pt_by_id[endo.id_pnt.value].id_ch
        # if main_channel is None:
        #     continue
        # Assign the channel’s endorheic id:
        model.dr_net_by_id[main_channel.value].id_endo = endo.id_eo

        # For each drainage point in the main channel:
        main_net = model.dr_net_by_id[main_channel.value]
        for dp_id in main_net.id_pnts.value:
            dp = model.dr_pt_by_id[dp_id]
            # Calculate downslope path length:
            # l1 is the channel length; l2 is the upstream length (upl) of the point.
            l1 = main_net.length
            l2 = dp.upl
            dp.dpl = l1 - l2
            # Also update the matrix (mat_id) if possible.
            i_curr, j_curr = dp.i, dp.j
            if dp.fldir.value is not None and dp.fldir.value > 0:
                # Retrieve the outflow point.
                out_dp = model.dr_pt_by_id.get(dp.fldir.value)
                if out_dp is not None:
                    i_out, j_out = out_dp.i, out_dp.j
                    # In Fortran, i_mat = i_curr*2-1 + (i_out - i_curr), etc.
                    i_mat = i_curr * 2 - 1 + (i_out - i_curr)
                    j_mat = j_curr * 2 - 1 + (j_out - j_curr)
                    # Update the corresponding cell in mat_id with the main channel id.
                    if model.mat_id[i_mat, j_mat] is not None:
                        model.mat_id[i_mat, j_mat].id_pnt = main_net.id_ch.value

        # For the current channel, set the downslope path identifier:
        main_net.n_path = 1
        main_net.id_path.value = [main_net.id_ch.value]
        
        # Now call the recursive routine to update tributaries.
        up_recurs(main_channel.value, model)
        
    
    finish_time = time.process_time()
    elapsed_time = finish_time - start_time
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = elapsed_time % 60
    print(f"Elapsed time (dpl): {hours}h {minutes}m {seconds:.2f}s")


def up_recurs(curr, model):
    """
    Recursively processes tributary channels for downslope path length updates.
    For each tributary channel of the channel identified by 'curr', this function:
      - Updates the downslope path length (dpl) for its drainage points,
      - Updates the associated matrix (mat_id),
      - Propagates the cumulative path (id_path) from the parent channel,
      - Sets the endorheic id to be the same as the parent's,
      - Recursively calls itself on the tributary channel.
    
    Parameters
    ----------
    curr : int
        The channel id (id_ch) of the current (parent) drainage network.
    model : object
        The hydrological model object containing dr_pt, dr_net, etc.
    """
    curr_net = model.dr_net_by_id[curr]
    # Loop over each tributary channel of the current channel.
    # In Fortran: do i=1, dr_net(curr)%n_jun, where dr_net(curr)%id_in holds tributary channels.
    for trib in curr_net.id_in.value:
        # Here, trib is the channel id of the tributary.
        trib_net = model.dr_net_by_id[trib]
        # For each drainage point in the tributary channel except the last one:
        for dp_id in trib_net.id_pnts.value[:-1]:
            dp = model.dr_pt_by_id[dp_id]
            # l1 is the length of the tributary channel.
            l1 = trib_net.length
            l2 = dp.upl
            # l3 is the downslope length at the end point of the tributary channel.
            end_dp = model.dr_pt_by_id[trib_net.id_end_pt.value]
            l3 = end_dp.dpl
            dp.dpl = l1 - l2 + l3
            i_curr, j_curr = dp.i, dp.j
            if dp.fldir.value is not None and dp.fldir.value > 0:
                out_dp = model.dr_pt_by_id.get(dp.fldir)
                if out_dp is not None:
                    i_out, j_out = out_dp.i, out_dp.j
                    i_mat = i_curr * 2 - 1 + (i_out - i_curr)
                    j_mat = j_curr * 2 - 1 + (j_out - j_curr)
                    if model.mat_id[i_mat, j_mat] is not None:
                        model.mat_id[i_mat, j_mat].id_pnt = trib_net.id_ch.value
        # Update the tributary channel’s path:
        trib_net.n_path = curr_net.n_path + 1
        # Build new id_path: first element is its own id_ch, then append parent's path.
        trib_net.id_path.value = [trib_net.id_ch.value] + curr_net.id_path.value
        # Propagate the endorheic id from the parent.
        trib_net.id_endo = curr_net.id_endo
        # Recursively process tributary's tributaries.
        up_recurs(trib, model)

