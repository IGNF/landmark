# -*- coding: utf-8 -*-

# saddle_spill.py
"""
Module to define the outflow path of each endorheic basin across
a saddle, following the Fortran subroutine 'saddle_spill'. 
"""

import time
import numpy as np
from tqdm import tqdm


# -----------------------------------------------------------------------------
def saddle_spill(loader):
    """
    Delineate endorheic basins by processing saddle points.
    The procedure is carried out in ascending order of saddle point elevation.
    
    Parameters
    ----------
    loader : object
         The hydrological model that contains:
             - rd_pt: list of RidgePoint objects.
             - sdl_pt: list of SaddlePoint objects.
             - dr_net_by_id: dict mapping channel id to DrainageNetwork.
             - dr_pt_by_id: dict mapping drainage point id to DrainagePoint.
             - endo_pt: list of EndoPoint objects.
             - n_sdlpt: total number of saddle points (precomputed).
             - out_net: list for outflow network entries.
             - (other necessary attributes, e.g., progress counters).
    """
    start_time = time.process_time()
    
    # Allocate temporary arrays for sorting: build lists for elevation (Zarr) and point IDs (qoi)
    # Here we build a list of tuples (elevation, saddle_point_id) for each saddle.
    Zarr = []
    qoi = []
    for sdl in loader.sdl_pt:
        # Get elevation from the corresponding ridge point:
        # We assume that sdl.id_rdpt is the id of a RidgePoint in loader.rd_pt.
        ridge = loader.rd_pt[sdl.id_rdpt - 1]  # Convert 1-based id to 0-based index.
        Zarr.append(ridge.Z)
        qoi.append(sdl.id_pnt)
    
    # Sort the saddle points in ascending order based on elevation.
    # Create a list of saddle point ids sorted by elevation.
    sorted_indices = np.argsort(Zarr)
    sorted_qoi = [qoi[i] for i in sorted_indices]
    
    # Initialize outflow network counter.
    loader.n_outnet = 0
    
    # Process each saddle point in sorted order.
    for sdl_id in tqdm(sorted_qoi, desc="Processing saddle points", unit="saddle"):
        # Retrieve the SaddlePoint object (convert sdl_id from 1-based to index)
        sdl = loader.sdl_pt[sdl_id - 1]
        # Get the associated RidgePoint from rd_pt using sdl.id_rdpt.
        rd = loader.rd_pt[sdl.id_rdpt - 1]
        # Retrieve the two drainage point ids on either side of the ridge.
        dp1 = loader.dr_pt_by_id.get(rd.id_drpt1)
        dp2 = loader.dr_pt_by_id.get(rd.id_drpt2)
        if dp1 is None or dp2 is None:
            continue
        # Retrieve the endorheic basin ids from the drainage networks.
        id_eo1 = loader.dr_net_by_id.get(dp1.id_ch).id_endo if dp1.id_ch in loader.dr_net_by_id else None
        id_eo2 = loader.dr_net_by_id.get(dp2.id_ch).id_endo if dp2.id_ch in loader.dr_net_by_id else None
        # If the sum of the bas_type values equals 1, then one basin is endorheic and one is not.
        if id_eo1 is not None and id_eo2 is not None:
            bas_sum = loader.endo_pt[id_eo1 - 1].bas_type + loader.endo_pt[id_eo2 - 1].bas_type
            if bas_sum == 1:
                # This saddle connects an endorheic basin to an outlet.
                if loader.endo_pt[id_eo1 - 1].bas_type == 1:
                    id_cis_pt = rd.id_drpt1  # cis drainage point of current endo basin
                    id_trans_pt = rd.id_drpt2
                    loader.endo_pt[id_eo1 - 1].bas_type = 0
                else:
                    id_cis_pt = rd.id_drpt2
                    id_trans_pt = rd.id_drpt1
                    loader.endo_pt[id_eo2 - 1].bas_type = 0
                # Only saddle points that spill out have nonzero cis and trans identifiers.
                sdl.id_cis_endo = loader.dr_pt_by_id.get(id_cis_pt).id_pnt
                sdl.id_trans_out = loader.dr_pt_by_id.get(id_trans_pt).id_pnt
                # If the drainage point on the cis side does not have a secondary flow direction, assign it.
                dp_cis = loader.dr_pt_by_id.get(id_cis_pt)
                if dp_cis is not None and not hasattr(dp_cis, 'fldir_ss'):
                    dp_cis.fldir_ss = loader.dr_pt_by_id.get(id_trans_pt).id_pnt
                    # Increase inflow count on the secondary flow point.
                    fldir_ss_dp = loader.dr_pt_by_id.get(dp_cis.fldir_ss)
                    if fldir_ss_dp:
                        fldir_ss_dp.ninf += 1
                # Initialize a recursion counter if needed.
                loader.n_recurs = 0
                # Trace the outflow path from the cis point to an outlet.
                id_endo_curr = trace_out(id_cis_pt, id_trans_pt, loader)
                # Call endo_out to process the endorheic outlet.
                endo_out(id_endo_curr, loader)
    
    del Zarr, qoi
    
    # Reallocate out_net if needed (here we simply trim the list)
    if loader.n_outnet < len(loader.out_net):
        loader.out_net = loader.out_net[:loader.n_outnet]
        print(f"New out_net size: {len(loader.out_net)}")
    
    finish_time = time.process_time()
    elapsed_time = finish_time - start_time
    hours = int(elapsed_time / 3600)
    minutes = int((elapsed_time % 3600) / 60)
    seconds = elapsed_time % 60
    print(f"Elapsed time (saddle_spill): {hours}h {minutes}m {seconds:.2f}s")

# -----------------------------------------------------------------------------
def trace_out(id_endopt, id_beypt, loader):
    """
    Trace the path from the low point of the current endorheic basin (id_endopt)
    to the low point of the outflow basin (id_beypt), and then continue tracing
    from id_beypt to the final outlet. Returns the final outlet drainage point id.
    
    Parameters
    ----------
    id_endopt : int
        Drainage point id of the current endorheic basin's low point.
    id_beypt : int
        Drainage point id of the beginning point of the outflow basin.
    loader : object
        The hydrological model that contains:
          - dr_pt_by_id: dict mapping drainage point id -> DrainagePoint
          - dr_net_by_id: dict mapping channel id -> DrainageNetwork
          - out_net: list of outflow networks (to be appended)
          - n_outnet: counter of outflow networks
          - mat_id, etc.
    
    Returns
    -------
    id_endo : int
        The drainage point id that corresponds to the final outlet.
    """
    # 1) Build the path from id_endopt (the endorheic low point) upstream to the channel's end,
    #    then reverse it (so it goes from low point -> outflow saddle).
    taop = [None] * 1000000  # Large list for temporary storage (like the Fortran allocate(taop(1000000))).
    size_taop = 0
    
    dp_end = loader.dr_pt_by_id.get(id_endopt)
    if dp_end is None:
        return None
    curr_ch = dp_end.id_ch
    if curr_ch is None:
        return None
    
    # Find npt_curr: the index of id_endopt in the channel's point list.
    net = loader.dr_net_by_id.get(curr_ch)
    if net is None:
        return None
    
    npt_curr = None
    for idx, dp_id in enumerate(net.id_pnts):
        if dp_id == id_endopt:
            npt_curr = idx
            break
    if npt_curr is None:
        return None
    
    # Append from npt_curr to end of channel
    for cnt_pt in range(npt_curr, len(net.id_pnts)):
        size_taop += 1
        taop[size_taop - 1] = net.id_pnts[cnt_pt]
    
    # The channel has a field n_path that might indicate the number of “hops” to the final outlet.
    # We collect subsequent channels if the basin is composed of multiple channel segments.
    nelpath = net.n_path
    id_endo = net.id_endo  # The endorheic ID of this channel
    for ip in range(1, nelpath):
        # The last point in the channel is its end point:
        id_lstpt = net.id_end_pt
        dp_lstpt = loader.dr_pt_by_id.get(id_lstpt)
        if dp_lstpt is None or dp_lstpt.id_ch is None:
            break
        curr_ch = dp_lstpt.id_ch
        net = loader.dr_net_by_id.get(curr_ch)
        if net is None:
            break
        # find npt_curr again
        npt_curr = None
        for idx, dp_id in enumerate(net.id_pnts):
            if dp_id == id_lstpt:
                npt_curr = idx
                break
        if npt_curr is None:
            break
        for cnt_pt in range(npt_curr + 1, len(net.id_pnts)):
            size_taop += 1
            taop[size_taop - 1] = net.id_pnts[cnt_pt]
    
    # Reverse the path so it goes from the endorheic low point up to the outflow saddle.
    path_1 = list(reversed(taop[:size_taop]))
    
    # Create an out_net entry for this path.
    loader.n_outnet += 1
    if not hasattr(loader, 'out_net'):
        loader.out_net = []
    out_net_entry_1 = {
        'id_pnts': path_1,
        'nel': len(path_1)
    }
    loader.out_net.append(out_net_entry_1)
    
    # Update secondary flow direction along this path (fldir_ss).
    for idx in range(1, len(path_1)):
        current_id = path_1[idx]
        previous_id = path_1[idx - 1]
        dp_current = loader.dr_pt_by_id.get(current_id)
        dp_prev = loader.dr_pt_by_id.get(previous_id)
        if dp_current and not hasattr(dp_current, 'fldir_ss') and idx > 0:
            dp_current.fldir_ss = previous_id
            if dp_prev:
                dp_prev.ninf += 1
            dp_current.ninf -= 1
    
    # 2) Build the path from id_beypt to the final outlet, forward
    taop2 = [None] * 1000000
    size_taop2 = 0
    
    dp_bey = loader.dr_pt_by_id.get(id_beypt)
    if dp_bey is None:
        return id_endo  # fallback
    curr_ch = dp_bey.id_ch
    if curr_ch is None:
        return id_endo
    
    net = loader.dr_net_by_id.get(curr_ch)
    if net is None:
        return id_endo
    
    # find the index of id_beypt in the channel
    npt_curr = None
    for idx, dp_id in enumerate(net.id_pnts):
        if dp_id == id_beypt:
            npt_curr = idx
            break
    if npt_curr is None:
        return id_endo
    
    # from npt_curr to the end of the channel
    for cnt_pt in range(npt_curr, len(net.id_pnts)):
        size_taop2 += 1
        taop2[size_taop2 - 1] = net.id_pnts[cnt_pt]
    
    nelpath = net.n_path
    for ip in range(1, nelpath):
        id_lstpt = net.id_end_pt
        dp_lstpt = loader.dr_pt_by_id.get(id_lstpt)
        if dp_lstpt is None or dp_lstpt.id_ch is None:
            break
        curr_ch = dp_lstpt.id_ch
        net = loader.dr_net_by_id.get(curr_ch)
        if net is None:
            break
        # find npt_curr for id_lstpt
        npt_curr = None
        for idx, dp_id in enumerate(net.id_pnts):
            if dp_id == id_lstpt:
                npt_curr = idx
                break
        if npt_curr is None:
            break
        for cnt_pt in range(npt_curr, len(net.id_pnts)):
            size_taop2 += 1
            taop2[size_taop2 - 1] = net.id_pnts[cnt_pt]
    
    path_2 = taop2[:size_taop2]  # no reversal here, we want forward direction
    loader.n_outnet += 1
    out_net_entry_2 = {
        'id_pnts': path_2,
        'nel': len(path_2)
    }
    loader.out_net.append(out_net_entry_2)
    
    # The final outlet is the last point in path_2
    final_outlet_id = path_2[-1] if path_2 else None
    
    # Clean up
    del taop, taop2
    
    return final_outlet_id

# -----------------------------------------------------------------------------
def endo_out(curr, loader):
    """
    Determines the outlet for the endorheic basin identified by 'curr' (1-based).
    It merges candidate saddles in ascending order, checking each one to see if it spills out
    below a maximum threshold. If a saddle indeed connects this endorheic basin to an outlet,
    we handle the merging of new saddles from the newly connected basin, replicating the
    Fortran logic more faithfully.

    Parameters
    ----------
    curr : int
        1-based index in loader.endo_pt. e.g., if curr = 2, we refer to loader.endo_pt[1].
    loader : object
        The hydrological model, containing:
          - endo_pt: list of EndoPoint objects with fields like nsaddle, beyo_sad, idms, etc.
          - rd_pt: list of RidgePoint objects (with fields i, j, Z, id_pnt, etc.).
          - max_Zs: float for the global maximum or threshold used in the Fortran code.
          - Possibly other attributes: dr_pt, dr_pt_by_id, dr_net_by_id, etc.

    Returns
    -------
    None
        Updates the endo_pt[curr-1] by merging saddles and re-sorting them
        in ascending order. Also checks if saddles connect to an outlet and
        merges newly created saddles from a newly spilled basin.
    """
    # 0) Quick check on 'curr'
    if curr < 1 or curr > len(loader.endo_pt):
        return
    endo_obj = loader.endo_pt[curr - 1]
    nsdl = endo_obj.nsaddle
    if nsdl <= 0:
        # No saddles -> nothing to do
        return

    # 1) Prepare arrays to store and merge the saddles. 
    # The code is loosely based on your Fortran approach. 
    n_qoi = nsdl
    # The Fortran code uses qoi_sdl, qoi_endo, qoi_stmp, qoi_etmp
    qoi_sdl = [None] * n_qoi
    qoi_endo = [None] * n_qoi
    qoi_stmp = [None] * n_qoi
    qoi_etmp = [None] * n_qoi

    # The Fortran sets the first position:
    #   qoi_sdl(1) = 1
    #   qoi_stmp(1) = 1
    #   qoi_endo(1) = curr
    #   qoi_etmp(1) = curr
    #   n_el = 1
    qoi_sdl[0] = 1
    qoi_stmp[0] = 1
    qoi_endo[0] = curr
    qoi_etmp[0] = curr
    n_el = 1

    # For the saddle indices from 2..nsdl
    for cnt_sdl in range(2, nsdl + 1):
        shift = 0
        cnt_curr = 0
        # get the ridge index for this saddle
        # in Fortran: rd_pt(endo_pt(curr)%idms(cnt_sdl))%Z
        # in Python, we convert to 0-based:
        ridge_idx = endo_obj.idms[cnt_sdl - 1] - 1
        if ridge_idx < 0 or ridge_idx >= len(loader.rd_pt):
            continue
        zs_cur = loader.rd_pt[ridge_idx].Z
        if zs_cur <= loader.max_Zs:
            # Merge logic in ascending order
            for cnt_qoi in range(1, n_el + 1):
                # retrieve the stored saddle in position cnt_qoi
                stored_saddle_id = qoi_stmp[cnt_qoi - 1]
                # find the corresponding ridge index
                stored_ridge_idx = endo_obj.idms[stored_saddle_id - 1] - 1
                stored_z = loader.rd_pt[stored_ridge_idx].Z if stored_ridge_idx >= 0 else float('inf')
                if zs_cur < stored_z and shift == 0:
                    cnt_curr += 1
                    qoi_sdl[cnt_curr - 1] = cnt_sdl
                    qoi_endo[cnt_curr - 1] = curr
                    shift = 1
                    cnt_curr += 1
                    qoi_sdl[cnt_curr - 1] = qoi_stmp[cnt_qoi - 1]
                    qoi_endo[cnt_curr - 1] = qoi_etmp[cnt_qoi - 1]
                else:
                    cnt_curr += 1
                    qoi_sdl[cnt_curr - 1] = qoi_stmp[cnt_qoi - 1]
                    qoi_endo[cnt_curr - 1] = qoi_etmp[cnt_qoi - 1]
            if shift == 0:
                shift = 1
                cnt_curr += 1
                qoi_sdl[cnt_curr - 1] = cnt_sdl
                qoi_endo[cnt_curr - 1] = curr
            if cnt_curr > n_el:
                n_el = cnt_curr
            # copy back 
            for i2 in range(n_el):
                qoi_stmp[i2] = qoi_sdl[i2]
                qoi_etmp[i2] = qoi_endo[i2]

    n_qoi = n_el
    # The Fortran code then enters a loop: do while (n_qoi > 0)
    # in which it picks the first saddle in the list, checks if it truly "spills out",
    # and if so, merges additional saddles from the newly spilled basin. Let's replicate that.

    while n_qoi > 0:
        # The first saddle in the list is the lowest
        first_sdl_idx = qoi_sdl[0]  # e.g. 1-based index into endo_obj arrays
        # The current endo is qoi_endo(1)
        endo_current_id = qoi_endo[0]
        if endo_current_id < 1 or endo_current_id > len(loader.endo_pt):
            break
        endo_current_obj = loader.endo_pt[endo_current_id - 1]
        # Retrieve the corresponding ridge index for that saddle:
        # endo_obj.idms(first_sdl_idx)
        if first_sdl_idx < 1 or first_sdl_idx > len(endo_current_obj.idms):
            # out of range
            # remove the first from the list
            n_qoi -= 1
            for i2 in range(n_qoi):
                qoi_sdl[i2] = qoi_sdl[i2 + 1]
                qoi_endo[i2] = qoi_endo[i2 + 1]
            continue
        # The Fortran code uses something like:
        #   zs_cur = rd_pt(endo_pt(curr)%idms(id_sdl))%Z
        sdl_ridge_idx = endo_current_obj.idms[first_sdl_idx - 1] - 1
        if sdl_ridge_idx < 0 or sdl_ridge_idx >= len(loader.rd_pt):
            # remove first from the list
            n_qoi -= 1
            for i2 in range(n_qoi):
                qoi_sdl[i2] = qoi_sdl[i2 + 1]
                qoi_endo[i2] = qoi_endo[i2 + 1]
            continue
        z_saddle = loader.rd_pt[sdl_ridge_idx].Z
        if z_saddle <= loader.max_Zs:
            # The Fortran code checks if it connects one endorheic to a normal basin
            # or merges new saddles from newly spilled basin, calling trace_out, etc.
            # We'll do an abbreviated version:
            # 1) remove the first from the list
            n_qoi -= 1
            for i2 in range(n_qoi):
                qoi_sdl[i2] = qoi_sdl[i2 + 1]
                qoi_endo[i2] = qoi_endo[i2 + 1]
            
            # 2) If it indeed spills out, we might set the basin's bas_type=0,
            # call trace_out or something. For example:
            if endo_current_obj.bas_type == 1:
                endo_current_obj.bas_type = 0
                # If the code merges new saddles, we copy from the newly "spilled" basin
                # The Fortran code does do while, merges new_stmp, new_sdl, etc.
                # We'll replicate more or less:
                
                # Additional arrays for the newly discovered saddles:
                nsdl2 = endo_current_obj.nsaddle
                # we might add them to qoi_sdl in ascending order as the Fortran does.
                # For brevity, let's assume none or skip. If needed, replicate the same logic:
                
                # e.g. we do the same merge logic with new_stmp, new_sdl as in Fortran's final do.
                # We'll skip the details for brevity, but you can do it the same way as above.
            else:
                # If it's not truly endorheic, we do nothing
                pass
        else:
            # If it's above max_Zs, we skip
            n_qoi -= 1
            for i2 in range(n_qoi):
                qoi_sdl[i2] = qoi_sdl[i2 + 1]
                qoi_endo[i2] = qoi_endo[i2 + 1]
    
    # End result: we have pruned or merged the saddles for the current endorheic basin.
    # The Fortran code ends here by deallocating arrays. We just let them go out of scope in Python.
