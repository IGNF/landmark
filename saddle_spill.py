# -*- coding: utf-8 -*-

# saddle_spill.py
"""
This module define the outflow path of each endorheic basin across
the saddle that connect it with the lower outflow basin.
The procedure is carried out in saddle point elevation ascending order 
"""

from tqdm import tqdm


def saddle_spill(model):
    """Delineate endorheic basins by processing saddle points.
    The procedure is carried out in ascending order of saddle point elevation.
    
    Parameters
    ----------
    model : object
         The hydrological model that contains:
             - rd_pt: list of RidgePoint objects.
             - sdl_pt: list of SaddlePoint objects.
             - dr_pt : list of DrainagePoint objects.
             - dr_net : list of DrainageNetwork objects.
             - endo_pt: list of EndoPoint objects.
             - n_sdlpt: total number of saddle points (precomputed).
             - out_net: list for outflow network entries.
    """
    
    # Allocate temporary arrays for sorting: build lists for elevation (Zarr) and point IDs (qoi)
    model.out_net = []
    model.n_outnet = 0

    # Build the list 'qoi' by sorting Saddle_point in ascending Z
    # qoi will store the 'id_pnt' for each SaddlePoint
    qoi = [sp.id_pnt for sp in sorted(model.sdl_pt, key=lambda sp: sp.Z)]



    # Process each saddle point in sorted order.
    for cnt_sdl in tqdm(range(len(qoi))):
        id_sdl = qoi[cnt_sdl]
        sp = model.sdl_pt[id_sdl-1]
        rp = model.rd_pt[sp.id_rdpt-1]
        id_eo1 = model.dr_net[model.dr_pt[rp.id_drpt1.value-1].id_ch.value-1].id_endo.value
        id_eo2 = model.dr_net[model.dr_pt[rp.id_drpt2.value-1].id_ch.value-1].id_endo.value
        

        if model.l_endo_pt[id_eo1-1].bas_type + model.l_endo_pt[id_eo2-1].bas_type == 1:
            #This saddle point connect and endorheic basin to an outlet

        
            max_Zs = rp.Z
            if model.l_endo_pt[id_eo1-1].bas_type == 1:
                id_cis_pt = rp.id_drpt1.value #id cis point of the current endo basin
                id_trans_pt = rp.id_drpt2.value #id trans point of the current endo basin
                model.l_endo_pt[id_eo1-1].bas_type = 0
            else: #that means (endo_pt(id_eo2)%bas_type == 1)
                id_cis_pt = rp.id_drpt2.value #id cis point of the current endo basin
                id_trans_pt = rp.id_drpt1.value #id trans point of the current endo basin
                model.l_endo_pt[id_eo2-1].bas_type = 0
            # Only saddle points that spill out have id_cis_endo and in_trans_endo gt 0
            sp.id_cis_endo = model.dr_pt[id_cis_pt-1].id_pnt
            sp.id_trans_out = model.dr_pt[id_trans_pt-1].id_pnt
            if model.dr_pt[id_cis_pt-1].fldir_ss.value == None:
                model.dr_pt[id_cis_pt-1].fldir_ss = model.dr_pt[id_trans_pt-1].id_pnt
                model.dr_pt[model.dr_pt[id_cis_pt-1].fldir_ss.value-1].ninf += 1
                

            
            id_endo_curr = trace_out(id_cis_pt, id_trans_pt, model)
            endo_out(id_endo_curr, max_Zs, model)
    
    del(model.out_net)

            
def trace_out(id_endopt, id_beypt, model):
    """Trace the drainage path from an endorheic point to its outflow into a connected basin,
    then trace the flow path downstream to the final outlet.

    This function updates the model's outflow networks and adjusts subsurface flow directions (`fldir_ss`)
    for each drainage point along the traced paths.

    Parameters
    ----------
    id_endopt : int
        Drainage point ID of the endorheic basin's lowest point.
    id_beypt : int
        Drainage point ID at the beginning of the outflow basin (saddle outflow).
    model : object
        Hydrological model containing:
        - dr_pt : list of DrainagePoint
        - dr_net : list of DrainageNetwork
        - out_net : list
            Stores the list of outflow paths.
        - n_outnet : int
            Number of outflow networks currently stored.
        - mat_id : ndarray
            DEM grid representation.

    Returns
    -------
    id_endo : int
        Drainage point ID corresponding to the final outlet.
    """

    taop = []  # Accumulate point IDs forming the outflow path

    # ▸ Phase 1: From endorheic basin low point to its outlet
    curr_ch = model.dr_pt[id_endopt - 1].id_ch.value
    npt_curr = model.dr_net[curr_ch - 1].id_pnts.value.index(id_endopt)

    # Append remaining points in the current channel starting at id_endopt
    taop.extend(model.dr_net[curr_ch - 1].id_pnts.value[npt_curr:])

    id_endo = model.dr_net[curr_ch - 1].id_endo.value  # Remember endorheic basin ID

    # Trace downstream path across the hierarchical channels
    nelpath = model.dr_net[curr_ch - 1].n_path
    for _ in range(1, nelpath):
        id_lstpt = model.dr_net[curr_ch - 1].id_end_pt.value
        curr_ch = model.dr_pt[id_lstpt - 1].id_ch.value
        npt_curr = model.dr_net[curr_ch - 1].id_pnts.value.index(id_lstpt)
        taop.extend(model.dr_net[curr_ch - 1].id_pnts.value[npt_curr + 1:])

    # Add reversed path to model.out_net (to go upstream → downstream)
    new_outflow = {
        "id_pnts": list(reversed(taop)),
        "nel": len(taop)
    }
    model.out_net.append(new_outflow)
    model.n_outnet += 1

    # ▸ Update fldir_ss for the reversed path
    for cnt_taop in range(len(taop), 1, -1):
        dp = model.dr_pt[taop[cnt_taop - 1] - 1]
        if dp.fldir_ss.value is None:
            dp.fldir_ss.value = taop[cnt_taop - 2]
            model.dr_pt[dp.fldir_ss.value - 1].ninf += 1
            dp.ninf -= 1

    # ▸ Phase 2: From the saddle outflow point to the true outlet
    taop = []

    curr_ch = model.dr_pt[id_beypt - 1].id_ch.value
    npt_curr = model.dr_net[curr_ch - 1].id_pnts.value.index(id_beypt)
    taop.extend(model.dr_net[curr_ch - 1].id_pnts.value[npt_curr:])

    nelpath = model.dr_net[curr_ch - 1].n_path
    for _ in range(1, nelpath):
        id_lstpt = model.dr_net[curr_ch - 1].id_end_pt.value
        curr_ch = model.dr_pt[id_lstpt - 1].id_ch.value
        npt_curr = model.dr_net[curr_ch - 1].id_pnts.value.index(id_lstpt)
        taop.extend(model.dr_net[curr_ch - 1].id_pnts.value[npt_curr + 1:])

    new_outflow = {
        "id_pnts": taop,
        "nel": len(taop)
    }
    model.out_net.append(new_outflow)
    model.n_outnet += 1

    return id_endo

    
                
def endo_out(curr, max_Zs, model):
    """Handles overflow detection and routing for endorheic basins.

    For a given endorheic basin (index `curr`), this function:
    - Sorts all associated saddle points by elevation.
    - Tests each saddle: does it spill over below the threshold elevation (`max_Zs`)?
    - If so, traces the spill path using `trace_out`, merges saddles from the newly reached basin,
      and re-processes them recursively.

    Parameters
    ----------
    curr : int
        1-based index of the endorheic basin in model.l_endo_pt.
    max_Zs : float
        Elevation threshold above which a saddle is not considered spillable.
    model : object
        Hydrological model containing ridge points, drainage networks, endorheic basins, and saddles.

    Returns
    -------
    None
        Updates the model in place (spill routing, subsurface flow direction, saddle priorities).
    """

    # Initialize local storage for saddle sorting
    nsdl = model.l_endo_pt[curr - 1].nsaddle
    n_qoi = nsdl

    qoi_sdl = [0] * n_qoi  # Indices of saddle points
    qoi_stmp = [0] * n_qoi
    qoi_endo = [0] * n_qoi  # Corresponding basin indices
    qoi_endo[0] = curr
    qoi_etmp = [0] * n_qoi
    qoi_etmp[0] = curr
    n_el = 1

    # ▸ Sort all saddle points of current endorheic basin
    for cnt_sdl in range(1, nsdl):
        shift = 0
        cnt_curr = 0
        zs_cur = model.rd_pt[model.l_endo_pt[curr - 1].idms.value[cnt_sdl] - 1].Z

        if zs_cur <= max_Zs:
            for cnt_qoi in range(n_el):
                ref_z = model.rd_pt[model.l_endo_pt[qoi_etmp[cnt_qoi] - 1].idms.value[qoi_stmp[cnt_qoi]] - 1].Z
                if zs_cur < ref_z and shift == 0:
                    cnt_curr += 1
                    qoi_sdl[cnt_curr - 1] = cnt_sdl
                    qoi_endo[cnt_curr - 1] = curr
                    shift = 1
                    cnt_curr += 1
                    qoi_sdl[cnt_curr - 1] = qoi_stmp[cnt_qoi]
                    qoi_endo[cnt_curr - 1] = qoi_etmp[cnt_qoi]
                else:
                    cnt_curr += 1
                    qoi_sdl[cnt_curr - 1] = qoi_stmp[cnt_qoi]
                    qoi_endo[cnt_curr - 1] = qoi_etmp[cnt_qoi]

            if shift == 0:
                cnt_curr += 1
                qoi_sdl[cnt_curr - 1] = cnt_sdl
                qoi_endo[cnt_curr - 1] = curr

        if cnt_curr > n_el:
            n_el = cnt_curr

        qoi_stmp[:n_el] = qoi_sdl[:n_el]
        qoi_etmp[:n_el] = qoi_endo[:n_el]

    qoi_sdl[:n_el] = qoi_stmp[:n_el]
    qoi_endo[:n_el] = qoi_etmp[:n_el]
    n_qoi = n_el

    # ▸ Main processing loop: evaluate each saddle in order
    while n_qoi > 0:
        id_sdl = qoi_sdl[0]
        curr = qoi_endo[0]
        rp = model.rd_pt[model.l_endo_pt[curr - 1].idms.value[id_sdl] - 1]
        zs_cur = rp.Z

        # Remove the current entry from the top of the priority queue
        qoi_stmp[:n_qoi - 1] = qoi_sdl[1:n_qoi]
        qoi_etmp[:n_qoi - 1] = qoi_endo[1:n_qoi]
        n_el = n_qoi - 1

        if zs_cur <= max_Zs:
            id_eo1 = curr
            id_eo2 = model.l_endo_pt[curr - 1].beyo_sad.value[id_sdl]

            # Check whether one of the two connected basins is still endorheic
            if model.l_endo_pt[id_eo1 - 1].bas_type + model.l_endo_pt[id_eo2 - 1].bas_type == 1:
                # Convert endorheic basin into open (bas_type = 0)
                if model.l_endo_pt[id_eo1 - 1].bas_type == 1:
                    model.l_endo_pt[id_eo1 - 1].bas_type = 0
                    if model.dr_net[model.dr_pt[rp.id_drpt1.value - 1].id_ch.value - 1].id_endo.value == id_eo1:
                        id_cis_pt = rp.id_drpt1.value
                        id_trans_pt = rp.id_drpt2.value
                    else:
                        id_cis_pt = rp.id_drpt2.value
                        id_trans_pt = rp.id_drpt1.value
                else:
                    model.l_endo_pt[id_eo2 - 1].bas_type = 0
                    if model.dr_net[model.dr_pt[rp.id_drpt1.value - 1].id_ch.value - 1].id_endo.value == id_eo2:
                        id_cis_pt = rp.id_drpt1.value
                        id_trans_pt = rp.id_drpt2.value
                    else:
                        id_cis_pt = rp.id_drpt2.value
                        id_trans_pt = rp.id_drpt1.value

                # Record connection information in saddle point
                model.sdl_pt[rp.id_sdl - 1].id_cis_endo = model.dr_pt[id_cis_pt - 1].id_pnt
                model.sdl_pt[rp.id_sdl - 1].id_trans_out = model.dr_pt[id_trans_pt - 1].id_pnt

                # Set subsurface direction if not already set
                if model.dr_pt[id_cis_pt - 1].fldir_ss.value is None:
                    model.dr_pt[id_cis_pt - 1].fldir_ss = model.dr_pt[id_trans_pt - 1].id_pnt
                    model.dr_pt[model.dr_pt[id_cis_pt - 1].fldir_ss.value - 1].ninf += 1

                # Trace and process the new outflow from this basin
                id_endo_curr = trace_out(id_cis_pt, id_trans_pt, model)
                curr = id_endo_curr

                # Recollect saddle candidates from this new basin
                nsdl = model.l_endo_pt[curr - 1].nsaddle
                new_sdl = [0] * nsdl
                new_stmp = [0] * nsdl
                new_endo = [0] * nsdl
                new_endo[0] = curr
                new_etmp = [0] * nsdl
                new_etmp[0] = curr
                n_el_new = 1

                for cnt_sdl in range(1, nsdl):
                    shift = 0
                    cnt_curr = 0
                    zs_cur = model.rd_pt[model.l_endo_pt[curr - 1].idms.value[cnt_sdl] - 1].Z
                    if zs_cur <= max_Zs:
                        for cnt_qoi in range(n_el_new):
                            ref_z = model.rd_pt[model.l_endo_pt[new_etmp[cnt_qoi] - 1].idms.value[new_stmp[cnt_qoi]] - 1].Z
                            if zs_cur < ref_z and shift == 0:
                                cnt_curr += 1
                                new_sdl[cnt_curr - 1] = cnt_sdl
                                new_endo[cnt_curr - 1] = curr
                                shift = 1
                                cnt_curr += 1
                                new_sdl[cnt_curr - 1] = new_stmp[cnt_qoi]
                                new_endo[cnt_curr - 1] = new_etmp[cnt_qoi]
                            else:
                                cnt_curr += 1
                                new_sdl[cnt_curr - 1] = new_stmp[cnt_qoi]
                                new_endo[cnt_curr - 1] = new_etmp[cnt_qoi]
                        if shift == 0:
                            cnt_curr += 1
                            new_sdl[cnt_curr - 1] = cnt_sdl
                            new_endo[cnt_curr - 1] = curr
                    if cnt_curr > n_el_new:
                        n_el_new = cnt_curr
                    new_stmp[:n_el_new] = new_sdl[:n_el_new]
                    new_etmp[:n_el_new] = new_endo[:n_el_new]

                # Merge newly collected saddle queue with remaining old one
                cnt_new = n_el_new
                new_sdl[:n_el_new] = new_stmp[:n_el_new]
                new_endo[:n_el_new] = new_etmp[:n_el_new]

                cnt_curr = 0
                mem_new = 0
                qoi_sdl = [0] * (n_el + cnt_new)
                qoi_endo = [0] * (n_el + cnt_new)
                for cnt_qoi in range(n_el):
                    for cnt_sdl in range(mem_new, cnt_new):
                        zs_cur = model.rd_pt[model.l_endo_pt[new_endo[cnt_sdl] - 1].idms.value[new_sdl[cnt_sdl]] - 1].Z
                        ref_z = model.rd_pt[model.l_endo_pt[qoi_etmp[cnt_qoi] - 1].idms.value[qoi_stmp[cnt_qoi]] - 1].Z
                        if zs_cur < ref_z:
                            cnt_curr += 1
                            qoi_sdl[cnt_curr - 1] = new_sdl[cnt_sdl]
                            qoi_endo[cnt_curr - 1] = new_endo[cnt_sdl]
                            mem_new += 1
                    cnt_curr += 1
                    qoi_sdl[cnt_curr - 1] = qoi_stmp[cnt_qoi]
                    qoi_endo[cnt_curr - 1] = qoi_etmp[cnt_qoi]
                if mem_new < cnt_new:
                    for cnt_sdl in range(mem_new, cnt_new):
                        cnt_curr += 1
                        qoi_sdl[cnt_curr - 1] = new_sdl[cnt_sdl]
                        qoi_endo[cnt_curr - 1] = new_endo[cnt_sdl]

                n_el = cnt_curr
                qoi_stmp[:n_el] = qoi_sdl[:n_el]
                qoi_etmp[:n_el] = qoi_endo[:n_el]
                n_qoi = n_el

        # Update the global queue after iteration
        qoi_sdl[:n_el] = qoi_stmp[:n_el]
        qoi_endo[:n_el] = qoi_etmp[:n_el]
        n_qoi = n_el
    
