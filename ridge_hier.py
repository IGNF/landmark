# -*- coding: utf-8 -*-
"""The subroutine RIDGE HIERarchization define an ordered ridge network
where the main rdigelines starts from lower and end to
higher ridges points
"""

from tqdm import tqdm

from data_structures import RidgeNetwork

def ridge_hier(model):
    """
    Builds a hierarchical ridgeline network starting from low elevation ridge points.

    Each ridgeline is traced from lower to higher ridge points using `rdl_trace()`.
    Only points with a single active neighbor and not yet part of a ridgeline
    are used as starting points.

    Parameters
    ----------
    model : object
        Hydrological model containing:
        - rd_pt : list of RidgePoint
        - rd_net : list to be filled with RidgeNetwork objects
        - delta_x : float, grid spacing
        - type_of_landscape : int, controls merge logic (natural/floodplain)
    """

    # Initialize junction count (used for recursive tracing)
    for rp in model.rd_pt:
        rp.n_jun = rp.nen  # Store initial number of neighbors

    # Sort ridge points by elevation (Z), ascending
    qoi = [rp.id_pnt for rp in sorted(model.rd_pt, key=lambda rp: rp.Z)]

    n_rdnet = 0  # Counter for ridgeline IDs
    model.rd_net = []

    # Loop over all ridge points (from lowest to highest)
    for cnt_qoi in tqdm(range(len(qoi))):
        id_nxpt_first = qoi[cnt_qoi]
        curr_pt_nen = model.rd_pt[id_nxpt_first - 1].nen
        Zcurr = model.rd_pt[id_nxpt_first - 1].Z

        # Start new ridgeline only from ridge points with one neighbor
        if curr_pt_nen == 1 and model.rd_pt[id_nxpt_first - 1].Z <= Zcurr:
            flag = 1
            while flag == 1:
                flag = 0
                id_nxpt_first, flag_back, n_rdnet = rdl_trace(id_nxpt_first, n_rdnet, Zcurr, model)
                if flag_back == 1:
                    flag = 1



def rdl_trace(id_rd, n_rdnet, Zcurr, model):
    """
    Recursive function that traces a ridgeline starting from a given ridge point.

    The procedure stops if:
    - the next point has more than one neighbor still available, or
    - the elevation of the next point exceeds Zcurr (initial starting Z).

    Parameters
    ----------
    id_rd : int
        Starting ridge point ID.
    n_rdnet : int
        Current ridgeline counter.
    Zcurr : float
        Elevation threshold for tracing (initial elevation).
    model : object
        Contains the ridge point list and existing ridgelines.

    Returns
    -------
    id_nxpt : int
        Next ridge point ID.
    flag_new : int
        Whether to continue tracing recursively (1) or stop (0).
    n_rdnet : int
        Updated count of ridgelines.
    """

    flag_new = 0
    mem_junc = 0
    curr_rp = model.rd_pt[id_rd - 1]

    if curr_rp.nrdl == 0:
        # ▸ Starting a new ridgeline
        n_rdnet += 1
        new_rl = RidgeNetwork(n_rdnet)
        new_rl.nel = 1
        new_rl.id_pnts = [id_rd]
        if curr_rp.junc == 1:
            new_rl.jun_el = 1

        # Get the only neighbor
        id_nxpt = curr_rp.id_neigh[0]
        nxt_rp = model.rd_pt[id_nxpt - 1]

        # Update new ridgeline
        new_rl.id_pnts.append(id_nxpt)
        new_rl.nel = 2
        new_rl.nrdpt_down = 2

        if nxt_rp.junc == 1 and new_rl.jun_el == 0:
            new_rl.jun_el = 2

        # Compute segment length
        i_start, j_start = curr_rp.i, curr_rp.j
        i_end, j_end = nxt_rp.i, nxt_rp.j
        new_rl.length = model.delta_x * 0.5 * ((i_start - i_end)**2 + (j_start - j_end)**2)**0.5

        # Mark current point as processed
        curr_rp.nen = 0
        curr_rp.id_neigh[0] = None
        curr_rp.nrdl = 1
        curr_rp.id_rdl[0] = new_rl.id_rdl

        # Update next point
        nxt_rp.id_neigh = [n for n in nxt_rp.id_neigh if n != id_rd]
        nxt_rp.nen = len(nxt_rp.id_neigh)
        nxt_rp.n_ptsa += curr_rp.n_ptsa
        nxt_rp.nrdl += 1
        nxt_rp.id_rdl[nxt_rp.nrdl - 1] = new_rl.id_rdl

        # If neighbor still has one connection and is below Zcurr, we continue
        if nxt_rp.nen == 1 and nxt_rp.Z <= Zcurr:
            flag_new = 1

        model.rd_net.append(new_rl)
        return id_nxpt, flag_new, n_rdnet

    else:
        # ▸ Point already belongs to a ridgeline — connect it to the dominant one

        # Choose the ridgeline with the largest upslope area or lowest average Z
        id_mrl = curr_rp.id_rdl[0]
        max_Astart = model.rd_pt[model.rd_net[id_mrl - 1].id_pnts[0] - 1].A_in
        mem_junc = model.rd_net[id_mrl - 1].jun_el

        # Compute mean elevation
        Zsum = sum(model.rd_pt[pid - 1].Z for pid in model.rd_net[id_mrl - 1].id_pnts)
        max_Zmean = Zsum / model.rd_net[id_mrl - 1].nel

        # Try to replace if other ridgelines have larger upstream area
        for cnt_el in range(1, curr_rp.nrdl):
            id_crl = curr_rp.id_rdl[cnt_el]
            net = model.rd_net[id_crl - 1]
            Zmean_curr = sum(model.rd_pt[pid - 1].Z for pid in net.id_pnts) / net.nel
            if model.type_of_landscape == 0:
                if model.rd_pt[net.id_pnts[0] - 1].A_in > max_Astart:
                    max_Astart = model.rd_pt[net.id_pnts[0] - 1].A_in
                    max_Zmean = Zmean_curr
                    id_mrl = id_crl

        id_mz = id_mrl  # Chosen ridgeline to extend
        id_nxpt = curr_rp.id_neigh[0]

        # Update ridgeline with new point
        model.rd_net[id_mz - 1].Zmean = max_Zmean
        model.rd_net[id_mz - 1].id_pnts.append(id_nxpt)
        model.rd_net[id_mz - 1].nel += 1
        model.rd_net[id_mz - 1].nrdpt_down += 1

        if mem_junc > 0 and model.rd_net[id_mz - 1].jun_el == 0:
            model.rd_net[id_mz - 1].jun_el = model.rd_net[id_mz - 1].nel - 1

        # Update length
        i_start = model.rd_pt[model.rd_net[id_mz - 1].id_pnts[-2] - 1].i
        j_start = model.rd_pt[model.rd_net[id_mz - 1].id_pnts[-2] - 1].j
        i_end = model.rd_pt[model.rd_net[id_mz - 1].id_pnts[-1] - 1].i
        j_end = model.rd_pt[model.rd_net[id_mz - 1].id_pnts[-1] - 1].j
        model.rd_net[id_mz - 1].length += model.delta_x * 0.5 * ((i_start - i_end)**2 + (j_start - j_end)**2)**0.5

        # Update current point
        curr_rp.nen = 0
        curr_rp.id_neigh[0] = None
        curr_rp.nrdl += 1
        curr_rp.id_rdl[curr_rp.nrdl - 1] = model.rd_net[id_mz - 1].id_rdl

        # Update next point
        model.rd_pt[id_nxpt - 1].id_neigh = [n for n in model.rd_pt[id_nxpt - 1].id_neigh if n != id_rd]
        model.rd_pt[id_nxpt - 1].nen = len(model.rd_pt[id_nxpt - 1].id_neigh)
        model.rd_pt[id_nxpt - 1].n_ptsa += model.rd_pt[id_rd - 1].n_ptsa
        model.rd_pt[id_nxpt - 1].id_rdl[model.rd_pt[id_nxpt - 1].nrdl] = model.rd_net[id_mz - 1].id_rdl
        model.rd_pt[id_nxpt - 1].nrdl += 1

        if model.rd_pt[id_nxpt - 1].nen == 1 and model.rd_pt[id_nxpt - 1].Z <= Zcurr:
            flag_new = 1

        return id_nxpt, flag_new, n_rdnet
        

        
        
        
    
    
