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
        if curr_pt_nen == 1 and model.rd_pt[id_nxpt_first-1].Z <= Zcurr:
            flag = 1
            while flag == 1:
                flag = 0
                id_nexpt_new, flag_back, n_rdnet = rdl_trace(id_nxpt_first, n_rdnet, Zcurr, model)
                if flag_back == 1:
                    id_nxpt_first = id_nexpt_new
                    flag = flag_back



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
        new_rl.n_down = 0
        new_rl.Zmean_down = 0
        new_rl.nel = 2
        new_rl.nrdpt_down = 2
        new_rl.id_pnts.append(id_nxpt)

        if nxt_rp.junc == 1 and new_rl.jun_el == 0:
            new_rl.jun_el = 2

        i_start = model.rd_pt[id_rd-1].i
        j_start = model.rd_pt[id_rd-1].j
        i_end = model.rd_pt[id_nxpt-1].i
        j_end = model.rd_pt[id_nxpt-1].j
        new_rl.length = model.delta_x*0.5*((i_start-i_end)**2+(j_start-j_end)**2)**0.5
        curr_rp.nen = 0
        # if curr_rp.md == 1:
        if curr_rp.id_sdl != None:
            curr_rp.n_ptsa = 0 #splling saddle
        curr_rp.id_neigh[0] = None
        curr_rp.nrdl = 1
        curr_rp.id_rdl[0] = new_rl.id_rdl
        
        # Keep only neighbors not equal to id_rd
        nxt_rp.id_neigh = [
        neighbor for neighbor in nxt_rp.id_neigh 
        if neighbor != id_rd
        ]
        
        # Update neighbor count (nen)
        nxt_rp.nen = len(nxt_rp.id_neigh)
        
        # Update n_ptsa
        nxt_rp.n_ptsa += curr_rp.n_ptsa
        nxt_rp.nrdl += 1
        nxt_rp.id_rdl[nxt_rp.nrdl-1] = new_rl.id_rdl
        curr_pt_nen = nxt_rp.nen
        if curr_pt_nen == 1 and nxt_rp.Z <= Zcurr:
            flag_new = 1
        
        model.rd_net.append(new_rl)
        
        return id_nxpt, flag_new, n_rdnet

    else:
        # ▸ Point already belongs to a ridgeline — connect it to the dominant one

        max_Astart = model.rd_pt[model.rd_net[curr_rp.id_rdl[0]-1].id_pnts[0]-1].A_in
        
        mem_junc = model.rd_net[curr_rp.id_rdl[0]-1].jun_el
        id_mrl = curr_rp.id_rdl[0] #id ridgeline max mutual distance, first assignement
        
        #Calculate the average altitude, it is valid for the plain areas.
        Zsum = 0
        Zmean_down = 0
        n_down = 0
        nrdpt_down = 0
        for cnt_pt in range(model.rd_net[id_mrl-1].nel):
            Zsum += model.rd_pt[model.rd_net[id_mrl-1].id_pnts[cnt_pt]-1].Z
        max_Zmean = Zsum / model.rd_net[id_mrl-1].nel
        
        if curr_rp.junc == 1 and model.rd_net[id_mrl-1].jun_el == 0 :
            model.rd_net[id_mrl-1].jun_el = model.rd_net[id_mrl-1].nel
        
        for cnt_el in range(1, curr_rp.nrdl):
            id_crl = curr_rp.id_rdl[cnt_el] #id current ridgeline mutual distance
            Zsum = 0
            for cnt_pt in range(model.rd_net[id_crl-1].nel):
                Zsum += model.rd_pt[model.rd_net[id_crl-1].id_pnts[cnt_pt]-1].Z
            Zmean_curr = Zsum / model.rd_net[id_crl-1].nel
            
            if model.rd_net[curr_rp.id_rdl[cnt_el]-1].jun_el > 0:
                mem_junc = model.rd_net[curr_rp.id_rdl[cnt_el]-1].jun_el
            
            #Natural basin (zero), Flood Plane (one)
            if model.type_of_landscape == 0:
                #natural basins
                if model.rd_pt[model.rd_net[curr_rp.id_rdl[cnt_el]-1].id_pnts[0]-1].A_in > max_Astart:
                    max_Astart = model.rd_pt[model.rd_net[curr_rp.id_rdl[cnt_el]-1].id_pnts[0]-1].A_in
                    max_Zmean = Zmean_curr
                    id_mrl = id_crl
                
        id_mz = id_mrl
        #At this point %id_neigh have to be only one element
        id_nxpt = curr_rp.id_neigh[0] #id next point. The current rd_pnt have to be only one neighbour point
        model.rd_net[id_mz-1].Zmean = max_Zmean #it updates every time I add a point and recalculate Zmean
        if model.rd_net[id_mz-1].n_down + n_down > 0:
            model.rd_net[id_mz-1].Zmean_down = ((model.rd_net[id_mz-1].Zmean_down * model.rd_net[id_mz-1].n_down) + (Zmean_down * n_down)) /\
                (model.rd_net[id_mz-1].n_down + n_down) #it updates every time I add a point and recalculate Zmean
            model.rd_net[id_mz-1].n_down += n_down
            model.rd_net[id_mz-1].nrdpt_down += nrdpt_down
        
        model.rd_net[id_mz-1].id_pnts.append(id_nxpt)
        model.rd_net[id_mz-1].nel += 1
        model.rd_net[id_mz-1].nrdpt_down += 1
        if mem_junc > 0 and model.rd_net[id_mz-1].jun_el == 0:
            model.rd_net[id_mz-1].jun_el = model.rd_net[id_mz-1].nel - 1
        i_start = model.rd_pt[model.rd_net[id_mz-1].id_pnts[model.rd_net[id_mz-1].nel-2]-1].i
        j_start = model.rd_pt[model.rd_net[id_mz-1].id_pnts[model.rd_net[id_mz-1].nel-2]-1].j
        i_end = model.rd_pt[model.rd_net[id_mz-1].id_pnts[model.rd_net[id_mz-1].nel-1]-1].i
        j_end = model.rd_pt[model.rd_net[id_mz-1].id_pnts[model.rd_net[id_mz-1].nel-1]-1].j
        model.rd_net[id_mz-1].length += model.delta_x*0.5*((i_start-i_end)**2+(j_start-j_end)**2)**0.5
        
        curr_rp.nen = 0
        curr_rp.id_neigh[0] = None
        curr_rp.nrdl += 1
        curr_rp.id_rdl[curr_rp.nrdl-1] = model.rd_net[id_mz-1].id_rdl
        
        model.rd_pt[id_nxpt-1].id_neigh = [
            neigh for neigh in model.rd_pt[id_nxpt-1].id_neigh if neigh != id_rd
        ]
        
        model.rd_pt[id_nxpt-1].nen = len(model.rd_pt[id_nxpt-1].id_neigh)
        model.rd_pt[id_nxpt-1].n_ptsa += model.rd_pt[id_rd-1].n_ptsa
        
        model.rd_pt[id_nxpt-1].id_rdl[model.rd_pt[id_nxpt-1].nrdl] = model.rd_net[id_mz-1].id_rdl

        model.rd_pt[id_nxpt-1].nrdl += 1
        
        curr_pt_nen = model.rd_pt[id_nxpt-1].nen
        if curr_pt_nen == 1 and model.rd_pt[id_nxpt-1].Z <= Zcurr:
            flag_new = 1

        return id_nxpt, flag_new, n_rdnet
        

        
        
        
    
    
