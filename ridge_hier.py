# -*- coding: utf-8 -*-
"""The subroutine RIDGE HIERarchization define an ordered ridge network
where the main rdigelines starts from lower and end to
higher ridges points
"""

from tqdm import tqdm

from data_structures import RidgeNetwork

def ridge_hier(model):
    """
    The subroutine RIDGE HIERarchization define an ordered ridge network
    where the main rdigelines starts from lower and end to
    higher ridges points
    """
    
    for rp in model.rd_pt:
        rp.n_jun = rp.nen
    
    # Build the list 'qoi' by sorting Ridge_point in ascending Z
    # qoi will store the 'id_pnt' for each RidgePoint
    qoi = [rp.id_pnt for rp in sorted(model.rd_pt, key=lambda rp: rp.Z)]
    
    n_rdnet=0  #id ridgelines in rd_net
    model.rd_net = []
    for cnt_qoi in tqdm(range(len(qoi))):
        id_nxpt_first = qoi[cnt_qoi]
        curr_pt_nen = model.rd_pt[id_nxpt_first-1].nen
        Zcurr = model.rd_pt[id_nxpt_first-1].Z
        if curr_pt_nen == 1 and model.rd_pt[id_nxpt_first-1].Z <= Zcurr:
            flag = 1
            while flag == 1:
                flag = 0
                id_nexpt_new, flag_back, n_rdnet = rdl_trace(id_nxpt_first, n_rdnet, Zcurr, model)
                if flag_back == 1:
                    id_nxpt_first = id_nexpt_new
                    flag = flag_back


def rdl_trace(id_rd, n_rdnet, Zcurr, model):
    """RiDgeLine TRACEr
    
    The subruotine trace recursively rdgelines. The procedure stop when 
    meet a ridge point with more then 1 neighbour still disconnected or the
    elevation is more then Zcurr
    The recursivity needs to fill ridgelines depressions and saddles
    """
    flag_new = 0
    mem_junc = 0
    curr_rp = model.rd_pt[id_rd-1]
    
    if curr_rp.nrdl == 0:
        # if curr_rp.id_pnt == 35745:
        #     print(f"passage du ridge point {curr_rp.id_pnt} dans le if  avec une liste id_neigh = {curr_rp.id_neigh}")

        #first point of new ridgeline
        n_rdnet += 1
        new_rl = RidgeNetwork(n_rdnet)
        new_rl.nel = 1
        new_rl.id_pnts = [id_rd]
        if curr_rp.junc == 1:
            new_rl.jun_el = 1
        id_nxpt = curr_rp.id_neigh[0] #id next point. The current rd_pnt have to be only one neighbour point
        nxt_rp = model.rd_pt[id_nxpt-1]
        new_rl.n_down = 0
        new_rl.Zmean_down = 0
        new_rl.nel = 2
        new_rl.nrdpt_down = 2
        new_rl.id_pnts.append(id_nxpt)
        
        if nxt_rp.junc == 1  and new_rl.jun_el == 0:
            new_rl.jun_el = 2
        
        i_start = model.rd_pt[id_rd-1].i
        j_start = model.rd_pt[id_rd-1].j
        i_end = model.rd_pt[id_nxpt-1].i
        j_end = model.rd_pt[id_nxpt-1].j
        new_rl.length = model.delta_x*0.5*((i_start-i_end)**2+(j_start-j_end)**2)**0.5
        curr_rp.nen = 0
        if curr_rp.md == 1:
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
        
        return id_nxpt,flag_new, n_rdnet
    
    else:
        # min_Zstart = model.rd_pt[model.rd_net[curr_rp.id_rdl[0]-1].id_pnts[0]-1].Z
        # print("curr_rp.id_rdl[0] : ", curr_rp.id_rdl[0])
        # print("curr_rp.id_rdl : ", curr_rp.id_rdl)
        # print("model.rd_net[curr_rp.id_rdl[0]-1].id_pnts[0] : ", model.rd_net[curr_rp.id_rdl[0]-1].id_pnts[0])
        
        max_Astart = model.rd_pt[model.rd_net[curr_rp.id_rdl[0]-1].id_pnts[0]-1].A_in
        
        mem_junc = model.rd_net[curr_rp.id_rdl[0]-1].jun_el
        id_mrl = curr_rp.id_rdl[0] #id ridgeline max mutual distance, first assignement
        # prev_rdpt = model.rd_net[id_mrl-1].nel-1 #penultimate ridge point on max ridgeline
        # max_md = model.rd_pt[model.rd_net[id_mrl-1].id_pnts[prev_rdpt-1]-1].md #max mutual distance first assignement
        
        #Calculate the average altitude, it is valid for the plain areas.
        Zsum = 0
        Zmean_down = 0
        n_down = 0
        nrdpt_down = 0
        for cnt_pt in range(model.rd_net[id_mrl-1].nel):
            Zsum += model.rd_pt[model.rd_net[id_mrl-1].id_pnts[cnt_pt]-1].Z
        max_Zmean = Zsum / model.rd_net[id_mrl-1].nel
        # n_down1 = model.rd_net[id_mrl-1].n_down
        # Zmean1 = (Zsum - model.rd_pt[model.rd_net[id_mrl-1].id_pnts[model.rd_net[id_mrl-1].nel-1]-1].Z) / (model.rd_net[id_mrl-1].nel-1)
        # max_Zmean is the average of the quotas that make up the current ridge
        # I also have to summarize the average of the quotas of the afferent ridges
        # I only need it if the current ridge is not the one that continues
        # in this case the ridge becomes afferent to the one that continues and I have to add max_Zmean with Zmean_down
        
        # n_down1 += 1

        # if n_down1 > 0:
        #     # Zmean_down1 = ((model.rd_net[id_mrl-1].Zmean_down*n_down1)+Zmean1)/(n_down1 +1)
        #     n_down1 += 1
        #     # nrdpt_down1 = model.rd_net[id_mrl-1].nrdpt_down - 1
        # else:
        #     # Zmean_down1 = Zmean1
        #     n_down1 += 1
        #     # nrdpt_down1 = model.rd_net[id_mrl-1].nrdpt_down - 1
        
        if curr_rp.junc == 1 and model.rd_net[id_mrl-1].jun_el == 0 :
            model.rd_net[id_mrl-1].jun_el = model.rd_net[id_mrl-1].nel
        
        for cnt_el in range(1, curr_rp.nrdl):
            #replaced the min_Zstart criterion with max_mutual_distance
            id_crl = curr_rp.id_rdl[cnt_el] #id current ridgeline mutual distance
            # prev_rdpt = model.rd_net[id_crl-1].nel - 1 #penultimate ridge point belonging to the current ridgeline
            # cur_md = model.rd_pt[model.rd_net[id_crl-1].id_pnts[prev_rdpt-1]-1].md #current mutual distance 
            #I calculate the average altitude, it is valid for the plain areas.
            Zsum = 0
            # flg_switch = 0
            for cnt_pt in range(model.rd_net[id_crl-1].nel):
                Zsum += model.rd_pt[model.rd_net[id_crl-1].id_pnts[cnt_pt]-1].Z
            Zmean_curr = Zsum / model.rd_net[id_crl-1].nel
            # Here too I have to calculate Zmean_down
            # which is ready in case I need it
            # I remove the last point from the mean
            # Zmean2 = (Zsum - model.rd_pt[model.rd_net[id_crl-1].id_pnts[model.rd_net[id_crl-1].nel-1]-1].Z) / (model.rd_net[id_crl-1].nel-1)
            # n_down2 = model.rd_net[id_crl-1].n_down
            
            # n_down2 += 1

            # if n_down2 > 0:
            #     # Zmean_down2 = ((model.rd_net[id_crl-1].Zmean_down * n_down2) + Zmean2) / (n_down2 + 1)
            #     n_down2 += 1
            #     # nrdpt_down2 = model.rd_net[id_crl-1].nrdpt_down - 1
            # else:
            #     # Zmean_down2 = Zmean2
            #     n_down2 += 1
            #     # nrdpt_down2 = model.rd_net[id_crl-1].nrdpt_down - 1
            
            if model.rd_net[curr_rp.id_rdl[cnt_el]-1].jun_el > 0:
                mem_junc = model.rd_net[curr_rp.id_rdl[cnt_el]-1].jun_el
            
            #Natural basin (zero), Flood Plane (one)
            if model.type_of_landscape == 0:
                #natural basins
                if model.rd_pt[model.rd_net[curr_rp.id_rdl[cnt_el]-1].id_pnts[0]-1].A_in > max_Astart:
                    max_Astart = model.rd_pt[model.rd_net[curr_rp.id_rdl[cnt_el]-1].id_pnts[0]-1].A_in
                    max_Zmean = Zmean_curr
                    # max_md = cur_md
                    id_mrl = id_crl
                
        id_mz = id_mrl
        #At this point %id_neigh have to be only one element
        id_nxpt = curr_rp.id_neigh[0] #id next point. The current rd_pnt have to be only one neighbour point
        # if curr_rp.id_pnt == 35745:
        #     print(f"passage du ridge point {curr_rp.id_pnt} avec une liste id_neigh = {curr_rp.id_neigh}")
        # if id_nxpt == None:
        #     print(f"\nid_nxpt est none pour le point {curr_rp.id_pnt}")
        #     print(f"Pour le ridge point {curr_rp.id_pnt} la liste id_neigh = {curr_rp.id_neigh}")
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
        
        # Mise à jour de la liste en excluant l'élément avec id_rd
        model.rd_pt[id_nxpt-1].id_neigh = [
            neigh for neigh in model.rd_pt[id_nxpt-1].id_neigh if neigh != id_rd
        ]
        
        model.rd_pt[id_nxpt-1].nen = len(model.rd_pt[id_nxpt-1].id_neigh)
        model.rd_pt[id_nxpt-1].n_ptsa += model.rd_pt[id_rd-1].n_ptsa
        
        # model.rd_pt[id_nxpt-1].id_rdl.append(model.rd_net[id_mz-1].id_rdl)
        model.rd_pt[id_nxpt-1].id_rdl[model.rd_pt[id_nxpt-1].nrdl] = model.rd_net[id_mz-1].id_rdl

        model.rd_pt[id_nxpt-1].nrdl += 1
        
        curr_pt_nen = model.rd_pt[id_nxpt-1].nen
        if curr_pt_nen == 1 and model.rd_pt[id_nxpt-1].Z <= Zcurr:
            flag_new = 1
    
        return id_nxpt,flag_new, n_rdnet
        

        
        
        
    
    
