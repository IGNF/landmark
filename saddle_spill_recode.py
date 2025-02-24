# -*- coding: utf-8 -*-

# saddle_spill.py
"""
This module define the outflow path of each endorheic basin across
the saddle that connect it with the lower outflow basin.
The procedure is carried out in saddle point elevation ascending order 
"""

import time
import numpy as np
from tqdm import tqdm


# -----------------------------------------------------------------------------
def saddle_spill(model):
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
    
    # Allocate temporary arrays for sorting: build lists for elevation (Zarr) and point IDs (qoi)
    model.out_net = []
    model.n_outnet = 0

    # Build the list 'qoi' by sorting Saddle_point in ascending Z
    # qoi will store the 'id_pnt' for each SaddlePoint
    qoi = [sp.id_pnt for sp in sorted(model.sdl_pt, key=lambda sp: sp.Z)]


    # Process each saddle point in sorted order.
    for cnt_sdl in tqdm(range(len(model.sdl_pt))):
        id_sdl = qoi[cnt_sdl]
        sp = model.sdl_pt[id_sdl-1]
        rp = model.rd_pt[sp.id_rdpt-1]
        id_eo1 = model.l_dr_net[model.dr_pt[rp.id_drpt1.value-1].id_ch.value-1].id_endo.value
        id_eo2 = model.l_dr_net[model.dr_pt[rp.id_drpt2.value-1].id_ch.value-1].id_endo.value
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
    
    for pt in model.sdl_pt:
        print(f"id_pnt : {pt.id_pnt}, id_cis_endo : {pt.id_cis_endo.value}, id_trans_endo : {pt.id_trans_out.value}")

            

# -----------------------------------------------------------------------------
def trace_out(id_endopt, id_beypt, model):
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
    
    taop = []

    # Trouver `curr_ch` et `npt_curr`
    curr_ch = model.dr_pt[id_endopt-1].id_ch.value
    npt_curr = model.l_dr_net[curr_ch-1].id_pnts.value.index(id_endopt) #Index de id_endopt dans la liste id_pnts du canal courant
    
    # Ajouter les points restants de `l_dr_net(curr_ch)` dans `taop`
    taop.extend(model.l_dr_net[curr_ch-1].id_pnts.value[npt_curr:])
    
    id_endo = model.l_dr_net[curr_ch-1].id_endo.value
    # Ajouter les chemins supplémentaires
    nelpath = model.l_dr_net[curr_ch-1].n_path
    #from the first junction to the endoreich basin lowest point 
    for _ in range(1, nelpath):
        id_lstpt = model.l_dr_net[curr_ch-1].id_end_pt.value
        curr_ch = model.dr_pt[id_lstpt-1].id_ch.value
        npt_curr = model.l_dr_net[curr_ch-1].id_pnts.value.index(id_lstpt)
        taop.extend(model.l_dr_net[curr_ch-1].id_pnts.value[npt_curr + 1:])  # Éviter le doublon du dernier point
        
    
    # Ajouter le chemin inversé à `out_net`
    new_outflow = {
        "id_pnts": list(reversed(taop)),
        "nel": len(taop)
    }
    model.out_net.append(new_outflow)
    model.n_outnet += 1
   
    
    # Mise à jour de `fldir_ss`
    for cnt_taop in range(len(taop), 1, -1): 
        dp = model.dr_pt[taop[cnt_taop-1]-1]
        if dp.fldir_ss.value is None:
            # Mise à jour de fldir_ss avec la valeur du point précédent dans `taop
            dp.fldir_ss.value = taop[cnt_taop-2]
            # Incrémentation du nombre d'affluents (`ninf`) du point précédent
            model.dr_pt[dp.fldir_ss.value-1].ninf += 1
            # Décrémentation du `ninf` du point actuel
            dp.ninf -= 1



    # #Code Fortran : 
    # # Boucle en ordre inverse (Fortran `do cnt_taop=size_taop,1,-1`)
    # for cnt_taop in range(len(taop), 0, -1):
    #     # Condition: Vérifier si `fldir_ss` est NULL et que `cnt_taop > 1`
    #     if model.dr_pt[taop[cnt_taop - 1]-1].fldir_ss is None and cnt_taop > 1:
    #         # Mise à jour de fldir_ss avec la valeur du point précédent dans `taop`
    #         model.dr_pt[taop[cnt_taop - 1]-1].fldir_ss = taop[cnt_taop - 2]
    
    #         # Incrémentation du nombre d'affluents (`ninf`) du point précédent
    #         model.dr_pt[model.dr_pt[taop[cnt_taop - 1]-1].fldir_ss-1].ninf += 1
    
    #         # Décrémentation du `ninf` du point actuel
    #         model.dr_pt[taop[cnt_taop - 1]-1].ninf -= 1


    # Path from the current endorheic saddle point to the free basin outlet point:
    taop = []

    # Trouver `curr_ch` et `npt_curr`
    curr_ch = model.dr_pt[id_beypt-1].id_ch.value
    npt_curr = model.l_dr_net[curr_ch-1].id_pnts.value.index(id_beypt) #Index de id_beypt dans la liste id_pnts du canal courant
    
    # Ajouter les points restants de `l_dr_net(curr_ch)` dans `taop`
    taop.extend(model.l_dr_net[curr_ch-1].id_pnts.value[npt_curr:])
    
    
    # Ajouter les chemins supplémentaires
    nelpath = model.l_dr_net[curr_ch-1].n_path
    #from the first junction to the endoreich basin lowest point 
    for _ in range(1, nelpath):
        id_lstpt = model.l_dr_net[curr_ch-1].id_end_pt.value
        curr_ch = model.dr_pt[id_lstpt-1].id_ch.value
        npt_curr = model.l_dr_net[curr_ch-1].id_pnts.value.index(id_lstpt)
        taop.extend(model.l_dr_net[curr_ch-1].id_pnts.value[npt_curr+1:])
        
        
    # Ajouter le chemin inversé à `out_net`
    new_outflow = {
        "id_pnts": taop,
        "nel": len(taop)
    }
    model.out_net.append(new_outflow)
    model.n_outnet += 1

    return id_endo


def endo_out(curr, max_Zs, model):
    """
    Determines the outlet for the endorheic basin identified by 'curr' (1-based).
    It merges candidate saddles in ascending order, checking each one to see if it spills out
    below a maximum threshold. If a saddle indeed connects this endorheic basin to an outlet,

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
    nsdl = model.l_endo_pt[curr-1].nsaddle
    n_qoi = nsdl
    
    qoi_sdl = []
    qoi_stmp = [0]
    qoi_endo = []
    qoi_etmp = [curr]
    n_el = 1
    
    #current endorheic basin saddles are considered 
    for cnt_sdl in range (1, nsdl): #the first one is already put in the tmp 
        shift = 0
        cnt_curr = 0
        zs_cur = model.rd_pt[model.l_endo_pt[curr-1].idms.value[cnt_sdl]-1].Z
        if zs_cur <= max_Zs:
            for cnt_qoi in range(n_el):

                if zs_cur < model.rd_pt[model.l_endo_pt[qoi_etmp[cnt_qoi]-1].idms.value[qoi_stmp[cnt_qoi]]-1].Z and shift ==0 :
                    # the current saddle has a lower quota than the one recorded in the current position
                    # so I insert that and then add the quota already recorded which moves
                    cnt_curr += 1
                    qoi_sdl.append(cnt_sdl)
                    qoi_endo.append(curr)
                    shift = 1
                    cnt_curr += 1
                    qoi_sdl.append(qoi_stmp[cnt_qoi])
                    qoi_endo.append(qoi_etmp[cnt_qoi])
                else :
                    #the current saddle has a higher height, the saddle that was already in the current position remains
                    cnt_curr += 1
                    qoi_sdl.append(qoi_stmp[cnt_qoi])
                    qoi_endo.append(qoi_etmp[cnt_qoi])
            # If shift is null it means that I have not inserted the saddle I have to insert anywhere because evidently
            # the height is greater than the heights of the saddles already present. I therefore put it last
            if shift == 0:
                shift = 1
                cnt_curr += 1
                qoi_sdl.append(cnt_sdl)
                qoi_endo.append(curr)
        if cnt_curr > n_el:
            n_el = cnt_curr
        
        if qoi_sdl == []:
            qoi_sdl = [0]
            qoi_endo = [curr]
            qoi_stmp[:n_el] = qoi_sdl[:n_el]
            qoi_etmp[:n_el] = qoi_endo[:n_el]
            qoi_sdl = []
            qoi_endo = []
            
        else:            
            qoi_stmp[:n_el] = qoi_sdl[:n_el]
            qoi_etmp[:n_el] = qoi_endo[:n_el]
    
    qoi_sdl[:n_el] = qoi_stmp[:n_el]
    qoi_endo[:n_el] = qoi_etmp[:n_el]
    
    n_qoi = n_el
    while n_qoi > 0:
        cnt_sdl = 1 #the first saddle is the one that currently has the lowest share
        id_sdl = qoi_sdl[cnt_sdl-1]
        curr = qoi_endo[cnt_sdl-1] #id of the current elements in endo
        rp = model.rd_pt[model.l_endo_pt[curr-1].idms.value[id_sdl]-1]
        zs_cur = rp.Z
        nsdl = model.l_endo_pt[curr-1].nsaddle
        qoi_stmp = qoi_sdl[1:n_qoi+1]
        qoi_etmp = qoi_endo[1:n_qoi+1]
        n_el = n_qoi -1
        if zs_cur <= max_Zs:
            id_eo1 = curr
            id_eo2 = model.l_endo_pt[curr-1].beyo_sad.value[id_sdl-1] #!!!!!!! pas certain du id_sdl-1
            if model.l_endo_pt[id_eo1-1].bas_type + model.l_endo_pt[id_eo2-1].bas_type == 1:
                #This saddle point connect and endorheic basin to an outlet
                if model.l_endo_pt[id_eo1-1].bas_type == 1:
                    model.l_endo_pt[id_eo1-1].bas_type = 0
                    if model.l_dr_net[model.dr_pt[rp.id_drpt1.value-1].id_ch.value-1].id_endo.value == id_eo1:
                        id_cis_pt = rp.id_drpt1.value #id cis point of the current endo basin
                        id_trans_pt = rp.id_drpt2.value #id trans point of the current endo basin
                    else:
                        id_cis_pt = rp.id_drpt2.value
                        id_trans_pt = rp.id_drpt1.value
                else: #that means (endo_pt(id_eo2)%bas_type == 1)
                    model.l_endo_pt[id_eo2-1].bas_type = 0
                    if model.l_dr_net[model.dr_pt[rp.id_drpt1.value-1].id_ch.value-1].id_endo.value == id_eo2:
                        id_cis_pt = rp.id_drpt1.value
                        id_trans_pt = rp.id_drpt2.value
                    else:
                        id_cis_pt = rp.id_drpt2.value
                        id_trans_pt = rp.id_drpt1.value
                #Only saddle points that spill out have id_cis_endo and id_trans_endo gt 0
                model.sdl_pt[rp.id_sdl-1].id_cis_endo = model.dr_pt[id_cis_pt-1].id_pnt
                model.sdl_pt[rp.id_sdl-1].id_trans_out = model.dr_pt[id_trans_pt-1].id_pnt
                if model.dr_pt[id_cis_pt-1].fldir_ss.value == None:
                    model.dr_pt[id_cis_pt-1].fldir_ss = model.dr_pt[id_trans_pt-1].id_pnt
                    model.dr_pt[model.dr_pt[id_cis_pt-1].fldir_ss.value-1].ninf += 1
                    
                id_endo_curr = trace_out(id_cis_pt, id_trans_pt, model) #trace the path from curent endorheic to outlet
                curr = id_endo_curr
                nsdl = model.l_endo_pt[curr-1].nsaddle
                cnt_new = 0
                new_sdl = []
                new_stmp = [0]
                new_endo = []
                new_etmp = [curr]
                n_el_new = 1
                #consider all the saddles of the current endorheic basin
                for cnt_sdl in range (1,nsdl):
                    shift = 0
                    cnt_curr = 0
                    zs_cur = model.rd_pt[model.l_endo_pt[curr-1].idms.value[cnt_sdl]-1].Z
                    if zs_cur <= max_Zs:
                        for cnt_qoi in range(n_el_new):
                            if zs_cur < model.rd_pt[model.l_endo_pt[new_etmp[cnt_qoi]-1].idms.value[new_stmp[cnt_qoi]]-1].Z and shift ==0 :
                                #the current saddle has a lower quota than the one recorded in the current position
                                #so I insert that and then add the quota already recorded which moves
                                cnt_curr += 1
                                new_sdl.append(cnt_sdl)
                                new_endo.append(curr)
                                shift = 1
                                cnt_sdl += 1
                                new_sdl.append(new_stmp[cnt_qoi])
                                new_endo.append(new_etmp[cnt_qoi])
                            else:
                                #the current saddle has a higher height, the saddle that was already in the current position remains
                                cnt_curr += 1
                                new_sdl.append(new_stmp[cnt_qoi])
                                new_endo.append(new_etmp[cnt_qoi])
                        # If shift is null it means that I have not inserted the saddle I have to insert anywhere because evidently
                        # the height is greater than the heights of the saddles already present. I therefore put it last
                        if shift == 0:
                            shift = 1
                            cnt_curr += 1
                            new_sdl.append(cnt_sdl)
                            new_endo.append(curr)
                    
                    
                    if cnt_curr > n_el_new:
                        n_el_new = cnt_curr
                        
                                        
                    if new_sdl == []:
                        new_sdl = [0]
                        new_endo = [curr]
                        new_stmp[:n_el_new] = new_sdl[:n_el_new]
                        new_etmp[:n_el_new] = new_endo[:n_el_new]
                        new_sdl = []
                        new_endo = []
                        
                    else:            
                        new_stmp[:n_el_new] = new_sdl[:n_el_new]
                        new_etmp[:n_el_new] = new_endo[:n_el_new]

                
                cnt_new = n_el_new
                new_sdl[:n_el_new] = new_stmp[:n_el_new]
                new_endo[:n_el_new] = new_etmp[:n_el_new]
                
                cnt_curr = 0
                mem_new = 0
                qoi_sdl = []
                qoi_endo = []
                for cnt_qoi in range(n_el): #loop on the elements of the main vector
                    for cnt_sdl in range(mem_new, cnt_new):
                        zs_cur = model.rd_pt[model.l_endo_pt[new_endo[cnt_sdl]-1].idms.value[new_sdl[cnt_sdl]]-1].Z
                        if zs_cur < model.rd_pt[model.l_endo_pt[qoi_etmp[cnt_qoi]-1].idms.value[qoi_stmp[cnt_qoi]]-1].Z:
                            cnt_curr += 1
                            qoi_sdl.append(new_sdl[cnt_sdl])
                            qoi_endo.append(new_endo[cnt_sdl])
                            mem_new += 1 #I count how many I've already written
                    cnt_curr  += 1
                    qoi_sdl.append(qoi_stmp[cnt_qoi])
                    qoi_endo.append(qoi_etmp[cnt_qoi])
                if mem_new < cnt_new:
                    qoi_sdl.extend(new_sdl[mem_new:cnt_new])
                    qoi_endo.extend(new_endo[mem_new:cnt_new])
                
                n_el = cnt_curr
                qoi_stmp = qoi_sdl[:n_el]
                qoi_etmp = qoi_endo[:n_el]
        
        qoi_sdl[:n_el] = qoi_stmp[:n_el]
        qoi_endo[:n_el] = qoi_etmp[:n_el]
        n_qoi = n_el
    
                
                
                    
                            
                            
    
                    
                
                

                                
 
                        
    
    
    
