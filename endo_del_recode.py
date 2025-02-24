# -*- coding: utf-8 -*-
"""
ENDOrheic basins DELineation
"""

# endo_dl.py
from tqdm import tqdm

from data_structures import SaddlePoint


def endo_del(model):
    """
    Delineates endorheic basins by processing ridge points whose mutual distance is equal to a very high value.
    For each ridge point (rd_pt) whose md attribute equals a large value (1e10),
    the function calls find_saddle to find or create saddle points.
    """

    # Initialize the list of saddle points.
    model.sdl_pt = []  # List of SaddlePoint objects.
    cnt_sdl = 0       # Counter for saddle points.

    for i in tqdm(range(len(model.rd_pt))):
        rp = model.rd_pt[i]
        if rp.md == 1e10:
            # Retrieve id_eo1 and id_eo2:
            # For id_eo1: access the DrainagePoint corresponding to rd.id_drpt1, then get its channel, then endo id.
            id_eo1 = model.l_dr_net[model.dr_pt[rp.id_drpt1.value-1].id_ch.value-1].id_endo.value
            id_eo2 = model.l_dr_net[model.dr_pt[rp.id_drpt2.value-1].id_ch.value-1].id_endo.value
            sdl_mem = 1
            cnt_sdl = find_saddle(rp.id_pnt, id_eo1, id_eo2, sdl_mem, model, cnt_sdl)
            sdl_mem = 0
            cnt_sdl = find_saddle(rp.id_pnt, id_eo2, id_eo1, sdl_mem, model, cnt_sdl)
    
    for sd_pt in model.sdl_pt:
        sd_pt.i = model.rd_pt[sd_pt.id_rdpt-1].i
        sd_pt.j = model.rd_pt[sd_pt.id_rdpt-1].j
        sd_pt.Z = model.rd_pt[sd_pt.id_rdpt-1].Z


            
def find_saddle(cnt_rdpt, id_cis, id_trans, sdl_mem, model, cnt_sdl):
    """FIND minimum SADDLE point for each ridge separating two endorheic basins.
    Searches (or creates) a minimal saddle point for the ridge point whose identifier
    is cnt_rdpt. For the provided "cis" (id_cis) and "trans" (id_trans) sides,
    the function updates the corresponding endo_pt structure.
    
    Parameters
    ----------
    cnt_rdpt : int
        Idof the ridge point.
    id_cis : int
        Identifier (endos) of the basin on the "cis" side.
    id_trans : int
        Identifier (endos) of the basin on the "trans" side.
    sdl_mem : int
        Flag indicating whether to create a new saddle point (1) or just update.
    model : object
        The hydrological model containing the necessary lists.
    cnt_sdl : int
        Current counter of saddle points.
        
    Returns
    -------
    cnt_sdl : int
        The updated saddle point counter after processing.
    """
    ep_cis = model.l_endo_pt[id_cis-1]
    if ep_cis.nsaddle == 0:
        ep_cis.nsaddle = 1
        ep_cis.beyo_sad.append(id_trans) #id endo other side of the saddle
        ep_cis.idms.append(cnt_rdpt) #id of rd_pt that is also the min saddle
        
        #New saddle point
        if sdl_mem == 1:
            cnt_sdl += 1
            # Create a new saddle point.
            new_sdl = SaddlePoint(id_pnt=cnt_sdl)
            new_sdl.id_rdpt = model.rd_pt[cnt_rdpt-1].id_pnt
            model.rd_pt[cnt_rdpt-1].id_sdl = new_sdl.id_pnt
            model.sdl_pt.append(new_sdl)

    else:
        flag = 0
        for cnt_en_sdl in range(model.l_endo_pt[id_cis-1].nsaddle):
            if ep_cis.beyo_sad.value[cnt_en_sdl-1] == id_trans:
                #already stored a saddle point between the two current endorheic basins
                flag = 1
                #It verify if the current is lower then the stored
                if model.rd_pt[cnt_rdpt-1].Z < model.rd_pt[ep_cis.idms.value[cnt_en_sdl-1]-1].Z:
                    if sdl_mem == 1:
                        model.rd_pt[cnt_rdpt-1].id_sdl = model.rd_pt[ep_cis.idms.value[cnt_en_sdl-1]-1].id_sdl
                        # print("\nmodel.rd_pt[ep_cis.idms.value[cnt_en_sdl-1]] :", model.rd_pt[ep_cis.idms.value[cnt_en_sdl-1]-1])
                        model.rd_pt[ep_cis.idms.value[cnt_en_sdl-1]-1].id_sdl = None
                        ep_cis.idms.value[cnt_en_sdl-1] = cnt_rdpt
                        #saddle point
                        curr_sdl = model.rd_pt[cnt_rdpt-1].id_sdl
                        # print("\ncnt_rdpt :", cnt_rdpt)

                        # print("\ncurr_sdl :", curr_sdl)
                        model.sdl_pt[curr_sdl-1].id_rdpt = model.rd_pt[cnt_rdpt-1].id_pnt
                        flag = 1
                    else:
                        ep_cis.idms.value[cnt_en_sdl-1] = cnt_rdpt
                        flag = 1
        
        if flag == 0:
            # No saddle point between the two current endorheic basins already stored
            # Add a new saddle point
            ep_cis.nsaddle += 1
            ep_cis.beyo_sad.append(id_trans)
            ep_cis.idms.append(cnt_rdpt)
            
            #New saddle point
            if sdl_mem == 1:
                cnt_sdl += 1
                new_sdl = SaddlePoint(id_pnt=cnt_sdl)
                new_sdl.id_rdpt=model.rd_pt[cnt_rdpt-1].id_pnt
                model.sdl_pt.append(new_sdl)
                model.rd_pt[cnt_rdpt-1].id_sdl = new_sdl.id_pnt
    return cnt_sdl

        
