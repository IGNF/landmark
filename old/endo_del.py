# -*- coding: utf-8 -*-
"""
ENDOrheic basins DELineation
"""

# endo_dl.py
import time
from tqdm import tqdm

from data_structures import EndoPoint, SaddlePoint

def endo_del(model):
    """
    Delineates endorheic basins by processing ridge points whose mutual distance is equal to a very high value.
    For each ridge point (rd_pt) whose md attribute equals a large value (e.g., 10e10),
    the function calls find_saddle to find or create saddle points.
    
    Updates model.endo_pt and model.sdl_pt and displays the processing time.
    """
    start_time = time.process_time()
    
    # Initialize the list of saddle points.
    model.sdl_pt = []  # List of SaddlePoint objects.
    cnt_sdl = 0       # Counter for saddle points.
    
    # Process all ridge points. (model.rd_pt is a list of RidgePoint)
    # Assume model.n_rdpnt is the total number of ridge points.
    for cnt_rdpt in tqdm(range(len(model.rd_pt))):  # Using a 0-indexed loop.
        rd = model.rd_pt[cnt_rdpt]  # rd is a RidgePoint object
        # If the mutual distance md is a very high value (e.g., 10e10)
        if rd.md == 10e10:
            # Retrieve id_eo1 and id_eo2:
            # For id_eo1: access the DrainagePoint corresponding to rd.id_drpt1, then get its channel, then endo id.
            dp1 = model.dr_pt_by_id.get(rd.id_drpt1)
            dp2 = model.dr_pt_by_id.get(rd.id_drpt2)
            if dp1 is None or dp2 is None:
                continue
            # Access networks via dp1 and dp2's channels.
            id_eo1 = model.dr_net_by_id.get(dp1.id_ch).id_endo if dp1.id_ch in model.dr_net_by_id else None
            id_eo2 = model.dr_net_by_id.get(dp2.id_ch).id_endo if dp2.id_ch in model.dr_net_by_id else None
            if id_eo1 is None or id_eo2 is None:
                continue
            sdl_mem = 1
            cnt_sdl = find_saddle(cnt_rdpt, id_eo1, id_eo2, sdl_mem, model, cnt_sdl)
            sdl_mem = 0
            cnt_sdl = find_saddle(cnt_rdpt, id_eo2, id_eo1, sdl_mem, model, cnt_sdl)
        
    
    model.n_sdlpt = cnt_sdl
    print(f"New number of saddle points: {len(model.sdl_pt)}")
    
    finish_time = time.process_time()
    elapsed_time = finish_time - start_time
    sph = 3600.0
    ms = 60.0
    hours = int(elapsed_time // sph)
    minutes = int((elapsed_time % sph) // ms)
    seconds = elapsed_time % ms
    print(f"Elapsed time (endo_del): {hours}h {minutes}m {seconds:.2f}s")

def find_saddle(cnt_rdpt, id_cis, id_trans, sdl_mem, model, cnt_sdl):
    """
    Searches (or creates) a minimal saddle point for the ridge point whose identifier
    is rd_pt[cnt_rdpt]. For the provided "cis" (id_cis) and "trans" (id_trans) sides,
    the function updates the corresponding endo_pt structure.
    
    Parameters
    ----------
    cnt_rdpt : int
        Index (0-indexed) of the ridge point in model.rd_pt.
    id_cis : int
        Identifier (endos) of the basin on the "cis" side.
    id_trans : int
        Identifier (endos) of the basin on the "trans" side.
    sdl_mem : int
        Flag indicating whether to create a new saddle point (1) or just update.
    model : object
        The hydrological model containing the necessary lists and dictionaries.
    cnt_sdl : int
        Current counter of saddle points.
        
    Returns
    -------
    cnt_sdl : int
        The updated saddle point counter after processing.
    """
    # Access the EndoPoint object corresponding to id_cis.
    # Check if the EndoPoint corresponding to id_cis is in model.endo_pt.
    endo_index = id_cis - 1  # Convert from 1-based to 0-based indexing
    
    if endo_index < 0:
        # If id_cis is less than 1, this is abnormal -> exit.
        return cnt_sdl

    # If endo_index is beyond the existing endo_pt list, create missing elements.
    while endo_index >= len(model.endo_pt):
        # Create a new EndoPoint with a new identifier endo_eo_val
        new_id_eo = len(model.endo_pt) + 1
        new_endo = EndoPoint(
            id_eo=new_id_eo,
            id_pnt=0,    # Use 0 or an nonexistent drainage point id
            bas_type=0,  # Default value
            nsaddle=0,
            beyo_sad=[],
            idms=[]
        )
        model.endo_pt.append(new_endo)

    # Now endo_index is necessarily valid in the list.
    endo_obj = model.endo_pt[endo_index]
    
    if endo_obj.nsaddle == 0:
        # No saddle has been stored for this basin yet.
        endo_obj.beyo_sad = [id_trans]  # Store id_trans as the other side of the saddle.
        endo_obj.idms = [cnt_rdpt + 1]   # Store cnt_rdpt+1 as an identifier (to match 1-based indexing)
        endo_obj.nsaddle = 1
        if sdl_mem == 1:
            cnt_sdl += 1
            # Create a new saddle point.
            new_sdl = SaddlePoint(id_pnt=cnt_sdl, id_rdpt=model.rd_pt[cnt_rdpt].id_pnt)
            model.sdl_pt.append(new_sdl)
            # Update the ridge point: associate its id_sdl with the new saddle.
            model.rd_pt[cnt_rdpt].id_sdl = new_sdl.id_pnt
        return cnt_sdl
    else:
        flag = 0
        # Traverse the saddle points already stored for this basin.
        for cnt_en_sdl in range(endo_obj.nsaddle):
            if endo_obj.beyo_sad[cnt_en_sdl] == id_trans:
                flag = 1
                # Check if the current ridge point is lower than the already stored one.
                stored_rdpt_id = endo_obj.idms[cnt_en_sdl]
                stored_rd = model.rd_pt[stored_rdpt_id - 1]  # Convert to 0-indexed
                current_rd = model.rd_pt[cnt_rdpt]
                if current_rd.Z < stored_rd.Z:
                    if sdl_mem == 1:
                        current_rd.id_sdl = stored_rd.id_sdl
                        stored_rd.id_sdl = None
                        endo_obj.idms[cnt_en_sdl] = cnt_rdpt + 1
                    else:
                        endo_obj.idms[cnt_en_sdl] = cnt_rdpt + 1
                break
        if flag == 0:
            endo_obj.nsaddle += 1
            endo_obj.beyo_sad.append(id_trans)
            endo_obj.idms.append(cnt_rdpt + 1)
            if sdl_mem == 1:
                cnt_sdl += 1
                new_sdl = SaddlePoint(id_pnt=cnt_sdl, id_rdpt=model.rd_pt[cnt_rdpt].id_pnt)
                model.sdl_pt.append(new_sdl)
                model.rd_pt[cnt_rdpt].id_sdl = new_sdl.id_pnt
        return cnt_sdl
