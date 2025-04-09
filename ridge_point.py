# -*- coding: utf-8 -*-
"""
The subroutine Ridge Point define the relationship between ridge points
as they are ready to be used to extract an ordered ridge network.
Modify fldir and id_ch
"""


from tqdm import tqdm
import numpy as np

from mutual_dist import md
from data_structures import RidgePoint, DrainagePoint, IDPointer
from thal_net_hso_length import thal_net_hso_length
from dpl_ss import dpl_ss
from a_endo import a_endo
from junc_update import junc_update



def find_ridge_neighbors(model):
    """
    Finds neighboring ridge points for each ridge point in the model using NumPy vectorized operations.
    
    Parameters
    ----------
    model : object
        Model containing:
        - mat_id: NumPy array of RidgePoint objects (or None)
        - rd_pt: List of RidgePoint objects
        - N, M: Number of rows and columns of the DEM.
    """
    
    n_rdpnt = len(model.rd_pt)
    n_rdpnt_old=n_rdpnt
    i_max, j_max = model.mat_id.shape
    i_max, j_max = i_max-1, j_max-1

     # Create a boolean mask where True indicates a RidgePoint
    ridge_mask = np.vectorize(lambda x: isinstance(x, RidgePoint))(model.mat_id)
    drainage_mask = np.vectorize(lambda x: isinstance(x, DrainagePoint))(model.mat_id)
    river_mask = np.vectorize(lambda x: isinstance(x, IDPointer))(model.mat_id)

    associated_mask = ridge_mask + drainage_mask + river_mask

    # Shifted arrays to find direct neighbors
    up = np.roll(ridge_mask, shift=1, axis=0) #Il faut décaler dans le sens opposé
    up[0,:] = False
    down = np.roll(ridge_mask, shift=-1, axis=0) #Il faut décaler dans le sens opposé
    down[-1,:] = False
    left = np.roll(ridge_mask, shift=1, axis=1)  
    left[:,0] = False
    right = np.roll(ridge_mask, shift=-1, axis=1)  
    right[:,-1] = False
    


    # Iterate over ridge points and assign neighbors based on the computed neighbor mask
    print("Find the neighbors")
    for cnt_rdpt in tqdm(range(len(model.rd_pt))):
        rp = model.rd_pt[cnt_rdpt]
        i = rp.i
        j = rp.j
        
        if up[i,j]: #up
            rp.id_neigh.append(model.mat_id[i-1,j].id_pnt)
            rp.nen += 1
        if down[i,j]: #down
            rp.id_neigh.append(model.mat_id[i+1,j].id_pnt)
            rp.nen += 1
        if left[i,j]: #left
            rp.id_neigh.append(model.mat_id[i,j-1].id_pnt)
            rp.nen += 1
        if right[i,j]: #right
            rp.id_neigh.append(model.mat_id[i,j+1].id_pnt)
            rp.nen += 1
        
        
        #bottom left
        if 0 <= i+1 <= i_max and 0 <= j-1 <= j_max:
            if ridge_mask[i+1, j-1] and associated_mask[i+1,j] and associated_mask[i,j-1]:
                if not(ridge_mask[i+1, j] or ridge_mask[i, j-1]):
                    # means that in one of the two corners of the square there is another ridge point for which a diagonal would be cut through the right angle
                    # No need to trace diagonal segment
                    if drainage_mask[i+1, j]:
                        irr1 = i-1
                        icc1 = j-2
                        irr2 = i+1
                        icc2 = j
                    if drainage_mask[i, j-1]:
                        irr1 = i+2
                        icc1 = j+1
                        irr2 = i
                        icc2 = j-1
                    if 0 <= irr1 <= i_max and 0 <= icc1 <= j_max:
                        if associated_mask[irr1, icc1] and associated_mask[irr2, icc2]:
                            if not ridge_mask[irr1, icc1]:
                                if (model.mat_id[irr1, icc1].fldir.value != model.mat_id[irr2, icc2].id_pnt.value 
                                    and model.mat_id[irr2, icc2].fldir.value != model.mat_id[irr1, icc1].id_pnt.value):
                                    rp.id_neigh.append(model.mat_id[i+1,j-1].id_pnt)
                                    rp.nen += 1
        
        #bottom right
        if 0 <= i+1 <= i_max and 0 <= j+1 <= j_max:
            if ridge_mask[i+1, j+1] and associated_mask[i+1,j] and associated_mask[i,j+1]:
                if not(ridge_mask[i+1, j] or ridge_mask[i, j+1]):
                    # No need to trace diagonal segment
                    if drainage_mask[i+1, j]:
                        irr1 = i-1
                        icc1 = j+2
                        irr2 = i+1
                        icc2 = j
                    if drainage_mask[i, j+1]:
                        irr1 = i+2
                        icc1 = j-1
                        irr2 = i
                        icc2 = j+1
                    if 0 <= irr1 <= i_max and 0 <= icc1 <= j_max:
                        if associated_mask[irr1, icc1] and associated_mask[irr2, icc2]:
                            if not ridge_mask[irr1, icc1]:
                                if (model.mat_id[irr1, icc1].fldir.value != model.mat_id[irr2, icc2].id_pnt.value 
                                    and model.mat_id[irr2, icc2].fldir.value != model.mat_id[irr1, icc1].id_pnt.value):
                                    rp.id_neigh.append(model.mat_id[i+1,j+1].id_pnt)
                                    rp.nen += 1
                                
        #top right
        if 0 <= i-1 <= i_max and 0 <= j+1 <= j_max:
            if ridge_mask[i-1, j+1] and associated_mask[i-1,j] and associated_mask[i,j+1]:
                if not(ridge_mask[i-1, j] or ridge_mask[i, j+1]):
                    # No need to trace diagonal segment
                    if drainage_mask[i-1, j]:
                        irr1 = i+1
                        icc1 = j+2
                        irr2 = i-1
                        icc2 = j
                    if drainage_mask[i, j+1]:
                        irr1 = i-2
                        icc1 = j-1
                        irr2 = i
                        icc2 = j+1
                    if 0 <= irr1 <= i_max and 0 <= icc1 <= j_max:
                        if associated_mask[irr1, icc1] and associated_mask[irr2, icc2]:
                            if not ridge_mask[irr1, icc1]:
                                if (model.mat_id[irr1, icc1].fldir.value != model.mat_id[irr2, icc2].id_pnt.value 
                                    and model.mat_id[irr2, icc2].fldir.value != model.mat_id[irr1, icc1].id_pnt.value):
                                    rp.id_neigh.append(model.mat_id[i-1,j+1].id_pnt)
                                    rp.nen += 1


        #top left
        if 0 <= i-1 <= i_max and 0 <= j-1 <= j_max:
            if ridge_mask[i-1, j-1] and associated_mask[i-1,j] and associated_mask[i,j-1]:
                if not(ridge_mask[i-1, j] or ridge_mask[i, j-1]):
                    # No need to trace diagonal segment
                    if drainage_mask[i-1, j]:
                        irr1 = i+1
                        icc1 = j-2
                        irr2 = i-1
                        icc2 = j
                    if drainage_mask[i, j-1]:
                        irr1 = i-2
                        icc1 = j+1
                        irr2 = i
                        icc2 = j-1
                    if 0 <= irr1 <= i_max and 0 <= icc1 <= j_max:
                        if associated_mask[irr1, icc1] and associated_mask[irr2, icc2]:
                            if not ridge_mask[irr1, icc1]:
                                if (model.mat_id[irr1, icc1].fldir.value != model.mat_id[irr2, icc2].id_pnt.value 
                                    and model.mat_id[irr2, icc2].fldir.value != model.mat_id[irr1, icc1].id_pnt.value):
                                    rp.id_neigh.append(model.mat_id[i-1,j-1].id_pnt)
                                    rp.nen += 1

                
            
    print("Update saddle neighbor")
    for cnt_sdl in tqdm(range(len(model.sdl_pt))):
        irs = model.sdl_pt[cnt_sdl].id_rdpt
        if model.rd_pt[irs-1].md == 1e10 and model.rd_pt[irs-1].id_sdl != None:
            if model.sdl_pt[cnt_sdl].id_cis_endo.value != None:
                cnt_nen1 = 0
                tmp_neigh1 = []
                for cnt_el in range(model.rd_pt[irs-1].nen):
                    id_ne = model.rd_pt[irs-1].id_neigh[cnt_el]
                    if model.rd_pt[id_ne-1].md == 1e10 and model.rd_pt[id_ne-1].id_sdl != None:
                        if model.sdl_pt[model.rd_pt[id_ne-1].id_sdl-1].id_cis_endo.value != None:
                            # it means that the neighbor point is a saddle point too
                            # the reference to the neighbor saddle in the current saddle is deleted
                            cnt_nen2 = 0
                            tmp_neigh2 = []
                            # reference to the current saddle in the neighbor saddle deleted
                            for i in range(model.rd_pt[id_ne-1].nen):
                                if model.rd_pt[id_ne-1].id_neigh[i] != irs:
                                    cnt_nen2 += 1
                                    tmp_neigh2.append(model.rd_pt[id_ne-1].id_neigh[i])
                            model.rd_pt[id_ne-1].nen = cnt_nen2
                            model.rd_pt[id_ne-1].id_neigh = tmp_neigh2
                        else:
                            cnt_nen1=cnt_nen1+1
                            tmp_neigh1.append(model.rd_pt[irs-1].id_neigh[cnt_el])
                    else:
                       cnt_nen1 += 1
                       tmp_neigh1.append(model.rd_pt[irs-1].id_neigh[cnt_el])
                model.rd_pt[irs-1].nen = cnt_nen1
                model.rd_pt[irs-1].id_neigh = tmp_neigh1
    
    cnt_sdl_2pt = 0
    
    print("splitting endorheic saddle points in two different coincident points")
    for cnt_sdl in tqdm(range(len(model.sdl_pt))):
        # splitting endorheic saddle points in two different coincident points
        # with one neighbor 
        irs = model.sdl_pt[cnt_sdl].id_rdpt
        rp_irs = model.rd_pt[irs-1]
        if rp_irs.nen > 0:
            if rp_irs.md == 1e10 and model.sdl_pt[cnt_sdl].id_cis_endo.value != None:
                # the if below was added  to account for
                # saddles with more than two adjacent points that need to be
                # simplified by leaving only two points and the excess
                # connecting them to the adjacent ones
                #      o                   o         Or any of the other similar combinations
                #      |                 /           
                #   o--o--o    ----->   o--o--o
                #      |                     /
                #      o                   o
                if rp_irs.nen > 2:
                    # Endorheic saddle with more than 2 neighbor points
                    cnt_sdl_2pt += 1
                    tmp_elab = rp_irs.id_neigh
                    for cnt_el in range(rp_irs.nen):
                        if tmp_elab[cnt_el] != None:
                            id_ne = rp_irs.id_neigh[cnt_el]
                            rr = model.rd_pt[id_ne-1].i
                            cr = model.rd_pt[id_ne-1].j
                            for cnt_el2 in range(rp_irs.nen):#looking for another neighbour where moving the segment
                                if tmp_elab[cnt_el2] != None:
                                    id_ne2 = rp_irs.id_neigh[cnt_el2]
                                    rr2 = model.rd_pt[id_ne2-1].i
                                    cr2 = model.rd_pt[id_ne2-1].j
                                    d=((rr-rr2)**2+(cr-cr2)**2)**0.5
                                    # if the two ridge_id2 are the same it means that the exit channel of the endorheic basin passes through the middle
                                    rr_sdl = rp_irs.i
                                    cr_sdl = rp_irs.j
                                    rr4=rr_sdl-(rr_sdl-rr+rr_sdl-rr2)
                                    cr4=cr_sdl-(cr_sdl-cr+cr_sdl-cr2)
                                    if (d == 2**0.5):
                                        if not(model.sdl_pt[cnt_sdl].id_cis_endo.value == model.mat_id[rr4, cr4].id_pnt.value #channel spilling from saddle 
                                            or model.sdl_pt[cnt_sdl].id_trans_out.value == model.mat_id[rr4, cr4].id_pnt.value):
                                            #no channel between points (rr,cr) and (rr2,cr2) 
                                            tmp_elab[cnt_el2] = None
                                            tmp_elab2 = model.rd_pt[id_ne-1].id_neigh
                                            model.rd_pt[id_ne-1].nen += 1
                                            model.rd_pt[id_ne-1].id_neigh = tmp_elab2
                                            model.rd_pt[id_ne-1].id_neigh.append(id_ne2)
                                            for id2_ne in range(model.rd_pt[id_ne2-1].nen):
                                                if model.rd_pt[id_ne2-1].id_neigh[id2_ne] == irs:
                                                    model.rd_pt[id_ne2-1].id_neigh[id2_ne] = id_ne
                    
                    cnt = 0
                    tmp_elab1 = []
                    for cnt_el3 in range(rp_irs.nen):
                        if tmp_elab[cnt_el3] != None:
                            cnt += 1
                            tmp_elab1.append(rp_irs.id_neigh[cnt_el3])
                    rp_irs.nen = cnt
                    rp_irs.id_neigh = tmp_elab1
                
                #splitting all saddle, those with 2 neighbhour points and those just reduced  
                for cnt_el in range(1,rp_irs.nen):
                    n_rdpnt += 1
                    model.sdl_pt[cnt_sdl].id_rdpt2 = n_rdpnt
                    new_rp = RidgePoint(rp_irs.i, rp_irs.j, rp_irs.Z, rp_irs.id_pnt)
                    new_rp.md = rp_irs.md
                    new_rp.id_drpt1 = rp_irs.id_drpt1
                    new_rp.id_drpt2 = rp_irs.id_drpt2
                    new_rp.nen = 1
                    new_rp.id_neigh = [rp_irs.id_neigh[cnt_el]]
                    new_rp.id_sdl = rp_irs.id_sdl
                    model.rd_pt.append(new_rp)
                rp_irs.nen = 1
                rp.id_neigh = rp.id_neigh[:1]
    
    #id correction
    print("Id correction")
    for cnt_rdpt in tqdm(range(n_rdpnt_old,n_rdpnt)):
        rp = model.rd_pt[cnt_rdpt]
        for cnt_el in range(rp.nen):
            id_ne = rp.id_neigh[cnt_el]
            # id_org = model.sdl_pt[rp.id_sdl-1].id_rdpt2.value
            id_ne2 = rp.id_neigh[cnt_el]
            for cnt_el2 in range(model.rd_pt[id_ne2-1].nen):
                if model.rd_pt[id_ne2-1].id_neigh[cnt_el2] == rp.id_pnt:
                    rp.id_pnt = model.sdl_pt[rp.id_sdl-1].id_rdpt2
                    model.rd_pt[id_ne2-1].id_neigh[cnt_el2] = rp.id_pnt
                    
    if model.main_channel_choice == 2:
        thal_net_hso_length(model, river_mask)
    
    dpl_ss(model)
    a_endo(model)
    junc_update(model)
    
    #Mutual distance update
    print("Update ridge points mutual distance")
    for rp in tqdm(model.rd_pt, desc="Calculating md for ridge points"):
        if model.dr_pt[rp.id_drpt1.value-1].dpl > 0 and model.dr_pt[rp.id_drpt2.value-1].dpl > 0:
            rp.md = md(model.dr_pt[rp.id_drpt1.value-1], model.dr_pt[rp.id_drpt2.value-1], model)
        
        rp.A_in = max(model.dr_pt[rp.id_drpt1.value-1].A_in, model.dr_pt[rp.id_drpt2.value-1].A_in)
        rp.A_in_min = min(model.dr_pt[rp.id_drpt1.value-1].A_in, model.dr_pt[rp.id_drpt2.value-1].A_in)

        
    
                


                            

  

        
        
        
        