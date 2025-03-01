# -*- coding: utf-8 -*-
"""
 THALweg NETwork Horton Stream Order LENGTH method construction

 The subroutine thal_net_hso_length define the thalweg network. If necessary
 the dowslope path that spill out from saddles of endorheic basins is traced.
 At each junction the main channel of the thalweg newtwork is defined as 
 the channel with the maximum upslope length.
"""

from tqdm import tqdm
import numpy as np


from data_structures import IDPointer, ListPointer, DrainageNetwork, DrainagePointInflow


def thal_net_hso_length(model, river_mask):
    
    # Build the list 'qoi' by sorting DrainagePoints in ascending Z
    # qoi will store the 'id_pnt' for each DrainagePoint
    qoi = [dp.id_pnt.value for dp in sorted(model.dr_pt, key=lambda dp: dp.Z)]
    
    #Reinitialize id_ch, ninf, inflow, Linflow
    for dp in model.dr_pt:
        dp.id_ch = IDPointer()
        dp.inflow = ListPointer()
        dp.Linflow = ListPointer()
        
    
    #Reinitialize dr_net
    model.dr_net = []
    # On recrée la liste avec des objets vides possédant les bons attributs
    dr_net = [DrainageNetwork() for _ in range(len(model.dr_pt))]
    
    dr_pt_in = [DrainagePointInflow() for _ in range(len(model.dr_pt))]

    # Reinitialize rivers
    model.mat_id = np.where(river_mask, None, model.mat_id)
    
    inet = 0
    
    for id_dr in tqdm(qoi[::-1], desc="Building drainage networks HSO length", unit="point"):
        # model.dr_pt[id_dr-1] = model.dr_pt[id_dr-1]
        # id_dr = model.dr_pt[id_dr-1].id_pnt #Pointer in this version
        i_curr = model.dr_pt[id_dr-1].i
        j_curr = model.dr_pt[id_dr-1].j
        if model.dr_pt[id_dr-1].ninf == 0 and model.dr_pt[id_dr-1].id_endo.value >= 0:

            if model.dr_pt[id_dr-1].fldir.value != None or model.dr_pt[id_dr-1].fldir_ss.value != None:
                if model.dr_pt[id_dr-1].fldir_ss.value != None:
                    curr_fldir = model.dr_pt[id_dr-1].fldir_ss
                    model.dr_pt[id_dr-1].fldir = model.dr_pt[id_dr-1].fldir_ss
                else:
                    curr_fldir = model.dr_pt[id_dr-1].fldir
                max_Z = model.dr_pt[id_dr-1].Z
                if curr_fldir.value != None:

                    # model.dr_pt[curr_fldir.value-1] = model.dr_pt[curr_fldir.value-1]
                    #drainage network 
                    if model.dr_pt[id_dr-1].A_in > 0 :

                        # channel points
                        # all upslope points are already considered 
                        curr_ch = model.dr_pt[id_dr-1].id_ch.value
                        dr_net[curr_ch-1].nel += 1
                        dr_net[curr_ch-1].sso = dr_net[curr_ch-1].hso
                        dr_net[curr_ch-1].id_pnts.append(curr_fldir.value)
                        dr_net[curr_ch-1].id_end_pt = curr_fldir #id of the end point
                        dr_net[curr_ch-1].length += model.delta_x*((i_curr-model.dr_pt[curr_fldir.value-1].i)**2 + (j_curr-model.dr_pt[curr_fldir.value-1].j)**2)**0.5
                        i_mat = i_curr*2 + (model.dr_pt[curr_fldir.value-1].i - i_curr)
                        j_mat = j_curr*2 + (model.dr_pt[curr_fldir.value-1].j - j_curr)
                        model.mat_id[i_mat, j_mat] = dr_net[curr_ch-1].id_ch
                        
                        if model.dr_pt[curr_fldir.value-1].A_in > 0:
                            # This means that it has already been processed in dr_net_ss and the upl value has been updated
                            # if it has already been processed then dr_pt_in is already allocated
                            dr_pt_in[curr_fldir.value-1].ninf += 1
                            dr_pt_in[curr_fldir.value-1].inflow.append(model.dr_pt[id_dr-1].id_ch.value) #Maintenant on stocke dans inflow les canaux, plus les points...
                            if dr_net[curr_ch-1].sso == dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso:
                                # upgrading order of the channel with greater length
                                if dr_net[curr_ch-1].length > model.dr_pt[curr_fldir.value-1].upl:
                                    #feature point update as the main channel is the current
                                    dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].hso = dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso
                                    dr_net[curr_ch-1].hso += 1
                                    model.dr_pt[curr_fldir.value-1].upl = dr_net[curr_ch-1].length
                                    model.dr_pt[curr_fldir.value-1].id_ch = dr_net[curr_ch-1].id_ch
                                    model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in+1
                                    model.dr_pt[curr_fldir.value-1].ninf -= 1
                                    dr_net[curr_ch-1].id_ch_out = model.dr_pt[id_dr-1].id_ch
                                    if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                        dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                                else:
                                    dr_net[curr_ch-1].hso = dr_net[curr_ch-1].sso
                                    dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].hso += 1
                                    model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in + 1
                                    model.dr_pt[curr_fldir.value-1].ninf -= 1
                                    if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                        dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                            
                            if dr_net[curr_ch-1].sso > dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso:                            
                                # upgrading order of the channel with greater length                  
                                # feature point update as the main channel is the current
                                dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].hso = dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso
                                model.dr_pt[curr_fldir.value-1].upl = dr_net[curr_ch-1].length
                                model.dr_pt[curr_fldir.value-1].id_ch = dr_net[curr_ch-1].id_ch
                                model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in+1
                                model.dr_pt[curr_fldir.value-1].ninf -= 1
                                dr_net[curr_ch-1].id_ch_out = model.dr_pt[id_dr-1].id_ch
                                if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                                    
                            if dr_net[curr_ch-1].sso < dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso:                            
                                model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in+1
                                model.dr_pt[curr_fldir.value-1].ninf -= 1
                                if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                        
                        else:
                            model.dr_pt[curr_fldir.value-1].upl = dr_net[curr_ch-1].length
                            model.dr_pt[curr_fldir.value-1].id_ch = dr_net[curr_ch-1].id_ch
                            model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in+1
                            model.dr_pt[curr_fldir.value-1].ninf -= 1
                            dr_net[curr_ch-1].id_ch_out = model.dr_pt[id_dr-1].id_ch
                            dr_pt_in[curr_fldir.value-1].ninf = 1
                            dr_pt_in[curr_fldir.value-1].inflow.append(model.dr_pt[id_dr-1].id_ch.value)
                            if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                    
                    else:
                        # (dr_pt(id_dr)%A_in == 0)
                        # CHANNEL HEAD
                        # new element (channel) in drainage network
                        inet += 1
                        # dr_net[inet-1] = dr_net[inet-1]
                        if dr_net[inet-1].id_ch.value == 0 or dr_net[inet-1].id_ch.value == None:
                            dr_net[inet-1].id_ch = IDPointer(inet)
                            model.dr_pt[id_dr-1].id_ch = dr_net[inet-1].id_ch
                            dr_net[inet-1].nel = 2
                            dr_net[inet-1].id_pnts = ListPointer()
                            dr_net[inet-1].id_pnts.append(model.dr_pt[id_dr-1].id_pnt.value)
                            dr_net[inet-1].id_start_pt = model.dr_pt[id_dr-1].id_pnt #id of the start point
                            dr_net[inet-1].id_ch_out = dr_net[inet-1].id_ch #first assignement
                            dr_net[inet-1].id_pnts.append(curr_fldir.value)
                            dr_net[inet-1].id_end_pt = curr_fldir
                            dr_net[inet-1].length = model.delta_x*((i_curr-model.dr_pt[curr_fldir.value-1].i)**2 + (j_curr-model.dr_pt[curr_fldir.value-1].j)**2)**0.5
                            dr_net[inet-1].sso = 1
                            dr_net[inet-1].hso = 1
                            i_mat = i_curr*2 + (model.dr_pt[curr_fldir.value-1].i - i_curr)
                            j_mat = j_curr*2 + (model.dr_pt[curr_fldir.value-1].j - j_curr)
                            model.mat_id[i_mat, j_mat] = dr_net[inet-1].id_ch
                            dr_pt_in[id_dr-1].ninf = 1
                            dr_pt_in[id_dr-1].inflow.append(dr_net[inet-1].id_ch.value)
                            if model.dr_pt[curr_fldir.value-1].A_in > 0:
                                dr_pt_in[curr_fldir.value-1].ninf += 1
                                dr_pt_in[curr_fldir.value-1].inflow.append(dr_net[inet-1].id_ch.value)
                                if dr_net[inet-1].sso == dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso:
                                    #upgrading order of the channel with greater length
                                    if dr_net[inet-1].length > model.dr_pt[curr_fldir.value-1].upl:
                                        #feature point update as the main channel is the current
                                        dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].hso = dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso
                                        dr_net[inet-1].hso += 1
                                        model.dr_pt[curr_fldir.value-1].upl = dr_net[inet-1].length
                                        model.dr_pt[curr_fldir.value-1].id_ch = dr_net[inet-1].id_ch
                                        model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in+1
                                        model.dr_pt[curr_fldir.value-1].ninf -= 1
                                        if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                            dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                                    else:
                                        dr_net[inet-1].hso = dr_net[inet-1].sso
                                        if dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].hso == dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso:
                                            dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].hso += 1
                                        model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in+1
                                        model.dr_pt[curr_fldir.value-1].ninf -= 1
                                        if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                            dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                                if dr_net[inet-1].sso > dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso:
                                    # upgrading order of the channel with greater length                  
                                    # feature point update as the main channel is the current
                                    model.dr_pt[curr_fldir.value-1].upl = dr_net[inet-1].length
                                    model.dr_pt[curr_fldir.value-1].id_ch = dr_net[inet-1].id_ch
                                    model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in+1
                                    model.dr_pt[curr_fldir.value-1].ninf -= 1
                                    if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                        dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                                        
                                if dr_net[inet-1].sso < dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso:
                                    model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in+1
                                    model.dr_pt[curr_fldir.value-1].ninf -= 1
                                    if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                        dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                            else:
                                #downslope point pointed by curr_fldir never processed
                                model.dr_pt[curr_fldir.value-1].upl = dr_net[inet-1].length
                                model.dr_pt[curr_fldir.value-1].id_ch = dr_net[inet-1].id_ch
                                model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in+1
                                model.dr_pt[curr_fldir.value-1].ninf -= 1
                                dr_pt_in[curr_fldir.value-1].ninf = 1
                                dr_pt_in[curr_fldir.value-1].inflow.append(dr_net[inet-1].id_ch.value)
                                if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                
                model.dr_pt[id_dr-1].id_endo.value = -1
    
    model.dr_net = dr_net[:inet]

                                
def dwnslp_hso(model, id_dr, dr_net, dr_pt_in, max_Z):
    """
    Recursive subroutine to trace dowslope path that spill out from saddles
    of endorheic basins coherent with Horton stream order

    """
    # model.dr_pt[id_dr-1] = model.dr_pt[id_dr-1]
    i_curr = model.dr_pt[id_dr-1].i
    j_curr = model.dr_pt[id_dr-1].j
    
    if model.dr_pt[id_dr-1].ninf == 0 and model.dr_pt[id_dr-1].id_endo.value >= 0:
        if model.dr_pt[id_dr-1].fldir.value != None or model.dr_pt[id_dr-1].fldir_ss.value != None:
            if model.dr_pt[id_dr-1].fldir_ss.value != None:
                curr_fldir = model.dr_pt[id_dr-1].fldir_ss
                model.dr_pt[id_dr-1].fldir = model.dr_pt[id_dr-1].fldir_ss
            else:
                curr_fldir = model.dr_pt[id_dr-1].fldir
        
            if curr_fldir.value != None:
                # model.dr_pt[curr_fldir.value-1] = model.dr_pt[curr_fldir.value-1]
                #drainage network 
                if model.dr_pt[id_dr-1].A_in > 0 :
                    # channel points
                    # all upslope points are already considered 
                    curr_ch = model.dr_pt[id_dr-1].id_ch.value
                    dr_net[curr_ch-1].nel += 1
                    dr_net[curr_ch-1].sso = dr_net[curr_ch-1].hso
                    dr_net[curr_ch-1].id_pnts.append(curr_fldir.value)
                    dr_net[curr_ch-1].id_end_pt = curr_fldir #id of the end point
                    dr_net[curr_ch-1].length += model.delta_x*((i_curr-model.dr_pt[curr_fldir.value-1].i)**2 + (j_curr-model.dr_pt[curr_fldir.value-1].j)**2)**0.5
                    i_mat = i_curr*2 + (model.dr_pt[curr_fldir.value-1].i - i_curr)
                    j_mat = j_curr*2 + (model.dr_pt[curr_fldir.value-1].j - j_curr)
                    model.mat_id[i_mat, j_mat] = dr_net[curr_ch-1].id_ch
                    
                    if model.dr_pt[curr_fldir.value-1].A_in > 0:
                        # This means that it has already been processed in dr_net_ss and the upl value has been updated
                        dr_pt_in[curr_fldir.value-1].ninf += 1
                        dr_pt_in[curr_fldir.value-1].inflow.append(dr_net[curr_ch-1].id_ch.value) 
                        if dr_net[curr_ch-1].sso == dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso:
                            # upgrading order of the channel with greater length
                            if dr_net[curr_ch-1].length > model.dr_pt[curr_fldir.value-1].upl:
                                #feature point update as the main channel is the current
                                dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].hso = dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso
                                dr_net[curr_ch-1].hso += 1
                                model.dr_pt[curr_fldir.value-1].upl = dr_net[curr_ch-1].length
                                model.dr_pt[curr_fldir.value-1].id_ch = dr_net[curr_ch-1].id_ch
                                model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in+1
                                model.dr_pt[curr_fldir.value-1].ninf -= 1
                                dr_net[curr_ch-1].id_ch_out = model.dr_pt[curr_fldir.value-1].id_ch
                                if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                            else:
                                dr_net[curr_ch-1].hso = dr_net[curr_ch-1].sso
                                dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].hso += 1
                                model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in + 1
                                model.dr_pt[curr_fldir.value-1].ninf -= 1
                                if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                        
                        if dr_net[curr_ch-1].sso > dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso:                            
                            # upgrading order of the channel with greater length                  
                            # feature point update as the main channel is the current
                            dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].hso = dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso
                            model.dr_pt[curr_fldir.value-1].upl = dr_net[curr_ch-1].length
                            model.dr_pt[curr_fldir.value-1].id_ch = dr_net[curr_ch-1].id_ch
                            model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in+1
                            model.dr_pt[curr_fldir.value-1].ninf -= 1
                            print("hello model.dr_pt[id_dr-1].id_ch")
                            dr_net[curr_ch-1].id_ch_out = model.dr_pt[id_dr-1].id_ch
                            if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                                
                        if dr_net[curr_ch-1].sso < dr_net[model.dr_pt[curr_fldir.value-1].id_ch.value-1].sso:                            
                            model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in+1
                            model.dr_pt[curr_fldir.value-1].ninf -= 1
                            if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                                dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                    
                    else:
                        model.dr_pt[curr_fldir.value-1].upl = dr_net[curr_ch-1].length
                        model.dr_pt[curr_fldir.value-1].id_ch = dr_net[curr_ch-1].id_ch
                        model.dr_pt[curr_fldir.value-1].A_in += model.dr_pt[id_dr-1].A_in+1
                        model.dr_pt[curr_fldir.value-1].ninf -= 1
                        dr_net[curr_ch-1].id_ch_out = model.dr_pt[curr_fldir.value-1].id_ch
                        #has not already been processed then alloco dr_pt_in
                        dr_pt_in[curr_fldir.value-1].ninf = 1
                        dr_pt_in[curr_fldir.value-1].inflow.append(dr_net[curr_ch-1].id_ch.value)
                        if model.dr_pt[curr_fldir.value-1].Z >= max_Z:
                            dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
            model.dr_pt[id_dr-1].id_endo.value = -1
        
    return dr_net, dr_pt_in
                                



                                


                                           

                                        
    

                            
                            




                                        
                                    
                                
                    
                
            

