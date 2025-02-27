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
        dp = model.dr_pt[id_dr-1]
        id_dr = dp.id_pnt #Pointer in this version
        i_curr = dp.i
        j_curr = dp.j
        if dp.ninf == 0 and dp.id_endo.value >= 0:

            if dp.fldir.value != None or dp.fldir_ss.value != None:

                if dp.fldir_ss.value != None:
                    curr_fldir = dp.fldir_ss
                    dp.fldir = dp.fldir_ss
                else:
                    curr_fldir = dp.fldir
                max_Z = dp.Z
                if curr_fldir.value != None:
                    print("hello")

                    dp_fldir = model.dr_pt[curr_fldir.value-1]
                    #drainage network 
                    if dp.A_in > 0 :
                        # channel points
                        # all upslope points are already considered 
                        curr_ch = dp.id_ch.value
                        dr_net[curr_ch-1].net += 1
                        dr_net[curr_ch-1].sso = dr_net[curr_ch].hso
                        dr_net[curr_ch-1].append(curr_fldir.value)
                        dr_net[curr_ch-1].id_end_pt = curr_fldir #id of the end point
                        dr_net[curr_ch-1].length += model.delta_x*((i_curr-dp_fldir.i)**2 + (j_curr-dp_fldir.j)**2)**0.5
                        i_mat = i_curr*2 + (dp_fldir.i - i_curr)
                        j_mat = j_curr*2 + (dp.fldir.j - j_curr)
                        model.mat_id[i_mat, j_mat] = dr_net[curr_ch-1].id_ch
                        
                        if dp_fldir.A_in > 0:
                            # This means that it has already been processed in dr_net_ss and the upl value has been updated
                            # if it has already been processed then dr_pt_in is already allocated
                            dr_pt_in[curr_fldir.value-1].ninf += 1
                            dr_pt_in[curr_fldir.value-1].inflow.append(dp.id_ch.value) #Maintenant on stocke dans inflow les canaux, plus les points...
                            if dr_net[curr_ch-1].sso == dr_net[dp_fldir.id_ch.value-1].sso:
                                # upgrading order of the channel with greater length
                                if dr_net[curr_ch-1].length > dp_fldir.upl:
                                    #feature point update as the main channel is the current
                                    dr_net[dp_fldir.id_ch.value-1].hso = dr_net[dp_fldir.id_ch.value-1].sso
                                    dr_net[curr_ch-1].hso += 1
                                    dp_fldir.upl = dr_net[curr_ch-1].length
                                    dp_fldir.id_ch = dr_net[curr_ch-1].id_ch
                                    dp_fldir.A_in += dp.A_in+1
                                    dp_fldir.ninf -= 1
                                    dr_net[curr_ch-1].id_ch_out = dp.id_ch
                                    if dp_fldir.Z >= max_Z:
                                        dr_net = dwnslp_hso(model, curr_fldir.value)
                                else:
                                    dr_net[curr_ch-1].hso = dr_net[curr_ch].sso
                                    dr_net[dp_fldir.id_ch-1].hso += 1
                                    dp_fldir.A_in += dp.A_in + 1
                                    dp_fldir.nin -= 1
                                    if dp_fldir >= max_Z:
                                        dr_net = dwnslp_hso(model, curr_fldir.value)
                            
                            if dr_net[curr_ch-1].sso > dr_net[dp_fldir.id_ch.value-1].sso:                            
                                # upgrading order of the channel with greater length                  
                                # feature point update as the main channel is the current
                                dr_net[dp_fldir.id_ch.value-1].hso = dr_net[dp_fldir.id_ch.value-1].sso
                                dp_fldir.upl = dr_net[curr_ch-1].length
                                dp_fldir.id_ch = dr_net[curr_ch-1].id_ch
                                dp_fldir.A_in += dp.A_in+1
                                dp_fldir.ninf -= 1
                                dr_net[curr_ch-1].id_ch_out = dp.id_ch
                                if dp_fldir.Z >= max_Z:
                                    dr_net = dwnslp_hso(model, curr_fldir.value)
                                    
                            if dr_net[curr_ch-1].sso < dr_net[dp_fldir.id_ch.value-1].sso:                            
                                dp_fldir.A_in += dp.A_in+1
                                dp_fldir.ninf -= 1
                                dr_net[curr_ch-1].id_ch_out = dp.id_ch
                                if dp_fldir.Z >= max_Z:
                                    dr_net = dwnslp_hso(model, curr_fldir.value)
                        
                        else:
                            dp_fldir.upl = dr_net[curr_ch-1].length
                            dp_fldir.id_ch = dr_net[curr_ch-1].id_ch
                            dp_fldir.A_in += dp.A_in+1
                            dp_fldir.ninf -= 1
                            dr_net[curr_ch-1].id_ch_out = dp.id_ch
                            dr_pt_in[curr_fldir.value-1].ninf = 1
                            dr_pt_in[curr_fldir.value-1].inflow.append(dp.id_ch.value)
                            if dp_fldir.Z >= max_Z:
                                dr_net = dwnslp_hso(model, curr_fldir.value)
                    
                    else:
                        # (dr_pt(id_dr)%A_in == 0)
                        # CHANNEL HEAD
                        # new element (channel) in drainage network
                        inet += 1
                        dr_net_inet = dr_net[inet-1]
                        if dr_net_inet.id_ch.value == 0:
                            dr_net_inet.id_ch = IDPointer(inet)
                            dp.id_ch = dr_net_inet.id_ch
                            dr_net_inet.nel = 2
                            dr_net_inet.id_pnts = ListPointer()
                            dr_net_inet.id_pnts.append(dp.id_pnt.value)
                            dr_net.id_start_pt = dp.id_pnt #id of the start point
                            dr_net_inet.id_ch_out = dr_net_inet.id_ch #first assignement
                            dr_net_inet.id_pnts.append(curr_fldir.value)
                            dr_net_inet.id_end_pt = curr_fldir
                            dr_net_inet.length = model.delta_x*((i_curr-dp_fldir.i)**2 + (j_curr-dp_fldir.j)**2)**0.5
                            dr_net_inet.sso = 1
                            dr_net_inet.hso = 1
                            i_mat = i_curr*2 + (dp_fldir.i - i_curr)
                            j_mat = j_curr*2 + (dp.fldir.j - j_curr)
                            model.mat_id[i_mat, j_mat] = dr_net[inet-1].id_ch
                            dr_pt_in[id_dr-1].ninf = 1
                            dr_pt_in[id_dr-1].inflow.append(dr_net_inet.id_ch.value)
                            
                            if dp_fldir.A_in > 0:
                                dr_pt_in[curr_fldir.value-1].ninf += 1
                                dr_pt_in[curr_fldir.value-1].inflow.append(dr_net_inet.id_ch.value)
                                if dr_net_inet.sso == dr_net[dp_fldir.id_ch.value-1].sso:
                                    #upgrading order of the channel with greater length
                                    if dr_net_inet.length > dp_fldir.upl:
                                        #feature point update as the main channel is the current
                                        dr_net[dp_fldir.id_ch.value-1].hso = dr_net[dp_fldir.id_ch.value-1].sso
                                        dr_net_inet.hso += 1
                                        dp_fldir.upl = dr_net_inet.length
                                        dp_fldir.id_ch = dr_net_inet.id_ch
                                        dp_fldir.A_in += dp.A_in+1
                                        dp_fldir.ninf -= 1
                                        if dp_fldir.Z >= max_Z:
                                            dr_net = dwnslp_hso(model, curr_fldir.value)
                                    else:
                                        dr_net_inet.hso = dr_net_inet.sso
                                        if dr_net[dp_fldir.id_ch.value-1].hso == dr_net[dp_fldir.id_ch.value-1].sso:
                                            dr_net[dp_fldir.id_ch.value-1].hso += 1
                                        dp_fldir.A_in += dp.A_in+1
                                        dp_fldir.ninf -= 1
                                        if dp_fldir.Z >= max_Z:
                                            dr_net = dwnslp_hso(model, curr_fldir.value)
                                if dr_net_inet.sso > dr_net[dp_fldir.id_ch.value-1].sso:
                                    # upgrading order of the channel with greater length                  
                                    # feature point update as the main channel is the current
                                    dp_fldir.upl = dr_net_inet.length
                                    dp_fldir.id_ch = dr_net_inet.id_ch
                                    dp_fldir.A_in += dp.A_in+1
                                    dp_fldir.ninf -= 1
                                    if dp_fldir.Z >= max_Z:
                                        dr_net = dwnslp_hso(model, curr_fldir.value)
                                        
                                if dr_net_inet.sso < dr_net[dp_fldir.id_ch.value-1].sso:
                                    dp_fldir.A_in += dp.A_in+1
                                    dp_fldir.ninf -= 1
                                    if dp_fldir.Z >= max_Z:
                                        dr_net = dwnslp_hso(model, curr_fldir.value)
                            else:
                                #downslope point pointed by curr_fldir never processed
                                dp_fldir.upl = dr_net_inet.length
                                dp_fldir.id_ch = dr_net_inet.id_ch
                                dp_fldir.A_in += dp.A_in+1
                                dp_fldir.ninf -= 1
                                dr_pt_in[curr_fldir.value-1].ninf = 1
                                dr_pt_in[curr_fldir.value-1].inflow.append(dr_net_inet.id_ch.value)
                                if dp_fldir.Z >= max_Z:
                                    dr_net = dwnslp_hso(model, curr_fldir.value)
                
                dp.id_endo.value = IDPointer(-1)
    
    model.dr_net = dr_net[:inet]

                                
def dwnslp_hso(model, id_dr):
    pass
                                



                                


                                           

                                        
    

                            
                            




                                        
                                    
                                
                    
                
            

