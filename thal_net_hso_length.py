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
    #!!!!!!!!Ces 3 étapes sont très très longues
    #Reinitialize id_ch, inflow, Linflow
    print("Reset id_ch, inflow, Linflow")
    for dp in tqdm(model.dr_pt):
        dp.id_ch = IDPointer()
        dp.inflow = ListPointer()
        dp.Linflow = ListPointer()
        
    
    #Reinitialize dr_net
    model.dr_net = []
    # print("Drainage network reset")
    # On recrée la liste avec des objets vides possédant les bons attributs
    dr_net = [DrainageNetwork() for _ in tqdm(range(len(model.dr_pt)), desc="Drainage network reset")]
    # dr_net = [DrainageNetwork() for _ in range(len(model.dr_net)*2)]

    
    # print("Temporary drainage points")
    dr_pt_in = [DrainagePointInflow() for _ in tqdm(range(len(model.dr_pt)), desc="Temporary drainage points")]

    # Reinitialize rivers
    model.mat_id = np.where(river_mask, None, model.mat_id)
    
    inet = 0
    
    # print("\n----avant la boucle thal_net_hso_length------ ")
    # print("Valeur de ninf pour le point 24020 : ", model.dr_pt[24019].ninf)
    # print("Valeur de id_endo pour le point 24020 : ", model.dr_pt[24019].id_endo.value)
    # print("Valeur de fldir pour le point 24020 : ", model.dr_pt[24019].fldir.value)
    # print("Valeur de fldir_ss pour le point 24020 : ", model.dr_pt[24019].fldir_ss.value)

    
    for id_dr in tqdm(qoi[::-1], desc="Building drainage networks HSO length", unit="point"):
        dp = model.dr_pt[id_dr-1]
        # id_dr = dp.id_pnt #Pointer in this version
        i_curr = dp.i
        j_curr = dp.j
        # if dp.id_pnt.value == 24020:
        #     print("\n-----dans la boucle thal_net_hso_length------ ")
        #     print("Valeur de ninf pour le point 24020 : ", dp.ninf)
        #     print("Valeur de id_endo pour le point 24020 : ", dp.id_endo.value)
        #     print("Valeur de fldir pour le point 24020 : ", dp.fldir.value)
        #     print("Valeur de fldir_ss pour le point 24020 : ", dp.fldir_ss.value)


        if dp.ninf == 0 and dp.id_endo.value >= 0:
            if dp.fldir.value != None or dp.fldir_ss.value != None:
                if dp.fldir_ss.value != None:
                    curr_fldir = dp.fldir_ss
                    dp.fldir = dp.fldir_ss
                    
                else:
                    curr_fldir = dp.fldir
                max_Z = dp.Z
                if curr_fldir.value != None:
                    # if curr_fldir.value == 24020 :
                    #     print("\n-------curr_fldir = 24020-----------")
                    #     print("Traitement du point : ", dp.id_pnt.value)
                    #     print("ninf de 24020 à l'entrée", model.dr_pt[curr_fldir.value-1].ninf)
                    dp_fldir = model.dr_pt[curr_fldir.value-1]
                    #drainage network 
                    if dp.A_in > 0 :

                        # channel points
                        # all upslope points are already considered 
                        curr_ch = dp.id_ch
                        dr_net[curr_ch.value-1].nel += 1
                        dr_net[curr_ch.value-1].sso = dr_net[curr_ch.value-1].hso
                        dr_net[curr_ch.value-1].id_pnts.append(curr_fldir.value)
                        dr_net[curr_ch.value-1].id_end_pt = curr_fldir #id of the end point
                        dr_net[curr_ch.value-1].length += model.delta_x*((i_curr-dp_fldir.i)**2 + (j_curr-dp_fldir.j)**2)**0.5
                        i_mat = i_curr*2 + (dp_fldir.i - i_curr)
                        j_mat = j_curr*2 + (dp_fldir.j - j_curr)
                        model.mat_id[i_mat, j_mat] = dr_net[curr_ch.value-1].id_ch
                        
                        if dp_fldir.A_in > 0:
                            # This means that it has already been processed in dr_net_ss and the upl value has been updated
                            # if it has already been processed then dr_pt_in is already allocated
                            dr_pt_in[curr_fldir.value-1].ninf += 1
                            dr_pt_in[curr_fldir.value-1].inflow.append(dp.id_ch.value) #Maintenant on stocke dans inflow les canaux, plus les points...
                            if dr_net[curr_ch.value-1].sso == dr_net[dp_fldir.id_ch.value-1].sso:
                                # upgrading order of the channel with greater length
                                if dr_net[curr_ch.value-1].length > dp_fldir.upl:
                                    #feature point update as the main channel is the current
                                    dr_net[dp_fldir.id_ch.value-1].hso = dr_net[dp_fldir.id_ch.value-1].sso
                                    dr_net[curr_ch.value-1].hso += 1
                                    dp_fldir.upl = dr_net[curr_ch.value-1].length
                                    dp_fldir.id_ch = dr_net[curr_ch.value-1].id_ch
                                    dp_fldir.A_in += dp.A_in+1
                                    dp_fldir.ninf -= 1
                                    dr_net[curr_ch.value-1].id_ch_out = dp.id_ch
                                    if dp_fldir.Z >= max_Z:
                                        dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                                else:
                                    dr_net[curr_ch.value-1].hso = dr_net[curr_ch.value-1].sso
                                    dr_net[dp_fldir.id_ch.value-1].hso += 1
                                    dp_fldir.A_in += dp.A_in + 1
                                    dp_fldir.ninf -= 1
                                    if dp_fldir.Z >= max_Z:
                                        dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                            
                            if dr_net[curr_ch.value-1].sso > dr_net[dp_fldir.id_ch.value-1].sso:                            
                                # upgrading order of the channel with greater length                  
                                # feature point update as the main channel is the current
                                dr_net[dp_fldir.id_ch.value-1].hso = dr_net[dp_fldir.id_ch.value-1].sso
                                dp_fldir.upl = dr_net[curr_ch.value-1].length
                                dp_fldir.id_ch = dr_net[curr_ch.value-1].id_ch
                                dp_fldir.A_in += dp.A_in+1
                                dp_fldir.ninf -= 1
                                dr_net[curr_ch.value-1].id_ch_out = dp.id_ch
                                if dp_fldir.Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                                    
                            if dr_net[curr_ch.value-1].sso < dr_net[dp_fldir.id_ch.value-1].sso:                            
                                dp_fldir.A_in += dp.A_in+1
                                dp_fldir.ninf -= 1
                                if dp_fldir.Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                        
                        else:
                            dp_fldir.upl = dr_net[curr_ch.value-1].length
                            dp_fldir.id_ch = dr_net[curr_ch.value-1].id_ch
                            dp_fldir.A_in += dp.A_in+1
                            dp_fldir.ninf -= 1
                            dr_net[curr_ch.value-1].id_ch_out = dp.id_ch
                            dr_pt_in[curr_fldir.value-1].ninf = 1
                            dr_pt_in[curr_fldir.value-1].inflow.append(dp.id_ch.value)
                            if dp_fldir.Z >= max_Z:
                                dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                    
                    else:
                        # (dr_pt(id_dr)%A_in == 0)
                        # CHANNEL HEAD
                        # new element (channel) in drainage network
                        inet += 1
                        dr_net_inet = dr_net[inet-1]
                        if dr_net_inet.id_ch.value == 0 or dr_net_inet.id_ch.value == None:
                            dr_net_inet.id_ch = IDPointer(inet)
                            dp.id_ch = dr_net_inet.id_ch
                            dr_net_inet.nel = 2
                            dr_net_inet.id_pnts = ListPointer()
                            dr_net_inet.id_pnts.append(dp.id_pnt.value)
                            dr_net_inet.id_start_pt = dp.id_pnt #id of the start point
                            dr_net_inet.id_ch_out = dr_net_inet.id_ch #first assignement
                            dr_net_inet.id_pnts.append(curr_fldir.value)
                            dr_net_inet.id_end_pt = curr_fldir
                            dr_net_inet.length = model.delta_x*((i_curr-dp_fldir.i)**2 + (j_curr-dp_fldir.j)**2)**0.5
                            dr_net_inet.sso = 1
                            dr_net_inet.hso = 1
                            i_mat = i_curr*2 + (dp_fldir.i - i_curr)
                            j_mat = j_curr*2 + (dp_fldir.j - j_curr)
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
                                            dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                                    else:
                                        dr_net_inet.hso = dr_net_inet.sso
                                        if dr_net[dp_fldir.id_ch.value-1].hso == dr_net[dp_fldir.id_ch.value-1].sso:
                                            dr_net[dp_fldir.id_ch.value-1].hso += 1
                                        dp_fldir.A_in += dp.A_in+1
                                        dp_fldir.ninf -= 1
                                        if dp_fldir.Z >= max_Z:
                                            dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                                if dr_net_inet.sso > dr_net[dp_fldir.id_ch.value-1].sso:
                                    # upgrading order of the channel with greater length                  
                                    # feature point update as the main channel is the current
                                    dp_fldir.upl = dr_net_inet.length
                                    dp_fldir.id_ch = dr_net_inet.id_ch
                                    dp_fldir.A_in += dp.A_in+1
                                    dp_fldir.ninf -= 1
                                    if dp_fldir.Z >= max_Z:
                                        dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                                        
                                if dr_net_inet.sso < dr_net[dp_fldir.id_ch.value-1].sso:
                                    dp_fldir.A_in += dp.A_in+1
                                    dp_fldir.ninf -= 1
                                    if dp_fldir.Z >= max_Z:
                                        dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                            else:
                                #downslope point pointed by curr_fldir never processed
                                dp_fldir.upl = dr_net_inet.length
                                dp_fldir.id_ch = dr_net_inet.id_ch
                                dp_fldir.A_in += dp.A_in+1
                                dp_fldir.ninf -= 1
                                dr_pt_in[curr_fldir.value-1].ninf = 1
                                dr_pt_in[curr_fldir.value-1].inflow.append(dr_net_inet.id_ch.value)
                                if dp_fldir.Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                # if curr_fldir.value == 24020 :
                #     print("ninf de 24020 à la sortie", model.dr_pt[curr_fldir.value-1].ninf)

                
                dp.id_endo.value = -1
    
    model.dr_net = dr_net[:inet]
    model.dr_pt_in = dr_pt_in
    
    #update ch_out
    for curr_net in model.dr_net:
        curr_net.id_ch_out = model.dr_pt[curr_net.id_end_pt.value-1].id_ch


                                
def dwnslp_hso(model, id_dr, dr_net, dr_pt_in, max_Z):
    """
    Recursive subroutine to trace dowslope path that spill out from saddles
    of endorheic basins coherent with Horton stream order

    """
    dp = model.dr_pt[id_dr-1]
    i_curr = dp.i
    j_curr = dp.j
    
    if dp.ninf == 0 and dp.id_endo.value >= 0:
        if dp.fldir.value != None or dp.fldir_ss.value != None:
            if dp.fldir_ss.value != None:
                curr_fldir = dp.fldir_ss
                dp.fldir = dp.fldir_ss
            else:
                curr_fldir = dp.fldir
        
            if curr_fldir.value != None:
                # if curr_fldir == 24020 :
                #     print("\n-------curr_fldir = 24020-----------")

                dp_fldir = model.dr_pt[curr_fldir.value-1]
                #drainage network 
                if dp.A_in > 0 :
                    # channel points
                    # all upslope points are already considered 
                    curr_ch = dp.id_ch
                    dr_net[curr_ch.value-1].nel += 1
                    dr_net[curr_ch.value-1].sso = dr_net[curr_ch.value-1].hso
                    dr_net[curr_ch.value-1].id_pnts.append(curr_fldir.value)
                    dr_net[curr_ch.value-1].id_end_pt = curr_fldir #id of the end point
                    dr_net[curr_ch.value-1].length += model.delta_x*((i_curr-dp_fldir.i)**2 + (j_curr-dp_fldir.j)**2)**0.5
                    i_mat = i_curr*2 + (dp_fldir.i - i_curr)
                    j_mat = j_curr*2 + (dp_fldir.j - j_curr)
                    model.mat_id[i_mat, j_mat] = dr_net[curr_ch.value-1].id_ch
                    
                    if dp_fldir.A_in > 0:
                        # This means that it has already been processed in dr_net_ss and the upl value has been updated
                        dr_pt_in[curr_fldir.value-1].ninf += 1
                        dr_pt_in[curr_fldir.value-1].inflow.append(dr_net[curr_ch.value-1].id_ch.value) 
                        if dr_net[curr_ch.value-1].sso == dr_net[dp_fldir.id_ch.value-1].sso:
                            # upgrading order of the channel with greater length
                            if dr_net[curr_ch.value-1].length > dp_fldir.upl:
                                #feature point update as the main channel is the current
                                dr_net[dp_fldir.id_ch.value-1].hso = dr_net[dp_fldir.id_ch.value-1].sso
                                dr_net[curr_ch.value-1].hso += 1
                                dp_fldir.upl = dr_net[curr_ch.value-1].length
                                dp_fldir.id_ch = dr_net[curr_ch.value-1].id_ch
                                dp_fldir.A_in += dp.A_in+1
                                dp_fldir.ninf -= 1
                                dr_net[curr_ch.value-1].id_ch_out = dp_fldir.id_ch
                                if dp_fldir.Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                            else:
                                dr_net[curr_ch.value-1].hso = dr_net[curr_ch.value-1].sso
                                dr_net[dp_fldir.id_ch.value-1].hso += 1
                                dp_fldir.A_in += dp.A_in + 1
                                dp_fldir.ninf -= 1
                                if dp_fldir.Z >= max_Z:
                                    dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                        
                        if dr_net[curr_ch.value-1].sso > dr_net[dp_fldir.id_ch.value-1].sso:                            
                            # upgrading order of the channel with greater length                  
                            # feature point update as the main channel is the current
                            dr_net[dp_fldir.id_ch.value-1].hso = dr_net[dp_fldir.id_ch.value-1].sso
                            dp_fldir.upl = dr_net[curr_ch.value-1].length
                            dp_fldir.id_ch = dr_net[curr_ch.value-1].id_ch
                            dp_fldir.A_in += dp.A_in+1
                            dp_fldir.ninf -= 1
                            dr_net[curr_ch.value-1].id_ch_out = dp.id_ch
                            if dp_fldir.Z >= max_Z:
                                dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                                
                        if dr_net[curr_ch.value-1].sso < dr_net[dp_fldir.id_ch.value-1].sso:                            
                            dp_fldir.A_in += dp.A_in+1
                            dp_fldir.ninf -= 1
                            if dp_fldir.Z >= max_Z:
                                dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
                    
                    else:
                        dp_fldir.upl = dr_net[curr_ch.value-1].length
                        dp_fldir.id_ch = dr_net[curr_ch.value-1].id_ch
                        dp_fldir.A_in += dp.A_in+1
                        dp_fldir.ninf -= 1
                        dr_net[curr_ch.value-1].id_ch_out = dp_fldir.id_ch
                        dr_pt_in[curr_fldir.value-1].ninf = 1
                        dr_pt_in[curr_fldir.value-1].inflow.append(dr_net[curr_ch.value-1].id_ch.value)
                        if dp_fldir.Z >= max_Z:
                            dr_net, dr_pt_in = dwnslp_hso(model, curr_fldir.value, dr_net, dr_pt_in, max_Z)
            dp.id_endo.value = -1
        
    return dr_net, dr_pt_in
                                



                                


                                           

                                        
    

                            
                            




                                        
                                    
                                
                    
                
            

