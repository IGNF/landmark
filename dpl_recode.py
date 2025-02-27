# -*- coding: utf-8 -*-
"""
Transcript in python of the original dpl.f90

! The subroutine Downslope Path Legth calculates the legth of the path 
! between each DTM cell and the outflow point even if the basin is 
! endorheic 
"""

from tqdm import tqdm
import time
import sys
# sys.setrecursionlimit(10000)  # Augmenter la limite (ex: 5000 appels récursifs)



def dpl(model):
    """
    Calculates the downslope path length (dpl) for each drainage point and updates the
    drainage network (dr_net) information. This routine works over the endorheic points
    (model.endo_pt) and then, for each channel, computes the downslope path length for each
    drainage point, and then recursively updates tributary channels.
    
    Parameters
    ----------
    model : object
        The hydrological model object that contains at least:
           - dr_pt: list of DrainagePoint objects.
           - dr_pt_by_id: dict mapping id_pnt -> DrainagePoint.
           - dr_net: list of DrainageNetwork objects.
           - dr_net_by_id: dict mapping channel id (id_ch) -> DrainageNetwork.
           - mat_id: 2D numpy array (of size (N*2-1, M*2-1)) containing DrainagePoint objects.
           - endo_pt: list of EndoPoint objects.
           - N, M, delta_x, delta_y, nodata: global parameters.
    """
    start_time = time.process_time()
    
    # Process each endorheic point.
    for endo in tqdm(model.l_endo_pt, desc="Processing endorheic channels", unit="channel"):
        # Determine the main channel for this endorheic basin:
        # main is obtained from the DrainagePoint corresponding to the endo_pt.
        # (En Fortran : main => dr_pt(endo_pt(cnt_endo)%id_pnt)%id_ch)
        main = model.dr_pt[endo.id_pnt.value-1].id_ch
        # Assign the channel’s endorheic id:
        model.dr_net[main.value-1].id_endo = endo.id_eo

        # For each drainage point in the main channel:
            
        main_net = model.dr_net[main.value-1]
        l1 = main_net.length
        for dp_id in main_net.id_pnts.value:
            dp = model.dr_pt[dp_id-1]
            l2 = dp.upl
            dp.dpl = l1 - l2
            
            i_curr = dp.i
            j_curr = dp.j
            
        
            if dp.fldir.value is not None and dp.fldir.value > 0:
                # Retrieve the outflow point.
                out_dp = model.dr_pt[dp.fldir.value-1]
                i_out, j_out = out_dp.i, out_dp.j
                i_mat = i_curr * 2  + (i_out - i_curr)
                j_mat = j_curr * 2  + (j_out - j_curr)
                                
                # Update the corresponding cell in mat_id with the main channel id.
                # if model.mat_id[i_mat, j_mat] is None:
                #     print("\ni_mat", i_mat)
                #     print("j_mat", j_mat)
                # if model.mat_id[i_mat, j_mat] is not None:
                #     model.mat_id[i_mat, j_mat].id_pnt = main_net.id_ch
                
                model.mat_id[i_mat, j_mat] = main_net.id_ch


        # For the current channel, set the downslope path identifier:
        main_net.n_path = 1
        main_net.id_path.append(main_net.id_ch.value)
        
        # Now call the recursive routine to update tributaries.
        up_recurs_recode(main.value, model)
        
    
    finish_time = time.process_time()
    elapsed_time = finish_time - start_time
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = elapsed_time % 60
    print(f"Elapsed time (dpl): {hours}h {minutes}m {seconds:.2f}s")




def up_recurs_recode(curr, model, visited=None):
    """
    Recursively processes tributary channels for downslope path length updates.
    For each tributary channel of the channel identified by 'curr', this function:
      - Updates the downslope path length (dpl) for its drainage points,
      - Updates the associated matrix (mat_id),
      - Propagates the cumulative path (id_path) from the parent channel,
      - Sets the endorheic id to be the same as the parent's,
      - Recursively calls itself on the tributary channel.
    
    Parameters
    ----------
    curr : int
        The channel id (id_ch) of the current (parent) drainage network.
    model : object
        The hydrological model object containing dr_pt, dr_net, etc.
    """
    # if visited is None:
    #     visited = set()

    # # Si on est déjà passé par ce canal, on ne le traite plus
    # if curr in visited:
    #     return
    # visited.add(curr)

    
    for i in range(model.dr_net[curr-1].n_jun):
        in_curr = model.dr_net[curr-1].id_in.value[i-1]
        for cnt_pt in range(model.dr_net[in_curr.value-1].nel-1): #last pnt belongs to main channel
        #Marque tous les points sur les lignes de flux pour ne pas apparaitre comme points de crête pour la suite.
            dp = model.dr_pt[model.dr_net[in_curr.value-1].id_pnts.value[cnt_pt]-1]
            # l1 is the length of the tributary channel.
            l1 = model.dr_net[in_curr.value-1].length
            l2 = dp.upl
            l3 = model.dr_pt[model.dr_net[in_curr.value-1].id_end_pt.value-1].dpl
            model.dr_pt[model.dr_net[in_curr.value-1].id_pnts.value[cnt_pt]-1].dpl = l1-l2+l3
            i_curr = dp.i
            j_curr = dp.j
            # print(i_curr)
            if dp.fldir.value is not None and dp.fldir.value > 0:
                out_dp = model.dr_pt[dp.fldir.value-1]
                i_out, j_out = out_dp.i, out_dp.j
                i_mat = i_curr * 2 + (i_out - i_curr)
                j_mat = j_curr * 2 + (j_out - j_curr)
                # print("\ni_mat", i_mat)
                # print("j_mat", j_mat)

                model.mat_id[i_mat, j_mat] = model.dr_net[in_curr.value-1].id_ch
        
        # Update the tributary channel’s path:
        model.dr_net[in_curr.value-1].n_path = model.dr_net[curr-1].n_path+1
        # Build new id_path: first element is its own id_ch, then append parent's path.
        n_path = model.dr_net[in_curr.value-1].n_path
        # model.dr_net[in_curr.value-1].id_path = []
        model.dr_net[in_curr.value-1].id_path.append(model.dr_net[in_curr.value-1].id_ch.value)
        for cnt_in in range(2,n_path+1):
            model.dr_net[in_curr.value-1].id_path.append(model.dr_net[curr-1].id_path[cnt_in-2])
        # Propagate the endorheic id from the parent.
        model.dr_net[in_curr.value-1].id_endo = model.dr_net[curr-1].id_endo
        
        up_recurs_recode(in_curr.value, model, visited)
        
        
        

        
            

            
