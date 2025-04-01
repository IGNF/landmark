# -*- coding: utf-8 -*-
"""
The subroutine VOLume and ARea allow to calculate for all the elements
of the thalweg network ordered with the Horton Stream Order,
the upslope area and the volume eventually added with the depression
filling procedure.
"""

import numpy as np


def vol_ar_hso(model):
    id_ch_rs = model.dr_net[0].id_path[model.dr_net[0].n_path-1]
    n_hso = model.dr_net[id_ch_rs-1].hso
    
    hso_arr = [i for i in range(1, n_hso+1)]
    cnt_ci = [0]*n_hso
    cnt_mis = [0]*n_hso
    A_arr = [0]*n_hso
    V_arr = [0]*n_hso
    V_mat = np.zeros((n_hso, n_hso))
    n_omega = np.zeros((n_hso, n_hso))
    ch_fv = [0] * len(model.dr_net)
    
    #Cumulative values calculation
    for curr_hso in range(n_hso):
        for cnt_el in range(len(model.dr_net)):
            curr_net = model.dr_net[cnt_el]
            if curr_net.hso == curr_hso:
                el_fv = 0
                curr_net.id_ch_out = model.dr_pt[curr_net.id_end_pt.value-1].id_ch
                
                A_arr[curr_hso] += model.dr_pt[curr_net.id_pnts.value[-2]-1].A_in + 1
                cnt_ci[curr_hso] += 1
                for id_dr in curr_net.id_pnt.value[:-1]:
                    i_curr = model.dr_pt[id_dr-1].i
                    j_curr = model.dr_pt[id_dr-1].j
                    

                    
    
    