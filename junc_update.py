# -*- coding: utf-8 -*-
"""
Update junc value to delineate bassins
"""


from tqdm import tqdm
import numpy as np

from data_structures import RidgePoint



def junc_update(model):
    ncell_th = model.a_out_threshold/model.delta_x/model.delta_y
    ridge_mask = np.vectorize(lambda x: isinstance(x, RidgePoint))(model.mat_id)

    for net in tqdm(model.dr_net):
        if net.hso >= model.hso_th:
            n_pnt = net.nel
            if ncell_th > 0:
                n_pnt_A_th = 0
                for cnt_pnt in range(n_pnt):
                    if model.dr_pt[net.id_pnts.value[cnt_pnt]-1].A_in + 1 >= ncell_th:
                        n_pnt_A_th += 1
            
            else:
                n_pnt_A_th = n_pnt
            
            if n_pnt_A_th > 1 :
                i_last_dtm = model.dr_pt[net.id_pnts.value[-1]-1].i
                j_last_dtm = model.dr_pt[net.id_pnts.value[-1]-1].j
                i_slast_dtm = model.dr_pt[net.id_pnts.value[-2]-1].i
                j_slast_dtm = model.dr_pt[net.id_pnts.value[-2]-1].j
                
                i_passo_dtm = i_slast_dtm - i_last_dtm
                j_passo_dtm = j_slast_dtm - j_last_dtm
                
                if i_passo_dtm == 0:
                    passo_mem = 0
                    i_passo = 1
                    j_passo = j_passo_dtm
                    i_last = i_last_dtm *2 - 1
                    i_slast = i_slast_dtm * 2 +1
                    j_last=j_last_dtm*2 + j_passo
                    j_slast=j_slast_dtm*2
                    set=0
                    ir=i_last           
                    for ic in range (j_last, j_slast+j_passo, j_passo):
                        if set == 0:
                            if ir>=0 and ir <= model.mat_id.shape[0]-1 and ic>=0 and ic <= model.mat_id.shape[1]-1:
                                if ridge_mask[ir, ic]:
                                    model.mat_id[ir, ic].junc = 1
                                    set = 1
                    
                    set = 0
                    ir = i_slast
                    for ic in range (j_last, j_slast+j_passo, j_passo):
                        if set == 0:
                            if ir>=0 and ir <= model.mat_id.shape[0]-1 and ic>=0 and ic <= model.mat_id.shape[1]-1:
                                if ridge_mask[ir, ic]:
                                    model.mat_id[ir, ic].junc = 1
                                    set = 1
                    
                else :
                    passo_mem=i_passo_dtm
                    i_passo=i_passo_dtm
                    j_passo=j_passo_dtm
                    i_last=i_last_dtm*2+i_passo
                    i_slast=i_slast_dtm*2
                    j_last=j_last_dtm*2+j_passo
                    j_slast=j_slast_dtm*2
                
                if j_passo_dtm == 0:
                    passo_mem = 0
                    i_passo = i_passo_dtm
                    j_passo = 1
                    i_last = i_last_dtm*2 + i_passo
                    i_slast = i_slast_dtm*2
                    j_last = j_last_dtm*2 - 1
                    j_slast = j_slast_dtm*2 + 1
                    set = 0
                    ic=j_last
                    for ir in range(i_last, i_slast+i_passo, i_passo):
                        if set == 0:
                            if ir>=0 and ir <= model.mat_id.shape[0]-1 and ic>=0 and ic <= model.mat_id.shape[1]-1:
                                if ridge_mask[ir, ic]:
                                    model.mat_id[ir, ic].junc = 1
                                    set = 1
                    
                    set = 0
                    ic = j_slast
                    for ir in range(i_last, i_slast+i_passo, i_passo):
                        if set == 0:
                            if ir>=0 and ir <= model.mat_id.shape[0]-1 and ic>=0 and ic <= model.mat_id.shape[1]-1:
                                if ridge_mask[ir, ic]:
                                    model.mat_id[ir, ic].junc = 1
                                    set = 1
                
                else:
                    passo_mem=i_passo_dtm
                    i_passo=i_passo_dtm
                    j_passo=j_passo_dtm
                    i_last=i_last_dtm*2 + i_passo
                    i_slast=i_slast_dtm*2
                    j_last=j_last_dtm*2 + j_passo
                    j_slast=j_slast_dtm*2
                
                if passo_mem != 0:
                    # diagonal direction dr_net element
                    # nodes placed on the vertical
                    for ir in range(i_last - i_passo, i_slast + 2*i_passo, i_passo):
                        ic = j_last + j_passo
                        if ir>=0 and ir <= model.mat_id.shape[0]-1 and ic>=0 and ic <= model.mat_id.shape[1]-1:
                            if ridge_mask[ir, ic]:
                                model.mat_id[ir, ic].junc = 1
                    # nodes placed on the horizontal
                    for ic in range(j_last - j_passo, j_slast + 2*j_passo, j_passo):
                        ir=i_last + i_passo
                        if ir>=0 and ir <= model.mat_id.shape[0]-1 and ic>=0 and ic <= model.mat_id.shape[1]-1:
                            if ridge_mask[ir, ic]:
                                model.mat_id[ir, ic].junc = 1
                            
                   
    for sp in model.sdl_pt:
        if sp.id_rdpt2 > 0 and sp.A_endo > ncell_th :
            model.rd_pt[sp.id_rdpt-1].junc = 1
            model.rd_pt[sp.id_rdpt2-1].junc = 1



            
        