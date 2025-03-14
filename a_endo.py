# -*- coding: utf-8 -*-
"""
The subroutine A_ENDO calculates the area of the portion of upslope
basin of each endorehic saddle spill that is less then OR equal to
the elevation of the spill
"""


import numpy as np

from data_structures import IDPointer, DrainagePoint


def a_endo(model):
    # river_mask = np.vectorize(lambda x: isinstance(x, IDPointer))(model.mat_id)
    drainage_mask = np.vectorize(lambda x: isinstance(x, DrainagePoint))(model.mat_id)


    for sp in model.sdl_pt:
        if sp.id_cis_endo.value != None:
            Z_saddle = model.rd_pt[sp.id_rdpt-1].Z
            cnt_fifo = 1
            pos_fifo = 1
            fifo = [sp.id_cis_endo.value]
            while cnt_fifo >= pos_fifo:
                i_curr = model.dr_pt[fifo[pos_fifo-1]-1].i
                j_curr = model.dr_pt[fifo[pos_fifo-1]-1].j
                
                for ir_dtm in range(i_curr-1, i_curr+2):
                    for ic_dtm in range(j_curr-1, j_curr+2):
                        ir = ir_dtm*2
                        ic = ic_dtm*2
                        if ir>=0 and ir< model.N*2-1 and ic >= 0 and ic < model.M*2-1:
                            if drainage_mask[ir, ic]:
                                if model.mat_id[ir, ic].fldir.value == fifo[pos_fifo-1]:
                                    if model.mat_id[ir, ic].Z <= Z_saddle:
                                        cnt_fifo += 1
                                        fifo.append(model.mat_id[ir, ic].id_pnt.value)
                pos_fifo += 1
            sp.A_endo = cnt_fifo
    

                