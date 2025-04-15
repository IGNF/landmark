# -*- coding: utf-8 -*-
"""
Update junc value to delineate bassins
"""


from tqdm import tqdm
import numpy as np

from data_structures import RidgePoint



def junc_update(model):
    """
    Identifies and marks junction points in the ridgeline network.

    A junction is defined as a ridge cell that connects two significant ridgelines,
    typically at a saddle or confluence area. This routine updates the `junc` attribute
    of RidgePoint objects stored in `model.mat_id` and `model.rd_pt`.

    The method inspects the last two drainage points of each network element,
    deduces the direction, and searches for intersecting ridge points in the corresponding
    neighborhood to mark as junctions.

    Parameters
    ----------
    model : object
        Hydrological model containing:
        - dr_net : list of DrainageNetwork
        - mat_id : 2D grid of RidgePoint or DrainagePoint
        - delta_x, delta_y : float
        - hso_th : int, minimum Horton Stream Order threshold
        - a_out_threshold : float
        - rd_pt : list of RidgePoint
        - sdl_pt : list of SaddlePoint
    """

    # Compute the area threshold in cell units
    ncell_th = model.a_out_threshold / model.delta_x / model.delta_y

    # Create a mask for identifying RidgePoints in the grid
    ridge_mask = np.vectorize(lambda x: isinstance(x, RidgePoint))(model.mat_id)

    for net in tqdm(model.dr_net):
        if net.hso >= model.hso_th:
            n_pnt = net.nel

            # Count how many points exceed the area threshold
            if ncell_th > 0:
                n_pnt_A_th = 0
                for cnt_pnt in range(n_pnt):
                    if model.dr_pt[net.id_pnts[cnt_pnt] - 1].A_in + 1 >= ncell_th:
                        n_pnt_A_th += 1
            else:
                n_pnt_A_th = n_pnt

            if n_pnt_A_th > 1:
                # Get coordinates of last and second-to-last drainage points
                i_last_dtm = model.dr_pt[net.id_pnts[-1] - 1].i
                j_last_dtm = model.dr_pt[net.id_pnts[-1] - 1].j
                i_slast_dtm = model.dr_pt[net.id_pnts[-2] - 1].i
                j_slast_dtm = model.dr_pt[net.id_pnts[-2] - 1].j

                i_passo_dtm = i_slast_dtm - i_last_dtm
                j_passo_dtm = j_slast_dtm - j_last_dtm

                # ▸ Case: Horizontal direction
                if i_passo_dtm == 0:
                    passo_mem = 0
                    i_passo = 1
                    j_passo = j_passo_dtm
                    i_last = i_last_dtm * 2 - 1
                    i_slast = i_slast_dtm * 2 + 1
                    j_last = j_last_dtm * 2 + j_passo
                    j_slast = j_slast_dtm * 2

                    set = 0
                    ir = i_last
                    for ic in range(j_last, j_slast + j_passo, j_passo):
                        if set == 0 and 0 <= ir < model.mat_id.shape[0] and 0 <= ic < model.mat_id.shape[1]:
                            if ridge_mask[ir, ic]:
                                model.mat_id[ir, ic].junc = 1
                                set = 1

                    set = 0
                    ir = i_slast
                    for ic in range(j_last, j_slast + j_passo, j_passo):
                        if set == 0 and 0 <= ir < model.mat_id.shape[0] and 0 <= ic < model.mat_id.shape[1]:
                            if ridge_mask[ir, ic]:
                                model.mat_id[ir, ic].junc = 1
                                set = 1

                else:
                    # ▸ Case: Diagonal or vertical direction
                    passo_mem = i_passo_dtm
                    i_passo = i_passo_dtm
                    j_passo = j_passo_dtm
                    i_last = i_last_dtm * 2 + i_passo
                    i_slast = i_slast_dtm * 2
                    j_last = j_last_dtm * 2 + j_passo
                    j_slast = j_slast_dtm * 2

                # ▸ Case: Vertical
                if j_passo_dtm == 0:
                    passo_mem = 0
                    i_passo = i_passo_dtm
                    j_passo = 1
                    i_last = i_last_dtm * 2 + i_passo
                    i_slast = i_slast_dtm * 2
                    j_last = j_last_dtm * 2 - 1
                    j_slast = j_slast_dtm * 2 + 1

                    set = 0
                    ic = j_last
                    for ir in range(i_last, i_slast + i_passo, i_passo):
                        if set == 0 and 0 <= ir < model.mat_id.shape[0] and 0 <= ic < model.mat_id.shape[1]:
                            if ridge_mask[ir, ic]:
                                model.mat_id[ir, ic].junc = 1
                                set = 1

                    set = 0
                    ic = j_slast
                    for ir in range(i_last, i_slast + i_passo, i_passo):
                        if set == 0 and 0 <= ir < model.mat_id.shape[0] and 0 <= ic < model.mat_id.shape[1]:
                            if ridge_mask[ir, ic]:
                                model.mat_id[ir, ic].junc = 1
                                set = 1

                else:
                    passo_mem = i_passo_dtm
                    i_passo = i_passo_dtm
                    j_passo = j_passo_dtm
                    i_last = i_last_dtm * 2 + i_passo
                    i_slast = i_slast_dtm * 2
                    j_last = j_last_dtm * 2 + j_passo
                    j_slast = j_slast_dtm * 2

                # ▸ Diagonal ridgeline segment: check for intersections
                if passo_mem != 0:
                    # Search vertical neighborhood
                    for ir in range(i_last - i_passo, i_slast + 2 * i_passo, i_passo):
                        ic = j_last + j_passo
                        if 0 <= ir < model.mat_id.shape[0] and 0 <= ic < model.mat_id.shape[1]:
                            if ridge_mask[ir, ic]:
                                model.mat_id[ir, ic].junc = 1
                    # Search horizontal neighborhood
                    for ic in range(j_last - j_passo, j_slast + 2 * j_passo, j_passo):
                        ir = i_last + i_passo
                        if 0 <= ir < model.mat_id.shape[0] and 0 <= ic < model.mat_id.shape[1]:
                            if ridge_mask[ir, ic]:
                                model.mat_id[ir, ic].junc = 1

    # ▸ Final: mark junctions from saddle points if spread area is large enough
    for sp in model.sdl_pt:
        if sp.id_rdpt2 > 0 and sp.A_endo > ncell_th:
            model.rd_pt[sp.id_rdpt - 1].junc = 1
            model.rd_pt[sp.id_rdpt2 - 1].junc = 1
