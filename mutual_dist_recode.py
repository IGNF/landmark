# -*- coding: utf-8 -*-
"""
Module to calculate the mutual distance (or "rejunction path") for each ridge point.
This code translates the Fortran subroutine 'mutual_dist' (and its helper 'md')
into Python. It assumes that the hydrological model (model) contains the following:
  - model.N, model.M: dimensions of the DEM.
  - model.mat_id: 2D numpy array (size (2*N-1, 2*M-1)) containing DrainagePoint objects.
  - model.dr_pt_by_id: dict mapping drainage point id (dp.id_pnt) to the DrainagePoint object.
  - model.dr_net_by_id: dict mapping channel id (net.id_ch) to the DrainageNetwork object.
  - model.dr_pt: list of all DrainagePoint objects.
  - model.dr_net: list of all DrainageNetwork objects.
"""

import time
from tqdm import tqdm

from data_structures import RidgePoint


def mutual_dist(model):
    """
    Computes the mutual distances (the "rejunction paths") for each ridge point.
    It creates a list of ridge points (rd_pt) and fills in their attributes based on the 
    drainage points in the matrix model.mat_id.
    
    The procedure is performed separately for:
      - Vertical cardinal directions,
      - Horizontal cardinal directions,
      - Diagonal directions.
    
    The function updates model.rd_pt (a list of RidgePoint objects) and prints timing information.
    
    Parameters
    ----------
    model : object
         The hydrological model object containing:
            - N, M: dimensions (used to iterate over the mat_id grid of size (2*N-1, 2*M-1))
            - mat_id: 2D numpy array with DrainagePoint objects at positions (i*2-1, j*2-1)
            - dr_pt_by_id: dict mapping drainage point id -> DrainagePoint
            - dr_net_by_id: dict mapping channel id -> DrainageNetwork
            - dr_pt: list of DrainagePoint objects.
    """

    cnt_rdpt = 0  # counter for ridge points
    rd_pt = []
    N2 = model.N * 2 - 1
    M2 = model.M * 2 - 1
    
    # cardinal direction vertical
    # o    drainage point (rr-1,cr)
    # +    ridge point  (rr,cr)
    # o    drainage point (rr+1,cr)
    
    for rr in tqdm(range(1, N2, 2)): #Row ridge
        for cr in range (0, M2,2):#Column ridge
            if model.mat_id[rr, cr] is not None:
                continue

            if (model.mat_id[rr - 1, cr] is not None) and (model.mat_id[rr + 1, cr] is not None):
                dp1 = model.mat_id[rr - 1, cr]
                dp2 = model.mat_id[rr + 1, cr]
                if dp1.Z < 0 or dp2.Z < 0:
                    continue
                else:
                    # Call md to compute mutual distance.
                    mtldst = md(dp1, dp2, model)
                    cnt_rdpt += 1
                    # Create a new ridge point.
                    # For the elevation, take the maximum of dp1.Z and dp2.Z.
                    new_rd = RidgePoint(
                        i=rr,
                        j=cr,
                        Z=max(dp1.Z, dp2.Z),
                        id_pnt=cnt_rdpt,
                    )
                    new_rd.md=mtldst

                    rd_pt.append(new_rd)
                    # Update mat_id: assign the ridge point to the cell.
                    model.mat_id[rr, cr] = new_rd
                    id_eo1 = model.l_dr_net[dp1.id_ch.value-1].id_endo.value
                    id_eo2 = model.l_dr_net[dp2.id_ch.value-1].id_endo.value
                    if (id_eo1 <= id_eo2):
                        new_rd.id_drpt1 = dp1.id_pnt
                        new_rd.id_drpt2 = dp2.id_pnt
                    else:
                        new_rd.id_drpt1 = dp2.id_pnt
                        new_rd.id_drpt2 = dp1.id_pnt



    # cardinal direction horizontal
    #  drainage point (rr-1,cr)
    #  |  /ridge point  (rr,cr)
    # 
    #  o + o
    # 
    #       \drainage point (rr+1,cr)
    for rr in tqdm(range(0, N2, 2)): #Row ridge
        for cr in range(1, M2, 2): #Column ridge
            if model.mat_id[rr, cr] is not None:
                continue

            if (model.mat_id[rr , cr - 1] is not None) and (model.mat_id[rr , cr + 1] is not None):
                dp1 = model.mat_id[rr, cr - 1]
                dp2 = model.mat_id[rr, cr + 1]
                if dp1.Z < 0 or dp2.Z < 0:
                    continue
                else:
                    mtldst = md(dp1, dp2, model)
                    cnt_rdpt += 1
                    # Create a new ridge point.
                    # For the elevation, take the maximum of dp1.Z and dp2.Z.
                    new_rd = RidgePoint(
                        i=rr,
                        j=cr,
                        Z=max(dp1.Z, dp2.Z),
                        id_pnt=cnt_rdpt,
                    )
                    new_rd.md=mtldst

                    rd_pt.append(new_rd)
                    # Update mat_id: assign the ridge point to the cell.
                    model.mat_id[rr, cr] = new_rd
                    id_eo1 = model.l_dr_net[dp1.id_ch.value-1].id_endo.value
                    id_eo2 = model.l_dr_net[dp2.id_ch.value-1].id_endo.value
                    if (id_eo1 <= id_eo2):
                        new_rd.id_drpt1 = dp1.id_pnt
                        new_rd.id_drpt2 = dp2.id_pnt
                    else:
                        new_rd.id_drpt1 = dp2.id_pnt
                        new_rd.id_drpt2 = dp1.id_pnt


    # diagonal directions 
    #  drainage point (rr-1,cr)
    #  |  /ridge point  (rr,cr)
    #  o   o
    #    + 
    #  o   o
    #       \drainage point (rr+1,cr)
    # 
    for rr in tqdm(range(1, N2,2)):
        for cr in range(1, M2, 2):
            if model.mat_id[rr, cr] is not None:
                continue

            # Check that all four diagonal neighbors exist.
            if (model.mat_id[rr - 1, cr - 1] is not None and
                model.mat_id[rr + 1, cr + 1] is not None and
                model.mat_id[rr + 1, cr - 1] is not None and
                model.mat_id[rr - 1, cr + 1] is not None):
                dp1 = model.mat_id[rr - 1, cr - 1]
                dp2 = model.mat_id[rr + 1, cr + 1]
                dp3 = model.mat_id[rr + 1, cr - 1]
                dp4 = model.mat_id[rr - 1, cr + 1]
                if (dp1 is None or dp2 is None or dp3 is None or dp4 is None or
                    dp1.Z < 0 or dp2.Z < 0 or dp3.Z < 0 or dp4.Z < 0):
                    continue
                else:
                    mtldst12 = md(dp1, dp2, model)
                    mtldst34 = md(dp3, dp4, model)  # Note: the Fortran code seems to call md with the same arguments twice.
                    mtldst_ar = [mtldst12, mtldst34]
                    cnt_rdpt += 1
                    new_rd = RidgePoint(
                        i=rr,
                        j=cr,
                        Z=max(dp1.Z, dp2.Z, dp3.Z, dp4.Z),
                        id_pnt=cnt_rdpt
                    )
                    new_rd.md=max(mtldst_ar)

                    rd_pt.append(new_rd)
                    model.mat_id[rr, cr] = new_rd
                    id_max = 0 if mtldst_ar[0] >= mtldst_ar[1] else 1
                    if id_max == 0:
                        id_eo1 = model.l_dr_net[dp1.id_ch.value-1].id_endo.value
                        id_eo2 = model.l_dr_net[dp1.id_ch.value-1].id_endo.value
                        if id_eo1 <= id_eo2:
                            new_rd.id_drpt1 = dp1.id_pnt
                            new_rd.id_drpt2 = dp2.id_pnt
                        else:
                            new_rd.id_drpt1 = dp2.id_pnt
                            new_rd.id_drpt2 = dp1.id_pnt
                    else:
                        id_eo3 = model.l_dr_net[dp3.id_ch.value-1].id_endo.value
                        id_eo4 = model.l_dr_net[dp4.id_ch.value-1].id_endo.value
                        if id_eo3 <= id_eo4:
                            new_rd.id_drpt1 = dp3.id_pnt
                            new_rd.id_drpt2 = dp4.id_pnt
                        else:
                            new_rd.id_drpt1 = dp4.id_pnt
                            new_rd.id_drpt2 = dp3.id_pnt

    model.rd_pt = rd_pt


def md(dp1, dp2, model):
    """
    Computes the mutual distance (mutdist) between two drainage points (dp1 and dp2)
    using the channel path information stored in their respective drainage networks.
    
    Parameters
    ----------
    dp1, dp2 : DrainagePoint
         The two drainage points (identified by their id_pnt) whose mutual distance is to be computed.
    model : object
         The hydrological model object containing:
           - dr_pt_by_id: dict mapping id_pnt -> DrainagePoint
           - dr_net_by_id: dict mapping channel id -> DrainageNetwork
           
    Returns
    -------
    mutdist : float
         The computed mutual distance. A very high value (e.g. 1e11) is returned if the two points belong to different basins.
    """
    curr_ch1 = dp1.id_ch
    curr_ch2 = dp2.id_ch
    # Access the drainage networks.
    net1 = model.l_dr_net[curr_ch1.value-1]
    net2 = model.l_dr_net[curr_ch2.value-1]
    if net1 is None or net2 is None:
        return 1e10


    n_path1 = net1.n_path
    n_path2 = net2.n_path
    min_n_pth = min(n_path1, n_path2)

    n_com=0  #Number of channels belonging to both path
    
    # Loop from min_n_pth down to 1.
    for cnt_path in range(min_n_pth-1, 0, -1): #!!!!!!!!!!!!!!!!! J'ai un doute sur le min_n_pth-1
        curr_cnt1 = cnt_path + n_path1 - min_n_pth  
        curr_cnt2 = cnt_path + n_path2 - min_n_pth
        # print("\nn_path1 : ", n_path1)
        # print("\nn_path2 : ", n_path2)

        
        # print("\ncurr_cnt1 : ", curr_cnt1)
        # print("\ncurr_cnt2 : ", curr_cnt2)

        # print(f'\nnet1.id_path : {[i.value for i in net1.id_path]}')
        # print(f'\nnet2.id_path : {[i.value for i in net2.id_path]}')
        if net1.id_path[curr_cnt1 - 1] == net2.id_path[curr_cnt2 - 1]:
            n_com += 1
    
    if n_com == 0:
        mutdist = 1e10
    
    else:
        if (n_path1 == n_path2) and (net1.id_path[0] == net2.id_path[0]):
            #Both points belogns to the same stream
            mutdist = abs(dp1.upl - dp2.upl)
        else:
            if (n_path1 - n_com > 0) and (n_path2 - n_com > 0):
                # Get the channel at position (n_path - n_com)
                jun1_channel_id = net1.id_path[n_path1 - n_com - 1].value
                jun2_channel_id = net2.id_path[n_path2 - n_com - 1].value
                #id of the first point of the path #1 belonging to the same path after the rejunction
                jun1 = model.l_dr_net[jun1_channel_id-1].id_end_pt.value
                #d of the first point of the path #2 belonging to the same path after the rejunction
                jun2 = model.l_dr_net[jun2_channel_id-1].id_end_pt.value
            else:
                #the two points belongs to the same basin
                if (n_path1 == min_n_pth) and (n_path1 != n_path2):
                    jun1 = dp1.id_pnt.value
                else:
                    jun1 = model.l_dr_net[jun1_channel_id-1].id_end_pt.value
                if (n_path2 == min_n_pth) and (n_path1 != n_path2):
                    jun2 = dp2.id_pnt.value
                else:
                    jun2 = model.l_dr_net[jun2_channel_id-1].id_end_pt.value
            if (jun1 is None) or (jun2 is None):
                mutdist = 1e11
            else:
                end_pt1 = model.dr_pt[jun1-1]
                end_pt2 = model.dr_pt[jun2-1]
                mutdist = dp1.dpl + dp2.dpl - 2 * min(end_pt1.dpl, end_pt2.dpl)
    return mutdist


