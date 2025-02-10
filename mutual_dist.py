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


#---------------------------------------------------------------------------
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
    # Retrieve the channel id for each drainage point.
    curr_ch1 = dp1.id_ch
    curr_ch2 = dp2.id_ch
    # Access the drainage networks.
    net1 = model.dr_net_by_id.get(curr_ch1)
    net2 = model.dr_net_by_id.get(curr_ch2)
    if net1 is None or net2 is None:
        return 1e11
    n_path1 = net1.n_path
    n_path2 = net2.n_path
    min_n_pth = min(n_path1, n_path2)
    n_com = 0  # number of channels common to both paths
    
    # Loop from min_n_pth down to 1.
    for cnt_path in range(min_n_pth, 0, -1):
        # Adjust indices for 0-indexing.
        curr_cnt1 = cnt_path + n_path1 - min_n_pth  # 1-indexed equivalent
        curr_cnt2 = cnt_path + n_path2 - min_n_pth
        # In Python, list indices start at 0.
        if net1.id_path[curr_cnt1 - 1] == net2.id_path[curr_cnt2 - 1]:
            n_com += 1

    if n_com == 0:
        mutdist = 1e11
    else:
        if (n_path1 == n_path2) and (net1.id_path[0] == net2.id_path[0]):
            mutdist = abs(dp1.upl - dp2.upl)
        else:
            if (n_path1 - n_com > 0) and (n_path2 - n_com > 0):
                # Get the channel at position (n_path - n_com)
                jun1_channel_id = net1.id_path[n_path1 - n_com - 1]
                jun2_channel_id = net2.id_path[n_path2 - n_com - 1]
                jun1 = model.dr_pt_by_id.get(model.dr_net_by_id[jun1_channel_id].id_end_pt)
                jun2 = model.dr_pt_by_id.get(model.dr_net_by_id[jun2_channel_id].id_end_pt)
            else:
                if (n_path1 == min_n_pth) and (n_path1 != n_path2):
                    jun1 = dp1
                else:
                    jun1 = model.dr_pt_by_id.get(model.dr_net_by_id[net1.id_path[n_path1 - n_com - 1]].id_end_pt)
                if (n_path2 == min_n_pth) and (n_path1 != n_path2):
                    jun2 = dp2
                else:
                    jun2 = model.dr_pt_by_id.get(model.dr_net_by_id[net2.id_path[n_path2 - n_com - 1]].id_end_pt)
            if (jun1 is None) or (jun2 is None):
                mutdist = 1e11
            else:
                mutdist = dp1.dpl + dp2.dpl - 2 * min(jun1.dpl, jun2.dpl)
    return mutdist

#---------------------------------------------------------------------------
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
    start_time = time.process_time()
    cnt_rdpt = 0  # counter for ridge points
    # Allocate a list to hold ridge points; size maximum is 2*N*2*M.
    rd_pt = []
    N2 = model.N * 2 - 1
    M2 = model.M * 2 - 1
    
    # --- Vertical cardinal direction ---
    # For rows from 2 to (N*2-1) step 2, columns from 1 to (M*2-1) step 2.
    for rr in tqdm(range(1, N2, 2)):   # In Fortran, rows: 2,4,...,N*2-1; here 0-index: 1,3,...
        for cr in range(0, M2, 2):  # Fortran: 1,3,..., M*2-1; here 0-index: 0,2,...
            # If the current cell in mat_id does not have an associated drainage point,
            # it is considered a ridge point.
            if model.mat_id[rr, cr] is not None:
                continue
            # Check if the cells above and below are associated.
            if (model.mat_id[rr - 1, cr] is not None) and (model.mat_id[rr + 1, cr] is not None):
                pos1 = model.mat_id[rr - 1, cr].id_pnt
                pos2 = model.mat_id[rr + 1, cr].id_pnt
                dp1 = model.dr_pt_by_id.get(pos1)
                dp2 = model.dr_pt_by_id.get(pos2)
                if (dp1 is None or dp2 is None) or (dp1.Z < 0 or dp2.Z < 0):
                    continue
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
                    md=mtldst,
                    nen=0
                )
                rd_pt.append(new_rd)
                # Update mat_id: assign the ridge point to the cell.
                model.mat_id[rr, cr] = new_rd
                # For the drainage points associated with this ridge point, assign their ids.
                # Use the smaller id first.
                if dp1.Z <= dp2.Z:
                    new_rd.id_drpt1 = dp1.id_pnt
                    new_rd.id_drpt2 = dp2.id_pnt
                else:
                    new_rd.id_drpt1 = dp2.id_pnt
                    new_rd.id_drpt2 = dp1.id_pnt
    # --- Horizontal cardinal direction ---
    for rr in tqdm(range(0, N2, 2)):   # rows: 1,3,... in Fortran -> here indices 0,2,...
        for cr in range(1, M2, 2):  # columns: 2,4,... -> here 1,3,...
            if model.mat_id[rr, cr] is not None:
                continue
            if (model.mat_id[rr, cr - 1] is not None) and (model.mat_id[rr, cr + 1] is not None):
                pos1 = model.mat_id[rr, cr - 1].id_pnt
                pos2 = model.mat_id[rr, cr + 1].id_pnt
                dp1 = model.dr_pt_by_id.get(pos1)
                dp2 = model.dr_pt_by_id.get(pos2)
                if (dp1 is None or dp2 is None) or (dp1.Z < 0 or dp2.Z < 0):
                    continue
                mtldst = md(dp1, dp2, model)
                cnt_rdpt += 1
                new_rd = RidgePoint(
                    i=rr,
                    j=cr,
                    Z=max(dp1.Z, dp2.Z),
                    id_pnt=cnt_rdpt,
                    md=mtldst,
                    nen=0
                )
                rd_pt.append(new_rd)
                model.mat_id[rr, cr] = new_rd
                if dp1.Z <= dp2.Z:
                    new_rd.id_drpt1 = dp1.id_pnt
                    new_rd.id_drpt2 = dp2.id_pnt
                else:
                    new_rd.id_drpt1 = dp2.id_pnt
                    new_rd.id_drpt2 = dp1.id_pnt
    # --- Diagonal directions ---
    for rr in tqdm(range(1, N2 - 1, 2)):  # rows: 2 to N*2-2 in Fortran -> here indices 1 to N2-2 (step 2)
        for cr in range(1, M2 - 1, 2):  # columns: similarly
            if model.mat_id[rr, cr] is not None:
                continue
            # Check that all four diagonal neighbors exist.
            if (model.mat_id[rr - 1, cr - 1] is not None and
                model.mat_id[rr + 1, cr + 1] is not None and
                model.mat_id[rr + 1, cr - 1] is not None and
                model.mat_id[rr - 1, cr + 1] is not None):
                pos1 = model.mat_id[rr - 1, cr - 1].id_pnt
                pos2 = model.mat_id[rr + 1, cr + 1].id_pnt
                pos3 = model.mat_id[rr + 1, cr - 1].id_pnt
                pos4 = model.mat_id[rr - 1, cr + 1].id_pnt
                dp1 = model.dr_pt_by_id.get(pos1)
                dp2 = model.dr_pt_by_id.get(pos2)
                dp3 = model.dr_pt_by_id.get(pos3)
                dp4 = model.dr_pt_by_id.get(pos4)
                if (dp1 is None or dp2 is None or dp3 is None or dp4 is None or
                    dp1.Z < 0 or dp2.Z < 0 or dp3.Z < 0 or dp4.Z < 0):
                    continue
                # For the diagonal case, the Fortran code calls md twice to get two mutual distances,
                # here we assume mtldst12 and mtldst34.
                mtldst12 = md(dp1, dp2, model)
                mtldst34 = md(dp1, dp2, model)  # Note: the Fortran code seems to call md with the same arguments twice.
                mtldst_ar = [mtldst12, mtldst34]
                cnt_rdpt += 1
                new_rd = RidgePoint(
                    i=rr,
                    j=cr,
                    Z=max(dp1.Z, dp2.Z, dp3.Z, dp4.Z),
                    id_pnt=cnt_rdpt,
                    md=max(mtldst_ar),
                    nen=0
                )
                rd_pt.append(new_rd)
                model.mat_id[rr, cr] = new_rd
                id_max = 0 if mtldst_ar[0] >= mtldst_ar[1] else 1
                if id_max == 0:
                    id_eo1 = model.dr_net_by_id.get(dp1.id_ch).id_endo
                    id_eo2 = model.dr_net_by_id.get(dp2.id_ch).id_endo
                    if id_eo1 <= id_eo2:
                        new_rd.id_drpt1 = dp1.id_pnt
                        new_rd.id_drpt2 = dp2.id_pnt
                    else:
                        new_rd.id_drpt1 = dp2.id_pnt
                        new_rd.id_drpt2 = dp1.id_pnt
                else:
                    id_eo1 = model.dr_net_by_id.get(dp1.id_ch).id_endo
                    id_eo2 = model.dr_net_by_id.get(dp2.id_ch).id_endo
                    if id_eo1 <= id_eo2:
                        new_rd.id_drpt1 = dp1.id_pnt
                        new_rd.id_drpt2 = dp2.id_pnt
                    else:
                        new_rd.id_drpt1 = dp2.id_pnt
                        new_rd.id_drpt2 = dp1.id_pnt
    # Save the ridge points list into model.
    model.rd_pt = rd_pt
    model.n_rdpnt = cnt_rdpt
    
    print(f"Number of ridge points: {cnt_rdpt}")
    
    # Timing information.
    finish_time = time.process_time()
    elapsed_time = finish_time - start_time
    sph = 3600.0
    ms = 60.0
    hours = int(elapsed_time // sph)
    minutes = int((elapsed_time % sph) // ms)
    seconds = elapsed_time % ms
    print(f"Mutual distance elapsed time: {hours}h {minutes}m {seconds:.2f}s")
