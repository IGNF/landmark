# -*- coding: utf-8 -*-
"""
Calculates the Downslope Path Length (DPL) for each drainage point, even in endorheic basins.
The path length is computed from each cell to the outflow point of its associated drainage basin.
"""

from tqdm import tqdm
import numpy as np

from data_structures import IDPointer


def dpl(model):
    """Computes the downslope path length (DPL) for all drainage points.

    This function processes each endorheic point to determine the corresponding main channel,
    and calculates the path length from each cell to the outflow for that basin. Tributary
    networks are recursively updated.

    Parameters
    ----------
    model : object
        Hydrological model containing at least the following attributes:

        - dr_pt : list of DrainagePoint
            All drainage points in the DEM.
        - dr_net : list of DrainageNetwork
            Drainage channels (each as a network).
        - mat_id : np.ndarray
            2D matrix for visualization and ID storage (size: (N*2-1, M*2-1)).
        - l_endo_pt : list of EndoPoint
            List of endorheic points.
        - N, M : int
            Dimensions of the DEM grid.
        - delta_x, delta_y : float
            Grid spacing.
        - nodata : float
            Nodata value used in the DEM.

    Returns
    -------
    None
    """
    # Process each endorheic point.
    for endo in tqdm(model.l_endo_pt, desc="Processing endorheic channels", unit="channel"):
        # Determine the main channel for this endorheic basin:
        # main is obtained from the DrainagePoint corresponding to the endo_pt.
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
                                                
                model.mat_id[i_mat, j_mat] = main_net.id_ch

        # For the current channel, set the downslope path identifier:
        main_net.n_path = 1
        main_net.id_path.append(main_net.id_ch.value)
        
        # Now call the recursive routine to update tributaries.
        up_recurs(main.value, model)
        

        


def up_recurs(curr, model):
    """Recursively updates tributary drainage networks with downslope path length.

    For each tributary channel:
    - Computes the downstream path length for each of its drainage points.
    - Updates the path matrix (mat_id).
    - Propagates the flow path and endorheic basin ID.

    Parameters
    ----------
    curr : int
        ID of the current (parent) channel.
    model : object
        Hydrological model object (same as in `dpl`).

    Returns
    -------
    None
    """
    for i in range(model.dr_net[curr-1].n_jun):
        in_curr = model.dr_net[curr-1].id_in.value[i-1]
        for cnt_pt in range(model.dr_net[in_curr.value-1].nel-1): #last pnt belongs to main channel
            dp = model.dr_pt[model.dr_net[in_curr.value-1].id_pnts.value[cnt_pt]-1]
            # l1 is the length of the tributary channel.
            l1 = model.dr_net[in_curr.value-1].length
            l2 = dp.upl
            l3 = model.dr_pt[model.dr_net[in_curr.value-1].id_end_pt.value-1].dpl
            model.dr_pt[model.dr_net[in_curr.value-1].id_pnts.value[cnt_pt]-1].dpl = l1-l2+l3
            i_curr = dp.i
            j_curr = dp.j
            if dp.fldir.value is not None and dp.fldir.value > 0:
                out_dp = model.dr_pt[dp.fldir.value-1]
                i_out, j_out = out_dp.i, out_dp.j
                i_mat = i_curr * 2 + (i_out - i_curr)
                j_mat = j_curr * 2 + (j_out - j_curr)

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
        
        up_recurs(in_curr.value, model)
        
        
        

        
            

            
