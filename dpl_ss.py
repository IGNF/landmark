# -*- coding: utf-8 -*-
"""
 Downslope Path Length Saddle Spill

The subroutine Downslope Path Legth calculates the length of the path 
 between each DTM cell and the outflow point even if the basin is 
 endorheic. This task is crried out after the definition of saddle spill path
"""

from tqdm import tqdm



def dpl_ss(model):
    """Computes downslope path length (dpl) for all drainage points connected via
    subsurface flow (`fldir_ss`) starting from endorheic basin outlets.

    This is the subsurface complement to `dpl()` using the upstream inflow structure
    stored in `model.dr_pt_in`.

    Parameters
    ----------
    model : object
        Hydrological model containing:
        - l_endo_pt : list of EndoPoint
        - dr_pt : list of DrainagePoint
        - dr_net : list of DrainageNetwork
        - dr_pt_in : list of DrainagePointInflow
    """
    for ep in tqdm(model.l_endo_pt):
        if model.dr_pt[ep.id_pnt-1].upl != 0:

            if model.dr_net[model.dr_pt[ep.id_pnt-1].id_ch-1].id_end_pt == model.dr_pt[ep.id_pnt-1].id_pnt:
                main = model.dr_pt[ep.id_pnt-1].id_ch
                model.dr_net[main-1].id_endo = ep.id_eo
                model.dr_net[main-1].n_path = 1
                model.dr_net[main-1].id_path = [model.dr_net[main-1].id_ch]
                                
                for in_curr in model.dr_pt_in[ep.id_pnt-1].inflow:
                    if in_curr != main:
                        model.dr_net[in_curr-1].n_path = model.dr_net[main-1].n_path + 1
                        model.dr_net[in_curr-1].id_path = [model.dr_net[in_curr-1].id_ch]
                        model.dr_net[in_curr-1].id_path.extend(model.dr_net[main-1].id_path[:model.dr_net[in_curr-1].n_path-1])
                    
                    up_recurs_ss(in_curr, model)
    
    del(model.dr_pt_in)
    
    
def up_recurs_ss(curr, model):
    """Recursively updates downslope path length for all channels upstream of `curr`,
    following the `dr_pt_in` inflow graph (subsurface).

    Parameters
    ----------
    curr : int
        ID of the current drainage network being updated.
    model : object
        Hydrological model (see above).
    """
    end_el = model.dr_net[curr - 1].nel
    for cnt_pt in range(end_el - 1, 0, -1):
        id_dr = model.dr_net[curr - 1].id_pnts[cnt_pt - 1]
        l1 = model.dr_net[curr - 1].length
        l2 = model.dr_pt[id_dr - 1].upl
        l3 = model.dr_pt[model.dr_net[curr - 1].id_pnts[end_el - 1] - 1].dpl

        model.dr_pt[id_dr - 1].dpl = l1 - l2 + l3

        if model.dr_pt_in[id_dr - 1].ninf > 1:
            for in_curr in model.dr_pt_in[id_dr - 1].inflow:
                if in_curr != curr:
                    inflow_net = model.dr_net[in_curr - 1]
                    inflow_net.n_path = model.dr_net[curr - 1].n_path + 1
                    inflow_net.id_path = [inflow_net.id_ch]
                    inflow_net.id_path.extend(model.dr_net[curr - 1].id_path[:inflow_net.n_path - 1])

                    up_recurs_ss(in_curr, model)
            
        
    
                        