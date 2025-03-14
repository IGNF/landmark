# -*- coding: utf-8 -*-
"""
 Downslope Path Length Saddle Spill

The subroutine Downslope Path Legth calculates the legth of the path 
 between each DTM cell and the outflow point even if the basin is 
 endorheic. This task is crried out after the definition of saddle spill path
"""



def dpl_ss(model):
    for ep in model.l_endo_pt:
        if model.dr_pt[ep.id_pnt.value-1].upl != 0:

            if model.dr_net[model.dr_pt[ep.id_pnt.value-1].id_ch.value-1].id_end_pt.value == model.dr_pt[ep.id_pnt.value-1].id_pnt.value:
                main = model.dr_pt[ep.id_pnt.value-1].id_ch.value
                model.dr_net[main-1].id_endo = ep.id_eo
                model.dr_net[main-1].n_path = 1
                model.dr_net[main-1].id_path = [model.dr_net[main-1].id_ch.value]
                                
                for in_curr in model.dr_pt_in[ep.id_pnt.value-1].inflow:
                    if in_curr != main:
                        model.dr_net[in_curr-1].n_path = model.dr_net[main-1].n_path + 1
                        model.dr_net[in_curr-1].id_path = [model.dr_net[in_curr-1].id_ch.value]
                        model.dr_net[in_curr-1].id_path.extend(model.dr_net[main-1].id_path[:model.dr_net[in_curr-1].n_path-1])
                    
                    up_recurs_ss(in_curr, model)
    
    
    
def up_recurs_ss(curr, model):
    # print("up_recurs_ss")
    end_el = model.dr_net[curr-1].nel
    for cnt_pt in range(end_el-1, 0, -1):
        id_dr = model.dr_net[curr-1].id_pnts.value[cnt_pt-1]
        l1 = model.dr_net[curr-1].length
        l2 = model.dr_pt[id_dr-1].upl
        l3 = model.dr_pt[model.dr_net[curr-1].id_pnts.value[end_el-1]-1].dpl
        model.dr_pt[id_dr-1].dpl = l1 - l2 + l3
        
        if model.dr_pt_in[id_dr-1].ninf > 1:
            for in_curr in model.dr_pt_in[id_dr-1].inflow:
                if in_curr != curr:
                    model.dr_net[in_curr-1].n_path = model.dr_net[curr-1].n_path + 1
                    model.dr_net[in_curr-1].id_path = [model.dr_net[in_curr-1].id_ch.value]
                    model.dr_net[in_curr-1].id_path.extend(model.dr_net[curr-1].id_path[:model.dr_net[in_curr-1].n_path-1])
                    up_recurs_ss(in_curr, model)
    
            
        
    
                        