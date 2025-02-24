# -*- coding: utf-8 -*-
"""
The subroutine Ridge Point define the relationship between ridge points
as they are ready to be used to extract an ordered ridge network.
Modify fldir
"""


from tqdm import tqdm

from mutual_dist_recode import md



def ridge_point(model):
    n_rdpnt = len(model.rd_pt)
    n_rdpnt_old=n_rdpnt
    
    for cnt_rdpt in tqdm(range(n_rdpnt)):
        rr = model.rd_pt[cnt_rdpt].i
        cr = model.rd_pt[cnt_rdpt].j
        nen = 0
        taonp = []
        
        #  CARDINAL POINTS
        # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
        #  .  o  .   -1
        #  .  +  .   rr
        #  .  o  .   +1
        # 
        #  -1 cr +1
        
        


