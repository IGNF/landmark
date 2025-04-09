# -*- coding: utf-8 -*-
"""
The subroutine A_ENDO calculates the area of the portion of upslope
basin of each endorehic saddle spill that is less then OR equal to
the elevation of the spill
"""


import numpy as np

from data_structures import DrainagePoint


def a_endo(model):
    """
    Computes the internal contributing area (A_endo) of endorheic basins
    for each saddle point.

    This function uses a flood-fill (FIFO-style) approach to explore all drainage
    points upslope of the saddle, but only within the same endorheic basin.
    A point is included if:
        - it is a drainage point (in `mat_id`),
        - it flows into the current point (via `fldir`),
        - its elevation is lower than or equal to the saddle elevation.

    Parameters
    ----------
    model : object
        Hydrological model containing:
        - mat_id : 2D ndarray with DrainagePoint or RidgePoint.
        - N, M : int, grid size.
        - dr_pt : list of DrainagePoint
        - sdl_pt : list of SaddlePoint
        - rd_pt : list of RidgePoint (contains Z of saddle via id_rdpt)

    Returns
    -------
    None
        Updates each SaddlePoint's `A_endo` attribute with the cell count.
    """

    # Create a boolean mask indicating which grid cells are DrainagePoints
    drainage_mask = np.vectorize(lambda x: isinstance(x, DrainagePoint))(model.mat_id)

    # Iterate over each saddle point
    for sp in model.sdl_pt:
        if sp.id_cis_endo.value is not None:

            # Get elevation of the saddle (threshold for contributing points)
            Z_saddle = model.rd_pt[sp.id_rdpt - 1].Z

            # Initialize FIFO queue (manual implementation)
            cnt_fifo = 1
            pos_fifo = 1
            fifo = [sp.id_cis_endo.value]  # Start from the point just upstream of the saddle

            # Begin flood-fill loop
            while cnt_fifo >= pos_fifo:
                id_curr = fifo[pos_fifo - 1]
                i_curr = model.dr_pt[id_curr - 1].i
                j_curr = model.dr_pt[id_curr - 1].j

                # Explore 3x3 neighborhood around current point
                for ir_dtm in range(i_curr - 1, i_curr + 2):
                    for ic_dtm in range(j_curr - 1, j_curr + 2):
                        ir = ir_dtm * 2
                        ic = ic_dtm * 2

                        if 0 <= ir < model.N * 2 - 1 and 0 <= ic < model.M * 2 - 1:
                            if drainage_mask[ir, ic]:
                                neighbor = model.mat_id[ir, ic]
                                # Check if neighbor flows to current cell and is below saddle
                                if neighbor.fldir.value == id_curr and neighbor.Z <= Z_saddle:
                                    cnt_fifo += 1
                                    fifo.append(neighbor.id_pnt.value)

                pos_fifo += 1

            # Save the area (number of contributing cells)
            sp.A_endo = cnt_fifo
