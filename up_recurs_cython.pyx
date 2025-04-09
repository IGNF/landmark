# up_recurs_cython.pyx
# Cython implementation of the up_recurs function

from libc.stdlib cimport malloc, free
import numpy as np
cimport numpy as np



cpdef void up_recurs_cython(int curr, object model, np.ndarray mat_val):
    """
    Optimized Cython version of up_recurs for updating downslope path lengths.

    Parameters
    ----------
    curr : int
        ID of the current (parent) channel.
    model : object
        The hydrological model object containing drainage networks and points.
    mat_val : ndarray of int or object
        Matrix to store flow line identifiers temporarily (same shape as mat_id).
    """
    cdef int i, cnt_pt, cnt_in
    cdef int in_curr, i_curr, j_curr, i_out, j_out
    cdef double l1, l2, l3
    cdef object dp, out_dp

    net_parent = model.dr_net[curr - 1]
    for i in range(net_parent.n_jun):
        in_curr = net_parent.id_in.value[i - 1].value
        net_child = model.dr_net[in_curr - 1]

        for cnt_pt in range(net_child.nel - 1):
            dp = model.dr_pt[net_child.id_pnts.value[cnt_pt] - 1]
            l1 = net_child.length
            l2 = dp.upl
            l3 = model.dr_pt[net_child.id_end_pt.value - 1].dpl
            dp.dpl = l1 - l2 + l3

            i_curr = dp.i
            j_curr = dp.j

            if dp.fldir.value is not None and dp.fldir.value > 0:
                out_dp = model.dr_pt[dp.fldir.value - 1]
                i_out = out_dp.i
                j_out = out_dp.j
                i_mat = i_curr * 2 + (i_out - i_curr)
                j_mat = j_curr * 2 + (j_out - j_curr)
                mat_val[i_mat, j_mat] = net_child.id_ch.value

        net_child.n_path = net_parent.n_path + 1
        net_child.id_path.append(net_child.id_ch.value)

        for cnt_in in range(2, net_child.n_path + 1):
            net_child.id_path.append(net_parent.id_path[cnt_in - 2])

        net_child.id_endo = net_parent.id_endo

        up_recurs_cython(in_curr, model, mat_val)
