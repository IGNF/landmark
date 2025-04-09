from libc.stdlib cimport malloc, free
import numpy as np
cimport numpy as np


cdef extern from "math.h":
    double sin(double x)
    double sqrt(double x)
    double fabs(double x)
    double atan2(double y, double x)
    double M_PI

cdef double EPSILON = 1e-7  # Approx equivalent to np.finfo(np.float32).eps

cpdef tuple facet(double e0, double e1, double e2, double delta_x, double delta_y):
    """
    Cython-accelerated version of the facet function.

    Parameters
    ----------
    e0 : float
        Elevation at the central point.
    e1 : float
        Elevation at the first neighbor point.
    e2 : float
        Elevation at the second neighbor point.
    delta_x : float
        Grid spacing in the x-direction.
    delta_y : float
        Grid spacing in the y-direction.

    Returns
    -------
    tuple of float
        r : float
            Computed aspect in radians.
        s_max_facet : float
            Maximum slope of the facet (positive downward).
    """
    cdef double s1, s2, r_val, sp, sd, s_max_facet

    s1 = (e0 - e1) / delta_x
    s2 = (e1 - e2) / delta_x

    if fabs(s1) < EPSILON:
        if s2 >= 0.0:
            r_val = M_PI / 2.0
        else:
            r_val = -M_PI / 2.0
    else:
        r_val = atan2(s2, s1)

    sp = sqrt(s1 * s1 + s2 * s2)
    sd = (e0 - e2) / sqrt(delta_x * delta_x + delta_y * delta_y)

    if 0.0 <= r_val <= M_PI / 4.0 and s1 >= 0.0:
        s_max_facet = sp
    elif s1 > sd:
        s_max_facet = s1
        r_val = 0.0
    else:
        s_max_facet = sd
        r_val = M_PI / 4.0

    return r_val, s_max_facet



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



cpdef void up_recurs_ss_cython(
    int curr,
    object dr_net,
    object dr_pt,
    object dr_pt_in
):
    """
    Cython-optimized version of up_recurs_ss to compute downslope path lengths
    in the subsurface drainage graph.

    Parameters
    ----------
    curr : int
        Current drainage network ID (1-based).
    dr_net : list of DrainageNetwork
        Drainage networks containing lengths, paths, etc.
    dr_pt : list of DrainagePoint
        Drainage points with upl, dpl.
    dr_pt_in : list of DrainagePointInflow
        Tracks inflows and their multiplicity.
    """
    cdef int cnt_pt, id_dr, end_el, in_curr
    cdef float l1, l2, l3

    end_el = dr_net[curr - 1].nel

    for cnt_pt in range(end_el - 1, 0, -1):
        id_dr = dr_net[curr - 1].id_pnts.value[cnt_pt - 1]
        l1 = dr_net[curr - 1].length
        l2 = dr_pt[id_dr - 1].upl
        l3 = dr_pt[dr_net[curr - 1].id_pnts.value[end_el - 1] - 1].dpl

        dr_pt[id_dr - 1].dpl = l1 - l2 + l3

        if dr_pt_in[id_dr - 1].ninf > 1:
            for in_curr in dr_pt_in[id_dr - 1].inflow:
                if in_curr != curr:
                    dr_net[in_curr - 1].n_path = dr_net[curr - 1].n_path + 1
                    dr_net[in_curr - 1].id_path = [dr_net[in_curr - 1].id_ch.value]
                    dr_net[in_curr - 1].id_path.extend(
                        dr_net[curr - 1].id_path[:dr_net[in_curr - 1].n_path - 1]
                    )
                    up_recurs_ss_cython(in_curr, dr_net, dr_pt, dr_pt_in)


# d8_ltd_cython.pyx
# Cython implementation of the D8-LTD flow direction algorithm

from libc.math cimport sin, sqrt, atan2, M_PI
cdef double pi = M_PI


cpdef tuple d8_ltd_cython(
    int i, int j, float sumdev_in,
    object[:, :] mat_id, int N, int M,
    float delta_x, float delta_y,
    float nodata
):
    """
    Compute D8-LTD flow direction and cumulative deviation for a drainage point.

    Parameters
    ----------
    i, j : int
        Coordinates of the drainage point.
    sumdev_in : float
        Cumulative deviation upstream.
    mat_id : 2D object array
        Grid containing DrainagePoint objects.
    N, M : int
        Grid size.
    delta_x, delta_y : float
        Grid resolution.
    nodata : float
        Nodata threshold.

    Returns
    -------
    (i_out, j_out, ndfl, sumdev) : tuple
        Outgoing flow indices, validity flag, and updated deviation.
    """
    cdef int l = 0
    cdef int ndfl = 1
    cdef int i_out = -1, j_out = -1
    cdef float sumdev = 0.0
    cdef float e[9]
    cdef int ii, jj, idx
    cdef object neighbor
    cdef float r_mx, sigma, dev_1, dev_2 


    # Load 3x3 elevation window
    for jj in range(j - 1, j + 2):
        for ii in range(i + 1, i - 2, -1):
            if ii < 0 or ii >= N or jj < 0 or jj >= M:
                ndfl = 0
                e[l] = 0.0
            else:
                neighbor = mat_id[ii * 2, jj * 2]
                if neighbor is None or neighbor.Z <= nodata:
                    ndfl = 0
                    e[l] = 0.0
                else:
                    e[l] = neighbor.Z
            l += 1

    cdef float e0 = e[4]  # Center cell

    # Triangle arrays
    cdef float e1_fmax[9]
    cdef float e2_fmax[9]
    cdef float r_max[9]
    cdef float s_max[9]
    cdef float sigma_arr[9]

    cdef int i_out1_arr[9]
    cdef int j_out1_arr[9]
    cdef int i_out2_arr[9]
    cdef int j_out2_arr[9]


    # Triangle definitions
    cdef int triangles[8][3]
    triangles[0][:] = [1, 0, 0]
    triangles[1][:] = [1, 2, 1]
    triangles[2][:] = [5, 2, 2]
    triangles[3][:] = [5, 8, 3]
    triangles[4][:] = [7, 8, 5]
    triangles[5][:] = [7, 6, 6]
    triangles[6][:] = [3, 6, 7]
    triangles[7][:] = [3, 0, 8]

    for idx in range(8):
        a, b, slot = triangles[idx]
        if abs(e[a] * e[b]) > EPSILON:
            r_val, s_val = facet(e0, e[a], e[b], delta_x, delta_y)
            e1_fmax[slot] = e[a]
            e2_fmax[slot] = e[b]
            r_max[slot] = r_val
            s_max[slot] = s_val

            if slot == 0:
                i_out1_arr[slot], j_out1_arr[slot] = i, j - 1
                i_out2_arr[slot], j_out2_arr[slot] = i + 1, j - 1
                sigma_arr[slot] = 1.0
            elif slot == 1:
                i_out1_arr[slot], j_out1_arr[slot] = i, j - 1
                i_out2_arr[slot], j_out2_arr[slot] = i - 1, j - 1
                sigma_arr[slot] = -1.0
            elif slot == 2:
                i_out1_arr[slot], j_out1_arr[slot] = i - 1, j
                i_out2_arr[slot], j_out2_arr[slot] = i - 1, j - 1
                sigma_arr[slot] = 1.0
            elif slot == 3:
                i_out1_arr[slot], j_out1_arr[slot] = i - 1, j
                i_out2_arr[slot], j_out2_arr[slot] = i - 1, j + 1
                sigma_arr[slot] = -1.0
            elif slot == 5:
                i_out1_arr[slot], j_out1_arr[slot] = i, j + 1
                i_out2_arr[slot], j_out2_arr[slot] = i - 1, j + 1
                sigma_arr[slot] = 1.0
            elif slot == 6:
                i_out1_arr[slot], j_out1_arr[slot] = i, j + 1
                i_out2_arr[slot], j_out2_arr[slot] = i + 1, j + 1
                sigma_arr[slot] = -1.0
            elif slot == 7:
                i_out1_arr[slot], j_out1_arr[slot] = i + 1, j
                i_out2_arr[slot], j_out2_arr[slot] = i + 1, j + 1
                sigma_arr[slot] = 1.0
            elif slot == 8:
                i_out1_arr[slot], j_out1_arr[slot] = i + 1, j
                i_out2_arr[slot], j_out2_arr[slot] = i + 1, j - 1
                sigma_arr[slot] = -1.0
        else:
            s_max[slot] = -1.0

    # Find maximum slope
    cdef float s_mx = -1.0
    cdef int id_mx = -1
    for idx in range(9):
        if s_max[idx] > s_mx:
            s_mx = s_max[idx]
            id_mx = idx

    if s_mx > 0:
        r_mx = r_max[id_mx]
        sigma = sigma_arr[id_mx]
        dev_1 = delta_x * sin(r_mx)
        dev_2 = delta_x * sqrt(2.0) * sin(pi / 4.0 - r_mx)

        if sigma == 1:
            dev_2 = -dev_2
        else:
            dev_1 = -dev_1

        sum1 = sumdev_in + dev_1
        sum2 = sumdev_in + dev_2
        a1 = abs(sum1)
        a2 = abs(sum2)

        if abs(a1 - a2) < EPSILON and e0 > e1_fmax[id_mx]:
            sumdev = sum1
            i_out, j_out = i_out1_arr[id_mx], j_out1_arr[id_mx]
        elif a1 < a2 and e0 > e1_fmax[id_mx]:
            sumdev = sum1
            i_out, j_out = i_out1_arr[id_mx], j_out1_arr[id_mx]
        elif a1 > a2 or e0 > e2_fmax[id_mx]:
            sumdev = sum2
            i_out, j_out = i_out2_arr[id_mx], j_out2_arr[id_mx]

    return i_out, j_out, ndfl, sumdev


