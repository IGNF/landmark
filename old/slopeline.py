# -*- coding: utf-8 -*-
"""
Module for calculating drainage directions using the D8-LTD method and defining slopelines.
This module builds on the previously loaded DEM and header data.
"""

#Intern imports
import numpy as np
from math import sqrt, sin, atan, pi
import time
from tqdm import tqdm

#Extern imports
from data_structures import DrainageNetwork, EndoPoint


# Small epsilon for floating point comparisons
EPSILON = np.finfo(float).eps


def facet(e0, e1, e2, delta_x, delta_y):
    """
    Calculates the aspect (r) and the maximum slope (s_max_facet)
    within a given triangle (facet), following the Fortran FACET subroutine.
    
    Parameters:
      e0, e1, e2 : float
         Elevations at the three vertices.
      delta_x, delta_y : float
         Grid spacings (assumed equal in this case).
    
    Returns:
      (r, s_max_facet) : tuple of float
         r: computed aspect (in radians),
         s_max_facet: the slope (positive downward).
    """
    # Use pi from math
    pi_val = pi
    s1 = (e0 - e1) / delta_x
    s2 = (e1 - e2) / delta_x
    if abs(s1) < EPSILON:
        if s2 >= 0.0:
            r = pi_val / 2.0
        else:
            r = -pi_val / 2.0
    else:
        r = atan(s2 / s1)
    sp = sqrt(s1**2 + s2**2)
    sd = (e0 - e2) / sqrt(delta_x**2 + delta_y**2)
    if (r >= 0.0 and r <= pi_val/4.0 and s1 >= 0.0):
        s_max_facet = sp
    else:
        if s1 > sd:
            s_max_facet = s1
            r = 0.0
        else:
            s_max_facet = sd
            r = pi_val / 4.0
    return r, s_max_facet

class SlopelineMixin:
    def d8_ltd(self, dp):
        """
        Implements the D8-LTD method to determine the drainage direction for a point.
        
        Parameters
        ----------
        dp : DrainagePoint
            The drainage point for which to compute the outflow.
        
        Returns
        -------
        (i_out, j_out, ndfl, sumdev) : tuple
            i_out, j_out : grid indices (in drainage point coordinates) of the selected outflow point.
            ndfl : integer flag (1 = valid, 0 = boundary/outflow)
            sumdev : cumulative deviation (updated)
        """
        # Initialize outputs
        ndfl = 1
        sumdev = 0.0
        i_out = 0
        j_out = 0

        # Current point indices (0-indexed)
        i = dp.i
        j = dp.j
        sumdev_in = dp.sumdev  # initial cumulative deviation

        # Allocate an array to store elevations of the 9 cells (neighbors including center)
        e = np.zeros(9)
        l = 0

        # Loop over neighbors:
        # In Fortran, the loop is: do jj = j-1, j+1 and do ii = i+1, i-1, -1.
        # In Python, we use a nested loop with ii descending from i+1 down to i-1.
        for jj in range(j - 1, j + 2):
            for ii in range(i + 1, i - 2, -1):  # descending: i+1, i, i-1
                l += 1
                # Boundary check (Python indices: valid if 0 <= index < N or M)
                if ii < 0 or ii >= self.N or jj < 0 or jj >= self.M:
                    ndfl = 0
                    continue  # Skip this neighbor
                # Retrieve the neighbor from mat_id (points stored at (ii*2, jj*2))
                neighbor = self.mat_id[ii * 2, jj * 2]
                if neighbor is None:
                    ndfl = 0
                    continue
                curr_elev = neighbor.Z
                if curr_elev <= self.nodata:
                    ndfl = 0
                    continue
                e[l - 1] = curr_elev
        # e[4] corresponds to the center cell (Fortran e(5))
        e0 = e[4]

        # Process 8 triangles (facets) to compute maximum slope on each facet.
        num_triangles = 8
        e1_fmax = np.zeros(num_triangles)
        e2_fmax = np.zeros(num_triangles)
        r_max_arr = np.zeros(num_triangles)
        s_max_arr = np.zeros(num_triangles)
        i_out1_arr = np.zeros(num_triangles, dtype=int)
        j_out1_arr = np.zeros(num_triangles, dtype=int)
        i_out2_arr = np.zeros(num_triangles, dtype=int)
        j_out2_arr = np.zeros(num_triangles, dtype=int)
        sigma_arr = np.zeros(num_triangles)

        # Triangle 021 (using Fortran indices: e(2) & e(1) → Python: e[1] & e[0])
        e1 = e[1]
        e2 = e[0]
        if abs(e1 * e2) > EPSILON:
            r_val, s_max_facet = facet(e0, e1, e2, self.delta_x, self.delta_y)
            e1_fmax[0] = e1
            e2_fmax[0] = e2
            r_max_arr[0] = r_val
            s_max_arr[0] = s_max_facet
            i_out1_arr[0] = i
            j_out1_arr[0] = j - 1
            i_out2_arr[0] = i + 1
            j_out2_arr[0] = j - 1
            sigma_arr[0] = 1.0

        # Triangle 023 (e(2) & e(3) → Python: e[1] & e[2])
        e1 = e[1]
        e2 = e[2]
        if abs(e1 * e2) > EPSILON:
            r_val, s_max_facet = facet(e0, e1, e2, self.delta_x, self.delta_y)
            e1_fmax[1] = e1
            e2_fmax[1] = e2
            r_max_arr[1] = r_val
            s_max_arr[1] = s_max_facet
            i_out1_arr[1] = i
            j_out1_arr[1] = j - 1
            i_out2_arr[1] = i - 1
            j_out2_arr[1] = j - 1
            sigma_arr[1] = -1.0

        # Triangle 063 (e(6) & e(3) → Python: e[5] & e[2])
        e1 = e[5]
        e2 = e[2]
        if abs(e1 * e2) > EPSILON:
            r_val, s_max_facet = facet(e0, e1, e2, self.delta_x, self.delta_y)
            e1_fmax[2] = e1
            e2_fmax[2] = e2
            r_max_arr[2] = r_val
            s_max_arr[2] = s_max_facet
            i_out1_arr[2] = i - 1
            j_out1_arr[2] = j
            i_out2_arr[2] = i - 1
            j_out2_arr[2] = j - 1
            sigma_arr[2] = 1.0

        # Triangle 069 (e(6) & e(9) → Python: e[5] & e[8])
        e1 = e[5]
        e2 = e[8]
        if abs(e1 * e2) > EPSILON:
            r_val, s_max_facet = facet(e0, e1, e2, self.delta_x, self.delta_y)
            e1_fmax[3] = e1
            e2_fmax[3] = e2
            r_max_arr[3] = r_val
            s_max_arr[3] = s_max_facet
            i_out1_arr[3] = i - 1
            j_out1_arr[3] = j
            i_out2_arr[3] = i - 1
            j_out2_arr[3] = j + 1
            sigma_arr[3] = -1.0

        # Triangle 089 (e(8) & e(9) → Python: e[7] & e[8])
        e1 = e[7]
        e2 = e[8]
        if abs(e1 * e2) > EPSILON:
            r_val, s_max_facet = facet(e0, e1, e2, self.delta_x, self.delta_y)
            e1_fmax[4] = e1
            e2_fmax[4] = e2
            r_max_arr[4] = r_val
            s_max_arr[4] = s_max_facet
            i_out1_arr[4] = i
            j_out1_arr[4] = j + 1
            i_out2_arr[4] = i - 1
            j_out2_arr[4] = j + 1
            sigma_arr[4] = 1.0

        # Triangle 087 (e(8) & e(7) → Python: e[7] & e[6])
        e1 = e[7]
        e2 = e[6]
        if abs(e1 * e2) > EPSILON:
            r_val, s_max_facet = facet(e0, e1, e2, self.delta_x, self.delta_y)
            e1_fmax[5] = e1
            e2_fmax[5] = e2
            r_max_arr[5] = r_val
            s_max_arr[5] = s_max_facet
            i_out1_arr[5] = i
            j_out1_arr[5] = j + 1
            i_out2_arr[5] = i + 1
            j_out2_arr[5] = j + 1
            sigma_arr[5] = -1.0

        # Triangle 047 (e(4) & e(7) → Python: e[3] & e[6])
        e1 = e[3]
        e2 = e[6]
        if abs(e1 * e2) > EPSILON:
            r_val, s_max_facet = facet(e0, e1, e2, self.delta_x, self.delta_y)
            e1_fmax[6] = e1
            e2_fmax[6] = e2
            r_max_arr[6] = r_val
            s_max_arr[6] = s_max_facet
            i_out1_arr[6] = i + 1
            j_out1_arr[6] = j
            i_out2_arr[6] = i + 1
            j_out2_arr[6] = j + 1
            sigma_arr[6] = 1.0

        # Triangle 041 (e(4) & e(1) → Python: e[3] & e[0])
        e1 = e[3]
        e2 = e[0]
        if abs(e1 * e2) > EPSILON:
            r_val, s_max_facet = facet(e0, e1, e2, self.delta_x, self.delta_y)
            e1_fmax[7] = e1
            e2_fmax[7] = e2
            r_max_arr[7] = r_val
            s_max_arr[7] = s_max_facet
            i_out1_arr[7] = i + 1
            j_out1_arr[7] = j
            i_out2_arr[7] = i + 1
            j_out2_arr[7] = j - 1
            sigma_arr[7] = -1.0

        # Select the triangle with the maximum slope
        s_mx = np.max(s_max_arr)
        id_mx = np.argmax(s_max_arr)
        e1_fmx_val = e1_fmax[id_mx]
        e2_fmx_val = e2_fmax[id_mx]
        r_mx = r_max_arr[id_mx]
        i_out1_mx = i_out1_arr[id_mx]
        j_out1_mx = j_out1_arr[id_mx]
        i_out2_mx = i_out2_arr[id_mx]
        j_out2_mx = j_out2_arr[id_mx]
        sigma_mx = sigma_arr[id_mx]

        if s_mx > 0.0:
            dev_1 = self.delta_x * sin(r_mx)
            dev_2 = self.delta_x * np.sqrt(2.0) * sin(pi/4.0 - r_mx)
            if sigma_mx == 1:
                dev_2 = -dev_2
            else:
                dev_1 = -dev_1
            sumdev_1 = sumdev_in + dev_1
            sumdev_2 = sumdev_in + dev_2
            a1 = abs(sumdev_1)
            a2 = abs(sumdev_2)
            if abs(a1 - a2) < EPSILON and e0 > e1_fmx_val:
                sumdev = sumdev_1
                i_out = i_out1_mx
                j_out = j_out1_mx
            else:
                if a1 < a2 and e0 > e1_fmx_val:
                    sumdev = sumdev_1
                    i_out = i_out1_mx
                    j_out = j_out1_mx
                elif a1 > a2 or e0 > e2_fmx_val:
                    sumdev = sumdev_2
                    i_out = i_out2_mx
                    j_out = j_out2_mx
        return i_out, j_out, ndfl, sumdev

    def calculate_slopelines(self):
        """
        Calculates drainage directions and defines slopelines using a D8-LTD approach.
        This method iterates over the drainage points in descending order of elevation.
        It updates each drainage point with its outflow (fdir) and builds drainage networks
        (channels) as well as endorheic point lists.
        """
        start_time = time.process_time()

        # Initialize drainage network and endorheic lists and counters
        self.dr_net = []
        self.dr_net_by_id = {}
        self.endo_pt = []
        self.channel_count = 0
        self.endorheic_count = 0

        # Iterate over drainage points in descending order (highest first)
        for dp in tqdm(reversed(self.dr_pt), desc="Calculating slopelines", unit="point"):

            if dp.ninf > 0:
                # Process channel points (with inflows)
                if dp.Linflow:
                    max_val = max(dp.Linflow)
                    max_index = dp.Linflow.index(max_val)
                    dp.upl = max_val
                    # Retrieve the inflow point corresponding to maximum upstream length
                    id_up_max = dp.inflow[max_index]
                    up_dp = next((p for p in self.dr_pt if p.id_pnt == id_up_max), None)
                    if up_dp is not None:
                        dp.id_ch = up_dp.id_ch
                        id_main = dp.id_ch
                    else:
                        dp.id_ch = None
                        id_main = None
                for inflow_id in dp.inflow:
                    inflow_dp = next((p for p in self.dr_pt if p.id_pnt == inflow_id), None)
                    if inflow_dp is None:
                        continue
                    curr_ch = inflow_dp.id_ch
                    if curr_ch not in self.dr_net_by_id:
                        net = DrainageNetwork(
                            id_ch=curr_ch,
                            id_pnts=[inflow_dp.id_pnt],
                            id_start_pt=inflow_dp.id_pnt,
                            id_end_pt=inflow_dp.id_pnt,
                            length=0.0,
                            id_ch_out=curr_ch,
                            n_jun=0,
                            id_in=[]
                        )
                        self.dr_net.append(net)
                        self.dr_net_by_id[curr_ch] = net
                    net = self.dr_net_by_id[curr_ch]
                    net.id_pnts.append(dp.id_pnt)
                    if dp.Linflow:
                        net.length = dp.Linflow[dp.inflow.index(inflow_id)]
                    net.id_ch_out = dp.id_ch
                    net.id_end_pt = dp.id_pnt
                    if id_main is not None and curr_ch != id_main:
                        main_net = self.dr_net_by_id.get(id_main)
                        if main_net is not None:
                            main_net.id_in.append(curr_ch)
            else:
                # Process channel heads (no inflows)
                self.channel_count += 1
                dp.id_ch = self.channel_count
                net = DrainageNetwork(
                    id_ch=self.channel_count,
                    id_pnts=[dp.id_pnt],
                    id_start_pt=dp.id_pnt,
                    id_end_pt=dp.id_pnt,
                    length=0.0,
                    id_ch_out=self.channel_count,
                    n_jun=0,
                    id_in=[]
                )
                self.dr_net.append(net)
                self.dr_net_by_id[self.channel_count] = net

            # Call the updated D8-LTD function to determine the outflow for the current point.
            i_out, j_out, ndfl, sumdev = self.d8_ltd(dp)
            if i_out is not None and j_out is not None and i_out != 0 and j_out != 0:
                out_dp = self.mat_id[i_out * 2, j_out * 2]
                if out_dp is not None:
                    dp.fdir = out_dp.id_pnt
                    out_dp.ninf += 1
                    out_dp.inflow.append(dp.id_pnt)
                    distance = sqrt(((dp.i - out_dp.i) * self.delta_x)**2 +
                                    ((dp.j - out_dp.j) * self.delta_y)**2)
                    out_dp.Linflow.append(dp.upl + distance)
                    out_dp.Sinflow.append(sumdev)
            else:
                # Endorheic point processing
                self.endorheic_count += 1
                dp.id_endo = self.endorheic_count
                new_endo = EndoPoint(
                    id_eo=self.endorheic_count,
                    id_pnt=dp.id_pnt,
                    bas_type=ndfl,
                    nsaddle=0
                )
                self.endo_pt.append(new_endo)

        elapsed_time = time.process_time() - start_time
        hours = int(elapsed_time // 3600)
        minutes = int((elapsed_time % 3600) // 60)
        seconds = elapsed_time % 60
        print(f"Elapsed time: {hours}h {minutes}m {seconds:.2f}s")

