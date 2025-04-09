# -*- coding: utf-8 -*-
"""
D8_LTD algoritm from the original D8_LTD.f90 Fortran code
"""

import numpy as np
from math import sin, sqrt, pi, atan2
from typing import Tuple, Optional

from hydro_utils_cython import facet

# A small epsilon for floating point comparisons.
EPSILON = np.finfo(np.float32).eps 

# def facet(e0: float, e1: float, e2: float, delta_x: float, delta_y: float) -> Tuple[float, float]:
#     """
#     Calculate the aspect and maximum slope for a triangular facet.

#     Parameters
#     ----------
#     e0 : float
#         Elevation at the central point.
#     e1 : float
#         Elevation at the first neighbor point.
#     e2 : float
#         Elevation at the second neighbor point.
#     delta_x : float
#         Grid spacing in the x-direction.
#     delta_y : float
#         Grid spacing in the y-direction.

#     Returns
#     -------
#     tuple of float
#         r : float
#             Computed aspect in radians.
#         s_max_facet : float
#             Maximum slope of the facet (positive downward).
#     """    
#     s1 = (e0 - e1) / delta_x
#     s2 = (e1 - e2) / delta_x
#     if abs(s1) < EPSILON:
#         if s2 >= 0.0:
#             r_val = pi / 2.0
#         else:
#             r_val = -pi / 2.0
#     else:
#         r_val = atan2(s2, s1)
#     sp = sqrt(s1**2 + s2**2)
#     sd = (e0 - e2) / sqrt(delta_x**2 + delta_y**2)
#     if (r_val >= 0.0 and r_val <= pi/4.0 and s1 >= 0.0):
#         s_max_facet = sp
#     else:
#         if s1 > sd:
#             s_max_facet = s1
#             r_val = 0.0
#         else:
#             s_max_facet = sd
#             r_val = pi/4.0
#     return r_val, s_max_facet



class SlopelineMixin:
    def d8_ltd(self, dp) -> Tuple[Optional[int], Optional[int], int, float]:
        """Compute the D8 flow direction and cumulative deviation from a drainage point.

        Parameters
        ----------
        dp : DrainagePoint
            The drainage point for which to compute the direction. Should include:
            - i, j : int
                Grid coordinates in the DEM.
            - sumdev : float
                Initial cumulative deviation.

        Returns
        -------
        tuple
            i_out : int or None
                Row index of selected outflow cell.
            j_out : int or None
                Column index of selected outflow cell.
            ndfl : int
                Validity flag (1 if all data valid, 0 otherwise).
            sumdev : float
                Updated cumulative deviation for the selected flow path.
        """
        # Get current point coordinates and initial cumulative deviation.
        i = dp.i
        j = dp.j
        sumdev_in = dp.sumdev  # This is the upstream (or initial) cumulative deviation.

        
        # Initialize outputs.
        ndfl = 1
        i_out = None
        j_out = None
        sumdev = 0
        
        
        # Prepare an array (list) to store the elevations in a 3×3 window around (i, j).
        e = [0.0] * 9
        l = 0
        # Loop over neighbors: For jj = j-1 to j+1; for ii = i+1 down to i-1.
        for jj in range(j - 1, j + 2):
            for ii in range(i + 1, i - 2, -1):  # Descending: i+1, i, i-1.
                l += 1
                # Boundary check: valid drainage point indices are 0 <= index < self.N (or self.M).
                if ii < 0 or ii >= self.N:
                    ndfl = 0
                    continue
                if jj < 0 or jj >= self.M:
                    ndfl = 0
                    continue
                # Retrieve the neighbor from self.mat_id.
                # In our reading process, a drainage point for cell (ii, jj) is stored at index (ii*2, jj*2).
                neighbor = self.mat_id[ii * 2, jj * 2]
                if neighbor is None:
                    ndfl = 0
                    continue
                curr_elev = neighbor.Z
                if curr_elev <= self.nodata:
                    ndfl = 0
                    continue
                e[l - 1] = curr_elev
        
        # The center cell's elevation.
        e0 = e[4]
        
        # print('e = ', e)
        
        # Prepare lists for the triangle (facet) calculations.
        num_triangles = 9  # We define 8 triangles (facets) around the center + center.
        e1_fmax = [0.0] * num_triangles
        e2_fmax = [0.0] * num_triangles
        r_max = [0.0] * num_triangles
        s_max = [self.nodata] * num_triangles  # Initialize slopes with nodata.
        i_out1_arr = [0] * num_triangles
        j_out1_arr = [0] * num_triangles
        i_out2_arr = [0] * num_triangles
        j_out2_arr = [0] * num_triangles
        sigma_arr = [0.0] * num_triangles
        
        # For each triangle, we use specific pairs of neighbor elevations.
        # Triangle 021: uses e[1] and e[0] (Fortran indices 2 and 1).
        if abs(e[1] * e[0]) > EPSILON:
            r_val, s_max_facet = facet(e0, e[1], e[0], self.delta_x, self.delta_y)
            e1_fmax[0] = e[1]
            e2_fmax[0] = e[0]
            r_max[0] = r_val
            s_max[0] = s_max_facet
            i_out1_arr[0] = i       # Candidate outflow: same row.
            j_out1_arr[0] = j - 1   # and column one less.
            i_out2_arr[0] = i + 1   # Alternative: next row.
            j_out2_arr[0] = j - 1
            sigma_arr[0] = 1.0
        # Triangle 023: uses e[1] and e[2] (Fortran indices 2 and 3).
        if abs(e[1] * e[2]) > EPSILON:
            r_val, s_max_facet = facet(e0, e[1], e[2], self.delta_x, self.delta_y)
            e1_fmax[1] = e[1]
            e2_fmax[1] = e[2]
            r_max[1] = r_val
            s_max[1] = s_max_facet
            i_out1_arr[1] = i
            j_out1_arr[1] = j - 1
            i_out2_arr[1] = i - 1
            j_out2_arr[1] = j - 1
            sigma_arr[1] = -1.0
        # Triangle 063: uses e[5] and e[2] (Fortran indices 6 and 3).
        if abs(e[5] * e[2]) > EPSILON:
            r_val, s_max_facet = facet(e0, e[5], e[2], self.delta_x, self.delta_y)
            e1_fmax[2] = e[5]
            e2_fmax[2] = e[2]
            r_max[2] = r_val
            s_max[2] = s_max_facet
            i_out1_arr[2] = i - 1
            j_out1_arr[2] = j
            i_out2_arr[2] = i - 1
            j_out2_arr[2] = j - 1
            sigma_arr[2] = 1.0
        # Triangle 069: uses e[5] and e[8] (Fortran indices 6 and 9).
        if abs(e[5] * e[8]) > EPSILON:
            r_val, s_max_facet = facet(e0, e[5], e[8], self.delta_x, self.delta_y)
            e1_fmax[3] = e[5]
            e2_fmax[3] = e[8]
            r_max[3] = r_val
            s_max[3] = s_max_facet
            i_out1_arr[3] = i - 1
            j_out1_arr[3] = j
            i_out2_arr[3] = i - 1
            j_out2_arr[3] = j + 1
            sigma_arr[3] = -1.0
        # Triangle 089: uses e[7] and e[8] (Fortran indices 8 and 9).
        if abs(e[7] * e[8]) > EPSILON:
            r_val, s_max_facet = facet(e0, e[7], e[8], self.delta_x, self.delta_y)
            e1_fmax[5] = e[7]
            e2_fmax[5] = e[8]
            r_max[5] = r_val
            s_max[5] = s_max_facet
            i_out1_arr[5] = i
            j_out1_arr[5] = j + 1
            i_out2_arr[5] = i - 1
            j_out2_arr[5] = j + 1
            sigma_arr[5] = 1.0
        # Triangle 087: uses e[7] and e[6] (Fortran indices 8 and 7).
        if abs(e[7] * e[6]) > EPSILON:
            r_val, s_max_facet = facet(e0, e[7], e[6], self.delta_x, self.delta_y)
            # if dp.id_pnt.value == 3474 :
            #     print('e0 = ', e0)
            #     print('e1 = ', e[7])
            #     print('e2 = ', e[6])
                
            e1_fmax[6] = e[7]
            e2_fmax[6] = e[6]
            r_max[6] = r_val
            s_max[6] = s_max_facet
            i_out1_arr[6] = i
            j_out1_arr[6] = j + 1
            i_out2_arr[6] = i + 1
            j_out2_arr[6] = j + 1
            sigma_arr[6] = -1.0
        # Triangle 047: uses e[3] and e[6] (Fortran indices 4 and 7).
        if abs(e[3] * e[6]) > EPSILON:
            r_val, s_max_facet = facet(e0, e[3], e[6], self.delta_x, self.delta_y)
            e1_fmax[7] = e[3]
            e2_fmax[7] = e[6]
            r_max[7] = r_val
            s_max[7] = s_max_facet
            i_out1_arr[7] = i + 1
            j_out1_arr[7] = j
            i_out2_arr[7] = i + 1
            j_out2_arr[7] = j + 1
            sigma_arr[7] = 1.0
        # Triangle 041: uses e[3] and e[0] (Fortran indices 4 and 1).
        if abs(e[3] * e[0]) > EPSILON:
            r_val, s_max_facet = facet(e0, e[3], e[0], self.delta_x, self.delta_y)
            e1_fmax[8] = e[3]
            e2_fmax[8] = e[0]
            r_max[8] = r_val
            s_max[8] = s_max_facet
            i_out1_arr[8] = i + 1
            j_out1_arr[8] = j
            i_out2_arr[8] = i + 1
            j_out2_arr[8] = j - 1
            sigma_arr[8] = -1.0
            
        
        # Select the triangle (facet) with the maximum slope.
        s_mx = max(s_max)
        # id_mx = s_max.index(s_mx)
        id_mx = np.argmax(s_max)
        e1_fmx_val = e1_fmax[id_mx]
        e2_fmx_val = e2_fmax[id_mx]
        r_mx = r_max[id_mx]
        i_out1_mx = i_out1_arr[id_mx]
        j_out1_mx = j_out1_arr[id_mx]
        i_out2_mx = i_out2_arr[id_mx]
        j_out2_mx = j_out2_arr[id_mx]
        sigma_mx = sigma_arr[id_mx]
                    
        
        # If the maximum slope is positive, compute deviations and decide the final output.
        if s_mx > 0.0:
            dev_1 = self.delta_x * sin(r_mx)
            dev_2 = self.delta_x * sqrt(2.0) * sin(pi/4.0 - r_mx)
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
