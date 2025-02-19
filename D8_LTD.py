# -*- coding: utf-8 -*-
"""
D8_LTD algoritm from the original D8_LTD.f90 Fortran code
"""

import numpy as np
from math import sin, atan, sqrt, pi

# A small epsilon for floating point comparisons.
EPSILON = 1e-6

def facet(e0, e1, e2, delta_x, delta_y):
    """
    Calculates the aspect (r) and the maximum slope (s_max_facet)
    within a given triangle (facet), following the Fortran FACET subroutine.
    
    Parameters:
      e0, e1, e2 : float
          Elevations at the three vertices.
      delta_x, delta_y : float
          Grid spacings.
    
    Returns:
      (r, s_max_facet) : tuple of floats
          r: computed aspect (in radians),
          s_max_facet: the slope (positive downward) of the facet.
    """
    s1 = (e0 - e1) / delta_x
    s2 = (e1 - e2) / delta_y
    if abs(s1) < EPSILON:
        if s2 >= 0.0:
            r_val = pi / 2.0
        else:
            r_val = -pi / 2.0
    else:
        r_val = atan(s2 / s1)
    sp = sqrt(s1**2 + s2**2)
    sd = (e0 - e2) / sqrt(delta_x**2 + delta_y**2)
    if (r_val >= 0.0 and r_val <= pi/4.0 and s1 >= 0.0):
        s_max_facet = sp
    else:
        if s1 > sd:
            s_max_facet = s1
            r_val = 0.0
        else:
            s_max_facet = sd
            r_val = pi/4.0
    return r_val, s_max_facet

class SlopelineMixin:
    def d8_ltd(self, dp):
        """
        Computes the drainage direction for a given drainage point using a method
        similar to the Fortran D8_LTD subroutine.
        
        Parameters:
          dp : DrainagePoint
              The drainage point (object) for which to compute the direction.
              It is assumed that dp has attributes:
                 - i, j: the row and column indices (0-indexed) in the DEM,
                 - sumdev: the cumulative deviation (initially set, e.g., to 0),
              And that the following attributes are available in self:
                 - self.N, self.M: dimensions of the DEM (number of drainage points in each direction),
                 - self.nodata: the no‑data value,
                 - self.delta_x, self.delta_y: grid spacings,
                 - self.mat_id: a 2D numpy array of objects (DrainagePoint) 
                   stored at positions (i*2, j*2) for each cell.
        
        Returns:
          (i_out, j_out, ndfl, sumdev) where:
            - i_out, j_out: the indices (0-indexed) of the selected outflow neighbor,
            - ndfl: integer flag (1 if valid, 0 if a boundary or invalid cell was encountered),
            - sumdev: the cumulative deviation computed for this point.
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
        
        if dp.id_pnt.value == 2802:
            print("=" * 40)
            print(f"DEBUG - Values for id_dr = {dp.id_pnt.value}")
            
            print("s_max values:")
            print(s_max)
            
            print("r_max values:")
            print(r_max)
            
            print(f"id_mx = {id_mx}")
            print(f"s_mx = {s_mx}")
            
            print(f"e1_fmx = {e1_fmx_val}, e2_fmx = {e2_fmx_val}")
            print(f"r_mx = {r_mx}")
            
            print(f"i_out1_mx = {i_out1_mx}, j_out1_mx = {j_out1_mx}")
            print(f"i_out2_mx = {i_out2_mx}, j_out2_mx = {j_out2_mx}")
            
            print(f"sigma_mx = {sigma_mx}")
            print("=" * 40)
            
        
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
        
            if dp.id_pnt.value == 2802:
                print("=" * 25)
                print(f"DEBUG - Computed values for id_pnt = {dp.id_pnt.value}")
                print(f"sumdev_in = {sumdev_in:.6f}")
                print(f"dev_1 = {dev_1:.6f}, dev_2 = {dev_2:.6f}")
                print(f"sumdev_1 = {sumdev_1:.6f}, sumdev_2 = {sumdev_2:.6f}")
                print(f"i_out = {i_out}, j_out = {j_out}")
                print("=" * 25)
        return i_out, j_out, ndfl, sumdev
