# -*- coding: utf-8 -*-
"""
Transcript in python of the original dpl.f90

! The subroutine Downslope Path Legth calculates the legth of the path 
! between each DTM cell and the outflow point even if the basin is 
! endorheic 
"""

import numpy as np
from tqdm import tqdm

class DPL:
    def __init__(self, loader):
        """
        Initializes the DPL class with the data from the loader object.
        Parameters
        ----------
        loader : LoadData
            An instance of the LoadData class containing the drainage points and grid.
        """
        self.dr_pt = loader.dr_pt  # List of drainage points
        self.mat_id = loader.mat_id  # 2D grid of drainage points
        self.delta_x = loader.delta_x
        self.delta_y = loader.delta_y

    def calculate_downslope_length(self):
        """
        Calculates the downslope path length for each drainage point.
        """
        for point in tqdm(self.dr_pt, desc="Calculating downslope path length", unit="point"):
            self.compute_path_length(point)

    def compute_path_length(self, point):
        """
        Recursively calculates the downslope path length for a given point.
        Parameters
        ----------
        point : dict
            The drainage point for which to calculate the downslope path length.
        """
        if point['dpl'] > 0:
            # Path length already computed
            return point['dpl']

        if point['fldir'] is None:
            # No flow direction, end of the path
            point['dpl'] = 0
            return 0

        # Get the downstream point
        downstream_point = next((p for p in self.dr_pt if p['id_pnt'] == point['fldir']), None)
        if downstream_point is None:
            point['dpl'] = 0
            return 0

        # Recursive calculation
        distance = self.compute_distance(point, downstream_point)
        point['dpl'] = distance + self.compute_path_length(downstream_point)
        return point['dpl']

    def compute_distance(self, point, downstream_point):
        """
        Computes the distance between two points.
        Parameters
        ----------
        point : dict
            The starting drainage point.
        downstream_point : dict
            The downstream drainage point.

        Returns
        -------
        float
            The distance between the two points.
        """
        dx = abs(point['i'] - downstream_point['i']) * self.delta_x
        dy = abs(point['j'] - downstream_point['j']) * self.delta_y
        return np.sqrt(dx**2 + dy**2)

# Example usage
if __name__ == "__main__":
    from load_data_translation import LoadData
    from slopeline_translation import SlopeLine

    # Load the data
    loader = LoadData()
    loader.process('header.dat', 'dtm.dat')

    # Perform slopeline calculation
    slopeline = SlopeLine(loader)
    slopeline.calculate_slopeline()

    # Perform downslope path length calculation
    dpl = DPL(loader)
    dpl.calculate_downslope_length()

    print("Drainage points after downslope path length calculation:")
    for point in loader.dr_pt[:5]:
        print(point)
