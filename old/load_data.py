# -*- coding: utf-8 -*-
"""
Transcript in python of the original load_data.f90

! The subroutine Load Data loads dtm data and header.
! It sorts also elevations in ascending order.
"""

import numpy as np
from tqdm import tqdm

class LoadData:
    def __init__(self):
        self.delta_x = None
        self.delta_y = None
        self.N = None
        self.M = None
        self.xllcorner = None
        self.yllcorner = None
        self.nodata = None
        self.dr_pt = []
        self.mat_id = None

    def read_header(self, header_file):
        """
        Reads the header file and extracts metadata about the grid.
        Parameters
        ----------
        header_file : str
            Path to the header file (e.g., 'header.dat').
        """
        try:
            with open(header_file, 'r') as f:
                # Skip lines until "Grid" is found
                for line in f:
                    if "Grid" in line:
                        break
                    
                # Start parsing from the current line
                self.delta_x = float(line.split('=')[-1].strip())
                self.delta_y = float(next(f).split('=')[-1].strip())
                self.N = int(next(f).split('=')[-1].strip())
                self.M = int(next(f).split('=')[-1].strip())
                self.xllcorner = float(next(f).split('=')[-1].strip())
                self.yllcorner = float(next(f).split('=')[-1].strip())
                self.nodata = float(next(f).split('=')[-1].strip())
        except Exception as e:
            raise RuntimeError(f"Error reading header file: {e}")

    def load_dtm(self, dtm_file):
        """
        Reads the DTM file and populates drainage points and the grid.
        Parameters
        ----------
        dtm_file : str
            Path to the DTM file (e.g., 'dtm.dat').
        """
        try:
            # Initialize structures
            self.mat_id = np.full((self.N * 2 - 1, self.M * 2 - 1), None)
            with open(dtm_file, 'r') as f:
                for i in tqdm(range(self.N)):
                    row = list(map(float, f.readline().strip().split()))
                    for j, elevation in enumerate(row):
                        if elevation > self.nodata:
                            point = {
                                'i': i + 1,
                                'j': j + 1,
                                'Z': elevation,
                                'id_pnt': len(self.dr_pt) + 1,
                                'fldir': None,
                                'A_in': 0,
                                'upl': 0,
                                'dpl': 0,
                                'id_endo': 0
                            }
                            self.dr_pt.append(point)
                            self.mat_id[i * 2, j * 2] = point
        except Exception as e:
            raise RuntimeError(f"Error reading DTM file: {e}")

    def sort_points_by_elevation(self):
        """
        Sorts the drainage points by elevation in ascending order.
        """
        self.dr_pt = list(tqdm(sorted(self.dr_pt, key=lambda x: x['Z']), desc="Sorting points by elevation", unit="point"))

    def process(self, header_file, dtm_file):
        """
        Executes the entire loading process: reads header, loads DTM, and sorts points.
        Parameters
        ----------
        header_file : str
            Path to the header file (e.g., 'header.dat').
        dtm_file : str
            Path to the DTM file (e.g., 'dtm.dat').
        """
        self.read_header(header_file)
        self.load_dtm(dtm_file)
        self.sort_points_by_elevation()

if __name__ == "__main__":
    
    header_path  = "../../out_scripts/cordevole_extrait_coord_header.dat"
    dat_dem_path = "../../out_scripts/cordevole_extrait_coord.dat"
    
    
    loader = LoadData()
    loader.process(header_path, dat_dem_path)

    print("Metadata:")
    print(f"delta_x: {loader.delta_x}, delta_y: {loader.delta_y}")
    print(f"Grid size: {loader.N}x{loader.M}")
    print(f"NoData value: {loader.nodata}")

    print("First 5 drainage points:")
    for point in loader.dr_pt[:5]:
        print(point)
