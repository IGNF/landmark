# -*- coding: utf-8 -*-
"""
Transcript in python of the original load_data.f90

! The subroutine Load Data loads dtm data and header.
! It sorts also elevations in ascending order.
"""

#Extern imports
import numpy as np
import time
import os
from tqdm import tqdm


#intern imports
from data_structures import DrainagePoint


class LoadData:
    def __init__(self):
        # Header parameters
        self.delta_x = None
        self.delta_y = None
        self.N = None
        self.M = None
        self.xllcorner = None
        self.yllcorner = None
        self.nodata = None
        self.header_options = {}       # To store the 7 main header parameters
        self.landmarks_options = {}    # Options related to landmarks
        self.saving_options = {}       # Saving options

        # Data structures for the DEM
        self.dr_pt = []  # List of drainage points (DrainagePoint objects)
        # Matrix for direct access to points (size: (N*2-1, M*2-1))
        # Points will be placed at indices (i*2, j*2) to allow for interpolation/neighborhood processing.
        self.mat_id = None
        self.dem = None #A supprimer si inutile

        # Optionally: list of river sections, if present
        self.rs = []

    def read_header(self, header_file):
        """
        Reads the header file and extracts the parameters.
        Assumes the file is organized in three sections:
         - HEADER (7 lines)
         - LANDMARKS options (8 lines)
         - SAVING options (7 lines)
        Separator lines (dashes) and comment lines are ignored.
        """
        try:
            with open(header_file, 'r') as f:
                lines = f.readlines()
        except Exception as e:
            raise RuntimeError(f"Error opening header file '{header_file}': {e}")

        # Clean lines: remove unnecessary spaces and lines with only dashes.
        cleaned_lines = []
        for line in lines:
            l = line.strip()
            if l == "" or set(l) == {"-"}:
                continue
            cleaned_lines.append(l)

        # Use a state machine to distribute lines into sections
        current_section = None
        sections = {"HEADER": [], "LANDMARKS": [], "SAVING": []}
        for l in cleaned_lines:
            if l.startswith('!'):
                # Identify section based on the comment (case-insensitive)
                if "HEADER" in l.upper():
                    current_section = "HEADER"
                elif "LANDMARKS" in l.upper():
                    current_section = "LANDMARKS"
                elif "SAVING" in l.upper():
                    current_section = "SAVING"
                else:
                    current_section = None
                continue
            if current_section:
                sections[current_section].append(l)

        # Verify sections content
        if len(sections["HEADER"]) < 7:
            raise ValueError("Incomplete HEADER section in header file.")
        if len(sections["LANDMARKS"]) < 8:
            raise ValueError("Incomplete LANDMARKS options section in header file.")
        if len(sections["SAVING"]) < 7:
            raise ValueError("Incomplete SAVING options section in header file.")

        # Read HEADER section
        try:
            self.delta_x    = float(sections["HEADER"][0].split('=')[-1].strip())
            self.delta_y    = float(sections["HEADER"][1].split('=')[-1].strip())
            self.N          = int(sections["HEADER"][2].split('=')[-1].strip())
            self.M          = int(sections["HEADER"][3].split('=')[-1].strip())
            self.xllcorner  = float(sections["HEADER"][4].split('=')[-1].strip())
            self.yllcorner  = float(sections["HEADER"][5].split('=')[-1].strip())
            self.nodata     = float(sections["HEADER"][6].split('=')[-1].strip())
        except Exception as e:
            raise ValueError(f"Error reading HEADER parameters: {e}")

        self.header_options = {
            "delta_x": self.delta_x,
            "delta_y": self.delta_y,
            "N": self.N,
            "M": self.M,
            "xllcorner": self.xllcorner,
            "yllcorner": self.yllcorner,
            "nodata": self.nodata
        }

        # Read LANDMARKS options section
        try:
            self.landmarks_options = {
                "bas_flp": int(sections["LANDMARKS"][0].split('=')[-1].strip()),
                "length_th": float(sections["LANDMARKS"][1].split('=')[-1].strip()),
                "height_th": float(sections["LANDMARKS"][2].split('=')[-1].strip()),
                "slope_th": float(sections["LANDMARKS"][3].split('=')[-1].strip()),
                "A_th": float(sections["LANDMARKS"][4].split('=')[-1].strip()),
                "hso_th": int(sections["LANDMARKS"][5].split('=')[-1].strip()),
                "n_rs": int(sections["LANDMARKS"][6].split('=')[-1].strip()),
                "al": int(sections["LANDMARKS"][7].split('=')[-1].strip())
            }
            # Compute additional value if needed (ncell_th = A_th / (delta_x*delta_y))
            self.landmarks_options["ncell_th"] = self.landmarks_options["A_th"] / (self.delta_x * self.delta_y)
        except Exception as e:
            raise ValueError(f"Error reading LANDMARKS parameters: {e}")

        # Read SAVING options section
        try:
            self.saving_options = {
                "drpt_save": int(sections["SAVING"][0].split('=')[-1].strip()),
                "drnet_save": int(sections["SAVING"][1].split('=')[-1].strip()),
                "endo_save": int(sections["SAVING"][2].split('=')[-1].strip()),
                "saddle_save": int(sections["SAVING"][3].split('=')[-1].strip()),
                "spill_save": int(sections["SAVING"][4].split('=')[-1].strip()),
                "rdpt_save": int(sections["SAVING"][5].split('=')[-1].strip()),
                "rdnet_save": int(sections["SAVING"][6].split('=')[-1].strip())
            }
        except Exception as e:
            raise ValueError(f"Error reading SAVING options: {e}")

    def load_dtm(self, dtm_file):
        """
        Reads the DEM file (dtm_file).
        For each line, reads M floating-point values and creates a DrainagePoint object
        when the elevation is greater than nodata.
        Points are stored in self.dr_pt and placed in self.mat_id.
        """
        # Allocate the matrix for pointers
        self.mat_id = np.full((self.N * 2 - 1, self.M * 2 - 1), None, dtype=object)
        self.dem = np.full((self.N, self.M), None, dtype=float)
        start_time = time.process_time()

        try:
            with open(dtm_file, 'r') as f:
                # Use tqdm for the outer loop over rows
                for id_i in tqdm(range(self.N), desc="Loading DEM", unit="row"):
                    line = f.readline()
                    if not line:
                        raise ValueError(f"Insufficient data in DEM file at row {id_i+1}.")
                    try:
                        # Read a row containing M floating-point values
                        row_data = list(map(float, line.strip().split()))
                    except Exception as e:
                        raise ValueError(f"Error reading row {id_i+1}: {e}")
                    if len(row_data) < self.M:
                        raise ValueError(f"Row {id_i+1} has less than {self.M} columns.")

                    for id_j in range(self.M):
                        elevation = row_data[id_j]
                        if elevation > self.nodata:
                            # Create a DrainagePoint (id_pnt starts at 1)
                            dp = DrainagePoint(i=id_i, j=id_j, Z=elevation, id_pnt=len(self.dr_pt) + 1)
                            self.dr_pt.append(dp)
                            # Place the point in the matrix.
                            self.mat_id[id_i * 2, id_j * 2] = dp
                            self.dem[id_i, id_j] = dp.Z
        except Exception as e:
            raise RuntimeError(f"Error reading DEM file '{dtm_file}': {e}")

        finish_time = time.process_time()
        elapsed_time = finish_time - start_time
        print(f"DEM loading completed in {elapsed_time:.2f} seconds.")

    def load_rs(self, rs_file="rs.dat"):
        """
        Reads the river sections file (rs_file) if the number of sections (n_rs) > 0.
        The file is expected to have three comment lines followed by n_rs lines of data.
        Each data line should contain 2 values (rs_x and rs_y).
        """
        if self.landmarks_options.get("n_rs", 0) > 0:
            try:
                with open(rs_file, 'r') as f:
                    # Skip the first 3 lines (comments)
                    for _ in range(3):
                        f.readline()
                    # Use tqdm to track reading each river section
                    for cnt in tqdm(range(self.landmarks_options["n_rs"]), desc="Loading river sections", unit="section"):
                        line = f.readline()
                        if not line:
                            raise ValueError(f"Insufficient data in {rs_file} at section {cnt+1}.")
                        try:
                            rs_x, rs_y = map(float, line.strip().split())
                            self.rs.append((rs_x, rs_y))
                        except Exception as e:
                            raise ValueError(f"Error reading coordinates for section {cnt+1} in {rs_file}: {e}")
                print(f"River sections loaded, count: {len(self.rs)}.")
            except Exception as e:
                raise RuntimeError(f"Error reading river sections file '{rs_file}': {e}")
        else:
            print("No river sections to load (n_rs = 0).")

    def sort_points_by_elevation(self):
        """
        Sorts the list of drainage points by elevation (ascending order).
        This sort is essential for further processing (e.g., flow calculations).
        """
        self.dr_pt.sort(key=lambda dp: dp.Z)

    def process(self, header_file, dtm_file, rs_file=None):
        """
        Executes the full loading process:
          - Reads the header file
          - Loads the DEM file
          - Sorts the drainage points by elevation
          - Optionally, loads the river sections file
        """
        print("Reading header file...")
        self.read_header(header_file)
        print("Loading DEM file...")
        self.load_dtm(dtm_file)
        print("Sorting drainage points by elevation...")
        self.sort_points_by_elevation()
        if rs_file is not None and os.path.exists(rs_file):
            print("Reading river sections file...")
            self.load_rs(rs_file)
        else:
            print("No river sections file provided or file not found.")

    # def get_neighbors(self, dp):
    #     """
    #     Returns a list of neighboring points for a given drainage point dp in the grid.
    #     We consider the 8-connected neighborhood, taking into account that
    #     drainage points are stored at indices (i*2, j*2) in self.mat_id.
        
    #     Parameters
    #     ----------
    #     dp : DrainagePoint
    #         The drainage point for which to find neighbors.
        
    #     Returns
    #     -------
    #     List[DrainagePoint]
    #         List of found neighbors (can be empty).
    #     """
    #     # Convert the point's coordinates to matrix indices
    #     grid_i = dp.i * 2
    #     grid_j = dp.j * 2
    #     neighbors = []
    #     # Offsets for 8-connected neighborhood (steps of 2 units)
    #     offsets = [(-2, 0), (2, 0), (0, -2), (0, 2),
    #                (-2, -2), (-2, 2), (2, -2), (2, 2)]
    #     nrows, ncols = self.mat_id.shape
    #     for di, dj in offsets:
    #         ni, nj = grid_i + di, grid_j + dj
    #         if 0 <= ni < nrows and 0 <= nj < ncols:
    #             neighbor = self.mat_id[ni, nj]
    #             if neighbor is not None:
    #                 neighbors.append(neighbor)
    #     return neighbors
    
    
if __name__ == "__main__":
    # Paths to input files
    # header_path  = "../../cordevole/cordevole_extrait_coord_header.dat"
    # dtm_path = "../../cordevole/cordevole_extrait_coord.dat"
    
    header_path  = "../../out_scripts/cordevole_extrait_minimum_header.dat"
    dtm_path = "../../out_scripts/cordevole_extrait_minimum.dat"

    rs_path = "rs.dat"           # Optionnel : chemin vers rs.dat



    loader = LoadData()
    loader.process(header_path, dtm_path, rs_file=rs_path)

    print("\n--- Paramètres HEADER ---")
    for key, value in loader.header_options.items():
        print(f"{key} : {value}")

    print("\n--- Options LANDMARKS ---")
    for key, value in loader.landmarks_options.items():
        print(f"{key} : {value}")

    print("\n--- Options SAVING ---")
    for key, value in loader.saving_options.items():
        print(f"{key} : {value}")

    print(f"\nNombre total de points de drainage traités : {len(loader.dr_pt)}\n")

# Display a few points and their neighbors    
    nb_affichage = 5
    for dp in loader.dr_pt[:nb_affichage]:
        print(f"Point id {dp.id_pnt} à la position (i={dp.i}, j={dp.j}) avec élévation {dp.Z}")