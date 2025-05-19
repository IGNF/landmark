# -*- coding: utf-8 -*-
"""
Moretti, G., & Orlandini, S. (2023).
Thalweg and ridge network extraction from unaltered topographic data as a
basis for terrain partitioning. Journal of Geophysical Research: Earth Surface,
128, e2022JF006943.
https://doi.org/10.1029/2022JF006943
"""

import os
import argparse

from landmark.main_landmark_geotiff import landmark_processing

def main():
    parser = argparse.ArgumentParser(
        description="LANDMARK hydrological and geomorphological processing pipeline",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("mnt_path", type=str, help="Path to the input GeoTIFF DEM")
    parser.add_argument("--output_dir", type=str, default="../../outputs", help="Directory to save output files")
    parser.add_argument("--a_spread", type=float, default=1e5, help="Threshold for ridge dispersion area (in m^2)")
    parser.add_argument("--a_out", type=float, default=1e5, help="Threshold for drainage area (in m^2)")
    parser.add_argument("--hso_th", type=int, default=5, help="Horton stream order threshold")
    parser.add_argument("--curvature_slope", action="store_true", help="Enable slope/curvature computation")
    parser.add_argument("--n_pts_calc_slope", type=int, default=5, help="Number of points to compute slope")
    parser.add_argument("--no_data_values", nargs="*", type=float, default=[-9999, 0], help="List of NoData values")

    args = parser.parse_args()

    # Argument validation
    if not os.path.isfile(args.mnt_path):
        raise FileNotFoundError(f"Input file not found: {args.mnt_path}")
    if args.a_spread < 0:
        raise ValueError("a_spread must be non-negative.")
    if args.a_out < 0:
        raise ValueError("a_out must be non-negative.")
    if args.hso_th < 0:
        raise ValueError("hso_th must be non-negative.")
    if args.n_pts_calc_slope < 3:
        raise ValueError("n_pts_calc_slope must be at least 3.")
    if not isinstance(args.no_data_values, list) or not all(isinstance(v, (int, float)) for v in args.no_data_values):
        raise ValueError("no_data_values must be a list of numbers.")

    landmark_processing(
        mnt_path=args.mnt_path,
        output_dir=args.output_dir,
        a_spread=args.a_spread,
        a_out=args.a_out,
        hso_th=args.hso_th,
        curvature_slope=args.curvature_slope,
        n_pts_calc_slope=args.n_pts_calc_slope,
        no_data_values=args.no_data_values,
    )

if __name__ == "__main__":
    main()
