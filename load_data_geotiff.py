# -*- coding: utf-8 -*-
"""
Geotiff loader for Landmark algorithm
Include export methods
"""

#Extern import
import numpy as np
import rasterio
from tqdm import tqdm
import geopandas as gpd
from shapely.geometry import LineString, Point

#intern imports
from data_structures import DrainagePoint


class LoadData:
    """Class to load a DEM from a GeoTIFF file using rasterio.
    Stores the DEM data and its metadata for further processing.
    """
    def __init__(self):
        # self.dem = None  # DEM as a NumPy array
        self.transform = None  # Affine transform for georeferencing
        self.crs = None  # Coordinate reference system
        self.nodata = None  # NoData value
        self.delta_x = None # Grid spacing along the x-direction
        self.delta_y = None # Grid spacing along the y-direction
        self.N = None # Number of row, DEM size along the y-direction
        self.M = None # Number of column, DEM size along the x-direction
        self.xllcorner = None # X low left corner coordinate
        self.yllcorner = None # Y low left corner coordinate
        self.nodata = None # NoData value
        
        # Data structures for the DEM
        self.dr_pt = []  # List of drainage points (DrainagePoint objects)
        # Matrix for direct access to points (size: (N*2-1, M*2-1))
        # Points will be placed at indices (i*2, j*2) to allow for interpolation/neighborhood processing.
        self.mat_id = None



    def read_geotiff(self, geotiff_file: str):
        """Reads a DEM from a GeoTIFF file and stores its metadata.

        Parameters
        ----------
        geotiff_file : str
            Path to the GeoTIFF file.
        """
        try:
            with rasterio.open(geotiff_file) as src:
                dem = src.read(1)  # Read the first band
                self.transform = src.transform  # Store georeferencing transform
                self.crs = src.crs  # Store coordinate reference system
                self.nodata = -9999.00
                
                # Replace NoData values with NaN for easier processing
                dem = np.where(dem == self.nodata, np.nan, dem)
                
                self.delta_x = self.transform[0]
                self.delta_y = -self.transform[4]
                self.N, self.M = dem.shape
                self.xllcorner = self.transform[2]
                self.yllcorner = self.transform[5] - self.delta_y * self.N
                self.dr_pt = np.full((self.N * self.M), None, dtype =object)
                self.mat_id = np.full((self.N * 2 - 1, self.M * 2 - 1), None, dtype=object)

                
                for id_i in tqdm(range(self.N), desc="Loading DEM", unit="row"):
                    for id_j in range(self.M):
                        elevation = dem[id_i, id_j]
                        # Create a DrainagePoint (id_pnt starts at 1)
                        dp = DrainagePoint(
                            i=id_i,
                            j=id_j,
                            Z=float(elevation),
                            id_pnt=(id_i * self.M + id_j) + 1
                            )
                        self.dr_pt[id_i * self.M + id_j] = dp
                        # Place the point in the matrix.
                        self.mat_id[id_i * 2, id_j * 2] = dp


                # Build the list 'qoi' by sorting DrainagePoints in ascending Z
                # qoi will store the 'id_pnt' for each DrainagePoint
                # self.qoi = [dp.id_pnt.value for dp in sorted(self.dr_pt, key=lambda dp: dp.Z)]

                
        except Exception as e:
            raise RuntimeError(f"Error reading GeoTIFF file '{geotiff_file}': {e}")
            
 
    
    def export_slopelines_to_shapefile(self, output_shapefile: str) -> None:
        """Export the drainage network (slopelines) as a Shapefile.

        Only networks whose stream order (`hso`) is greater than or equal to a threshold
        (`self.hso_th`) are included. Each slopeline is converted to a LineString geometry,
        and a set of associated attributes is computed for each feature.

        Parameters
        ----------
        output_shapefile : str
            Path to the output Shapefile (.shp) file.

        Returns
        -------
        None

        Notes
        -----
        The output Shapefile contains one feature per complete slopeline (drainage path).
        The attribute table includes the following fields:

        +---------------+---------+--------------------------------------------------------------+
        | Field         | Type    | Description                                                  |
        +===============+=========+==============================================================+
        | id_ch         | int     | ID of the current slopeline (channel).                       |
        | nel           | int     | Number of elements (points) in the slopeline.                |
        | id_pnts       | str     | List of point IDs forming the slopeline (as string).         |
        | id_start_pt   | int     | ID of the first point in the slopeline.                      |
        | id_end_pt     | int     | ID of the last point in the slopeline.                       |
        | length        | float   | Total length of the slopeline (in spatial units).            |
        | id_ch_out     | int     | ID of the downstream slopeline.                              |
        | id_ch_main    | int     | ID of the main downstream slopeline in the hierarchy.        |
        | A_out         | int     | Contributing area at the second-to-last point (in map units²).|
        | n_jun         | int     | Number of upstream junctions flowing into this slopeline.    |
        | id_in         | list    | List of IDs of incoming slopelines (if any).                 |
        | n_path        | int     | Length of the downstream path in number of segments.         |
        | id_path       | str     | List of slopeline IDs along the downstream path.             |
        | id_endo       | int     | ID of the endorheic point (if applicable).                   |
        | sso           | int     | Stream Segment Order.                                        |
        | hso           | int     | Horton Stream Order.                                         |
        +---------------+---------+--------------------------------------------------------------+
        """
        lines = []
        attributes = []
        a_unit = self.delta_x * self.delta_y  # Area of one cell (assumed square)

        for net in tqdm(self.dr_net):
            if net.hso >= self.hso_th:
                coords = []
                for pnt_id in net.id_pnts.value:
                    dp = self.dr_pt[pnt_id - 1]
                    if dp:
                        # Transform row/column indices to spatial coordinates
                        x, y = self.transform * (dp.j, dp.i)
                        # Adjust to center of the cell
                        coords.append((x + self.delta_x / 2, y - self.delta_y / 2))

                if len(coords) > 1:
                    lines.append(LineString(coords))

                    # Prepare attributes dictionary for the current line
                    attributes.append({
                        "id_ch": int(net.id_ch.value) if net.id_ch.value is not None else 0,
                        "nel": int(net.nel),
                        "id_pnts": str([int(p) if p is not None else 0 for p in net.id_pnts.value])[:254],
                        "id_start_pt": int(net.id_start_pt.value) if net.id_start_pt.value is not None else 0,
                        "id_end_pt": int(net.id_end_pt.value) if net.id_end_pt.value is not None else 0,
                        "length": float(net.length),
                        "id_ch_out": int(net.id_ch_out.value) if net.id_ch_out.value is not None else 0,
                        "id_ch_main": int(net.id_path[-1]),
                        "A_out": int(self.dr_pt[net.id_pnts.value[-2] - 1].A_in * a_unit),
                        "n_jun": int(net.n_jun),
                        "id_in": [int(i.value) if i.value is not None else 0 for i in net.id_in.value]
                                if net.id_in is not None else [],
                        "n_path": int(net.n_path),
                        "id_path": str([int(i) if i is not None else 0 for i in net.id_path]
                                    if net.id_path is not None else [])[:254],
                        "id_endo": int(net.id_endo.value) if net.id_endo.value is not None else 0,
                        "sso": int(net.sso) if net.sso is not None else 0,
                        "hso": int(net.hso) if net.hso is not None else 0
                    })

        # Create GeoDataFrame and export to a Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=lines, crs=self.crs)
        gdf.to_file(output_shapefile)
        
    
    def export_slopelines_single_element_to_shapefile(self, output_shapefile: str) -> None:
        """Export the drainage network (slopelines) as individual segments into a Shapefile.

        Each segment between two consecutive drainage points is exported as a separate LineString,
        provided it passes the stream order (`hso`) and contributing area (`A_out`) thresholds.

        Parameters
        ----------
        output_shapefile : str
            Path to the output Shapefile (.shp) file.

        Returns
        -------
        None

        Notes
        -----
        The output Shapefile contains one feature per elementary segment of the slopeline network.
        The attribute table includes the following fields:

        +--------------+---------+--------------------------------------------------------------+
        | Field        | Type    | Description                                                  |
        +==============+=========+==============================================================+
        | id_ch        | int     | ID of the current slopeline segment.                         |
        | start_pnt    | int     | ID of the start point of the segment.                        |
        | end_pnt      | int     | ID of the end point of the segment.                          |
        | id_ch_out    | int     | ID of the downstream slopeline segment.                      |
        | id_ch_main   | int     | ID of the main downstream slopeline in the hierarchy.        |
        | A_out        | int     | Contributing area at the start point (in map units²).        |
        | length       | float   | Length of the segment (in spatial units, e.g. meters).       |
        | sso          | int     | Stream Segment Order.                                        |
        | hso          | int     | Horton Stream Order.                                         |
        +--------------+---------+--------------------------------------------------------------+
        """        
        lines = []
        attributes = []
        a_unit = self.delta_x * self.delta_y  # Area of one raster cell

        for net in tqdm(self.dr_net):
            if net.hso >= self.hso_th:
                points_ids = net.id_pnts.value

                for idx in range(len(points_ids) - 1):
                    dp1 = self.dr_pt[points_ids[idx] - 1]
                    dp2 = self.dr_pt[points_ids[idx + 1] - 1]

                    if dp1 and dp2:
                        a_out = dp1.A_in * a_unit
                        if a_out >= self.a_out_threshold:
                            # Compute coordinates from grid indices
                            x1, y1 = self.transform * (dp1.j, dp1.i)
                            x2, y2 = self.transform * (dp2.j, dp2.i)

                            # Adjust coordinates to cell centers
                            line = LineString([
                                (x1 + self.delta_x / 2, y1 - self.delta_y / 2),
                                (x2 + self.delta_x / 2, y2 - self.delta_y / 2)
                            ])

                            lines.append(line)

                            # Collect attributes for the current segment
                            attributes.append({
                                "id_ch": int(net.id_ch.value) if net.id_ch.value is not None else 0,
                                "start_pnt": int(dp1.id_pnt.value) if dp1.id_pnt.value is not None else 0,
                                "end_pnt": int(dp2.id_pnt.value) if dp2.id_pnt.value is not None else 0,
                                "id_ch_out": int(net.id_ch_out.value) if net.id_ch_out.value is not None else 0,
                                "id_ch_main": int(net.id_path[-1]),
                                "A_out": int(dp1.A_in * a_unit),
                                "length": float(line.length),
                                "sso": int(net.sso) if net.sso is not None else 0,
                                "hso": int(net.hso) if net.hso is not None else 0
                            })

        # Create GeoDataFrame and export to a Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=lines, crs=self.crs)
        gdf.to_file(output_shapefile)

        
        
    def export_drainage_point(self, output_shapefile: str) -> None:
        """Export all drainage points as a Shapefile.

        Each drainage point is represented as a Point geometry, with a comprehensive
        set of hydrological and topological attributes.

        Parameters
        ----------
        output_shapefile : str
            Path to the output Shapefile (.shp) file.

        Returns
        -------
        None

        Notes
        -----
        The output Shapefile contains one feature per drainage point. The attribute table
        includes the following fields:

        +------------+---------+--------------------------------------------------------------+
        | Field      | Type    | Description                                                  |
        +============+=========+==============================================================+
        | id_pnt     | int     | Unique identifier of the drainage point.                    |
        | i          | int     | Row index in the DEM grid.                                  |
        | j          | int     | Column index in the DEM grid.                               |
        | Z          | float   | Elevation of the point.                                     |
        | fldir      | int     | Flow direction (point id) from this point.                   |
        | fldir_ss   | int     | Secondary flow direction (for debug).                                  |
        | A_in       | int     | Contributing area (in number of cells).                     |
        | upl        | float   | Upland distance from the ridge to the point.                |
        | dpl        | float   | Downslope distance from the point to the outlet.            |
        | sumdev     | float   | Sum of elevation differences along the flow path.           |
        | id_endo    | int     | ID of the associated endorheic basin (if any).              |
        | ninf       | int     | Number of inflowing cells.                                  |
        | inflow     | str     | String-encoded list of inflowing cell IDs.                  |
        | Linflow    | str     | String-encoded list of inflowing line segment IDs.          |
        | Sinflow    | str     | String-encoded list of inflowing slope segment IDs.         |
        | id_ch      | int     | ID of the slopeline (channel) to which this point belongs.  |
        +------------+---------+--------------------------------------------------------------+
        """
        points = []
        attributes = []

        for dp in tqdm(self.dr_pt):
            # Convert grid indices to spatial coordinates
            x, y = self.transform * (dp.j, dp.i)
            # Adjust to center of the cell
            points.append(Point(x + self.delta_x / 2, y - self.delta_y / 2))

            # Prepare attribute dictionary for this drainage point
            attributes.append({
                "id_pnt": int(dp.id_pnt.value) if dp.id_pnt.value is not None else 0,
                "i": int(dp.i),
                "j": int(dp.j),
                "Z": float(dp.Z),
                "fldir": int(dp.fldir.value) if dp.fldir.value is not None else 0,
                "fldir_ss": int(dp.fldir_ss.value) if dp.fldir_ss.value is not None else 0,
                "A_in": int(dp.A_in),
                "upl": float(dp.upl),
                "dpl": float(dp.dpl),
                "sumdev": float(dp.sumdev),
                "id_endo": int(dp.id_endo.value) if dp.id_endo.value is not None else 0,
                "ninf": int(dp.ninf) if dp.ninf is not None else 0,
                "inflow": str(dp.inflow.value) if dp.inflow.value is not None else "NULL",
                "Linflow": str(dp.Linflow.value) if dp.Linflow.value is not None else "NULL",
                "Sinflow": str(dp.Sinflow.value) if dp.Sinflow.value is not None else "NULL",
                "id_ch": int(dp.id_ch.value) if dp.id_ch.value is not None else 0
            })

        # Create GeoDataFrame and export to a Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=points, crs=self.crs)
        gdf.to_file(output_shapefile)
            
            
    def export_ridge_point(self, output_shapefile: str) -> None:
        """Export ridge points as a Shapefile.

        Each ridge point corresponds to a cell located along the ridge network.
        These points are stored with their hydrological and topological properties.

        Parameters
        ----------
        output_shapefile : str
            Path to the output Shapefile (.shp) file.

        Returns
        -------
        None

        Notes
        -----
        The output Shapefile contains one feature per ridge point. The attribute table includes:

        +-------------+---------+------------------------------------------------------------+
        | Field       | Type    | Description                                                |
        +=============+=========+============================================================+
        | id_pnt      | int     | Ridge point ID.                                            |
        | i, j        | int     | Row and column indices in the grid (multiplied by 0.5).    |
        | Z           | float   | Elevation of the point.                                    |
        | md          | float   | Minimum distance to a drainage point (capped).             |
        | id_drpt1    | int     | ID of the first associated drainage point.                 |
        | id_drpt2    | int     | ID of the second associated drainage point.                |
        | A_in        | int     | Local contributing area.                                   |
        | A_in_min    | int     | Minimum contributing area between associated drains.       |
        | nen         | int     | Number of end points connected.                            |
        | n_jun       | int     | Number of ridge junctions.                                 |
        | id_neigh    | str     | Neighboring ridge point IDs.                               |
        | id_sdl      | int     | Associated saddle point ID.                                |
        | nrdl        | int     | Number of ridge lines.                                     |
        | id_rdl      | str     | List of ridge line IDs.                                    |
        | junc        | int     | Junction flag.                                             |
        | n_ptsa      | float   | Approximated area upstream of this point (in map units²). |
        +-------------+---------+------------------------------------------------------------+
        """        
        points = []
        attributes = []
        
        for rp in tqdm(self.rd_pt):
            x, y = self.transform * (rp.j / 2, rp.i / 2)
            points.append(Point(x + self.delta_x / 2, y - self.delta_y / 2))
    
            # Prepare attribute dictionary for this drainage point
            attributes.append({
                "id_pnt": int(rp.id_pnt) if rp.id_pnt is not None else 0,
                "i": int(rp.i),
                "j": int(rp.j),
                "Z": float(rp.Z),
                "md": min(float(rp.md), 99999999) if rp.md is not None else 0.0,  # Shapefile value limit
                "id_drpt1": int(rp.id_drpt1.value) if rp.id_drpt1 and rp.id_drpt1.value is not None else 0,
                "id_drpt2": int(rp.id_drpt2.value) if rp.id_drpt2 and rp.id_drpt2.value is not None else 0,
                "A_in": int(rp.A_in) if rp.A_in is not None else 0,
                "A_in_min": int(rp.A_in_min) if rp.A_in_min is not None else 0,
                "nen": int(rp.nen) if rp.nen is not None else 0,  
                "n_jun": int(rp.n_jun) if rp.n_jun is not None else 0,
                "id_neigh": str(rp.id_neigh) if rp.id_neigh is not None else 0,
                "id_sdl": int(rp.id_sdl) if rp.id_sdl is not None else 0,
                "nrdl": int(rp.nrdl) if rp.nrdl is not None else 0,
                "id_rdl": str(rp.id_rdl) if rp.id_rdl is not None else 0,
                "junc": int(rp.junc) if rp.junc is not None else 0,
                "n_ptsa": float(rp.n_ptsa*0.5*self.delta_x*self.delta_y ) if rp.n_ptsa is not None else 0
            })
    
        # Create GeoDataFrame and export to a Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=points, crs=self.crs)
        gdf.to_file(output_shapefile)
   
   
    
    def export_saddle_points(self, output_shapefile: str) -> None:
        """Export saddle points as a Shapefile.

        Saddle points represent topographic thresholds between basins.
        Only saddle points linked to an endorheic basin are exported.

        Parameters
        ----------
        output_shapefile : str
            Path to the output Shapefile (.shp) file.

        Returns
        -------
        None

        Notes
        -----
        The output Shapefile contains one feature per saddle point. The attribute table includes:

        +----------------+---------+-------------------------------------------------------------+
        | Field          | Type    | Description                                                 |
        +================+=========+=============================================================+
        | id_pnt         | int     | Saddle point ID.                                            |
        | id_rdpt        | int     | First associated ridge point ID.                            |
        | id_rdpt2       | int     | Second associated ridge point ID.                           |
        | Z              | float   | Elevation of the saddle point.                              |
        | id_cis_endo    | int     | ID of the cis-endorheic system.                             |
        | id_trans_out   | int     | ID of the downstream transfer ridge segment.                |
        | A_endo         | float   | Endorheic contributing area linked to the saddle point.      |
        +----------------+---------+-------------------------------------------------------------+
        """        
        points = []
        attributes = []
        
        for sp in tqdm(self.sdl_pt): 
            if sp.id_cis_endo.value != None :
                x, y = self.transform * (sp.j / 2, sp.i / 2)
                points.append(Point(x + self.delta_x / 2, y - self.delta_y / 2))
        
                attributes.append({
                    "id_pnt": int(sp.id_pnt) if sp.id_pnt is not None else 0,
                    "id_rdpt": int(sp.id_rdpt) if sp.id_rdpt is not None else 0,
                    "id_rdpt2": int(sp.id_rdpt2) if sp.id_rdpt2 is not None else 0,
                    "Z": float(sp.Z),
                    "id_cis_endo": int(sp.id_cis_endo.value) if sp.id_cis_endo and sp.id_cis_endo.value is not None else 0,
                    "id_trans_out": int(sp.id_trans_out.value) if sp.id_trans_out and sp.id_trans_out.value is not None else 0,
                    "A_endo": float(sp.A_endo) if sp.A_endo is not None else 0.0
                })
    
        gdf = gpd.GeoDataFrame(attributes, geometry=points, crs=self.crs)
        gdf.to_file(output_shapefile) 


    def export_endo_points(self, output_shapefile: str) -> None:
        """Export endorheic basin points as a Shapefile.

        Endorheic points are drainage endpoints where water does not reach the sea.
        The points are derived from the drainage points list.

        Parameters
        ----------
        output_shapefile : str
            Path to the output Shapefile (.shp) file.

        Returns
        -------
        None

        Notes
        -----
        The output Shapefile contains one feature per endorheic basin point. The attribute table includes:

        +-------------+---------+-------------------------------------------------------------+
        | Field       | Type    | Description                                                 |
        +=============+=========+=============================================================+
        | id_eo       | int     | Endorheic point ID.                                         |
        | id_pnt      | int     | ID of the associated drainage point.                        |
        | bas_type    | int     | Type of the basin (endorheic or exorheic).           |
        | nsaddle     | int     | Number of saddle points linked to this endorheic basin.     |
        | beyo_sad    | str     | List of saddle points beyond this basin.                    |
        | idms        | str     | ID of the main sink or depression system.                   |
        +-------------+---------+-------------------------------------------------------------+
        """        
        points = []
        attributes = []
        
        for ep in tqdm(self.l_endo_pt):
            dp = self.dr_pt[ep.id_pnt.value-1]
    
            x, y = self.transform * (dp.j, dp.i)
            points.append(Point(x + self.delta_x / 2, y - self.delta_y / 2))
            
    
            attributes.append({
                "id_eo": int(ep.id_eo.value) if ep.id_eo.value is not None else 0,
                "id_pnt": int(ep.id_pnt.value) if ep.id_pnt.value is not None else 0,
                "bas_type": int(ep.bas_type) if ep.bas_type is not None else 0,
                "nsaddle": int(ep.nsaddle) if ep.nsaddle is not None else 0,
                "beyo_sad": str(ep.beyo_sad.value) if ep.beyo_sad.value is not None else "NULL",
                "idms": str(ep.idms.value) if ep.idms.value is not None else "NULL",
            })
    
        gdf = gpd.GeoDataFrame(attributes, geometry=points, crs=self.crs)
        gdf.to_file(output_shapefile) 
        
        
    def export_ridgelines_to_shapefile(self, output_shapefile: str) -> None:
        """Export ridgeline networks as a Shapefile.

        Only ridgelines corresponding to watersheds that pass a stream order threshold (`jun_el > 0`)
        are included. Each ridgeline is exported as a LineString geometry with associated attributes.

        Parameters
        ----------
        output_shapefile : str
            Path to the output Shapefile (.shp) file.

        Returns
        -------
        None

        Notes
        -----
        The output Shapefile contains one feature per ridgeline. The attribute table includes:

        +-------------+---------+--------------------------------------------------------------+
        | Field       | Type    | Description                                                  |
        +=============+=========+==============================================================+
        | id_rdl      | int     | Ridgeline ID.                                                |
        | nel         | int     | Number of elements (points) in the ridgeline.                |
        | id_pnts     | str     | List of ridge point IDs (as a string).                       |
        | length      | float   | Total length of the ridgeline (in spatial units).            |
        | Zmean       | float   | Mean elevation along the ridgeline.                          |
        | nrdpt_down  | int     | Number of downstream ridge points.                           |
        | n_down      | int     | Number of ridgelines downstream.                             |
        | Zmean_down  | float   | Mean elevation of downstream ridgelines.                     |
        | Z_diff      | float   | Elevation difference (Zmean - Zmean_down).                   |
        | A_in        | float   | Contributing area at the upstream point (in map units²).     |
        | A_in_min    | float   | Minimum contributing area among connected upstream points.   |
        | jun_el      | int     | Index of the first junction point in the ridgeline.          |
        +-------------+---------+--------------------------------------------------------------+

        Remarks
        -------
        The full ridgeline is used (`range(0, net.nel)`) to avoid breaks at saddle points,
        which may be problematic when calculating curvature. Using `range(net.jun_el, net.nel)`
        would mimic the original Fortran behavior, but could result in visual discontinuities.
        """
        lines = []
        attributes = []

        for net in tqdm(self.rd_net):
            if net.jun_el > 0:  # Only include ridgelines of filtered watersheds
                coords = []

                # NOTE: Starting from 0 (not from jun_el) avoids visual breaks at saddle points.
                for cnt_pnt in range(0, net.nel):
                    pnt_id = net.id_pnts[cnt_pnt]
                    rp = self.rd_pt[pnt_id - 1]  # RidgePoint
                    if rp:
                        # Convert to spatial coordinates and center the point in the cell
                        x, y = self.transform * (rp.j / 2, rp.i / 2)
                        coords.append((x + self.delta_x / 2, y - self.delta_y / 2))

                if len(coords) > 1:
                    lines.append(LineString(coords))

                    attributes.append({
                        "id_rdl": int(net.id_rdl),
                        "nel": int(net.nel),
                        "id_pnts": str(net.id_pnts)[:254],
                        "length": float(net.length),
                        "Zmean": float(net.Zmean),
                        "nrdpt_down": int(net.nrdpt_down),
                        "n_down": int(net.n_down),
                        "Zmean_down": float(net.Zmean_down),
                        "Z_diff": float(net.Zmean - net.Zmean_down) if net.Zmean_down > 0 else 0,
                        "A_in": float((self.rd_pt[net.id_pnts[0] - 1].A_in + 1) * self.delta_x * self.delta_y),
                        "A_in_min": float((self.rd_pt[net.id_pnts[0] - 1].A_in_min + 1) * self.delta_x * self.delta_y),
                        "jun_el": int(net.jun_el)
                    })

        # Create GeoDataFrame and export to a Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=lines, crs=self.crs)
        gdf.to_file(output_shapefile)



    def export_ridgelines_single_element_to_shapefile(self, output_shapefile: str) -> None:
        """Export ridgeline networks as individual segments into a Shapefile.

        Each elementary segment between two consecutive ridge points is exported as a LineString,
        provided the spread area (`A_spread`) exceeds a given threshold. This allows filtering
        of minor ridges and focuses on morphologically significant features.

        Parameters
        ----------
        output_shapefile : str
            Path to the output Shapefile (.shp) file.

        Returns
        -------
        None

        Notes
        -----
        The output Shapefile contains one feature per ridgeline segment. The attribute table includes:

        +--------------+---------+--------------------------------------------------------------+
        | Field        | Type    | Description                                                  |
        +==============+=========+==============================================================+
        | id_rdl       | int     | Ridgeline ID to which the segment belongs.                   |
        | start_pnt    | int     | ID of the start ridge point of the segment.                  |
        | end_pnt      | int     | ID of the end ridge point of the segment.                    |
        | length       | float   | Segment length (in spatial units, e.g., meters).             |
        | jun_el       | int     | Index of the first junction in the full ridgeline.           |
        | A_spread     | float   | Approximated spread area covered by the segment (in m²).     |
        +--------------+---------+--------------------------------------------------------------+

        Remarks
        -------
        Starting the iteration at index 0 (rather than `jun_el`) avoids gaps at saddle points,
        which is useful for morphological analyses such as curvature computations.
        """
        lines = []
        attributes = []

        for net in tqdm(self.rd_net):
            if net.jun_el > 0:
                points_ids = net.id_pnts

                # NOTE: Full range used to preserve geometry continuity at saddle points
                for idx in range(0, net.nel - 1):
                    rp1 = self.rd_pt[points_ids[idx] - 1]
                    rp2 = self.rd_pt[points_ids[idx + 1] - 1]

                    if rp1 and rp2:
                        # Estimate spread area between the two ridge points
                        A_last = (1 / rp2.n_jun) * self.delta_x * self.delta_y
                        a_spread = (rp1.n_ptsa * 0.5 * self.delta_x * self.delta_y) + A_last

                        if a_spread >= self.a_spread_threshold:
                            # Compute coordinates of the two endpoints
                            x1, y1 = self.transform * (rp1.j / 2, rp1.i / 2)
                            x2, y2 = self.transform * (rp2.j / 2, rp2.i / 2)

                            # Build LineString segment centered in each raster cell
                            line = LineString([
                                (x1 + self.delta_x / 2, y1 - self.delta_y / 2),
                                (x2 + self.delta_x / 2, y2 - self.delta_y / 2)
                            ])

                            lines.append(line)

                            # Store segment attributes
                            attributes.append({
                                "id_rdl": int(net.id_rdl),
                                "start_pnt": int(rp1.id_pnt),
                                "end_pnt": int(rp2.id_pnt),
                                "length": float(line.length),
                                "jun_el": int(net.jun_el),
                                "A_spread": float(a_spread)
                            })

        # Create GeoDataFrame and export to a Shapefile
        gdf = gpd.GeoDataFrame(attributes, geometry=lines, crs=self.crs)
        gdf.to_file(output_shapefile)


if __name__ == "__main__":
    dtm_path = "./dtm.tif"
    
    
    loader = LoadData()
    loader.read_geotiff(dtm_path)
    metadata = loader.get_metadata()
    # print(metadata)
