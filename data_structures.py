# -*- coding: utf-8 -*-
"""
Centralizes data structures and global parameters for hydrological analysis.
From pnt_net.f90 Fortran code
"""

from dataclasses import dataclass, field
from typing import List, Optional

# Data Structures

@dataclass
class DrainagePoint:
    """
    Represents a drainage point in the DEM.
    
    Attributes:
      i, j: int
          Row and column indices in the DEM (0-indexed).
      Z: float
          Elevation at the point.
      id_pnt: int
          Unique identifier of the drainage point.
      fldir: int
          Identifier of the downstream (outflow) point.
      fldir_ss: int
          flow direction endorheic saddle spill, index of outflow point.
      A_in: int
          Inflow area in number of cells.
      upl: float
          Maximum upstream length (used in flow accumulation).
      dpl: float
          Downslope path length.
      sumdev: float
          Cumulative deviation (used in flow direction calculations).
      id_endo: int
          Identifier of the endorheic basin, if applicable.
      ninf: int
          Number of inflows (number of upstream points).
      inflow: list
          List of identifiers of inflow points.
      Linflow: list
          List of upstream flow path lengths associated with each inflow.
      Sinflow: list
          List of deviations (or slopes) associated with each inflow.
      id_ch: int
          Identifier of the channel (or network) to which this point belongs.
    """
    i: int
    j: int
    Z: float
    id_pnt: int
    fldir: int = None
    fldir_ss: int = None
    A_in: int = 0
    upl: float = 0.0
    dpl: float = 0.0
    sumdev: float = 0.0
    id_endo: int = None
    ninf: int = 0
    inflow: list = field(default_factory=list)
    Linflow: list = field(default_factory=list)
    Sinflow: list = field(default_factory=list)
    id_ch: int = None

        
@dataclass
class DrainageNetwork:
    """
    Represents a drainage network (or channel) composed of several drainage points.
    
    Attributes:
      id_ch: int
          Identifier of the channel.
      nel: int
          Number of points of the channel
      id_pnts: list
          List of identifiers of drainage points belonging to the channel.
      id_start_pt: int
          Identifier of the head point of the channel.
      id_end_pt: int
          Identifier of the terminal point of the channel.
      length: float
          Total length of the channel.
      id_ch_out: int
          Identifier of the main channel (for outflow direction).
      n_jun: int
          Number of tributaries (number of upstream connected channels).
      id_in: list
          List of identifiers of tributary channels.
      n_path: int
          Number of flow paths (optional, for additional analyses).
      id_path: int
          id of the downslope channel path.
      id_endo: int
          id of endorheic basin that the current channel belongs.
      sso: int
          Strhaler Stream Order for hso calculation, at the end sso=hso.
      hso: int
          Horton Stream Order
    """
    id_ch: int
    nel: int 
    id_pnts: list = field(default_factory=list)
    id_start_pt: int = None
    id_end_pt: int = None
    length: float = 0.0
    id_ch_out: int = None
    n_jun: int = 0
    id_in: list = field(default_factory=list)
    n_path: int = 0
    id_path: int = None
    id_endo: int = None
    sso: int = None
    hso: int = None

@dataclass
class EndoPoint:
    """
    Represents an endorheic (closed-basin) point or an outflow point for endorheic basins.
    
    Attributes:
      id_eo: int
          Unique identifier of the endorheic point.
      id_pnt: int
          Identifier of the associated drainage point.
      bas_type: int
          Basin type (e.g., 0 for outflow, 1 for endorheic).
      nsaddle: int
          Number of saddles (if applicable).
      beyo_sad: list (optional)
          List of identifiers of points located on the other side of the saddle (if needed).
      idms: list (optional)
          List of identifiers of saddle minima.
    """
    id_eo: int
    id_pnt: int
    bas_type: int
    nsaddle: int = 0
    beyo_sad: list = field(default_factory=list)
    idms: list = field(default_factory=list)
    
    
@dataclass
class RidgePoint:
    """
    Ridge point class, equivalent to `rd_pt_type` in Fortran.
    Represents a point on the ridgeline network with associated attributes.
    """
    i: int  # Row index in mat_id (not necessarily the same as DEM grid indices)
    j: int  # Column index in mat_id
    Z: float  # Elevation of the ridge point
    id_pnt: int  # Unique identifier of this ridge point
    
    md: Optional[float] = None  # Mutual distance between two drainage points
    id_drpt1: Optional[int] = None  # ID of the first neighboring drainage point
    id_drpt2: Optional[int] = None  # ID of the second neighboring drainage point
    A_in: Optional[int] = None  # Maximum inflow area between the two drainage points
    A_in_min: Optional[int] = None  # Minimum inflow area of the two drainage points
    
    nen: int = 0  # Number of neighboring points
    n_jun: int = 0  # Number of junctions (equals `nen` in Fortran)
    id_neigh: List[int] = field(default_factory=list)  # List of neighboring ridge points
    id_sdl: Optional[int] = None  # ID of the corresponding saddle point
    
    nrdl: int = 0  # Number of ridgeline segments this point belongs to
    id_rdl: List[int] = field(default_factory=list)  # List of ridgeline IDs associated with this point
    junc: int = 0  # Junction flag (0 = no junction, 1 = junction)
    n_ptsa: int = 1  # Number of ridge points defining the spread area (related to DEM resolution)


@dataclass
class SaddlePoint:
    """
    Saddle point class, equivalent to `sdl_pt_type` in Fortran.
    Represents a saddle point connecting different drainage areas.
    """
    id_pnt: int  # Unique identifier of the saddle point

    id_rdpt: Optional[int] = None  # ID of the corresponding ridge point (linked to `rd_pt.id_pnt`)
    id_rdpt2: int = 0  # ID of the secondary ridge point if the saddle is split

    id_cis_endo: Optional[int] = None  # ID of the endorheic basin on one side of the saddle
    id_trans_out: Optional[int] = None  # ID of the outflow point for the other side of the saddle

    A_endo: Optional[int] = None  # Area of the endorheic basin affected by the saddle
    


@dataclass
class RidgeNetwork:
    """
    Ridge network class, equivalent to `rd_net_type` in Fortran.
    Represents a ridgeline segment with associated attributes.
    """
    id_rdl: int  # Unique identifier of the ridgeline segment
    
    length: float = 0.0  # Length of the ridgeline segment
    nel: int = 0  # Number of points in the ridgeline

    id_pnts: List[int] = field(default_factory=list)  # List of point IDs in the ridgeline

    Zmean: float = 0.0  # Mean elevation of the ridgeline
    nrdpt_down: int = 0  # Number of ridge points connected, including the current ridge
    n_down: int = 0  # Number of ridges connected to this ridge segment
    Zmean_down: float = 0.0  # Mean elevation of the ridges connected to this one

    jun_el: int = 0  # Index of the ridgeline element that serves as a junction
    
    