# -*- coding: utf-8 -*-
"""
Centralizes data structures and global parameters for hydrological analysis.
From pnt_net.f90 Fortran code
"""

from dataclasses import dataclass, field


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
      fdir: int
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
    fdir: int = None
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
    i: int              # Row index in mat_id (not necessarily same as DEM cell indices)
    j: int              # Column index in mat_id
    Z: float            # Elevation (taken as the maximum of the two drainage points)
    id_pnt: int         # Unique identifier for this ridge point
    md: float = None    # Mutual distance (to be computed)
    id_sdl: int = None  # (For saddle point; here we just nullify it)
    id_drpt1: int = None  # id of one neighboring drainage point
    id_drpt2: int = None  # id of the other neighboring drainage point
    nen: int = 0        # Neighborhood number (set to 0 here)



class SaddlePoint:
    def __init__(self, id_pnt, id_rdpt):
        self.id_pnt = id_pnt        # Unique saddle point id
        self.id_rdpt = id_rdpt      # Associated RidgePoint id
        self.id_rdpt2 = 0           # Secondary ridge point id (default 0)
        self.id_cis_endo = None     # Endorheic basin id on the cis side
        self.id_trans_out = None    # Drainage point id on the trans (outflow) side

    
    ################ A compléter avec RidgelineNetwork #########
