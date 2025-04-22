# -*- coding: utf-8 -*-
"""
Centralizes data structures and global parameters for hydrological analysis.
Originally translated and adapted from Fortran (pnt_net.f90).
"""


# --- Domain-specific structures ---

class DrainagePoint:
    """Represents a point in the drainage system.

    Attributes
    ----------
    i, j : int
        Row and column indices in the DEM.
    Z : float
        Elevation of the point.
    id_pnt : int
        Unique identifier of the drainage point.
    fldir : int
        Flow direction index (point id).
    fldir_ss : int
        Secondary flow direction index.
    A_in : int
        Inflow area in number of cells.
    upl : float
        Upland path length.
    dpl : float
        Downslope path length.
    sumdev : float
        Cumulative elevation deviation.
    id_endo : int
        Associated endorheic basin ID.
    ninf : int
        Number of inflow directions.
    inflow, Linflow, Sinflow : List
        IDs, lengths, and slopes of inflow segments.
    id_ch : int
        Associated channel ID.
    """
    def __init__(self, i, j, Z, id_pnt):
        self.i = i
        self.j = j
        self.Z = Z
        self.id_pnt = id_pnt
        self.fldir = None
        self.fldir_ss = None
        self.A_in = 0
        self.upl = 0.0
        self.dpl = 0.0
        self.sumdev = 0.0
        self.id_endo = 0
        self.ninf = 0
        self.inflow = []
        self.Linflow = []
        self.Sinflow = []
        self.id_ch = None
        
    def reset_flow_data(self):
        """Resets flow-related attributes (id_ch, inflow, Linflow) to initial state.
    
        This method clears any existing hydrological path information while preserving
        object references (pointers are not replaced, only their content is reset).
        """
        self.id_ch = None
        self.inflow.clear()
        self.Linflow.clear()

    def __repr__(self):
        attributes = {key: value
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"


class DrainagePointInflow:
    """Represents inflow data for a drainage point during network construction.

    Attributes
    ----------
    ninf : int
        Number of inflow directions.
    inflow : list of int
        IDs of the inflow points.
    """
    def __init__(self):
        self.ninf = 0
        self.inflow = []

    def __repr__(self):
        attributes = {key: value
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"


class DrainageNetwork:
    """Represents a drainage network or channel composed of drainage points.

    Attributes
    ----------
    id_ch : int
        Identifier of the channel.
    nel : int
        Number of points.
    id_pnts : List
        List of drainage point IDs.
    id_start_pt, id_end_pt : int
        IDs of the start and end points of the channel.
    length : float
        Total channel length.
    id_ch_out : int
        Downstream channel ID.
    n_jun : int
        Number of upstream tributaries.
    id_in : List
        List of IDs of upstream channels.
    n_path : int
        Number of segments in downstream path.
    id_path : list of int
        IDs of channels in the downstream path.
    id_endo : int
        Endorheic basin ID if any.
    sso, hso : int
        Stream segment and Horton orders.
    """
    def __init__(self, id_ch=None, nel=0):
        self.id_ch = id_ch
        self.nel = nel
        self.id_pnts = []
        self.id_start_pt = None
        self.id_end_pt = None
        self.length = 0.0
        self.id_ch_out = None
        self.n_jun = 0
        self.id_in = []
        self.n_path = 0
        self.id_path = []
        self.id_endo = 0
        self.sso = None
        self.hso = None

    def __repr__(self):
        attributes = {key: value
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"

        

class EndoPoint:
    """Represents an endorheic (closed-basin) point or an outflow point for endorheic basins.
    
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
      idms: list
          List of identifiers of saddle minima.
    """
    def __init__(self):
        self.id_eo = None 
        self.id_pnt = None 
        self.bas_type = None
        self.nsaddle = None
        self.beyo_sad = []
        self.idms = []
    
    def __repr__(self):
        attributes = {key: value
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"
    
    

class RidgePoint:
    """Represents a point on the ridgeline network.

    Attributes
    ----------
    i, j : int
        Indices in the working matrix (may differ from DEM grid).
    Z : float
        Elevation.
    id_pnt : int
        Unique ridge point ID.
    md : float
        Minimum mutual distance between associated drainage points.
    id_drpt1, id_drpt2 : int
        IDs of neighboring drainage points.
    A_in, A_in_min : float
        Max and min inflow areas from drainage neighbors.
    nen, n_jun : int
        Number of neighbors and junctions.
    id_neigh : list
        IDs of neighboring ridge points.
    id_sdl : int
        ID of associated saddle point.
    nrdl : int
        Number of ridge segments connected.
    id_rdl : list
        IDs of connected ridge segments.
    junc : int
        Junction flag (0: no junction, 1: junction).
    n_ptsa : float
        Number of points used for area estimation.
    """
    def __init__(self, i, j, Z, id_pnt):
        self.i = i
        self.j = j
        self.Z = Z
        self.id_pnt = id_pnt
        self.md = None
        self.id_drpt1 = None
        self.id_drpt2 = None
        self.A_in = None
        self.A_in_min = None
        self.nen = 0
        self.n_jun = 0
        self.id_neigh = []
        self.id_sdl = None
        self.nrdl = 0
        self.id_rdl = [None, None, None, None, None, None, None, None]
        self.junc = 0
        self.n_ptsa = 1

    def __repr__(self):
        attributes = {key: value
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"


class SaddlePoint:
    """Represents a saddle point connecting ridges and drainage basins.

    Attributes
    ----------
    id_pnt : int
        Saddle point ID.
    id_rdpt, id_rdpt2 : int
        IDs of associated ridge points.
    id_cis_endo : IDPointer
        Endorheic basin ID (cis-side).
    id_trans_out : IDPointer
        Downstream outflow ID (trans-side).
    A_endo : float
        Associated endorheic area.
    """
    def __init__(self, id_pnt):
        self.id_pnt = id_pnt
        self.id_rdpt = 0
        self.id_rdpt2 = 0
        self.id_cis_endo = None
        self.id_trans_out = None
        self.A_endo = None

    def __repr__(self):
        attributes = {key: value
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"



class OutNetwork:
    """Represents an outflow network composed of drainage point IDs.

    Attributes
    ----------
    nel : int
        Number of points.
    id_pnts : List
        List of drainage point IDs.
    """
    def __init__(self, nel):
        self.nel = nel
        self.id_pnts = []

    def __repr__(self):
        attributes = {key: value
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"



class RidgeNetwork:
    """Represents a segment of a ridgeline in the topographic network.

    Attributes
    ----------
    id_rdl : int
        Ridgeline ID.
    length : float
        Total length of the ridgeline.
    nel : int
        Number of points in the ridgeline.
    id_pnts : list of int
        List of ridge point IDs.
    Zmean : float
        Mean elevation.
    nrdpt_down : int
        Number of connected downstream ridge points.
    n_down : int
        Number of downstream ridge segments.
    Zmean_down : float
        Mean elevation of connected downstream segments.
    jun_el : int
        Index of the junction point within the segment.
    """
    def __init__(self, id_rdl):
        self.id_rdl = id_rdl
        self.length = 0.0
        self.nel = 0
        self.id_pnts = []
        self.Zmean = 0.0
        self.nrdpt_down = 0
        self.n_down = 0
        self.Zmean_down = 0.0
        self.jun_el = 0

    def __repr__(self):
        attributes = {key: value
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"
    
    