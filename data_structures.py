# -*- coding: utf-8 -*-
"""
Centralizes data structures and global parameters for hydrological analysis.
From pnt_net.f90 Fortran code
"""


# Data Structures

class IDPointer:
    """Simulates a pointer in Python by encapsulating a mutable value."""
    def __init__(self, value=None):
        self.value = value


class ListPointer:
    """Simulates a pointer to a list in Python by encapsulating a mutable list."""
    def __init__(self, value=None):
        if value is None:
            value = []
        self.value = value
        
    def append(self, item):
        """Appends an item to the list."""
        self.value.append(item)



class DrainagePoint:
    """
    Class representing a drainage point, with attributes simulating pointers in Fortran.

    Attributes:
      i, j: int
          Row and column indices in the DEM.
      Z: float
          Elevation of the point.
      id_pnt: int
          Unique identifier of the drainage point.
      fldir: IDPointer
          Identifier of the downstream point.
      fldir_ss: IDPointer
          Identifier of the downstream point in endorheic basins.
      A_in: int
          Inflow area in number of cells.
      upl: float
          Upslope path length.
      dpl: float
          Downslope path length.
      sumdev: float
          Cumulative deviation (used in flow direction calculations).
      id_endo: int
          Identifier of the endorheic basin, if applicable.
      ninf: int
          Number of inflow points.
      inflow: list
          List of identifiers of inflow points.
      Linflow: list
          List of upstream flow path lengths for each inflow.
      Sinflow: list
          List of deviations (or slopes) associated with each inflow.
      id_ch: IDPointer
          Identifier of the channel (or network) to which this point belongs.
    """
    def __init__(self, i, j, Z, id_pnt):
        self.i = i
        self.j = j
        self.Z = Z
        self.id_pnt = IDPointer(id_pnt)
        self.fldir = IDPointer(None)      # Simulates a pointer
        self.fldir_ss = IDPointer(None)   # Simulates a pointer
        self.A_in = 0
        self.upl = 0.0
        self.dpl = 0.0
        self.sumdev = 0.0
        self.id_endo = IDPointer(None)
        self.ninf = 0
        self.inflow = ListPointer()  # Simulates a pointer to a list
        self.Linflow = ListPointer()
        self.Sinflow = ListPointer()
        self.id_ch = IDPointer(None)  # Simulates a pointer

    def __repr__(self):
        attributes = {key: (value.value if isinstance(value, (IDPointer, ListPointer)) else value)
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"
    
    
class DrainagePointInflow:
    """
    Class needed in drainage network saddle spill construction.
    Attributes:
      ninf: int
          Number of inflow points.
      inflow: list
          List of identifiers of inflow points.
    """
    
    def __init__(self):
        self.ninf = 0
        self.inflow = []

    def __repr__(self):
        attributes = {key: (value.value if isinstance(value, (IDPointer, ListPointer)) else value)
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"

        

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
    def __init__(self, id_ch, nel):
        self.id_ch = IDPointer(id_ch)  # Simulates a pointer
        self.nel = nel
        self.id_pnts = ListPointer()  # Simulates a pointer to a list of points
        self.id_start_pt = IDPointer(None)  # Simulates a pointer
        self.id_end_pt = IDPointer(None)  # Simulates a pointer
        self.length = 0.0
        self.id_ch_out = IDPointer(None)  # Simulates a pointer
        self.n_jun = 0
        self.id_in = ListPointer()  # Simulates a pointer to a list
        self.n_path = 0
        self.id_path = []
        self.id_endo = IDPointer(None)  # Simulates a pointer
        self.sso = None
        self.hso = None
        
    def __repr__(self):
        attributes = {key: (value.value if isinstance(value, (IDPointer, ListPointer)) else value)
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"
        

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
    def __init__(self):
        self.id_eo = IDPointer()  # Simulates a pointer
        self.id_pnt = IDPointer()  # Simulates a pointer
        self.bas_type = None
        self.nsaddle = None
        self.beyo_sad = ListPointer()  # Simulates a pointer to a list
        self.idms = ListPointer()  # Simulates a pointer to a list
    
    def __repr__(self):
        attributes = {key: (value.value if isinstance(value, (IDPointer, ListPointer)) else value)
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"
    
    

class RidgePoint:
    """
    Ridge point class, equivalent to `rd_pt_type` in Fortran.
    Represents a point on the ridgeline network with associated attributes.
    """
    def __init__(self, i, j, Z, id_pnt):
        self.i = i  # Row index in mat_id (not necessarily the same as DEM grid indices)
        self.j = j  # Column index in mat_id
        self.Z = Z  # Elevation of the ridge point
        self.id_pnt = id_pnt # Unique identifier of this ridge point
        
        self.md = None  # Mutual distance between two drainage points
        self.id_drpt1 = IDPointer(None) # ID of the first neighboring drainage point
        self.id_drpt2 = IDPointer(None)  # ID of the second neighboring drainage point
        self.A_in = None  # Maximum inflow area between the two drainage points
        self.A_in_min = None  # Minimum inflow area of the two drainage points
        
        self.nen = 0  # Number of neighboring points
        self.n_jun = 0  # Number of junctions (equals `nen` in Fortran)
        self.id_neigh = ListPointer()  # List of neighboring ridge points
        self.id_sdl = None  # ID of the corresponding saddle point
        
        self.nrdl = 0  # Number of ridgeline segments this point belongs to
        self.id_rdl = ListPointer()  # List of ridgeline IDs associated with this point
        self.junc = 0  # Junction flag (0 = no junction, 1 = junction)
        self.n_ptsa = 1  # Number of ridge points defining the spread area (related to DEM resolution)

    def __repr__(self):
        attributes = {key: (value.value if isinstance(value, (IDPointer, ListPointer)) else value)
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"


class SaddlePoint:
    """
    Saddle point class, equivalent to `sdl_pt_type` in Fortran.
    Represents a saddle point connecting different drainage areas.
    """
    def __init__(self, id_pnt):
        self.id_pnt = id_pnt  # Unique identifier of the saddle point
    
        self.id_rdpt = 0  # ID of the corresponding ridge point (linked to `rd_pt.id_pnt`)
        self.id_rdpt2 = 0 # ID of the secondary ridge point if the saddle is split
    
        self.id_cis_endo = IDPointer(None)  # ID of the endorheic basin on one side of the saddle
        self.id_trans_out = IDPointer(None)  # ID of the outflow point for the other side of the saddle
    
        self.A_endo = None  # Area of the endorheic basin affected by the saddle
    
    def __repr__(self):
        attributes = {key: (value.value if isinstance(value, (IDPointer, ListPointer)) else value)
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"



class OutNetwork:
    """
    Class representing an outflow network, with attributes simulating pointers in Fortran.
    """
    def __init__(self, nel):
        self.nel = nel # Number of points of of the channel
        self.id_pnts = ListPointer()  # Each points to the point in the dr_pt

    def __repr__(self):
        attributes = {key: (value.value if isinstance(value, (IDPointer, ListPointer)) else value)
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"



class RidgeNetwork:
    """
    Ridge network class, equivalent to `rd_net_type` in Fortran.
    Represents a ridgeline segment with associated attributes.
    """
    def __init__(self, id_rdl):
        self.id_rdl = IDPointer(id_rdl)  # Unique identifier of the ridgeline segment
        self.length = 0.0  # Length of the ridgeline segment
        self.nel = 0  # Number of points in the ridgeline
        self.id_pnts = ListPointer()  # List of point IDs in the ridgeline
        self.Zmean = 0.0  # Mean elevation of the ridgeline
        self.nrdpt_down = 0  # Number of ridge points connected, including the current ridge
        self.n_down = 0  # Number of ridges connected to this ridge segment
        self.Zmean_down = 0.0  # Mean elevation of the ridges connected to this one
        self.jun_el = 0  # Index of the ridgeline element that serves as a junction
        
    def __repr__(self):
        attributes = {key: (value.value if isinstance(value, (IDPointer, ListPointer)) else value)
                      for key, value in vars(self).items()}
        return f"{self.__class__.__name__}({attributes})"
    
    