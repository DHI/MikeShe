########## MIKE SHE plugin for open bore holes ##########
# Subject:      Simulate open bore holes, calculate exchange flows between SZ layers and the bore hole and the resulting head in the bore hole
# Usage:        Attach to an MIKE SHE model as a plugin, adapt the few hard coded global variables (after imports)
# Dependencies: mikeio
# author:       uha@dhigroup.com
# date:         10/2022
#
#########################################################
#
# Setup:
# - gauging bore holes with small diameter
# - open at the top, no source/sink (pumping) applied
# - connected to sz layers (leakage coefficient specified)
# 
# Desired Results:
# - Exchange flow between each SZ layer and the bore hole
# - Resulting head in the bore hole
#
# Calculating flow: For simplicity we'll use
#   Q = (h_bh - h_i) * A_i * l_i
#     h_bh: head in the bore hole
#     h_i: head in layer i
#     A_i: bore hole exchange area in layer i
#     l_i: leakage in layer i
#   It is a simplification because this equation is for parallel flow lines (straight isobars), while around a well the flow paths are radial (circular isobars).
#   A common equation for wells was formulated by Dupuit-Thiem (actually 2 equations, 1 for free and one for confined aquifers, this is the one for free):
#     Q = Pi * k_f * (h_2^2 - h_1^2)/ln(r_2/r_1)
#       k_f: hydr. conductivity
#       r_1, r_2: distance locations 1 and 2 from bore hole
#       r_1, r_2: distance locations 1 and 2 from bore hole
#   In order to use this equation we'd need 2 points outside the well with known head and distance to the well. Potentially we could use neighbouring grid cells, or there may be other formulations allowing to use just a single head, or maybe it does not even make a huge difference - this could be investigated in the future.
#
# Reasoning: The flow in the bore hole is very fast compared to the SZ flow. If the pressure distribution in the SZ changes,
#            the head in the bore hole will almost instantaneously reach a new equilibrium so that the sum of all exchange
#            flows will be 0.
#
# The sum of all exchange flow bore hole-layer should be 0:
# 0 = SUM_0n[(h_bh - h_i) * A_i * l_i]
#   SUM_0n: Sum for all i in [0..n]
#
# Resolve to h_bh:
# 0 = SUM_0n[h_bh * A_i * l_i - h_i * A_i * l_i]
# 0 = SUM_0n[h_bh * A_i * l_i] - SUM_0n[h_i * A_i * l_i]
# 0 = h_bh * SUM_0n[A_i * l_i] - SUM_0n[h_i * A_i * l_i]
# SUM_0n[h_i * A_i * l_i] = h_bh * SUM_0n[A_i * l_i]
#________________________________________________________
#          SUM_0n[h_i * A_i * l_i]
#   h_bh = -----------------------
#            SUM_0n[A_i * l_i]
#________________________________________________________

#########################################################
#
# Leakage coefficient, conductivity and flow distance
# 
# The leakage coefficient can be interpreted as an abstraction of the hydr. conductivity k_f and the flow distance (let's call it df).
# leakge = k_f / df.
# It is not trivial to chose df.
# For simplicity let's assume this: The head in a cell is an average of the sub-cell distribution. If water is going into or out of the bore hole a cone will form around the cell center=bore hole location. The average cell head will represent the actual value everywhere on a circle around the center where the circle covers half the area of the cell. Such a circle has a radius of about 0.4 times the cell size, so let's use this for the flow distance (see variable flow_distance_cell_ratio).
#
#########################################################
#
# Current Limitations and Omissions
#
# UZ inclusion:
#   Due to a bug in MIKE SHE the plugin in the current form does not work with MIKE versions before Rel2023 if UZ is included and does have a shorter time step than SZ.
#
# Exchange area: 
#   Currently the exchange area is static, assuming that flow is through the entire depth of each layer. This is of course not always the case, especially for the top layer. The exchange area should be updated dynamically where the layer is not fully submerged, e. g. using the average head of bore hole and layer (limited to within the layer thickness).
#
# Driving head:
#   In the current formulation the driving head will be calculated using the bore hole head even if it is below the bottom of the layer. However in this situation instead the elevation of the bottom of the layer should be used.
#
# Head calculation in bore holes:
#   The advantage of the current determination of the head in the bore hole is that it is a closed equation. To solve the inaccuaracies described above (see "Exchange area" and "Driving head") an iterative approach would be needed. This has been implemented in an earlier version of this script and could be restored if needed.
#
# Timing:
#   Due to the explicit coupling of the bore hole to the SZ (the bore hole exchange and the SZ solver are executed one after the other) the exchange flows are based on the heads of the _previous_ SZ time step. The only way to change this would be to integrate the bore hole exchange directly into the SZ solver. Usually we don't even do this for components inside mshe. If a better temporal resolution is required, the time step length should be reduced.
#
# Number of wells:
#   Currently only one bore hole is supported with the location and diameter hard coded in the script. To support more bore holes we'd just have to use a list of bore holes and iterate over that list where required. It may make sense to read the required values from a csv file.
#
# Open length of bore hole:
#   Currently each bore hole will extend from the bottom to the top of the model. A limited depth and a potential coating in the top layer(s) would have to be implemented if needed.
#
# Custom result folder:
#   Not currently supported - not provided by MShePy, result folder needed for reading SZ layer bottoms.
#
# Memory:
#   Currently all data for the bore hole dfs0 file is collected over the entire simulation to be written at the very end. This may or may not become a problem with many bore holes, short SZ time steps and long simulations.
#

import MShePy as ms
import math
import mikeio
import os
import numpy as np
import textwrap as tw

# hard coded for bore hole 1:
bh1_name = "Bore Hole #1"
bh1_coords = [24.0, 26.0] # coordinates in model coordinate system
d_bh = 0.15  # borehole diameter (m)

flow_distance_cell_ratio = 0.4

class BoreHole:
  def __init__(self, name, x, y, d, zlyg, leak):
    self.i, self.j = ms.wm.gridCoordToCell(x, y)
    if not ms.wm.gridIsInternal(self.i, self.j):
      raise ValueError("Cannot create borehole {0} outside the model area!".format(name))
    self.name = name
    self.n_lay = len(zlyg) - 1
    self.leak = leak # bore hole to SZ exchange leakage coefficient
    self.zly = np.flip(zlyg[0:-1]) # bottom of layers - strip ground surface and change from bottom up to top down
    self.dz_ly = [zlyg[k] - zlyg[k - 1] for k in range(self.n_lay, 0 , -1)] # thickness of layers, from bottom up to top down
    self.d_bh = d # bore hole diameter
    self.ex_area = [dz * math.pi * self.d_bh for dz in self.dz_ly] # bore hole to SZ contact area
    self.heads = []
    self.times = []

  def __repr__(self):
    return __str__()
      
  def __str__(self):
    x, y = ms.wm.gridCellToCoord(self.i, self.j)
    return tw.dedent(f"""\
      Bore hole '{self.name}':
        (x,y):  ({x}, {y})
        (i,j):  ({self.i}, {self.j})
        zly:     [{', '.join(f'{x:9.4g}' for x in self.zly)}]
        dz_ly:   [{', '.join(f'{x:9.4g}' for x in self.dz_ly)}]
        leak:    [{', '.join(f'{x:9.4g}' for x in self.leak)}]
        ex_area: [{', '.join(f'{x:9.4g}' for x in self.ex_area)}]""")
    
def read_layer_bottoms():
  global setup_name
  # For now assume default result folder is being used. For extra safety we could check for a custom result folder in the she file, target "SYSTEM":
  #[SYSTEM]
  #   UseCustomResultFolder = true
  #   CustomResultFolder = |.\MyFolder|
  # For extra extra safety we'd have to check if the she file has been modified after preprocessing.
  fn = os.path.join(setup_dir, "{0}.she - Result Files/{0}_PreProcessed_3DSZ.dfs3".format(setup_name))
  item = "Lower level of computational layers in the saturated zone"
  zly = mikeio.read(fn, items=item)[0].values[0] # [0]: item, [0]: time step

  # other than mshe datasets zly is [z,y,x] and z is bottom up!
  return np.transpose(zly) # (z, y, x) to (x, y, z)

# Initializations
def postEnterSimulator():
  global bh_1
  global flow_distance_cell_ratio
  global bh1_name
  global bh1_coords
  global d_bh
  global setup_dir
  global setup_name
  # general values:
  dx = ms.wm.gridCellToCoord(1, 0)[0] - ms.wm.gridCellToCoord(0, 0)[0] # cell size
  _, _, z_ground = ms.wm.getValues(ms.paramTypes.DEM_Z)
  zg = np.array(z_ground[:])
  
  # retrieve other values for bore hole:
  i, j = ms.wm.gridCoordToCell(bh1_coords[0], bh1_coords[1])
  if not ms.wm.gridIsInternal(i, j):
    raise ValueError("Cannot create borehole at {0} which is outside the model area!".format(bh1_coords))
    
  setup_dir, setup_name = os.path.split(ms.wm.getSheFilePath())
  setup_name, _ = os.path.splitext(setup_name)

  # Currently bottom of sz layer elevation not available from MShePy - read thickness directly from PP'ed file:  
  zly = read_layer_bottoms()
  zgt = zg[..., None] # For ground elevations [x, y] add z dimension (with extent 1) for layer thickness calculation
  zlyg = np.concatenate([zly, zgt], 2) # Layer bottoms and ground elevation combined
  _, _, kh = ms.wm.getValues(ms.paramTypes.SZ_K_HOR) # Assuming k_f is static! Potentially it could be modified during run time via the api.
  leak = [khl / (flow_distance_cell_ratio * dx) for khl in kh[i, j]]
  bh_1 = BoreHole(bh1_name, bh1_coords[0], bh1_coords[1], d_bh, zlyg[i,j], leak)
  ms.wm.print(bh_1)

def preTimeStep():
  global bh_1
  _, sz_time, sz_heads = ms.wm.getValues(ms.paramTypes.SZ_HEAD) # sz_time is from previous time step, as are the heads!
  if(sz_time is None): # Not an SZ time step
    return

  ms.wm.print(f"\r\n\r\n########### {sz_time} ########################")
  szSource = ms.dataset(ms.paramTypes.SZ_SOURCE)
  dividend = 0
  divisor  = 0
  for i in range(bh_1.n_lay):
    if bh_1.leak[i] != 0:
      dividend += sz_heads[bh_1.i, bh_1.j, i] * bh_1.ex_area[i] * bh_1.leak[i]
      divisor  +=                           bh_1.ex_area[i] * bh_1.leak[i]

  if divisor > 0: # leakage to _any_ layer?
    bh_head = dividend / divisor
    for i in range(bh_1.n_lay):
      if bh_1.leak[i] != 0:
        # ex-flow for this layer/bh. Driving head _not_ limited to bottom of layer.
        qex = (bh_head - sz_heads[bh_1.i, bh_1.j, i]) * bh_1.ex_area[i] * bh_1.leak[i]
        szSource[bh_1.i, bh_1.j, i] = qex
        if bh_head < bh_1.zly[i]:
          msg =  "WARNING: Head in bore hole \"{:}\" has fallen to {:7.3f} m, below the bottom of layer no. {:} ({:7.3f} m) at {:}!".format(bh_1.name, bh_head, i, bh_1.zly[i], sz_time)
          msg += " The flow from this layer into the bore hole will be overestimated."
          # If this is an issue then an iterative approach for finding the solution is required, like in the 1st version of this plugin!
          ms.wm.log(msg)
  else:
    bh_head = float("NaN")
  
  bh_1.heads.append(bh_head)
  bh_1.times.append(sz_time)
  
  ms.wm.print(f"bh_head:   {bh_head:10.3f} m")
  l_per_m3 = 1000
  for head, flux in zip(sz_heads[bh_1.i, bh_1.j], szSource[bh_1.i, bh_1.j]):
    ms.wm.print(f"head: {head:10.3f} m, flux {flux * l_per_m3:10.3f} l/s")
  ms.wm.print("")
  ms.wm.setValues(szSource)

def preLeaveSimulator():
  global bh_1
  global setup_name
  preTimeStep() # capture end of last time step
  title = "Open bore hole heads"
  items = [mikeio.ItemInfo("Head {0}".format(bh_1.name), itemtype=mikeio.EUMType.Water_Level)]
  fname = os.path.join(setup_dir, "{0}.she - Result Files/{0}_BoreHoleHeads.dfs0".format(setup_name))

  dfs = mikeio.Dfs0()
  dfs.write(
    filename=fname,
    data=[bh_1.heads],
    datetimes=bh_1.times,
    items=items,
    title=title,
  )
