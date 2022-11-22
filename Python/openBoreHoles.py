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
#   In order to use this equation we'd need 2 points outside the well with known head and distance to the well. Potentially we could use neighbouring grid cells, 
#     or there may be other formulations allowing to use just a single head, or maybe it does not even make a huge difference - this could be investigated 
#     in the future.
#
# Reasoning: The flow in the bore hole is very fast compared to the SZ flow. If the pressure distribution in the SZ changes,
#            the head in the bore hole will almost instantaneously reach a new equilibrium so that the sum of all exchange
#            flows will be 0.
#
# Two different approaches are implemented to find the head:
#   1.: A closed equation with some limitations regarding accuracy, but quick to program and to execute
#   2.: An iterative solution - more correct, but less efficient
# Both methods can be used, change the source code below: Activate the desired method call in method update_flux_head. Note: The iterative 
#   method calls the closed solution for a starting point.
#
# The closed solution:
# The sum of all exchange flows bore hole-layer should be 0:
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
#
# The iterative solution:
# In each iteration the sum of all exchange flows bore hole-layer is calculated - the target for this is 0. The head for the next iteration will
# be calculated using the secant method. When abs value of the sum of flows is below the threshold (see constant EPS) the solution will be accepted.

#########################################################
#
# Leakage coefficient, conductivity and flow distance
# 
# The leakage coefficient can be interpreted as an abstraction of the hydr. conductivity k_f and the flow distance (let's call it df).
# leakge = k_f / df.
# It is not trivial to chose df.
# For simplicity let's assume this: The head in a cell is an average of the sub-cell distribution. If water is going into or out of the bore hole a cone will form around the cell center=bore hole location. The average cell head will represent the actual value everywhere on a circle around the center where the circle covers half the area of the cell. Such a circle has a radius of about 0.4 times the cell size, so let's use this for the flow distance (see variable FLOW_DISTANCE_CELL_RATIO).
#
#########################################################
#
# Current Limitations and Omissions
#
#### General:
# UZ inclusion:
#   Due to a bug in MIKE SHE the plugin in the current form does not work with MIKE versions before Rel2023 if UZ is included and does have a shorter time step than SZ.
#
# Timing:
#   Due to the explicit coupling of the bore hole to the SZ (the bore hole exchange and the SZ solver are executed one after the other) the exchange flows are 
#   based on the heads of the _previous_ SZ time step. The only way to change this would be to integrate the bore hole exchange directly into the SZ solver. 
#   Usually we don't even do this for components inside mshe. If a better temporal resolution is required, the time step length should be reduced.
#
# Custom result folder:
#   Not currently supported - not provided by MShePy, result folder needed for reading SZ layer bottoms.
#
# Memory:
#   Currently all data for the bore hole dfs0 file is collected over the entire simulation to be written at the very end. 
#   This may or may not become a problem with many bore holes, short SZ time steps and long simulations.
#
#### Limitations closed solution:
# Exchange area: 
#   The exchange area is static, assuming that flow is through the entire depth of each layer. 
#   This is of course not always the case, especially for the top layer. The exchange area should be updated dynamically 
#   where the layer is not fully submerged, e. g. using the average head of bore hole and layer (limited to within the layer thickness).
#
# Driving head:
#   The driving head will be calculated using the bore hole head even if it is below the bottom of the layer. 
#   However in this situation instead the elevation of the bottom of the layer should be used.
#  
# Open length of bore hole:
#   Each bore hole will extend from the bottom to the top of the model. A limited depth and a potential coating in the top layer(s) 
#   are not taken into account.
#

import MShePy as ms
import math
import os
import numpy as np
import textwrap as tw
import csv
import mikeio

COL_CNT_CFG = 6
FLOW_DISTANCE_CELL_RATIO = 0.4 # ()
MAX_ITER  = 10 # max number of iterations
EPS       = 1.0e-9 # (m^3/s) stop criterion, tolerance for remaining flow sum <> 0 of all layers<->bh
L_PER_M3  = 1000 # (l/m^3)
CFG_NAME = "openBoreHoles.txt"

class BoreHole:
  def __init__(self, name, x, y, d, bot, top, zlyg, leaks, sz_heads):
    self.i, self.j = ms.wm.gridCoordToCell(x, y)
    self.bot = bot
    self.top = top
    if not ms.wm.gridIsInternal(self.i, self.j):
      raise ValueError(f"Cannot create borehole \"{name}\" outside the model area!")
    self.name = name
    self.n_lay = len(zlyg) - 1
    self.leaks = leaks # bore hole to SZ exchange leakage coefficient
    self.zly = np.flip(zlyg[0:-1]) # bottom of layers - strip ground surface and change from bottom up to top down
    self.dz_ly = [zlyg[k] - zlyg[k - 1] for k in range(self.n_lay, 0 , -1)] # thickness of layers, from bottom up to top down
    self.d_bh = d # bore hole diameter
    self.ex_area = [dz * math.pi * self.d_bh for dz in self.dz_ly] # bore hole to SZ contact area (assuming bh from top to bottom of layer)
    self.heads = [] # result storing
    self.head = zlyg[-1] # initial head at ground surface for start of iterations only
    self.head_last = self.head
    nan = float("NaN")
    self.flo_sum = nan
    self.flo_sum_last = nan

  def __closed_solution(self, sz_heads, sz_time, silent = False):
    dividend = 0
    divisor  = 0
    for i in range(self.n_lay):
      if self.leaks[i] != 0:
        dividend += sz_heads[self.i, self.j, i] * self.ex_area[i] * self.leaks[i]
        divisor  +=                           self.ex_area[i] * self.leaks[i]

    if divisor > 0: # leakage to _any_ layer?
      self.head = dividend / divisor
      for i in range(self.n_lay):
        if self.leaks[i] != 0:
          # ex-flow for this layer/bh. Driving head _not_ limited to bottom of layer.
          qex = (self.head - sz_heads[self.i, self.j, i]) * self.ex_area[i] * self.leaks[i]
          sz_source[self.i, self.j, i] = qex
          if self.head < self.zly[i] and not silent:
            msg = f"WARNING: Head in bore hole \"{self.name}\" has fallen to {self.head:7.3f} m, below the bottom of layer no. {i} "\
                  f"({self.zly[i]:7.3f} m) at {sz_time}! The flow from this layer into the bore hole will be overestimated."
            ms.wm.log(msg)
    else:
      self.head = float("NaN")
    
  def __sum_flow(self, sz_heads):
    flo_sum = 0
    for ly in range(self.n_lay):
      if self.leaks[ly] == 0:
        continue

      ly_bot = self.zly[ly]
      if self.top < ly_bot: # If bore hole is all coated in this layer there cannot be any exchange
        continue

      ly_top = ly_bot + self.dz_ly[ly]
      if self.bot > ly_top: # If bore hole does not extend into this layer there cannot be any exchange
        continue

      sz_head = sz_heads[self.i, self.j, ly]

      # Find the top and bottom of contact between BH and SZ water.
      # Top of contact range cannot be above:
      ctc_top_bh = min(self.top, ly_top, max(self.bot, self.head)) # max(...): self.head should never be below bottom, but for extra safety...
      ctc_top_sz = min(self.top, ly_top, sz_head)

      ctc_bot = max(self.bot, ly_bot) # the same for SZ and BH!

      # Setting the contact range to half way between the WLs is really to account for head diff increasing from 0 to full head diff between WL of bh and WL of SZ
      contact_thick = (ctc_top_bh + ctc_top_sz) / 2 - ctc_bot
      contact_thick = max(contact_thick, 0)
      area_eff = self.ex_area[ly] * contact_thick / self.dz_ly[ly] # 0..1 fraction of full contact area
      bh_head_eff = max(self.head, ly_bot) # limit BH head to >= bottom of layer
      flo = (bh_head_eff - sz_head) * area_eff * self.leaks[ly]
      sz_source[self.i, self.j, ly] = flo
      flo_sum += flo
    return flo_sum
    
  def __iterative_solution(self, sz_heads, sz_time):
    ms.wm.print(" iteration              head         head_last          flow_sum      flo_sum_last")
    ms.wm.print("                         (m)               (m)             (l/s)             (l/s)")
    if math.isnan(self.flo_sum): # start of simulation - create starting point for iterations
      self.flo_sum_last = self.__sum_flow(sz_heads)
      self.__closed_solution(sz_heads, sz_time, silent=True)
    for its in range(MAX_ITER):
      self.flo_sum = self.__sum_flow(sz_heads)
      ms.wm.print("{:10g} {:17.7g} {:17.7g} {:17.4g} {:17.4g}".format(its, self.head, self.head_last, self.flo_sum * L_PER_M3, self.flo_sum_last * L_PER_M3))
      tmp = self.head
      if abs(self.flo_sum) < EPS: # solution good enough?
        its -= 1 # to save another EPS check for warning below
        self.flo_sum_last = self.flo_sum
        break
      self.head = self.head - (self.head - self.head_last) / (self.flo_sum - self.flo_sum_last) * (self.flo_sum) # secant method
      self.head = max(self.head, self.bot) # BH head below bottom end is invalid
      self.head_last = tmp
      self.flo_sum_last = self.flo_sum
    if its + 1 == MAX_ITER:
      msg = f"WARNING {sz_time}: Max. number of {MAX_ITER} iterations in bore hole \"{self.name}\", "\
            f"remaining net exchange flow of {self.flo_sum * L_PER_M3:7.3g} l/s "\
            f"is above threshold of {EPS * L_PER_M3:7.3g} l/s and will cause a water balance error! New head is {self.head:7.3f} m."
      ms.wm.log(msg)

  def update_flux_head(self, sz_heads, sz_time):
    # Activate either one of these solution methods:
    # self.__closed_solution(sz_heads, sz_time)
    self.__iterative_solution(sz_heads, sz_time)
    self.heads.append(self.head)

    ms.wm.print(f"bh_head \"{self.name}\":   {self.head:10.3f} m")
    l = 1
    for head, flux in zip(sz_heads[self.i, self.j], sz_source[self.i, self.j]):
      ms.wm.print(f"    GW layer {l:3} - head: {head:10.3f} m, flux {flux * L_PER_M3:10.3g} l/s")
      l += 1

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
        leak:    [{', '.join(f'{x:9.4g}' for x in self.leaks)}]
        ex_area: [{', '.join(f'{x:9.4g}' for x in self.ex_area)}]""")
    
def read_layer_bottoms():
  # For now assume default result folder is being used. For extra safety we could check for a custom result folder in the she file, target "SYSTEM":
  #[SYSTEM]
  #   UseCustomResultFolder = true
  #   CustomResultFolder = |.\MyFolder|
  # For extra extra safety we'd have to check if the she file has been modified after preprocessing.
  fn = os.path.join(setup_dir, f"{setup_name}.she - Result Files/{setup_name}_PreProcessed_3DSZ.dfs3")
  item = "Lower level of computational layers in the saturated zone"
  zly = mikeio.read(fn, items=item)[0].values[0] # [0]: item, [0]: time step

  # other than mshe datasets zly is [z,y,x] and z is bottom up!
  return np.transpose(zly) # (z, y, x) to (x, y, z)

# Initializations
def postEnterSimulator():
  global bhs
  global times
  global setup_dir
  global setup_name
  global sz_source
  # general values:
  bhs = []
  times = []
  dx = ms.wm.gridCellToCoord(1, 0)[0] - ms.wm.gridCellToCoord(0, 0)[0] # cell size
  _, _, z_ground = ms.wm.getValues(ms.paramTypes.DEM_Z)
  zg = np.array(z_ground[:])

  setup_dir, setup_name = os.path.split(ms.wm.getSheFilePath())
  setup_name, _ = os.path.splitext(setup_name)

  # Currently bottom of sz layer elevation not available from MShePy - read thickness directly from PP'ed file:
  zly = read_layer_bottoms()
  zgt = zg[..., None] # For ground elevations [x, y] add z dimension (with extent 1) for layer thickness calculation
  zlyg = np.concatenate([zly, zgt], 2) # Layer bottoms and ground elevation combined
  _, _, kh = ms.wm.getValues(ms.paramTypes.SZ_K_HOR) # Assuming k_f is static! Potentially it could be modified during run time via the api.

  cfg_path = os.path.join(setup_dir, CFG_NAME)
  if not os.path.exists(cfg_path):
    raise ValueError(f"Looking for bore hole cfg file \"{cfg_path}\", but it does not exist!")

  sz_source = ms.dataset(ms.paramTypes.SZ_SOURCE)
  _, sz_time, sz_heads = ms.wm.getValues(ms.paramTypes.SZ_HEAD)
  with open(cfg_path) as csvFile:
    cfg = csv.reader(csvFile, delimiter=';')
    i = 0
    for tokens in cfg:
      if len(tokens) == 0 or tokens[0].strip()[0] == '#':
        continue
      if not len(tokens) == COL_CNT_CFG:
        ms.wm.log("Expected bore hole configuration file format:")
        ms.wm.log("  name; x; y; d; bot; top")
        raise ValueError(f"Invalid bore hole configuration file {cfg_path}: expected {COL_CNT_CFG} columns but found {len(tokens)}.")
      i += 1
      if i > 1:
        # retrieve other values for bore hole:
        bh_name = tokens[0].strip()
        bh_x = float(tokens[1])
        bh_y = float(tokens[2])
        bh_d = float(tokens[3])
        bh_bot = float(tokens[4])
        bh_top = float(tokens[5])
        if bh_bot > bh_top:
          msg = f"Cannot create borehole \"{bh_name}\" with the top end below the bottom end!\n"
          raise ValueError(msg)
        i, j = ms.wm.gridCoordToCell(bh_x, bh_y)
        if not ms.wm.gridIsInternal(i, j):
          msg = f"Cannot create borehole \"{bh_name}\" outside the model area!\n"\
                f"  x={bh_x}, y={bh_y};  i={i}, j={j}"
          raise ValueError(msg)
        ms.wm.log(f"i={i}, j={j}")
        leaks = [khl / (FLOW_DISTANCE_CELL_RATIO * dx) for khl in kh[i, j]] # leakages in this SZ column
        bh = BoreHole(bh_name, bh_x, bh_y, bh_d, bh_bot, bh_top, zlyg[i,j], leaks, sz_heads)
        ms.wm.print(bh)
        bhs.append(bh)

def preTimeStep():
  _, sz_time, sz_heads = ms.wm.getValues(ms.paramTypes.SZ_HEAD) # sz_time and heads are from the end of the previous time step!
  if(sz_time is None): # Not an SZ time step
    return

  ms.wm.print(f"\r\n\r\n########### {sz_time} ########################")
  sz_source.value(0.0)
  for bh in bhs:
    bh.update_flux_head(sz_heads, sz_time)
    ms.wm.print("")
  times.append(sz_time)
  ms.wm.setValues(sz_source)

def preLeaveSimulator():
  preTimeStep() # capture end of last time step
  title = "Open bore hole heads"
  items = [mikeio.ItemInfo(f"Head \"{bh.name}\"", itemtype=mikeio.EUMType.Water_Level) for bh in bhs]
  fname = os.path.join(setup_dir, f"{setup_name}.she - Result Files/{setup_name}_BoreHoleHeads.dfs0")

  dfs = mikeio.Dfs0()
  dfs.write(
    filename = fname,
    data = [bh.heads for bh in bhs],
    datetimes = times,
    items = items,
    title = title,
  )
