########## MIKE SHE plugin for open bore holes ##########
# Subject:      Simulate open bore holes, calculate exchange flows between SZ layers and the bore hole and the resulting head in the bore hole
# Usage:        Create a config file, attach to a MIKE SHE model as a plugin, adapt the few hard coded global variables (after imports)
# Dependencies: mikeio (which requires: numpy, scipy - these are also used here)
# author:       uha@dhigroup.com
# date:         10/2022
#
#########################################################
#
# Configuration file
# The configuration file is a csv file using ';' as a separator. The first line contains
# headers, subsequent lines either data or comments denoted by '#'. The columns are:
#   name: An arbitrary name for the bore hole (must not contain ';')
#   x:    x-coordinate of bore hole (model coordinate system)
#   y:    y-coordinate of bore hole (model coordinate system)
#   d:    bore hole diameter (m)
#   bot:  Elevation of bore hole bottom (m)
#   top:  Elevation of bore hole top (m) (really the top of the filter for bore holes that have a coated top end)
# Decimal separator is '.'. Extra white space may be present anywhere and will be ignored.
# 
# Example:
# +----------------------------------------------------+
# |               name;    x;    y;    d;   bot;   top |
# | # Bore hole right in the center:                   |
# |   Center bore hole; 24.2; 26.1; 0.15; -80.3;  10.7 |
# +----------------------------------------------------+
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
# Both methods can be used, change the source code below: Activate the desired method call in method update_flux_head.
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
# Uses Brent's method from scipy. This is combining the bisection method, the secant method and inverse quadratic interpolation.
# In each iteration the sum of all exchange flows bore hole-layer is calculated - the target for this is 0.

#########################################################
#
# Leakage coefficient, conductivity and flow distance
#
# The leakage coefficient can be interpreted as an abstraction of the hydr. conductivity k_f and the flow distance (let's call it df).
# leakge = k_f / df.
# It is not trivial to chose df.
#
#   For simplicity let's assume this: The head in a cell is an average of the sub-cell distribution. If water is going into or out of the bore hole
# a cone will form around the cell center=bore hole location. The average cell head will represent the actual value everywhere on a circle around
# the center where the circle covers half the area of the cell. Such a circle has a radius of about 0.4 times the cell size, so let's use this for
# the flow distance (see variable FLOW_DISTANCE_CELL_RATIO).
#
#########################################################
#
# Current Limitations and Omissions
#
#### General:
# Timing:
#   Due to the explicit coupling of the bore hole to the SZ (the bore hole exchange and the SZ solver are executed one after the other) the exchange flows are
#   based on the heads of the _previous_ SZ time step. If a better temporal resolution is required, the time step length should be reduced.
#
# Custom result folder:
#   Not currently supported - not provided by MShePy, result folder needed for reading SZ layer bottoms.
#
# Memory:
#   Currently all data for the bore hole dfs0 file is collected over the entire simulation to be written at the very end.
#   This may or may not become a problem with many bore holes, short SZ time steps and long simulations.
# 
# Time varying SZ conductivities:
#   Currently the assumption is that SZ conductivities are constant throughout the simulation. If the conductivities 
#   are time varying instead the bore hole leakages will not be updated accordingly after initialization.
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
import csv
import math
import mikeio
import numpy as np
import os
from scipy import optimize
import textwrap as tw

COL_CNT_CFG = 6
FLOW_DISTANCE_CELL_RATIO = 0.4 # ()

# Max number of iterations in Brent's method. 
#   NOTE 1: This does NOT include bracket finding which takes 2 or more iterations!
#   NOTE 2: The solution may still be below the threshold, but selection of the best solution may be skipped!
MAX_ITER  = 10
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
    self.flow_cache = {} # For Brent's method we need to pre-calculate 2 starting points, cache these to avoid double-calculating
    self.iterations = [] # list of iterations per time step
    self.no_converge = 0 # count of total non-convergences
    self.its = 0         # counter for iterations within one time step

  def __closed_solution(self, sz_heads, sz_time, silent = False):
    dividend = 0
    divisor  = 0
    for i in range(self.n_lay):
      if self.leaks[i] != 0:
        dividend += sz_heads[self.i, self.j, i] * self.ex_area[i] * self.leaks[i]
        divisor  +=                               self.ex_area[i] * self.leaks[i]

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

  def __calc_flow(self, head, sz_heads):
    if head in self.flow_cache:
      return self.flow_cache[head]
    self.its += 1

    flow_sum = 0
    for ly in range(self.n_lay):
      if self.leaks[ly] == 0:
        continue

      ly_bot = self.zly[ly]
      if self.top < ly_bot: # If bore hole is all coated in this layer there cannot be any exchange
        continue

      ly_top = ly_bot + self.dz_ly[ly]
      if self.bot > ly_top: # If bore hole does not extend into this layer there cannot be any exchange
        continue

      sz_head_eff = max(sz_heads[self.i, self.j, ly], self.bot) # Driving head must not consider SZ head below bottom of BH

      # Find the top and bottom of contact between BH and SZ.
      # - ctc_bh is the length whith water in the bore hole in this layer
      # - ctc_sz is the depth whith water in this sz layer (in contact with bh)
      # Top of contact range cannot be above:
      ctc_top_bh = min(self.top, ly_top, max(self.bot, head)) # max(...): head should never be below bottom, but for extra safety...
      ctc_top_sz = min(self.top, ly_top, sz_head_eff)

      ctc_bot = max(self.bot, ly_bot) # the same for SZ and BH!

      ctc_bh = max(ctc_top_bh - ctc_bot, 0)
      ctc_sz = max(ctc_top_sz - ctc_bot, 0)

      # Setting the contact range to half way between the WLs is really to account for head diff increasing from 0 to full head diff between WL of bh and WL of SZ
      contact_thick = (ctc_bh + ctc_sz) / 2

      area_eff = self.ex_area[ly] * contact_thick / self.dz_ly[ly] # 0..1 fraction of full contact area
      bh_head_eff = max(head, ly_bot)  # Driving head must not consider BH head below bottom of layer
      flow = (bh_head_eff - sz_head_eff) * area_eff * self.leaks[ly]
      sz_source[self.i, self.j, ly] = flow
      flow_sum += flow
    ms.wm.print(f"{self.its:10g} {head:.9f} {flow_sum * L_PER_M3:.4g}")
    self.flow_cache[head] = flow_sum
    return flow_sum

  def __iterative_solution(self, sz_heads, sz_time):
    ms.wm.print(f"Bore hole '{self.name}'")
    ms.wm.print(f"======================================")
    ms.wm.print(" iteration   head         flow_sum")
    ms.wm.print("              (m)            (l/s)")
    self.its = 0

    flow_sum_0 = self.__calc_flow(self.head, sz_heads)

    # If the flow from bore hole into SZ is positive we need to reduce the head to reach 0 flow - and vice versa.
    # This should be ok even when going to/below the bottom of bore hole because this means SZ head is at/below bottom
    # of the bore hole and the resulting flow is 0. In this case we have a valid bracket: f(flow_sum_1) == 0
    if flow_sum_0 > 0:
      add = -1
    else:
      add = 1

    # Find bracket, i. e. 2 starting points where sign of flow differs (required for scipy brent)
    while True:
      add *= 2
      if abs(add) >= 1024: # A kilometer of head change and bracket still not found??? Something must be wrong, but I don't know what!
        msg = f"Bore hole '{self.name}': Unable to find a pair of head values where the resulting flow is pos for one and neg for the other."
        raise ValueError(msg)
      flow_sum_1 = self.__calc_flow(self.head + add, sz_heads)
      if (flow_sum_0 * flow_sum_1) <= 0: # bracket found?
        break
      self.head = self.head + add

    self.head, r = optimize.brentq(self.__calc_flow, self.head, self.head + add, rtol=EPS, args = sz_heads, maxiter = MAX_ITER, full_output=True, disp=False)
    self.iterations.append(self.its)

    flow_sum_0 = self.__calc_flow(self.head, sz_heads) # don't worry about performance doing this again, it's cached!
    ms.wm.print(f"accepted:  {self.head:.9f} {flow_sum_0 * L_PER_M3:.4g}")
    if not r.converged:
      self.no_converge += 1
      ms.wm.print(f"        *NO CONVERGENCE!*")
    ms.wm.print(f"--------------------------------------\r\n")

  def update_flux_head(self, sz_heads, sz_time):
    try:
      # Activate either one of these solution methods:
      # self.__closed_solution(sz_heads, sz_time)
      self.__iterative_solution(sz_heads, sz_time)
      self.heads.append(self.head)
    except:
      l = 1
      ms.wm.print(f"bh_head \"{self.name}\"")
      for head, flux in zip(sz_heads[self.i, self.j], sz_source[self.i, self.j]):
        ms.wm.print(f"    GW layer {l:3} - head: {head:10.3f} m, flux {flux * L_PER_M3:10.3g} l/s")
        l += 1
      raise

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

  # set floating point exception handling to raise an error for overflow, 0-divide and invalid operations
  np.seterr(divide='raise', over='raise', invalid='raise')

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
        j, k = ms.wm.gridCoordToCell(bh_x, bh_y)
        if not ms.wm.gridIsInternal(j, k):
          msg = f"Cannot create borehole \"{bh_name}\" outside the model area!\n"\
                f"  x={bh_x}, y={bh_y};  j={j}, k={k}"
          raise ValueError(msg)
        leaks = [khl / (FLOW_DISTANCE_CELL_RATIO * dx) for khl in kh[j, k]] # leakages in this SZ column
        bh = BoreHole(bh_name, bh_x, bh_y, bh_d, bh_bot, bh_top, zlyg[j,k], leaks, sz_heads)
        ms.wm.print(bh)
        bhs.append(bh)

def preTimeStep():
  _, sz_time, sz_heads = ms.wm.getValues(ms.paramTypes.SZ_HEAD) # sz_time and heads are from the end of the previous time step!
  if(sz_time is None): # Not an SZ time step
    return

  ms.wm.print(f"\r\n\r\n########### {sz_time} ########################")
  sz_source.value(0.0)
  for bh in bhs:
    bh.flow_cache.clear()
    bh.update_flux_head(sz_heads, sz_time)
    # ms.wm.print("")
  times.append(sz_time)
  ms.wm.setValues(sz_source)

def preLeaveSimulator():
  preTimeStep() # capture end of last time step
  title = "Open bore hole heads"
  items = [mikeio.ItemInfo(f"Head \"{bh.name}\"", itemtype=mikeio.EUMType.Water_Level) for bh in bhs]
  fname = os.path.join(setup_dir, f"{setup_name}.she - Result Files/{setup_name}_BoreHoleHeads.dfs0")

  ds = mikeio.Dataset(
    data = [bh.heads for bh in bhs],
    time = times,
    items = items
  )
  ds.to_dfs(fname, title=title)

  ms.wm.print("\r\nBore Hole Convergence Report\r\n============================")
  ms.wm.print("(Iteration count includes finding bracket and discarded iterations!)")
  for bh in bhs:
    ms.wm.print(f"Bore hole '{bh.name}'")
    itrs = np.array(bh.iterations)
    ms.wm.print(f"  Solution above threshold count:   {bh.no_converge} ({bh.no_converge / len(itrs):.2%})")
    ms.wm.print(f"  Average number of iterations:     {itrs.mean():.2}")
    ms.wm.print(f"  Maximum number of iterations:     {itrs.max ()}")
    ms.wm.print(f"  Standard deviation of iterations: {itrs.std ():.2}")
    ms.wm.print(f"  Total number of iterations:       {itrs.sum ()}")
  if any(bh.no_converge > 0 for bh in bhs):
    ms.wm.log("\r\nBore hole plugin: Solution threshold exceeded, see print log file for more details.")
