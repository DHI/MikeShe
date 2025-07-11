# Subject: This is a MIKE SHE plugin that writes SZ boundary flows of the current simulation to a dfs2 file,
#          one item per SZ-layer, each simulation time step.
# Usage:   Reference this file as a plugin in the MIKE SHE GUI and run the simulation.
# Dependencies:
#   - mikecore (required for lower level dfs access than mikeio provides)
#   - Python in a version supported by MIKE SHE (last tested on MIKE 2026 with python 3.13)
# \author dhi\uha
# \date 04/2020


import numpy as np
import os
from pathlib import Path
import sys

import MShePy as ms
from mikecore.eum         import eumUnit, eumItem, eumQuantity
from mikecore.DfsFactory  import DfsBuilder, DfsFactory, DataValueType
from mikecore.DfsFile     import DfsSimpleType


class Cell:
  ix = 0 # location of this cell
  iy = 0 # location of this cell
  in_n = False # Is the cell in this direction inside the model?
  in_e = False # Is the cell in this direction inside the model?
  in_s = False # Is the cell in this direction inside the model?
  in_w = False # Is the cell in this direction inside the model?


# globals
layers = []
bnd_cells = []
she_path = ""
sim_start = 0
dfs = None
data_x = None
data_y = None
nx = 0
ny = 0
nz = 0


def setup_dfs0():
  global she_path
  global dfs
  global data_x
  global data_y
  global nx
  global ny
  global nz
  global sim_start
  sim_start, sim_end = ms.wm.simPeriod()
  she_path = ms.wm.getSheFilePath()
  setup_name = Path(she_path).stem
  setup_name, _ = os.path.splitext(setup_name)
  filename = os.path.join(she_path + " - Result Files", setup_name + "_SzBndFluxes.dfs2")
  builder = DfsBuilder.Create(filename, "MSHE SZ boundary fluxes output per layer", 0)
  builder.SetDataType(1)
  factory = DfsFactory()

  projection, x_origin, y_origin, rotation = ms.wm.gridGeoInfo()
  builder.SetGeographicalProjection(factory.CreateProjectionGeoOrigin(projection, x_origin, y_origin, rotation))
  sim_start = ms.wm.currentTime()
  # note: time unit given here has to be used in WriteItemTimeStepNext
  axis = factory.CreateTemporalNonEqCalendarAxis(eumUnit.eumUsec, sim_start)
  builder.SetTemporalAxis(axis)
  builder.DeleteValueFloat = -1e-30
  (_, _, values) = ms.wm.getValues(ms.paramTypes.SZ_X_FLO) # just for the geometry
  (nx, ny, nz) = values.shape()
  (x0, y0) = ms.wm.gridCellToCoord(0, 0)
  (x1, y1) = ms.wm.gridCellToCoord(1, 1)
  data_x = np.ndarray(nx * ny, "single") # flux in x-direction for each cell
  data_y = np.ndarray(nx * ny, "single") # flux in y-direction for each cell

  # - Set "delete" values outside model area, not to be touched again
  # - Set 0.0 for internal cells, not to be touched again
  # (only boundary cell values will be set in each time step)
  for x in range(nx):
    for y in range(ny):
      if(not ms.wm.gridIsInModel(x, y)):
        data_x[x + y * nx] = builder.DeleteValueFloat
        data_y[x + y * nx] = builder.DeleteValueFloat
      else:
        data_x[x + y * nx] = 0.0
        data_y[x + y * nx] = 0.0

  dx = x1 - x0  # cell size, dx == dy
  axis = factory.CreateAxisEqD2(eumUnit.eumUmeter, nx, x0 - dx / 2, dx, ny, y0 - dx / 2, dx)
  item_builder = builder.CreateDynamicItemBuilder()
  for iz in range(nz):
    for xy in ['x', 'y']:
      item_builder.Set(f'Boundary inflow layer {iz + 1}, {xy}-direction', eumQuantity.Create(eumItem.eumIDischarge, eumUnit.eumUm3PerSec), DfsSimpleType.Float)
      item_builder.SetValueType(DataValueType.MeanStepBackward)
      item_builder.SetAxis(axis)
      builder.AddDynamicItem(item_builder.GetDynamicItemInfo())

  builder.CreateFile(filename)
  dfs = builder.GetFile()


def postEnterSimulator():
  global layers
  global bnd_cells
  setup_dfs0()
  for x in range(nx):
    for y in range(ny):
      in_model = ms.wm.gridIsInModel(x, y) # includes bnd
      internal = ms.wm.gridIsInternal(x, y) # does not include bnd
      if(in_model and not internal): # on boundary?
        c = Cell()
        c.iX = x
        c.iY = y
        c.in_n = ms.wm.gridIsInternal(x,     y + 1)
        c.in_e = ms.wm.gridIsInternal(x + 1, y    )
        c.in_s = ms.wm.gridIsInternal(x,     y - 1)
        c.in_w = ms.wm.gridIsInternal(x - 1, y    )
        bnd_cells.append(c)

  postTimeStep() # to catch initial values


def postTimeStep():
  global data_x
  global data_y
  global nx
  global ny
  global nz
  global sim_start
  global bnd_cells
  (startTime, _, x_flows) = ms.wm.getValues(ms.paramTypes.SZ_X_FLO)
  if(startTime is None):
    return # not an SZ time step
  (_, _, y_flows) = ms.wm.getValues(ms.paramTypes.SZ_Y_FLO)
  now = ms.wm.currentTime()
  dfsTime = (now - sim_start).total_seconds() # use same time unit as in CreateTemporalNonEqCalendarAxis
  for z in range(nz):
    for c in bnd_cells:
      x_flow = 0
      y_flow = 0

      if(c.in_n):
        y_flow = y_flows[c.iX,      c.iY,     z] # flow in pos y direction is _into_   model
      if(c.in_e):
        x_flow = x_flows[c.iX,      c.iY,     z] # flow in pos x direction is _into_   model
      if(c.in_s):
        y_flow -= y_flows[c.iX,     c.iY - 1, z] # flow in pos y direction is _out_of_ model
      if(c.in_w):
        x_flow -= x_flows[c.iX - 1, c.iY,     z] # flow in pos x direction is _out_of_ model
      data_x[c.iX + c.iY * nx] = x_flow
      data_y[c.iX + c.iY * nx] = y_flow
    dfs.WriteItemTimeStepNext(dfsTime, data_x)
    dfs.WriteItemTimeStepNext(dfsTime, data_y)
  dfs.Flush()


def preLeaveSimulator():
  global dfs
  dfs.Close()
