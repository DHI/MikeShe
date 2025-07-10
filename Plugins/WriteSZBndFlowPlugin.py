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

import MShePy
from mikecore.eum         import eumUnit, eumItem, eumQuantity
from mikecore.DfsFactory  import DfsBuilder, DfsFactory, DataValueType
from mikecore.DfsFile     import DfsSimpleType


class Cell:
  ix = 0 # location of this cell
  iy = 0 # location of this cell
  inN = False # Is the cell in this direction inside the model?
  inE = False # Is the cell in this direction inside the model?
  inS = False # Is the cell in this direction inside the model?
  inW = False # Is the cell in this direction inside the model?


# globals
layers = []
bndCells = []
shePath = ""
simStart = 0
dfs = None
dfsDataX = None
dfsDataY = None
nX = 0
nY = 0
nZ = 0


def setupDfs0():
  global shePath
  global dfs
  global dfsDataX
  global dfsDataY
  global nX
  global nY
  global nZ
  global simStart
  sim_start, sim_end = MShePy.wm.simPeriod()
  shePath = MShePy.wm.getSheFilePath()
  sheDir = os.path.dirname(shePath)
  setup_name = Path(shePath).stem
  setup_name, _ = os.path.splitext(setup_name)
  filename = os.path.join(shePath + " - Result Files", setup_name + "_SzBndFluxes.dfs2")
  builder = DfsBuilder.Create(filename, "MSHE SZ boundary fluxes output per layer", 0)
  builder.SetDataType(1)
  factory = DfsFactory()

  projection, x_origin, y_origin, rotation = MShePy.wm.gridGeoInfo()
  builder.SetGeographicalProjection(factory.CreateProjectionGeoOrigin(projection, x_origin, y_origin, rotation))
  simStart = MShePy.wm.currentTime()
  # note: time unit given here has to be used in WriteItemTimeStepNext
  axis = factory.CreateTemporalNonEqCalendarAxis(eumUnit.eumUsec, simStart)
  builder.SetTemporalAxis(axis)
  builder.DeleteValueFloat = -1e-30
  (startTime, endTime, values) = MShePy.wm.getValues(MShePy.paramTypes.SZ_X_FLO) # just for the geometry
  (nX, nY, nZ) = values.shape()
  (x0, y0) = MShePy.wm.gridCellToCoord(0, 0)
  (x1, y1) = MShePy.wm.gridCellToCoord(1, 1)
  dfsDataX = np.ndarray(nX * nY, "single") # flux in x-direction for each cell
  dfsDataY = np.ndarray(nX * nY, "single") # flux in y-direction for each cell

  # - Set "delete" values outside model area, not to be touched again
  # - Set 0.0 for internal cells, not to be touched again
  # (only boundary cell values will be set in each time step)
  for x in range(nX):
    for y in range(nY):
      if(not MShePy.wm.gridIsInModel(x, y)):
        dfsDataX[x + y * nX] = builder.DeleteValueFloat
        dfsDataY[x + y * nX] = builder.DeleteValueFloat
      else:
        dfsDataX[x + y * nX] = 0.0
        dfsDataY[x + y * nX] = 0.0

  dx = x1 - x0  # cell size, dx == dy
  axis = factory.CreateAxisEqD2(eumUnit.eumUmeter, nX, x0 - dx / 2, dx, nY, y0 - dx / 2, dx)
  itemBuilder = builder.CreateDynamicItemBuilder()
  for iz in range(nZ):
    for xy in ['x', 'y']:
      itemBuilder.Set(f'Boundary inflow layer {iz + 1}, {xy}-direction', eumQuantity.Create(eumItem.eumIDischarge, eumUnit.eumUm3PerSec), DfsSimpleType.Float)
      itemBuilder.SetValueType(DataValueType.MeanStepBackward)
      itemBuilder.SetAxis(axis)
      builder.AddDynamicItem(itemBuilder.GetDynamicItemInfo())

  builder.CreateFile(filename)
  dfs = builder.GetFile()


def postEnterSimulator():
  global layers
  global bndCells
  setupDfs0()
  for x in range(nX):
    for y in range(nY):
      inModel = MShePy.wm.gridIsInModel(x, y) # includes bnd
      internal = MShePy.wm.gridIsInternal(x, y) # does not include bnd
      if(inModel and not internal): # on boundary?
        c = Cell()
        c.iX = x
        c.iY = y
        c.inN = MShePy.wm.gridIsInternal(x,     y + 1)
        c.inE = MShePy.wm.gridIsInternal(x + 1, y    )
        c.inS = MShePy.wm.gridIsInternal(x,     y - 1)
        c.inW = MShePy.wm.gridIsInternal(x - 1, y    )
        bndCells.append(c)

  postTimeStep() # to catch initial values


def postTimeStep():
  global dfsDataX
  global dfsDataY
  global nX
  global nY
  global nZ
  global simStart
  global bndCells
  (startTime, _, xFlows) = MShePy.wm.getValues(MShePy.paramTypes.SZ_X_FLO)
  if(startTime is None):
    return # not an SZ time step
  (_, _, yFlows) = MShePy.wm.getValues(MShePy.paramTypes.SZ_Y_FLO)
  now = MShePy.wm.currentTime()
  dfsTime = (now - simStart).total_seconds() # use same time unit as in CreateTemporalNonEqCalendarAxis
  for z in range(nZ):
    for c in bndCells:
      xFlow = 0
      yFlow = 0

      if(c.inN):
        yFlow = yFlows[c.iX,      c.iY,     z] # flow in pos y direction is _into_   model
      if(c.inE):
        xFlow = xFlows[c.iX,      c.iY,     z] # flow in pos x direction is _into_   model
      if(c.inS):
        yFlow -= yFlows[c.iX,     c.iY - 1, z] # flow in pos y direction is _out_of_ model
      if(c.inW):
        xFlow -= xFlows[c.iX - 1, c.iY,     z] # flow in pos x direction is _out_of_ model
      dfsDataX[c.iX + c.iY * nX] = xFlow
      dfsDataY[c.iX + c.iY * nX] = yFlow
    dfs.WriteItemTimeStepNext(dfsTime, dfsDataX)
    dfs.WriteItemTimeStepNext(dfsTime, dfsDataY)
  dfs.Flush()


def preLeaveSimulator():
  global dfs
  dfs.Close()
