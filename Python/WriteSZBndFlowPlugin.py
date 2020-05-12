# Subject: This is a MIKE SHE plugin that writes SZ boundary flows of the current simulation to a dfs2 file,
#          one item per layer, each simulation time step.
# Usage:   
#   - Reference this file as a plugin in the MIKE SHE GUI and run the simulation.
# Limitations:
#   - The spatial reference cannot currently be retrieved from MIKE SHE, so the output dfs2 file will be set to "local coordinates".
#     If required the actual spatial reference has to be set manually after the simulation.
# Prerequisites:
#   - MIKE SDK
#   - Python in a version supported by MIKE SHE (tested on MIKE 2020 with python 3.7) (e.g. from https://www.python.org/downloads/)
#   - python package "pythonnet":
#     + open windows cmd (as administrator if python installed for all users)
#     + cd to the python installation folder/Scripts directory (e.g. "C:\Program Files\Python37\Scripts")
#     + type:
#       pip install pythonnet
# \author dhi\uha
# \date 04/2020

import sys
import os
import subprocess as sp
import MShePy

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
  import clr
  global simStart
  global simStart
  now = MShePy.wm.currentTime()
  clr.AddReference("DHI.Mike.Install, Version=1.0.0.0, Culture=neutral, PublicKeyToken=c513450b5d0bf0bf") # "fully qualified" name required!
  from DHI.Mike.Install import MikeImport, MikeProducts
  MikeImport.SetupLatest()
  clr.AddReference("DHI.Generic.MikeZero.DFS")
  clr.AddReference("DHI.Generic.MikeZero.EUM")
  clr.AddReference("System")
  import System
  from System import Array
  from DHI.Generic.MikeZero import eumUnit, eumItem, eumQuantity
  from DHI.Generic.MikeZero.DFS import DfsFactory, DfsBuilder, DfsSimpleType, DataValueType
  shePath = MShePy.wm.getSheFilePath()
  sheDir = os.path.dirname(shePath)
  filename = os.path.join(sheDir, 'BndFluxes.dfs2')
  builder = DfsBuilder.Create(filename, "MSHE SZ boundary fluxes output per layer", 0)
  builder.SetDataType(1)
  factory = DfsFactory()
  builder.SetGeographicalProjection(factory.CreateProjectionGeoOrigin("NON-UTM", 0, 0, 0))
  simStart = now
  nowSys = System.DateTime(now.year, now.month, now.day, now.hour, now.minute, now.second)
  # note: time unit given here has to be used in WriteItemTimeStepNext
  axis = factory.CreateTemporalNonEqCalendarAxis(eumUnit.eumUsec, nowSys)
  builder.SetTemporalAxis(axis)
  builder.DeleteValueFloat = -1e-30
  (startTime, endTime, values) = MShePy.wm.getValues(MShePy.paramTypes.SZ_X_FLO) # just for the geometry
  (nX, nY, nZ) = values.shape()
  (x0, y0) = MShePy.wm.gridCellToCoord(0, 0)
  (x1, y1) = MShePy.wm.gridCellToCoord(1, 1)
  dfsDataX = Array.CreateInstance(System.Single, nX * nY)
  dfsDataY = Array.CreateInstance(System.Single, nX * nY)
  for x in range(nX):
    for y in range(nY):
      if(not MShePy.wm.gridIsInModel(x, y)):
        dfsDataX[x + y * nX] = builder.DeleteValueFloat
        dfsDataY[x + y * nX] = builder.DeleteValueFloat
  dx = x1 - x0  # cell size, dx == dy
  axis = factory.CreateAxisEqD2(eumUnit.eumUmeter, nX, x0 - dx / 2, dx, nY, y0 - dx / 2, dx)
  itemBuilder = builder.CreateDynamicItemBuilder()
  itemBuilder.SetValueType(DataValueType.MeanStepBackward)
  itemBuilder.SetAxis(axis)
  for iz in range(nZ):
    for xy in ['x', 'y']:
      itemBuilder.Set('Boundary inflow layer {0}, {1}-direction'.format(iz + 1, xy), eumQuantity.Create(eumItem.eumIDischarge, eumUnit.eumUm3PerSec), DfsSimpleType.Float)
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
  (startTime, endTime, xFlows) = MShePy.wm.getValues(MShePy.paramTypes.SZ_X_FLO)
  if(startTime is None):
    return # not an SZ time step
  (startTime, endTime, yFlows) = MShePy.wm.getValues(MShePy.paramTypes.SZ_Y_FLO)
  now = MShePy.wm.currentTime()
  dfsTime = (now - simStart).total_seconds()
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

