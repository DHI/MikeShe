########## MIKE SHE plugin, prototype
# Subject: Simulate open bore holes, calculate exchange flows between SZ layers and the bore hole and the resulting head in the bore hole
# Usage:   Attach to an MIKE SHE model as a plugin, make sure the hard coded values match
# author:  uha@dhigroup.com
# date:    10/2022

import MShePy as m
import math

pipeHead = 10000 # initial value at sim start

def preTimeStep():
  global pipeHead
  eps = 0.0000001 # stop criterion, tolerance for remaining flow sum <> 0 of all layers<->pipe
  maxIts = 10
  pipeHeadTmp = pipeHead
  leak = [0, 8e-4, 0, 8e-4] # pipe leakage factor (1/s). TODO: use SZ-conductivity?
  dzLy = [10, 9.5, 0.5, 80] # layer thickness (m). TODO: retrieve actual values per layer
  dPipe = 0.15  # borehole diameter (m)
  pipeBaseArea = math.pi * pow(dPipe / 2, 2)
  exArea = [d * math.pi * dPipe for d in dzLy]
  (startTime, endTime, szHeads) = m.wm.getValues(m.paramTypes.SZ_HEAD)
  m.wm.print("\r\n\r\n########### {0} ########################".format(endTime))
  szSource = m.dataset(m.paramTypes.SZ_SOURCE)
  dtHrs = m.wm.nextTimeStep()[0]
  m.wm.print("dt: {0} min".format(dtHrs * 60))
  m.wm.print("iteration         pipeHeadTest      pipeHeadTestLast  floSum            floSumLast        ")
  pipeHeadTest = pipeHead
  pipeHeadLast = pipeHead
  pipeHeadTestLast = float("NaN")
  floSumLast       = float("NaN")
  for its in range(maxIts):
    floSum = 0
    for i in range(4):
      if leak[i] == 0:
        continue
      flo = (pipeHeadTest - szHeads[2,2,i]) * exArea[i] * leak[i]
      szSource[2,2,i] = flo
      floSum += flo
      
    m.wm.print("{:10g} {:17.10f} {:17.10f} {:17.10f} {:17.10f}".format(its + 1, pipeHeadTest, pipeHeadTestLast, floSum, floSumLast))
    tmp = pipeHeadTest
    if its == 0:
      # First 2 iterations needed to get a starting point for secant method. Could use values from last time step instead.
      if floSum > 0:
        dh = -0.01
      else:
        dh = 0.01
      pipeHeadTest += dh
    else:
      if abs(floSum) < eps:
        break
      pipeHeadTest = pipeHeadTest - (pipeHeadTest - pipeHeadTestLast) / (floSum - floSumLast) * (floSum) # secant method
    pipeHeadTestLast = tmp
    floSumLast = floSum
    # end iterations

  pipeHead = pipeHeadTest

  m.wm.print("SZ Head LB: {:10.3f} m".format(szHeads[2,2,1]))
  m.wm.print("pipeHead:   {:10.3f} m".format(pipeHead))
  m.wm.print("SZ Head LD: {:10.3f} m".format(szHeads[2,2,3]))
  m.wm.print("")
  m.wm.setValues(szSource)
  