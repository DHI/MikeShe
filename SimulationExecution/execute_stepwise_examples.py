# Subject:      This is a python script demonstrating how to run a MIKE She model using the MShePy API taking one or more time steps.
# Usage:        Execute as a standalone script. Current working directory must be the directory where this script is located.
# Dependencies: mikeio, MIKE Zero
#               In order to make the MShePy module (MIKE Zero) accessible either
#                 - add the MIKE Zero installation path/bin/x64 directory to the system/user variable PYTHONPATH (preferred) or
#                 - call sys.path.append("MIKE Zero installation path/bin/x64")
# Data:         DEMO_MODEL (see variable below)
# \author dhi\uha
# \date 06/2025

from datetime import datetime, timedelta
import mikeio
import MShePy as ms
import os
import subprocess as sp
import sys
DEMO_MODEL = "Data\\3x3_Box\\3x3_Box.she"


#######################################
# Helper functions
#######################################

def pp(setup):
  """Run the MIKE She preprocessor.
  As no python interface is available use the executable and run it as an external process.
  """
  # Get the directory where MShePy was loaded from and use PP from the same installation
  mz_dir = os.path.dirname(ms.__file__)
  pp_exe = os.path.join(mz_dir, "MShe_Preprocessor.exe")
  sp.run([pp_exe, setup])


def enable_plugin(in_path, out_path):
  """Enable using plugins, set the path to the python interpreter and the
  path to this file for plugins - save as copy.
  """
  # Use the same interpreter that is running this script - but we need the dll, not the exe
  py_dir = os.path.dirname(sys.executable) # path to python.exe
  py_dll = f"python{sys.version_info.major}{sys.version_info.minor}.dll"
  py_path = os.path.join(py_dir, py_dll)

  pfs = mikeio.PfsDocument(in_path, unique_keywords=False)
  pfs.MIKESHE_FLOWMODEL.SimSpec.ModelComp.Plugins = 1
  pfs.MIKESHE_FLOWMODEL.Plugins.PyResolve = 1
  # Set paths - special syntax for pfs paths using '|'
  pfs.MIKESHE_FLOWMODEL.Plugins.PyPath = f"|{py_path}|"
  pfs.MIKESHE_FLOWMODEL.Plugins.PluginFileList.PluginFile_1.FILE_NAME = f"|{__file__}|" # this file
  pfs.write(out_path) # save copy


#######################################
# Plugins
#######################################

def postTimeStep():
  """A MIKE She plugin function. No technical problem to put it in the same file as the code
  calling the MIKE She engine, however in a larger project it might be cleaner to separate it.
  """
  ms.wm.log(f"Message from plugin: Time step performed, time now: {ms.wm.currentTime()}")

#######################################
# Demo functions
#######################################

def exectue_time_steps(setup):
  """Demo: Run the MIKE She water movement engine by taking time steps via the python API."""
  p, e = os.path.splitext(setup)
  setup_plugin = p + "_plugin" + e
  # optional: Enable plugins
  enable_plugin(setup, setup_plugin)

  pp(setup_plugin)
  ms.wm.initialize(setup_plugin)
  performed, dt_hours, first_time = ms.wm.performTimeStep()
  if performed:
    ms.wm.log("Initial time step performed!")
    ms.wm.log(f"  Duration: {dt_hours} h")
    ms.wm.log(f"  New time: {first_time}")
  else:
    ms.wm.log("Failed to execute time step via python")
    return
  dt = timedelta(hours=dt_hours * 3.5)
  next_time = first_time + dt
  ms.wm.log(f"Requesting to run to: {next_time}")
  next_time = ms.wm.runToTime(next_time)
  ms.wm.log(f"Actual new date-time after calling runToTime: {next_time}")
  ms.wm.terminate(True)


if __name__ == "__main__":
  exectue_time_steps(DEMO_MODEL)
