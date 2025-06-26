# Subject:      This is a python script demonstrating how to run a MIKE She model using the MShePy API.
#               The simulation will be run as one without interaction between the time steps.
# Usage:        Execute as a standalone script. Current working directory must be the directory where this script is located.
# Dependencies: mikeio, MIKE Zero
#               In order to make the MShePy module (MIKE Zero) accessible either
#                 - add the MIKE Zero installation path/bin/x64 directory to the system/user variable PYTHONPATH (preferred) or
#                 - call sys.path.append("MIKE Zero installation path/bin/x64")
# Data:         DEMO_MODEL (see variable below)
# \author dhi\uha
# \date 06/2025

import shutil
import MShePy as ms
import mikeio
import os
import shutil
import subprocess as sp
from concurrent.futures import ThreadPoolExecutor, wait
from multiprocessing import Process
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


def wm(setup):
  """Run the MIKE She water movement engine via its python interface."""
  ms.wm.runAll(setup)


def pp_wm(setup):
  """Run the MIKE She preprocessor followed by the water movement engine."""
  pp(setup)
  wm(setup)


def make_prec_copy(in_path, out_path, new_prec):
  """Make a copy of a .she-file and set the value for global constant precipitation."""
  pfs = mikeio.PfsDocument(in_path, unique_keywords=False)
  pfs.MIKESHE_FLOWMODEL.Climate.PrecipitationRate.GLOBAL.FIXED_VALUE = new_prec
  pfs.write(out_path)


def prepare_and_run(setup, var, prec_val):
  path, extension = os.path.splitext(setup)
  setup_var = path + var + extension
  make_prec_copy(setup, setup_var, prec_val)
  print(f"starting {var}")
  proc = Process(target=pp_wm, args=(setup_var,))
  proc.start()
  # wait until this process has finished
  proc.join()


#######################################
# Demo functions
#######################################

def run_single(setup):
  """Demo: Run the MIKE She water movement engine for a single setup."""
  wm(setup)


def run_variants_serial(setup):
  """Demo: Run several variants of the same MIKE She simulation.
  Even though the simulations are run sequentially a process has to be spawned for each
  model run! This is because MIKE She is not reentrant, i. e. when run within a
  single process instances would share data and not be independant of each other.
  """
  p, e = os.path.splitext(setup)
  setup_v = p + "_Variant" + e
  shutil.copy(setup, setup_v)
  for p in [1,2,4]:
    pfs = mikeio.PfsDocument(setup_v, unique_keywords=False)
    pfs.MIKESHE_FLOWMODEL.Climate.PrecipitationRate.GLOBAL.FIXED_VALUE = p
    pfs.write(setup_v)
    proc = Process(target=pp_wm, args=(setup_v,))
    proc.start()
    # wait until this process has finished before starting the next one or returning
    proc.join()


def run_variants_parallel(setup):
  """Demo: Run several different MIKE She simulations in parallel.
  Note that it is not possible to run several instances of the same
  model (same .she file) in parallel, you need to make a copy.

  With many setups it will be more efficient to run multiple processes
  with single threaded models instead of using multithreading. The total
  number of processes should be controlled to avoid exhausting resources,
  see run_variants_parallel_pool demo.
  """
  path, extension = os.path.splitext(setup)
  procs = []
  for var, prec_val in [["_A", 1], ["_B", 2], ["_C", 4]]:
    setup_var = path + var + extension
    make_prec_copy(setup, setup_var, prec_val)
    proc = Process(target=pp_wm, args=(setup_var,))
    proc.start()
    procs.append(proc)

  # for each process wait until it has finished
  (p.join() for p in procs)


def run_variants_parallel_pool(setup):
  """Demo: Run several different MIKE She simulations in parallel.
  Limit the number of simulations run in parallel.
  Note that you cannot just use a process pool. Processes from the pool
  will be reused which does not work for the execution of MIKE She as
  it is not reentrant.
  Instead use a thread pool with the specified number of maximum threads (or simulations),
  and in each thread start a new process for each simulation.
  """
  path, extension = os.path.splitext(setup)
  job_info = [["_E", 3], ["_F", 8], ["_G", 10], ["_H", 12]]
  with ThreadPoolExecutor(max_workers=2) as executor:
    futures = [executor.submit(prepare_and_run, setup, var, prec_val) for var, prec_val in job_info]
    wait(futures)


if __name__ == "__main__":
  run_single(DEMO_MODEL)
  run_variants_serial(DEMO_MODEL)
  run_variants_parallel(DEMO_MODEL)
  run_variants_parallel_pool(DEMO_MODEL)
