########## MIKE SHE plugin: UZ debug output post-processing ##########
# Subject:      Reads the UZ debug output written by MIKE SHE (_2DUZ_UzCells.dfs2) after a simulation
#               and collapses the full time series into a single static 2D grid per item by averaging
#               over all time steps. The result shows spatially where the UZ solver struggled most
#               during the simulation (numerically expensive cells).
# Usage:        Reference this file as a plugin in the MIKE SHE GUI and run the simulation.
#               Output is written to the model result folder as "UZ_Debug.dfs2".
#               Requires "UZ debug output" to be enabled in the MIKE SHE simulation settings.
# Dependencies: mikeio (which requires numpy - also used here directly)
# author:
# date:
######################################################################

import os
import mikeio
import MShePy as ms
import numpy as np


# Names of the two UZ debug items as written by MIKE SHE into _2DUZ_UzCells.dfs2.
# These count how often the UZ solver had to reduce its time step or iterate per cell,
# and are indicators of numerical difficulty in the unsaturated zone.
item_name_uz_red = "UZ time step reduction count"  # high values -> solver instability
item_name_uz_its = "UZ iteration count"             # high values -> slow convergence / high compute cost
item_names = [item_name_uz_red, item_name_uz_its]

# globals
she_path = None  # full path to the .she file, retrieved at initialization


def postEnterSimulator():
  # Capture the .she file path early, before any output files exist.
  # leaveSimulator() needs it to construct the result file path.
  global she_path
  she_path = ms.wm.getSheFilePath()


def leaveSimulator():
  # Runs after all output files have been closed, so the dfs2 result file is fully written at this point.

  # Construct the path to the UZ debug result file using MIKE SHE's default result folder convention.
  _, she_file    = os.path.split(she_path)
  she_base_name, _ = os.path.splitext(she_file)
  she_res_dir    = she_path + " - Result Files"
  uz_res_file    = os.path.join(she_res_dir, she_base_name + "_2DUZ_UzCells.dfs2")
  os.chdir(she_res_dir)

  ds = mikeio.open(uz_res_file).read()

  # Build a NaN mask from the reduction count item to identify inactive/dry cells.
  # MIKE SHE writes NaN for cells outside the model domain or not part of the UZ computation.
  # We apply the same mask to all output items so inactive cells stay NaN (not 0) in the result.
  nan_mask = None
  if ds[item_name_uz_red] is not None:
    nan_mask = np.isnan(ds[item_name_uz_red][0].values)

  # For each requested item: average over all time steps to produce a single static 2D map.
  # nanmean ignores NaN (inactive cells), so only active cells contribute to the average.
  # The result represents the mean per-time-step count over the full simulation period.
  selected_items = []
  for name in item_names:
    try:
      item           = ds[name]
      averaged_values = np.nanmean(item.values, axis=0, keepdims=True)  # collapse time axis -> [1, ny, nx]
      if nan_mask is not None:
        averaged_values[0, nan_mask] = np.nan  # restore inactive cells
      selected_items.append(mikeio.DataArray(averaged_values, time=[ds.time[0]], item=item.item, geometry=item.geometry))
    except KeyError:
      print(f"Warning: Item '{name}' not found in input file.")

  if not selected_items:
    raise ValueError("None of the requested items were found in the input file.")

  # Write averaged result as a single-timestep dfs2 file for easy inspection in MIKE ZERO / MIKE VIEW.
  ds_averaged = mikeio.Dataset(selected_items)
  ds_averaged.to_dfs("UZ_Debug.dfs2")
