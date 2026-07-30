########## MIKE SHE plugin: SZ debug output post-processing ##########
# Subject:      Reads the SZ debug output written by MIKE SHE (_2DSZ.dfs2) after a simulation
#               and collapses the full time series of "Avg. no. of SZ iterations above threshold"
#               into two static 2D grids by reducing over all time steps:
#                 1. Mean: average per-time-step iteration count over the full simulation period
#                 2. Max:  worst-case iteration count observed in any single time step
#               Together these two maps show where the SZ solver consistently struggled (Mean)
#               and where peak numerical difficulty occurred (Max).
# Usage:        Reference this file as a plugin in the MIKE SHE GUI and run the simulation.
#               Output is written to the model result folder as "SZ_Debug.dfs2", containing
#               both items (Mean and Max) as a single-timestep dfs2 file.
#               Requires SZ debug output to be enabled in the MIKE SHE simulation settings.
# Dependencies: mikeio (which requires numpy - also used here directly)
# author:
# date:
######################################################################

import os
import mikeio
import MShePy as ms
import numpy as np


# Name of the SZ debug item as written by MIKE SHE into _2DSZ.dfs2.
# It counts how often the SZ solver exceeded its iteration threshold per cell and time step,
# and is an indicator of numerical difficulty in the saturated zone.
item_name_sz_its = "Avg. no. of SZ iterations above threshold"  # high values -> slow convergence / high compute cost

# globals
she_path = None  # full path to the .she file, retrieved at initialization


def postEnterSimulator():
  # Capture the .she file path early, before any output files exist.
  # leaveSimulator() needs it to construct the result file path.
  global she_path
  she_path = ms.wm.getSheFilePath()


def leaveSimulator():
  # Runs after all output files have been closed, so the dfs2 result file is fully written at this point.

  # Construct the path to the SZ debug result file using MIKE SHE's default result folder convention.
  _, she_file      = os.path.split(she_path)
  she_base_name, _ = os.path.splitext(she_file)
  she_res_dir      = she_path + " - Result Files"
  sz_res_file      = os.path.join(she_res_dir, she_base_name + "_2DSZ.dfs2")
  os.chdir(she_res_dir)

  ds   = mikeio.open(sz_res_file).read(items=[item_name_sz_its])
  item = ds[item_name_sz_its]

  # Build a NaN mask from the first time step to identify inactive cells.
  # MIKE SHE writes NaN for cells outside the model domain or not part of the SZ computation.
  # We apply the same mask to all output items so inactive cells stay NaN (not 0) in the result.
  nan_mask = np.isnan(item[0].values)

  # --- Item 1: Mean ---
  # Average over all time steps. Represents typical solver effort per cell across the simulation.
  # nanmean ignores NaN (inactive cells), so only active cells contribute to the average.
  mean_values = np.nanmean(item.values, axis=0, keepdims=True)  # collapse time axis -> [1, ny, nx]
  mean_values[0, nan_mask] = np.nan  # restore inactive cells
  mean_item = mikeio.DataArray(
    mean_values,
    time=[ds.time[0]],
    item=mikeio.ItemInfo(item_name_sz_its + " (Mean)", itemtype=item.item.type, unit=item.item.unit),
    geometry=item.geometry
  )

  # --- Item 2: Max ---
  # Maximum over all time steps. Identifies cells where peak numerical difficulty occurred,
  # even if only in a single time step - useful for spotting transient instability hotspots.
  # nanmax ignores NaN (inactive cells), so only active cells contribute to the maximum.
  max_values = np.nanmax(item.values, axis=0, keepdims=True)   # collapse time axis -> [1, ny, nx]
  max_values[0, nan_mask] = np.nan  # restore inactive cells
  max_item = mikeio.DataArray(
    max_values,
    time=[ds.time[0]],
    item=mikeio.ItemInfo(item_name_sz_its + " (Max)", itemtype=item.item.type, unit=item.item.unit),
    geometry=item.geometry
  )

  # Write both items as a single-timestep dfs2 file for easy inspection in MIKE ZERO / MIKE VIEW.
  ds_result = mikeio.Dataset([mean_item, max_item])
  ds_result.to_dfs("SZ_Debug.dfs2")
