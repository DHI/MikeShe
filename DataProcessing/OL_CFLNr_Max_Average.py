########## MIKE SHE plugin: Overland flow Courant number post-processing ##########
# Subject:      Reads the overland flow debug output written by MIKE SHE (_overland.dfs2) after
#               a simulation and collapses the full time series of the two OL Wave Courant number
#               items into static 2D grids by reducing over all time steps. For each input item,
#               two output grids are produced:
#                 - Mean: average Courant number per cell over the full simulation period
#                 - Max:  worst-case Courant number observed in any single time step
#               This yields 4 output items in total, covering both typical and peak conditions
#               for both the mean and max Courant numbers reported by MIKE SHE.
#               Courant numbers above 1.0 indicate cells where the explicit OL solver may become
#               unstable; the static maps help identify persistent and transient problem areas.
# Usage:        Reference this file as a plugin in the MIKE SHE GUI and run the simulation.
#               Output is written to the model result folder as "OL_Debug.dfs2", containing
#               all 4 items as a single-timestep dfs2 file.
#               Requires overland flow and OL debug output to be enabled in the MIKE SHE
#               simulation settings.
# Dependencies: mikeio (which requires numpy - also used here directly)
# author:
# date:
###################################################################################

import os
import mikeio
import MShePy as ms
import numpy as np


# Names of the two OL Courant number items as written by MIKE SHE into _overland.dfs2.
# Both are indicators of numerical stability in the explicit overland flow solver:
#   Mean: average Courant number across all OL sub-steps within a WM time step
#   Max:  maximum Courant number across all OL sub-steps within a WM time step
item_name_ol_mean = "Mean OL Wave Courant number (explicit OL)"  # persistent instability indicator
item_name_ol_max  = "Max OL Wave Courant number (explicit OL)"   # peak / worst-case instability indicator
item_names = [item_name_ol_mean, item_name_ol_max]

# globals
she_path = None  # full path to the .she file, retrieved at initialization


def postEnterSimulator():
  # Capture the .she file path early, before any output files exist.
  # leaveSimulator() needs it to construct the result file path.
  global she_path
  she_path = ms.wm.getSheFilePath()


def leaveSimulator():
  # Runs after all output files have been closed, so the dfs2 result file is fully written at this point.

  # Construct the path to the OL debug result file using MIKE SHE's default result folder convention.
  _, she_file      = os.path.split(she_path)
  she_base_name, _ = os.path.splitext(she_file)
  she_res_dir      = she_path + " - Result Files"
  ol_res_file      = os.path.join(she_res_dir, she_base_name + "_overland.dfs2")
  os.chdir(she_res_dir)

  ds = mikeio.open(ol_res_file).read(items=item_names)

  # Build a NaN mask from the first item's first time step to identify inactive cells.
  # MIKE SHE writes NaN for cells outside the model domain or not part of the OL computation.
  # We apply the same mask to all output items so inactive cells stay NaN (not 0) in the result.
  nan_mask = np.isnan(ds[item_name_ol_mean][0].values)

  # For each input item, compute both Mean and Max over all time steps.
  # This produces 4 output items in total, ordered as: Mean(Mean), Max(Mean), Mean(Max), Max(Max).
  output_items = []
  for name in item_names:
    item = ds[name]

    # --- Mean reduction ---
    # Average over all time steps. Represents the typical Courant number per cell across the simulation.
    # nanmean ignores NaN (inactive cells), so only active cells contribute to the average.
    mean_values = np.nanmean(item.values, axis=0, keepdims=True)  # collapse time axis -> [1, ny, nx]
    mean_values[0, nan_mask] = np.nan  # restore inactive cells
    output_items.append(mikeio.DataArray(
      mean_values,
      time=[ds.time[0]],
      item=mikeio.ItemInfo(name + " (Mean)", itemtype=item.item.type, unit=item.item.unit),
      geometry=item.geometry
    ))

    # --- Max reduction ---
    # Maximum over all time steps. Identifies cells where peak Courant numbers occurred,
    # even if only in a single time step - useful for spotting transient instability hotspots.
    # nanmax ignores NaN (inactive cells), so only active cells contribute to the maximum.
    max_values = np.nanmax(item.values, axis=0, keepdims=True)    # collapse time axis -> [1, ny, nx]
    max_values[0, nan_mask] = np.nan  # restore inactive cells
    output_items.append(mikeio.DataArray(
      max_values,
      time=[ds.time[0]],
      item=mikeio.ItemInfo(name + " (Max)", itemtype=item.item.type, unit=item.item.unit),
      geometry=item.geometry
    ))

  # Write all 4 items as a single-timestep dfs2 file for easy inspection in MIKE ZERO / MIKE VIEW.
  ds_result = mikeio.Dataset(output_items)
  ds_result.to_dfs("OL_Debug.dfs2")
