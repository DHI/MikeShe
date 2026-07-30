########## MIKE SHE plugin: SZ head below layer bottom diagnostic ##########
# Subject:      Computes the difference between the minimum head elevation in the first SZ layer
#               (over the full simulation) and the preprocessed bottom elevation of that layer:
#
#                 Result = Min(head elevation in SZ) - Lower level of first computational layer
#
#               Negative values indicate cells where the head temporarily fell below the bottom
#               of the first SZ layer during the simulation. This is a known numerical issue in
#               MIKE SHE that can cause mass balance errors and should be minimized by adjusting
#               layer geometry, time steps, or boundary conditions.
#               Two output items are written:
#                 1. Min head - Layer bottom:  the raw difference (m); negative = problem
#                 2. Min head:                 the minimum head over the simulation, for reference
# Usage:        Reference this file as a plugin in the MIKE SHE GUI and run the simulation.
#               Output is written to the model result folder as "SZ_HeadBelowBottom.dfs2".
# Dependencies: mikeio (which requires numpy - also used here directly)
# author:
# date:
############################################################################

import os
import mikeio
import MShePy as ms
import numpy as np


# Name of the SZ head item as written by MIKE SHE into _3DSZ.dfs3.
item_name_sz_head   = "head elevation in saturated zone"           # (m), simulation result

# Name of the layer bottom item in the preprocessed SZ file.
# Contains the lower boundary elevation of each computational SZ layer.
item_name_sz_bottom = "Lower level of computational layers in the saturated zone"  # (m), static

# globals
she_path = None  # full path to the .she file, retrieved at initialization


def postEnterSimulator():
  # Capture the .she file path early, before any output files exist.
  # leaveSimulator() needs it to construct the result file paths.
  global she_path
  she_path = ms.wm.getSheFilePath()


def leaveSimulator():
  # Runs after all output files have been closed, so all result files are fully written.

  # Construct paths to both the simulation result and the preprocessed file,
  # using MIKE SHE's default result folder convention.
  _, she_file      = os.path.split(she_path)
  she_base_name, _ = os.path.splitext(she_file)
  she_res_dir      = she_path + " - Result Files"
  sz_res_file      = os.path.join(she_res_dir, she_base_name + "_3DSZ.dfs3")
  preproc_file     = os.path.join(she_res_dir, she_base_name + "_PreProcessed_3DSZ.DFS3")
  os.chdir(she_res_dir)

  # --- Read minimum head over the full simulation ---
  ds     = mikeio.open(sz_res_file).read(items=[item_name_sz_head])
  item   = ds[item_name_sz_head]   # shape: (time, z, y, x) for a true 3D SZ output
  values = item.values

  # We are interested in the first SZ layer only, as that is the layer whose bottom
  # we compare against. Extract layer index 0 along the z-axis before any further processing
  # so that nan_mask and min_values are consistently 2D spatially (y, x).
  layer0_values = values[:, 0, :, :]   # (time, y, x)

  # Build a NaN mask from the first time step of layer 0 to identify inactive cells.
  # MIKE SHE writes NaN for cells outside the model domain or not part of the SZ computation.
  nan_mask = np.isnan(layer0_values[0])   # (y, x)

  # Minimum head in layer 0 over all time steps - this is the worst-case head level
  # in the first layer, which is the relevant quantity for the below-bottom diagnostic.
  min_values              = np.nanmin(layer0_values, axis=0, keepdims=True)  # (1, y, x)
  min_values[0, nan_mask] = np.nan

  # --- Read preprocessed bottom elevation of the first SZ layer ---
  preproc_ds   = mikeio.open(preproc_file).read(items=[item_name_sz_bottom])
  bottom_item  = preproc_ds[item_name_sz_bottom]

  # The preprocessed file contains all SZ layers as a dfs3.
  # Extract the first layer (index 0) and the first (only) time step to get a 2D (y, x) array.
  bottom_values               = bottom_item.isel(z=0).values[0]    # (y, x)
  bottom_values[nan_mask]     = np.nan
  bottom_values               = bottom_values[np.newaxis, ...]      # (1, y, x)

  # --- Compute difference: Min head - Layer bottom ---
  # Positive: head stayed above the layer bottom throughout the simulation (healthy)
  # Negative: head fell below the layer bottom at some point (problematic)
  diff_values              = min_values - bottom_values             # (1, y, x)
  diff_values[0, nan_mask] = np.nan

  # --- Build output items ---

  # 2D geometry for layer 0, used for all output DataArrays
  layer0_geometry = item.isel(z=0).geometry

  # Item 1: the raw difference - main diagnostic output
  diff_item = mikeio.DataArray(
    diff_values,
    time=[ds.time[0]],
    item=mikeio.ItemInfo("Min SZ head - Layer 1 bottom (negative = head below bottom)", itemtype=item.item.type, unit=item.item.unit),
    geometry=layer0_geometry
  )

  # Item 2: minimum head for reference, so the result file is self-contained
  min_item = mikeio.DataArray(
    min_values,
    time=[ds.time[0]],
    item=mikeio.ItemInfo(item_name_sz_head + " (Min)", itemtype=item.item.type, unit=item.item.unit),
    geometry=layer0_geometry
  )

  # Write both items as a single-timestep dfs2 file for easy inspection in MIKE ZERO / MIKE VIEW.
  ds_result = mikeio.Dataset([diff_item, min_item])
  ds_result.to_dfs("SZ_HeadBelowBottom.dfs2")
