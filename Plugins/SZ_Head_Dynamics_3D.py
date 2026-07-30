########## MIKE SHE plugin: SZ head dynamics and turning point analysis ##########
# Subject:      Reads the SZ head results written by MIKE SHE (_3DSZ.dfs3) after a simulation
#               and computes six static 2D grids from the full transient time series of
#               "head elevation in saturated zone", written as one dfs2 file per SZ layer:
#                 1. Time mean:              average head level per cell over the full simulation
#                 2. Mean absolute change:   mean of |delta head| between consecutive time steps,
#                                            indicating how dynamically active each cell is
#                 3. Turning point count:    number of times the head reverses direction per cell
#                                            (i.e. switches from rising to falling or vice versa),
#                                            indicating oscillatory / unstable behaviour
#                 4. Max:                    maximum head level observed in any single time step
#                 5. Min:                    minimum head level observed in any single time step
#                 6. Range (Max - Min):      total variation range per cell over the simulation
#               Zero differences are skipped when detecting turning points: a timestep with no
#               change does not interrupt a trend (same logic as the reference xarray implementation).
#               Since the SZ head item is 3D (time, z, y, x), each SZ layer is processed
#               independently and written to a separate dfs2 file:
#               SZ_HeadDynamics_Layer1.dfs2, SZ_HeadDynamics_Layer2.dfs2, ...
# Usage:        Reference this file as a plugin in the MIKE SHE GUI and run the simulation.
#               Output is written to the model result folder, one dfs2 file per SZ layer.
# Dependencies: mikeio (which requires numpy - also used here directly)
# author:
# date:
##################################################################################

import os
import mikeio
import MShePy as ms
import numpy as np


# Name of the SZ head item as written by MIKE SHE into _3DSZ.dfs3.
item_name_sz_head = "head elevation in saturated zone"  # (m), 3D item across all SZ layers

# globals
she_path = None  # full path to the .she file, retrieved at initialization


# ---------- helper functions ----------

def _ffill_axis0(arr):
  # Forward fill NaN values along the time axis (axis 0).
  # Used to propagate the last known sign direction through zero-difference timesteps.
  result = arr.copy()
  for t in range(1, result.shape[0]):
    nan_mask = np.isnan(result[t])
    result[t] = np.where(nan_mask, result[t - 1], result[t])
  return result


def _bfill_axis0(arr):
  # Backward fill NaN values along the time axis (axis 0).
  # Used to handle leading zeros at the start of the sign series (before any direction is established).
  result = arr.copy()
  for t in range(result.shape[0] - 2, -1, -1):
    nan_mask = np.isnan(result[t])
    result[t] = np.where(nan_mask, result[t + 1], result[t])
  return result


def _count_turning_points(values):
  # Count the number of turning points (direction reversals) per cell across all time steps.
  # A turning point is where the head switches from rising to falling or vice versa.
  # Shape of values: (time, y, x) for a single layer.

  # Step 1: first difference -> direction of change at each timestep
  diff = np.diff(values, axis=0).astype(float)        # shape: (time-1, y, x)

  # Step 2: sign of each difference: +1 = rising, -1 = falling, 0 = no change
  sgn = np.sign(diff)

  # Step 3: replace zeros with NaN so that flat periods are skipped during fill
  sgn[sgn == 0] = np.nan

  # Step 4: forward- then backward-fill NaN along time axis so each cell always has
  #         a defined direction, carrying the last known trend through flat periods
  sgn = _ffill_axis0(sgn)
  sgn = _bfill_axis0(sgn)

  # Step 5: a turning point occurs where the sign flips between consecutive timesteps
  #         product < 0 means the two signs are opposite (one +1, one -1)
  sign_change = (sgn[1:] * sgn[:-1]) < 0              # shape: (time-2, y, x)

  # Step 6: count turning points per cell, ignoring NaN (inactive cells)
  return np.nansum(sign_change.astype(float), axis=0) # shape: (y, x)


# ---------- plugin slots ----------

def postEnterSimulator():
  # Capture the .she file path early, before any output files exist.
  # leaveSimulator() needs it to construct the result file path.
  global she_path
  she_path = ms.wm.getSheFilePath()


def leaveSimulator():
  # Runs after all output files have been closed, so the result file is fully written at this point.

  # Construct the path to the SZ result file using MIKE SHE's default result folder convention.
  _, she_file      = os.path.split(she_path)
  she_base_name, _ = os.path.splitext(she_file)
  she_res_dir      = she_path + " - Result Files"
  sz_res_file      = os.path.join(she_res_dir, she_base_name + "_3DSZ.dfs3")
  os.chdir(she_res_dir)

  ds     = mikeio.open(sz_res_file).read(items=[item_name_sz_head])
  item   = ds[item_name_sz_head]   # shape: (time, z, y, x) for a true 3D SZ output
  values = item.values
  nz     = values.shape[1]         # number of SZ layers

  # Process each SZ layer independently and write one dfs2 file per layer.
  # Extracting per-layer slices (time, y, x) ensures nan_mask is always 2D (y, x),
  # avoiding indexing errors that arise when operating on the full 4D array.
  for iz in range(nz):
    layer_values   = values[:, iz, :, :]      # (time, y, x) - one layer
    layer_geometry = item.isel(z=iz).geometry # 2D geometry for this layer

    # Build a NaN mask from the first time step to identify inactive cells in this layer.
    # MIKE SHE writes NaN for cells outside the model domain or below the aquifer bottom.
    nan_mask = np.isnan(layer_values[0])      # (y, x)

    # --- Item 1: Time mean ---
    # Average head level per cell over the full simulation period.
    mean_values              = np.nanmean(layer_values, axis=0, keepdims=True)       # (1, y, x)
    mean_values[0, nan_mask] = np.nan
    mean_item = mikeio.DataArray(
      mean_values,
      time=[ds.time[0]],
      item=mikeio.ItemInfo(item_name_sz_head + " (Mean)", itemtype=item.item.type, unit=item.item.unit),
      geometry=layer_geometry
    )

    # --- Item 2: Mean absolute change per timestep ---
    # Mean of |delta head| between consecutive timesteps. High values indicate strongly
    # dynamic cells; near-zero values indicate cells that barely change over the simulation.
    diff                     = np.diff(layer_values, axis=0)                          # (time-1, y, x)
    dma_values               = np.nanmean(np.abs(diff), axis=0, keepdims=True)       # (1, y, x)
    dma_values[0, nan_mask]  = np.nan
    dma_item = mikeio.DataArray(
      dma_values,
      time=[ds.time[0]],
      item=mikeio.ItemInfo(item_name_sz_head + " (Mean abs. change per timestep)", itemtype=item.item.type, unit=item.item.unit),
      geometry=layer_geometry
    )

    # --- Item 3: Turning point count ---
    # Number of times the head changes direction per cell. High values indicate oscillatory
    # behaviour which may point to numerical instability or strong forcing variability.
    count_values              = _count_turning_points(layer_values)                   # (y, x)
    count_values              = count_values[np.newaxis, ...]                         # (1, y, x)
    count_values[0, nan_mask] = np.nan
    count_item = mikeio.DataArray(
      count_values,
      time=[ds.time[0]],
      item=mikeio.ItemInfo(item_name_sz_head + " (Turning point count)", itemtype=item.item.type, unit=item.item.unit),
      geometry=layer_geometry
    )

    # --- Item 4: Max ---
    # Maximum head level observed in any single time step per cell.
    max_values              = np.nanmax(layer_values, axis=0, keepdims=True)          # (1, y, x)
    max_values[0, nan_mask] = np.nan
    max_item = mikeio.DataArray(
      max_values,
      time=[ds.time[0]],
      item=mikeio.ItemInfo(item_name_sz_head + " (Max)", itemtype=item.item.type, unit=item.item.unit),
      geometry=layer_geometry
    )

    # --- Item 5: Min ---
    # Minimum head level observed in any single time step per cell.
    min_values              = np.nanmin(layer_values, axis=0, keepdims=True)          # (1, y, x)
    min_values[0, nan_mask] = np.nan
    min_item = mikeio.DataArray(
      min_values,
      time=[ds.time[0]],
      item=mikeio.ItemInfo(item_name_sz_head + " (Min)", itemtype=item.item.type, unit=item.item.unit),
      geometry=layer_geometry
    )

    # --- Item 6: Range (Max - Min) ---
    # Difference between the maximum and minimum head level per cell.
    # Represents the total variation experienced by each cell over the simulation.
    range_values              = max_values - min_values                               # (1, y, x)
    range_values[0, nan_mask] = np.nan
    range_item = mikeio.DataArray(
      range_values,
      time=[ds.time[0]],
      item=mikeio.ItemInfo(item_name_sz_head + " (Range)", itemtype=item.item.type, unit=item.item.unit),
      geometry=layer_geometry
    )

    # Write all 6 items for this layer as a single-timestep dfs2 file.
    ds_result = mikeio.Dataset([mean_item, dma_item, count_item, max_item, min_item, range_item])
    ds_result.to_dfs(f"SZ_HeadDynamics_Layer{iz + 1}.dfs2")
