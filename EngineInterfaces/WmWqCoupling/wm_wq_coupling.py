"""MIKE-SHE plugin to enable WM-WQ coupling.
Usage:
  - Set any extra transfer items that are required for your model
  - change LOG_LEVEL below if desired (0 for no logging, 3 for max)
  - attach this plugin to a MIKE-SHE model.
Transfer wm->wq:
  - Find a WM result item that is available as output via the python api
  - Create a WQ species with a telling name
  - Add data to extra_items list in the code below
  - ElementSets: Either equal or WM "BaseGrid" to WQ "WMUZ3DGrid"
"""

import sys
import os

# make sure the wm- and wq- python wrappers can be found, assuming they are in the same directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wq_wrapper as wq
import wm_wrapper as wm
import numpy as np

if "d.exe" in sys.executable:
  import MShePyd as ms
else:
  import MShePy as ms

LOG_LEVEL = 2 # 1-3 currently used
items_needed = []
sim_start = None
total_h_wq = 0
ok = True
nuz3 = 0
lluz = 0
z_top_uz = None
z_bot_uz = None
S_PER_H = 3600
H_PER_DAY = 24
S_PER_DAY = S_PER_H * H_PER_DAY
MM_PER_M = 1000


class Transfer:
  """WM-WQ item transfer definition, project specific"""
  def __init__(self, source_item: str, source_elem: str, target_item: str, target_elem: str, fac = 1):
    self.source_item = source_item
    self.source_elem = source_elem
    self.target_item = target_item
    self.target_elem = target_elem
    self.conversion_factor = fac

  def __str__(self):
    return f"{self.source_item} => {self.target_item}"

  def __repr__(self):
    return self.__str__()

"""List of Transfers from WM to WQ, project specific"""
transfers_wm2wq = []


"""List of Transfers from WQ to WM, project specific"""
transfers_wq2wm = []

transfers_sums      = [] # items that shall be summed up over all layers, found by naming convention
transfers_conc_sums = [] # species concentrations that shall be used to calculate mass, then sum up over layers and copy result to each layer, found by naming convention

# for now hard code a few names, may later treat those by a naming convention
copydown_names   = [] # items that shall be copied down from top to all layers (with same name)


def log(l: int, *args):
  """Write to mshe log depending on current LOG_LEVEL"""
  if l <= LOG_LEVEL:
    ms.wm.log(*args)


def transfer_wm2wq():
  """Transfer values:
  - WM values that are needed for WQ to run (flows etc.)
  - ExtraItems that are needed for a specific setup"""
  ok = False
  # wm-wq-coupling items, can be transferred verbatim
  for item, el_set in items_needed:
    start, end, values = wm.get_values(item, el_set)
    wq.set_values(item, el_set, values)

  # wm-wq-coupling extra items specific to a setup
  for transfer in transfers_wm2wq:
    _, t, wm_values = wm.get_values(transfer.source_item, transfer.source_elem)
    log(3, transfer)
    if transfer.source_elem == "BaseGrid" and transfer.target_elem == "WMUZ3DGrid":
      _, _, wq_values = wq.get_values(transfer.target_item, transfer.target_elem)
      # wq_values are internal cells only while base grid is 2d including boundary cells
      n = 0
      for i_ref in range(len(wm_values)):
        if wm.is_internal(i_ref):
          wq_values[-1, n] = wm_values[i_ref] * transfer.conversion_factor # -1: top layer
          n += 1
      wq.set_values(transfer.target_item, transfer.target_elem, wq_values)
    elif transfer.source_elem == transfer.target_elem:
      wq.set_values(transfer.target_item, transfer.target_elem, [v * transfer.conversion_factor for v in wm_values])
      log(3, [v * transfer.conversion_factor for v in wm_values])
    else:
      raise Exception(f"Transfer of extra items from WM '{transfer.source_elem}' to WQ '{transfer.target_elem}' not implemented!")
  ok = True


def transfer_wq2wm():
  """Transfer values:
  - WQ to WM feedback values"""
  ok = False
  # LAI: transfer from "fixed" in WQ to LAI in WM and take care of the different element set
  #      (nUz does _not_ include boundary cells, nRef does!)
  for tr in transfers_wq2wm:
    tStart, tEnd, wq_values = wq.get_values(tr.source_item, tr.source_elem)
    if not uz_values_available(tStart, tr.source_item):
      return
    _, _, wm_values = ms.wm.getValues(tr.target_item)
    n = 0
    shp = wm_values.shape()

    # transfer values from ecolab to WM, skip boundary cells as they are not included in UZ WQ cells
    for ix in range(shp[0]):
      for iy in range(shp[1]):
        if ms.wm.gridIsInternal(ix, iy):
          wm_values[ix, iy] = wq_values[-1, n] # use top layer, counting from bottom
          n += 1

    ms.wm.setValues(wm_values)
  ok = True


def uz_values_available(t_start: float, item: str) -> bool:
  """WQ fills UZ items (fixed and concentration, element set WMUZ3DGrid) only after a WQ time step
  that included a UZ transport step. After an OL-only step (OL AD time step shorter than the UZ AD
  time step) get_values returns zeros and t_start = -1; those zeros must not be processed or written
  back, they would overwrite the UZ state."""
  if t_start >= 0:
    return True
  log(3, f"No UZ step in the last WQ step, '{item}' not available - skipping UZ transfers")
  return False


def transfer_wq2wq():
  """Transfer values:
  - WQ to WQ (aggregation)"""

  tr = transfers_wq2wm[1] # RD
  tStart, _, rd_3d = wq.get_values(tr.source_item, tr.source_elem)
  if not uz_values_available(tStart, tr.source_item):
    return
  rd_2d = rd_3d[-1, :]  # only calculated in top layer (?), makes no sense to have different values over depth, copy down not done yet

  # Calculate the mass from concentration and water content per layer, then sum up across all layers and set that sum in each layer
  log(3, "sum up mass from layers...")
  _, _, uz_wc = ms.wm.getValues(ms.paramTypes.UZ_WC)
  log(3, "UZ_WC")
  log(3, uz_wc)
  log(3, uz_wc.value())
  shp = uz_wc.shape()

  for transfer in transfers_conc_sums:
    log(3, transfer)
    # get values as numpy array
    _, _, wq_values = wq.get_values(transfer.source_item, transfer.source_elem)
    log(3, '\r\nwq_values before summing:\r\n', wq_values)
    i_uz = 0

    # each 3d cell, skip boundary cells as they are not included in UZ WQ cells
    for ix in range(shp[0]):
      for iy in range(shp[1]):
        if ms.wm.gridIsInternal(ix, iy):
          rd = rd_2d[i_uz]
          z_ground = z_top_uz[-1, i_uz]

          # any roots inside a layer?
          has_roots = np.fromiter((rd > z_ground - z_top_uz[z, i_uz] for z in range(shp[2])), dtype=bool)

          # for convenience extract this column in the same order (bottom up) as the other data
          uz_wc_col = np.fromiter(uz_wc[ix, iy, ::-1], dtype=float)

          # sum up over all layers with roots in them (as requested by customer, interpolation, just black and white)
          sum_z = (uz_wc_col[has_roots] * wq_values[has_roots, i_uz]).sum()
          # if rd > 0 and sum_z > 0:
          #   breakpoint()
          wq_values[ has_roots, i_uz] = sum_z
          wq_values[~has_roots, i_uz] = 0
          i_uz += 1
    log(3, '\r\nwq_values after summing:\r\n', wq_values)
    wq.set_values(transfer.target_item, transfer.target_elem, wq_values)

  log(3, "sum up values from layers...")
  # Sum up values of all layers and assign to all layers
  for transfer in transfers_sums:
    log(3, transfer)
    # get values as numpy array
    _, _, wq_values = wq.get_values(transfer.source_item, transfer.source_elem)
    log(3, wq_values, '\r\n\r\n')
    value_sum = wq_values
    wq_values[:] = [sum(col) for col in zip(*wq_values)]
    log(3, wq_values, '\r\n\r\n')
    wq.set_values(transfer.target_item, transfer.target_elem, wq_values)

  # copy values from top layer to all layers below (makes no sense if there is just 1 layer)
  if lluz < 2:
    return

  log(3, "copy down...")
  for name in copydown_names:
    log(3, name)
    # get values as numpy array
    tStart, tEnd, wq_values = wq.get_values(f"{name}: UZ fixed(undef) (matrix phase)", "WMUZ3DGrid")
    log(4, type(wq_values))
    log(3, wq_values, '\r\n\r\n')
    wq_values[:-1, :] = wq_values[-1, :]
    log(4, type(wq_values))
    log(3, wq_values, '\r\n\r\n')
    wq.set_values(f"{name}: UZ fixed(undef) (matrix phase)", "WMUZ3DGrid", wq_values)


def postEnterSimulator():
  """After having intialized mshe-WM: Initialize WQ model and WM-WQ coupling"""
  global items_needed, sim_start, ok, nuz3, lluz, ol_drain, z_top_uz, z_bot_uz
  ok = False
  # breakpoint()
  she_path = ms.wm.getSheFilePath()
  sim_start, _ = ms.wm.simPeriod()
  ok = wq.init(she_path, wm_coup=True)
  if not ok:
    raise Exception("Error initializing WQ")
  ok = wq.openmi_init()
  if not ok:
    raise Exception("Error initializing WQ OpenMi")
  n_items = wq.get_input_exchange_item_count()
  item_names_nos = set()
  item_names_s = set()

  ol_drain = ms.dataset(ms.paramTypes.OLDR_IN_FLX)
  ol_drain.unit("mm/day")

  log(2, "")
  log(2, "Identify items for exchange (plugin)...")
  # - find "needed" items (WM-WQ-coupling)
  # - find items that might need to receive a sum up and be copied to layers
  for i in range(n_items):
    item_id, element_set, _, _, _, _, needed = wq.get_input_exchange_item(i)
    log(2, f"Item: {item_id}")
    if needed:
      items_needed.append([item_id, element_set])
    if element_set == "WMUZ3DGrid":
      if item_id[:2].lower() == "s_":
        item_names_s.add(item_id)

  log(2, f"\r\nItems needed for WM-WQ-coupling: {items_needed}")
  log(2, f"\r\nItem candidates for SUM: {item_names_s}")

  # Find items that might have to be summed up
  n_items = wq.get_output_exchange_item_count()
  for i in range(n_items):
    item_id, element_set, _, _, _, _ = wq.get_output_exchange_item(i)
    if element_set == "WMUZ3DGrid":
      if item_id[:2].lower() != "s_":
        item_names_nos.add(item_id)

  log(2, f"\r\nItem candidates for SUM input: {item_names_nos}")

  for s_item in item_names_s:
    species_name = s_item[2:s_item.find(":")]
    found = any(s.startswith(f"{species_name}: ") for s in item_names_nos)
    log(3, f"{species_name}: in item_names_nos? ", f"{found}!")
    if found:
      item_name = next(s for s in item_names_nos if s.startswith(f"{species_name}: "))
      transfer = Transfer(
                   source_item = item_name,
                   source_elem = "WMUZ3DGrid",
                   target_item = s_item,
                   target_elem = "WMUZ3DGrid"
                 )
      if "UZ fixed(undef)" in item_name:
        transfers_sums.append(transfer)
      elif "UZ concentration" in item_name:
        transfers_conc_sums.append(transfer)

  log(2, f"\r\nResulting names for plain SUM: {transfers_sums}")
  log(2, f"\r\nResulting names for conc-SUM:  {transfers_conc_sums}")

  all_items = items_needed + [[e.source_item, e.source_elem] for e in transfers_sums] + [[e.source_item, e.source_elem] for e in transfers_conc_sums]
  log(2, f"\r\nAll exchange items:         {all_items}")

  item_list = ";".join(sublist[0] for sublist in all_items)
  log(2, f"\r\nAll exchange items as list: {item_list}")
  wq.openmi_prepare(item_list)
  nuz2 = wq.get_element_count("WMUZ2DGrid")
  nuz3 = wq.get_element_count("WMUZ3DGrid")
  lluz = int(nuz3 / nuz2)

  z_top_uz, z_bot_uz = wm.get_wmuz_3d_gridinfo()
  log(3, "\r\nZ-Top:")
  log(3, z_top_uz[:])
  log(3, "\r\nZ-Bot:")
  log(3, z_bot_uz[:])
  postTimeStep()
  ok = True


def postTimeStep():
  """After a WM time step: Transfer values between WM and WQ and take WQ time steps to catch up with WM"""
  global total_h_wq, ok
  ok = False
  current_wm = ms.wm.currentTime()
  total_h_wm = (current_wm - sim_start).total_seconds() / S_PER_H
  log(3, f"New WM time: {total_h_wm} h, {current_wm}")
  # Note: transferring wm data must not be done more than once per wm time step, otherwise wq gets confused about time stepping!
  wq.set_wm_coupling_time_h(total_h_wm)
  transfer_wm2wq()
  while total_h_wm > total_h_wq:
    log(3, f"Running WQ from {total_h_wq} h ...")
    transfer_wq2wq()
    time_step_h, total_h_wq = wq.perform_time_step()
    log(3, f"... finishing at {total_h_wq} h")
  transfer_wq2wm()
  ok = True


def leaveSimulator():
  """WM Simulation completed, stop WQ"""
  wq.terminate(ok)
