"""Sample MIKE-SHE plugin for the WM-WQ coupling: the transfers of one setup. The coupling itself is
wm_wq_coupling.py next to this file; see its header for what it is and is not.

Usage:
  - copy this file, wm_wq_coupling.py and the two wrappers (wq_wrapper.py, wm_wrapper.py) into one folder
  - replace the transfers below by those of your setup
  - attach this file as plugin to the MIKE-SHE model."""

import os
import sys

plugin_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, plugin_dir)
import wm_wq_coupling as wmwq
from wm_wq_coupling import Transfer, MM_PER_M, S_PER_DAY, postEnterSimulator, postTimeStep, leaveSimulator

wmwq.LOG_LEVEL = 2 # 1-4 currently used

"""List of Transfers from WM to WQ, project specific"""
wmwq.transfers_wm2wq = [
  Transfer(
    source_item = "reference evapotranspiration (OpenMI)",
    source_elem = "BaseGrid",
    target_item = "RefEt: UZ fixed(undef) (matrix phase)",
    target_elem = "WMUZ3DGrid",
    fac = MM_PER_M * S_PER_DAY # convert "back" from m/s to mm/d
  ),
  Transfer(
    source_item = "actual evapotranspiration",
    source_elem = "BaseGrid",
    target_item = "EtA: UZ fixed(undef) (matrix phase)",
    target_elem = "WMUZ3DGrid",
    fac = MM_PER_M * S_PER_DAY # convert "back" from m/s to mm/d
  ),
  Transfer(
    source_item = "root water uptake",
    source_elem = "WMUZ3DGrid",
    target_item = "RootUp: UZ fixed(undef) (matrix phase)",
    target_elem = "WMUZ3DGrid",
    fac = MM_PER_M * S_PER_DAY # convert "back" from m/s to mm/d
  ),
  Transfer(
    source_item = "crop coefficient",
    source_elem = "BaseGrid",
    target_item = "Kc: UZ fixed(undef) (matrix phase)",
    target_elem = "WMUZ3DGrid",
    fac = 1
  )
]


"""List of Transfers from WQ to WM, project specific"""
wmwq.transfers_wq2wm = [
  Transfer( # LAI
    source_item = "LAI: UZ fixed(undef) (matrix phase)",
    source_elem = "WMUZ3DGrid",
    target_item = 12, # mshe "new" datatype, found in manual. We cannot use the enum ms.paramTypes.[...] before initializing the model.
    target_elem = "BaseGrid"
  ),
  Transfer( # Root depth (hard coded index 1 for this one in wm_wq_coupling.py!)
   source_item  = "RD: UZ fixed(undef) (matrix phase)",
   source_elem  = "WMUZ3DGrid",
   target_item  = 11, # mshe "new" datatype, found in manual. We cannot use the enum ms.paramTypes.[...] before initializing the model.
   target_elem  = "BaseGrid"
  )
]

# for now hard code a few names, may later treat those by a naming convention
wmwq.copydown_names   = ["LAI", "RD"] # items that shall be copied down from top to all layers (with same name)
