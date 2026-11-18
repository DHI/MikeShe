# WM-WQ coupling: running water quality from a water movement plugin

`wm_wq_coupling.py` steps a water quality (WQ) simulation along with the water movement (WM) simulation and
exchanges values between the two, through the wrappers in the parent folder.

**What this is and is not.** The module was developed for a specific setup, not as a general WM-WQ coupling
that works out of the box. It is intended as a starting point for developing other model specific couplings.
DHI does not support it, takes no responsibility for a coupling built from it or for the engines processing
it correctly in such a configuration. The documentation of how to create a setup may be incomplete and may
not name all requirements. The module and the engines may not guard against invalid use.

The original model included MIKE ECO Lab. The transfers inside WQ (the sum of a species over the UZ layers
into an `s_<name>` species, the root depth taken from the second WQ-to-WM transfer, the copy-down of the top
layer) are the conventions of that setup.

A plugin states the transfers of its setup on the module and re-exports the module's hooks, so that MIKE SHE
finds them in the plugin:

```python
import os, sys
sys.path.insert(0, "<path to this folder>")
import wm_wq_coupling as wmwq
from wm_wq_coupling import Transfer, MM_PER_M, S_PER_DAY, postEnterSimulator, postTimeStep, leaveSimulator

wmwq.LOG_LEVEL = 2 # 0 none, 3 everything
wmwq.transfers_wm2wq = [ # WM output item -> WQ species
  Transfer("actual evapotranspiration", "BaseGrid", "EtA: UZ fixed(undef) (matrix phase)", "WMUZ3DGrid",
           fac = MM_PER_M * S_PER_DAY)]
wmwq.transfers_wq2wm = [ # WQ species (top layer) -> WM item, by MIKE SHE parameter type number
  Transfer("LAI: UZ fixed(undef) (matrix phase)", "WMUZ3DGrid", 12, "BaseGrid")]
wmwq.copydown_names = ["LAI"] # species whose top-layer value is copied to the layers below
```

`sample_plugin.py` in this folder is a complete plugin of this kind, with the transfers of the setup the module
was developed for.

The module finds `wm_wrapper.py` and `wq_wrapper.py` next to itself or in its parent folder, so either copy this
folder together with the two wrappers into your setup, or keep the repository layout.

The WM engine must run in AD-coupled mode for this, otherwise it does not write the flow results WQ
reads (`.xfrf`). The switch is a flag file next to the setup: the path of the `.she` file with `.she`
replaced by `_AdCoup`, no extension, any or no content. For `MySetup.she` create an empty `MySetup_AdCoup`
in the same folder. Without it the plugin fails in `postEnterSimulator`.

The items the WM-WQ coupling itself needs are found and transferred automatically. A species named `s_<name>`
receives, in every UZ layer, the sum over the layers of the species `<name>`: a plain sum for a fixed item, a
water-content-weighted sum over the root zone for a concentration. The root depth is taken from the second
entry of `transfers_wq2wm`.
