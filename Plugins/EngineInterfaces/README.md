# Engine interfaces: water movement and water quality from a Python plugin

`MShePy` gives convenient Python access to the water movement (WM) simulation, but not to water quality (WQ).
Both engines also have C interfaces. The two modules in this folder wrap them with `ctypes`, so that a MIKE SHE
Python plugin can read and write WQ data, and WM data in the same way:

| Module | Engine library | Main functions |
|---|---|---|
| `wq_wrapper.py` | `MShe_wq.dll` | `init`, `openmi_init`, `get_input_exchange_item(_count)`, `get_output_exchange_item(_count)`, `openmi_prepare`, `get_values`, `set_values`, `perform_time_step`, `set_wm_coupling_time_h`, `terminate`, plus grid information |
| `wm_wrapper.py` | `MShe_wm.dll` | `get_values`, `get_element_count`, `is_internal`, `get_uz_dimensions`, `get_sz_layers`, `get_3d_gridinfo`, `get_wmuz_3d_gridinfo` |

Access through these wrappers is less convenient than through `MShePy`: items are addressed by their exchange item
ID and element set (for example `BaseGrid`, `SZ3DGrid`, `WMUZ2DGrid`, `WMUZ3DGrid`), and values come as flat arrays
in the engine's (Fortran) order. Use `MShePy` for water movement wherever it offers what you need.

## Requirements

- The modules load the engine libraries from the folder of the running executable, so they only work inside a
  simulation, i.e. imported by a plugin that MIKE SHE runs. They do not work from `python.exe`.
- The engine libraries must be from the same MIKE SHE release as the setup.

## Use

Put the two modules next to your plugin and import them:

```python
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wq_wrapper as wq
import wm_wrapper as wm
```

The available WQ items are listed by `get_input_exchange_item_count()` / `get_input_exchange_item(i)` and the
corresponding output functions; each item returns its ID, element set, description and EUM type and unit.
Pass the IDs you use, joined by `;`, to `openmi_prepare` before the first `get_values` / `set_values`.

A WQ simulation driven from a plugin follows the hooks of the plugin:

- `postEnterSimulator`: `wq.init(she_path, wm_coup=True)`, `wq.openmi_init()`, select the items, `wq.openmi_prepare(...)`
- `postTimeStep`: exchange values, `wq.set_wm_coupling_time_h(t)`, then `wq.perform_time_step()` until WQ has reached the WM time
- `leaveSimulator`: `wq.terminate(success)`
