# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A collection of **standalone Python scripts** for MIKE SHE (DHI's integrated hydrological modelling
system): engine plugins, data-processing utilities and MShePy API examples. There is no package, no
build step, no test suite and no dependency manifest — each `.py` file is self-contained and is either
executed directly or referenced as a plugin from a MIKE SHE `.she` setup.

Three directories, three different execution models:

| Directory | How the code runs |
|---|---|
| `Plugins/` | Loaded *by* the MIKE SHE engine during a simulation; the engine calls named module-level hook functions |
| `SimulationExecution/` | Standalone scripts that *drive* the MIKE SHE engine through `MShePy` |
| `DataProcessing/` | Standalone CLI/GUI tools that read or rewrite MIKE SHE files; no engine involved |

## Dependencies and environment

- **`MShePy`** ships with MIKE Zero, not pip. For standalone scripts it must be importable: add
  `<MIKE Zero install>/bin/x64` to `PYTHONPATH` (preferred) or `sys.path.append(...)`.
  Inside a plugin the engine has already made it available.
- Python version must be one supported by the installed MIKE release (MIKE 2026: Python 3.13).
  For plugins, the `.she` file points at the **`pythonXY.dll`**, not `python.exe`.
- pip packages used across scripts: `mikeio` (high-level dfs/pfs), `mikecore` (low-level dfs writing,
  used where mikeio is not granular enough), `numpy`, `scipy` (`openBoreHoles.py`), `pyshp`
  (`PtPathline*.py`), `tkinter` (optional, for `pfs_pack.py`'s file dialog).
- mikeio makes breaking API changes between majors and these scripts track the current release
  (3.0.1 at time of writing). Notably `Dataset` takes only DataArrays — the old
  `Dataset(data=…, time=…, items=…)` constructor is gone. Check signatures against the installed
  version rather than assuming, and record what you tested against in the file header.

## Running things

```powershell
# MShePy examples — cwd must be the script's own directory (DEMO_MODEL is a relative path)
cd SimulationExecution
python run_all_examples.py          # runAll: single, serial variants, parallel, pooled
python execute_stepwise_examples.py # initialize + performTimeStep + runToTime, with an inline plugin

# Pack a .she setup and everything it references
python DataProcessing/pfs_pack.py <setup.she> [-o OUT] [-a|-m] [-d] [-f]
# no args -> interactive (tkinter dialog); -d writes a staging dir instead of a zip
```

There is no lint config and no automated tests. Verification means running an actual MIKE SHE model:
preprocess, then run, then inspect the produced dfs file. The 3-cell demo setup in
`SimulationExecution/Data/3x3_Box/` is the fast smoke-test model.

The preprocessor has no Python API — run `MShe_Preprocessor.exe` as a subprocess, locating it via
`os.path.dirname(ms.__file__)` so it comes from the same MIKE installation as the loaded `MShePy`
(see `pp()` in both `SimulationExecution` scripts).

## Plugin architecture

A plugin is a module defining any subset of hook functions the engine calls at fixed points. Hooks
used in this repo, in call order:

- `postEnterSimulator()` — after engine init. Grid geometry and static parameters are available;
  result files are not. This is where setups read config, build state, and capture
  `ms.wm.getSheFilePath()` for later.
- `preTimeStep()` / `postTimeStep()` — per water-movement time step.
- `preLeaveSimulator()` — output files still **open**. Use for closing files you opened yourself
  (`write_sz_bnd_flow.py`) or writing dfs0 results accumulated in memory (`openBoreHoles.py`).
- `leaveSimulator()` — output files **closed**. This is the only safe place to read the engine's own
  result files with mikeio; all the post-processing plugins do their work here.

Cross-cutting conventions the plugins rely on:

- **Result-file paths** are derived, not configured: the folder is `<full path to .she> + " - Result Files"`
  (note: the `.she` extension stays in the folder name) and files are `<she stem>_3DSZ.dfs3`,
  `_2DSZ.dfs2`, `_overland.dfs2`, `_2DUZ_UzCells.dfs2`, `_PreProcessed_3DSZ.dfs3` (preprocessed static
  data such as layer bottoms). Custom result folders are **not** supported by these plugins.
- **Axis order differs between the two data sources.** MShePy datasets are `(x, y, z)`; dfs files read
  through mikeio are `(time, z, y, x)` and `z` runs bottom-up. Mixing them requires an explicit
  transpose/flip — see `read_layer_bottoms()` in `openBoreHoles.py`.
- **Inactive cells.** MIKE SHE writes NaN outside the model domain. The reduction plugins capture a
  `nan_mask` from the first time step, reduce with `np.nan*` functions, then write NaN back into the
  masked cells so inactive area does not become 0. Use `ms.wm.gridIsInModel()` / `gridIsInternal()`
  when working on the MShePy side instead.
- **Output.** Post-processing plugins emit a single-timestep dfs2 built from `mikeio.DataArray` /
  `mikeio.Dataset.to_dfs()`, inheriting `itemtype`, `unit` and `geometry` from the source item.
- **Feedback into the simulation** goes through `ms.wm.setValues(dataset)` after filling an
  `ms.dataset(ms.paramTypes...)` — `openBoreHoles.py` writes `SZ_SOURCE` this way.
- **Logging**: `ms.wm.log()` for messages the user should see in the MIKE SHE log; `ms.wm.print()` for
  verbose output to the print log file. Do not use bare `print()` in a plugin.

The `Plugins/` post-processing scripts (`OL_*`, `SZ_*`, `UZ_*`) are deliberate near-duplicates of one
another: same skeleton, different item names and reductions. When changing that skeleton, expect to
change it in several files rather than refactoring them into a shared module.

## Editing `.she` / pfs files programmatically

Use `mikeio.PfsDocument(path, unique_keywords=False)` — `unique_keywords=False` matters, MIKE SHE pfs
files legitimately repeat keywords. File-valued keywords must be written wrapped in pipes:
`f"|{path}|"`. `enable_plugin()` in `execute_stepwise_examples.py` shows the full set of keys needed to
turn plugins on (`SimSpec.ModelComp.Plugins`, `Plugins.PyResolve`, `Plugins.PyPath`,
`Plugins.PluginFileList.PluginFile_1.FILE_NAME`).

## Concurrency constraint

**MIKE SHE is not reentrant.** Two simulations cannot share a process, and a process cannot be reused
for a second simulation — which rules out `multiprocessing.Pool`. The working pattern is a
`ThreadPoolExecutor` bounding concurrency, with each thread spawning a fresh `Process` per simulation
(`run_variants_parallel_pool` in `run_all_examples.py`). Parallel variants also need distinct `.she`
copies; the same setup file cannot be run twice concurrently.

## Style

- New/maintained code uses **2-space indentation**; the older `DataProcessing/PtPathline*.py` and
  `ReadPtBin.py` use 4. Match the file you are in.
- Every script starts with a header comment block: `Subject`, `Usage`, `Dependencies`/`Requires`,
  optionally `Guarantee`/`Limitations`, then author and date. Keep it and keep it accurate — for
  plugins it is the only user documentation, and several files document known limitations there
  (e.g. `openBoreHoles.py` on explicit coupling lag, static exchange area and time-varying K).
- Comments in this codebase carry hydrological and numerical reasoning (why a formulation is a
  simplification, what a diagnostic means physically). That density is intentional; preserve it.
- `snake_case` for new code; some older files and MShePy call sites are `camelCase`.
