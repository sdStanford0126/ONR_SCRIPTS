# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Utility scripts for the ONR jet acoustics project — a computational fluid dynamics (CFD) study of supersonic rectangular jets (with and without port injection). The simulations run on the Anvil HPC cluster at Purdue. Post-processing scripts read probe data output by the solver and produce turbulence statistics and spectral plots.

## Running Scripts

Scripts are run directly as Python modules from within their directory (they use relative imports and hardcoded I/O paths):

```bash
# Boundary layer axial profile post-processing (MPI-parallel, runs on Anvil)
cd python_postproc
python BL_ax_prof_post.py

# Probe generation (generates .txt probe coordinate files for the solver)
cd python_probes
python write_BL_probes.py
python write_BL_interior_probes.py

# Older BL probes post-processing
cd python_postproc
python BL_probes_post.py

# Axial profile post-processing (reads/writes HDF5)
cd python_postproc
python ax_prof_post.py
```

MATLAB scripts in `matlab_scripts/` are run directly in MATLAB; see the [critical matlab scripts doc](https://docs.google.com/document/d/19NgOg48rY7gr2NTyHzAptX-TrWndukzRqh06U9prddA/edit?usp=sharing) for original file locations.

## Key Dependencies

- `numpy`, `matplotlib`, `scipy`, `h5py` — core scientific stack
- `mpi4py` — MPI parallelism (used in `BL_ax_prof_post.py`)
- `tqdm` — progress bars
- `concurrent.futures.ThreadPoolExecutor` — threaded I/O parallelism in `BL_ax_prof_post.py`
- LaTeX + TeX fonts (via `setPlotpref`) — required for publication-quality plots; the path `/Library/TeX/texbin` is macOS-specific and must exist or be updated for the current environment

## Code Architecture

### Data Flow

1. **Probe generation** (`python_probes/`) — writes plain-text `.txt` files of (x,y,z) coordinates that are fed to the CFD solver as probe locations.

2. **Solver output** — produces two file types per timestep:
   - `.pxyz` — position file: columns are `[probe_index, x, y, z]`; indices are **not** contiguous (sparse/reordered), so all readers reconstruct a dense array via `arr[ind] = data`.
   - `.pcd` — data file per timestep: columns are flow variables (u, v, w, p, rho or T depending on probe type).

3. **Post-processing** (`python_postproc/`) — reads `.pxyz` + time series of `.pcd` files, computes statistics, and saves `.png` plots (and intermediate HDF5 in `ax_prof_post.py`).

### Important Conventions

**Index reordering pattern** — used in every reader: the `.pxyz` file stores a sparse probe index in column 0. Data must be placed at `arr[ind]` before slicing by spatial coordinate:
```python
Pos = np.loadtxt(posName, skiprows=1)
ind = Pos[:, 0].astype(int)
X = np.zeros(Npts); X[ind] = Pos[:, 1]
```

**Geometry constants** (shared across probe writers and post-processors):
- `L = 2.56` — nozzle base length (non-dimensional)
- `LEVEL = 9` — mesh refinement level; `delta = L / 2**LEVEL` is the grid spacing
- Nozzle converging section: x ∈ [-1.0152632, 0], z_wall varies linearly from 0.42498737 to 0.5
- Major axis (y): wall at ±1.0; Minor axis (z): wall at ±0.5 (at nozzle exit)

**Non-dimensionalization**: All quantities are non-dimensional. Viscosity uses a power-law: `mu = mu_ref * (T/T_ref)^n` with `mu_ref=3.26788998e-6`, `T_ref=1.0`, `n=0.7`, `gamma=1.4`. Pressure is used as a proxy for temperature via ideal gas: `T = p/rho`.

**Wall shear / friction velocity**: computed via second-order one-sided finite differences at the wall (`-3f[0] + 4f[1] - f[2]) / (2*dy`), then `u_tau = sqrt(tau_w / rho_w)`. y+ is `y * u_tau / nu`.

**Streamwise averaging** in `plotTurbProf_C`: profiles are averaged over a band near the centerline (±y_lim=0.85, ±z_lim=0.35) rather than a single centerline point.

### Module Structure

| Path | Purpose |
|------|---------|
| `python_postproc/Universal_Subroutines.py` | Shared `setPlotpref()` — sets matplotlib rcParams for publication plots |
| `python_postproc/BL_ax_prof_post.py` | Main BL post-proc: threaded `.pcd` loading, turbulence stats, law-of-wall plots |
| `python_postproc/BL_ax_prof_avg_post.py` | Earlier/simpler version of BL axial profile post-proc |
| `python_postproc/BL_probes_post.py` | Line-probe post-proc at nozzle exit (4 lines: 2 y-lines, 2 z-lines) |
| `python_postproc/ax_prof_post.py` | Far-field axial profile post-proc; converts `.pcd` time series → HDF5 |
| `python_probes/write_BL_probes.py` | Generates 4-line BL probe coordinates at a given x station |
| `python_probes/write_BL_interior_probes.py` | Generates volumetric BL probes over the converging nozzle section |
| `matlab_scripts/post_processing/` | MATLAB PSD, SPOD, FWH, and surface-plot scripts |
| `matlab_scripts/probes/` | MATLAB probe-coordinate generators |
| `python_old/` | Superseded scripts; kept for reference |

### I/O Paths

All data paths are hardcoded to Anvil scratch: `/anvil/scratch/x-sdai/` or `/anvil/scratch/x-akhiln/`. Update the `data_dir`, `out_dir`, and related variables in `main()` when switching cases. Output directories must exist before running (scripts do not create them).
