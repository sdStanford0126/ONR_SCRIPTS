"""Define a structured LSA probe volume; preview only unless --write is set.

Run from any directory:
    python generate_lsa_volumetric_probes.py
    python generate_lsa_volumetric_probes.py --output probes.txt --write

The text file has one comment header and three columns: x y z. Rows advance
x fastest, then y, then z, all increasing. Assuming the sampling solver
preserves row order, restore sequential sampled values q with:
    Q = q.reshape((Nx, Ny, Nz), order="F")  # NumPy: Q[i, j, k]
    Q = reshape(q, [Nx, Ny, Nz]);          % MATLAB: Q(i, j, k)

Coordinates use h=1. Both transverse axes cover +/-8.8De, rounded outward
to multiples of 0.02h about zero. Only coordinate vectors are allocated in
preview mode. Writing streams bounded blocks rather than a full point cloud.
Requires NumPy.
"""

import argparse
from math import ceil
from pathlib import Path

import numpy as np


def coordinate_blocks(x, y, z, y_per_block=64):
    """Yield (n, 3) coordinate blocks in x-fastest order.

    Zero-based row i + Nx*j + Nx*Ny*k contains (x[i], y[j], z[k]).
    Blocks stop at each z plane and contain at most Nx*y_per_block rows.
    """
    if not isinstance(y_per_block, int) or y_per_block < 1:
        raise ValueError("y_per_block must be a positive integer.")
    nx = len(x)
    for zk in z:
        for start in range(0, len(y), y_per_block):
            ys = y[start:start + y_per_block]
            block = np.empty((nx * len(ys), 3), dtype=float)
            block[:, 0] = np.tile(x, len(ys))
            block[:, 1] = np.repeat(ys, nx)
            block[:, 2] = zk
            yield block


def generate_lsa_volumetric_probes(output_file=None, write_grid=False):
    """Return a definition dictionary and optionally write the probe file.

    Defaults to preview only. Existing output files are never overwritten.
    The returned x/y/z vectors define the structured array's three axes.
    """
    if not isinstance(write_grid, bool):
        raise TypeError("write_grid must be a boolean.")
    output_file = (Path(output_file) if output_file is not None else
                   Path(__file__).resolve().with_name("lsa_volumetric_probes.txt"))
    h = 1.0
    L = 2 * h
    dx, dy, dz = 0.05 * h, 0.02 * h, 0.02 * h
    De = 1.3 * (L * h) ** 0.625 / (L + h) ** 0.25
    requested_half_extent = 8.8 * De
    nx_intervals = round(35 * h / dx)
    ny_half = ceil(requested_half_extent / dy)
    nz_half = ceil(requested_half_extent / dz)
    x = np.arange(nx_intervals + 1) * dx
    y = np.arange(-ny_half, ny_half + 1) * dy
    z = np.arange(-nz_half, nz_half + 1) * dz
    Nx, Ny, Nz = len(x), len(y), len(z)
    point_count = Nx * Ny * Nz
    header = (
        f"# LSA volumetric probes; columns=x y z; h={h:.8f}; L={L:.8f}; "
        f"De={De:.12g}; dx={dx:.8f}; dy={dy:.8f}; dz={dz:.8f}; "
        f"x=[{x[0]:.8f},{x[-1]:.8f}]; y=[{y[0]:.8f},{y[-1]:.8f}]; "
        f"z=[{z[0]:.8f},{z[-1]:.8f}]; "
        "requested_yz_extent=+/-8.8De (rounded outward); "
        f"Nx={Nx}; Ny={Ny}; Nz={Nz}; N={point_count}; "
        "order=x-fastest,y-next,z-slowest; "
        "NumPy_reshape=(Nx,Ny,Nz),order=F; MATLAB_reshape=[Nx,Ny,Nz]"
    )
    # Exact byte count for ASCII %.8f columns, two spaces and an LF per row.
    estimated_bytes = (
        sum(len(f"{v:.8f}") for v in x) * Ny * Nz
        + sum(len(f"{v:.8f}") for v in y) * Nx * Nz
        + sum(len(f"{v:.8f}") for v in z) * Nx * Ny
        + 3 * point_count + len(header.encode("ascii")) + 1
    )
    grid = dict(
        h=h, L=L, De=De, dx=dx, dy=dy, dz=dz, x=x, y=y, z=z,
        Nx=Nx, Ny=Ny, Nz=Nz, point_count=point_count,
        requested_half_extent=requested_half_extent,
        estimated_bytes=estimated_bytes, output_file=output_file,
        header=header, order="x-fastest,y-next,z-slowest", written=False,
    )
    print(f"De/h = {De/h:.12g}; requested transverse half-extent/h = {requested_half_extent/h:.12g}")
    print(f"x/h: [{x[0]/h:g}, {x[-1]/h:g}]; "
          f"y/h: [{y[0]/h:g}, {y[-1]/h:g}]; z/h: [{z[0]/h:g}, {z[-1]/h:g}]")
    print(f"[Nx, Ny, Nz] = [{Nx}, {Ny}, {Nz}]; total probes = {point_count:,}")
    print(f"Text size: {estimated_bytes/1e9:.3f} GB ({estimated_bytes/2**30:.3f} GiB)")
    if not write_grid:
        print("Preview only: no grid file written.")
        return grid

    # Exclusive creation prevents accidental replacement. Context management
    # closes the file even if writing is interrupted or the disk fills.
    with output_file.open("x", encoding="ascii", newline="\n") as stream:
        stream.write(header + "\n")
        blocks_per_plane = ceil(Ny / 64)
        progress_interval = max(1, ceil(Nz / 20))
        for number, block in enumerate(coordinate_blocks(x, y, z), start=1):
            np.savetxt(stream, block, fmt="%.8f", delimiter=" ", newline="\n")
            if number % blocks_per_plane == 0:
                plane = number // blocks_per_plane
                if plane % progress_interval == 0 or plane == Nz:
                    print(f"Written z planes: {plane}/{Nz} ({100*plane/Nz:.1f}%)", flush=True)
    grid["written"] = True
    print(f"Saved {output_file}")
    return grid


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, help="Output path; defaults to lsa_volumetric_probes.txt beside this script.")
    parser.add_argument("--write", action="store_true", help="Write the grid file. Without this flag, only preview the definition.")
    args = parser.parse_args()
    try:
        generate_lsa_volumetric_probes(args.output, write_grid=args.write)
    except (OSError, ValueError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
