import argparse
import os
import numpy as np
import matplotlib.pyplot as plt


def build_line_probes(x_0,L, level):
    dx = L / (2 ** level)
    
    if x_0 < 0: #upstream location
        x_str = -1.01052632
        x_end = 0.0
        z_str = 0.42498737
        z_end = 0.5 
        slope = (z_end - z_str) /(x_end - x_str)
        zlim = z_str + slope*(x_0 - x_str)
        print("zlim is %f" % zlim)
    else:
        zlim = 0.5

    y = np.arange(-(1.0 - dx / 2.0), (1.0 - dx / 2.0) + dx * 0.5, dx)
    z = np.arange(-(zlim - dx / 2.0), (zlim - dx / 2.0) + dx * 0.5, dx)
    
    # Line definitions (in order provided by user)
    # First two lines along y, second two lines along z
    lines = [
        ("y", x_0, 0.0),   # x=0, z=0
        ("y", x_0, 0.45),  # x=0, z=0.45
        ("z", x_0, 0.0),   # x=0, y=0
        ("z", x_0, 0.95),  # x=0, y=0.95
    ]

    xyz_all = []
    for axis, x_val, other_val in lines:
        if axis == "y":
            xyz = np.zeros((y.size, 3))
            xyz[:, 0] = x_val
            xyz[:, 1] = y
            xyz[:, 2] = other_val
        else:
            xyz = np.zeros((z.size, 3))
            xyz[:, 0] = x_val
            xyz[:, 1] = other_val
            xyz[:, 2] = z
        xyz_all.append(xyz)

    y_n = y.size
    z_n = z.size
    return xyz_all, np.vstack(xyz_all), dx, y_n, z_n


def plot_probes(xyz_lines, dx):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
    for i, xyz in enumerate(xyz_lines):
        ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color=colors[i % len(colors)], linewidth=1.5)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(f"BL probes (dx={dx})")
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    plt.show()


def main():
    L = 2.56
    level = 9
    outdir = "/Users/steven/OneDrive/Stanford/ONR project/Simulations Utilities" 
    fname_fmt = "BL_probes_x_{:.2f}.txt"
    x0=0.04
    fname = fname_fmt.format(x0)
    xyz_lines, xyz, dx, y_n, z_n = build_line_probes(x0,L, level)
    header = f"x y z  # BL probes, dx={dx}, Ny={y_n}, Nz={z_n}"
    outpath = os.path.join(outdir, fname)
    np.savetxt(outpath, xyz, delimiter=" ", header=header)
    plot_probes(xyz_lines, dx)


if __name__ == "__main__":
    main()
