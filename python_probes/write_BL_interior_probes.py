import os
import numpy as np
import matplotlib.pyplot as plt

"""
build probes for BL asessment on axial planes on the converging section
"""
L = 2.56 #baselength
LEVEL=9 #mesh refinement level
delta = L/2**LEVEL

x_str = -1.0152632
x_end = 0.0
z_str = 0.42498737
z_end = 0.5
slope = (z_end-z_str)/(x_end-x_str)

def z_lim(x):
    if x > 0:
        raise ValueError("x must be smaller or equal to 0")
    else:
        z_lim = z_str + slope*(x - x_str)
    return z_lim

xyz_all = []

xs = np.linspace(-1,0,11)
y = np.arange(-(1.0 - delta / 2.0), (1.0 - delta / 2.0) + delta * 0.5, delta)

for x in xs:
    z_max = z_lim(x)
    z = np.arange(-(z_max - delta / 2.0), (z_max - delta / 2.0) + delta * 0.5, delta)

    yz_count = y.size * z.size
    xyz = np.empty((yz_count, 3), dtype=float)
    xyz[:, 0] = x
    xyz[:, 1] = np.repeat(y, z.size)
    xyz[:, 2] = np.tile(z, y.size)

    xyz_all.append(xyz)

xyz_comb = np.vstack(xyz_all)
print("yz_count:", yz_count)

outdir = "/Users/steven/OneDrive/Stanford/ONR project/Simulations Utilities" 
fname = "BL_interior_axprof.txt"
header = f"x y z  # BL probes, delta={delta}, Nlayer={xs.size}, Ny={y.size}, Nz=varied"
outpath = os.path.join(outdir,fname)
np.savetxt(outpath,xyz_comb,delimiter=" ", header=header)
#debug validation 
#fig=plt.figure()
#ax = fig.add_subplot(projection='3d')
#ax.scatter(xyz_comb[:,0], xyz_comb[:,1], xyz_comb[:,2])
#ax.set_aspect('equal')
#ax.view_init(0,-90)
#plt.show()
