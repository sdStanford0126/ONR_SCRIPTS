#THis is being used to post processing time averaged
#results for the BL_interior_probes (axial profiles)
import os 
import numpy as np
import matplotlib as plt

L = 2.56 #baselength
LEVEL=9 #mesh refinement level
delta = L/2**LEVEL

#Nz is in fact a function of x
x_str = -1.0152632
x_end = 0.0
z_str = 0.42498737
z_end = 0.5
slope = (z_end-z_str)/(x_end-x_str)

#old versoin
xs = np.linspace(-1,0,5)
#new version

y = np.arange(-(1.0 - delta / 2.0), (1.0 - delta / 2.0) + delta * 0.5, delta)




def z_lim(x):
    if x > 0:
        raise ValueError("x must be smaller or equal to 0")
    else:
        z_lim = z_str + slope*(x - x_str)
    return z_lim

def findNz(x,delta):
    z_max = z_lim(x)
    z = np.arange(-(z_max - delta / 2.0), (z_max - delta / 2.0) + delta * 0.5, delta)
    Nz = len(z)
    return Nz

Nx = xs.size
Ny = y.size
Nzs = [findNz(x,delta) for x in xs]

print("Nz_position by count",Nzs)
print("total probe pts:", np.sum(Ny*Nzs))

#TODO: need to implement extraction of each axial plane of probes, 
#first need to get index from .pxyz and record x,y,z locatiions (no read original probe file)
#then for a given x value, find the closest plane and extract the data from 
# that plane and organize them in the y,z plane


def extract_plane():
    pass


