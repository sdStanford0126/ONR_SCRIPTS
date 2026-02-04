import numpy as np
import matplotlib
matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as plt
import pdb; pdb.set_trace()
plt.rcParams.update({
    'axes.linewidth':     1.5,
    'xtick.major.width' : 1.5,
    'ytick.major.width' : 1.5,
    "font.family" : "serif",
    "text.latex.preamble" : r"\usepackage{amsmath} \usepackage{times} \usepackage{fourier} \usepackage{bm}",
    "font.size" : 14,
    "legend.fontsize" : 12})
import sys
import os
import scipy.io as io

"""
The purpose of this script is to read injecito ports probe data
"""


