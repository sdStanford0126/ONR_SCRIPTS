#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 23:30:18 2025

@author: steven
"""

def setPlotpref(fontsize=14):
    import os
    tex_lib = "/Library/TeX/texbin"
    os.environ["PATH"] += ":" + tex_lib
    import matplotlib.pyplot as plt
    plt.rcParams.update({
     'axes.linewidth':     1.5,
     'xtick.major.width' : 1.5,
     'ytick.major.width' : 1.5,
     "font.family" : "serif",
     "text.latex.preamble" : r"\usepackage{amsmath} \usepackage{times} \usepackage{fourier} \usepackage{bm}",
     "font.size" : 14,
     "legend.fontsize" : 12})