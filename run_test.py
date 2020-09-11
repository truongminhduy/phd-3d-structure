import numpy as np
import scipy as sc
import scipy.io as scio
import pylab as pl
import os
import sys
import scipy.signal as signal
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

np.set_printoptions(precision=3)
pi = np.pi
### MAIN 1 #######################################################################
path = ''
# os.system(path+'gmsh mesh_3D_line.geo')
os.system(path+'gmsh mesh_compute.geo')
# os.system(path+'gmsh mesh_3D_2sphere.geo')
# os.system(path+'gmsh mesh_3D_line_simple.geo')
# os.system(path+'gmsh mesh_plot.geo u.pos')
# os.system(path+'gmsh mesh_3D_plot.geo')