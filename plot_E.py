################################################################################
#### PLOT COMPARE ################################################
####################################################################!!!#########
import numpy as np
import scipy as sc
import scipy.io as scio
import pylab as plt
import pandas as pd
import sys
import seaborn as sns
from function_e import *
# print(sys.version) 
np.set_printoptions(precision=4)
import os
import sys
import matplotlib.pyplot as plt
plt.rc('font', size=13)
pi = np.pi

plt.rcParams['font.family'] = 'serif'
plt.style.use('seaborn-muted')
plt.rc('text', usetex=True)
plt.rc('font', **{'family' : "serif"})
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)
#plt.rcParams['text.latex.preamble'] = [r'\boldmath']
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 16
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['figure.titlesize'] = 18

# #### rational function (N/D)' ###############################################
# res = np.load('noPML/output.npz')['saveddic'].item()
res = np.load('output.npz')['saveddic'].item()

E2   = res['E2']
E2_all = res['E2_all']
E0_0 = res['E0_0']
E0_1 = res['E0_1']
E0_2 = res['E0_2']
cel = res['cel']
# num = res['num']
# den = res['den']
omegas    = res['omegas']
# pole_root = res['pole_root']
# epsi_root = res['epsi_root']

eigs = res['eigs']
# eigs0 = res['eigs0']
# eigs1 = res['eigs1']
# eigs2 = res['eigs2']
# eigs = np.concatenate((eigs0, eigs1))

# eigs = eigs0
plt.figure()
# ax1 = plt.subplot(2, 1, 1)
ax1 = plt.subplot(1, 1, 1)
ax1.plot(eigs.real, eigs.imag, 'C0o',markersize=3 ,mfc='none',label='Eigenmode')
plt.ylabel(r"$\omega.imag$")
# plt.xlabel(r"$\omega.real$")
# plt.plot(pole_root.imag,-pole_root.real,'C2+',label='Pole',  mew=3, ms=7)
# plt.plot(epsi_root.imag,-epsi_root.real,'C3+',label=r"$\omega_2$",  mew=3, ms=7)
plt.legend()
# plt.xlim(0, 130) 
# plt.xlim(-70, 130)  


# ##### PLOT #################################################################################################################
tab_colors = ['C0','C1','C2','C3','C4','C5','C6','C7','C8']
plt.figure(figsize=(10,9))
###------------------------------------------------------------------------------##
ax1 = plt.subplot(3, 1, 2)
ax1.plot(omegas,E2[:,0],'C18',label='Direct' , fillstyle='none',markersize=6)
ax1.plot(omegas,E2[:,3],'C3-.',label=r'$f_{\rho}=\lambda-\lambda_0$',linewidth=1.8)
ax1.plot(omegas,E2[:,1],'C0-',label=r'$f_{\rho}=1$', linewidth=1.8)
ax1.plot(omegas,E2[:,2],'C2:',label=r'$f_{\rho}=\lambda$',linewidth=1.8)
# ax1.plot(omegas,E2[:,2],'C2-.',label=r'$f_{\rho}=\lambda$',linewidth=1.8)

plt.ticklabel_format(axis='y', scilimits=(-4,-3))
plt.ylabel(r"$\vert E \vert$")
ax1.xaxis.grid() 
# plt.xlabel("omega")
plt.yscale('log')
plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=4)
# plt.legend()


# #-----------------------------------######
ax2 = plt.subplot(3, 1, 3)
ax2.plot(omegas,E0_2.real[:,0],'C18',label='Direct' , fillstyle='none',markersize=5)
ax2.plot(omegas,E0_2.real[:,1],'C0-',label=r'$f_{\rho}=1$', linewidth=1.6)
ax2.plot(omegas,E0_2.real[:,2],'C2:' ,label=r'$f_{\rho}=\lambda$',linewidth=1.8)
# ax2.plot(omegas,E0_2.real[:,3],'C3-.', label=r'$f_{\rho}=\lambda-\lambda_0$',linewidth=1.8)
plt.ylabel(r'$\Re(E_2)$')
plt.xlabel(r"$\omega$")
plt.yscale('log')
# plt.ticklabel_format(axis='y', scilimits=(-4,-3))
# plt.legend()
ax2.xaxis.grid() 

ax3 = plt.subplot(3, 1, 1)
ax3.plot(eigs.real, eigs.imag, 'C0o',markersize=4 ,mfc='none',label='Eigenmode')
plt.ylabel(r"$\omega.imag$")
# plt.xlabel(r"$\omega.real$")
# plt.plot(pole_root.imag,-pole_root.real,'C2+',label='Pole',  mew=3, ms=7)
# plt.plot(epsi_root.imag,-epsi_root.real,'C3+',label=r"$\omega_2$",  mew=3, ms=7)
plt.legend()
ax3.xaxis.grid()
ax3.xaxis.set_ticks_position('top')
ax3.xaxis.set_label_position('top') 
# ##------------------------------------------------------------------------------##

plt.show()
