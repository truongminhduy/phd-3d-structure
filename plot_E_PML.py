##################################################################
#### PLOT COMPARE ################################################
##################################################################
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
import matplotlib.gridspec as gridspec

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
res = np.load('PML/output.npz')['saveddic'].item()
# res = np.load('output.npz')['saveddic'].item()

# cel = res['cel']
# num = res['num']
# den = res['den']

# pole_root = res['pole_root']
# epsi_root = res['epsi_root']

eigs = res['eigs']

# plt.figure()
# ax1 = plt.subplot(1, 1, 1)
# ax1.plot(eigs.real, eigs.imag, 'C0o',markersize=3 ,mfc='none',label='Eigenmode')
# # for k in range(len(eigs.real)):
# #       plt.annotate(k, xy=(eigs.real[k],eigs.imag[k]), xytext=None, xycoords='data', textcoords='data', arrowprops=None)
# plt.ylabel(r"$\omega.imag$")
# # plt.xlabel(r"$\omega.real$")
# # plt.plot(pole_root.imag,-pole_root.real,'C2+',label='Pole',  mew=3, ms=7)
# # plt.plot(epsi_root.imag,-epsi_root.real,'C3+',label=r"$\omega_2$",  mew=3, ms=7)
# plt.legend()
# # plt.xlim(0, 130) 
# # plt.xlim(-70, 130)  
# ##### PLOT #################################################################################################################
# E2   = res['E2']
E2   = res['E2_allt']
# E0_0 = res['E0_0']
# E0_2 = res['E0_2']
E0_2 = res['E0_2t']
omegas = res['omegas']
tab_colors = ['C0','C1','C2','C3','C4','C5','C6','C7','C8']
# plt.figure(figsize=(12,15))
plt.figure(figsize=(10,10))
# plt.figure()

gs = gridspec.GridSpec(4, 1, height_ratios=[4,0.2,6,6]) 
###------------------------------------------------------------------------------##
ax1 = plt.subplot(gs[2])
ax1.plot(omegas,E2[:,0],'C18',label='Direct' , fillstyle='none',markersize=5)
# ax1.plot(omegas,E2[:,3],'C3-.',label=r'$f_{\rho}=\lambda^2$',linewidth=1.8)
ax1.plot(omegas,E2[:,1],'C0-',label=r'$f_{\rho}=1$', linewidth=1.6)
ax1.plot(omegas,E2[:,2],'C2-.',label=r'$f_{\rho}=\lambda$',linewidth=1.8)
# ax1.plot(omegas,E2[:,3],'C3-.',label=r'$f_{\rho}=\lambda^2$',linewidth=1.8)
plt.ticklabel_format(axis='y', scilimits=(-4,-3))
plt.ylabel(r"$\vert E \vert$")

plt.yscale('log')
# ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
#           ncol=3, fancybox=True, shadow=True)
plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=4)
# plt.legend()
# plt.xlim(24, 52)  
plt.xlim(20, 65) 
ax1.xaxis.grid()
# plt.axvline(x=33,color='C4')
# ##-----------------------------------######
ax2 = plt.subplot(gs[3])
ax2.plot(omegas,E0_2.real[:,0],'C18',label='Direct' , fillstyle='none',markersize=5)
# ax2.plot(omegas,E0_2.real[:,3],'C3-.',label=r'$f_{\rho}=\lambda^2$',linewidth=1.8)
ax2.plot(omegas,E0_2.real[:,1],'C0-',label=r'$f_{\rho}=1$', linewidth=1.6)
ax2.plot(omegas,E0_2.real[:,2],'C2-.',label=r'$f_{\rho)}=\lambda$',linewidth=1.8)
# ax2.plot(omegas,E0_2.real[:,3],'C3-.',label=r'$f_{\rho}=\lambda^2$',linewidth=1.8)
plt.ylabel(r'$\vert\Re(E_P)\vert$')
plt.xlabel(r"$\omega$")
plt.yscale('log')
# plt.ticklabel_format(axis='y', scilimits=(-4,-3))
plt.xlim(20, 65) 
ax2.xaxis.grid() 
# plt.legend()
# plt.axvline(x=33,color='C4')

ax3 = plt.subplot(gs[0])
ax3.plot(eigs.real, eigs.imag, 'C0o',markersize=4 ,mfc='none',label='Eigenmode')
plt.ylabel(r"$\Im(\omega)$")
plt.xlabel(r"$\Re(\omega)$")
# plt.plot(pole_root.imag,-pole_root.real,'C2+',label='Pole',  mew=3, ms=7)
# plt.plot(epsi_root.imag,-epsi_root.real,'C3+',label=r"$\omega_2$",  mew=3, ms=7)
plt.legend()
# plt.xlim(24, 52)  
plt.xlim(20, 65) 
ax3.xaxis.grid()
ax3.xaxis.set_ticks_position('top')
ax3.xaxis.set_label_position('top') 
# plt.axvline(x=33,color='C4')
# ##------------------------------------------------------------------------------##

plt.show()
