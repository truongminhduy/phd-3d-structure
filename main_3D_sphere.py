##############################################################################################################################
###### SPHERE IN BOX #########################################################################################################
##############################################################################################################################
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as linalg
import scipy.signal as signal
import matplotlib.gridspec as gridspec
from function_e import *
import os
import sys
np.set_printoptions(precision=4)
# os.environ['OPENBLAS_NUM_THREADS'] = '4'
import time
start = time.time()
###### MAIN PROGRAM ##########################################################################################################
## Mau's constants
pi    = np.pi
micro = 1e-6
femto = 1e-15
peta  = 1e15
cc    = 299792458
mu_0  = pi*4e-7
eps_0 = 1./(mu_0*cc**2)
##--- PARAMETER ----------------#########################
scale_light = 1e8
micro       = 1e-6
scale_frequ = scale_light/micro # 1e14
cel = cc/scale_light
mu0 = mu_0  *scale_light
ep0 = eps_0 *scale_light

e_scale = 2e-1
r_ellipse_2 = 1*e_scale

xS , yS , zS  = 1.2 *e_scale, 1.2 *e_scale, 1.2 *e_scale
xD1, yD1, zD1 =-1.2 *e_scale,-1.2 *e_scale,-1.2 *e_scale
xD2, yD2, zD2 = 1.2 *e_scale,-1.2 *e_scale, 1.2 *e_scale

r_source =  10e-2*e_scale

lambda_m = 0.25*e_scale
D_pml_in = 3*e_scale
pml_size = 4*e_scale

epsr = 16+1j
epsilon_oo = 3.36174
omega_p    = 1.3388e16/scale_frequ
gamma      = 7.07592e13/scale_frequ

eps_mil = 1
E_or_H  = 0			# electric or magnetic
PML_or_not = 0
# vector0vs1  = 0	# scalar or vector
type_source = 0 	# line source or plane wave
rati_or_poly = 1 	# rational or polynomial eigenvalue problem??
neig = 50

par_gmsh_getdp = open('parameters_gmsh_getdp.dat', 'w')
par_gmsh_getdp.write('neig = %d;\n'%(neig))
par_gmsh_getdp.write('xS          = %3.3e;\n'%(xS))
par_gmsh_getdp.write('yS          = %3.3e;\n'%(yS))
par_gmsh_getdp.write('zS          = %3.3e;\n'%(zS))
par_gmsh_getdp.write('xD1         = %3.3e;\n'%(xD1))
par_gmsh_getdp.write('yD1         = %3.3e;\n'%(yD1))
par_gmsh_getdp.write('zD1         = %3.3e;\n'%(zD1))
par_gmsh_getdp.write('xD2         = %3.3e;\n'%(xD2))
par_gmsh_getdp.write('yD2         = %3.3e;\n'%(yD2))
par_gmsh_getdp.write('zD2         = %3.3e;\n'%(zD2))
# par_gmsh_getdp.write('r_ellipse_1 = %3.3e;\n'%(r_ellipse_1))
par_gmsh_getdp.write('r_ellipse_2 = %3.3e;\n'%(r_ellipse_2))
# par_gmsh_getdp.write('r_ellipse_3 = %3.3e;\n'%(r_ellipse_3))
par_gmsh_getdp.write('r_source    = %3.3e;\n'%(r_source))
par_gmsh_getdp.write('lambda_m    = %3.3e;\n'%(lambda_m))
par_gmsh_getdp.write('D_pml_in    = %3.3e;\n'%(D_pml_in))
par_gmsh_getdp.write('pml_size    = %3.3e;\n'%(pml_size))
par_gmsh_getdp.write('cel         = %3.3e;\n'%(cel))
par_gmsh_getdp.write('mu0         = %3.3e;\n'%(mu0))
par_gmsh_getdp.write('epsilon0    = %3.3e;\n'%(ep0))
par_gmsh_getdp.write('eps_mil_re = %3.3e;\n'%(eps_mil.real))
par_gmsh_getdp.write('eps_mil_im = %3.3e;\n'%(eps_mil.imag))
par_gmsh_getdp.write('eps_diff_re = %3.3e;\n'%(epsr.real ))
par_gmsh_getdp.write('eps_diff_im = %3.3e;\n'%(epsr.imag ))
par_gmsh_getdp.write('PML         = 0.5;\n')
# par_gmsh_getdp.write('PML         = 0;\n')
# par_gmsh_getdp.write('angle       = Pi/10;\n')
par_gmsh_getdp.write('E_or_H      = %d;\n'%(E_or_H))
par_gmsh_getdp.write('PML_or_not  = %d;\n'%(PML_or_not))
par_gmsh_getdp.write('type_source = %d;\n'%(type_source))
par_gmsh_getdp.write('epsilon_oo  = %3.5e;\n'%(epsilon_oo))
par_gmsh_getdp.write('omega_p     = %3.5e;\n'%(omega_p))
par_gmsh_getdp.write('gamma       = %3.5e;\n'%(gamma))
par_gmsh_getdp.write('zone1       = 1;\n')
par_gmsh_getdp.write('zone2       = 30;\n')
par_gmsh_getdp.close()

gmsh_path  = ''
getdp_path = ''
if PML_or_not == 1:
	mesh_filename = 'mesh_3D_line.msh'
	mesh_geo      = 'mesh_3D_line.geo'
else:
	mesh_filename = 'mesh_3D_line_simple.msh'
	mesh_geo      = 'mesh_3D_line_simple.geo'
os.system(gmsh_path+'gmsh '+mesh_geo+' -3 -o '+mesh_filename)

###### SPECTRAL PROBLEM ############################################################################
slepc_options_pep = ' -pep_max_it 400 -pep_target_magnitude'
slepc_options_rational = ' -nep_max_it 200 -nep_type nleigs -nep_rational -nep_target_magnitude -petsc_prealloc 100 \
					-rg_interval_endpoints -50,50,-70,70'# -v 0'
					# ' -nep_max_it 100 -nep_type nleigs -nep_rational -nep_target_magnitude -rg_interval_endpoints -50,9,-30,50'
spec_getdp = getdp_path+'getdp spectral_3D.pro -pre Projection -msh '+mesh_filename+' -cal -pos postop -slepc'# -v 0' 
os.system(spec_getdp + slepc_options_rational)
eigs = np.loadtxt('EigenValuesReal.txt',usecols=[5]) + 1j*np.loadtxt('EigenValuesImag.txt',usecols=[5])
neig = len(eigs)					

plt.figure()
vpr = np.loadtxt('EigenValuesReal.txt',usecols=[5])
vpi = np.loadtxt('EigenValuesImag.txt',usecols=[5])
plt.plot(vpr, vpi, 'C0o',markersize=3 ,mfc='none')
for k in range(len(vpr)):
      plt.annotate(k, xy=(vpr[k],vpi[k]), xytext=None, xycoords='data', textcoords='data', arrowprops=None)
plt.title("Eigenvalues")
plt.xlabel("Re")
plt.ylabel("Im")

if rati_or_poly == 0:
	os.system('rm norm1.txt norm2.txt')
	norm_getdp = getdp_path+'getdp norm.pro -pre Projection -msh '+mesh_filename+' -res spectral_3D.res -pos postop_norm -v 0'
	os.system(norm_getdp)
	norm1 = np.loadtxt('norm1.txt',usecols=[5]) + 1j*np.loadtxt('norm1.txt',usecols=[6] )  
	norm2 = np.loadtxt('norm2.txt',usecols=[5]) + 1j*np.loadtxt('norm2.txt',usecols=[6] ) 
	d1 = 1/cel**2 * ( 2*epsilon_oo*eigs - 1j*gamma*omega_p**2/(eigs+1j*gamma)**2 ) * norm1
	d2 = 2/cel**2 * eigs * norm2 
	dd = d1+d2
else:
	os.system('rm norm1.txt norm2.txt')
	norm_getdp = getdp_path+'getdp norm.pro -pre Projection -msh '+mesh_filename+' -res spectral_3D.res -pos postop_norm -v 0'
	os.system(norm_getdp)
	norm1 = np.loadtxt('norm1.txt',usecols=[5]) + 1j*np.loadtxt('norm1.txt',usecols=[6] )  
	norm2 = np.loadtxt('norm2.txt',usecols=[5]) + 1j*np.loadtxt('norm2.txt',usecols=[6] )
	norm3 = np.loadtxt('norm3.txt',usecols=[5]) + 1j*np.loadtxt('norm3.txt',usecols=[6] )  
	d1 = 1/cel**2 * (3*epsilon_oo*eigs**2 + 2*epsilon_oo*gamma*eigs*1j - omega_p**2) * norm1
	d2 = 1/cel**2 * (3*eigs**2 + 2*gamma*1j*eigs)* norm2 
	d3 = -1* norm3 
	dd1 = d1+d2+d3
##################################################################################################################################
N_frequency = 10
E0_0 = np.zeros((N_frequency,4),dtype=complex)
E0_1 = np.zeros((N_frequency,4),dtype=complex)
E0_2 = np.zeros((N_frequency,4),dtype=complex)
E2   = np.zeros((N_frequency,4))
E2_all = np.zeros((N_frequency,4))
omegas = np.linspace(25,50,N_frequency)
for i,omega0 in enumerate(omegas) :
	lambda0 = 2.*pi*cel/omega0	
	# omega0 = 2.*pi*cel/lambda0
	epsr = epsilon_oo - omega_p**2/(omega0**2+1j*omega0*gamma) 
	# print(epsr)
	par_gmsh_getdp = open('parameters_gmsh_getdp.dat', 'a')
	par_gmsh_getdp.write('lambda0 = %3.3e;\n'%(lambda0))
	par_gmsh_getdp.write('neig = %d;\n'%(neig))
	par_gmsh_getdp.write('eps_diff_re = %3.3e;\n'%(epsr.real ))
	par_gmsh_getdp.write('eps_diff_im = %3.3e;\n'%(epsr.imag ))
	par_gmsh_getdp.close()
	############ direct ##################################################################
	direct_getdp = getdp_path+'getdp direct_3D.pro -pre Scattering -msh '+mesh_filename+' -cal -pos postop_scat -v 0' 
	os.system(direct_getdp)
	E2[i,0] = np.loadtxt('E2.txt')[1]
	E2_all[i,0] = np.loadtxt('E2_all.txt')[1]
	E0_txt, E1_txt, E2_txt = np.loadtxt('E0_0.txt'), np.loadtxt('E0_1.txt'), np.loadtxt('E0_2.txt')
	E0_0[i,0] = np.sqrt( E0_txt[5]**2 + E0_txt[6]**2 + E0_txt[7]**2 ) + 1j*np.sqrt( E0_txt[8]**2 + E0_txt[9]**2 + E0_txt[10]**2 )
	E0_1[i,0] = np.sqrt( E1_txt[5]**2 + E1_txt[6]**2 + E1_txt[7]**2 ) + 1j*np.sqrt( E1_txt[8]**2 + E1_txt[9]**2 + E1_txt[10]**2 )
	E0_2[i,0] = np.sqrt( E2_txt[5]**2 + E2_txt[6]**2 + E2_txt[7]**2 ) + 1j*np.sqrt( E2_txt[8]**2 + E2_txt[9]**2 + E2_txt[10]**2 )
	#### COMPUTE PROJECTION PROBLEM ######################################################### 
	#### norm ------------------------------------#####
	os.system('rm Jns.txt')
	Jn_getdp = getdp_path+'getdp Jn.pro -pre Projection -msh '+mesh_filename+' -res spectral_3D.res -pos postop_Jn -v 0'
	os.system(Jn_getdp)
	Jn  = np.loadtxt('Jns.txt',usecols=[5]) + 1j*np.loadtxt('Jns.txt',usecols=[6] ) 
	### Pns  #########################
	if rati_or_poly == 0:
		Pns0 = Jn/(omega0-eigs)/dd
		Pns1 = Jn/(omega0-eigs)/dd*(eigs/omega0)
		Pns2 = Jn/(omega0-eigs)/dd*((eigs-33)/(omega0-33))
	else:
		Pns0 = Jn*(eigs+1j*gamma)/(omega0-eigs)/dd1*(eigs/omega0)**0
		Pns1 = Jn*(eigs+1j*gamma)/(omega0-eigs)/dd1*(eigs/omega0)**1
		Pns2 = Jn*(eigs+1j*gamma)/(omega0-eigs)/dd1*(eigs/omega0)**2
	#####################################################################	
	file_Pns = open('Pns.dat', 'w')
	for k, Pn in enumerate(Pns0):
	    file_Pns.write('Pns_re_%d = %3.5e;\n'%(k,Pn.real))
	    file_Pns.write('Pns_im_%d = %3.5e;\n'%(k,Pn.imag))
	file_Pns.close()
	project_getdp = getdp_path+'getdp projection_3D.pro -pre Projection -msh  '+mesh_filename+' -res spectral_3D.res -pos  postop_modal -bin -v 0'
	os.system(project_getdp)
	os.system('rm E2.txt E2_all.txt E0_0.txt E0_1.txt E0_2.txt')
	projection_coeff_getdp = getdp_path+'getdp projection_coeff_3D.pro -pre Projection -msh  '+mesh_filename+' -cal -pos postop_modal -v 0' 
	os.system(projection_coeff_getdp)
	E2[i,1] = np.loadtxt('E2.txt')[1]
	E2_all[i,1] = np.loadtxt('E2_all.txt')[1]
	E0_txt, E1_txt, E2_txt = np.loadtxt('E0_0.txt'), np.loadtxt('E0_1.txt'), np.loadtxt('E0_2.txt')
	E0_0[i,1] = np.sqrt( E0_txt[5]**2 + E0_txt[6]**2 + E0_txt[7]**2 ) + 1j*np.sqrt( E0_txt[8]**2 + E0_txt[9]**2 + E0_txt[10]**2 )
	E0_1[i,1] = np.sqrt( E1_txt[5]**2 + E1_txt[6]**2 + E1_txt[7]**2 ) + 1j*np.sqrt( E1_txt[8]**2 + E1_txt[9]**2 + E1_txt[10]**2 )
	E0_2[i,1] = np.sqrt( E2_txt[5]**2 + E2_txt[6]**2 + E2_txt[7]**2 ) + 1j*np.sqrt( E2_txt[8]**2 + E2_txt[9]**2 + E2_txt[10]**2 )
	########2----------------------------------#####
	file_Pns = open('Pns.dat', 'w')
	for k, Pn in enumerate(Pns1):
	    file_Pns.write('Pns_re_%d = %3.5e;\n'%(k,Pn.real))
	    file_Pns.write('Pns_im_%d = %3.5e;\n'%(k,Pn.imag))
	file_Pns.close()
	project_getdp = getdp_path+'getdp projection_3D.pro -pre Projection -msh  '+mesh_filename+' -res spectral_3D.res -pos  postop_modal -bin -v 0'
	os.system(project_getdp)
	os.system('rm E2.txt E2_all.txt E0_0.txt E0_1.txt E0_2.txt')
	projection_coeff_getdp = getdp_path+'getdp projection_coeff_3D.pro -pre Projection -msh  '+mesh_filename+' -cal -pos postop_modal -v 0' 
	os.system(projection_coeff_getdp)
	E2[i,2] = np.loadtxt('E2.txt')[1]
	E2_all[i,2] = np.loadtxt('E2_all.txt')[1]
	E0_txt, E1_txt, E2_txt = np.loadtxt('E0_0.txt'), np.loadtxt('E0_1.txt'), np.loadtxt('E0_2.txt')
	E0_0[i,2] = np.sqrt( E0_txt[5]**2 + E0_txt[6]**2 + E0_txt[7]**2 ) + 1j*np.sqrt( E0_txt[8]**2 + E0_txt[9]**2 + E0_txt[10]**2 )
	E0_1[i,2] = np.sqrt( E1_txt[5]**2 + E1_txt[6]**2 + E1_txt[7]**2 ) + 1j*np.sqrt( E1_txt[8]**2 + E1_txt[9]**2 + E1_txt[10]**2 )
	E0_2[i,2] = np.sqrt( E2_txt[5]**2 + E2_txt[6]**2 + E2_txt[7]**2 ) + 1j*np.sqrt( E2_txt[8]**2 + E2_txt[9]**2 + E2_txt[10]**2 )
	########3----------------------------------#####
	file_Pns = open('Pns.dat', 'w')
	for k, Pn in enumerate(Pns2):
	    file_Pns.write('Pns_re_%d = %3.5e;\n'%(k,Pn.real))
	    file_Pns.write('Pns_im_%d = %3.5e;\n'%(k,Pn.imag))
	file_Pns.close()
	project_getdp = getdp_path+'getdp projection_3D.pro -pre Projection -msh  '+mesh_filename+' -res spectral_3D.res -pos  postop_modal -bin -v 0'
	os.system(project_getdp)
	os.system('rm E2.txt E2_all.txt E0_0.txt E0_1.txt E0_2.txt')
	projection_coeff_getdp = getdp_path+'getdp projection_coeff_3D.pro -pre Projection -msh  '+mesh_filename+' -cal -pos postop_modal -v 0' 
	os.system(projection_coeff_getdp)
	E2[i,3] = np.loadtxt('E2.txt')[1]
	E2_all[i,3] = np.loadtxt('E2_all.txt')[1]
	E0_txt, E1_txt, E2_txt = np.loadtxt('E0_0.txt'), np.loadtxt('E0_1.txt'), np.loadtxt('E0_2.txt')
	E0_0[i,3] = np.sqrt( E0_txt[5]**2 + E0_txt[6]**2 + E0_txt[7]**2 ) + 1j*np.sqrt( E0_txt[8]**2 + E0_txt[9]**2 + E0_txt[10]**2 )
	E0_1[i,3] = np.sqrt( E1_txt[5]**2 + E1_txt[6]**2 + E1_txt[7]**2 ) + 1j*np.sqrt( E1_txt[8]**2 + E1_txt[9]**2 + E1_txt[10]**2 )
	E0_2[i,3] = np.sqrt( E2_txt[5]**2 + E2_txt[6]**2 + E2_txt[7]**2 ) + 1j*np.sqrt( E2_txt[8]**2 + E2_txt[9]**2 + E2_txt[10]**2 )

# os.system(gmsh_path+'gmsh '+mesh_geo+' us.pos up.pos')
# os.system(gmsh_path+'gmsh '+mesh_geo+' us.pos up.pos -merge compare.geo')
dumpglob2npz('output.npz',globals()) 
plt.show()