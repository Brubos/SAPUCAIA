#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 16:08:48 2023

@author: humberto.junior and Bruno.souza
"""

#%% PACKAGES:

from optlnls.plot import plot_beam  #  functionalities for optical beam plotting
import matplotlib.pyplot as plt     #  Python plotting library
import optlnls.shadow               #  library related to optical simulations
import numpy as np                  #  library for numerical computation
import h5py                         #  library for handling HDF5 files

#%% FUNCTIONS: 
#%%     gaussian_beam function
def gaussian_beam(z, s0, z0, beta):
    return s0*np.sqrt(1 + ((z-z0)/beta)**2)

#%%     fit gaussian_beam function
def fit_gaussian_beam(z, s, auto_guess=True, p0=(1,1,1), n_points=0):
    
    # Packages:
    
    from scipy.optimize import curve_fit
    
    
    # Automatic Guess:
    
    if(auto_guess):
        
        arg_min = np.argmin(s)
        
        s0_guess = s[arg_min]
        
        z0_guess = z[arg_min]
        
        p0 = [s0_guess, z0_guess, 0.1]
    
    
    # Fit:
    
    popt, pcov = curve_fit(gaussian_beam, z, s, p0=p0, maxfev=10000) # z[np.argmin(s)]
    
    s0, z0, beta = popt
    
    
    # New arrays: 
    
    if(not(n_points)): n_points = len(z)    
    
    z_new = np.linspace(z[0], z[-1], n_points)
    
    s_new = gaussian_beam(z_new, s0, z0, beta)
    
    
    # Calculate R2:

    residuals = s - gaussian_beam(z, s0, z0, beta)

    ss_res = np.sum(residuals**2)
        
    ss_tot = np.sum((s-np.mean(s))**2)

    r_squared = 1 - (ss_res / ss_tot)
    
    
    # Prints:
    
    print('\n')    
    
    print('R² = %.6f' %(r_squared))
    
    print('s0 = %.3f' %s0)
    
    print('z0 = %.3f' %z0)
    
    print('Beta = %.3f' %beta)
    
    
    return z_new, s_new, popt

#%%     read_caustic function    
def read_caustic(filename, write_attributes=False, plot=False, plot2D=False, 
                 print_minimum=False, cmap='viridis', figprefix=''):
    
    with h5py.File(filename, 'r+') as f:
    
        # g = f['datasets']
    
        dset_names = list(f.keys())
        
        center_shadow = np.zeros((len(dset_names)-2, 2), dtype=float)
        center = np.zeros((len(dset_names)-2, 2), dtype=float)
        rms = np.zeros((len(dset_names)-2, 2), dtype=float)
        fwhm = np.zeros((len(dset_names)-2, 2), dtype=float)
        fwhm_shadow = np.zeros((len(dset_names)-2, 2), dtype=float)
        
        ###### READ DATA #######################
        zOffset = f.attrs['zOffset']        
        zStart = f.attrs['zStart']
        zFin = f.attrs['zFin']
        nz = f.attrs['nz']
                
        #if(plot2D): 
        xStart = f[dset_names[2]].attrs['xStart']
        xFin = f[dset_names[2]].attrs['xFin']
        nx = f[dset_names[2]].attrs['nx']
        yStart = f[dset_names[2]].attrs['yStart']
        yFin = f[dset_names[2]].attrs['yFin']
        ny = f[dset_names[2]].attrs['ny']
            
        histoH = np.zeros((nx, nz))
        histoV = np.zeros((ny, nz))
        
        z_points = np.linspace(zStart, zFin, nz) + zOffset
        
        for i, dset in enumerate(dset_names[2:]):
            #dset_keys = list(f[dset].attrs.keys())
        
            center_shadow[i,0] = f[dset].attrs['center_h_shadow']
            center_shadow[i,1] = f[dset].attrs['center_v_shadow']
            center[i,0] = f[dset].attrs['mean_h']
            center[i,1] = f[dset].attrs['mean_v']
            rms[i,0] = f[dset].attrs['rms_h']
            rms[i,1] = f[dset].attrs['rms_v']
            fwhm[i,0] = f[dset].attrs['fwhm_h'][0]
            fwhm[i,1] = f[dset].attrs['fwhm_v'][0]
            fwhm_shadow[i,0] = f[dset].attrs['fwhm_h_shadow']
            fwhm_shadow[i,1] = f[dset].attrs['fwhm_v_shadow']
            
            #if(plot2D):
            histo2D = np.array(f[dset])
            histoH[:,i] = histo2D.sum(axis=1)
            histoV[:,i] = histo2D.sum(axis=0)
                
    #### FIND MINIMUMS AND ITS Z POSITIONS

    rms_min = [np.min(rms[:,0]), np.min(rms[:,1])]
    fwhm_min = [np.min(fwhm[:,0]), np.min(fwhm[:,1])]
    fwhm_shadow_min = [np.min(fwhm_shadow[:,0]), np.min(fwhm_shadow[:,1])]

    rms_min_z=np.array([z_points[np.abs(rms[:,0]-rms_min[0]).argmin()],
                        z_points[np.abs(rms[:,1]-rms_min[1]).argmin()]])

    fwhm_min_z=np.array([z_points[np.abs(fwhm[:,0]-fwhm_min[0]).argmin()],
                         z_points[np.abs(fwhm[:,1]-fwhm_min[1]).argmin()]])

    fwhm_shadow_min_z=np.array([z_points[np.abs(fwhm_shadow[:,0]-fwhm_shadow_min[0]).argmin()],
                                z_points[np.abs(fwhm_shadow[:,1]-fwhm_shadow_min[1]).argmin()]])

    center_rms = np.array([center[:,0][np.abs(z_points-rms_min_z[0]).argmin()],
                           center[:,1][np.abs(z_points-rms_min_z[1]).argmin()]])

    center_fwhm = np.array([center[:,0][np.abs(z_points-fwhm_min_z[0]).argmin()],
                            center[:,1][np.abs(z_points-fwhm_min_z[1]).argmin()]])

    center_fwhm_shadow = np.array([center[:,0][np.abs(z_points-fwhm_shadow_min_z[0]).argmin()],
                                   center[:,1][np.abs(z_points-fwhm_shadow_min_z[1]).argmin()]])
 
    
    outdict = {'xStart': xStart,
               'xFin': xFin,
               'nx': nx,
               'yStart': yStart,
               'yFin': yFin,
               'ny': ny,
               'zStart': zStart,
               'zFin': zFin,
               'nz': nz,
               'zOffset': zOffset,
               'center_h_array': center[:,0], 
               'center_v_array': center[:,1],
               'center_shadow_h_array': center_shadow[:,0], 
               'center_shadow_v_array': center_shadow[:,1],
               'rms_h_array': rms[:,0], 
               'rms_v_array': rms[:,1],
               'fwhm_h_array': fwhm[:,0], 
               'fwhm_v_array': fwhm[:,1],
               'fwhm_shadow_h_array': fwhm_shadow[:,0], 
               'fwhm_shadow_v_array': fwhm_shadow[:,1],
               'rms_min_h': rms_min[0],
               'rms_min_v': rms_min[1],
               'fwhm_min_h': fwhm_min[0],
               'fwhm_min_v': fwhm_min[1],
               'fwhm_shadow_min_h': fwhm_shadow_min[0],
               'fwhm_shadow_min_v': fwhm_shadow_min[1],
               'z_rms_min_h': rms_min_z[0],
               'z_rms_min_v': rms_min_z[1],
               'z_fwhm_min_h': fwhm_min_z[0],
               'z_fwhm_min_v': fwhm_min_z[1],
               'z_fwhm_shadow_min_h': fwhm_shadow_min_z[0],
               'z_fwhm_shadow_min_v': fwhm_shadow_min_z[1],
               'center_rms_h': center_rms[0],
               'center_rms_v': center_rms[1],
               'center_fwhm_h': center_fwhm[0],
               'center_fwhm_v': center_fwhm[1],
               'center_fwhm_shadow_h': center_fwhm_shadow[0],
               'center_fwhm_shadow_v': center_fwhm_shadow[1],
               }
    
    if(write_attributes):
        with h5py.File(filename, 'a') as f:
            for key in list(outdict.keys()):
                f.attrs[key] = outdict[key]
                
            f.create_dataset('histoXZ', data=histoH, dtype=np.float, compression="gzip")
            f.create_dataset('histoYZ', data=histoV, dtype=np.float, compression="gzip")
            
    if(print_minimum):
        print('\n   ****** \n' + '   Z min (rms-hor): {0:.3e}'.format(rms_min_z[0]))
        print('   Z min (rms-vert): {0:.3e}\n   ******'.format(rms_min_z[1]))
        
    if(plot):
        
        plt.figure()
        plt.title('rms')
        plt.plot(z_points, rms[:,0], label='rms_h')
        plt.plot(z_points, rms[:,1], label='rms_v')
        plt.legend()
        plt.minorticks_on()
        plt.grid(which='both', alpha=0.2)    
    
        plt.figure()
        plt.title('fwhm')
        plt.plot(z_points, fwhm[:,0], label='fwhm_h')
        plt.plot(z_points, fwhm[:,1], label='fwhm_v')
        # plt.plot(z_points, fwhm_shadow[:,0], label='fwhm_h_shadow')
        # plt.plot(z_points, fwhm_shadow[:,1], label='fwhm_v_shadow')
        plt.legend()
        plt.minorticks_on()
        plt.grid(which='both', alpha=0.2)    
        
        plt.figure()
        plt.title('center')
        plt.plot(z_points, center[:,0], label='center_h')
        # plt.plot(z_points, center_shadow[:,0], label='center_h_shadow')
        plt.legend()
        plt.minorticks_on()
        plt.grid(which='both', alpha=0.2)    
        
        plt.figure()
        plt.title('center')
        plt.plot(z_points, center[:,1], label='center_v')
        # plt.plot(z_points, center_shadow[:,1], label='center_v_shadow')
        plt.legend()
        plt.minorticks_on()
        plt.grid(which='both', alpha=0.2)    
            
        plt.show()
        
    if(plot2D):
        
        extHZ = [zStart+zOffset, zFin+zOffset, xStart, xFin]
        extVZ = [zStart+zOffset, zFin+zOffset, yStart, yFin]
        
        plt.figure(figsize=(6,2))
        plt.subplots_adjust(0.13, 0.22, 0.97, 0.95)
        # plt.title('XZ')
        plt.imshow(histoH, extent=extHZ, origin='lower', aspect='auto', cmap=cmap)
        plt.xlabel('Z [mm]')
        plt.ylabel('Horizontal [mm]')
        plt.minorticks_on()
        plt.tick_params(which='both', axis='both', top=True, right=True)
        if(figprefix != ''):
            plt.savefig(figprefix + '_XZ.png', dpi=300)
    
        plt.figure(figsize=(6,2))
        plt.subplots_adjust(0.13, 0.22, 0.97, 0.95)    
        # plt.title('YZ')
        plt.imshow(histoV, extent=extVZ, origin='lower', aspect='auto', cmap=cmap)
        plt.xlabel('Z [mm]')
        plt.ylabel('Vertical [mm]')
        plt.minorticks_on()
        plt.tick_params(which='both', axis='both', top=True, right=True)
        if(figprefix != ''):
            plt.savefig(figprefix + '_YZ.png', dpi=300)
        

    
    return histoH, histoV, outdict

#%% DATA: 
   
d=31000  #distance of the mirror    
   
filename = 'SPU_8007eV_PAPU_002.h5'

histo, histoV, outdict = read_caustic(filename, write_attributes=False, plot=False, plot2D=False, print_minimum=False, cmap='viridis', figprefix='')

xi = outdict['xStart']
xf = outdict['xFin']
nx = outdict['nx']

x = np.linspace(xi, xf, nx)

yi = outdict['yStart']
yf = outdict['yFin']
ny = outdict['ny']

y = np.linspace(yi, yf, ny)

zi = outdict['zStart']
zf = outdict['zFin']
nz = outdict['nz']

z = np.linspace(zi, zf, nz) + d

# plt.imshow(histo)                       # show the caustic


#%% BEAM SIZE DATA:

    
fwhm_x, fwhm_y = np.genfromtxt('SPU_caustic_PAPU_31-61m_FWHM.txt', unpack=True, comments='#', usecols=(-2,-1))
fwtm_x, fwtm_y = np.genfromtxt('SPU_caustic_PAPU_31-61m_FWTM.txt', unpack=True, comments='#', usecols=(-2,-1))

fwhm_x_fit, fwhm_y_fit =  np.genfromtxt('SPU_caustic_PAPU_31-61m_FWHM.txt', unpack=True, comments='#', usecols=(1,2))
distance = np.genfromtxt('SPU_caustic_PAPU_31-61m_FWHM.txt', unpack=True, comments='#', usecols=(0))
    
#%% DEPTH OF FOCUS (DOF):

# Calculating the increase in 10% of the minimum fwhm
threshold_x = 1.1 * min(fwhm_x)  # determines the threshold_x
threshold_y = 1.1 * min(fwhm_y)  # determines the threshold_y

# Two lists with True in the values within the DOF
dof_x = fwhm_x <= threshold_x    
dof_y = fwhm_y <= threshold_y

# Horizontal
# Find the index of the first True from left to right 
ih0 = None
for index, value in enumerate(dof_x):
    if value:
        ih0 = index
        break

# Find the index of the first True from right to left
ih1 = None
for index, value in reversed(list(enumerate(dof_x))):
    if value:
        ih1 = index
        break
    
# Vertical
# Find the index of the first True from left to right 
iv0 = None
for index, value in enumerate(dof_y):
    if value:
        iv0 = index
        break

# Find the index of the first True from right to left
iv1 = None
for index, value in reversed(list(enumerate(dof_y))):
    if value:
        iv1 = index
        break


dof0_h = z[ih0] # horizontal depth of focus
dof1_h = z[ih1] # horizontal depth of focus 

dof0_v = z[iv0] # vertical depth of focus 
dof1_v = z[iv1] # vertical depth of focus 



#%% PLOT:


## Fit:
    
z_new, fwhm_h_new, popt_fwhm_h = fit_gaussian_beam(z, fwhm_x, auto_guess=True)
z_new, fwhm_v_new, popt_fwhm_v = fit_gaussian_beam(z, fwhm_y, auto_guess=True)

z_new, fwtm_h_new, popt_fwtm_h = fit_gaussian_beam(z, fwtm_x, auto_guess=True)
z_new, fwtm_v_new, popt_fwtm_v = fit_gaussian_beam(z, fwtm_y, auto_guess=True)


## Horizontal:
    
fs = 11
lw = 2

savefig = True

legend_fwhm = r'$FWHM_0 \sqrt{1+\left(\frac{z-z_0}{\beta}\right)^2}$'
legend_fwtm = r'$FWTM_0 \sqrt{1+\left(\frac{z-z_0}{\beta}\right)^2}$'

fig = plt.figure(figsize=(10,9))
ax1 = plt.subplot(211)
ax2 = plt.subplot(212, sharex = ax1)



ax1.imshow(histo, origin='lower', extent=(z[0]/1000, z[-1]/1000, -900, 900), interpolation = 'sinc', vmin=0, vmax=np.max(histo), aspect='auto')
ax1.plot(z_new/1000, -1*fwhm_h_new/2, '--', color='white', alpha=1.0, linewidth=lw, label='FWHM')
ax1.plot(z_new/1000, fwhm_h_new/2, '--', color='white', alpha=1.0, linewidth=lw)
ax1.plot(z_new/1000, -1*fwtm_h_new/2, '--', color='wheat', alpha=1.0, linewidth=lw, label='FWTM')
ax1.plot(z_new/1000, fwtm_h_new/2, '--', color='wheat', alpha=1.0, linewidth=lw)
# Plot D.O.F (Depth of Focus)
ax1.axvline(x=dof0_h/1000, color='red',alpha=0.5, linestyle='--',label="Depth of focus")
ax1.axvline(x=dof1_h/1000, color='red',alpha=0.5, linestyle='--')

ax2.plot(z_new/1000, fwhm_h_new, '-', color='C0', alpha=1.0, label=legend_fwhm) #label='Fit'
ax2.plot(z_new/1000, fwtm_h_new, '-', color='C1', alpha=1.0, label=legend_fwtm) #label='Fit' 


plt.text(0.8, 0.55, r'$z_0 = {0:.1f} \ m$'.format(np.abs(popt_fwhm_h[1]/1000))+'\n'+r'$FWHM_0 = {0:.1f} \ \mu m$'.format(np.abs(popt_fwhm_h[0]))+'\n'+r'$\beta = {0:.2f} \ m$'.format(np.abs(popt_fwhm_h[2]/1000)), 
          horizontalalignment='left', verticalalignment='bottom', fontsize = 10.0, color='C0', transform=plt.gca().transAxes)  

plt.text(0.8, 0.35, r'$z_0 = {0:.1f} \ m$'.format(np.abs(popt_fwtm_h[1]/1000))+'\n'+r'$FWTM_0 = {0:.1f} \ \mu m$'.format(np.abs(popt_fwtm_h[0]))+'\n'+r'$\beta = {0:.2f} \ m$'.format(np.abs(popt_fwtm_h[2]/1000)), 
          horizontalalignment='left', verticalalignment='bottom', fontsize = 10.0, color='C1', transform=plt.gca().transAxes)  

ax1.legend()
ax2.legend()

ax1.set_ylabel('Beam Size [$\mu$m]', fontsize=fs)
ax2.set_ylabel('Beam Size [$\mu$m]', fontsize=fs)

ax2.set_xlabel('Distance [m]', fontsize=fs) 

ax1.set_title('Horizontal Caustic', fontsize=fs+2) 

ax1.set_xlim(31, 61)
ax1.set_ylim(-900, 900)
ax1.minorticks_on()
ax1.tick_params(which='both', axis='both', direction='out', right=False, top=False, labelsize=fs)
ax1.set_xticks(np.linspace(31,61, round((61-31)/1+1)))
ax1.grid(which='both', alpha=0.2)

ax2.set_xlim(31, 61)
ax2.set_ylim(0, 1600)
ax2.minorticks_on()
ax2.tick_params(which='both', axis='both', direction='out', right=False, top=False, labelsize=fs)
ax2.set_xticks(np.linspace(31,61, round((61-31)/1+1)))
ax2.grid(which='both', alpha=0.2)

plt.tight_layout()
if(savefig): plt.savefig('Horizontal_caustic_31-61m_PAPU.png', dpi=500)


## Vertical:

fig = plt.figure(figsize=(10,9))
ax1 = plt.subplot(211)
ax2 = plt.subplot(212, sharex = ax1)

ax1.imshow(histoV, origin='lower', extent=(z[0]/1000, z[-1]/1000, -900, 900), interpolation = 'sinc', vmin=0, vmax=np.max(histoV), aspect='auto')
    
ax1.plot(z_new/1000, -1*fwhm_v_new/2, '--', color='white', alpha=1.0, linewidth=lw, label='FWHM')
ax1.plot(z_new/1000, fwhm_v_new/2, '--', color='white', alpha=1.0, linewidth=lw)
ax1.plot(z_new/1000, -1*fwtm_v_new/2, '--', color='wheat', alpha=1.0, linewidth=lw, label='FWTM')
ax1.plot(z_new/1000, fwtm_v_new/2, '--', color='wheat', alpha=1.0, linewidth=lw)
# Plot D.O.F (Depth of Focus)
ax1.axvline(x=dof0_v/1000, color='red',alpha=0.5, linestyle='--',label="Depth of focus")
ax1.axvline(x=dof1_v/1000, color='red',alpha=0.5, linestyle='--')

ax2.plot(z_new/1000, fwhm_v_new, '-', color='C0', alpha=1.0, label=legend_fwhm) #label='Fit'
ax2.plot(z_new/1000, fwtm_v_new, '-', color='C1', alpha=1.0, label=legend_fwtm) #label='Fit' 


plt.text(0.8, 0.55, r'$z_0 = {0:.1f} \ m$'.format(np.abs(popt_fwhm_v[1]/1000))+'\n'+r'$FWHM_0 = {0:.1f} \ \mu m$'.format(np.abs(popt_fwhm_v[0]))+'\n'+r'$\beta = {0:.2f} \ m$'.format(np.abs(popt_fwhm_v[2]/1000)), 
          horizontalalignment='left', verticalalignment='bottom', fontsize = 10.0, color='C0', transform=plt.gca().transAxes)  # 0.04, 0.05

plt.text(0.8, 0.35, r'$z_0 = {0:.1f} \ m$'.format(np.abs(popt_fwtm_v[1]/1000))+'\n'+r'$FWTM_0 = {0:.1f} \ \mu m$'.format(np.abs(popt_fwtm_v[0]))+'\n'+r'$\beta = {0:.2f} \ m$'.format(np.abs(popt_fwtm_v[2]/1000)), 
          horizontalalignment='left', verticalalignment='bottom', fontsize = 10.0, color='C1', transform=plt.gca().transAxes)  # 0.04, 0.05

ax1.legend()
ax2.legend()

ax1.set_ylabel('Beam Size [$\mu$m]', fontsize=fs)
ax2.set_ylabel('Beam Size [$\mu$m]', fontsize=fs)

ax2.set_xlabel('Distance [m]', fontsize=fs) 

ax1.set_title('Vertical Caustic', fontsize=fs+2) 

ax1.set_xlim(31, 61)
ax1.set_ylim(-550, 550)
ax1.minorticks_on()
ax1.tick_params(which='both', axis='both', direction='out', right=False, top=False, labelsize=fs)
ax1.set_xticks(np.linspace(31,61, round((61-31)/1+1)))
ax1.grid(which='both', alpha=0.2)




ax2.set_xlim(31, 61)
ax2.set_ylim(0, 1600)
ax2.minorticks_on()
ax2.tick_params(which='both', axis='both', direction='out', right=False, top=False, labelsize=fs)
ax2.set_xticks(np.linspace(31,61, round((61-31)/1+1)))
ax2.grid(which='both', alpha=0.2)

plt.tight_layout()
if(savefig): plt.savefig('Vertical_caustic_31-61m_PAPU.png', dpi=500)   