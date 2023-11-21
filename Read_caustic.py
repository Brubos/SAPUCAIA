#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 16:08:48 2023

@author: humberto.junior and Bruno.souza
"""

#%% Packages:

from optlnls.plot import plot_beam       # plot the beam information
import matplotlib.pyplot as plt          # library for creating charts and plots 
import optlnls.shadow                    # library related to optical simulations
import numpy as np                       # scientific computing library 
import h5py                              # library for handling HDF5 files


#%% Functions:
    
def gaussian_beam(z, s0, z0, beta):
    return s0*np.sqrt(1 + ((z-z0)/beta)**2)


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
        
    # if(plot2D == 'log'):
        
    #     xc_min_except_0 = np.min(histoH[histoH>0])
    #     histoH[histoH<=0.0] = xc_min_except_0/2.0
        
    #     plt.figure()
    #     plt.title('XZ')
    #     plt.imshow(histoH, extent=[zStart, zFin, xStart, xFin], origin='lower', aspect='auto',
    #                norm=LogNorm(vmin=xc_min_except_0/2.0, vmax=np.max(histoH)))
    #     plt.xlabel('Z')
    #     plt.ylabel('Horizontal')

    #     yc_min_except_0 = np.min(histoV[histoV>0])
    #     histoV[histoV<=0.0] = yc_min_except_0/2.0

    #     plt.figure()
    #     plt.title('YZ')
    #     plt.imshow(histoV, extent=[zStart, zFin, yStart, yFin], origin='lower', aspect='auto',
    #                norm=LogNorm(vmin=xc_min_except_0/2.0, vmax=np.max(histoV)))
    #     plt.xlabel('Z')
    #     plt.ylabel('Vertical')
    
    
    return histoH, histoV, outdict


#%% Data: 
    
filename = 'SPU_8007eV_PAPU.h5'
d = 31000   # distance from source

histoH, histoV, outdict = read_caustic(filename, write_attributes=False, plot=False, plot2D=False, print_minimum=False, cmap='viridis', figprefix='')

xi = outdict['xStart']
xf = outdict['xFin']
nx = outdict['nx']

x_array = np.linspace(xi, xf, nx)

yi = outdict['yStart']
yf = outdict['yFin']
ny = outdict['ny']

y_array = np.linspace(yi, yf, ny)

zi = outdict['zStart']
zf = outdict['zFin']
nz = outdict['nz']

z_array = np.linspace(zi, zf, nz) + d       

# plt.imshow(histoH)                    # show the caustic


#%% Run:

plot_FWHM = True
plot_FWTM = True    

fwhm_h_fit_list = []; fwhm_v_fit_list = []; fwhm_h_sli_list = []; fwhm_v_sli_list = [];
fwtm_h_fit_list = []; fwtm_v_fit_list = []; fwtm_h_sli_list = []; fwtm_v_sli_list = [];

with h5py.File(filename, 'r+') as f:
    
    dset_names = list(f.keys())
    
    for i in range(len(dset_names[2:])):
        
        print('Calculating for d = {} m'.format(d))
        dset = dset_names[i+2]
        
        beam = np.array(f[dset])
        
        beam_2D = np.zeros((ny+1, nx+1))
        
        beam_2D[0,1:] = x_array
        beam_2D[1:,0] = y_array
        beam_2D[1:,1:] = beam.transpose() # Não deveria ter o .transpose(). É um bug da nossa ferramenta "caustic" do Oasys.
            
        sr = 1.5 # size range
        
        d = z_array[i]
        
        title = 'd = %.d mm' %d
        
        if(plot_FWHM):
            
            outfilename_fwhm = '' #'SPU_beam_size_'+str(int(d))+'mm_FWHM.png'
            
            outputs_fwhm = plot_beam(beam_2D, outfilename=outfilename_fwhm, cut=0, textA=1, textB=5, textC=6, textD=10, fitType=3, 
                                     xlabel='Hor.', ylabel='Ver.', plot_title=title, unitFactor=1e3, fwhm_threshold=0.5, x_range=1, y_range=1, 
                                     x_range_min=-sr, x_range_max=sr, y_range_min=-sr, y_range_max=sr, cmap='viridis', zero_pad_x=1, zero_pad_y=5, 
                                     export_slices=0, zlabel='ph/s/100mA/$\mu$m', show_plot=False)
            
            fwhm_h_fit, fwhm_v_fit = outputs_fwhm['fit_fwhm_x'], outputs_fwhm['fit_fwhm_z']
            fwhm_h_sli, fwhm_v_sli = outputs_fwhm['fwhm_x'], outputs_fwhm['fwhm_z']
            
            fwhm_h_fit_list.append(fwhm_h_fit); fwhm_v_fit_list.append(fwhm_v_fit); 
            fwhm_h_sli_list.append(fwhm_h_sli); fwhm_v_sli_list.append(fwhm_v_sli);
            
        if(plot_FWTM):
            
            outfilename_fwtm = '' #'SPU_beam_size_'+str(int(d))+'mm_FWTM.png'
            
            outputs_fwtm = plot_beam(beam_2D, outfilename=outfilename_fwtm, cut=0, textA=1, textB=5, textC=6, textD=10, fitType=3, 
                                     xlabel='Hor.', ylabel='Ver.', plot_title=title, unitFactor=1e3, fwhm_threshold=0.1, x_range=1, y_range=1, 
                                     x_range_min=-sr, x_range_max=sr, y_range_min=-sr, y_range_max=sr, cmap='viridis', zero_pad_x=1, zero_pad_y=5, 
                                     export_slices=0, zlabel='ph/s/100mA/$\mu$m', show_plot=False)
            
            fwtm_h_fit, fwtm_v_fit = outputs_fwtm['fit_fwhm_x'], outputs_fwtm['fit_fwhm_z']
            fwtm_h_sli, fwtm_v_sli = outputs_fwtm['fwhm_x'], outputs_fwtm['fwhm_z']
            
            fwtm_h_fit_list.append(fwtm_h_fit); fwtm_v_fit_list.append(fwtm_v_fit); 
            fwtm_h_sli_list.append(fwtm_h_sli); fwtm_v_sli_list.append(fwtm_v_sli);
            

#%% Plot:
    
savefig = True

legend_fwhm = r'$FWHM_0 \sqrt{1+\left(\frac{z-z_0}{\beta}\right)^2}$'
legend_fwtm = r'$FWTM_0 \sqrt{1+\left(\frac{z-z_0}{\beta}\right)^2}$'

if(plot_FWHM):
    
    z_new, fwhm_h_new, popt_fwhm_h = fit_gaussian_beam(z_array, fwhm_h_fit_list, auto_guess=True)
    z_new, fwhm_v_new, popt_fwhm_v = fit_gaussian_beam(z_array, fwhm_v_fit_list, auto_guess=True)
    
    plt.figure(figsize=(4.5,3))
    plt.subplots_adjust(0.15, 0.15, 0.95, 0.95)
    plt.scatter(z_array/1000, fwhm_h_fit_list, color='C0', alpha=0.8, label='Shadow')
    plt.plot(z_new/1000, fwhm_h_new, '--', color='black', alpha=0.8, label=legend_fwhm) #label='Fit'
    # plt.title('Horizontal Caustic')
    # plt.ylabel('Beam Size FWHM [$\mu$m]')
    # plt.xlabel('Distance [m]')
    plt.title('Cáustica Horizontal')
    plt.ylabel('Tamanho do feixe FWHM [$\mu$m]')
    plt.xlabel('Distância [m]')
    plt.xlim(30, 60)
    plt.ylim(80, 480)
    plt.minorticks_on()
    plt.tick_params(which='both', axis='both', direction='in', right=True, top=True)
    plt.grid(which='both', alpha=0.2)
    plt.legend(fontsize=9) #10
    plt.text(0.04, 0.05, r'$z_0 = {0:.1f} \ m$'.format(np.abs(popt_fwhm_h[1]/1000))+'\n'+r'$FWHM_0 = {0:.1f} \ \mu m$'.format(np.abs(popt_fwhm_h[0]))+'\n'+r'$\beta = {0:.2f} \ m$'.format(np.abs(popt_fwhm_h[2]/1000)), 
             horizontalalignment='left', verticalalignment='bottom', fontsize = 8.0, color='C0', transform=plt.gca().transAxes)  # 0.04, 0.05
    plt.tight_layout()
    if(savefig): plt.savefig('Horizontal_Caustic_FWHM.png', dpi=600)
    plt.show()
    
    plt.figure(figsize=(4.5,3))
    plt.subplots_adjust(0.15, 0.15, 0.95, 0.95)
    plt.scatter(z_array/1000, fwhm_v_fit_list, color='C0', alpha=0.8, label='Shadow')
    plt.plot(z_new/1000, fwhm_v_new, '--', color='black', alpha=0.8, label=legend_fwhm) #label='Fit'
    # plt.title('Vertical Caustic')
    # plt.ylabel('Beam Size FWHM [$\mu$m]')
    # plt.xlabel('Distance [m]')
    plt.title('Cáustica Vertical')
    plt.ylabel('Tamanho do feixe FWHM [$\mu$m]')
    plt.xlabel('Distância [m]')
    plt.xlim(30, 60)
    plt.ylim(-20, 680)
    plt.minorticks_on()
    plt.tick_params(which='both', axis='both', direction='in', right=True, top=True)
    plt.grid(which='both', alpha=0.2)
    plt.legend(fontsize=9) #10
    plt.text(0.04, 0.05, r'$z_0 = {0:.1f} \ m$'.format(np.abs(popt_fwhm_v[1]/1000))+'\n'+r'$FWHM_0 = {0:.1f} \ \mu m$'.format(np.abs(popt_fwhm_v[0]))+'\n'+r'$\beta = {0:.2f} \ m$'.format(np.abs(popt_fwhm_v[2]/1000)), 
             horizontalalignment='left', verticalalignment='bottom', fontsize = 8.0, color='C0', transform=plt.gca().transAxes)  # 0.04, 0.05
    plt.tight_layout()
    if(savefig): plt.savefig('Vertical_Caustic_FWHM.png', dpi=600)
    plt.show()  
    

if(plot_FWTM):
    
    z_new, fwtm_h_new, popt_fwtm_h = fit_gaussian_beam(z_array, fwtm_h_fit_list, auto_guess=True)
    z_new, fwtm_v_new, popt_fwtm_v = fit_gaussian_beam(z_array, fwtm_v_fit_list, auto_guess=True)
    
    plt.figure(figsize=(4.5,3))
    plt.subplots_adjust(0.15, 0.15, 0.95, 0.95)
    plt.scatter(z_array/1000, fwtm_h_fit_list, color='C0', alpha=0.8, label='Shadow')
    plt.plot(z_new/1000, fwtm_h_new, '--', color='black', alpha=0.8, label=legend_fwtm) #label='Fit'
    # plt.title('Horizontal Caustic')
    # plt.ylabel('Beam Size FWTM [$\mu$m]')
    # plt.xlabel('Distance [m]')
    plt.title('Cáustica Horizontal')
    plt.ylabel('Tamanho do feixe FWTM [$\mu$m]')
    plt.xlabel('Distância [m]')
    plt.xlim(30, 60)
    plt.ylim(150, 900)
    plt.minorticks_on()
    plt.tick_params(which='both', axis='both', direction='in', right=True, top=True)
    plt.grid(which='both', alpha=0.2)
    plt.legend(fontsize=9) #10
    plt.text(0.04, 0.05, r'$z_0 = {0:.1f} \ m$'.format(np.abs(popt_fwtm_h[1]/1000))+'\n'+r'$FWTM_0 = {0:.1f} \ \mu m$'.format(np.abs(popt_fwtm_h[0]))+'\n'+r'$\beta = {0:.2f} \ m$'.format(np.abs(popt_fwtm_h[2]/1000)), 
             horizontalalignment='left', verticalalignment='bottom', fontsize = 8.0, color='C0', transform=plt.gca().transAxes)  # 0.04, 0.05
    plt.tight_layout()
    if(savefig): plt.savefig('Horizontal_Caustic_FWTM.png', dpi=600)
    plt.show()
    
    plt.figure(figsize=(4.5,3))
    plt.subplots_adjust(0.15, 0.15, 0.95, 0.95)
    plt.scatter(z_array/1000, fwtm_v_fit_list, color='C0', alpha=0.8, label='Shadow')
    plt.plot(z_new/1000, fwtm_v_new, '--', color='black', alpha=0.8, label=legend_fwtm) #label='Fit'
    # plt.title('Vertical Caustic')
    # plt.ylabel('Beam Size FWTM [$\mu$m]')
    # plt.xlabel('Distance [m]')
    plt.title('Cáustica Vertical')
    plt.ylabel('Tamanho do feixe FWTM [$\mu$m]')
    plt.xlabel('Distância [m]')
    plt.xlim(30, 60)
    plt.ylim(-50, 1100)
    plt.minorticks_on()
    plt.tick_params(which='both', axis='both', direction='in', right=True, top=True)
    plt.grid(which='both', alpha=0.2)
    plt.legend(fontsize=9) #10
    plt.text(0.04, 0.05, r'$z_0 = {0:.1f} \ m$'.format(np.abs(popt_fwtm_v[1]/1000))+'\n'+r'$FWTM_0 = {0:.1f} \ \mu m$'.format(np.abs(popt_fwtm_v[0]))+'\n'+r'$\beta = {0:.2f} \ m$'.format(np.abs(popt_fwtm_v[2]/1000)), 
             horizontalalignment='left', verticalalignment='bottom', fontsize = 8.0, color='C0', transform=plt.gca().transAxes)  # 0.04, 0.05
    plt.tight_layout()
    if(savefig): plt.savefig('Vertical_Caustic_FWTM.png', dpi=600)
    plt.show()  
    


#%% .txt:

if(plot_FWHM):
    
    filename_FWHM = 'SPU_caustic_PAPU_31-61m_FWHM.txt'
    
    data = np.array([z_array, fwhm_h_fit_list, fwhm_v_fit_list, fwhm_h_sli_list, fwhm_v_sli_list, fwhm_h_new, fwhm_v_new])
    
    np.savetxt(filename_FWHM, data.transpose(), '%.4f', delimiter='\t', header='Distance[mm]\tFWHM_X_Fit[um]\tFWHM_Y_Fit[um]\tFWHM_X_Slice[um]\tFWHM_Y_Slice[um]\tFWHM_X_Gaussian\tFWHM_Y_Gaussian')

if(plot_FWTM):
    
    filename_FWTM = 'SPU_caustic_PAPU_31-61m_FWTM.txt'
    
    data = np.array([z_array, fwtm_h_fit_list, fwtm_v_fit_list, fwtm_h_sli_list, fwtm_v_sli_list, fwtm_h_new, fwtm_v_new])
    
    np.savetxt(filename_FWTM, data.transpose(), '%.4f', delimiter='\t', header='Distance[mm]\tFWTM_X_Fit[um]\tFWTM_Y_Fit[um]\tFWTM_X_Slice[um]\tFWTM_Y_Slice[um]\tFWTM_X_Gaussian\tFWTM_Y_Gaussian')            
        
