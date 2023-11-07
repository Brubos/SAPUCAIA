#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 14:59:21 2023

@author: bruno.souza
"""

########### Simulate Sapucaia Beamline ###########
# this script obtain: 
    # the  beamline transmission;
    # the beam size and its divergence;
    # the total flux at the source and after the beamline;
    # the total power at the source and after the beamline;


## Packages:

import time                                    # package to work with time
import numpy as np                             # scientific computing library 
import sys                                     # library that allows indicating a path of a folder
sys.path.insert(0, '/home/ABTLUS/bruno.souza/miniconda3/envs/py38/lib/python3.8/site-packages/oasys_srw/')   
from scipy.signal import savgol_filter         # apply a Savitzky-Golay filter to an array
from optlnls.shadow import calc_und_flux       # used for calculate the flux undulator
from run_sapucaia import run_Sapucaia          # run sapucaia beamline
from optlnls.source import get_k, und_source   # obtain the K and the undulator source parameters (size and divergence)
from optlnls.importing import read_shadow_beam # used for read the shadow beam
from optlnls.math import get_fwhm              # obtain the fwhm from shadow beam
from optlnls.plot import plot_beam             # plot the beam information
import matplotlib.pyplot as plt                # library for creating charts and plots  


## Fixed Parameters:

eBeamEnergy = 3.0                                # storage ring energy [GeV]
e_spread = 0.00084                               # energy spread 
per = 22.0e-3                                    # undulator period [m]
L = 2.244                                        # undulator length [m]  
k_max = 1.71                                     # maximum deflection parameter
accept_hor,accept_ver = 66e-06, 66e-06           # beamline acceptance. Also used for calculating the undulator flux [rad]
source_beam = [65.3e-6, 3.0e-6, 3.8e-6, 1.4e-6]  # electron beam parameters [Sx [m], Sy [m], Sx' [rad], Sy'[rad]]


## User Defined parameters:

energy_array =np.linspace(5000, 18000, 2) # photon beam energy. Also used for calculating h, K, B [eV] 
atomic_plane = 'Si111'                    # crystal type. 'Si111' or 'Si311'
n_rays = 1000000                          # number of rays for shadow 
nbins = 151                               # bins for the undulator flux
distance_from_mirror=16850                # distance from mirror [mm]


plot_figure = True   # plots beam size and divergence
show_plots =  False  # show beam size and divergence plots 
plot_transm = True   # plots beamline transmission
show_transm = False  # show transmission plots


## Automatic Calculated parameters:

e_beam_size_X, e_beam_size_Z, e_beam_div_X, e_beam_div_Z = source_beam
betaX = e_beam_size_X/e_beam_div_X  
betaZ = e_beam_size_Z/e_beam_div_Z
emmX = e_beam_size_X*e_beam_div_X 
emmZ = e_beam_size_Z*e_beam_div_Z       


## Lists to store data:
    
flux_list = []; power_list = []; resolution_list = []; resolution_list_filt = []; harmonic_list = [];  

if(plot_figure):
    
    fwhm_h_size_fit_list = []; fwhm_v_size_fit_list = []; fwhm_h_size_sli_list = []; fwhm_v_size_sli_list = [];
    fwhm_h_div_fit_list = []; fwhm_v_div_fit_list = []; fwhm_h_div_sli_list = []; fwhm_v_div_sli_list = [];
    

## Run:

t0 = time.time()

for energy in energy_array:
    
    print('\n'+'E = %.1f eV' % energy)
    # print('\n')
    
    # Calculates beam size and divergence at source:
    
    h, K, B = get_k(per, 'max', energy, k_max)    

    size_X, div_X = und_source(emmX, betaX, e_spread, L, per, energy, h)
    size_Z, div_Z = und_source(emmZ, betaZ, e_spread, L, per, energy, h)
    
    delta_E =(11/13)*(energy/1000) - (55/13) + 3    
    
    harmonic_list.append(h)
    

    
    # Run Shadow and calculates flux:
    
    beam = run_Sapucaia(n_rays=n_rays, energy=energy, delta_E=delta_E, sig_h=size_X/1000, sig_v=size_Z/1000, div_h=1e-06*div_X, div_v=1e-06*div_Z, d_image=distance_from_mirror, 
                    atomic_plane=atomic_plane, hor_accept=accept_hor, ver_accept=accept_ver)
    
    
    outputs = calc_und_flux(beam=beam, nbins=nbins, eBeamEnergy=eBeamEnergy, eSpread=e_spread, current=0.1,
                            und_per=per, und_length=L, B=B, min_harmonic=1, max_harmonic=(h+2),
                            source_beam=source_beam, show_plots=False, accept_hor=accept_hor, accept_ver=accept_ver)
    
    flux_s, flux_b = outputs['total flux at source'], outputs['total flux propagated']
    
    power_s, power_b = outputs['total power at source'], outputs['total power propagated']
    
    E_b, T_E = outputs['energy array'] , outputs['transmission array']
    
    T_E_filt = savgol_filter(T_E, 15, 3)
    
    resol = get_fwhm(E_b, T_E)[0]/energy
    
    resol_filt = get_fwhm(E_b, T_E_filt)[0]/energy
    
    flux_list.append(flux_b); power_list.append(power_b)
    
    resolution_list.append(resol); resolution_list_filt.append(resol_filt);
    
    
    # Plot Transmission:
        
    if(plot_transm):
        
        fig = plt.figure()
        plt.plot(E_b, T_E, label='data')
        plt.plot(E_b, T_E_filt, label='filter')
        plt.xlabel('Energy [eV]', fontsize=13)
        plt.ylabel('Transmission', fontsize=13)
        plt.legend(fontsize=11)
        plt.yscale('linear')
        plt.xscale('linear')
        plt.minorticks_on()
        plt.tick_params(which='both', axis='both', direction='in', right=True, top=True, labelsize=12)
        plt.grid(which='both', alpha=0.2)    
        plt.tight_layout()
        plt.savefig('transmission_'+atomic_plane+'_'+str(int(energy))+'eV.png')
        if(not(show_transm)): plt.close(fig)
    
    
    # Plots beam size and divergence at sample:
    
    if(plot_figure):
    
        beam2D_size = read_shadow_beam(beam, x_column_index=3, y_column_index=1, nbins_x=150, nbins_y=150, nolost=1, 
                                        ref=23, zeroPadding=4.5, gaussian_filter=0)
        
        beam2D_div = read_shadow_beam(beam, x_column_index=6, y_column_index=4, nbins_x=150, nbins_y=150, nolost=1, 
                                      ref=23, zeroPadding=1.2, gaussian_filter=0)
        
        if((energy/1000) < 10):
            title = 'SPU - %.1f keV            F: %.2E ph/s/100mA' %(energy/1000, flux_b)
        else:
            title = 'SPU - %.1f keV           F: %.2E ph/s/100mA' %(energy/1000, flux_b)
        
        
        sr = 250e-3 # size range
        dr = 50e-6  # dive range
                
        
        outputs_beam = plot_beam(beam2D_size, outfilename='SPU_beam_size_'+atomic_plane+'_'+str(int(energy))+'eV.png',
                                  cut=0, textA=1, textB=5, textC=6, textD=10, fitType=3, xlabel='Hor.', ylabel='Ver.', plot_title=title, unitFactor=1e3,
                                  fwhm_threshold=0.5, x_range=1, y_range=1, x_range_min=-sr, x_range_max=sr, y_range_min=-sr, y_range_max=sr, cmap='viridis',
                                  integral=flux_b, zero_pad_x=1, zero_pad_y=5, export_slices=0, zlabel='ph/s/100mA/$\mu$m', show_plot=show_plots)
        
        outputs_div = plot_beam(beam2D_div, outfilename='SPU_beam_div_'+atomic_plane+'_'+str(int(energy))+'eV.png', cut=0,
                                textA=1, textB=5, textC=6, textD=10, fitType=3, xlabel='Hor.', ylabel='Ver.', plot_title=title, unitFactor=1e6,
                                fwhm_threshold=0.1, x_range=1, y_range=1, x_range_min=-dr, x_range_max=dr, y_range_min=-dr, y_range_max=dr,
                                cmap='viridis', integral=flux_b, units='$\mu$rad', zero_pad_x=1, zero_pad_y=1, export_slices=0, 
                                zlabel='ph/s/100mA/$\mu$rad', show_plot=show_plots)
    
    
        fwhm_h_size_fit, fwhm_v_size_fit = outputs_beam['fit_fwhm_x'], outputs_beam['fit_fwhm_z']
        fwhm_h_size_sli, fwhm_v_size_sli = outputs_beam['fwhm_x'], outputs_beam['fwhm_z']
        
        fwhm_h_div_fit, fwhm_v_div_fit = outputs_div['fit_fwhm_x'], outputs_div['fit_fwhm_z']
        fwhm_h_div_sli, fwhm_v_div_sli = outputs_div['fwhm_x'], outputs_div['fwhm_z']
        
        fwhm_h_div_fit_list.append(fwhm_h_div_fit); fwhm_v_div_fit_list.append(fwhm_v_div_fit); 
        fwhm_h_div_sli_list.append(fwhm_h_div_sli); fwhm_v_div_sli_list.append(fwhm_v_div_sli);
        
        fwhm_h_size_fit_list.append(fwhm_h_size_fit); fwhm_v_size_fit_list.append(fwhm_v_size_fit); 
        fwhm_h_size_sli_list.append(fwhm_h_size_sli); fwhm_v_size_sli_list.append(fwhm_v_size_sli);
    

        
## Writing Flux .txt file:

filename = 'SPU_flux_'+atomic_plane+'.txt'

data = np.array([energy_array, flux_list, power_list, resolution_list, resolution_list_filt, harmonic_list])

np.savetxt(filename, data.transpose(), '%.6E', delimiter='\t', header='Energy[eV]\tFlux[ph/s/100mA]\tPower[W/100mA]\tResolution\tResolution_Filtered\tHarmonic')


## Writing Size .txt file:
    
if(plot_figure):
    
    filename = 'SPU_size_'+atomic_plane+'.txt'
    
    data = np.array([energy_array, fwhm_h_size_fit_list, fwhm_v_size_fit_list, fwhm_h_size_sli_list, fwhm_v_size_sli_list])
    
    np.savetxt(filename, data.transpose(), '%.3f', delimiter='\t', header='Energy[eV]\tFWHM_X_Fit[\u03bcm]\tFWHM_Y_Fit[\u03bcm]\tFWHM_X_Slice[\u03bcm]\tFWHM_Y_Slice[\u03bcm]')


## Writing Divergence .txt file:
    
if(plot_figure):
    
    filename = 'SPU_divergence_'+atomic_plane+'.txt'
    
    data = np.array([energy_array, fwhm_h_div_fit_list, fwhm_v_div_fit_list, fwhm_h_div_sli_list, fwhm_v_div_sli_list])
    
    np.savetxt(filename, data.transpose(), '%.3f', delimiter='\t', header='Energy[eV]\tFWHM_X_Fit[\u03bcrad]\tFWHM_Y_Fit[\u03bcrad]\tFWHM_X_Slice[\u03bcrad]\tFWHM_Y_Slice[\u03bcrad]')


print('Total time = %.2f s' %(time.time()-t0))

