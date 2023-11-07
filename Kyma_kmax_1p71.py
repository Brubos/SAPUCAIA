#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 11:01:32 2023

@author: bruno.souza
"""

########### Calculate R.M.S beam informations ###########
# this script obtain: 
    # a graph of RMS beam size;
    # a graph of RMS beam divergence;
    # a file.txt with informations (E, h, K, B, σx, σz, σ'x, σ'z)				 
    


# Packages:  
from optlnls.source import get_k, und_source  # obtain the K and the undulator source parameters (size and divergence)
import numpy as np                            # scientific computing library
import matplotlib.pyplot as plt               # library for creating charts and plots  

# Parameters:
per = 22e-03                                      # undulator period [m]
k_max = 1.71                                      # maximum deflection parameter
e_spread = 0.084/100                              # energy spread
N_per = 51                                        # number of periods 
L = N_per*per                                     # undulator length [m]  
e_beam_size_h, e_beam_size_v = [64.7e-6, 5.1e-6]  # electron beam size [ Sx [m], Sy [m] ]
e_beam_div_h, e_beam_div_v = [3.8e-6, 1.4e-6]     # electron beam divergence [ Sx' [rad], Sy'[rad] ]
Eo = 1000                                         # initial energy [eV]
Ef = 37000                                        # final energy   [eV]
steps=300                                         # energy steps   [eV]
n = round((Ef-Eo)/(steps)+1)                      # point numbers           
energies=np.linspace(Eo,Ef,n)                     # create energy vector

# Automatically calculated parameters:
emm_h = (e_beam_size_h)*(e_beam_div_h) 
emm_v = (e_beam_size_v)*(e_beam_div_v)

beta_h = (e_beam_size_h)/(e_beam_div_h) 
beta_v = (e_beam_size_v)/(e_beam_div_v)

# Create seven null vectors
vetor_h=(np.zeros(shape=n))   # h - harmonic number [dimensionless]
vetor_K=np.zeros(shape=n)     # K - deflection parameter [dimensionless]
vetor_B=np.zeros(shape=n)     # B - magnetic Field [T]
vetor_szh=np.zeros(shape=n)   # horizontal beam size [mm]
vetor_szv=np.zeros(shape=n)   # vertical beam size  [mm]
vetor_divh=np.zeros(shape=n)  # horizontal divergence [μrad]
vetor_divv=np.zeros(shape=n)  # vertical divergence [μrad]

# Create an unitary vector
uni=(np.ones(shape=n))
uni_normalize=uni*(1e6)

for i in range(len(energies)):
    
    energy = energies[i] 
    
    # Calculate harmonic number and K:
    h, K, B = get_k(Period=per, what_harmonic='max', Energy=energy, k_ext=k_max)
    
    # Calculate beam size and divergence:
    size_h, div_h = und_source(emittance=emm_h, beta=beta_h, e_spread=e_spread, und_length=L, und_period=per, 
                            ph_energy=energy, harmonic=h) # [um], [urad]
    size_v, div_v = und_source(emittance=emm_v, beta=beta_v, e_spread=e_spread, und_length=L, und_period=per, 
                            ph_energy=energy, harmonic=h) # [um], [urad]
    
    # Save the values in a vector    
    vetor_h[i]=h
    vetor_K[i]=K
    vetor_B[i]=B  
    vetor_szh[i]=size_h
    vetor_szv[i]=size_v
    vetor_divh[i]=div_h
    vetor_divv[i]=div_v
    
vetor_data=np.array([energies,vetor_h,vetor_K,vetor_B,vetor_szh,vetor_szv,vetor_divh,vetor_divv]).transpose()


# Save .txt
filename = 'Sapucaia_sizes.txt'
np.savetxt(filename, vetor_data, fmt='%.6f', delimiter='\t', header="Energy (eV)\tHarmonic\tK\t\tB(T)\t\t\u03C3x\t\t\u03C3z\t\t\u03C3x'\t\t\u03C3z'")

# PLOT RMS BEAM SIZE
plt.figure()
plt.title("KYMA 2.244 m / high-$\u03B2$ SA")
plt.xlabel("Energy [keV]")
plt.ylabel("RMS Beam Size [$\mu$m]")
plt.xlim(0, 38)
plt.ylim(0, 80)
plt.plot(energies/1000,e_beam_size_h*uni_normalize, ':', label='Horizontal / e-beam')
plt.plot(energies/1000,e_beam_size_v*uni_normalize, ':', label='Vertical / e-beam')
plt.plot(energies/1000,vetor_szh, '.', label='Vertical / photon beam')
plt.plot(energies/1000,vetor_szv, '.', label='Horizontal / photon beam')
plt.legend(loc="center right")

# Save png 
plt.savefig('RMS_beam_size', dpi=600)

# PLOT RMS BEAM DIVERGENCE
plt.figure()
plt.title("KYMA 2.244 m / high-$\u03B2$ SA")
plt.xlabel("Energy [keV]")
plt.ylabel("RMS Beam Divergence [$\mu$rad]")
plt.xlim(0, 38)
plt.ylim(0, 20)
plt.plot(energies/1000,e_beam_div_h*uni_normalize, ':', label='Horizontal / e-beam')
plt.plot(energies/1000,e_beam_div_v*uni_normalize, ':', label='Vertical / e-beam')
plt.plot(energies/1000,vetor_divh, '.', label='Horizontal / photon beam')
plt.plot(energies/1000,vetor_divv, '.', label='Vertical / photon beam')
plt.legend()

# Save png 
plt.savefig('RMS_beam_divergence', dpi=600)

