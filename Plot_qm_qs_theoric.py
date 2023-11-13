#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 20:01:05 2023

@author: bruno
"""

import matplotlib.pyplot as plt
import numpy as np


#%% Parameters
Rm = 7908.179411                         # Meridional radius [m]
Rs = 0.096874                            # Sagittal radius   [m]
p = 31                                   # Distance source-optic element [m]
gamma_list = []                          # Parameter for calculating the meridional radius [adim]
beta_list = []                           # Parameter for calculating the sagittal radius   [adim]
qm_list = []                             # q meridional list [m]
qs_list = []                             # q sagittal list   [m]
pitch = np.linspace(2.5, 4.5, 21)        # Mirror angle [mrad]
delta_pitch = np.linspace(-1, 1, 21)     # Variation of mirror angle [mrad]


#%% Functions

def beta(Rs,pitch,p):
    for angle in pitch:
        beta=(Rs/(2*angle*p))*1000
        beta_list.append(beta)
    return beta_list

def gamma(Rm,pitch,p):
    for angle in pitch:
        gamma=((Rm*angle)/(2*p))/1000
        gamma_list.append(gamma)
    return gamma_list

def qm(gamma_list,p):
    for gamma in gamma_list:
        qm=(gamma*p)/(1-gamma)
        qm_list.append(qm+p)
    return qm_list

def qs(betta_list,p):
    for beta in beta_list:
        qs=(beta*p)/(1-beta)
        qs_list.append(qs+p)
    return qs_list


def plot_graphs(delta_qm, delta_qs, pitch):
    # Plotting delta_qm as a function of pitch
    plt.figure(figsize=(12, 5))

    plt.subplot(1, 2, 1)
    plt.plot(pitch, delta_qm, marker='o', linestyle='-', color='b')
    plt.title(r'$q_m = \frac{\gamma p}{(1-\gamma)}$' +r'        '+ r'$\gamma = \frac{R_m \theta}{2p}$', pad=20, fontsize=18)
    plt.xlabel('ΔPitch [mrad]')
    plt.ylabel('Δq Meridional')
    plt.grid(True)

    # Plotting delta_qs as a function of pitch
    plt.subplot(1, 2, 2)
    plt.plot(pitch, delta_qs, marker='o', linestyle='-', color='r')
    plt.title(r'$q_s = \frac{\beta p}{(1-\beta)}$' +r'        '+ r'$\beta = \frac{R_s}{2\theta p}$', pad=20, fontsize=18)
    plt.xlabel('ΔPitch [mrad]')
    plt.ylabel('Δq Sagittal')
    plt.grid(True)

    # Adjusting the layout to avoid label overlap
    plt.tight_layout()

    # Saving the plots with high quality
    plt.savefig('Δqm_Δqs_pitch', dpi=1000)

    # Displaying the plots
    plt.show()

#%% Calculating

beta_list=beta(Rs,pitch,p)
gamma_list=gamma(Rm,pitch,p)
qm_list=np.array(qm(gamma_list,p))
qs_list=np.array(qs(beta_list,p))
delta_qm=qm_list - 56
delta_qs=qs_list - 56
plot_graphs(delta_qm, delta_qs, delta_pitch)
