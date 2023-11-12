#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 14:41:02 2023

@author: bruno.souza
"""
########### Allows you to image the beams of a misalignment vector  ###########
# this script obtain:
    # the focal sensibility


#%% PACKAGES:
import numpy as np                                                           # Scientific computing library 
import Shadow                                                                # Used for shadow simulations
from optlnls.importing import read_shadow_beam                               # Used for read the shadow beam
from Plot_beam_modified import plot_beam                                     # Plot the beam information
from optlnls.shadow import run_shadow_caustic, read_caustic,get_good_ranges  # Custom functions related to Shadow software
import matplotlib.pyplot as plt                                              # Library for creating charts and plots 
import pandas as pd                                                          # Library for data manipulation and analysis

#%% PARAMETERS:
n_rays= 1000000                               # number of rays used in the simulation
misalig_mirror=np.linspace(2.5,4.5,21)      # mirror misalignment array;

# caustic parameters
z0 = 0           # starting z
zf = 30000       # final z
nz = 11          # number of points


# Set the default values for nbins_x and nbins_y
nbins_x = 400
nbins_y = 400


#%% Lists to store data:

# Initialize a dictionary to store data
data = {
    'misalignment': [],
    'fwhm_h_sli': [],
    'fwhm_v_sli': [],
    'mean_pos_h': [],
    'mean_pos_v': [],
    'z_fwhm_min_h': [],
    'z_fwhm_min_v': [],
    'z_rms_min_h':[],
    'z_rms_min_v':[],}

#%% Running the SAPUCAIA beamline to the DCM
# RUN SAPUCAIA:

# write (1) or not (0) SHADOW files start.xx end.xx star.xx
iwrite = 0

#
# initialize shadow3 source (oe0) and beam
#
beam = Shadow.Beam()
oe0 = Shadow.Source()
oe1 = Shadow.OE()
oe2 = Shadow.OE()
oe3 = Shadow.OE()
oe4 = Shadow.OE()

#
# Define variables. See meaning of variables in: 
#  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml 
#  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
#

oe0.FDISTR = 3
# oe0.FILE_BOUND = b'/home/ABTLUS/bruno.souza/GITHUB/SAPUCAIA/SPU_optimize_source_66x66urad2.txt'  # CNPEM
oe0.FILE_BOUND = b'/home/bruno/GITHUB/SAPUCAIA/SPU_optimize_source_66x66urad2.txt'             # Home
oe0.F_BOUND_SOUR = 2
oe0.F_PHOT = 0
oe0.HDIV1 = 0.0
oe0.HDIV2 = 0.0
oe0.IDO_VX = 0
oe0.IDO_VZ = 0
oe0.IDO_X_S = 0
oe0.IDO_Y_S = 0
oe0.IDO_Z_S = 0
oe0.ISTAR1 = 5676561
oe0.NPOINT = n_rays
oe0.PH1 = 8007.0
oe0.SIGDIX = 6.3e-06
oe0.SIGDIZ = 1.17e-05
oe0.SIGMAX = 0.0647
oe0.SIGMAZ = 0.0051
oe0.VDIV1 = 0.0
oe0.VDIV2 = 0.0

oe1.DUMMY = 0.1
oe1.FWRITE = 3
oe1.F_REFRAC = 2
oe1.F_SCREEN = 1
oe1.I_SLIT = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
oe1.K_SLIT = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
oe1.N_SCREEN = 1
oe1.RX_SLIT = np.array([3.1125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
oe1.RZ_SLIT = np.array([3.1125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
oe1.T_IMAGE = 0.0
oe1.T_INCIDENCE = 0.0
oe1.T_REFLECTION = 180.0
oe1.T_SOURCE = 20750.0

oe2.DUMMY = 0.1
# oe2.FILE_REFL = b'/home/ABTLUS/bruno.souza/Oasys/Si111.dat'   # CNPEM
oe2.FILE_REFL = b'/home/bruno/Oasys/Si111.dat'              # HOME
oe2.FWRITE = 1
oe2.F_CENTRAL = 1
oe2.F_CRYSTAL = 1
oe2.PHOT_CENT = 8007.0
oe2.R_LAMBDA = 5000.0
oe2.T_IMAGE = 37.608
oe2.T_INCIDENCE = 75.7024111803
oe2.T_REFLECTION = 75.7024111803
oe2.T_SOURCE = 8250.0

oe3.ALPHA = 180.0
oe3.DUMMY = 0.1
# oe3.FILE_REFL = b'/home/ABTLUS/bruno.souza/Oasys/Si111.dat'  # CNPEM
oe3.FILE_REFL = b'/home/bruno/Oasys/Si111.dat'              # HOME
oe3.FWRITE = 1
oe3.F_CENTRAL = 1
oe3.F_CRYSTAL = 1
oe3.PHOT_CENT = 8007.0
oe3.R_LAMBDA = 5000.0
oe3.T_IMAGE = 1962.392
oe3.T_INCIDENCE = 75.7024111803
oe3.T_REFLECTION = 75.7024111803
oe3.T_SOURCE = 0.0
 
#Run SHADOW to create the source

if iwrite:
    oe0.write("start.00")

beam.genSource(oe0)

if iwrite:
    oe0.write("end.00")
    beam.write("begin.dat")

############################################################### Add Bandwidth:

#A2EV = 50676.89919462:

codata_h = np.array(6.62606957e-34)     # Planck constant

codata_ec = np.array(1.602176565e-19)   # elementary charge

codata_c = np.array(299792458.0)        # speed of light

A2EV = 2.0*np.pi/(codata_h*codata_c/codata_ec*1e2)
 
#Shadow beam:

E_old = beam.getshonecol(11, nolost=1)  # energy column # beam.rays[:,10]/A2EV

E0 = E_old[0]
 
#Energy Bandwidth: 

delta_E = 6 #[eV]
 
#New beam:

beam_new = beam.duplicate()

E_new = E0 + delta_E*( np.random.random_sample((len(E_old),)) - 0.5 )

beam_new.rays[:,10] = E_new*A2EV
 
beam = beam_new.duplicate()
 
#Prints:

print('\n')

print('Energy (old): \n', E_old)

print('Energy (new): \n', beam_new.getshonecol(11, nolost=1))

print('\n')
 
print('New Energy Limits:')

print('E min = ', np.min(beam_new.getshonecol(11, nolost=1)))

print('E max = ', np.max(beam_new.getshonecol(11, nolost=1)))

print('\n')

##############################################################################

#
# Run optical element 1
#
print("    Running optical element: %d"%(1))
if iwrite:
    oe1.write("start.01")

beam.traceOE(oe1,1)

if iwrite:
    oe1.write("end.01")
    beam.write("star.01")


#
# Run optical element 2
#
print("    Running optical element: %d"%(2))
if iwrite:
    oe2.write("start.02")

beam.traceOE(oe2,2)

if iwrite:
    oe2.write("end.02")
    beam.write("star.02")


#
# Run optical element 3
#
print("    Running optical element: %d"%(3))
if iwrite:
    oe3.write("start.03")

beam.traceOE(oe3,3)

if iwrite:
    oe3.write("end.03")
    beam.write("star.03")


#%% RUNNING ONLY MIRROR FOR DIFFERENT MISALIGNMENT VALUES 

for misalig in misalig_mirror:    
    
    print("####################################################################")  
    print(misalig)
    print("####################################################################") 
    
    beam_copy = beam.duplicate()   # make a copy of the beam
    
    #
    # Run optical element 4 - Mirror
    #
    oe4.ALPHA = 90.0
    oe4.DUMMY = 0.1
    oe4.FHIT_C = 1 
    # oe4.FILE_REFL = b'/home/ABTLUS/bruno.souza/Oasys/Rh.dat'                                          # CNPEM
    # oe4.FILE_RIP = b'/home/ABTLUS/bruno.souza/GITHUB/SAPUCAIA/SPU_total_deformation_300mm_sh.dat'     # CNPEM
    oe4.FILE_REFL = b'/home/bruno/Oasys/Rh.dat'                                          # HOME
    oe4.FILE_RIP = b'/home/bruno/GITHUB/SAPUCAIA/SPU_total_deformation_300mm_sh.dat'     # HOME
    oe4.FMIRR = 3
    oe4.FWRITE = 1
    oe4.F_DEFAULT = 0
    oe4.F_G_S = 2
    oe4.F_MOVE = 1
    oe4.F_REFLEC = 1
    oe4.F_RIPPLE = 1
    oe4.RLEN1 = 150.0
    oe4.RLEN2 = 150.0
    oe4.RWIDX1 = 2.5
    oe4.RWIDX2 = 2.5
    oe4.SIMAG = 25000.0
    oe4.SSOUR = 31000.0
    oe4.THETA = 90-np.degrees(misalig*1e-3)
    oe4.T_IMAGE = 0.0
    oe4.T_INCIDENCE =  90-np.degrees(misalig*1e-3)
    oe4.T_REFLECTION = 90-np.degrees(misalig*1e-3)
    oe4.T_SOURCE = 0.0
    
    #
    # Run optical element 4
    #
    print("    Running optical element: %d"%(4))
    if iwrite:
        oe4.write("start.04")
    
    beam_copy.traceOE(oe4,4)
    
    if iwrite:
        oe4.write("end.04")
        beam_copy.write("star.04")
   
      
   #%% Read Shadow beam
    beam2D = read_shadow_beam(beam_copy, x_column_index=3, y_column_index=1, nbins_x=nbins_x, nbins_y=nbins_y, nolost=1, ref=23, zeroPadding=2, gaussian_filter=0)
        
    #%% Plot Beam
    
    zero_pad_x=0
    zero_pad_y=0
    plot_range_x=3000
    plot_range_y=3000
    
    outputs = plot_beam(beam2D, show_plot=False,outfilename='',outfileext='png',cut=0,textA=1,textB=5,textC=2,fitType=3,cmap='viridis',plot_title='',
                        zero_pad_x=zero_pad_x, zero_pad_y=zero_pad_y,
                        x_range = plot_range_x, y_range = plot_range_y,
                        x_range_min= -(plot_range_x/2)/1000,
                        x_range_max=  (plot_range_x/2)/1000,
                        y_range_min= -(plot_range_y/2)/1000,
                        y_range_max=  (plot_range_x/2)/1000)
    
    #%% Run Caustic
    

    print("####################################################################")
    print('Running caustic...')
    print("####################################################################") 


    Caustic = "Caustic (z0="+str(z0)+", zf="+str(zf)+", nz="+str(nz)+')'
    
    # Get the good range for x and y
    goodRange = get_good_ranges(beam=beam_copy, zStart=z0, 
                                zFin=zf, colh=3, colv=1)
    
    # Ru caustic
    run_shadow_caustic(filename=Caustic, beam=beam_copy, zStart=z0, zFin=zf, nz=nz, zOffset=0, colh=3, colv=1, colref=23, nbinsh=200, nbinsv=200,
                        xrange=goodRange[0:2], 
                        yrange=goodRange[2:4])

    histo_h, histo_v, caustic_dict = read_caustic(Caustic,plot=False,plot2D=False,cmap='viridis',figprefix='',print_minimum=False)
    
    #%% Store the outputs   
    data['misalignment'].append(misalig)
    data['fwhm_h_sli'].append(outputs['fwhm_x'])
    data['fwhm_v_sli'].append(outputs['fwhm_z'])
    data['mean_pos_h'].append(outputs['mean_x'])
    data['mean_pos_v'].append(outputs['mean_z'])
    data['z_fwhm_min_h'].append(caustic_dict['z_fwhm_min_h']*(1e-3)+31)   # sum 31000 due position of the mirror
    data['z_fwhm_min_v'].append(caustic_dict['z_fwhm_min_v']*(1e-3)+31)   # sum 31000 due position of the mirror
    data['z_rms_min_h'].append(caustic_dict['center_rms_h'])
    data['z_rms_min_v'].append(caustic_dict['center_rms_v'])   
         
#%% Save the data

# Original header list
header_list = [' Misalignment [mrad]', 
               'z_fwhm_min_h [m]',
               'z_fwhm_min_v [m]']
               # 'z_rms_min_h []',
               # 'z_rms_min_v []',]

# Custom formatting function to center-align header
def format_center_align_header(value):
    return f'{value:^20}'  # You can adjust the width as needed

# Apply the formatting function to each header
formatted_headers = [format_center_align_header(header) for header in header_list]

# Create a DataFrame with the formatted headers
data2 = pd.DataFrame({
    formatted_headers[0]: data['misalignment'],
    formatted_headers[1]: data['z_fwhm_min_h'],
    formatted_headers[2]: data['z_fwhm_min_v'],
    # formatted_headers[3]: data['z_rms_min_h'],  
    # formatted_headers[4]: data['z_rms_min_v']
})

# Specify the file name
file_name = 'teste'

# Custom formatting function to center-align data
def format_center_align_data(value):
    return f'{value:^20}'  # You can adjust the width as needed

# Apply the formatting function to each data column
for column in data2.columns:
    data2[column] = data2[column].apply(format_center_align_data)

# Save the DataFrame to a .txt file with a tab separator
data2.to_csv(file_name, sep='\t', index=False, header=False)

# Save the header separately
with open(file_name, 'r') as f:
    content = f.read()

with open(file_name, 'w') as f:
    f.write('\t'.join(formatted_headers) + '\n')
    f.write(content)

print(f'Data saved to {file_name}')
