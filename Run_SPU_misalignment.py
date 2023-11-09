#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 14:41:02 2023

@author: bruno.souza
"""
########### Allows you to image the beams of a misalignment vector  ###########
# this script obtain:
    # the misalignment of the beam
    # the values of FWHM slices and fit
    # the values of mean position and mean position fit
    

#%% PACKAGES:
import numpy as np                                                           # Scientific computing library 
import Shadow                                                                # Used for shadow simulations
from optlnls.importing import read_shadow_beam                               # Used for read the shadow beam
from Plot_beam_modified import plot_beam                                     # Plot the beam information
from optlnls.shadow import run_shadow_caustic, read_caustic,get_good_ranges  # Custom functions related to Shadow software
import matplotlib.pyplot as plt                                              # Library for creating charts and plots 
import pandas as pd                                                          # Library for data manipulation and analysis

#%% PARAMETERS:
n_rays= 1000000                            # number of rays used in the simulation
deg_f='Ry'                                 # the degree of freedom in Sirius frame;
unit ='µrad'                               # µrad, mrad for rotations; [µm], [mm] for translations
misalig_array=[54] #np.linspace(-90,90,11)      # misalignment array;

position = 16970                           # Position of interest

# Define a dictionary to map positions to devices
position_to_device = {
    841:   'DVF3',
    13750: 'DVF4',
    16970: 'SMP',
    22000: 'Hor_f',
    25000: 'Ver_f'    }

# Set the default values for nbins_x and nbins_y
nbins_x = 200
nbins_y = 200

# Check if the position is in the dictionary
if position in position_to_device:
    device = position_to_device[position]
    
    # Update nbins_x based on the device
    if device == 'DVF3':
        nbins_x = 90        # better parameters for this case

#%% Lists to store data:
    
# Initialize a dictionary to store data
data = {
    'misalignment': [],
    'fwhm_h_sli': [],
    'fwhm_v_sli': [],
    'mean_pos_h': [],
    'mean_pos_v': [],
    'fwhm_h_fit': [],
    'fwhm_v_fit': [],
    'mean_pos_h_fit': [],
    'mean_pos_v_fit': [],   }

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
oe0.FILE_BOUND = b'/home/ABTLUS/bruno.souza/GITHUB/SAPUCAIA/SPU_optimize_source_66x66urad2.txt'  # CNPEM
# oe0.FILE_BOUND = b'/home/bruno/GITHUB/SAPUCAIA/SPU_optimize_source_66x66urad2.txt'                 # Home
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
oe2.FILE_REFL = b'/home/ABTLUS/bruno.souza/Oasys/Si111.dat'   # CNPEM
# oe2.FILE_REFL = b'/home/bruno/Oasys/Si111.dat'              # HOME
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
oe3.FILE_REFL = b'/home/ABTLUS/bruno.souza/Oasys/Si111.dat'  # CNPEM
# oe3.FILE_REFL = b'/home/bruno/Oasys/Si111.dat'              # HOME
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

for misalig in misalig_array:    
    
    print("####################################################################")
    
    if deg_f in ['Rx', 'Ry', 'Rz']:
        if unit == 'µrad': print('For variation of {} µrad'.format(misalig))
        if unit == 'mrad': print('For variation of {} mrad'.format(misalig))
    
    if deg_f in ['Tx', 'Ty', 'Tz']:
        if unit == 'mm': print('For variation of {} mm'.format(round(misalig,3))) # round the number
        if unit == 'µm': print('For variation of {} µm'.format(misalig))
        
    print("####################################################################") 
    
    beam_copy = beam.duplicate()   # make a copy of the beam
    
    #
    # Run optical element 4 - Mirror
    #
    oe4.ALPHA = 90.0
    oe4.DUMMY = 0.1
    oe4.FHIT_C = 1 
    oe4.FILE_REFL = b'/home/ABTLUS/bruno.souza/Oasys/Rh.dat'                                          # CNPEM
    oe4.FILE_RIP = b'/home/ABTLUS/bruno.souza/GITHUB/SAPUCAIA/SPU_total_deformation_300mm_sh.dat'     # CNPEM
    # oe4.FILE_REFL = b'/home/bruno/Oasys/Rh.dat'                                          # HOME
    # oe4.FILE_RIP = b'/home/bruno/GITHUB/SAPUCAIA/SPU_total_deformation_300mm_sh.dat'     # HOME
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
    oe4.THETA = 89.7994647717
    oe4.T_IMAGE = 0.0
    oe4.T_INCIDENCE =  89.7994647717
    oe4.T_REFLECTION = 89.7994647717
    oe4.T_SOURCE = 0.0

    
    #%% Assigning misalignment values:
    
   # Rotation angle mapping to variables
    rotation_mapping = {
    "Rx": "Z_ROT",
    "Ry": "X_ROT",
    "Rz": "Y_ROT"      }
     
    rotation_factors = {
    ('Rx', 'µrad'): -1 * np.rad2deg(misalig * 1e-6),
    ('Ry', 'µrad'): -1 * np.rad2deg(misalig * 1e-6),
    ('Rz', 'µrad'):  1 * np.rad2deg(misalig * 1e-6),                                    
    ('Rx', 'mrad'): -1 * np.rad2deg(misalig * 1e-3),
    ('Ry', 'mrad'): -1 * np.rad2deg(misalig * 1e-3),
    ('Rz', 'mrad'):  1 * np.rad2deg(misalig * 1e-3)
    }
       
    
    # Translation mapping to variables
    translation_mapping = {
    "Tx": "OFFZ",
    "Ty": "OFFX",
    "Tz": "OFFY"          }
    
    translation_factors = {
    ('Tx', 'mm'): -1 * misalig,
    ('Ty', 'mm'): -1 * misalig,
    ('Tz', 'mm'):  1 * misalig,
    ('Tx', 'µm'): -1 * misalig * (1e-3),
    ('Ty', 'µm'): -1 * round(misalig, 3) * 1e-3,
    ('Tz', 'µm'):  1 * misalig * 1e-3             }
       
    key = (deg_f,unit)
    
    if deg_f in rotation_mapping:
        rotation_variable = rotation_mapping[deg_f]
        value = rotation_factors[key]
        setattr(oe4, rotation_variable, value)
       
    if deg_f in translation_mapping:
        translation_variable = translation_mapping[deg_f]
        value = translation_factors[key]
        setattr(oe4, translation_variable, value)
  
    print("####################################################################") 
    print(oe4.Z_ROT)
    print("####################################################################") 

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
 
    
    # Choose the visualization element based on the position              
    device = position_to_device.get(position, 'Unknown Device')
    
        
    beam_copy.retrace(position)                        # Propagating the beam to the position
   
 
    #%% Read Shadow beam
    beam2D = read_shadow_beam(beam_copy, x_column_index=3, y_column_index=1, nbins_x=nbins_x, nbins_y=nbins_y, nolost=1, ref=23, zeroPadding=2, gaussian_filter=0)
         
    #%% Adjustable the scale and plotting
    
    #%% Rotations
    if deg_f == 'Rx': 
        plot_range_x, plot_range_y = 400, 400      # the total range of the plot in [µm]    
        zero_pad_x, zero_pad_y = 4, 4              # the zeros      
        
    if deg_f == 'Ry': 
        plot_range_x, plot_range_y = 7000, 7000    # the total range of the plot in [µm] 
        zero_pad_x, zero_pad_y = 12, 6              # the zeros      
        
    if deg_f == 'Rz': 
        plot_range_x, plot_range_y = 10400, 10400  #  the total range of the plot in [µm]   
        zero_pad_x, zero_pad_y = 25, 25            # the zeros
        
    #%% Translations    
    if deg_f == 'Tx': 
        plot_range_x, plot_range_y = 450, 450      # the total range of the plot in [µm]    
        zero_pad_x, zero_pad_y = 2, 2              # the zeros      
        
    if deg_f == 'Ty': 
        plot_range_x, plot_range_y = 9000, 9000    # the total range of the plot in [µm] 
        zero_pad_x, zero_pad_y = 10, 40              # the zeros      
        
    if deg_f == 'Tz': 
        plot_range_x, plot_range_y = 460, 460  # the total range of the plot in [µm]   
        zero_pad_x, zero_pad_y = 5, 5            # the zeros
        
    #%% Filename
    
    # Rotations
    if deg_f in ['Rx','Ry','Rz']: 
        if unit =='µrad': filename = 'SPU '+ device + ' ' + deg_f +'=%.0fµrad'%(misalig)
        if unit =='mrad': filename = 'SPU '+ device + ' ' + deg_f +'=%.0fmrad'%(misalig)
        
    # Translations
    if deg_f in ['Tx','Ty','Tz']:
        if unit =='mm': filename = "SPU "+ device + ' ' + deg_f +"=%.1fmm"%(misalig)
        if unit =='µm': filename = 'SPU '+ device + ' ' + deg_f +'=%.0fµm'%(misalig)
        
    #%% Plot Beam
    outputs = plot_beam(beam2D, show_plot=True,outfilename=filename+'.png',outfileext='png',cut=0,textA=1,textB=5,textC=2,fitType=3,cmap='viridis',plot_title=filename,
                        zero_pad_x=zero_pad_x, zero_pad_y=zero_pad_y,
                        x_range = plot_range_x, y_range = plot_range_y,
                        x_range_min= -(plot_range_x/2)/1000,
                        x_range_max=  (plot_range_x/2)/1000,
                        y_range_min= -(plot_range_y/2)/1000,
                        y_range_max=  (plot_range_x/2)/1000)
   
    #%% Store the outputs    
 
    data['misalignment'].append(misalig)
    data['fwhm_h_sli'].append(outputs['fwhm_x'])
    data['fwhm_v_sli'].append(outputs['fwhm_z'])
    data['mean_pos_h'].append(outputs['mean_x'])
    data['mean_pos_v'].append(outputs['mean_z'])
    data['fwhm_h_fit'].append(outputs['fit_fwhm_x'])
    data['fwhm_v_fit'].append(outputs['fit_fwhm_z'])
    data['mean_pos_h_fit'].append(outputs['fit_mean_x'])
    data['mean_pos_v_fit'].append(outputs['fit_mean_z'])

    
#%% Save the data

# Titles for the columns
header_list = [deg_f + '_Misalignment' + '[' + unit + ']',
               "X FWHM [µm]",
               "Y FWHM [µm]",
               "X mean position",
               "Y mean position",
               "X FWHM fit [µm]",
               "Y FWHM fit [µm]",
               "X mean fit",
               "Y mean fit"]

# Create a DataFrame with your data
data2 = {
    header_list[0]: data['misalignment'],
    header_list[1]: data['fwhm_h_sli'],
    header_list[2]: data['fwhm_v_sli'],
    header_list[3]: data['mean_pos_h'],
    header_list[4]: data['mean_pos_v'],
    header_list[5]: data['fwhm_h_fit'],
    header_list[6]: data['fwhm_v_fit'],
    header_list[7]: data['mean_pos_h_fit'],
    header_list[8]: data['mean_pos_v_fit']
}

df = pd.DataFrame(data2)

# Specify the file name
file_name = filename

# Custom formatting function to center-align data
def format_center_align_data(value):
    return f'{value:^20}'  # You can adjust the width as needed

# Apply the formatting function to each data cell
df = df.applymap(format_center_align_data)

# Save the DataFrame to a .txt file with a tab separator, including header
df.to_csv(file_name, sep='\t', index=False)

# Reopen the file for modification
with open(file_name, 'r') as f:
    content = f.readlines()

# Modify the header line to center-align titles
content[0] = '\t'.join([f'{title:^20}' for title in header_list]) + '\n'

# Write the modified content back to the file
with open(file_name, 'w') as f:
    f.writelines(content)

print(f'Data saved to {file_name}')