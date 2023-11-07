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
from plot_beam_modified import plot_beam                                     # Plot the beam information
from optlnls.shadow import run_shadow_caustic, read_caustic,get_good_ranges  # Custom functions related to Shadow software
import matplotlib.pyplot as plt                                              # Library for creating charts and plots 
import pandas as pd                                                          # Library for data manipulation and analysis

#%% PARAMETERS:
n_rays= 1000000                            # number of rays used in the simulation
deg_f='Rx'                                 # the degree of freedom in Sirius frame;
unit ='µrad'                               # µrad, mrad for rotations; [µm], [mm] for translations
misalig_array=np.linspace(-450,450,3)      # misalignment array;

position = 13750                           # Position of interest

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
fwhm_h_fit_list = []; fwhm_v_fit_list = []; 
fwhm_h_sli_list = []; fwhm_v_sli_list = [];

z_fwhm_min_h_list=[]; z_fwhm_min_v_list=[];
z_rms_min_h_list= []; z_rms_min_v_list= [];

misalignment_list=[];

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
oe0.FILE_BOUND = b'/home/bruno/GITHUB/SAPUCAIA/SPU_optimize_source_66x66urad2.txt'
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
oe2.FILE_REFL = b'/home/bruno/Oasys/Si111.dat'
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
oe3.FILE_REFL = b'/home/bruno/Oasys/Si111.dat'
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
 
 
#
#run optical element 1
#
print("    Running optical element: %d"%(1))
if iwrite:
    oe1.write("start.01")
 
beam.traceOE(oe1,1)
 
if iwrite:
    oe1.write("end.01")
    beam.write("star.01")
 
 
#
#run optical element 2
#
print("    Running optical element: %d"%(2))
if iwrite:
    oe2.write("start.02")
 
beam.traceOE(oe2,2)
 
if iwrite:
    oe2.write("end.02")
    beam.write("star.02")
     
 
#
#run optical element 3 
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
        if unit == 'mm': print('For variation of {} mm'.format(misalig))
        if unit == 'µm': print('For variation of {} µm'.format(misalig))
        
    print("####################################################################") 
    
    beam_copy = beam.duplicate()   # make a copy of the beam 01
     
    #
    #run optical element 4 - Mirror
    #
    oe4.ALPHA = 90.0
    oe4.DUMMY = 0.1
    oe4.FHIT_C = 1
    oe4.FILE_REFL = b'/home/bruno/Oasys/Rh.dat'
    oe4.FILE_RIP = b'/home/bruno/GITHUB/SAPUCAIA/SPU_total_deformation_300mm_sh.dat'
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
    oe4.T_INCIDENCE = 89.7994647717
    oe4.T_REFLECTION = 89.7994647717
    oe4.T_SOURCE = 0.0
    
    # Assigning misalignment values:
    
    # Rotations [µrad]
    if(deg_f == 'Rx') and unit =='µrad': oe4.Z_ROT = (  -1*np.rad2deg(misalig*(1e-6))  )
    if(deg_f == 'Ry') and unit =='µrad': oe4.X_ROT = (  -1*np.rad2deg(misalig*(1e-6))  )
    if(deg_f == 'Rz') and unit =='µrad': oe4.Y_ROT = (  -1*np.rad2deg(misalig*(1e-6))  )
    # Rotations [mrad]
    if(deg_f == 'Rx') and unit =='mrad': oe4.Z_ROT = (  -1*np.rad2deg(misalig)*(1e-3)  )
    if(deg_f == 'Ry') and unit =='mrad': oe4.X_ROT = (  -1*np.rad2deg(misalig)*(1e-3)  )
    if(deg_f == 'Rz') and unit =='mrad': oe4.Y_ROT = (  -1*np.rad2deg(misalig)*(1e-3)  )
  
    # Translations [mm]   
    if(deg_f == 'Tx') and unit =='mm': oe4.OFFZ = -1*misalig
    if(deg_f == 'Ty') and unit =='mm': oe4.OFFX = -1*misalig
    if(deg_f == 'Tz') and unit =='mm': oe4.OFFY = misalig
    # Translations [µm]  
    if(deg_f == 'Tx') and unit =='µm': oe4.OFFZ = (-1*misalig)*(1e-3)
    if(deg_f == 'Ty') and unit =='µm': oe4.OFFX = (-1*misalig)*(1e-3)
    if(deg_f == 'Tz') and unit =='µm': oe4.OFFY = (misalig)*(1e-3)
    
    
    #
    #run optical element 4
    #
    print("    Running optical element: %d"%(4))
    if iwrite:
        oe4.write("start.04")
    
    
    beam_copy.traceOE(oe4,4)
    
    if iwrite:
        oe4.write("end.04")
        beam_copy.write("star.04")
 
    
    # Choosing the visualization element
    if(position == 841):   device='DVF3 '
    if(position == 13750): device='DVF4 '
    if(position == 16970): device='SMP '
    if(position == 22000): device='Hor.f '
    if(position == 25000): device='Ver.f '
    
    beamcaustic = beam_copy.duplicate()                # Beam after running the mirror
    beam_copy.retrace(position)                        # Beam propagate to a specific element
   
 
    #%% Read Shadow beam
    beam2D = read_shadow_beam(beam_copy, x_column_index=3, y_column_index=1, nbins_x=nbins_x, nbins_y=nbins_y, nolost=1, ref=23, zeroPadding=2, gaussian_filter=0)
        
    
    #%% Adjustable the scale and plotting
    
    # Rotations
    if deg_f == 'Rx': 
        plot_range_x, plot_range_y = 600, 600      # the total range of the plot in [µm]    
        zero_pad_x, zero_pad_y = 2, 2              # the zeros      
        
    if deg_f == 'Ry': 
        plot_range_x, plot_range_y = 6000, 6000    # the total range of the plot in [µm] 
        zero_pad_x, zero_pad_y = 8, 4              # the zeros      
        
    if deg_f == 'Rz': 
        plot_range_x, plot_range_y = 17000, 17000  #  the total range of the plot in [µm]   
        zero_pad_x, zero_pad_y = 14, 14            # the zeros
        
    # Translations    
    if deg_f == 'Tx': 
        plot_range_x, plot_range_y = 800, 800      # the total range of the plot in [µm]    
        zero_pad_x, zero_pad_y = 2, 2              # the zeros      
        
    if deg_f == 'Ty': 
        plot_range_x, plot_range_y = 4600, 4600    # the total range of the plot in [µm] 
        zero_pad_x, zero_pad_y = 6, 6              # the zeros      
        
    if deg_f == 'Tz': 
        plot_range_x, plot_range_y = 800, 800  # the total range of the plot in [µm]   
        zero_pad_x, zero_pad_y = 5, 5            # the zeros
        
    #%% Filename
    
    # Rotations
    if deg_f in ['Rx','Ry','Rz']: 
        if unit =='µrad': filename = 'SPU '+ device + deg_f +'=%.0fµrad'%(misalig)
        if unit =='mrad': filename = 'SPU '+ device + deg_f +'=%.0fmrad'%(misalig)
        
    # Translations
    if deg_f in ['Tx','Ty','Tz']:
        if unit =='mm': filename = 'SPU '+ device + deg_f +'=%.0fmm'%(misalig)
        if unit =='µm': filename = 'SPU '+ device + deg_f +'=%.0fµm'%(misalig)
        
       #%% Plot Beam
    outputs = plot_beam(beam2D, show_plot=False,outfilename=filename,outfileext='png',cut=0,textA=1,textB=5,textC=2,fitType=3,cmap='viridis',plot_title=filename,
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


    z0 = 0           # starting z
    zf = 30000       # final z
    nz = 31          # number of points

    Caustic = "Caustic (z0="+str(z0)+", zf="+str(zf)+", nz="+str(nz)+')'
    
    # Get the good range for x and y
    goodRange = get_good_ranges(beam=beamcaustic, zStart=z0, 
                                zFin=zf, colh=3, colv=1)
    
    # Ru caustic
    run_shadow_caustic(filename=Caustic, beam=beamcaustic, zStart=z0, zFin=zf, nz=nz, zOffset=0, colh=3, colv=1, colref=23, nbinsh=200, nbinsv=200,
                        xrange=[-2,2], 
                        yrange=[-2,2])
                       
                       # xrange=goodRange[0:2], 
                       # yrange=goodRange[2:4])

    histo_h, histo_v, caustic_dict = read_caustic(Caustic,plot=False,plot2D=False,cmap='viridis',figprefix='')
    
    z_fwhm_min_h = caustic_dict['z_fwhm_min_h'] 
    z_fwhm_min_v = caustic_dict['z_fwhm_min_v']
    z_rms_min_h = caustic_dict['center_rms_h']
    z_rms_min_v = caustic_dict['center_rms_v']
    
   
    #%% Store the outputs    
 
    # FOCUS 
    z_rms_min_h,z_rms_min_v = caustic_dict['center_rms_h'], caustic_dict['center_rms_v']
    z_fwhm_min_h_list.append(z_fwhm_min_h); z_fwhm_min_v_list.append(z_fwhm_min_v);
    z_rms_min_h_list.append(z_rms_min_h); z_rms_min_v_list.append(z_rms_min_v);
    
    # FWHM
    fwhm_h_sli, fwhm_v_sli = outputs['fwhm_x'], outputs['fwhm_z']
    fwhm_h_fit, fwhm_v_fit = outputs['fit_fwhm_x'], outputs['fit_fwhm_z']
    
    # FWHM slice and fit
    fwhm_h_sli_list.append(fwhm_h_sli); fwhm_v_sli_list.append(fwhm_v_sli);
    fwhm_h_fit_list.append(fwhm_h_fit); fwhm_v_fit_list.append(fwhm_v_fit);  
    
    # Misalignments
    misalignment_list.append(misalig)
    
#%% Plot the focus sensibility
#%% Rotations

############################ Rx - Yall ############################
if(deg_f == 'Rx'): 
   
    # Create a figure with two subplots and increase the vertical gap between them
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'hspace': 0.4})

    # Plot the horizontal focus position on the first subplot
    ax1.plot(misalignment_list, z_rms_min_h_list, 'b', label='Horizontal Focus Position')
    ax1.set_ylabel('Horizontal focus position')

    # Plot the vertical focus position on the second subplot
    ax2.plot(misalignment_list, z_rms_min_v_list, 'r', label='Vertical Focus Position')
    ax2.set_xlabel(deg_f + ' [' + unit + '] - Yall')
    ax2.set_ylabel('Vertical focus position')

    # Enable minor ticks for both subplots
    ax1.minorticks_on()
    ax2.minorticks_on()

    # Show the grid for both subplots
    ax1.grid(True)
    ax2.grid(True)

    # Set a title for the entire figure
    plt.suptitle('Focal Sensitivity')

    plt.tight_layout()  # Adjust subplot spacing
    plt.show()
    
############################ Ry - Pitch ############################    
if(deg_f == 'Ry') and caustic:
    plt.plot(misalignment_list,z_rms_min_h_list, '.', color='black', label='data')
    plt.xlabel(deg_f+' ['+unit+'] - Pitch') 
    plt.grid(True)
    plt.ylabel('Horizontal focus position')
    plt.title('Focal Sensitivity')

############################ Rz - Roll ############################    
if(deg_f == 'Rz') and caustic:
    plt.plot(misalignment_list,z_rms_min_v_list,  '.', color='black', label='data')
    plt.xlabel(deg_f+' ['+unit+'] - Roll') 
    plt.grid(True)
    plt.ylabel('Vertical focus position')
    plt.title('Focal Sensitivity')

#%% Translations


############################### Tx ###############################
if(deg_f == 'Tx') and caustic: 
   
    # Create a figure with two subplots and increase the vertical gap between them
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'hspace': 0.4})

    # Plot the horizontal focus position on the first subplot
    ax1.plot(misalignment_list, z_rms_min_h_list, 'b', label='Horizontal Focus Position')
    ax1.set_ylabel('Horizontal focus position')

    # Plot the vertical focus position on the second subplot
    ax2.plot(misalignment_list, z_rms_min_v_list, 'r', label='Vertical Focus Position')
    ax2.set_xlabel(deg_f + ' [' + unit + '] - Yall')
    ax2.set_ylabel('Vertical focus position')

    # Enable minor ticks for both subplots
    ax1.minorticks_on()
    ax2.minorticks_on()

    # Show the grid for both subplots
    ax1.grid(True)
    ax2.grid(True)

    # Set a title for the entire figure
    plt.suptitle('Focal Sensitivity')

    plt.tight_layout()  # Adjust subplot spacing
    plt.show()


############################### Ty ###############################
if(deg_f == 'Ty') and caustic: 
   
   # Create a figure with two subplots and increase the vertical gap between them
   fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'hspace': 0.4})

   # Plot the horizontal focus position on the first subplot
   ax1.plot(misalignment_list, z_rms_min_h_list, 'b', label='Horizontal Focus Position')
   ax1.set_ylabel('Horizontal focus position')

   # Plot the vertical focus position on the second subplot
   ax2.plot(misalignment_list, z_rms_min_v_list, 'r', label='Vertical Focus Position')
   ax2.set_xlabel(deg_f + ' [' + unit + '] - Yall')
   ax2.set_ylabel('Vertical focus position')

   # Enable minor ticks for both subplots
   ax1.minorticks_on()
   ax2.minorticks_on()

   # Show the grid for both subplots
   ax1.grid(True)
   ax2.grid(True)

   # Set a title for the entire figure
   plt.suptitle('Focal Sensitivity')

   plt.tight_layout()  # Adjust subplot spacing
   plt.show()
    
############################### Tz ###############################
if(deg_f == 'Tz') and caustic: 
    
   # Create a figure with two subplots and increase the vertical gap between them
   fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'hspace': 0.4})

   # Plot the horizontal focus position on the first subplot
   ax1.plot(misalignment_list, z_rms_min_h_list, 'b', label='Horizontal Focus Position')
   ax1.set_ylabel('Horizontal focus position')

   # Plot the vertical focus position on the second subplot
   ax2.plot(misalignment_list, z_rms_min_v_list, 'r', label='Vertical Focus Position')
   ax2.set_xlabel(deg_f + ' [' + unit + '] - Yall')
   ax2.set_ylabel('Vertical focus position')

   # Enable minor ticks for both subplots
   ax1.minorticks_on()
   ax2.minorticks_on()

   # Show the grid for both subplots
   ax1.grid(True)
   ax2.grid(True)

   # Set a title for the entire figure
   plt.suptitle('Focal Sensitivity')

   plt.tight_layout()  # Adjust subplot spacing
   plt.show()
          


#%% Save the data

#%% Caustic

if caustic:
    # Original header list
    header_list = [deg_f+' Misalignment '+'['+unit+']', 
                   'z_fwhm_min_h',
                   'z_fwhm_min_v',
                   'z_rms_min_h []',
                   'z_rms_min_v []',]
    
    # Custom formatting function to center-align header
    def format_center_align_header(value):
        return f'{value:^20}'  # You can adjust the width as needed
    
    # Apply the formatting function to each header
    formatted_headers = [format_center_align_header(header) for header in header_list]
    
    # Create a DataFrame with the formatted headers
    data = pd.DataFrame({
        formatted_headers[0]: misalignment_list,
        formatted_headers[1]: z_fwhm_min_h_list,
        formatted_headers[2]: z_fwhm_min_v_list
        formatted_headers[3]: z_rms_min_h_list,  
        formatted_headers[4]: z_rms_min_v_list
    })
    
    # Specify the file name
    file_name = filename
    
    # Custom formatting function to center-align data
    def format_center_align_data(value):
        return f'{value:^20}'  # You can adjust the width as needed
    
    # Apply the formatting function to each data column
    for column in data.columns:
        data[column] = data[column].apply(format_center_align_data)
    
    # Save the DataFrame to a .txt file with a tab separator
    data.to_csv(file_name, sep='\t', index=False, header=False)
    
    # Save the header separately
    with open(file_name, 'r') as f:
        content = f.read()
    
    with open(file_name, 'w') as f:
        f.write('\t'.join(formatted_headers) + '\n')
        f.write(content)
    
    print(f'Data saved to {file_name}')
