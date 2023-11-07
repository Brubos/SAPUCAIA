#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 14:11:04 2023

@author: bruno.souza
"""

########### Run Sapucaia Beamline ###########

### Packages:
    
import numpy as np  # scientific computing library 
import Shadow       # used for shadow simulations


### Function to run Sapucaia:

def run_Sapucaia(n_rays=1000000, energy=12000, delta_E=10, hor_accept=66e-6, ver_accept=66e-6, sig_h=0.065347,
                 sig_v=0.003893, div_h=1.0112e-05, div_v=9.405e-06, d_image=16850, atomic_plane="Si111"):
    
    '''
    
   Run Sapucaia Beamline 
   
   Parameters:
       
       - n_rays: number of rays in the simulation (float);
       - energy: energy [eV] (float, array or list);
       - delta_E: energy variation [eV] (float); 
       - hor_accept: mirror horizontal acceptability [rad] (float); 
       - ver_accept: mirror vertical acceptability [rad] (float);
       - sig_h: sigma X [mm] (float);  
       - sig_v: sigma Z [mm] (float); 
       - div_h: horizontal sigma divergence [rad] (float);  
       - div_v: vertical sigma divergence [rad] (float);  
       - d_image: image plane distance from mirror [mm] (float); 
       - atomic_plane: filename used in the simulation (string); 
    
    Returns:
     
       - R: a file with the characteristics of the beam after passing through the optical elements of the simulation (Shadow extension);
       
     '''
     
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
    
    
    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0
    
    
    # initialize shadow3 source (oe0) and beam
    # oe - optical elements
    beam = Shadow.Beam()
    oe0  = Shadow.Source()
    oe1  = Shadow.OE()
    oe2  = Shadow.OE()
    oe3  = Shadow.OE()
    oe4  = Shadow.OE()
    
    
    # Define variables. See meaning of variables in: 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    
    
    oe0.FDISTR = 3
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = hor_accept/2 #3.3e-05
    oe0.HDIV2 = hor_accept/2 #3.3e-05
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 5676561
    oe0.NPOINT = n_rays #1000000
    oe0.PH1 = energy - delta_E/2 #11995.0
    oe0.PH2 = energy + delta_E/2 #12005.0
    oe0.SIGDIX = div_h        #1.0112e-05  
    oe0.SIGDIZ = div_v        #9.405e-06
    oe0.SIGMAX = sig_h        #0.065347
    oe0.SIGMAZ = sig_v        #0.003893
    oe0.VDIV1 = ver_accept/2 #3.3e-05
    oe0.VDIV2 = ver_accept/2 #3.3e-05
    
    oe1.DUMMY = 0.1
    oe1.FWRITE = 3
    oe1.F_REFRAC = 2
    oe1.F_SCREEN = 1
    oe1.I_SLIT = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe1.N_SCREEN = 1
    oe1.RX_SLIT = np.array([hor_accept*26000, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe1.RZ_SLIT = np.array([ver_accept*26000, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 0.0
    oe1.T_REFLECTION = 180.0
    oe1.T_SOURCE = 26000.0
    
    oe2.DUMMY = 0.1
    oe2.FILE_REFL = bytes('/home/ABTLUS/bruno.souza/Oasys/'+atomic_plane+'.dat', 'utf-8') #b'/home/ABTLUS/bruno.souza/Oasys/Si111.dat'
    oe2.FWRITE = 3
    oe2.F_CENTRAL = 1
    oe2.F_CRYSTAL = 1
    oe2.PHOT_CENT = energy #12000.0
    oe2.T_IMAGE = 50.0
    oe2.T_INCIDENCE = 0.0
    oe2.T_REFLECTION = 0.0
    oe2.T_SOURCE = 2950.0
    
    oe3.ALPHA = 180.0
    oe3.DUMMY = 0.1
    oe3.FILE_REFL = bytes('/home/ABTLUS/bruno.souza/Oasys/'+atomic_plane+'.dat', 'utf-8') #b'/home/ABTLUS/bruno.souza/Oasys/Si111.dat'
    oe3.FWRITE = 3
    oe3.F_CENTRAL = 1
    oe3.F_CRYSTAL = 1
    oe3.PHOT_CENT = energy #12000.0
    oe3.T_IMAGE = -50.0
    oe3.T_INCIDENCE = 0.0
    oe3.T_REFLECTION = 0.0
    oe3.T_SOURCE = 50.0
    
    oe4.ALPHA = 270.0
    oe4.DUMMY = 0.1
    oe4.FHIT_C = 1
    oe4.FILE_REFL = b'/home/ABTLUS/bruno.souza/Oasys/Rh.dat'
    oe4.FMIRR = 3
    oe4.FWRITE = 1
    oe4.F_DEFAULT = 0
    oe4.F_REFLEC = 1
    oe4.RLEN1 = 150.0
    oe4.RLEN2 = 150.0
    oe4.RWIDX1 = 2.5
    oe4.RWIDX2 = 2.5
    oe4.SIMAG = 25000.0
    oe4.SSOUR = 31000.0
    oe4.THETA = 89.7994647717
    oe4.T_IMAGE = d_image # 16850.0 # in relation to the mirror
    oe4.T_INCIDENCE = 89.7994647717
    oe4.T_REFLECTION = 89.7994647717
    oe4.T_SOURCE = 2000.0
    
    
    
    # Run SHADOW to create the source
    
    if iwrite:
        oe0.write("start.00")
    
    beam.genSource(oe0)
    
    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")
    
    
    
    #run optical element 1
    
    
    print("    Running optical element: %d"%(1))
    if iwrite:
        oe1.write("start.01")
    
    beam.traceOE(oe1,1)
    
    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")
    
    
    
    #run optical element 2
    
    print("    Running optical element: %d"%(2))
    if iwrite:
        oe2.write("start.02")
    
    beam.traceOE(oe2,2)
    
    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")
    
    

    #run optical element 3
    
    print("    Running optical element: %d"%(3))
    if iwrite:
        oe3.write("start.03")
    
    beam.traceOE(oe3,3)
    
    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")
       
    
 
    #run optical element 4
  
    print("    Running optical element: %d"%(4))
    if iwrite:
        oe4.write("start.04")
    
    beam.traceOE(oe4,4)
    
    if iwrite:
        oe4.write("end.04")
        beam.write("star.04")
    
    
    Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")
    Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    
    return beam
