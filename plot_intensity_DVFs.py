
from optlnls.importing import read_srw_int
from optlnls.plot import plot_beam
import matplotlib.pyplot as plt
import glob

# PLOT MULTI-E:

if(1): # cut at (0,0)

    filelist = glob.glob('*_mE.dat')
    filelist.sort(key = lambda x: x.split('_DVF')[-1][0])
    
    ur = 1
    
    cmap= 'viridis' # 'jet'# 'plasma' # 'viridis'
    
    title_list = ['PAPU: 8 keV @ DVF 2      Ideal Field',
                  'PAPU: 8 keV @ DVF 2      Measured Field',
                  'PAPU: 8 keV @ DVF 3      Ideal Field',
                  'PAPU: 8 keV @ DVF 3      Measured Field',
                  'PAPU: 8 keV @ DVF 4      Ideal Field',
                  'PAPU: 8 keV @ DVF 4      Measured Field']
    
    rw_list = [0.8, 0.8, 0.8, 0.8, 0.4, 0.4]
    
    for i in range(len(filelist)):
    
        filename = filelist[i]
        
        title = title_list[i]   
        
        rw = rw_list[i]
        
        beam = read_srw_int(filename)
        beam = beam[0]
        
        beam_corrected = beam[1:,1:].copy()
        beam[1:,1:] = beam_corrected[-1::-1]
        
        prefix = filename[:-4]
        
        pic_filename = prefix + '_cut.png'
        
        plot_beam(beam2D=beam, outfilename=pic_filename, cut=0, textA=1, textB=5, textC=2, textD=13, x_range=ur, y_range=ur,
    			  cmap=cmap, x_range_min=-rw, x_range_max=rw, y_range_min=-rw, y_range_max=rw,
    			  fitType=0, plot_title=title, zero_pad_x=1, zero_pad_y=1, export_slices=False)    


if(0): # cut at (0,0)
    
    filename = 'SPU_PAPU_8keV_DVF2_Measured_Field_mE.dat' # 'SPU_PAPU_8keV_48m_Ideal_Field_mE.dat' # 'SPU_PAPU_8keV_48m_Measured_Field_mE.dat'
    rw = 1.6
    
    title = 'PAPU: 8 keV @ DVF 2'        
    
    beam = read_srw_int(filename)
    beam = beam[0]
    
    beam_corrected = beam[1:,1:].copy()
    beam[1:,1:] = beam_corrected[-1::-1]
    
    prefix = filename[:-4]
    
    pic_filename = prefix + '_cut.png'
    
    cmap= 'jet' # 'jet'# 'plasma' # 'viridis'
    
    ur=1
    
    plot_beam(beam2D=beam, outfilename=pic_filename, cut=2, textA=1, textB=5, textC=6, textD=13, x_range=ur, y_range=ur,
			  cmap=cmap, x_range_min=-rw, x_range_max=rw, y_range_min=-rw, y_range_max=rw,
			  fitType=0, plot_title=title, zero_pad_x=1, zero_pad_y=1, export_slices=False)
    
    
if(0): # integrated
    
    filename = 'SPU_PAPU_8keV_DVF2_Measured_Field_mE.dat' # 'SPU_PAPU_8keV_48m_Ideal_Field_mE.dat'
    rw = 1.6
    
    title = 'PAPU: 8 keV @ DVF 2'        
    
    beam = read_srw_int(filename)
    beam = beam[0]
    
    beam_corrected = beam[1:,1:].copy()
    beam[1:,1:] = beam_corrected[-1::-1]
    
    prefix = filename[:-4]
    
    pic_filename = prefix + '_integrated.png'
    
    cmap= 'jet' # 'jet'# 'plasma' # 'viridis'
    
    ur=1
    
    plot_beam(beam2D=beam, outfilename=pic_filename, cut=0, textA=1, textB=5, textC=6, textD=13, x_range=ur, y_range=ur,
			  cmap=cmap, x_range_min=-rw, x_range_max=rw, y_range_min=-rw, y_range_max=rw,
			  fitType=0, plot_title=title, zero_pad_x=1, zero_pad_y=1, export_slices=False)


# Teste:

if(0):
    
    # beam_corrected = beam[1:,1:].copy()
    
    # beam[1:,1:] = beam_corrected[-1::-1]
    
    plt.figure()
    plt.imshow(beam[1:,1:], extent=[beam[0,1:][0], beam[0,1:][-1], beam[1:,0][0], beam[1:,0][-1]], cmap='jet')
    



