# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 18:27:20 2023

@author: sofia
"""

import sys;
import struct
import numpy as np
import math
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
ifile = "InSitu/GPRdata/IS3_GPR0-5.txt" 
#when only 2 shots were taken, the metal plate was not used, meaning there is uss and uis.

f = open(ifile, 'r')
data = np.genfromtxt(f, dtype=float)

window_size = 100

depth=data[:,0]-3 # cut the surface to be at the first peak
depth=depth/10
#graph has been cut so that y=0 corresponds with the surface." (Sørens bachelor)
amp=data[:,1] #unit: dB

uss_depth=depth[:50000]
uss_amp=amp[:50000]
uis_depth=depth[50001:]
uis_amp=amp[50001:]


moving_averages_uss_depth = []
moving_averages_uss_amp = []
moving_averages_uis_depth = []
moving_averages_uis_amp = []

i=0 
  
# Loop through the array t o 
#consider every window of size 3 
while i < len(uss_depth) - window_size + 1: 
  
    # Calculate the average of current window 
    window_average = round(np.sum(uss_depth[ 
      i:i+window_size]) / window_size, 2) 
      
    # Store the average of current 
    # window in moving average list 
    moving_averages_uss_depth.append(window_average) 
      
    # Shift window to right by one position 
    i += 1

i=0
while i < len(uss_amp) - window_size + 1: 
  
    # Calculate the average of current window 
    window_average = round(np.sum(uss_amp[ 
      i:i+window_size]) / window_size, 2) 
      
    # Store the average of current 
    # window in moving average list 
    moving_averages_uss_amp.append(window_average) 
      
    # Shift window to right by one position 
    i += 1
  
    
i=0 
  
# Loop through the array t o 
#consider every window of size 3 
while i < len(uis_depth) - window_size + 1: 
  
    # Calculate the average of current window 
    window_average = round(np.sum(uis_depth[ 
      i:i+window_size]) / window_size, 2) 
      
    # Store the average of current 
    # window in moving average list 
    moving_averages_uis_depth.append(window_average) 
      
    # Shift window to right by one position 
    i += 1

i=0
while i < len(uis_amp) - window_size + 1: 
  
    # Calculate the average of current window 
    window_average = round(np.sum(uis_amp[ 
      i:i+window_size]) / window_size, 2) 
      
    # Store the average of current 
    # window in moving average list 
    moving_averages_uis_amp.append(window_average) 
      
    # Shift window to right by one position 
    i += 1
    

# undisturbed_snow_surface = uss
# metal_plate_on_snow_surface = mpss
# metal plate on ice surface = mpis
# undisturbed ice surface = uis. Here a pit was dug of at least 1x1 m.
# metal plate on ice surface (no snow) = mpisns


# #Inserting the specific profile IS3-GPR-0-2
# fine=0.013
# depth_hoare=0.013+0.184
# fine_01=0.013+0.208
# depth_hoare_01 = 0.013 + 0.292
# fine_02=0.013+0.32
# ice_lense=0.013+0.41
# fine_03=0.013+0.411
# ice_lense_01=0.0130+0.435
# depth_hoare_02=0.013+0.436
# ice_lense_02=0.013+0.46
# depth_hoare_03=0.013+0.461
# medium=0.013+0.50
# depth_hoare_04=0.013+0.92
# fine_04=0.013+0.94
# depth_hoare_05=0.013+0.98
# fine_05=0.013+1.01
# depth_hoare_06=0.013+1.24
# coarse=0.013+1.28
# depth_hoare_07=0.013+1.65
# coarse_01=0.013+1.69
# ice_lense_03=0.013+2.146
# coarse_02=0.013+2.147
# ice_lense_04=0.013+2.67
# coarse_03=0.013+2.675


# plt.figure()
# plt.plot(moving_averages_uss_amp, moving_averages_uss_depth, label='undisturbed snow surface', color='darkmagenta')
# #plt.plot(moving_averages_uis_amp, moving_averages_uis_depth, label='undisturbed ice surface', color='cadetblue')
# plt.axhline(fine, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare, color = 'k', linewidth = 1)
# plt.axhline(fine_01, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare_01, color = 'k', linewidth = 1)
# plt.axhline(fine_02, color = 'k', linewidth = 1)
# plt.axhline(ice_lense, color = 'indianred', alpha=0.6, linewidth = 2)
# plt.axhline(fine_03+0.05,  color = 'indianred', alpha=0.6, linewidth = 2)
# plt.axhline(ice_lense_01, color = 'indianred', alpha=0.6, linewidth = 2)
# plt.axhline(depth_hoare_02, color = 'k', linewidth = 1)
# plt.axhline(ice_lense_02, color = 'indianred', alpha=0.6, linewidth = 2)
# plt.axhline(depth_hoare_03,  color = 'k', linewidth = 1)
# plt.axhline(medium, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare_04, color = 'k', linewidth = 1)
# plt.axhline(fine_04, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare_05, color = 'k', linewidth = 1)
# plt.axhline(fine_05, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare_06, color = 'k', linewidth = 1)
# plt.axhline(coarse, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare_07, color = 'k', linewidth = 1)
# plt.axhline(coarse_01, color = 'k', linewidth = 1)
# plt.axhline(ice_lense_03, color = 'indianred', alpha=0.6, linewidth = 2)
# plt.axhline(coarse_02, color = 'k', linewidth = 1)
# plt.axhline(ice_lense_03, color = 'indianred', alpha=0.6, linewidth = 2)
# plt.axhline(coarse_03, color = 'k', linewidth = 1)


# plt.text(-79, fine, 'fine', size=5)
# plt.text(-79, depth_hoare, 'depth hoare', size=5)
# plt.text(-79, fine_01+0.03, 'fine', size=5)
# plt.text(-79, depth_hoare_01, 'depth hoare, larger crystals', size=5)
# plt.text(-79, fine_02+0.03, 'fine', size=5)
# plt.text(-79, ice_lense, 'ice lense', size=5, color = 'indianred')
# plt.text(-79+6, fine_03, 'fine', size=5)
# plt.text(-79+9, ice_lense_01, 'ice lense', size=5, color = 'indianred')
# plt.text(-79+15, depth_hoare_02, 'depth hoare', size=5)
# plt.text(-79+24, ice_lense_02, 'ice lense', size=5,color = 'indianred')
# plt.text(-79, depth_hoare_03+0.02, 'depth hoare', size=5)
# plt.text(-79, medium+0.04, 'medium', size=5)
# plt.text(-79, depth_hoare_04-0.03, 'depth hoare, airy gaps', size=5)
# plt.text(-79, fine_04, 'fine, firm powder', size=5)
# plt.text(-79+9, depth_hoare_05, 'depth hoare, airy gaps', size=5)
# plt.text(-79, fine_05+0.02, 'fine, firm powder', size=5)
# plt.text(-79, depth_hoare_06-0.01, 'depth hoare, crystal air gap', size=5)
# plt.text(-79, coarse, 'coarse, firm, granular ice crystals', size=5)
# plt.text(-79, depth_hoare_07, 'depth hoare, crystal air gap', size=5)
# plt.text(-79, coarse_01+0.03, 'coarse, firm', size=5)
# plt.text(-79, ice_lense_03, 'ice lense', size=5,color = 'indianred')
# plt.text(-79, coarse_02+0.05, 'coarse, firm', size=5)
# plt.text(-79, ice_lense_03, 'ice lense, very solid', size=5,color = 'indianred')
# plt.text(-79, coarse_03, 'coarse, firm', size=5)


# plt.legend(loc='best')
# plt.ylim([2.83,0])
# plt.xlabel('amplitude [dB]')
# plt.ylabel('depth (m)')
# plt.savefig('Figures/InSitu/IS3_TRANS_GPR0-2', bbox_inches='tight', dpi=200)
# plt.show()





# #Inserting the specific profile IS3-GPR-150-2
# fine=0.013
# ice_lense=0.013+0.30
# fine_01=0.013+0.301
# fine_02 = 0.013 + 0.37
# depth_hoare=0.013+0.60
# fine_03=0.013+0.62
# depth_hoare_01=0.013+1.0
# fine_04=0.0130+1.02
# depth_hoare_02=0.013+1.14
# coarse=0.013+1.17
# ice_lense_02=0.013+1.45
# coarse_01=0.013+1.46
# depth_hoare_03=0.013+2.20
# ice_lense_03=0.013+2.50
# depth_hoare_04=0.013+2.51

# plt.figure()
# plt.plot(moving_averages_uss_amp, moving_averages_uss_depth, label='undisturbed snow surface', color='darkmagenta')
# #plt.plot(moving_averages_uis_amp, moving_averages_uis_depth, label='undisturbed ice surface', color='cadetblue')
# plt.axhline(fine, color = 'k', linewidth = 1)
# plt.axhline(ice_lense, color = 'indianred', alpha=0.6, linewidth = 2)
# plt.axhline(fine_01, color = 'k', linewidth = 1)
# plt.axhline(fine_02, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare, color = 'k', linewidth = 1)
# plt.axhline(fine_03, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare_01, color = 'k', linewidth = 1)
# plt.axhline(fine_04, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare_02, color = 'k', linewidth = 1)
# plt.axhline(coarse, color = 'k', linewidth = 1)
# plt.axhline(ice_lense_02, color = 'indianred', alpha=0.6, linewidth = 2)
# plt.axhline(coarse_01, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare_03, color = 'k', linewidth = 1)
# plt.axhline(ice_lense_03, color = 'indianred', alpha=0.6, linewidth = 2)
# plt.axhline(depth_hoare_04, color = 'k', linewidth = 1)

# plt.text(-79, fine, 'fine', size=8)
# plt.text(-79, ice_lense-0.15, 'ice lense', color = 'indianred', size=8)
# plt.text(-79, fine_01, 'fine, lighter crystals', size=8)
# plt.text(-79, fine_02, 'fine, firm snow', size=8)
# plt.text(-79, depth_hoare, 'depth hoare', size=8)
# plt.text(-79, fine_03+0.05, 'fine', size=8)
# plt.text(-79, depth_hoare_01, 'depth hoare', size=8)
# plt.text(-79, fine_04 + 0.06, 'fine', size=8)
# plt.text(-79, depth_hoare_02 + 0.01, 'depth hoare', size=8)
# plt.text(-79, coarse + 0.05, 'coarse', size=8)
# plt.text(-79, ice_lense_02, 'ice lense', color = 'indianred', size=8)
# plt.text(-79, coarse_01+0.07, 'coarse', size=8)
# plt.text(-79, depth_hoare_03, 'depth hoare', size=8)
# plt.text(-79, ice_lense_03, 'ice lense', color = 'indianred', size=8)
# plt.text(-79, depth_hoare_04+0.08, 'depth hoare', size=8)


# plt.legend(loc='best')
# plt.ylim([2.82,0])
# plt.xlabel('amplitude [dB]')
# plt.ylabel('depth (m)')
# plt.savefig('Figures/InSitu/IS3_TRANS_GPR150-2', bbox_inches='tight', dpi=200)
# plt.show()






# #Inserting the specific profile IS3-GPR-300-2
# fine=0.013
# medium=0.013+0.13
# fine_01=0.013+0.15
# coarse = 0.013 + 0.285
# coarse_01=0.013+0.38
# depth_hoare=0.013+0.48
# medium_01=0.013+0.505
# coarse_02=0.013+0.75
# depth_hoare_02=0.013+0.95
# slush=0.0130+1.45
# plt.figure()
# plt.plot(moving_averages_uss_amp, moving_averages_uss_depth, label='undisturbed snow surface', color='darkmagenta')
# #plt.plot(moving_averages_uis_amp, moving_averages_uis_depth, label='undisturbed ice surface', color='cadetblue')
# plt.axhline(fine, color = 'k', linewidth = 1)
# plt.axhline(medium, color = 'k', linewidth = 1)
# plt.axhline(fine_01, color = 'k', linewidth = 1)
# plt.axhline(coarse, color = 'k', linewidth = 1)
# plt.axhline(coarse_01, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare, color = 'k', linewidth = 1)
# plt.axhline(medium_01, color = 'k', linewidth = 1)
# plt.axhline(coarse_02, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare_02, color = 'k', linewidth = 1)
# plt.axhline(slush, color = 'k', linewidth = 1)
# plt.text(-80, fine+0.03, 'fine', size=8)
# plt.text(-80, medium, 'medium', size=8)
# plt.text(-80, fine_01+0.03, 'fine', size=8)
# plt.text(-80, coarse, 'coarse', size=8)
# plt.text(-80, coarse_01, 'coarse, more compact', size=8)
# plt.text(-80, depth_hoare-0.01, 'depth hoare', size=8)
# plt.text(-80, medium_01+0.01, 'medium', size=8)
# plt.text(-80, coarse_02, 'coarse',size=8)
# plt.text(-80, depth_hoare_02, 'depth hoare', size=8)
# plt.text(-80, slush + 0.01, 'slush', size=8)
# plt.legend(loc='best')
# plt.ylim([1.45,0])
# plt.xlabel('amplitude [dB]')
# plt.ylabel('depth (m)')
# plt.savefig('Figures/InSitu/IS3_TRANS_GPR300-2', bbox_inches='tight', dpi=200)
# plt.show()



# #Inserting the specific profile GPR-1-343-5
# fresh_snow=0.013
# fine_firm_snow=0.013+0.02
# hoare=0.013+0.25
# medium = 0.013 + 0.30
# depth_hoare=0.013+0.38
# coarse=0.013+0.42
# very_coarse=0.013+0.745
# depth_hoare_01=0.013+0.87
# depth_hoare_02=0.013+1.10
# slush=0.0130+1.20
# plt.figure()
# plt.plot(moving_averages_uss_amp, moving_averages_uss_depth, label='undisturbed snow surface', color='darkmagenta')
# #plt.plot(moving_averages_uis_amp, moving_averages_uis_depth, label='undisturbed ice surface', color='cadetblue')
# plt.axhline(fresh_snow, color = 'k', linewidth = 1)
# plt.axhline(fine_firm_snow, color = 'k', linewidth = 1)
# plt.axhline(hoare, color = 'k', linewidth = 1)
# plt.axhline(medium, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare, color = 'k', linewidth = 1)
# plt.axhline(coarse, color = 'k', linewidth = 1)
# plt.axhline(very_coarse, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare_01, color = 'k', linewidth = 1)
# plt.axhline(depth_hoare_02, color = 'k', linewidth = 1)
# plt.axhline(slush, color = 'k', linewidth = 1)
# plt.text(-80, fresh_snow-0.004, 'fresh snow', size=8)
# plt.text(-80, fine_firm_snow+0.018, 'fine, firm snow', size=8)
# plt.text(-80, hoare, 'depth hoare', size=8)
# plt.text(-80, medium, 'firm, larger crystals', size=8)
# plt.text(-80, depth_hoare, 'firm, coarser granules', size=8)
# plt.text(-80, coarse, 'coarse', size=8)
# plt.text(-80, very_coarse, 'very coarse', size=8)
# plt.text(-80, depth_hoare_01, 'depth hoare',size=8)
# plt.text(-80, depth_hoare_02 - 0.01, 'depth hoare, larger crystals', size=8)
# plt.text(-80, slush + 0.01, 'slush', size=8)
# plt.legend(loc='best')
# plt.ylim([1.23,0])
# plt.xlabel('amplitude [dB]')
# plt.ylabel('depth (m)')
# plt.savefig('Figures/InSitu/IS3_TRANS_GPR343-5', bbox_inches='tight', dpi=200)
# plt.show()

#Inserting the specific profile IS3_GPR0-5

fresh_snow=0.013
medium=0.013+0.03
coarse=0.013+0.255
coarsening=0.013+0.325
coarse_01 = 0.013 + 0.79
very_coarse=0.013+0.85
ice=0.013+1.01

plt.figure()
plt.plot(moving_averages_uss_amp, moving_averages_uss_depth, label='undisturbed snow surface', color='darkmagenta')
#plt.plot(moving_averages_uis_amp, moving_averages_uis_depth, label='undisturbed ice surface', color='cadetblue')
plt.axhline(fresh_snow, color = 'k', linewidth = 1)
plt.axhline(medium, color = 'k', linewidth = 1)
plt.axhline(coarse, color = 'k', linewidth = 1)
plt.axhline(coarsening, color = 'k', linewidth = 1)
plt.axhline(coarse_01, color = 'k', linewidth = 1)
plt.axhline(very_coarse, color = 'k', linewidth = 1)
plt.axhline(ice, color = 'k', linewidth = 1)

plt.text(-80, fresh_snow, 'fresh snow', size=8)
plt.text(-80, medium, 'medium', size=8)
plt.text(-80, coarse, 'coarse', size=8)
plt.text(-80, coarsening, 'coarsening', size=8)
plt.text(-80, coarse_01, 'coarse', size=8)
plt.text(-80, very_coarse, 'very coarse', size=8)
plt.text(-80, ice, 'ice', size=8)

plt.legend(loc='best')
plt.ylim([1.02,0])
plt.xlabel('amplitude [dB]')
plt.ylabel('depth (m)')
plt.savefig('Figures/InSitu/IS3_GPR0-5', bbox_inches='tight', dpi=200)
plt.show()



