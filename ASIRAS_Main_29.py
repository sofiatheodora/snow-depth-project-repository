#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: Sofia Thirslund
Uses importer from Henriette Skourup
re-tracks using tfmra provided by Alessandro Di Bella

"""

from ASIRAS_Importer_29 import *

TrmfaThreshold=0.5

from tfmra_py3 import tfmra_retracker
from ocog_retracker import OCOGretracker
import sys;
import struct
import numpy as np
import math
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
#import Echogram
plt.style.use('ggplot')
c_A=299702547 

# IS_3_Lat = -76.7767
# IS_3_Lon = -52.2400


#Plotting selected Lon Lat
plt.plot(Lon,Lat)
plt.title('flown route')
plt.savefig('Figures/122/LonLatRoute.png',dpi=170)
plt.show()


# #Picking data closest to the line
from scipy import spatial

r2=[-52.202,-76.7358]
r1=[-52.191,-76.7368]

coord_Lon=np.expand_dims(Lon, axis=0)
coord_Lat=np.expand_dims(Lat, axis=0)

coord=np.concatenate((coord_Lon.T,coord_Lat.T),axis=1)
print(np.shape(coord))

distance_r2,index_r2 = spatial.KDTree(coord).query(r2)
print(distance_r2)
print(index_r2)

distance_r1,index_r1 = spatial.KDTree(coord).query(r1)
print(distance_r1)
print(index_r1)

data_start=index_r2-5
data_end=index_r1+5
data_length=data_end-data_start

selected_Lon=Lon[data_start:data_end]
selected_Lat=Lat[data_start:data_end]
selected_Hwgs84=Hwgs84[data_start:data_end]

#Plotting selected Lon Lat
plt.plot(Lon,Lat, color='k')
plt.plot(selected_Lon,selected_Lat, color='r', linewidth=5)
plt.title('Data Selection')
plt.xlabel('Lon (\N{DEGREE SIGN})')
plt.ylabel('Lat (\N{DEGREE SIGN})')
plt.savefig('Figures/0011/DataSelectionOnFull.png',dpi=170)
plt.show()



#Plotting Hwgs as a colormap 
# Apply a fancy colormap to the figure
cmap = plt.get_cmap('plasma')
plt.set_cmap(cmap)
plt.scatter(selected_Lon, selected_Lat, c=selected_Hwgs84)
cbar=plt.colorbar()
cbar.set_label('Hwgs84 (m)')
plt.title('Selected bit of flown route with elevation')
plt.xlabel('Lon (\N{DEGREE SIGN})')
plt.ylabel('Lat (\N{DEGREE SIGN})')
plt.savefig('Figures/0011/selected Hwgs as a colormap',dpi=180)
plt.show()


#Plotting selected Lon Lat
plt.plot(selected_Lon,selected_Lat)
plt.title('Data Selection')
plt.savefig('Figures/0011/DataSelection.png',dpi=170)
plt.show()




print(data_length)

#Plotting selected data H_ocog (from ESA)
#Plotting Hocog as a function of time. Something is not right here.
selected_samples_var=np.linspace(0, int(data_length)-1, int(data_length))
plt.plot(selected_samples_var,Hocog[data_start:data_end],color='k')
plt.title('Hocog')
plt.xlabel('sample point')
plt.ylabel('Hocog (m)')
plt.title('ESA Hocog not roll sorted')
plt.savefig('Figures/0011/ESA Hocog not roll sorted.png', dpi=200)
plt.show()


g = range(1,np.size(PWF[0,:])+1)
print(g)

#Calculating H tfmra (original) for selection


selected_PWF=PWF[data_start:data_end]
selected_WindowDelay=WindowDelay[data_start:data_end]

H_trmfa_selected=np.ones(data_length)

gate_rtck_array=np.ones(data_length)
for x in range(data_length):
    gate_rtck, i_peak, AmpPmax, i_peakMax  = tfmra_retracker(selected_PWF[x,:],TrmfaThreshold, .8, 5)
    if gate_rtck != math.nan:
        H_trmfa_selected[x]=selected_Hwgs84[x]-((gate_rtck-BinRef)*BinWidth+(selected_WindowDelay[x]*c_A)/2)
    gate_rtck_array[x]=gate_rtck

#################################################################################3
# plot single waveforms     
#################################################################################3
plt.figure()
plt.plot(g, selected_PWF[round(99),:].transpose())
#plt.plot(g, gate_rtck)
plt.title('Waveform')
plt.savefig('Figures/122/ASIRAS_SCALED_'+str(Time[data_length])+'.png',dpi=200)        
#################################################################################3


#Plotting roll on selection
selected_Hocog=Hocog[data_start:data_end]
selected_Roll=Roll[data_start:data_end]

plt.figure()
plt.plot(selected_samples_var,selected_Roll,color='k')
plt.xlabel('samples')
plt.ylabel('Roll angle (\N{DEGREE SIGN})')
plt.savefig('Figures/0011/RollonSelection.png')
plt.show()

from tfmra_py3 import tfmra_retracker
color_list=['powderblue','lightskyblue','deepskyblue','blue','darkblue','k','darkred','red','salmon','mistyrose']
#color_list=['k','k','k','k','k','k','k','k','k','k']
aval=[0.3,0.4,0.52,0.65,0.8,1,0.8,0.6,0.4,0.3]
BinWidthSnow=0.0618

#################################################################################3
# plot  waveforms GPR_0
#################################################################################3
begin=data_length-10

for x in range(10):
    plt.figure()
    plt.plot(g, selected_PWF[begin+round(x),:].transpose())
    plt.title('Waveform' + str(x))
    plt.savefig('Figures/0011/GPR_0/GPR_0_wf_'+str(x)+'.png',dpi=200)        

gate=tfmra_retracker(selected_PWF[begin+5,:],TrmfaThreshold, .9, 5)[0]
print(gate)


sd=2.83
gate_end=gate + sd/BinWidthSnow
print(gate_end)

y=BinWidthSnow*np.linspace(0,(round(gate_end)-round(gate)),round(sd/BinWidthSnow))
print(y)
plt.figure()

for x in range(10):
    plt.plot(selected_PWF[begin+round(x),round(gate):round(gate_end)].transpose(),y, color=color_list[x] ,alpha=aval[x],label=str(x))
plt.ylabel('depth (m)')
plt.xlabel('amplitude')
plt.ylim([sd,0])
plt.savefig('Figures/0011/GPR_0/GPR_0_wfs', dpi=200)
plt.show()



#################################################################################3

#################################################################################3
# plot  waveforms GPR_150
#################################################################################3
begin=round((data_length/3)*2)-5
for x in range(10):
    plt.figure()
    plt.plot(g, selected_PWF[begin+round(x),:].transpose())
    plt.title('Waveform' + str(x))
    plt.savefig('Figures/0011/GPR_150/GPR_150_wf_'+str(x)+'.png',dpi=200)   

gate=tfmra_retracker(selected_PWF[begin+5,:],TrmfaThreshold, .9, 5)[0]
print(gate)


sd=2.83
gate_end=gate + sd/BinWidthSnow
print(gate_end)


y=BinWidthSnow*np.linspace(0,(round(gate_end)-round(gate)),len(selected_PWF[begin+round(x),round(gate):round(gate_end)]))
plt.figure()
for x in range(10):
    plt.plot(selected_PWF[begin+round(x),round(gate):round(gate_end)].transpose(),y, color=color_list[x] ,alpha=aval[x],label=str(x))
plt.ylabel('depth (m)')
plt.xlabel('amplitude')
plt.ylim([sd,0])
plt.savefig('Figures/0011/GPR_150/GPR_150_wfs', dpi=200)
plt.show()
     
#################################################################################3

#################################################################################3
# plot  waveforms GPR_300
#################################################################################3
begin=round((data_length/3)*1)-5
for x in range(10):
    plt.figure()
    plt.plot(g, selected_PWF[begin+round(x),:].transpose())
    plt.title('Waveform' + str(x))
    plt.savefig('Figures/0011/GPR_300/GPR_300_wf_'+str(x)+'.png',dpi=200)     
    
gate=tfmra_retracker(selected_PWF[begin+5,:],TrmfaThreshold, .9, 5)[0]
print(gate)

    

sd=1.45
gate_end=gate + sd/BinWidthSnow
y=BinWidthSnow*np.linspace(0,(round(gate_end)-round(gate)),len(selected_PWF[begin+round(x),round(gate):round(gate_end)]))

plt.figure()
print(gate_end)
plt.figure()
for x in range(10):
    plt.plot(selected_PWF[begin+round(x),round(gate):round(gate_end)].transpose(),y, color=color_list[x] ,alpha=aval[x],label=str(x))
plt.ylabel('depth (m)')
plt.xlabel('amplitude')
plt.ylim([sd,0])
plt.savefig('Figures/0011/GPR_300/GPR_300_wfs', dpi=200)
plt.show()    
#################################################################################3

#################################################################################3
# plot  waveforms GPR_343
#################################################################################3

begin=0
for x in range(10):
    plt.figure()
    plt.plot(g, selected_PWF[begin+round(x),:].transpose())
    plt.title('Waveform' + str(x))
    plt.savefig('Figures/0011/GPR_343/GPR_343_wf_'+str(x)+'.png',dpi=200)   

gate=tfmra_retracker(selected_PWF[begin+5,:],TrmfaThreshold, .9, 5)[0]
print(gate)

sd=1.23
gate_end=gate + sd/BinWidthSnow
y=BinWidthSnow*np.linspace(0,(round(gate_end)-round(gate)),len(selected_PWF[begin+round(x),round(gate):round(gate_end)]))

print(gate_end)
plt.figure()
for x in range(10):
    plt.plot(selected_PWF[begin+round(x),round(gate):round(gate_end)].transpose(),y, color=color_list[x] ,alpha=aval[x],label=str(x))
plt.ylabel('depth(m)')
plt.xlabel('amplitude')
plt.ylim([sd,0])
plt.savefig('Figures/0011/GPR_343/GPR_343_wfs', dpi=200)
plt.show()       
#################################################################################3





print(data_length)
print(len(selected_Roll))
values_below = (selected_Roll < -1.5).sum()
print(values_below)
values_above=(selected_Roll > 1.5).sum()
print(values_above)
print(len(selected_Roll)-values_below-values_above)
print(2704-247-546)
print('percent of discarded data on selection due to roll: ' + str((((values_below+values_above))/len(selected_Roll))*100) + '%' )

Roll_mask=(selected_Roll<=1.5) & (selected_Roll>=-1.5)


#Roll sorting
Roll_sorted_selected_Hocog_ESA = selected_Hocog[Roll_mask]
Roll_sorted_selected_Htrmfa = H_trmfa_selected[Roll_mask]
Roll_sorted_selected_PWF = selected_PWF[Roll_mask]
Roll_sorted_selected_WindowDelay=selected_WindowDelay[Roll_mask]
Roll_sorted_selected_Hwgs84=selected_Hwgs84[Roll_mask]
Roll_sorted_selected_Lon=selected_Lon[Roll_mask]
Roll_sorted_selected_Lat=selected_Lat[Roll_mask]
Roll_sorted_selected_gate_rtck_array=gate_rtck_array[Roll_mask]


#I'll plot Hocog ESA and trmfa together
k=len(Roll_sorted_selected_Htrmfa) #length of roll sorted selection
X=np.linspace(0, int(k)-1, int(k))

plt.plot(X,Roll_sorted_selected_Hocog_ESA, color='lightsteelblue', label = 'h ocog (ESA)')
plt.plot(X,Roll_sorted_selected_Htrmfa, color='darkblue', label = 'h tfmra 50%')
plt.legend()
plt.xlabel('samples')
plt.ylabel('h (WGS84) [m]')
plt.savefig('Figures/0011/H trmfa roll sorted'+'.png',dpi=200)
plt.show() 


#Loop to plot every 1 selected roll sorted waveforms starting at 500
# k=len(Roll_sorted_selected_PWF)   #data_length after roll sort
# for x in range(0, 500, 1):
#     plt.figure()
#     plt.plot(g, Roll_sorted_selected_PWF[x,:].transpose())
#     plt.ylim(0,180)
#     #plt.plot(g, gate_rtck)
#     plt.title('Waveform' + str(x))
#     plt.savefig('Figures/0011/waveforms/ASIRAS_SCALED_'+ str(x) +'.png',dpi=200)  







