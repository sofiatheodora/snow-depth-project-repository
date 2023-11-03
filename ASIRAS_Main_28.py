#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: Sofia Thirslund
Uses importer from Henriette Skourup
re-tracks using tfmra provided by Alessandro Di Bella

"""

from ASIRAS_Importer_28 import *

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



#Plotting selected Lon Lat
plt.plot(Lon,Lat)
plt.title('flown route')
plt.savefig('Figures/122/LonLatRoute.png',dpi=170)
plt.show()


#Picking data closest to the line
from scipy import spatial

b2=[-52.205,-76.738]
b1=[-52.195,-76.7394]

coord_Lon=np.expand_dims(Lon, axis=0)
coord_Lat=np.expand_dims(Lat, axis=0)

coord=np.concatenate((coord_Lon.T,coord_Lat.T),axis=1)
print(np.shape(coord))

distance_b2,index_b2 = spatial.KDTree(coord).query(b2)
print(distance_b2)
print(index_b2)

distance_b1,index_b1 = spatial.KDTree(coord).query(b1)
print(distance_b1)
print(index_b1)

data_start=index_b2
data_end=index_b1
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
plt.savefig('Figures/122/DataSelectionOnFull.png',bbox_inches='tight',dpi=170)
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
plt.savefig('Figures/122/selected Hwgs as a colormap',dpi=180)
plt.show()


#Plotting selcted Lon Lat
plt.plot(selected_Lon,selected_Lat)
plt.title('Data Selection - near line b2-b1')
plt.savefig('Figures/122/DataSelection.png',dpi=170)
plt.show()




#Plotting selected data H_ocog (from ESA)
#Plotting Hocog as a function of time. Something is not right here.
selected_samples_var=np.linspace(0, int(data_length)-1, int(data_length))
plt.plot(selected_samples_var,Hocog[data_start:data_end],color='k')
plt.title('Hocog')
plt.xlabel('sample point')
plt.ylabel('Hocog (m)')
plt.title('ESA Hocog not roll sorted')
plt.savefig('Figures/122/ESA Hocog not roll sorted.png', dpi=200)
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
plt.savefig('Figures/122/RollonSelection.png')
plt.show()

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
print(np.sum(Roll_mask))


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
plt.savefig('Figures/122/H trmfa roll sorted'+'.png',dpi=200)
plt.show() 


#Adding corrections to Htrmfa

#Finding lead
from tfmra_py3 import tfmra_retracker
lead_index=data_start+46
plt.plot(g, PWF[lead_index,:].transpose())

gate_rtck, i_peak, AmpPmax, i_peakMax  = tfmra_retracker(PWF[lead_index,:],TrmfaThreshold, .8, 5)
H_trmfa_lead=Hwgs84[lead_index]-((gate_rtck-BinRef)*BinWidth+(WindowDelay[lead_index]*c_A)/2)
print(H_trmfa_lead-(-12.7))


Lon_lead=Lon[lead_index]
Lat_lead=Lat[lead_index]
print(Lon_lead)
print(Lat_lead)



WindowDelayRange=c_A*Roll_sorted_selected_WindowDelay

ALS_ASIRAS_offset=-(H_trmfa_lead-(-12.7))

H_trmfa_corr=H_trmfa_selected+ALS_ASIRAS_offset

f = open("ALS/b_line_ALS" , 'r')
b_line_ALS=np.genfromtxt(f, dtype=float)


ALS_ASIRAS_dif=H_trmfa_corr-b_line_ALS

plt.figure()
plt.plot(X,H_trmfa_corr, color='k', label = 'H trmfa')
plt.plot(X,b_line_ALS, color='g', label = 'ALS')
plt.title('trmfa (50%) vs ALS elevation along b line')
plt.legend()
plt.xlabel('sample number')
plt.ylabel('h (WGS 84)')
plt.savefig('Figures/122/trmfa (50) vs ALS elevation along b line'+'.png',dpi=200)
plt.show() 





# Caluclate Htrmfa_30_corr and Htrmfa_70_corr

H_trmfa_30=np.ones(data_length)
TrmfaThreshold=0.3
for x in range(data_length):
    gate_rtck, i_peak, AmpPmax, i_peakMax  = tfmra_retracker(selected_PWF[x,:],TrmfaThreshold, .8, 5)
    if gate_rtck != math.nan:
        H_trmfa_30[x]=selected_Hwgs84[x]-((gate_rtck-BinRef)*BinWidth+(selected_WindowDelay[x]*c_A)/2)
 
H_trmfa_30_corr=H_trmfa_30+ALS_ASIRAS_offset

H_trmfa_70=np.ones(data_length)
TrmfaThreshold=0.7
for x in range(data_length):
    gate_rtck, i_peak, AmpPmax, i_peakMax  = tfmra_retracker(selected_PWF[x,:],TrmfaThreshold, .8, 5)
    if gate_rtck != math.nan:
        H_trmfa_70[x]=selected_Hwgs84[x]-((gate_rtck-BinRef)*BinWidth+(selected_WindowDelay[x]*c_A)/2)
 
H_trmfa_70_corr=H_trmfa_70+ALS_ASIRAS_offset

plt.plot(X,b_line_ALS, color='g', label = 'ALS')
plt.plot(X,H_trmfa_30_corr,color='r',label='h trmfa 30%')
plt.plot(X,H_trmfa_corr, color='k', label = 'h trmfa 50%')
plt.plot(X,H_trmfa_70_corr,color='b',label='h trmfa 70%')
plt.legend()
plt.xlabel('sample number')
plt.ylabel('h (WGS84) [m]')
plt.savefig('Figures/122/trmfa (30)(50)(70) vs ALS elevation along b line'+'.png',dpi=200)
plt.show() 







#Create echogram

g_128=range(1,round(np.size(PWF[0,:])/2)+1)
order=np.ones((128,k))
for x in X:
    order[:,round(x)]=x
    
    
plt.figure()
cmap = plt.get_cmap('plasma')
plt.set_cmap(cmap)
for x in X:
    plt.scatter(order[:,round(x)],g_128,c=Roll_sorted_selected_PWF[round(x),:128], marker='s') 
cbar=plt.colorbar(orientation='vertical')
cbar.set_label('Waveform amplitude')
plt.title('Echogram along b-line')
plt.xlabel('sample number')
plt.ylabel('bin number')
plt.savefig('Figures/122/echogram_b_line'+'.png',dpi=200) 




#Loop to plot all selected roll sorted waveforms
k=len(Roll_sorted_selected_PWF)   #data_length after roll sort
print(k)
# for x in range(0, k, 1):
#     plt.figure()
#     plt.plot(g, Roll_sorted_selected_PWF[x,:].transpose())
#     plt.ylim(0,60)
#     #plt.plot(g, gate_rtck)
#     plt.title('Waveform' + str(x))
#     plt.savefig('Figures/122/line_b waveforms/ASIRAS_SCALED_'+ str(x) +'.png',dpi=200)  





#doublepeak waveform
TrmfaThreshold=0.5

index_dp_1=99

print(Roll_sorted_selected_Lon[index_dp_1])
print(Roll_sorted_selected_Lat[index_dp_1])

#Plot roll sorted Lon Lat
plt.figure()
plt.plot(Roll_sorted_selected_Lon,Roll_sorted_selected_Lat)
plt.scatter(Roll_sorted_selected_Lon[index_dp_1],Roll_sorted_selected_Lat[index_dp_1])
plt.title('Location of doublepeak 1')
plt.savefig('Figures/122/Location of doublepeak.png')
plt.show()

#Plot known waveform
plt.figure()
plt.plot(g, Roll_sorted_selected_PWF[index_dp_1,:].transpose()) #Her har jeg bare brugt PWF i stedet for P_wf?
plt.ylim(0,25)
#plt.plot(g, gate_rtck)
plt.title('Waveform of double peak 1')
plt.savefig('Figures/122/DoublePeak 1 wave form.png',dpi=200)   


#RETRACKING the known double peak with threshold...
TrmfaThreshold=0.5
gate_rtck, i_peak, AmpPmax, i_peakMax  = tfmra_retracker(Roll_sorted_selected_PWF[index_dp_1,:],TrmfaThreshold, .8, 5)
print(gate_rtck)
print(Roll_sorted_selected_PWF[index_dp_1,40])
print(Roll_sorted_selected_PWF[index_dp_1,41])

print('leading edge position gate retracked by trmfa ' + str(gate_rtck))
print('index of retracked peak using tfrma ' + str(i_peak))
print('Max Amplitude of P trmfa ' + str(i_peakMax)) #equal to AmpPmax

#print(selected_Lon) #used to put in manually to ReadALS
#print(selected_Lat)

# ...and OCOG Retracker
LEP, g, COG, A_OCOG, W_OCOG = OCOGretracker(Roll_sorted_selected_PWF[index_dp_1,:], 0)
     #   ocogWidth[pp] = W_OCOG
     #   ocogRT[pp] = Hwgs84[1110]-((LEP-BinRef)*BinWidth+WindowDelayRange/2)
print('leading edge position using ocog ' + str(LEP))



WindowDelayRange_dp_1=Roll_sorted_selected_WindowDelay[index_dp_1]*c_A 
print('wdr ' + str(WindowDelayRange_dp_1/2))  
print('Hwgs ' + str(Roll_sorted_selected_Hwgs84[index_dp_1])) 
print("why is Hwgs lower than wdr?")

#trmfa (original)
H_trmfa_dp = Roll_sorted_selected_Hwgs84[index_dp_1]-((gate_rtck-BinRef)*BinWidth+WindowDelayRange_dp_1/2)
print('calculated trmfa H ' + str(H_trmfa_dp))

#(ocog)
ocog_height_doublepeak = Roll_sorted_selected_Hwgs84[index_dp_1]-((LEP-BinRef)*BinWidth+WindowDelayRange_dp_1/2)
print('calculated H_ocog ' + str(ocog_height_doublepeak))

ESA_ocog_doublepeak=Roll_sorted_selected_Hocog_ESA[index_dp_1]
print('ESA H_ocog ' + str(ESA_ocog_doublepeak))


#ALS data (dp = wf index_dp_1)
closest_coord=[-52.194689, -76.736812] #must be adjusted to dp_1
ALS_dp=-13.19  #must be adjusted to dp_1

diff=H_trmfa_dp-ALS_dp
print(diff) #ok, let's see if the off-set between the two fits this diff?

# Offset between ALS-ASIRAS () elevations: 3.94 m (Skyblue overflight).


remaining_diff=diff+ALS_ASIRAS_offset
print('diff using trmfa50 ' + str(remaining_diff))

print('diff using ESA ocog ' +str( Roll_sorted_selected_Hocog_ESA[index_dp_1]-ALS_dp+ALS_ASIRAS_offset))

ASIRAS_dp=H_trmfa_dp+ALS_ASIRAS_offset
print(ASIRAS_dp)
#MSS er ift. til WGS84:
    
MSS_dp=-16.449 #check if closest to new dp

dp_ASIRAS_MSS=ASIRAS_dp-MSS_dp #ref to MSS instead of HWGS84

print(dp_ASIRAS_MSS)



#Comparing to inSitu GPR500
BinWidthSnow=0.0618
color_list=['powderblue','lightskyblue','deepskyblue','blue','darkblue','k','darkred','red','salmon','mistyrose']
#color_list=['k','k','k','k','k','k','k','k','k','k']
aval=[0.3,0.4,0.52,0.65,0.8,1,0.8,0.6,0.4,0.3]
#################################################################################3
# plot  waveforms GPR_500
#################################################################################3
begin=data_length-5
for x in range(10):
    plt.figure()
    plt.plot(g, PWF[index_b2+begin+round(x),:].transpose())
    plt.title('Waveform' + str(x))
    plt.savefig('Figures/122/GPR_500/GPR_500_wf_'+str(x)+'.png',dpi=200)        

gate = tfmra_retracker(PWF[index_b2+begin+4,:],0.5,0.6,5)[0]
print(gate)
sd=1.01
gate_end=gate + sd/BinWidthSnow
print(gate_end)
y=BinWidthSnow*np.linspace(0,(round(gate_end)-round(gate)),(round(gate_end)-round(gate)))

plt.figure()
for x in range(10):
    plt.plot(PWF[index_b2+begin+round(x),round(gate):round(gate_end)].transpose(),y, color=color_list[x] ,alpha=aval[x],label=str(x))
plt.ylabel('depth (m)')
plt.xlabel('amplitude')
plt.ylim([sd,0])
plt.savefig('Figures/122/GPR_500/GPR_500_wfs', dpi=200)
plt.show()



#(my trmfa)

from tfmra_py3_MyVersion import tfmra_retracker
from tfmra_py3_MyVersion import find_first_peak
from tfmra_py3_MyVersion import find_second_peak


#Remember to adjust search area
first_peak_i,first_peak_value=find_first_peak(Roll_sorted_selected_PWF[index_dp_1,:], 0.7*30)
print(first_peak_i)
print(first_peak_value)

second_peak_i,second_peak_value=find_second_peak(Roll_sorted_selected_PWF[index_dp_1,:], 0.2*30)
print(second_peak_i)
print(second_peak_value)

print(tfmra_retracker(Roll_sorted_selected_PWF[index_dp_1,:],TrmfaThreshold, .8, .2, 5))
treshold_01_01=tfmra_retracker(Roll_sorted_selected_PWF[index_dp_1,:],TrmfaThreshold, .8, .2, 5)[3]
treshold_01_02=tfmra_retracker(Roll_sorted_selected_PWF[index_dp_1,:],TrmfaThreshold, .8, .2, 5)[7]

gate_rtck_01_01=tfmra_retracker(Roll_sorted_selected_PWF[index_dp_1,:],TrmfaThreshold, .8, .2, 5)[0]

gate_rtck_01_02=tfmra_retracker(Roll_sorted_selected_PWF[index_dp_1,:],TrmfaThreshold, .8, .2, 5)[4]


WF_dp_1=Roll_sorted_selected_PWF[index_dp_1,:]

#Plot known doublepeak waveform with retracked gates.
plt.figure()
plt.plot(g,WF_dp_1.transpose(),color='k')
plt.axvline(x=gate_rtck_01_01,color='b', ls = ':')
plt.axvline(x=gate_rtck_01_02,color='r', ls = ':')
plt.axhline(y=treshold_01_01, color='b', ls = ':')
plt.axhline(y=treshold_01_02,color='r', ls = ':')
plt.xlim(30,60)
plt.savefig('Figures/122/DoublePeak 1 wave form with retracked gates.png',dpi=200)   
plt.show() 


#Double peak 2
index_dp_2=82

plt.figure()
plt.plot(g, Roll_sorted_selected_PWF[index_dp_2,:].transpose())
plt.ylim(0,40)
plt.title('Waveform of double peak 2')
plt.savefig('Figures/122/DoublePeak 2 wave form.png',dpi=200)   

print(find_first_peak(Roll_sorted_selected_PWF[index_dp_2,:], 0.8*30))
print(find_second_peak(Roll_sorted_selected_PWF[index_dp_2,:], 0.2*30))
print(tfmra_retracker(Roll_sorted_selected_PWF[index_dp_2,:],TrmfaThreshold, .8, .2, 5))
treshold_02_01=tfmra_retracker(Roll_sorted_selected_PWF[index_dp_2,:],TrmfaThreshold, .8, .2, 5)[3]
treshold_02_02=tfmra_retracker(Roll_sorted_selected_PWF[index_dp_2,:],TrmfaThreshold, .8, .2, 5)[7]

gate_rtck_02_01=tfmra_retracker(Roll_sorted_selected_PWF[index_dp_2,:],TrmfaThreshold, .8, .2, 5)[0]

gate_rtck_02_02=tfmra_retracker(Roll_sorted_selected_PWF[index_dp_2,:],TrmfaThreshold, .8, .2, 5)[4]


WF_dp_2=Roll_sorted_selected_PWF[index_dp_2,:]

#Plot known doublepeak waveform with retracked gates.
plt.figure()
plt.plot(g,WF_dp_2.transpose(),color='k')
plt.axvline(x=gate_rtck_02_01,color='b', ls = ':')
plt.axvline(x=gate_rtck_02_02,color='r', ls = ':')
plt.axhline(y=treshold_02_01, color='b', ls = ':')
plt.axhline(y=treshold_02_02,color='r', ls = ':')
plt.xlim(30,60)
plt.savefig('Figures/122/DoublePeak 2 wave form with retracked gates.png',dpi=200)   
plt.show()


#Double peak 3
index_dp_3=43

plt.figure()
plt.plot(g, Roll_sorted_selected_PWF[index_dp_3,:].transpose())
plt.ylim(0,40)
plt.title('Waveform of double peak 3')
plt.savefig('Figures/122/DoublePeak 3 wave form.png',dpi=200)   

print(find_first_peak(Roll_sorted_selected_PWF[index_dp_3,:], 0.8*25))
print(find_second_peak(Roll_sorted_selected_PWF[index_dp_3,:], 0.2*25))
print(tfmra_retracker(Roll_sorted_selected_PWF[index_dp_3,:],TrmfaThreshold, .8, .2, 5))
treshold_03_01=tfmra_retracker(Roll_sorted_selected_PWF[index_dp_3,:],TrmfaThreshold, .8, .2, 5)[3]

treshold_03_02=tfmra_retracker(Roll_sorted_selected_PWF[index_dp_3,:],TrmfaThreshold, .8, .2, 5)[7]

gate_rtck_03_01=tfmra_retracker(Roll_sorted_selected_PWF[index_dp_3,:],TrmfaThreshold, .8, .2, 5)[0]

gate_rtck_03_02=tfmra_retracker(Roll_sorted_selected_PWF[index_dp_3,:],TrmfaThreshold, .8, .2, 5)[4]


WF_dp_3=Roll_sorted_selected_PWF[index_dp_3,:]

#Plot known doublepeak waveform with retracked gates. 
plt.figure()
plt.plot(g,WF_dp_3.transpose(),color='k')
plt.axvline(x=gate_rtck_03_01,color='b', ls = ':')
plt.axvline(x=gate_rtck_03_02,color='r', ls = ':')
plt.axhline(y=treshold_03_01, color='b', ls = ':')
plt.axhline(y=treshold_03_02,color='r', ls = ':')
plt.xlim(0,100)
plt.xlim([30,60])
plt.savefig('Figures/122/DoublePeak 3 wave form with retracked gates.png',dpi=200)   
plt.show()


#plot 3 double peaks together
plt.figure()
plt.plot(g,WF_dp_1.transpose(),'--', color='k', label = '01')
plt.plot(g,WF_dp_2.transpose(),':', color='k', label = '02')
plt.plot(g,WF_dp_3.transpose(),color='k', label = '03')
plt.legend()
plt.xlim([30,60])
plt.xlabel('gates')
plt.ylabel('amplitude')
plt.savefig('Figures/122/3doublepeaks', dpi=200)
plt.show()


#Create echogram with (trmfa50%) retracked gates marked. Also dp retracked second gates marked.

g_128=range(1,round(np.size(PWF[0,:])/2)+1)
order=np.ones((128,k))
for x in X:
    order[:,round(x)]=x
    
    
plt.figure()
cmap = plt.get_cmap('plasma')
plt.set_cmap(cmap)

for x in X:
    plt.scatter(order[:,round(x)],g_128,c=PWF[index_b2+round(x),:128],  marker='s') 
cbar=plt.colorbar(orientation='vertical')
cbar.set_label('Waveform amplitude')
plt.plot(X,Roll_sorted_selected_gate_rtck_array,color='k', alpha=0.8)
plt.scatter(index_dp_3,gate_rtck_03_02, marker='x',color='c', label = 'a')
plt.scatter(index_dp_2,gate_rtck_02_02, marker='x',color='k', label = 'b')
plt.scatter(index_dp_1,gate_rtck_01_02, marker='x',color='white', label = 'c')
plt.legend()

plt.xlabel('sample number')
plt.ylabel('bin number')
plt.savefig('Figures/122/echogram_b_line_with_retracked_gates'+'.png',dpi=200) 

plt.show()



#Calculating snow depth
C_S=1.688E8
print('C_S =' + str(C_S))

BinWidthSnow=0.0618

#dp 01
Lat_dp_01=Roll_sorted_selected_Lat[index_dp_1]
Lon_dp_01=Roll_sorted_selected_Lon[index_dp_1]
dp_01_coord=[Lat_dp_01,Lon_dp_01]
print('dp 01 coord:' + str(dp_01_coord) )

SD_01=BinWidthSnow*(gate_rtck_01_02-gate_rtck_01_01)
print('Snow depth dp 01: ' + str(SD_01) + 'm')

#dp 02

Lat_dp_02=Roll_sorted_selected_Lat[index_dp_2]
Lon_dp_02=Roll_sorted_selected_Lon[index_dp_2]
dp_02_coord=[Lat_dp_02,Lon_dp_02]
print('dp 02 coord:' + str(dp_02_coord) )

SD_02=BinWidthSnow*(gate_rtck_02_02-gate_rtck_02_01)
print('Snow depth dp 02: ' + str(SD_02) + 'm')

#dp 03
Lat_dp_03=Roll_sorted_selected_Lat[index_dp_3]
Lon_dp_03=Roll_sorted_selected_Lon[index_dp_3]
dp_03_coord=[Lat_dp_03,Lon_dp_03]
print('dp 03 coord:' + str(dp_03_coord) )

SD_03=BinWidthSnow*(gate_rtck_03_02-gate_rtck_03_01)
print('Snow depth dp 03: ' + str(SD_03) + 'm')



#Finding data closest to transect 1
r2=[-52.202,-76.736]
r1=[-52.191,-76.7367]

distance_t1,index_t1 = spatial.KDTree(coord).query([-52.191,-76.7367])
print(distance_t1)
print(index_t1)
distance_t2,index_t2 = spatial.KDTree(coord).query(r2)
print(distance_t2)
print(index_t2)

pad=index_t2-index_t1
print(pad)
t1_Lon=Lon[index_t1:index_t1+pad]
t1_Lat=Lat[index_t1:index_t1+pad]
t1_Hwgs84=Hwgs84[index_t1:index_t1+pad]

t1_PWF=PWF[index_t1:index_t1+pad]
t1_WindowDelay=WindowDelay[index_t1:index_t1+pad]
t1_Hocog=Hocog[index_t1:index_t1+pad]
t1_Roll=Roll[index_t1:index_t1+pad]
data_length=pad
selected_samples_var=np.linspace(0, int(data_length)-1, int(data_length))
plt.figure()
plt.plot(selected_samples_var,t1_Roll,color='k')
plt.xlabel('samples')
plt.ylabel('Roll angle (\N{DEGREE SIGN})')
plt.savefig('Figures/122/Rollont1.png')
plt.show()


# for x in range(0, pad, 1):
#     plt.figure()
#     plt.plot(g, PWF[index_t1+x,:].transpose())
#     plt.ylim(0,60)
#     #plt.plot(g, gate_rtck)
#     plt.title('Waveform' + str(x))
#     plt.savefig('Figures/122/t1/ASIRAS_SCALED_'+ str(x) +'.png',dpi=200)  



#Plotting selected Lon Lat
plt.figure()
plt.plot(Lon,Lat, color='k')
plt.plot(t1_Lon,t1_Lat, color='r', linewidth=5)
plt.xlabel('Lon (\N{DEGREE SIGN})')
plt.ylabel('Lat (\N{DEGREE SIGN})')
plt.savefig('Figures/122/t1OnFull.png',bbox_inches='tight',dpi=170)
plt.show()


#Double peak 4
dp4=20
index_dp_4=index_t1

plt.figure()
plt.plot(g, PWF[index_dp_4,:].transpose())
plt.ylim(0,40)
plt.title('Waveform of double peak 4')
plt.savefig('Figures/122/DoublePeak 4 wave form.png',dpi=200)   

print(find_first_peak(PWF[index_dp_4,:], 0.8*25))
print(find_second_peak(PWF[index_dp_4,:], 0.2*25))
print(tfmra_retracker(PWF[index_dp_4,:],TrmfaThreshold, .8, .2, 5))
treshold_04_01=tfmra_retracker(PWF[index_dp_4,:],TrmfaThreshold, .8, .2, 5)[3]

treshold_04_02=tfmra_retracker(PWF[index_dp_4,:],TrmfaThreshold, .8, .2, 5)[7]

gate_rtck_04_01=tfmra_retracker(PWF[index_dp_4,:],TrmfaThreshold, .8, .2, 5)[0]

gate_rtck_04_02=tfmra_retracker(PWF[index_dp_4,:],TrmfaThreshold, .8, .2, 5)[4]


WF_dp_4=PWF[index_dp_4,:]

#Plot known doublepeak waveform with retracked gates. 
plt.figure()
plt.plot(g,WF_dp_4.transpose(),color='k')
plt.axvline(x=gate_rtck_04_01,color='b', ls = ':')
plt.axvline(x=gate_rtck_04_02,color='r', ls = ':')
plt.axhline(y=treshold_04_01, color='b', ls = ':')
plt.axhline(y=treshold_04_02,color='r', ls = ':')
plt.xlim(0,100)
plt.xlim([30,60])
plt.savefig('Figures/122/DoublePeak 4 wave form with retracked gates.png',dpi=200)   
plt.show()



#Double peak 5
dp5=108
index_dp_5=index_t1+dp5

plt.figure()
plt.plot(g, PWF[index_dp_5,:].transpose())
plt.ylim(0,40)
plt.title('Waveform of double peak 5')
plt.savefig('Figures/122/DoublePeak 5 wave form.png',dpi=200)   

print(find_first_peak(PWF[index_dp_5,:], 0.8*25))
print(find_second_peak(PWF[index_dp_5,:], 0.3*25))
print(tfmra_retracker(PWF[index_dp_5,:],TrmfaThreshold, .8, .2, 5))
treshold_05_01=tfmra_retracker(PWF[index_dp_5,:],TrmfaThreshold, .8, .2, 5)[3]

treshold_05_02=tfmra_retracker(PWF[index_dp_5,:],TrmfaThreshold, .8, .2, 5)[7]

gate_rtck_05_01=tfmra_retracker(PWF[index_dp_5,:],TrmfaThreshold, .8, .2, 5)[0]

gate_rtck_05_02=tfmra_retracker(PWF[index_dp_4,:],TrmfaThreshold, .8, .2, 5)[4]


WF_dp_5=PWF[index_dp_5,:]

#Plot known doublepeak waveform with retracked gates. 
plt.figure()
plt.plot(g,WF_dp_5.transpose(),color='k')
plt.axvline(x=gate_rtck_05_01,color='b', ls = ':')
plt.axvline(x=gate_rtck_05_02,color='r', ls = ':')
plt.axhline(y=treshold_05_01, color='b', ls = ':')
plt.axhline(y=treshold_05_02,color='r', ls = ':')
plt.xlim(0,100)
plt.xlim([30,60])
plt.savefig('Figures/122/DoublePeak 5 wave form with retracked gates.png',dpi=200)   
plt.show()


#dp 04
Lat_dp_04=Lat[index_dp_4]
Lon_dp_04=Lon[index_dp_4]
dp_04_coord=[Lat_dp_04,Lon_dp_04]
print('dp 04 coord:' + str(dp_04_coord) )

SD_04=BinWidthSnow*(gate_rtck_04_02-gate_rtck_04_01)
print('Snow depth dp 04: ' + str(SD_04) + 'm')

#dp 05
Lat_dp_05=Lat[index_dp_5]
Lon_dp_05=Lon[index_dp_5]
dp_05_coord=[Lat_dp_05,Lon_dp_05]
print('dp 05 coord:' + str(dp_05_coord) )

SD_05=BinWidthSnow*(gate_rtck_05_02-gate_rtck_05_01)
print('Snow depth dp 05: ' + str(SD_05) + 'm')

k=len(t1_PWF) #length of roll sorted selection
X=np.linspace(0, int(k)-1, int(k))

#echogram along transect 1


from tfmra_py3 import tfmra_retracker
TrmfaThreshold=0.5


gate_rtck_array=np.ones(data_length)
for x in range(data_length):
    gate_rtck  = tfmra_retracker(PWF[index_t1+x,:],TrmfaThreshold, .8, 5)[0]
    if gate_rtck != math.nan:
        gate_rtck_array[x]=gate_rtck


g_128=range(1,round(np.size(PWF[0,:])/2)+1)
order=np.ones((128,k))
for x in X:
    order[:,round(x)]=x
    
    
plt.figure()
cmap = plt.get_cmap('plasma')
plt.set_cmap(cmap)
for x in X:
    plt.scatter(order[:,round(x)],g_128,c=t1_PWF[round(x),:128], marker='s') 
cbar=plt.colorbar(orientation='vertical')
cbar.set_label('Waveform amplitude')
plt.plot(X,gate_rtck_array,color='k', alpha=0.8)
plt.scatter(dp4,gate_rtck_04_02, marker='x',color='pink', label = 'd')
plt.scatter(dp5,gate_rtck_05_02, marker='x',color='maroon', label = 'e')
plt.legend()

plt.xlabel('sample number')
plt.ylabel('bin number')
plt.savefig('Figures/122/echogram_transect_1_with_retracked_gates'+'.png',dpi=200) 

plt.show()








#MAKING MAP
IS_3_Lat = -76.7767
IS_3_Lon = -52.2400


x_val=[r2[0],r1[0]]
y_val=[r2[1],r1[1]]
x_valb=[b2[0],b1[0]]
y_valb=[b2[1],-76.7387]

SD_01=np.around(SD_01,2)
SD_02=np.around(SD_02,2)
SD_03=np.around(SD_03,2)
SD_04=np.around(SD_04,2)
SD_05=np.around(SD_05,2)
import matplotlib.path as mpath
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

plt.style.use('ggplot')

fig = plt.figure(figsize=[10, 5])
ax = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())

fig.subplots_adjust(bottom=0.05, top=0.95,
                      left=0.04, right=0.95, wspace=0.02)



  # Limit the map to -60 degrees latitude and below.
ax.set_extent([-52.188,-52.21, -76.730, -76.746,], ccrs.PlateCarree())
ax.stock_img()

ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.RIVERS)


gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.05, linestyle='--', zorder=0)


#gl.xformatter = LONGITUDE_FORMATTER
#gl.yformatter = LATITUDE_FORMATTER


plt.plot(Lon,Lat,
          color='k', linestyle = ':', linewidth = 1, alpha=0.9,
          transform=ccrs.Geodetic(), zorder=1
          )

# plt.plot(Roll_sorted_selected_Lon,Roll_sorted_selected_Lat, linewidth = 3,
#           color='k',
#           transform=ccrs.Geodetic(),
#           )
#plt.scatter(IS_3_Lon,IS_3_Lat, color = 'k',
#            transform=ccrs.Geodetic(), marker = 'x')

plt.plot(x_val, y_val, color = 'green',transform=ccrs.Geodetic(), zorder=2,label='transect 1')
plt.plot(x_valb, y_valb, color = 'red',transform=ccrs.Geodetic(),zorder=2, label='transect 2')

plt.scatter(dp_03_coord[1],dp_03_coord[0], color = 'c', s = 40,marker='x', zorder=3,
            transform=ccrs.Geodetic(), label='a: '  + str(SD_03) + ' m')
plt.scatter(dp_02_coord[1],dp_02_coord[0], color = 'k', s = 40,marker='x', zorder=3,
            transform=ccrs.Geodetic(), label='b: '  + str(SD_02) + ' m')
plt.scatter(dp_01_coord[1],dp_01_coord[0], color = 'white', s = 40, marker='x',zorder=3,
            transform=ccrs.Geodetic(),label='c: '  + str(SD_01) + ' m')

plt.scatter(dp_04_coord[1],dp_04_coord[0], color = 'pink', s = 40,marker='x', zorder=3,
            transform=ccrs.Geodetic(), label='d: '  + str(SD_04) + ' m')

plt.scatter(dp_05_coord[1],dp_05_coord[0], color = 'maroon', s = 40,marker='x', zorder=3,
            transform=ccrs.Geodetic(), label='e: '  + str(SD_05) + ' m')

plt.legend(bbox_to_anchor=(1.20, 1), loc='upper left', borderaxespad=0)
plt.savefig('Figures/122/Route_dps.png',bbox_inches='tight', dpi=300)

plt.show()

print(b1)


