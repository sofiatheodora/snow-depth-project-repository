#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Read karen 
Created on Thu Jun 15 10:00:51 2017

@author: Sofia Thirslund
based on: ReadKarenL1B by Henriette Skourup
"""

#Hwgs84 is called FlightAltitude

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt 
from mpl_toolkits.basemap import pyproj as pyproj
plt.style.use('ggplot')
import sys
#sys.path.append('../../PySub/')

from scipy.interpolate import LinearNDInterpolator as LNDI
from ocog_retracker import OCOGretracker
import os
import argparse
import datetime

def parse_arguments():
        parser = argparse.ArgumentParser(description='reading and retracking KAren data')
        parser.add_argument("-i","--inputfile", type=str, help="Input metasensing nc file")
        parser.add_argument("-o","--outputdir", type=str,default='', help="Please define the output directory ")
        args = parser.parse_args()
        return args
    

if __name__ == '__main__':
    
#     args = parse_arguments()
#     if args.inputfile:
#         print('Analyzing '+ args.inputfile)
#     else:   
#         sys.exit('No input file please specify "-i file.nc')
            
#     ifile =args.inputfile
#     odir  =  args.outputdir

#'KAREN/KAR_OPER_Level1b_20171228T185413_20171228T185501_levc.nc' is the file that contains line b
    ifile ='KAREN/KAR_OPER_Level1b_20171229T124335_20171229T124423_levc.nc'
    odir  =  'Figures/KAREN/'
    
    #try:
        
    ncfile = Dataset(ifile,'r')
    print(ncfile)
        
    Lat = ncfile.variables["latitude_ka"][:]
    Lon = ncfile.variables["longitude_ka"][:]
    #VerticalVel = ncfile.variables["VerticalVel"][:] #m/s
    Wf = ncfile.variables["hr_power_waveform_ka"][:]  #Watts
    phase = ncfile.variables["hr_phase_waveform_ka"][:]  #
    coh = ncfile.variables["hr_coh_waveform_ka"][:]  #Watts
    FlightAltitude = ncfile.variables["com_altitude_ka"][:] #m 
    Range = ncfile.variables["range"][:] #Range [m]
    Roll = ncfile.variables["off_nadir_roll_angle_ka"][:] #long_name = "Roll Angle" units = "arc degree [deg]";
    leap = 18.
    GPSTime = ncfile.variables["time_ka"][:] - leap # long_name = "Global Positioning System Time" units = "second [s]";
    baselineVer = ncfile.variables["BaselineVer"][:]
    print(baselineVer)
    
    pitch = ncfile.variables["off_nadir_pitch_angle_pf_ka"][:]
    #roll = ncfile.variables["off_nadir_roll_angle_ka"][:]
    yaw = ncfile.variables["off_nadir_yaw_angle_pf_ka"][:]
    ncfile.close()
    
    #gimp = Ftri(Lon,Lat)
    #c = 299792458. # m / s
    
    #travel = c*RD/2.
    
    Pollat0  = 90.
    Pollon0  = np.nanmean(Lon)
    Pollatts = np.nanmean(Lat)
    ell = 'WGS84'
    Proj=pyproj.Proj(proj='stere', ellps=ell, lat_0=Pollat0, lon_0=Pollon0, lat_ts=Pollatts,x_0=0,y_0=0,k_0=1.0) 
    x, y = Proj(Lon,Lat)
    dist = np.zeros_like(x)
    srf = np.zeros_like(x)
    ocog = np.zeros_like(x)
    ph = np.zeros_like(x)
    ch = np.zeros_like(x)
    ocogRT = np.zeros_like(x)
    RT = np.zeros_like(x)
    RTmin = np.zeros_like(x)
    for ii in range(len(x)):
        if ii >0:
            dist[ii] = np.sqrt((x[ii]-x[ii-1])**2+(y[ii]-y[ii-1])**2)
        srf[ii] = FlightAltitude[ii]-Range[Wf[ii,:]==max(Wf[ii,:])]
        try: 
            LEP, g, COG, A_OCOG, W_OCOG = OCOGretracker(Wf[ii,:], 100, 0)
            ocog[ii] = np.interp(LEP,np.arange(len(Wf[ii,:])), Range)
            ocogRT[ii] = FlightAltitude[ii]-ocog[ii]
            ph[ii] = np.interp(LEP,np.arange(len(Wf[ii,:])), phase[ii,:])
            ch[ii] = np.interp(LEP,np.arange(len(Wf[ii,:])), coh[ii,:])
            noise = 1000.*np.median(Wf[ii,:]) 
            maxtr = Range[Wf[ii,:]==max(Wf[ii,:])]
            mintr = Range[Wf[ii,:]>noise][0]
            trmsk = (Range>=mintr) & (Range<=maxtr)
            RT[ii]= FlightAltitude[ii]-np.interp(max(Wf[ii,:])/2.,Wf[ii,trmsk] , Range[trmsk])
            RTmin[ii]= FlightAltitude[ii]-mintr
        except: pass
        
    atdist = np.cumsum(dist)
     
    tmsk = (Range>min(Range)) & (Range<max(Range))
    
    
    fig=plt.figure()
    plt.subplot(311)
    Wf[Wf==-9999.]=np.nan
    sub=10
    TT,DD = np.meshgrid(Range[tmsk],atdist[:-1:sub])
    data = np.log(Wf[:-1:sub, tmsk].T)
    plt.pcolor(DD.T, TT.T, data, vmin=np.nanmin(data), vmax=np.nanmax(data))
    plt.plot(atdist, ocog)
    plt.gca().invert_yaxis()
    plt.xlabel('Distance along track [m]')
    plt.ylabel('Range [m]')
    plt.subplot(323)
    plt.plot(atdist[srf<(FlightAltitude-50.)], srf[srf<(FlightAltitude-50.)],label='Karen max')
    plt.plot(atdist[(RT<(FlightAltitude-50.)) & (RT!=0)], RT[(RT<(FlightAltitude-50.)) & (RT!=0)], 'b',label='Karen 50%')
    plt.plot(atdist[(RTmin<(FlightAltitude-50.)) & (RTmin!=0)], RTmin[(RTmin<(FlightAltitude-50.)) & (RTmin!=0)], 'c',label='Karen 0%')
    plt.plot(atdist[ocogRT<(FlightAltitude-50.)], ocogRT[ocogRT<(FlightAltitude-50.)], 'm',label='OCOG')
    #plt.plot(atdist,gimp, 'k', label='GIMP')
    plt.legend(loc='best')
    plt.xlabel('Distance along track [m]')
    plt.ylabel('Elevation [m]')
    plt.subplot(324)
    plt.scatter(Lon, Lat, s=20,c=atdist, lw=0)
    plt.savefig(odir + '/' +  'KarenTest01_4423.png',dpi=200)    
    
    # tmpax = plt.gca()
    # ylim = tmpax.get_ylim()
    # xlim = tmpax.get_xlim()
    # data = np.loadtxt('GR-IC-CA.SVA', delimiter='\n', dtype=str)
    # tmp = np.empty((len(data),5))
    # for ii in range(len(data)):
    #     t = map(float,data[ii].split())
    #     tmp[ii,:len(t)] = t 
            
    # ii =0;
    # while ii < len(data):
    #     lat = tmp[ii+1:ii+int(tmp[ii,0]),0]
    #     lon = tmp[ii+1:ii+int(tmp[ii,0]),1]
    #     ii = ii+int(tmp[ii,0])+1
    #     plt.plot(lon, lat,'k')
    # tmpax.set_xlim((np.floor(xlim[0]-0.5), np.ceil(xlim[1])+0.5))
    # tmpax.set_ylim((np.floor(ylim[0]-0.5), np.ceil(ylim[1])+0.5))

    cbar = plt.colorbar()
    cbar.set_label('Distance along track [m]', rotation=90,fontsize=12)
    plt.subplot(313)
    plt.plot(atdist,pitch, 'b', label='Pitch') 
    plt.plot(atdist,Roll, 'r', label='Roll') 
    plt.plot(atdist,yaw, 'k', label='yaw') 
    plt.legend(loc='best')
    fig.set_size_inches(30,20)
    plt.savefig(odir + '/KarenTest02_4423.png', bbox_inches='tight', pad_inches=0)
    plt.show()
    
    day = datetime.datetime(int(os.path.basename(ifile)[-39:-35]), int(os.path.basename(ifile)[-35:-33]), int(os.path.basename(ifile)[-33:-31]), 0, 0)
    doy = (day - datetime.datetime(2000,1,1)).days +1
    sec = (GPSTime/(3600.*24.)-np.floor( GPSTime/(3600.*24.)))*(24.*3600)
    
    
    np.savetxt(odir+os.path.basename(ifile)[:-3]+'.txt', np.column_stack((np.ones_like(GPSTime)*doy, Lat, Lon, sec,  ocogRT, Roll, pitch, yaw, ph, ch)),  \
               fmt='%i %4.7f %4.7f %4.5f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f', \
               header='DayNo*,   Latitude,   Longitude,   seconds of day [UTC],  OCOG re-tracked elevation,  roll angle (degrees),  Pitch angle (degrees),  Yaw angle (degrees), Phase, coherence', \
               comments='#Readout based on the OCOG RT by adia@space.dtu.dk and implemented by ssim@space.dtu.dk\n')
#    except:
#        print 'the file failed'


# #Plot single waveform:
    
#     g = range(1,np.size(Wf[5,:])+1)
#     plt.figure()
#     plt.plot(g, Wf[7,:].transpose())
#     plt.title('Waveform')
#     plt.savefig('Figures/KAREN/line_r waveforms/KARENtestWF' + '.png',dpi=200)  

# #Plot Lon Lat
# plt.plot(Lon,Lat)
# plt.show()


#Picking data closest to the line
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

data_start=index_r2-5 #end and start point of data selection should be switched if going from right to left
data_end=index_r1+5
data_length=data_end-data_start
print(data_length) #sample size

selected_Lon=Lon[data_start:data_end]
selected_Lat=Lat[data_start:data_end]
selected_FlightAltitude=FlightAltitude[data_start:data_end]
selected_atdist=atdist[data_start:data_end]
selected_Roll=Roll[data_start:data_end]
selected_pitch=pitch[data_start:data_end]
selected_yaw=yaw[data_start:data_end]
selected_Range=Range[data_start:data_end]
selected_Wf=Wf[data_start:data_end]
selected_Hocog_ESA=ocog[data_start:data_end]

#Plotting selected Lon Lat
plt.plot(Lon,Lat, color='k')
plt.plot(selected_Lon,selected_Lat, color='r', linewidth=5)
#plt.title('Data Selection')
plt.savefig('Figures/KAREN/DataSelectionOnFull_file_4423',dpi=170)
plt.show()


#Plotting roll pitch and yaw on selection



#cbar = plt.colorbar()
#cbar.set_label('Distance along track [m]', rotation=90,fontsize=12)
plt.plot(range(data_length),selected_pitch, 'b', label='Pitch') 
plt.plot(range(data_length),selected_Roll, 'r', label='Roll') 
plt.plot(range(data_length),selected_yaw, 'k', label='yaw') 
plt.legend(loc='best')
fig.set_size_inches(30,20)
plt.xlabel('Sample number')
plt.ylabel('Angle [\u00b0]')
plt.savefig('Figures/KAREN/rpy_on_selection_file_4423', bbox_inches='tight', pad_inches=0)
plt.show()





from tfmra_py3_KAREN import tfmra_retracker
TrmfaThreshold=0.5
color_list=['powderblue','lightskyblue','deepskyblue','blue','darkblue','k','darkred','red','salmon','mistyrose']
#color_list=['k','k','k','k','k','k','k','k','k','k']
aval=[0.3,0.4,0.52,0.65,0.8,1,0.8,0.6,0.4,0.3]
BinWidthSnow=0.165*0.563

#################################################################################3
# plot  waveforms GPR_0
#################################################################################3
begin=round(data_length)-10
g=range(len(selected_Wf[0,:]))
for x in range(10):
    plt.figure()
    plt.plot(g, selected_Wf[begin+round(x),:].transpose())
    plt.title('Waveform' + str(x))
    plt.savefig('Figures/KAREN/GPR_0/GPR_0_wf_'+str(x)+'.png',dpi=200)        

gate = tfmra_retracker(selected_Wf[begin+5,:],TrmfaThreshold,0.5,40)[0]
print(gate)

sd=2.83
gate_end=gate + sd/BinWidthSnow
print(gate_end)
y=BinWidthSnow*np.linspace(0,(round(gate_end)-round(gate)),round(gate_end)-round(gate))

plt.figure()
for x in range(10):
    plt.plot(selected_Wf[begin+round(x),round(gate):round(gate_end)].transpose(),y, color=color_list[x] ,alpha=aval[x],label=str(x))
plt.ylabel('depth(m)')
plt.xlabel('amplitude')
plt.ylim([sd,0])
plt.savefig('Figures/KAREN/GPR_0/GPR_0_wfs', dpi=200)
plt.show()

#################################################################################3
# plot  waveforms GPR_150
#################################################################################3
begin=round((data_length/3)*2)-5
for x in range(10):
    plt.figure()
    plt.plot(g, selected_Wf[begin+round(x),:].transpose())
    plt.title('Waveform' + str(x))
    plt.savefig('Figures/KAREN/GPR_150/GPR_150_wf_'+str(x)+'.png',dpi=200)        
    
gate = tfmra_retracker(selected_Wf[begin+5,:],TrmfaThreshold,0.5,40)[0]

sd=2.82
gate_end=gate + sd/BinWidthSnow
print(gate_end)
y=BinWidthSnow*np.linspace(0,(round(gate_end)-round(gate)),(round(gate_end)-round(gate)))

print(gate)
plt.figure()
for x in range(10):
    plt.plot(selected_Wf[begin+round(x),round(gate):round(gate_end)].transpose(),y, color=color_list[x] ,alpha=aval[x],label=str(x))
plt.ylabel('depth (m)')
plt.xlabel('amplitude')
plt.ylim([sd,0])
plt.savefig('Figures/KAREN/GPR_150/GPR_150_wfs', dpi=200)
plt.show()


#################################################################################3
# plot  waveforms GPR_300
#################################################################################3
begin=round((data_length/3)*1)-5
for x in range(10):
    plt.figure()
    plt.plot(g, selected_Wf[begin+round(x),:].transpose())
    plt.title('Waveform' + str(x))
    plt.savefig('Figures/KAREN/GPR_300/GPR_300_wf_'+str(x)+'.png',dpi=200)        

gate = tfmra_retracker(selected_Wf[begin+5,:],TrmfaThreshold,0.5,40)[0]
sd=1.45
gate_end=gate + sd/BinWidthSnow
print(gate_end)
y=BinWidthSnow*np.linspace(0,(round(gate_end)-round(gate)),(round(gate_end)-round(gate)))

plt.figure()
for x in range(10):
    plt.plot(selected_Wf[begin+round(x),round(gate):round(gate_end)].transpose(),y, color=color_list[x] ,alpha=aval[x],label=str(x))
plt.ylabel('depth (m)')
plt.xlabel('amplitude')
plt.ylim([sd,0])
plt.savefig('Figures/KAREN/GPR_300/GPR_300_wfs', dpi=200)
plt.show()

#################################################################################3
# plot  waveforms GPR_343
#################################################################################3
begin=0
for x in range(10):
    plt.figure()
    plt.plot(g, selected_Wf[begin+round(x),:].transpose())
    plt.title('Waveform' + str(x))
    plt.savefig('Figures/KAREN/GPR_343/GPR_343_wf_'+str(x)+'.png',dpi=200)        

gate = tfmra_retracker(selected_Wf[begin+5,:],TrmfaThreshold,0.5,40)[0]
sd=1.23
gate_end=gate + sd/BinWidthSnow
print(gate_end)
y=BinWidthSnow*np.linspace(0,(round(gate_end)-round(gate)),round(gate_end)-round(gate))

plt.figure()
for x in range(10):
    plt.plot(selected_Wf[begin+round(x),round(gate):round(gate_end)].transpose(),y, color=color_list[x] ,alpha=aval[x],label=str(x))
plt.ylabel('depth (m)')
plt.xlabel('amplitude')
plt.ylim([sd,0])
plt.savefig('Figures/KAREN/GPR_343/GPR_343_wfs', dpi=200)
plt.show()



values_below = (selected_Roll < -1.5).sum()
print(values_below)
values_above=(selected_Roll > 1.5).sum()
print(values_above)
print('percent of discarded data on selection due to roll: ' + str((((values_below+values_above))/len(selected_Roll))*100) + '%' )


#Roll sorting
Roll_mask=(selected_Roll<=1.5) & (selected_Roll>=-1.5)
Roll_sorted_selected_Hocog_ESA = selected_Hocog_ESA[Roll_mask]
Roll_sorted_selected_Wf = selected_Wf[Roll_mask]
Roll_sorted_selected_Range=selected_Range[Roll_mask]
Roll_sorted_selected_Lon=selected_Lon[Roll_mask]
Roll_sorted_selected_Lat=selected_Lat[Roll_mask]
Roll_sorted_selected_FlightAltitude=selected_FlightAltitude[Roll_mask]

k=len(Roll_sorted_selected_Hocog_ESA)

# for x in range(0, k, 1):
#      plt.figure()
#      plt.plot(g, Roll_sorted_selected_Wf[round(x),:].transpose())
#      #plt.plot(g, gate_rtck)
#      plt.title('Waveform' + str(x))
#      plt.ylabel('Power [W]')
#      plt.savefig('Figures/KAREN/line_r waveforms/KAREN'+ str(x) +'.png',dpi=200)  




    
    
    