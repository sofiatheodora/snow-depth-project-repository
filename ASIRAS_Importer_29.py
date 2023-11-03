# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 22:16:40 2023

@author: Sofia Thirslund
Based on: "ReadAsirasL1B_trmfa_py3" by Henriette Skourup, DTU Space, hsk@space.dtu.dk, created June 2018
"""

import sys;
import struct
import numpy as np


idir = 'ASIRAS/'
ifile = 'AS6OA00_ASIWL1B040320171229T124148_20171229T135011_0001'
TotSize = 21095319 # taken manually from ASCII header
SPHsize = 2512  # might be variable, taken manually from ASCII header
TimeStart = 0. 
TimeEnd = 86400.

Retracker = 'TRMFA50'

TrmfaThreshold = 0.01

#TAI to UTC offset
# 2017 37s
TAI2UTCoffset=37.

# Waveform parameters
BinWidth=0.1098
BinRef=128
WvfBinNumber = 256

    
ASCIIheaderSize = 891
MPHsize = 1247  
#SPHsize = 2512 # might be variable, taken manually from ASCII header
SPHsize = int(SPHsize) # might be variable, taken manually from ASCII header

HeaderSize = MPHsize + SPHsize
print(HeaderSize)


DataRecSize = 16660     
WvfGroupSize = 624     
GroupRepeatSize = 20
TotSize = int(TotSize) # taken manually from ASCII header

NumRec = int((TotSize - HeaderSize)/DataRecSize)

endian  = 3210

print('Number of records:',NumRec)

print("ASCII header size %s" % ASCIIheaderSize)
print("MPH header size %s" % MPHsize)
print("SPH header size %s" % SPHsize)
print("Header size total %s" % HeaderSize)
print("# of records %s" % NumRec)

#print "# records %s" % numrec
#print "min_lat, max_lat, min_lon, max_lon %s" % min_lat, max_lat, min_lon, max_lon

    
#Open datafile for later use
f = open(idir+ifile + ".DBL","rb") 

print(ifile+'.DBL')

f.seek(HeaderSize)


Number_used=int(GroupRepeatSize*NumRec)
print(GroupRepeatSize*NumRec)

lat = np.zeros(Number_used) * np.nan 
lon = np.zeros(Number_used) * np.nan 
TAI = np.zeros(Number_used) * np.nan 
Time = np.zeros(Number_used) * np.nan 
#Tsec = np.zeros(GroupRepeatSize*NumRec) * np.nan 
#Tmicsec = np.zeros(GroupRepeatSize*NumRec) * np.nan 
Lat = np.zeros(Number_used) * np.nan 
Lon = np.zeros(Number_used) * np.nan 
Hwgs84 = np.zeros(Number_used) * np.nan 
WindowDelay = np.zeros(Number_used) * np.nan 
Hocog = np.zeros(Number_used) * np.nan 
IRCorr = np.zeros(Number_used) * np.nan 
Roll = np.zeros(Number_used) * np.nan
Hftrma = np.zeros(Number_used) * np.nan
ocog  = np.zeros(Number_used) * np.nan
ocogRT  = np.zeros(Number_used) * np.nan
ocogWidth  = np.zeros(Number_used) * np.nan
Pmax  = np.zeros(Number_used) * np.nan
iPmax  = np.zeros(Number_used) * np.nan
print(np.shape(lat))


wf = np.zeros((Number_used,WvfBinNumber)) * np.nan 
A = np.zeros(Number_used) * np.nan 
B = np.zeros(Number_used) * np.nan 

PWF = np.zeros((Number_used,WvfBinNumber)) * np.nan 


print(np.shape(PWF))
print(NumRec)

# loop over several records 
nloop=0
oloop=0
mloop=0

for ii in range(1,NumRec+1):
#for ii in range(1,3):
    # Read Time Orbit Group
#    print 'Read Time Orbit Group'
    for n in range(1,GroupRepeatSize+1):
        byteTO = f.read(84)
##        frb=[]
#
        TimeOrbitData =  struct.unpack('>l L L l H H L L l l l l 3l 3l 3l L',byteTO)

#        print nloop
#        print TimeOrbitData[8]*10E-8,TimeOrbitData[9]*10E-8
#        print TimeOrbitData[0],TimeOrbitData[1],TimeOrbitData[2]*1E-3,TimeOrbitData[10]*1E-3
#        print >> output, "%12.6f %13.6f" % (TimeOrbitData[8]*1E-7,TimeOrbitData[9]*1E-7)
        TAI[nloop]=TimeOrbitData[0]
#        Tsec[nloop]=TimeOrbitData[1]
#        Tmicsec[nloop]=TimeOrbitData[2]*10E-7
        Time[nloop]=TimeOrbitData[1]+TimeOrbitData[2]*1E-6-TAI2UTCoffset
        Lat[nloop]=TimeOrbitData[8]*1E-7
        Lon[nloop]=TimeOrbitData[9]*1E-7
        Hwgs84[nloop]=TimeOrbitData[10]*1E-3
        
        nloop=nloop+1

## working :-)
#    # Read Measurement Group
# #   print 'Read Measurement Group'
##    print struct.calcsize('=q 17l 4h l 3H')
#    f.read(94*GroupRepeatSize)    

    for o in range(1,GroupRepeatSize+1):
        byteMG = f.read(94)    
        MGData = struct.unpack('>q 17l 4h l 3H',byteMG)
#   Window-delay, Hocog, Roll
#        print MGData[0]*1E-12,MGData[4]*1E-3,MGData[18]*1E-3
        WindowDelay[oloop]=MGData[0]*1E-12
        Hocog[oloop]=MGData[4]*1E-3
        IRCorr[oloop]=MGData[12]*1E-3
        Roll[oloop]=MGData[18]*1E-3
        oloop = oloop + 1

    # Read Correction Group
    # Zeroed for ASIRAS
    f.read(64)

    # Read Average Waveform Group
    # Zeroed for ASIRAS
    f.read(556)

    # Read Waveform Group


#    print 'Read waveform Group'
    #print struct.calcsize('=256H l l H H 50H')

    for m in range(1,GroupRepeatSize+1):
        byteWvf=f.read(624)
        mmloop=0
        WvfGroupData =  struct.unpack('>256H l l H H 50H',byteWvf)
        for mm in range(1,WvfBinNumber+1):
#            WvfGroupData[nloop,mm]=WvfGroupData[mm]

#             wf[mloop,mmloop]=WvfGroupData[mmloop]

             PWF[mloop,mmloop]=1E-9*(2**WvfGroupData[257])*WvfGroupData[256]*WvfGroupData[mmloop]
            
#             print mloop,mmloop,WvfGroupData[mmloop]
             mmloop=mmloop+1
#        A[mloop]=WvfGroupData[256]
#        B[mloop]=WvfGroupData[257]

        mloop=mloop+1