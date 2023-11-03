# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 12:16:09 2023

@author: Sofia Thirslund
"""

from tfmra_py3 import tfmra_retracker
import sys;
import struct
import numpy as np
import math
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


ifile = "ALS/ALS_20171228_1844_1x1.scn" 


f = open(ifile, 'r')
data = np.genfromtxt(f, dtype=float)
Lat=data[:,1]
Lon=data[:,2]
ALS=data[:,3]

data_start=0
data_end=len(Lat)-1
data_length=data_end-data_start


selected_Lon=Lon[data_start:data_end]
selected_Lat=Lat[data_start:data_end]
selected_ALS=ALS[data_start:data_end]


#Sorting for outliers

#ALS_mask=(selected_ALS<40) & (selected_ALS>-50)

#selected_ALS=selected_ALS[ALS_mask]
#selected_Lon=selected_Lon[ALS_mask]
#selected_Lat=selected_Lat[ALS_mask]


# print(selected_Lat)
# Lat_min=min(selected_Lat)
# Lat_max=max(selected_Lat)
# plt.plot(selected_Lon,selected_Lat)
# plt.ylim(Lat_min,Lat_max)
# plt.show()



# print(min(Lon))
# print(max(Lon))
# plt.plot(Lon,Lat)
# plt.xlim(min(Lon),max(Lon))


# #Plot data and data selection and wf 615
# plt.figure()
# #cmap = plt.get_cmap('plasma')
# #plt.set_cmap(cmap)
# plt.scatter(Lon,Lat)
# plt.scatter(selected_Lon, selected_Lat, c='b')
# #cbar=plt.colorbar()
# #cbar.set_label('Elevation (m)')
# plt.scatter(Lon_doublepeak,Lat_doublepeak, color = 'k')
# plt.title('ALS data selection')
# plt.savefig('Figures/ALS/ALS datapoints with doublepeak as black')
# plt.show()







#Find nearest value
from scipy import spatial

selected_Lon = np.expand_dims(selected_Lon, axis=0)
selected_Lat = np.expand_dims(selected_Lat, axis=0)

selected_coordinates=np.concatenate((selected_Lon.T,selected_Lat.T),axis=1)

# dp=[Lon_doublepeak,Lat_doublepeak]

# distance,index = spatial.KDTree(selected_coordinates).query(dp)
# print(distance)
# print(index)

# dp_closest_coord=[selected_Lon[0,index],selected_Lat[0,index]]
# print(dp_closest_coord) #this is very close to dp!
# print(dp)

# ALS_dp=selected_ALS[index]
# print(ALS_dp)



#Extract ALS along b line
f = open("ALS/b_line_coord_Lon_KAREN" , 'r')
b_line_Lon= np.genfromtxt(f, dtype=float)

f = open("ALS/b_line_coord_Lat_KAREN" , 'r')
b_line_Lat=np.genfromtxt(f, dtype=float)

b_line_Lon= np.expand_dims(b_line_Lon, axis=0)
b_line_Lat = np.expand_dims(b_line_Lat, axis=0)
b_line_coordinates=np.concatenate((b_line_Lon.T,b_line_Lat.T),axis=1)

print(np.shape(b_line_Lat)[1])
b_line_closest_coord=np.zeros((np.shape(b_line_Lat)[1],2))
distance=np.zeros((np.shape(b_line_Lat)[1]))
indexes=np.zeros((np.shape(b_line_Lat)[1]))

for ii in range(np.shape(b_line_Lat)[1]):

    distance[ii],index = spatial.KDTree(selected_coordinates).query(b_line_coordinates[ii])
    print(index)
    b_line_closest_coord[ii]=[selected_Lon[0,index],selected_Lat[0,index]]
    indexes[ii]=index

    
indexes=indexes.astype(int)
b_line_ALS=selected_ALS[indexes]
    

plt.figure()
plt.plot(b_line_closest_coord[:,0],b_line_closest_coord[:,1],color='g',label='ALS samples on b line')
plt.plot(b_line_Lon[0,:],b_line_Lat[0,:],color='r',label='KAREN samples on b line')
plt.xlabel('Lat')
plt.ylabel('Lon')
plt.legend()
plt.savefig('Figures/ALS/b line coordinates ALS vs KAREN'+'.png',dpi=200)
plt.show()

print(b_line_ALS)

# plt.plot(Lon,Lat)

# # #Plot selected data and line b
# from matplotlib.ticker import ScalarFormatter

# class ScalarFormatterClass(ScalarFormatter):
#     def _set_format(self):
#        self.format = "%1.2f"

# fig, ax = plt.subplots()
# ax.set_ylabel('Lat')
# ax.set_xlabel('Lon')
# cmap = plt.get_cmap('plasma')
# plt.set_cmap(cmap)
# plt.scatter(selected_Lon, selected_Lat, c=selected_ALS)
# cbar=plt.colorbar()
# cbar.set_label('Elevation (m)')
# #plt.scatter(b_line_closest_coord[:,0],b_line_closest_coord[:,1],color='g',label='ALS coord of b line')
# plt.title('Elevation (ALS)')
# plt.savefig('Figures/ALS/ALS elevation_29_2800')
# plt.show()


#Extract lead value
# Lon_lead=-52.2010797
# Lat_lead=-76.73824669
# lead=[Lon_lead,Lat_lead]


# distance,index = spatial.KDTree(selected_coordinates).query(lead)
# print(distance)
# print('index = ' +str(data_start+index))

# print('elevation ALS: ' +str(ALS[data_start+index]))
# print('Lon of lead: ' + str(Lon[data_start+index]))
# print('Lat of lead:' +str(Lat[data_start+index]))

