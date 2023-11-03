# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 22:45:17 2023

@author: Sofia Thirslund
"""

from ASIRAS_Importer_28 import *

#MAKING MAP

import matplotlib.path as mpath
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
plt.style.use('ggplot')

fig = plt.figure(figsize=[10, 5])
ax = fig.add_subplot(1, 2, 1, projection=ccrs.SouthPolarStereo())

fig.subplots_adjust(bottom=0.05, top=0.95,
                     left=0.04, right=0.95, wspace=0.02)


ax.gridlines()
 # Limit the map to -60 degrees latitude and below.
ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())

ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.RIVERS)

 # Compute a circle in axes coordinates, which we can use as a boundary
 # for the map. We can pan/zoom as much as we like - the boundary will be
 # permanently circular.
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)

ax.gridlines()

start_lon, start_lat = Lon[0], Lat[0]
end_lon, end_lat = Lon[len(Lon)-1],Lat[len(Lon)-1]

plt.plot(Lon,Lat,
         color='red',
         transform=ccrs.Geodetic(),
         )
plt.scatter([start_lon,end_lon],[start_lat,end_lat],marker='o',color='red',
            transform=ccrs.Geodetic())
plt.savefig('Figures/122/AntarcticaRouteSouthPolarStereo.png',bbox_inches='tight', dpi=300)

plt.show()

#Plotting Lon Lat
plt.plot(Lon,Lat)
plt.savefig('Figures/122/LonLatRoute.png',dpi=170)
plt.show()

#Plotting Lon lat with square map
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-80,20, -75, -60])
ax.stock_img()
ax.gridlines(draw_labels=True)
plt.plot(Lon,Lat,transform=ccrs.Geodetic())
plt.xlabel('Lon (deg)')
plt.ylabel('Lat (deg)')
plt.savefig('Figures/122/SquareMapRoutePlateCarre',dpi=180)
plt.show()


#Plotting Hwgs84 as a function of time
samples_var=np.linspace(1.0, Number_used, Number_used)
plt.plot(samples_var,Hwgs84)
plt.title('Height of airplane in WGS84')
plt.xlabel('sample point')
plt.ylabel('Hwgs84 (m)')
plt.savefig('Figures/122/HeightPlane.png',dpi=180)
plt.show()


#Plotting Hwgs as a colormap 
# Apply a fancy colormap to the figure
cmap = plt.get_cmap('plasma')
plt.set_cmap(cmap)
plt.scatter(Lon, Lat, c=Hwgs84)
cbar=plt.colorbar()
cbar.set_label('Hwgs84 (m)')
plt.title('Flown route')
plt.xlabel('Lon (deg)')
plt.ylabel('Lat (deg)')
plt.savefig('Figures/122/Hwgs as a colormap',bbox_inches="tight",dpi=200)
plt.show()

#Plotting route and selected route
data_start=9900
data_end=12604


selected_Lon=Lon[data_start:data_end]
selected_Lat=Lat[data_start:data_end]

plt.plot(Lon,Lat)
plt.plot(selected_Lon,selected_Lat, color = 'b', linewidth=3)
plt.title('Data Selection')
plt.xlabel('Lon (deg)')
plt.ylabel('Lat (deg)')
plt.savefig('Figures/122/LonLatRoutewSelected.png',dpi=200)
plt.show()

