#!/usr/bin/env python3
# 
from netCDF4 import Dataset
import netCDF4
import numpy as np
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from wrf import smooth2d
import sys

### CHANGE THIS ################################################################

tracking = True    # draw recent path on every particle
if tracking: 
    ltrajectory = 24 # the path length (in hours)

################################################################################
################################################################################
################################################################################

### FUNCTIONS ################################################################

def get_z_ind( ncid, zarray ):

    zavg = np.average(zarray)

    z = np.average(ncidz.variables["z"][0,:,:,:],axis=(1,2))/9.81

    idx = np.abs(z - zavg).argmin()

    return idx

################################################################################
################################################################################
################################################################################

filename = sys.argv[1]

ncidp = Dataset(filename,"r")
time  = ncidp.variables['time']
data  = ncidp.variables["points"][:]

pt_file = ncidp.out_file

zfilename = ncidp.pwd+ncidp.src_file
ncidz = Dataset(zfilename,"r")

try:
    lon2d = ncidz.variables['XLONG'][:]
    lat2d = ncidz.variables['XLAT'][:]
except:
    lat1d = ncidz.variables["latitude"][:]
    lon1d = ncidz.variables["longitude"][:]
    lon2d, lat2d = np.meshgrid(lon1d,lat1d)
    del lat1d, lon1d

ztime = ncidz.variables['time']
ptime = ncidp.variables['time']
ztime_convert = netCDF4.num2date(ztime[:], ztime.units, ztime.calendar)
ptime_convert = netCDF4.num2date(ptime[:], ptime.units, ptime.calendar)

if tracking:
    daystep = int((ltrajectory*3600.)/(ptime_convert[1]-ptime_convert[0]).total_seconds())

data  = ncidp.variables["points"][:]
ntime = data.shape[0]
npts  = data.shape[1]

level = get_z_ind( ncidz, data[0,:,2] )
print("determined vertical level = ",level)

size = 2.
llcrnrlat = np.max([np.min(data[:,:,1])-size,-90.])
urcrnrlat = np.min([np.max(data[:,:,1])+size,+90.])
llcrnrlon = np.max([np.min(data[:,:,2])-size,-180.])
urcrnrlon = np.min([np.max(data[:,:,2])+size,180. ])

# if(urcrnrlat <= 70):
#     m = Basemap(projection='cyl',llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,
#                 llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon,resolution='l')
# else:
#     avglon = np.arange(data[:,:,2].any())
#     m = Basemap(projection='npstere',boundinglat=np.max([np.min(data[:,:,1])-5,60.]),lon_0=avglon,resolution='l')

# Define the NAAD map:
m = Basemap(width=10000000,height=10000000,
            rsphere=(6378137.00,6356752.3142),\
            resolution='c',area_thresh=1000.,projection='lcc',\
            lat_1=45.,lat_2=45,lat_0=45,lon_0=-48.)

if ncidp.accuracy == 1: accuracy = "Accurate scheme"
else: accuracy = "Coarse scheme"

if ncidp.horizontal == 1: horizontal = "2D"
else: horizontal = "3D"

for it in range(0,ntime): # ntime

    print(ptime_convert[it], end='\r', flush=True)

    lons = data[it,:,0]
    lats = data[it,:,1]
    hgts = data[it,:,2]/1000
    
    itt = np.where(ztime_convert == ptime_convert[it])[0].tolist()
    if itt != []: 
        z = ncidz.variables["z"][itt,level,:,:]
        smooth_z = smooth2d(z[0,:,:], 6, cenweight=4)
        del z

    fig = plt.figure(figsize=(6,4), dpi=300)
    plt.rcParams.update({'font.size': 6}) 
    
    m.drawcoastlines(linewidth=0.1)
    m.drawparallels(np.arange(-90.,91.,5.) ,linestyle='--',linewidth=0.1)
    m.drawmeridians(np.arange(-180,180.,5.),linestyle='--',linewidth=0.1)
    # m.drawparallels(np.arange(0,80,15),labels=[False,True,False,False])
    m.drawmeridians(np.arange(-180,180.,5.),linestyle='--',linewidth=0.1)
    m.fillcontinents('tab:gray', alpha = 0.7)

    titlestrL = "%s"%(ptime_convert[it])
    titlestrR = "%4.1f%%"%((float(lons.count())/float(lons.shape[0])*100.))
    
    plt.title(titlestrL,loc='left', fontsize=6, y=0.98)
    plt.title(titlestrR,loc='right', fontsize=6, y=0.98)

    if tracking:
        cs = m.contour(lon2d, lat2d, smooth_z, levels=np.arange(0, 100001, 100),
            linestyles='-', colors="g", linewidths=0.1, latlon=True)
    else:
        cs = m.contour(lon2d, lat2d, smooth_z, levels=np.arange(0, 100001, 100),
            linestyles='-', colors="black", linewidths=0.5, latlon=True)
    
    cs = m.contourf(lon2d, lat2d, smooth_z, 10, # levels=np.arange(6000, 11001, 500), #levels=np.arange(6400, 10501, 100),
        alpha = 0.5, extend='min', cmap="YlGn_r", latlon=True)
    
    if tracking:
        ts = it - np.min([daystep,it])
        for ip in range(0,npts):
            x,y = m(data[ts:it,ip,0], data[ts:it,ip,1])
            plt.plot(x,y,'-', 
                color="black", linewidth=0.5)

    cs = m.scatter(lons, lats, c=hgts, 
                   vmin=0, vmax=4, cmap="jet",s=[0.6], latlon=True)

    clb = fig.colorbar(cs)
    clb.set_label('Particle height [km]', labelpad=-33, fontsize=6, y=0.5, rotation=90)

    figname = "grid_%s_%07d.png"%(pt_file,it)
    # plt.show()
    fig.savefig(figname)
    plt.close()

    del lons
    del lats
    del hgts
print("\n DONE!")

print("Create video:")
print("ffmpeg -framerate 10 -i grid_%s_%%07d.png \\"%(pt_file))
print("   -c:v libx264 -r 30 -pix_fmt yuv420p grid_%s.mp4"%(pt_file))
print("or gif animation:")
print("convert -delay 20 grid_%s_00*.png grid_%s.gif"%(pt_file,pt_file))
