#!/usr/bin/env python3
# 
from netCDF4 import Dataset
import netCDF4
import numpy as np # conda install -c anaconda numpy
from mpl_toolkits.basemap import Basemap # conda install -c conda-forge matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from wrf import smooth2d
import sys
import cmaps  # conda install -c conda-forge cmaps

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

    z = np.average(ncid.variables["z"][0,:,:,:],axis=(1,2))/9.81

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

horizontal = ncidp.horizontal # if TRUE than it is 3D simulation
if horizontal:
    levels = ncidz.variables['level']
    level_real = ncidp.getncattr('horizontal_level')
    level = np.argmin(np.abs(levels-level_real))
else:
    levels = data[0,:,2]
    level = get_z_ind( ncidz, levels )

level_min = np.min(levels)
level_max = np.max(levels)

print(f"Vertical level index {level} and value {levels[level]:.1f} for plot pressure contours")

size = 2.
llcrnrlat = np.max([np.min(data[:,:,0])-size,-90.])
urcrnrlat = np.min([np.max(data[:,:,0])+size,+90.])
llcrnrlon = np.max([np.min(data[:,:,1])-size,-180.])
urcrnrlon = np.min([np.max(data[:,:,1])+size,180. ])

# if(urcrnrlat <= 70):
#     m = Basemap(projection='cyl',llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,
#                 llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon,resolution='l')
# else:
#     avglon = np.arange(data[:,:,2].any())
#     m = Basemap(projection='npstere',boundinglat=np.max([np.min(data[:,:,1])-5,60.]),lon_0=avglon,resolution='l')

# Define the NAAD map:

m = Basemap(width=8000000,height=8500000,
            rsphere=(6378137.00,6356752.3142),\
            resolution='c',area_thresh=1000.,projection='lcc',\
            lat_1=45.,lat_2=45,lat_0=45,lon_0=-45.)

if ncidp.accuracy == 1: accuracy = "Accurate scheme"
else: accuracy = "Coarse scheme"

for it in range(0,ntime): # ntime

    print(ptime_convert[it], end='\r', flush=True)
    # print(ptime_convert[it])

    lons = data[it,:,0]
    lats = data[it,:,1]
    hgts = data[it,:,2]
    
    itt = np.where(ztime_convert == ptime_convert[it])[0].tolist()
    if itt != []: 
        z = ncidz.variables["z"][itt,level,:,:]/100
        # levs = np.arange(np.floor(np.quantile(z,0.05)), np.ceil(np.quantile(z,0.95)), 2)
        levs = np.arange(np.floor(np.min(z)), np.ceil(np.max(z)), 2)
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
        cs = m.contour(lon2d, lat2d, smooth_z, levels=levs,
            linestyles='-', colors="g", linewidths=0.1, latlon=True)
    else:
        cs = m.contour(lon2d, lat2d, smooth_z, levels=levs,
            linestyles='-', colors="black", linewidths=0.5, latlon=True)
    
    cs = m.contourf(lon2d, lat2d, smooth_z, levels=levs, #levels=np.arange(6400, 10501, 100),
        alpha = 0.5, extend='min', cmap="YlGn_r", latlon=True)
    
    if tracking:
        ts = it - np.min([daystep,it])
        for ip in range(0,npts):
            x,y = m(data[ts:it,ip,0], data[ts:it,ip,1])
            plt.plot(x,y,'-', 
                color="black", alpha=0.5, linewidth=0.2)

    cs = m.scatter(lons, lats, c=hgts, 
                   vmin=level_min, vmax=level_max, cmap=cmaps.BlueYellowRed, s=[0.6], latlon=True)

    clb = fig.colorbar(cs)
    clb.set_label('Particle height [km]', labelpad=-33, fontsize=5, y=0.5, rotation=90)

    figname = f"grid_{filename}_{it:07d}.png"
    # plt.show()
    fig.savefig(figname)
    plt.close()

    del lons
    del lats
    del hgts
print("\n DONE!")

print(f"Create video:")
print(f"ffmpeg -framerate 10 -i grid_{filename}_%07d.png \\"        )
print(f"   -c:v libx264 -r 30 -pix_fmt yuv420p grid_{filename}.mp4"  )
print(f"or gif animation:"                                          )
print(f"convert -delay 20 grid_{filename}_00*.png grid_{filename}.gif")
