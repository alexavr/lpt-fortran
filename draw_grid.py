#!/usr/bin/env python3
# 
import numpy as np # conda install -c anaconda numpy
from mpl_toolkits.basemap import Basemap # conda install -c conda-forge matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import sys
import cmaps  # conda install -c conda-forge cmaps
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.ndimage import gaussian_filter
import pandas as pd 

### CHANGE THIS ################################################################

tracking = False    # draw recent path on every particle
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

ds_trk = xr.open_dataset(filename)
ds_src = xr.open_dataset(ds_trk.src_file)

time_trk = [ pd.to_datetime(str(i.values)) for i in ds_trk.time]
time_src = [ pd.to_datetime(str(i.values)) for i in ds_src.time] 

data  = ds_trk.points

try:
    lon2d = ds_src.XLONG
    lat2d = ds_src.XLAT
except:
    lat1d = ds_src.latitude
    lon1d = ds_src.longitude
    lon2d, lat2d = np.meshgrid(lon1d,lat1d)
    del lat1d, lon1d


# get N timesteps in daystep (1 day)
if tracking: 
    daystep = int(ltrajectory*3600./(time_trk[1] - time_trk[0]).total_seconds())

ntime = len(time_trk)
npts  = len(data[0,:,0])

# Получаем уровни по вертикали.
# Если плоский запуск, то читаем вектор высот (для отрисовки) и актуальную высоту.
if ds_trk.horizontal:
    levels = ds_src.level
    level_real = ds_trk.horizontal_level
    level = np.argmin(np.abs(levels-level_real))
else:
    levels = data[0,:,2]
    level = np.abs(ds_src.z[0,:,:,:].mean(axis=(1,2)) - levels.mean() ).argmin()
    level_real = levels[level]



# For colorbar
# level_real = level_real/100
levs = np.arange(np.floor(ds_src.z.isel(level=level).min()/100), np.ceil(ds_src.z.isel(level=level).max()/100), 4)

# Main loop
for it in range(0,ntime): # ntime
 

    lons = data[it,:,0]
    lats = data[it,:,1]
    hgts = data[it,:,2]

    try:
        itt = time_src.index(time_trk[it])
        z = ds_src.z.isel(time=itt,level=level)/100
        smooth_z = gaussian_filter(z, sigma=1)
        del z
    except:
        itt = None

    plt.figure(figsize=(6,4), dpi=150)

    proj = ccrs.PlateCarree()

    # ax = plt.axes(projection=ccrs.LambertConformal(
    #     central_longitude=-45.0,central_latitude=50.0))

    ax = plt.axes(projection=ccrs.NearsidePerspective(
        central_longitude=-45, central_latitude=45, 
        satellite_height=39785831, 
        false_easting=0, false_northing=0, globe=None))
    


    ax.coastlines('110m',linewidth=0.2, alpha=1, color="black")
    ax.add_feature(cfeature.LAND, alpha=0.5)
    ax.add_feature(cfeature.OCEAN, alpha=0.2)

    gl = ax.gridlines(crs=proj,
                        # draw_labels=True, dms=True, 
                        # x_inline=False, y_inline=False,  
                        linewidth=0.5, linestyle=":", color='gray', alpha=0.5)
    gl.xlabels_top = False
    gl.ylabels_right = False
    # gl.xlines = False
    # gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
    # gl.xformatter = LONGITUDE_FORMATTER
    # gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 6} # , 'color': 'gray'
    gl.ylabel_style = {'size': 6} # , 'color': 'gray'
    gl.rotate_labels=0

    ax.set_extent([-85, -5, 5, 80], crs=ccrs.PlateCarree())

    # Рисуем поля геопотенциала
    zfil = plt.contourf(ds_src.XLONG, ds_src.XLAT, smooth_z, alpha=0.5, levels=levs, cmap="YlGn_r", transform=proj)
    zcnt = plt.contour(ds_src.XLONG, ds_src.XLAT, smooth_z, levels=levs, linestyles='-', colors="g", linewidths=0.1, transform=proj)
    # plt.contourf(smooth_z,transform=proj)
    
    cbar = plt.colorbar(zfil, shrink=0.95, pad=0.01)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label('Geopotential height [hPa]', labelpad=-26, fontsize=6, y=0.5, rotation=90)

    # Разбрасываем точки на текущий момент
    if ds_trk.horizontal:
        cs = plt.scatter(lons, lats, color="black", s=[0.6], transform=proj)
    else:
        cs = plt.scatter(lons, lats, c=hgts, 
                       vmin=level_min, vmax=level_max, cmap=cmaps.BlueYellowRed, s=[0.6], transform=proj)

    # Рисуем хвоста за последние ltrajectory шагов
    if tracking:
        ts = it - np.min([daystep,it])
        for ip in range(0,npts):
            plt.plot(data[ts:it,ip,0], data[ts:it,ip,1],'-', color="black", alpha=0.9, linewidth=0.2, transform=proj)
            # plt.text(data[it,   ip,0], data[   it,ip,1], ip, fontsize=5, transform=proj)

    # print(time_trk[it])
    titlestrL = f"{time_trk[it]}"
    titlestrR = f"{(float(lons.count())/float(lons.shape[0])*100.):4.1f}%"
    plt.title(titlestrL,loc='left', fontsize=6, y=0.98)
    plt.title(titlestrR,loc='right', fontsize=6, y=0.98)

    print(f"{time_trk[it]}   {titlestrR}", end='\r', flush=True)

    # plt.tight_layout()
    figname = f"./pics/grid_{filename}_{it:07d}.png"
    # plt.show()
    plt.savefig(figname, dpi=150)
    plt.close()

print(f"Create video:")
print(f"ffmpeg -framerate 10 -i grid_{filename}_%07d.png \\"        )
print(f"   -c:v libx264 -r 30 -pix_fmt yuv420p grid_{filename}.mp4"  )
print(f"or gif animation:"                                          )
print(f"convert -delay 20 grid_{filename}_00*.png grid_{filename}.gif")
