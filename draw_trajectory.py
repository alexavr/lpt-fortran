import numpy as np                  # conda install -y -c anaconda numpy
import matplotlib.pyplot as plt     # conda install -y -c conda-forge matplotlib
from netCDF4 import Dataset
import xarray as xr                 # conda install -y xarray dask netCDF4 bottleneck
import cartopy.crs as ccrs          # conda install -y -c conda-forge cartopy
import cartopy.feature as cfeature
import cmaps                        # conda install -y -c conda-forge cmaps
import sys

################################################################################

def get_time_ind( time, time_step_hours ):

    HH = np.asarray([t.hour   for t in time[:]])
    MM = np.asarray([t.minute for t in time[:]])
    SS = np.asarray([t.second for t in time[:]])

    time_ind = np.where( (HH==time_step_hours) & (MM == 0) & (SS==0) )

    return time_ind

################################################################################
################################################################################
################################################################################

filename = sys.argv[1]

ds_trk = xr.open_dataset(filename)
ds_src = xr.open_dataset(ds_trk.src_file)

horizontal = ds_trk.horizontal # TRUE if it's 2D simulation

levels = ds_src.level

level_min = levels.min()
level_max = levels.max()

npts = len(ds_trk.n)
ntime = len(ds_trk.time)

plt.figure(figsize=(6,4), dpi=150)

proj = ccrs.PlateCarree()

ax = plt.axes(projection=proj)
ax.coastlines('110m',linewidth=0.2, alpha=1, color="black")
ax.add_feature(cfeature.LAND, alpha=0.5)
ax.add_feature(cfeature.OCEAN, alpha=0.2)

gl = ax.gridlines(crs=proj,
                  draw_labels=True, 
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




for ip in range(0,npts):
    lons = ds_trk.points[:,ip,0]
    lats = ds_trk.points[:,ip,1]
    hgts = ds_trk.points[:,ip,2]

    plt.plot(lons, lats, alpha=1, transform=proj, color="tab:red")
    # plt.scatter(lons, lats, alpha=1, transform=proj, s=10, marker='x', color="black", zorder=10)
    plt.scatter(lons[::8], lats[::8], alpha=1, transform=proj, s=10, marker='x',  color="black", zorder=10)

plt.show()
