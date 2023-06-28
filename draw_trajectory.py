from netCDF4 import Dataset
import netCDF4
import numpy as np
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import sys
from datetime import *
import cmaps
import xarray as xr

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

ncid = Dataset(filename,"r")
time = ncid.variables['time']
data = ncid.variables["points"][:]

ds = xr.open_dataset(ncid.src_file)
# vunits = ds.u.vert_units

horizontal = ncid.horizontal # if TRUE than it is 3D simulation
try:
    levels = ds.level
    levels_name = ds.level.long_name
    levels_units = ds.level.units
except:
    levels = data[0,:,2]
    levels_name = "model_levels"
    levels_units = "number"

    
level_min = np.min(levels)
level_max = np.max(levels)




time_convert = netCDF4.num2date(time[:], time.units, time.calendar)
time_ind = get_time_ind(time_convert,0)

npts = data.shape[1]
ntime = data.shape[0]
lons_all = data[:,:,0]
lats_all = data[:,:,1]

size = 2.
llcrnrlat = np.max([np.min(lats_all)-size,-90.])
urcrnrlat = np.min([np.max(lats_all)+size,+90.])
llcrnrlon = np.max([np.min(lons_all)-size,-180.])
urcrnrlon = np.min([np.max(lons_all)+size,180. ])


if(urcrnrlat <= 90):
    m = Basemap(projection='cyl',llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,
                llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon,resolution='l')
else:
    avglon = np.arange(lons_all.any())
    m = Basemap(projection='npstere',boundinglat=np.max([llcrnrlat,60.]),lon_0=avglon,resolution='c')


plt.figure(figsize=(6,4), dpi=150)

cs = m.drawcoastlines(linewidth=0.5)
cs = m.drawparallels(np.arange(-90., 91., 5), linewidth=0.1)
cs = m.drawmeridians(np.arange(-180,180., 5), linewidth=0.1)
cs = m.fillcontinents('grey', alpha = 0.2)

# boxstr = "start time: %s\nend time  :%s"%(ncid.start_time,ncid.end_time) 
# titlestr = "No of particles = %d\ndata = %s"%(npts,src_name) 

# if ncid.accuracy == 1: accuracy = "Accurate scheme"
# else: accuracy = "Coarse scheme"

# if ncid.horizontal == 1: horizontal = "2D"
# else: horizontal = "3D"

# titlestr = "Data: %s\nstart time: %s\nend time:  %s\n%s (%s), total # of particles = %d"%(src_name,ncid.start_time,ncid.end_time,accuracy,horizontal,npts) 
# plt.title(titlestr,loc='left', fontsize=6)

for ip in range(0,npts):
    lons = data[:,ip,0]
    lats = data[:,ip,1]
    hgts = data[:,ip,2]

    
    # time_ip = time[np.where( lons.mask )]
    # time_dur = (lons.shape[0]-np.sum(lons.mask))*dtime/24 
    # print( "Time integrated: ",ip, time_dur )
    
#     textstr = textstr+"\nip: %d"%ip+" duration (w): %6.4f"%(time_dur)
     
    color = "tab:red"
    linestyle='-'
#     if(hgts[0] <= 0.5): 
# #         color = "tab:red"
#         linestyle='-'
#     if(hgts[0] > 0.5 and hgts[0] <= 1 ): 
# #         color = "tab:blue"
#         linestyle='--'
#     if(hgts[0] > 1 and hgts[0] <= 3 ): 
# #         color = "tab:green"
#         linestyle=':'
#     if(hgts[0] > 3  ): 
# #         color = "tab:orange"
#         linestyle='-.'
    
    # label = "%d m (%4.1f days)"%(int(hgts[0]),time_dur)
    x,y = m(lons, lats)
    plt.plot(x,y,'-',  color=color, linestyle=linestyle, linewidth=0.7, alpha=0.5)
    cs = m.scatter(lons[time_ind], lats[time_ind], c=hgts[time_ind], 
                   vmin=level_min, vmax=level_max, cmap=cmaps.WhiteBlueGreenYellowRed, marker='s',
                   s=[10.], latlon=True)

    del lons
    del lats
    del hgts
    # del time_ip

ax = plt.gca()

# props = dict(boxstyle='round', facecolor='white', alpha=0.5) # wheat
# ax.text(0.05, 0.1, boxstr, transform=ax.transAxes, fontsize=6,
#         verticalalignment='top', bbox=props)

# ax.legend(loc='lower right',fontsize=4)

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color="black",linestyle='-' , lw=1),
                Line2D([0], [0], color="black",linestyle='--', lw=1),
                Line2D([0], [0], color="black",linestyle=':' , lw=1),
                Line2D([0], [0], color="black",linestyle='-.', lw=1)]

# ax.legend(custom_lines, ['0.3 km', '1 km ', '3 km ', '5 km'],
#           loc='lower right',
#           title = 'start height',
#           title_fontsize = 8,
#           fontsize=8)
# # ax.legend(custom_lines, ['0.3 km', '1 km ', '3km '],loc='upper right',fontsize=8)

ax = plt.gca()
im = ax.imshow(np.arange(100).reshape((10,10)), vmin=level_min, vmax=level_max, cmap=cmaps.WhiteBlueGreenYellowRed)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax, label=f"{levels_name} ({levels_units})")
cbar.ax.tick_params(labelsize=5) 
plt.show()

# figname = "%s.png"%(filename)
# plt.savefig(figname)
# plt.close()

# props = dict(boxstyle='round', facecolor='white', alpha=0.5) # wheat
# ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=6,
#         verticalalignment='top', bbox=props)
