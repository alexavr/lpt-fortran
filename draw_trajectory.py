from netCDF4 import Dataset
import netCDF4
import numpy as np
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import sys
from datetime import *


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

time_convert = netCDF4.num2date(time[:], time.units, time.calendar)
time_ind = get_time_ind(time_convert,0)

npts = data.shape[1]
ntime = data.shape[0]
lons_all = data[:,:,0]
lats_all = data[:,:,1]

size = 2.
llcrnrlat = np.min(lats_all)-size # np.max([plat-size,-90.])
urcrnrlat = np.max(lats_all)+size # np.min([plat+size,+90.])
llcrnrlon = np.min(lons_all)-size # np.max([plon-size,0.  ])
urcrnrlon = np.max(lons_all)+size # np.min([plon+size,360.])

m = Basemap(projection='cyl',llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,
            llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon,resolution='l')
# m = Basemap(projection='cyl',llcrnrlat=10,urcrnrlat=90,
#             llcrnrlon=-100,urcrnrlon=-10,resolution='l')
# m = Basemap(projection='npstere',boundinglat=60,lon_0=-50,resolution='l')

plt.figure(figsize=(6,4), dpi=150)

cs = m.drawcoastlines(linewidth=0.5)
cs = m.drawparallels(np.arange(-90.,91.,1.),linewidth=0.1)
cs = m.drawmeridians(np.arange(-180,180.,1.),linewidth=0.1)
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
    hgts = data[:,ip,2]/1000
    
    # time_ip = time[np.where( lons.mask )]
    # time_dur = (lons.shape[0]-np.sum(lons.mask))*dtime/24 
    # print( "Time integrated: ",ip, time_dur )
    
#     textstr = textstr+"\nip: %d"%ip+" duration (w): %6.4f"%(time_dur)
     
    color = "red"
    linestyle=':'
    if(hgts[0] <= 0.5): 
#         color = "tab:red"
        linestyle='-'
    if(hgts[0] > 0.5 and hgts[0] <= 1 ): 
#         color = "tab:blue"
        linestyle='--'
    if(hgts[0] > 1 and hgts[0] <= 3 ): 
#         color = "tab:green"
        linestyle=':'
    if(hgts[0] > 3  ): 
#         color = "tab:orange"
        linestyle='-.'
    
    # label = "%d m (%4.1f days)"%(int(hgts[0]),time_dur)
    plt.plot(lons, lats,'-',  color=color, linestyle=linestyle, linewidth=1)
    cs = m.scatter(lons[time_ind], lats[time_ind], c=hgts[time_ind], 
                   vmin=0, vmax=6, cmap="jet", marker='s', 
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

ax.legend(custom_lines, ['0.3 km', '1 km ', '3 km ', '5 km'],
          loc='lower right',
          title = 'start height',
          title_fontsize = 8,
          fontsize=8)
# ax.legend(custom_lines, ['0.3 km', '1 km ', '3km '],loc='upper right',fontsize=8)

ax = plt.gca()
im = ax.imshow(np.arange(100).reshape((10,10)), vmin=0, vmax=6, cmap="jet")
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.show()

# props = dict(boxstyle='round', facecolor='white', alpha=0.5) # wheat
# ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=6,
#         verticalalignment='top', bbox=props)