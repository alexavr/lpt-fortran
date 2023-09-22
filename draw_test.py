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
from pathlib import Path
import wrf # conda install -c conda-forge wrf-python (https://anaconda.org/conda-forge/wrf-python)
from netCDF4 import Dataset


itime = 5
iz = 5

# filename = "../data/em_tropical_cyclone_15km.nc"
# f = Dataset(filename)
# u = wrf.getvar(f, "ua",itime)[iz,:,:]
# v = wrf.getvar(f, "va",itime)[iz,:,:]
# nx = wrf.extract_dim(f, "west_east")
# ny = wrf.extract_dim(f, "south_north")

filename = "em_tropical_cyclone_15km_pl.nc"
ds = xr.open_dataset(filename)
u = ds.u[itime,iz,:,:]
v = ds.v[itime,iz,:,:]
nx = u.shape[1]
ny = u.shape[0]

x = np.arange(0, nx)
y = np.arange(0, ny)
X, Y = np.meshgrid(x, y)

wndspd = np.sqrt(u**2+v**2)

print(u.shape)
print(v.shape)
print(X.shape)
print(Y.shape)

plt.figure(figsize=(5,5), dpi=150)
plt.contourf(wndspd, alpha=0.8)
plt.streamplot(X,Y,u,v,color='black',linewidth=0.5)

plt.show()