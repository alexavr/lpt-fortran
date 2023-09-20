#!/usr/bin/env python3

# conda activate data
# python prep_wrf_pres.py {filename}.nc

# conda install -c conda-forge wrf-python
# 1. ln -sf /storage/NAD/NAAD/v4/LoRes/2015/wrfout_d01_2015-01-* .
# 2. for i in $(ls wrfout_d01_*); do ./prep_wrf_isobaric.py $i; done
# 3. cdo mergetime _pt_wrfout_d01_* _WRF.nc
# 4. rm -rf wrfout_d01_* _pt_wrfout_d01_* _WRF.nc

import numpy as np
import xarray as xr
import pandas as pd
from netCDF4 import Dataset
import wrf # conda install -c conda-forge wrf-python (https://anaconda.org/conda-forge/wrf-python)
import sys
from pathlib import Path
import time
import f90nml # conda install -c conda-forge f90nml

if len(sys.argv) != 2:
	print("Usage:")
	print("       python prep_wrf_pres.py wrfout_FILENAME")
	exit()
else:
	src_filename = str(sys.argv[1])

def prep_var(ds, name):

	try:
		del ds.attrs['projection']
		del ds.attrs['coordinates']
	except:
		pass

	ds.name = name

	try:
		ds.level.attrs['description'] = 'pressure'
		ds.level.attrs['units'] = 'Pa'
		ds = ds.rename({'Time': 'time'})
	except:
		pass

	return ds

print(f"   Infile: {src_filename}")
file_name = Path(src_filename).stem
# file_name_split = file_name.split("_")
# outfile = f"_pt_{file_name_split[2]}_{file_name_split[3]}.nc"
outfile = f"{file_name}_pl.nc"
print(f"  Outfile: {outfile}")

with open('./lpt.nml') as nml_file:
	nml = f90nml.read(nml_file)
	levels = nml['wrf_prep']['levels']
	cartesian = nml['data']['cartesian_grid']

ds = xr.open_dataset(src_filename)
f = Dataset(src_filename)

print(" * Reading coordinates...")
lon2d = ds.XLONG.isel(Time=0)
lat2d = ds.XLAT.isel(Time=0)

print(" * Reading p...")
p = wrf.getvar(f, "p",wrf.ALL_TIMES)

print(" * Reading u...")
u = wrf.getvar(f, "ua",wrf.ALL_TIMES)
u = wrf.interplevel(u, p, levels)
u = prep_var(u, "u")
u.to_netcdf(outfile)
del u

print(" * Reading v...")
v = wrf.getvar(f, "va",wrf.ALL_TIMES)
v = wrf.interplevel(v, p, levels)
v = prep_var(v, "v")
v.to_netcdf(outfile, mode='a')
del v

print(" * Reading w...")
w = wrf.getvar(f, "wa",wrf.ALL_TIMES)
w = wrf.interplevel(w, p, levels)
w = prep_var(w, "w")
w.to_netcdf(outfile, mode='a')
del w

print(" * Reading geopotential...")
z = wrf.getvar(f, "geopotential",wrf.ALL_TIMES)
z = wrf.interplevel(z, p, levels)
z = prep_var(z, "z")
z.to_netcdf(outfile, mode='a')
del z

lon2d = prep_var(lon2d, "XLONG")
lat2d = prep_var(lat2d, "XLAT")

lon2d.to_netcdf(outfile, mode='a')
lat2d.to_netcdf(outfile, mode='a')
