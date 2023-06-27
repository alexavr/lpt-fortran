#!/usr/bin/env python3
# conda install -c conda-forge wrf-python
# 1. ln -sf /storage/NAD/NAAD/v4/LoRes/2015/wrfout_d01_2015-01-* .
# 2. for i in $(ls wrfout_d01_*); do ./prep_wrf.py $i; done
# 3. cdo mergetime _pt_wrfout_d01_* _WRF.nc
# 4. cdo settunits,hours _WRF.nc wrf_2015-01.nc
# 5. rm -rf wrfout_d01_* _pt_wrfout_d01_* _WRF.nc
import sys
from netCDF4 import Dataset
import numpy as np
import wrf as wrf

if len(sys.argv) != 2:
	print("Usage:")
	print("      ./prep_wrf.py wrfout_file")
	exit()

ifile = str(sys.argv[1])
ofile = "_pt_"+ifile

print(f"Input file:  {ifile}")
print(f"Output file: {ofile}")

incid = Dataset(ifile,"r")
times = wrf.extract_times(incid, None, method='cat', squeeze=True, cache=None, meta=True, do_xtime=True)

ntimes = times.shape[0]

var = incid.variables["U"][:]
u = wrf.destagger(var, 3, meta=False)
var = incid.variables["V"][:]
v = wrf.destagger(var, 2, meta=False)
var = incid.variables["W"][:]
w = wrf.destagger(var, 1, meta=False)
ph  = incid.variables["PH"][:]
phb = incid.variables["PHB"][:]
p  = incid.variables["P"][:]
pb = incid.variables["PB"][:]
# slp = wrf.getvar(incid, "slp", timeidx=wrf.ALL_TIMES)
hgt = incid.variables["HGT"][0,:,:]

# Since geopotential hight in WRF goes from mean Earth hight
# We need to substruct hgt to make it be from the real surface
# 1. No correction: 
# var = (ph+phb)
# z = wrf.destagger(var, 1, meta=False)
# 2. With correction: 
var = ph+phb
z = wrf.destagger(var, 1, meta=False)
p = p+pb
h = z/9.81
for it in range(0,z.shape[0]):
	for iz in range(0,z.shape[1]):
		z[it,iz,:,:] = z[it,iz,:,:] - hgt*9.81
		h[it,iz,:,:] = h[it,iz,:,:] - hgt

del var, ph, phb, hgt

xlat  = incid.variables["XLAT"][0,:,:]
xlong = incid.variables["XLONG"][0,:,:]

oncid = Dataset(ofile,'w')
idtime = oncid.createDimension('time', u.shape[0])
idz = oncid.createDimension('bottom_top', u.shape[1])
idy = oncid.createDimension('south_north', u.shape[2])
idx = oncid.createDimension('west_east', u.shape[3])

vtime = oncid.createVariable('time', np.float64, ('time',))
vtime.units = times.units
vtime.long_name = 'time'
vtime.calendar = 'gregorian'

vlat = oncid.createVariable('XLAT' , np.float32, ('south_north','west_east'))
vlon = oncid.createVariable('XLONG', np.float32, ('south_north','west_east'))
# vslp = oncid.createVariable('slp', np.float32, ('time','south_north','west_east'))
vu = oncid.createVariable('u', np.float32, ('time','bottom_top','south_north','west_east'))
vv = oncid.createVariable('v', np.float32, ('time','bottom_top','south_north','west_east'))
vw = oncid.createVariable('w', np.float32, ('time','bottom_top','south_north','west_east'))
vz = oncid.createVariable('z', np.float32, ('time','bottom_top','south_north','west_east'))
vh = oncid.createVariable('h', np.float32, ('time','bottom_top','south_north','west_east'))
vp = oncid.createVariable('p', np.float32, ('time','bottom_top','south_north','west_east'))

vtime[:] = times[:]
vlat[:] = xlat
vlon[:] = xlong
vu[:] = u
vv[:] = v
vw[:] = w
vz[:] = z
# vslp[:] = slp
vh[:] = h
vp[:] = p

del u, v, w, z, h, p, xlat, xlong, times
# del slp
