#!/usr/bin/env python3

# Интерполируем сырые данные WRF на уровни

# Подготовка окружения:
# conda activate data
# python prep_wrf_pres.py {filename}.nc


# Настройка конфигурации подготовки (в разделе 'wrf_prep' lpt.nml):
# 'type' - ставим тип уровней, на которые хотим интерполировать (pressure, height или ml)
# 'levels' - список уровней в размерности type (метры, Па или номера уровней)

# Запуск подготовки:
# 1. ln -sf /storage/NAD/NAAD/v4/LoRes/2015/wrfout_d01_2015-01-* . 	# линкуем сырые данные в текущую папку
# 2. for i in $(ls wrfout_d01_*); do ./prep_wrf.py $i; done			# каждый файл интерполируем на уровни высот (используем пакет python-wrf)
# 3. cdo mergetime _{interp_level}_wrfout_d01_* {RESULT}_hgt.nc 					# склеиваем полученные в п.2 файлы в единый nc файл
# 4. rm -rf wrfout_d01_* _wrfout_d01_* _WRF.nc 						# удаляем ненужное


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
	print("       python prep_wrf.py wrfout_FILENAME")
	exit()
else:
	src_filename = str(sys.argv[1])

def prep_var(ds, name, nml):

	try:
		del ds.attrs['projection']
		del ds.attrs['coordinates']
	except:
		pass

	ds.name = name

	try:
		ds.level.attrs['description'] = nml['wrf_prep']['type']

		if ( ds.level.attrs['description'] == "height" ):
			ds.level.attrs['units'] = 'm'
		elif (ds.level.attrs['description'] == "pressure"):
			ds.level.attrs['units'] = 'Pa'
		elif (ds.level.attrs['description'] == "ml"):
			ds.level.attrs['units'] = 'ml'
		else:
			print(f" type of vertical levels has to be 'pressure', 'height' or 'ml' only")
			exit()

		ds = ds.rename({'Time': 'time'})
	except:
		pass

	return ds

start_time = time.time()
print(f"START ############################################################")

with open('./lpt.nml') as nml_file:
	nml = f90nml.read(nml_file)
	levels = nml['wrf_prep']['levels']
	cartesian = nml['data']['cartesian_grid']
	ideal_case = nml['data']['ideal_case']

print(f"Infile: {src_filename}")
file_name = Path(src_filename).stem
# file_name_split = file_name.split("_")
# outfile = f"_pt_{file_name_split[2]}_{file_name_split[3]}.nc"
outfile = f"_{nml['wrf_prep']['type']}_{file_name}.nc"
print(f"Outfile: {outfile}")

ds = xr.open_dataset(src_filename)
f = Dataset(src_filename)

print(f" * Reading {nml['wrf_prep']['type']}...")
if nml['wrf_prep']['type'] == "ml":
	levs = nml['wrf_prep']['levels']
else:
	levs = wrf.getvar(f, nml['wrf_prep']['type'],wrf.ALL_TIMES)

print(" * Reading u...")
u = wrf.getvar(f, "ua",wrf.ALL_TIMES)
if not (nml['wrf_prep']['type'] == "ml"):
	u = wrf.interplevel(u, levs, levels)
u = prep_var(u, "u",nml)
u.to_netcdf(outfile)
del u

print(" * Reading v...")
v = wrf.getvar(f, "va",wrf.ALL_TIMES)
if not (nml['wrf_prep']['type'] == "ml"):
	v = wrf.interplevel(v, levs, levels)
v = prep_var(v, "v",nml)
v.to_netcdf(outfile, mode='a')
del v

print(" * Reading w...")
w = wrf.getvar(f, "wa",wrf.ALL_TIMES)
if not (nml['wrf_prep']['type'] == "ml"):
	w = wrf.interplevel(w, levs, levels)
w = prep_var(w, "w",nml)
w.to_netcdf(outfile, mode='a')
del w

print(" * Reading geopotential...")
z = wrf.getvar(f, "geopotential",wrf.ALL_TIMES)
if not (nml['wrf_prep']['type'] == "ml"):
	z = wrf.interplevel(z, levs, levels)
z = prep_var(z, "z",nml)
z.to_netcdf(outfile, mode='a')
del z



print(" * Coordinates...")
if ideal_case: print(" * * It's ideal case -> replacing coordinates with indexes of grid")

if ideal_case:
    nx = ds.XLAT.shape[2]
    ny = ds.XLAT.shape[1]
    x = np.arange(1, nx+1)
    y = np.arange(1, ny+1)
    X, Y = np.meshgrid(x, y)
    lon2d = xr.Dataset({'XLONG': (('south_north', 'west_east'), X)})
    lat2d = xr.Dataset({'XLAT': (('south_north', 'west_east'), Y)})
else:
	lon2d = ds.XLONG.isel(Time=0)
	lat2d = ds.XLAT.isel(Time=0)
	lon2d = prep_var(lon2d, "XLONG")
	lat2d = prep_var(lat2d, "XLAT")

lon2d.to_netcdf(outfile, mode='a')
lat2d.to_netcdf(outfile, mode='a')

f.close()

with Dataset(src_filename) as src, Dataset(outfile, "a") as dst:
	dst.setncatts(src.__dict__)

print(f"END   ############################################################ {(time.time() - start_time):.1}s")

