*Detailed description is coming soon...*

[TOC]

The Lagrangian Particle Tracing of weightless particles. The scheme is based on pure interpolation (no modeling involved!). Basic idea: take the current particle coordinates and predict the future location base on its current velocity. As simple as that. The scheme uses 3D velocity components and geopotential fields.

## Requirements

[NetCDF](https://www.unidata.ucar.edu/), [datetime-fortran](https://github.com/wavebitscientific/datetime-fortran).

## Building 

First `vi Makefile` :

```bash
NETCDF=/opt/netcdf4-serial     # to add path to your NetCDF
DATETIME=/opt/datetime-fortran # add path to your datetime-fortran libs
FC=ifort                       # your FORTRAN compiler
```

Normally you don't need to change other stuff.

Building:

```bash
make
```

Cleaning in case anything went wrong:

```bash
make clean
```

## Data reprocessing

The  source file has to be properly prepared. 

First, all data has to be merged into a single file. The file should have 3D variables for wind velocity components and geopotential named as *u*,*v*,*w*,*z* and coordinates (time and location). 

The time dimension has to be relative (`second/minutes/hours since YYYY-MM-DD HH:mm:SS`) with attributes: units and calendar (normally its "gregorian").

If your data is global (reanalyzes?) then vertical velocity (*w*) is expected to have Pa/s unit (as reanalyzes are normally provide) and coordinates are: *longitude* or *lon* and *latitude* or *lat* for horizontal location. For vertical coordinate: *levels* (in Pa). All has to be vector (1D).

If you date is regional (assuming the WRF output for now, but could be anything with little modifications) then coordinates should be names as: *XLONG* and *XLAT* and expected to be 2D (to moving domains allowed for now). 

All this preparation could be accomplished using the [cdo](https://code.mpimet.mpg.de/projects/cdo/) library.

## Configuration

The configuration is done using the file *namelist.pt*.

The domain could be regional (`regional = .TRUE.`) or global (periodic boundary, `regional = .FALSE.`). You need also set the time interval using `stime` and `etime` (this interval has to be in your source file!).

The particles could be initiated manually (`pt_grid = .FALSE.`) in ASCII file (see *points_test_Reykjavik.dat* as an example). Or evenly spaced allover the domain (`pt_grid = .TRUE.`), `pt_step ` is number of steps between particles in both directions, `pt_height` is elevation in meters. 

Accuracy of the computations. There are two schemes provided: 

1. Coarse scheme (`accuracy = .FALSE.`). Predicts particle location using time step from the source file. Fast but very inaccurate. Recommended only for test purposes. 
2. Fine scheme (`accuracy = .TRUE.`). Interpolates data between source file time step. This option controlled by `timestep` (in minutes).

## Running

```bash
./lpt.exe src_file.nc
```

where *src_file.nc* comes from the "Data reprocessing".



