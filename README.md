*Detailed description is coming soon...*

The lagrangian tracing of weightless particles. The scheme is based on pure interpolation (no modeling involved!). Basic idea: take the particle coordinates and predict the future location based on its current velocity. As simple as that. The scheme uses 3D velocity components and geopotential fields.

## Requirements

[NetCDF](https://www.unidata.ucar.edu/), [datetime-fortran](https://github.com/wavebitscientific/datetime-fortran).

## Building 

First `vi Makefile` :

```bash
NETCDF=/opt/netcdf4-hdf5       # path to your NetCDF
DATETIME=/opt/datetime-fortran # path to your datetime-fortran libs
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

## Data preprocessing

The source file should contain 3D velocity and geopotential fields. But it has to be properly prepared. 

First, all data has to be merged into a single file. The file should have 3D variables for wind velocity components and geopotential named as *u*,*v*,*w*,*z* and coordinates (time and location). 

The time dimension has to be relative (`second/minutes/hours since YYYY-MM-DD HH:mm:SS`) with attributes: units and calendar (normally its "gregorian").

If your data is global (reanalyzes?) then vertical velocity (*w*) is expected to have Pa/s unit (as reanalyzes are normally provide) and coordinates are: *longitude* or *lon* and *latitude* or *lat* for horizontal location. For vertical coordinate: *levels* (in Pa). All has to be vector (1D).

If you data is regional (assuming the WRF output for now, but could be anything with little modifications) then coordinates should be names as: *XLONG* and *XLAT* and expected to be 2D (moving domains are not allowed for now). 

All this preparation could be accomplished by using the [cdo](https://code.mpimet.mpg.de/projects/cdo/) library. For preprocessing of WRF output you might want to use *prep_wrf.py* script (works, but written in a very hurry). 

## Configuration

The configuration is done using the file *namelist.pt*.

The domain could be regional (`regional = .TRUE.`) or global (periodic boundary, `regional = .FALSE.`). You need also set the time interval using `stime` and `etime` (this interval has to be present in the source file!).

The particles could be initiated manually (`pt_grid = .FALSE.`) in ASCII file named *init_particles.dat* (see *init_particles_test_Reykjavik.dat* or *init_particles_test_NorthPole.dat* as an example). 

Or evenly spaced allover the domain (`pt_grid = .TRUE.`), `pt_step` is number of steps (grid pionts) between particles in both directions, `pt_height` is an elevation in meters (could be an array, but haven't tested it yet). 

The accuracy could be controlled. There are two schemes provided: 
1. Coarse scheme (`accuracy = .FALSE.`). Predicts particle location using time step from the source file. Fast but very inaccurate. Recommended only for test purposes. 
2. Fine scheme (`accuracy = .TRUE.`). Interpolates data between source file time step. This option controlled by `timestep` (in minutes). The less the better but becomes very computationally heavy. 

## Running

```bash
./lpt.exe src_file.nc
```
where *src_file.nc* comes from the "Data reprocessing". The output will be *lptf90_src_file.nc*

## Plans

- [ ] Implement backward tracing
- [ ] Pix Pole issue -- sometimes particles are tend to stuck at Pole for a longer then usual (run with init_particles_test_NorthPole.dat to see)
