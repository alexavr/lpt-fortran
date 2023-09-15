################################################
# Change this part
#NETCDF=/opt/netcdf-4
DATETIME=/opt/datetime-fortran-1.6.0
FC = ifort
################################################

################################################
# Do not edit the rest!!!
# (unless you know what you are doing)

# Program Name
PROG = lpt.exe

# Scr folder
VPATH = ./src/

# Compiler Flags
LINKER = $(FC) 
INCLUDE = -I$(NETCDF)/include -I$(DATETIME)/include 
LIBS = -qopenmp -L$(DATETIME)/lib -ldatetime -L$(NETCDF)/lib -lnetcdf -lnetcdff 
FFLAGS = -O3

# Object files
OBJS = module_globals.o module_io.o module_grid.o module_timetools.o lpt.o

LPT: $(PROG)

# Create the LPT
$(PROG): $(OBJS)
	@echo "--------------------------------------"
	@echo "Creating the executable for the LPT"
	@echo "--------------------------------------"
	$(LINKER) -o $(PROG) $(OBJS) $(LIBS)
	mv *.o *.mod $(VPATH)
# 	ln -sf $(VPATH)${PROG} .

%.o: %.f90
	@echo "--------------------------------------"
	@echo "Compiling the file $<"
	@echo "--------------------------------------"
	$(FC) -c $(FFLAGS) $(INCLUDE) $(LIBS) $<
	
# Clean up everything
clean:
	@echo "--------------------------------------"
	@echo "Cleaning everything up in LPT"
	@echo "--------------------------------------"
	rm -f $(VPATH)*.o $(VPATH)*.mod $(VPATH)*.exe $(PROG)

module_globals.o    : module_globals.f90
module_timetools.o  : module_timetools.f90 module_globals.o
module_io.o         : module_io.f90 module_globals.o
module_grid.o       : module_grid.f90 module_globals.o
lpt.o               : lpt.f90 module_io.o module_grid.o module_timetools.o
