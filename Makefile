################################################
# Change this part
NETCDF=/opt/netcdf4-serial
DATETIME=/opt/datetime-fortran-1.6.0
FC = ifort
################################################

################################################
# Do not edit part
# (unless you know what you are doing)

# Program Name
PROG = lpt.exe

# Scr folder
VPATH = ./src/

# Compiler Flags
FFLAGS = -I$(NETCDF)/include -I$(DATETIME)/include
FLINK = -O3 -L$(DATETIME)/lib -ldatetime -L$(NETCDF)/lib -lnetcdf -lnetcdff
LINKER = $(FC) -o

# Object files
OBJS = module_io.o module_globals.o module_grid.o module_timetools.o pt.o

model: $(PROG)

# Creates the model
$(PROG): $(OBJS)
	@echo "--------------------------------------"
	@echo "Creating the executable for the PT"
	@echo "--------------------------------------"
	$(LINKER) $(PROG) $(OBJS) $(FLINK)
	mv *.o *.mod $(VPATH)
	ln -sf $(VPATH)${PROG} .

%.o: %.f90
	@echo "--------------------------------------"
	@echo "Compiling the file $<"
	@echo "--------------------------------------"
	$(FC) -c $(FFLAGS) $<
	
# Cleans up everything
clean:
	@echo "--------------------------------------"
	@echo "Cleaning everything up in model"
	@echo "--------------------------------------"
	rm -f $(VPATH)*.o $(VPATH)*.mod $(VPATH)*.exe pt.exe

module_globals.o: module_globals.f90
module_timetools.o  : module_timetools.f90 module_globals.o
module_io.o     : module_io.f90 module_globals.o
module_grid.o   : module_grid.f90 module_globals.o
pt.o            : pt.f90 module_io.o module_grid.o module_timetools.o