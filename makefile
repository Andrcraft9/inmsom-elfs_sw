include Makefile.inc

SRCFIRSTSERVICE =	     \
	Service/read_write_parameters.f90

SRCCONTROL =	\
	Control/dyncall.f90      \
	Control/init_arrays.f90  \
	Control/init_pars.f90    \
	Control/output.f90

SRCSERVICE =	    \
    Service/grid_construction.f90     \
	Service/basin_parameters.f90      \
	Service/input_output_data.f90     \
	Service/rw_ctl_file.f90           \
	Service/time_tools.f90

SRCFUNCTION =	\
	Function/depth.f90   \
	Function/mixing.f90  \
	Function/vel_ssh.f90 \
	Function/shallow_water.f90

SRCMODULES = 	\
	Modules/mod_main_basin_pars.f90      \
	Modules/mod_rec_length.f90           \
	Modules/mod_hilbert_curve.f90            \
	Modules/mod_mpi_parallel_tools.f90   \
	Modules/mod_basin_grid.f90           \
	Modules/mod_ocean_variables.f90      \
	Modules/mod_time_integration.f90

all: inmsom clean

inmsom:
#order is important
	$(FC) -o inmsom $(FCFLAGS) $(SRCFIRSTSERVICE) $(SRCMODULES) $(SRCSERVICE) $(SRCFUNCTION) $(SRCCONTROL) inmsom_head.f90	
clean:
	$(RM) *.o *.mod
