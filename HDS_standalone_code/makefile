# Compiler settings
FC = gfortran
COMMON_FCFLAGS = -Wall -Wextra -Wconversion -fbacktrace -ffree-line-length-0
DEBUG_FCFLAGS = -g -O0
RELEASE_FCFLAGS = -O3

# Source files
hds_src= \
		type_HDS.f90 \
        HDS.f90

# model driver
hds_drv = main_HDS.f90

# Executable name
hds_exe = HDS_program
path_exe = ../bin/

# Main target
all: release compile link clean install

debug: COMMON_FCFLAGS += $(DEBUG_FCFLAGS)
debug: $(EXECUTABLE)

release: COMMON_FCFLAGS += $(RELEASE_FCFLAGS)
release: $(EXECUTABLE)

# compile
compile:
	$(FC) $(COMMON_FCFLAGS) -c $(hds_src) $(hds_drv) \
	$(INCLUDES)

# link routines
link:
	$(FC) -g *.o -o $(hds_exe)

clean:
	rm -f *.o
	rm -f *.mod

install:
	@mkdir -p $(path_exe)
	@mv $(hds_exe) $(path_exe)
	$(info $(hds_exe) successfully installed in $(path_exe))
