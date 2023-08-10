# Compiler settings
FC = gfortran
COMMON_FCFLAGS = -Wall -Wextra -Wconversion -fbacktrace -ffree-line-length-0
DEBUG_FCFLAGS = -g -O0
RELEASE_FCFLAGS = -O3

# Source files
SRC_MAIN = main_HDS.f90
SRC_MODULE = HDS.f90

# Object files
OBJ_MAIN = $(SRC_MAIN:.f90=.o)
OBJ_MODULE = $(SRC_MODULE:.f90=.o)

# Module file
MOD_FILE = hds.mod

# Executable name
EXECUTABLE = HDS_program

# Main target
all: release

debug: COMMON_FCFLAGS += $(DEBUG_FCFLAGS)
debug: $(EXECUTABLE)

release: COMMON_FCFLAGS += $(RELEASE_FCFLAGS)
release: $(EXECUTABLE)

$(EXECUTABLE): $(OBJ_MAIN) $(OBJ_MODULE)
	$(FC) $(COMMON_FCFLAGS) -o $@ $(OBJ_MAIN) $(OBJ_MODULE)
	rm -f $(OBJ_MAIN) $(OBJ_MODULE) $(MOD_FILE)

%.o: %.f90 $(MOD_FILE)
	$(FC) $(COMMON_FCFLAGS) -c $<

$(MOD_FILE): $(SRC_MODULE)
	$(FC) $(COMMON_FCFLAGS) -c $(SRC_MODULE)

clean:
	rm -f $(EXECUTABLE) $(OBJ_MAIN) $(OBJ_MODULE) $(MOD_FILE)
