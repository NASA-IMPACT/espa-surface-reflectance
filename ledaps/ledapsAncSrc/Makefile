##############################################################################
# Makefile for ncep and ozone dataset repackaging.
# 
# Note: In order for the netCDF4 NCEP libraries and includes to work correctly
# with the HDF4 libraries and includes, the HDF4 libraries need to be built
# with netCDF disabled.  --disable-netcdf
##############################################################################
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ../..
include $(TOP)/make.config

#-----------------------------------------------------------------------------
# Set up compile options
CC    = gcc
RM    = rm
EXTRA = -Wall $(EXTRA_OPTIONS)

# Define the source code and object files
SRC1 = convert_ozone.c
OBJ1 = $(SRC1:.c=.o)

SRC2 = ncep_reanalysis_surface_repackage.o
OBJ2 = $(SRC2:.c=.o)

ALL_OBJ = $(OBJ1) $(OBJ2)

# Define include paths
INCDIR  = -I. -I$(NCDF4INC) -I$(HDF5INC) -I$(HDFINC) -I$(JPEGINC) -I$(ZLIBINC)
NCFLAGS = $(EXTRA) $(INCDIR)

# Define the object libraries and paths
EXLIB = -L$(CURLLIB) -lcurl \
        -L$(IDNLIB) -lidn -ldl -lrt \
        -L$(SZIPLIB) -lsz \
        -L$(ZLIBLIB) -lz
HDF_EXLIB = -L$(HDFLIB) -lmfhdf -ldf \
            -L$(JPEGLIB) -ljpeg
NCDF_EXLIB = -L$(NCDF4LIB) -lnetcdf \
             -L$(HDF5LIB) -lhdf5_hl -lhdf5 -pthread
MATHLIB = -lm
LOADLIB1 = $(HDF_EXLIB) $(EXLIB) $(MATHLIB)
LOADLIB2 = $(HDF_EXLIB) $(NCDF_EXLIB) $(EXLIB) $(MATHLIB)

# Define scripts
SCRIPTS = updatencep.py updatetoms.py

# Define C executables
EXE1 = convert_ozone
EXE2 = ncep_repackage
ALL_EXE = $(EXE1) $(EXE2)

#-----------------------------------------------------------------------------
all: $(ALL_OBJ) $(ALL_EXE)

$(EXE1): $(OBJ1)
	$(CC) $(EXTRA) $(OBJ1) -o $(EXE1) $(LOADLIB1)


$(EXE2): $(OBJ2)
	$(CC) $(EXTRA) $(OBJ2) -o $(EXE2) $(LOADLIB2)

#-----------------------------------------------------------------------------
install:
	install -d $(link_path)
	install -d $(ledaps_bin_install_path)
	@for executable in $(ALL_EXE); do \
            cmd="install -m 755 $$executable $(ledaps_bin_install_path)"; \
            echo "$$cmd"; $$cmd || exit 1; \
            cmd="ln -sf $(ledaps_link_source_path)/$$executable $(link_path)/$$executable"; \
            echo "$$cmd"; $$cmd; \
        done
	@for script in $(SCRIPTS); do \
            cmd="install -m 755 $$script $(ledaps_bin_install_path)"; \
            echo "$$cmd"; $$cmd || exit 1; \
            cmd="ln -sf $(ledaps_link_source_path)/$$script $(link_path)/$$script"; \
            echo "$$cmd"; $$cmd; \
        done

#-----------------------------------------------------------------------------
clean:
	rm -f *.o $(EXE1) $(EXE2)

#-----------------------------------------------------------------------------
$(OBJ1): $(INC1)
$(OBJ2): $(INC2)

.c.o:
	$(CC) $(NCFLAGS) -c $< -o $@
