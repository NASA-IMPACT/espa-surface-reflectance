#-----------------------------------------------------------------------------
# Makefile
#
# For building lndsr.
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ../../../..
include $(TOP)/make.config

#-----------------------------------------------------------------------------
# Set up compile options
CC    = gcc
FC    = gfortran
RM    = rm
EXTRA = -Wall $(EXTRA_OPTIONS)
LNDPM = ../lndpm

# Define the include files
C_INC = ar.h clouds.h const.h date.h error.h grib.h input.h keyvalue.h \
        lndsr.h lut.h myhdf.h mystring.h output.h param.h prwv_input.h \
        read_grib_tools.h sixs_runs.h sr.h

# Define the source code and object files
C_SRC = \
        ar.c              \
        clouds.c          \
        date.c            \
        error.c           \
        grib.c            \
        input.c           \
        lndsr.c           \
        lut.c             \
        myhdf.c           \
        mystring.c        \
        output.c          \
        param.c           \
        prwv_input.c      \
        read_grib_tools.c \
        sixs_runs.c       \
        sr.c
C_OBJ = $(C_SRC:.c=.o)

F_SRC = \
        CHAND.f \
        CSALBR.f
F_OBJ = $(F_SRC:.f=.o)

ALL_OBJ = $(C_OBJ) $(F_OBJ)

# Define include paths
INCDIR  = -I. -I${LNDPM} -I$(ESPAINC) -I$(XML2INC)
HDF_INCDIR = -I$(JPEGINC) -I$(HDFINC) -I$(HDFEOS_GCTPINC)
NCFLAGS = $(EXTRA) $(INCDIR) $(HDF_INCDIR)

# Define the object libraries and paths
EXLIB = -L$(ESPALIB) -l_espa_raw_binary -l_espa_common -l_espa_format_conversion \
        -L$(XML2LIB) -lxml2 \
        -L$(LZMALIB) -llzma \
        -L$(ZLIBLIB) -lz
HDF_EXLIB = -L$(HDFLIB) -lmfhdf -ldf \
            -L$(SZIPLIB) -lsz \
            -L$(JPEGLIB) -ljpeg \
            -L$(HDFEOS_GCTPLIB) -lGctp
MATHLIB = -lm
LOADLIB = $(EXLIB) $(HDF_EXLIB) $(MATHLIB)

# Define C executables
EXE = lndsr

#-----------------------------------------------------------------------------
all: $(EXE)

$(EXE): $(ALL_OBJ)
	$(CC) $(EXTRA) -o $(EXE) $(ALL_OBJ) $(LOADLIB)

#-----------------------------------------------------------------------------
install:
	install -d $(link_path)
	install -d $(ledaps_bin_install_path)
	install -m 755 $(EXE) $(ledaps_bin_install_path)
	ln -sf $(ledaps_link_source_path)/$(EXE) $(link_path)/$(EXE)

#-----------------------------------------------------------------------------
clean:
	rm -f *.o $(EXE)

#-----------------------------------------------------------------------------
$(C_OBJ): $(C_SRC) $(C_INC)

.c.o:
	$(CC) $(NCFLAGS) -c $< -o $@

$(F_OBJ): $(F_SRC)

.f.o:
	$(FC) $(NCFLAGS) -c $< -o $@
