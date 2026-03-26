#-----------------------------------------------------------------------------
# Makefile for gapfilling the VIIRS auxiliary data
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ../../..
include $(TOP)/make.config

#-----------------------------------------------------------------------------
# Set up compile options
CC    = gcc
RM    = rm
EXTRA = -Wall $(EXTRA_OPTIONS)

# Define the include files
INC = gapfill_viirs_aux.h

# Define the source code and object files
SRC = get_args.c            \
      gapfill_viirs_aux.c
OBJ = $(SRC:.c=.o)

# Define include paths
INCDIR = -I. -I$(ESPAINC) -I$(XML2INC)
HDF5_INCDIR = -I$(HDF5INC) -I$(HDFEOS5_INC)
NCFLAGS  = $(EXTRA) $(INCDIR) $(HDF5_INCDIR)

# Define the object libraries
EXLIB = -L$(ESPALIB) -l_espa_common \
        -L$(XML2LIB) -lxml2 \
        -L$(LZMALIB) -llzma \
        -L$(SZIPLIB) -lsz \
        -L$(ZLIBLIB) -lz
HDF5_EXLIB = -L$(HDF5LIB) -lhdf5_hl -lhdf5 \
             -L$(JPEGLIB) -ljpeg
MATHLIB = -lm
LOADLIB = $(HDF5_EXLIB) $(EXLIB) $(MATHLIB)

# Define C executable s
EXE = gapfill_viirs_aux

#-----------------------------------------------------------------------------
all: $(EXE)

$(EXE): $(OBJ) $(INC)
	$(CC) $(EXTRA) -o $(EXE) $(OBJ) $(LOADLIB)

#-----------------------------------------------------------------------------
install:
	install -d $(link_path)
	install -d $(lasrc_bin_install_path)
	install -m 755 $(EXE) $(lasrc_bin_install_path)
	ln -sf $(lasrc_link_source_path)/$(EXE) $(link_path)/$(EXE) || exit 1

#-----------------------------------------------------------------------------
clean:
	$(RM) *.o $(EXE)

#-----------------------------------------------------------------------------
$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) -c $<
