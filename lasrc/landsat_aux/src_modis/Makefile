#-----------------------------------------------------------------------------
# Makefile for combine L8 auxiliary code
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
INC = combine_l8_aux_data.h

# Define the source code and object files
SRC = get_args.c            \
      combine_l8_aux_data.c
OBJ = $(SRC:.c=.o)

# Define include paths
INCDIR = -I. -I$(ESPAINC) -I$(XML2INC)
HDF_INCDIR = -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC)
NCFLAGS  = $(EXTRA) $(INCDIR) $(HDF_INCDIR)

# Define the object libraries
EXLIB = -L$(ESPALIB) -l_espa_common \
        -L$(XML2LIB) -lxml2 \
        -L$(LZMALIB) -llzma \
        -L$(SZIPLIB) -lsz \
        -L$(ZLIBLIB) -lz
HDF_EXLIB = -L$(HDFEOS_LIB) -lhdfeos \
            -L$(HDFLIB) -lmfhdf -ldf \
            -L$(JPEGLIB) -ljpeg \
            -L$(HDFEOS_GCTPLIB) -lGctp
MATHLIB = -lm
LOADLIB = $(HDF_EXLIB) $(EXLIB) $(MATHLIB)

# Define C executable s
EXE = combine_l8_aux_data

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
