#-----------------------------------------------------------------------------
# Makefile
#
# For building lndpm.
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ../../../..
include $(TOP)/make.config

#-----------------------------------------------------------------------------
# Set up compile options
CC    = gcc
RM    = rm
EXTRA = -Wall $(EXTRA_OPTIONS)

# Define the include files
INC = lndpm.h

# Define the source code and object files
SRC = lndpm.c
OBJ = $(SRC:.c=.o)

# Define the object libraries and paths
EXLIB = -L$(ESPALIB) -l_espa_raw_binary -l_espa_common \
        -L$(XML2LIB) -lxml2 \
        -L$(LZMALIB) -llzma \
        -L$(ZLIBLIB) -lz

MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

# Define include paths
INCDIR  = -I. -I$(ESPAINC) -I$(XML2INC)
NCFLAGS = $(EXTRA) $(INCDIR)

# Define C executables
EXE = lndpm

#-----------------------------------------------------------------------------
all: $(EXE)

$(EXE): $(OBJ) $(INC)
	$(CC) $(EXTRA) -o $(EXE) $(OBJ) $(LOADLIB)

#-----------------------------------------------------------------------------
install:
	install -d $(link_path)
	install -d $(ledaps_bin_install_path)
	install -m 755 $(EXE) $(ledaps_bin_install_path)
	ln -sf $(ledaps_link_source_path)/$(EXE) $(link_path)/$(EXE)


#-----------------------------------------------------------------------------
clean:
	$(RM) -f *.o $(EXE)

#-----------------------------------------------------------------------------
$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) -c $< -o $@
