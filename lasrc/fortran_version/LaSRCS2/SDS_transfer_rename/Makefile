##############################################################################
# Compile commands from original implementation
##############################################################################
# gcc -O0 -o SDS_transfer_rename SDS_transfer_rename.c -I/usr/local/TOOLKIT/SDPToolkit_gfortran44-5.2.17/hdf/linux64/hdf-4.2.5/include/ -L/usr/local/TOOLKIT/SDPToolkit_gfortran44-5.2.17/hdf/linux64/hdf-4.2.5/lib/ -lmfhdf -ldf -ljpeg -lz -lsz -lm
##############################################################################

CC = cc

INCDIR  = -I. -I$(JPEGINC) -I$(TIFFINC) -I$(GEOTIFF_INC) -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC)
LIBDIR  = -L. -L$(JPEGLIB) -L$(TIFFLIB) -L$(GEOTIFF_LIB) -L$(HDFLIB) -L$(HDFEOS_LIB) -L$(HDFEOS_GCTPLIB)

C_EXTRA   = -Wall -O2
NCFLAGS   = $(CFLAGS) $(C_EXTRA) $(INCDIR)

C_OBJECTS = SDS_transfer_rename.o

TARGET = SDS_transfer_rename

all: $(TARGET)

$(TARGET): $(C_OBJECTS)
	$(CC) $(C_EXTRA) -o $(TARGET) $(C_OBJECTS) $(INCDIR) $(LIBDIR) -lmfhdf -ldf -ljpeg -lz -ldl -lm

clean:
	rm -f $(TARGET) *.o

#
# Rules
#
.c.o:
	$(CC) $(NCFLAGS) $(INCDIR) -c $< -o $@

