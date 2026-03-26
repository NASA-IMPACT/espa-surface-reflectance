##############################################################################
# Compile commands from original implementation
##############################################################################
# gcc -O0 -o tiff2hdf *.c -I/data/cots/hdf/include -I/data/cots/libgeotiff-1.2.5-64bit/include -L/data/cots/libgeotiff-1.2.5-64bit/lib -ltiff -L/data/cots/hdf/lib -lmfhdf -ldf -ljpeg -L/data/cots/szip-2.1/lib -lz -lsz -lm -lproj
##############################################################################

CC = cc

INCDIR  = -I. -I$(JPEGINC) -I$(TIFFINC) -I$(GEOTIFF_INC) -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC)
LIBDIR  = -L. -L$(JPEGLIB) -L$(SZIPLIB) -L$(TIFFLIB) -L$(GEOTIFF_LIB) -L$(HDFLIB) -L$(HDFEOS_LIB) -L$(HDFEOS_GCTPLIB)

C_EXTRA   = -Wall -O2
NCFLAGS   = $(CFLAGS) $(C_EXTRA) $(INCDIR)

C_OBJECTS = tif2hdf.o TIFF_reader2.o

#TARGET = TIFF_reader2
TARGET = tif2hdf

all: $(TARGET)

$(TARGET): $(C_OBJECTS)
	$(CC) $(C_EXTRA) -o $(TARGET) $(C_OBJECTS) $(INCDIR) $(LIBDIR) -ltiff -lgeotiff -lmfhdf -ldf -ljpeg -lz -lsz -lm

clean:
	rm -f $(TARGET) *.o

#
# Rules
#
.c.o:
	$(CC) $(NCFLAGS) -c $< -o $@

