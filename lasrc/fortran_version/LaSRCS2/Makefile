##############################################################################
# Compile commands from original implementation
##############################################################################
# cc -O0 -c setprojinfo.c
# cc -O0 -c qa_index_routine.c
#
# gfortran -ffixed-line-length-132  -O -o $1 $1.f sub_deg2utm.f sub_utm2deg.f subaeroretv4.f LUT-MSI_subr_v4.f subaeroretwat.f qa_index_routine.o setprojinfo.o -I /usr/local/opt/TOOLKIT-5.2.17/hdf/linux64/hdf-4.2.5/include -L /usr/local/opt/TOOLKIT-5.2.17/hdf/linux64/hdf-4.2.5/lib -lmfhdf -ldf -ljpeg -lz -ldl -lm -lsz
##############################################################################

CC = cc
FC = gfortran

INCDIR  = -I. -I$(JPEGINC) -I$(TIFFINC) -I$(GEOTIFF_INC) -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC)
LIBDIR  = -L. -L$(JPEGLIB) -L$(TIFFLIB) -L$(GEOTIFF_LIB) -L$(HDFLIB) -L$(HDFEOS_LIB) -L$(HDFEOS_GCTPLIB)

C_EXTRA   = -Wall -O2
NCFLAGS   = $(CFLAGS) $(C_EXTRA) $(INCDIR)

F_EXTRA   = -Wall -O -ffixed-line-length-132 -cpp

C_OBJECTS = qa_index_routine.o setprojinfo.o

F_OBJECTS = sub_deg2utm.o sub_utm2deg.o subaeroretwat.o subaeroretv4.o LUT-MSI_subr_v4.o LaSRCS2AV3.5.5.o

TARGET = LaSRCS2AV3.5.5

all: $(TARGET)

$(TARGET): $(C_OBJECTS) $(F_OBJECTS)
	$(FC) $(F_EXTRA) -o $(TARGET) $(F_OBJECTS) $(C_OBJECTS) $(INCDIR) $(LIBDIR) -lmfhdf -ldf -ljpeg -lz -ldl -lm

clean:
	rm -f $(TARGET) *.o

#
# Rules
#
.c.o:
	$(CC) $(NCFLAGS) -I$(HDFINC) -c $< -o $@

.f.o:
	$(FC) $(F_EXTRA) -c $< -o $@

