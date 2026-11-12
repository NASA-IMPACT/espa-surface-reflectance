##############################################################################
# Compile commands from original implementation
##############################################################################
# gfortran -ffixed-line-length-132 -O -o jdoy jdoy.f 
##############################################################################

FC = gfortran

F_EXTRA   = -Wall -O2 -ffixed-line-length-132

F_OBJECTS = jdoy.o

TARGET = jdoy

all: $(TARGET)

$(TARGET): $(F_OBJECTS)
	$(FC) $(F_EXTRA) -o $(TARGET) $(F_OBJECTS)

clean:
	rm -f $(TARGET) *.o

#
# Rules
#
.f.o:
	$(FC) $(F_EXTRA) -c $< -o $@

