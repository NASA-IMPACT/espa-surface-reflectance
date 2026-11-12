#-----------------------------------------------------------------------------
# Makefile
#
# For building 6S
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

# Define the source code and object files
F_SRC = AATSR.f ABSTRA.f AEROSO.f AKTOOL.f ATMREF.f AVHRR.f BBM.f BDM.f \
        BRDFGRID.f CHAND.f CLEARW.f CSALBR.f DICA1.f DICA2.f DICA3.f \
        DISCOM.f DISCRE.f DUST.f ENVIRO.f EQUIVWL.f GAUSS.f GLI.f GOES.f \
        HAPKALBE.f HAPKBRDF.f HRV.f IAPIALBE.f IAPIBRDF.f IAPITOOLS.f \
        INTERP.f ISO.f KERNEL.f KERNELPOL.f LAKEW.f MAS.f MERIS.f METEO.f \
        METH1.f METH2.f METH3.f METH4.f METH5.f METH6.f MIDSUM.f MIDWIN.f \
        MIE.f MINNALBE.f MINNBRDF.f MOCA1.f MOCA2.f MOCA3.f MOCA4.f MOCA5.f \
        MOCA6.f MODIS.f MSS.f NIOX1.f MODISBRDF.f MODISALBE.f NIOX2.f \
        NIOX3.f NIOX4.f NIOX5.f NIOX6.f OCEA.f OCEAALBE.f OCEABRDF.f \
        OCEABRDFFAST.f OCEATOOLS.f ODA550.f ODRAYL.f OS.f OSPOL.f OXYG3.f \
        OXYG4.f OXYG5.f OXYG6.f OZON1.f PLANPOL.f POLDER.f POSGE.f POSGW.f \
        POSLAN.f POSMTO.f POSNOA.f POSSOL.f POSSPO.f POLGLIT.f POLNAD.f \
        PRESPLANE.f PRESSURE.f PRINT_ERROR.f RAHMALBE.f RAHMBRDF.f \
        ROUJALBE.f ROUJBRDF.f SAND.f SCATRA.f SEAWIFS.f SOLIRR.f SOOT.f \
        SPECINTERP.f SPLIE2.f SPLIN2.f SPLINE.f SPLINT.f STM.f SUBSUM.f \
        SUBWIN.f TM.f TROPIC.f TRUNCA.f US62.f VARSOL.f VEGETA.f VERSALBE.f \
        VERSBRDF.f VERSTOOLS.f WALTALBE.f WALTBRDF.f WATE.f WAVA1.f WAVA2.f \
        WAVA3.f WAVA4.f WAVA5.f WAVA6.f AEROPROF.f main.f
F_OBJ = $(F_SRC:.f=.o)

# Define include paths
INCDIR  = -I.
NCFLAGS = $(EXTRA) $(INCDIR)

# Define the object libraries and paths
MATHLIB = -lm
LOADLIB = $(MATHLIB)

# Define C executables
EXE = sixsV1.0B

#-----------------------------------------------------------------------------
all: $(EXE)

$(EXE): $(F_OBJ)
	$(FC) $(EXTRA) $(F_OBJ) -o $(EXE) $(LOADLIB)

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
$(F_OBJ): $(F_SRC)

.f.o:
	gfortran $(NCFLAGS) -c $< -o $@
