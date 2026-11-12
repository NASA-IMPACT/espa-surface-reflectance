#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing L8 SR.
#-----------------------------------------------------------------------------
.PHONY: all install clean

all:
	echo "make all in src..."; \
        ($(MAKE) -C src || exit 1)
	echo "make all in src_modis..."; \
        ($(MAKE) -C src_modis || exit 1)

install: all
	echo "make install in scripts..."; \
        ($(MAKE) -C scripts install || exit 1)
	echo "make install in src..."; \
        ($(MAKE) -C src install || exit 1)
	echo "make install in src_modis..."; \
        ($(MAKE) -C src_modis install || exit 1)

clean:
	echo "make clean in src..."; \
        ($(MAKE) -C src clean || exit 1)
	echo "make clean in src_modis..."; \
        ($(MAKE) -C src_modis clean || exit 1)

