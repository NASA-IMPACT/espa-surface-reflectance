#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing LaSRC.
#----------------------------------------------------------------------------
.PHONY: all install clean all-lasrc install-lasrc clean-lasrc all-lasrc-aux install-lasrc-aux clean-lasrc-aux

DIR_LaSRC = c_version
DIR_LaSRC_AUX = landsat_aux

#-----------------------------------------------------------------------------
all: all-lasrc all-lasrc-aux

install: install-lasrc install-lasrc-aux

clean: clean-lasrc clean-lasrc-aux

#-----------------------------------------------------------------------------
all-lasrc:
	echo "make all in $(DIR_LaSRC)..."; \
        ($(MAKE) -C $(DIR_LaSRC) || exit 1)

install-lasrc: all-lasrc
	echo "make install in $(DIR_LaSRC)..."; \
        ($(MAKE) -C $(DIR_LaSRC) install || exit 1)

clean-lasrc:
	echo "make clean in $(DIR_LaSRC)..."; \
        ($(MAKE) -C $(DIR_LaSRC) clean || exit 1)

#-----------------------------------------------------------------------------
all-lasrc-aux:
	echo "make all in $(DIR_LaSRC_AUX)..."; \
        ($(MAKE) -C $(DIR_LaSRC_AUX) || exit 1)

install-lasrc-aux: all-lasrc-aux
	echo "make install in $(DIR_LaSRC_AUX)..."; \
        ($(MAKE) -C $(DIR_LaSRC_AUX) install || exit 1)

clean-lasrc-aux:
	echo "make clean in $(DIR_LaSRC_AUX)..."; \
        ($(MAKE) -C $(DIR_LaSRC_AUX) clean || exit 1)

