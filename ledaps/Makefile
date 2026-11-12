#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing Ledaps.
#----------------------------------------------------------------------------
.PHONY: all install clean all-ledaps install-ledaps clean-ledaps all-ledaps-aux install-ledaps-aux clean-ledaps-aux

DIR_LEDAPS = ledapsSrc
DIR_LEDAPS_AUX = ledapsAncSrc

#-----------------------------------------------------------------------------
all: all-ledaps all-ledaps-aux

install: install-ledaps install-ledaps-aux

clean: clean-ledaps clean-ledaps-aux

#-----------------------------------------------------------------------------
all-ledaps:
	echo "make all in $(DIR_LEDAPS)..."; \
        ($(MAKE) -C $(DIR_LEDAPS) || exit 1)

install-ledaps: all-ledaps
	echo "make install in $(DIR_LEDAPS)..."; \
        ($(MAKE) -C $(DIR_LEDAPS) install || exit 1)

clean-ledaps:
	echo "make clean in $(DIR_LEDAPS)..."; \
        ($(MAKE) -C $(DIR_LEDAPS) clean || exit 1)

#-----------------------------------------------------------------------------
all-ledaps-aux:
	echo "make all in $(DIR_LEDAPS_AUX)..."; \
        ($(MAKE) -C $(DIR_LEDAPS_AUX) || exit 1)

install-ledaps-aux: all-ledaps-aux
	echo "make install in $(DIR_LEDAPS_AUX)..."; \
        ($(MAKE) -C $(DIR_LEDAPS_AUX) install || exit 1)

clean-ledaps-aux:
	echo "make clean in $(DIR_LEDAPS_AUX)..."; \
        ($(MAKE) -C $(DIR_LEDAPS_AUX) clean || exit 1)

