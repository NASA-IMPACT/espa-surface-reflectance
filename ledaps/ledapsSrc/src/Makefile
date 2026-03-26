#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing Ledaps.
#-----------------------------------------------------------------------------
.PHONY: all install clean

MODULES = lndpm lndcal lndsr 6sV-1.0B

all:
	@for module in $(MODULES); do \
		echo "make all in $$module..."; \
	    ($(MAKE) -C $$module || exit 1); \
	done

install:
	@for module in $(MODULES); do \
		echo "installing all in $$module..."; \
	    ($(MAKE) -C $$module install || exit 1); \
	done

clean:
	@for module in $(MODULES); do \
		echo "make clean in $$module..."; \
	    ($(MAKE) -C $$module clean || exit 1); \
	done
