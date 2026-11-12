#-----------------------------------------------------------------------------
# Makefile
#
# Project Name: surface reflectance
#-----------------------------------------------------------------------------
.PHONY: check-environment all install clean all-script install-script clean-script all-ledaps install-ledaps clean-ledaps all-ledaps-aux install-ledaps-aux clean-ledaps-aux all-lasrc install-lasrc clean-lasrc all-lasrc-aux install-lasrc-aux clean-lasrc-aux all-aux install-aux build tests docs deploy

include make.config

DIR_LEDAPS = ledaps
DIR_LaSRC = lasrc

#-----------------------------------------------------------------------------
all: all-script all-ledaps all-lasrc

install: check-environment install-script install-ledaps install-lasrc

clean: clean-script clean-ledaps clean-lasrc

#-----------------------------------------------------------------------------
# ESPA Standard Makefile targets.  Do not remove.
#
build: all

tests:

docs:

deploy: 

#-----------------------------------------------------------------------------
all-script:
	echo "make all in scripts"; \
        (cd scripts; $(MAKE) all);

install-script:
	echo "make install in scripts"; \
        (cd scripts; $(MAKE) install);

clean-script:
	echo "make clean in scripts"; \
        (cd scripts; $(MAKE) clean);

#-----------------------------------------------------------------------------
all-lasrc:
	echo "make all in $(DIR_LaSRC)"; \
        (cd $(DIR_LaSRC); $(MAKE) all-lasrc);

install-lasrc:
	echo "make install in $(DIR_LaSRC)"; \
        (cd $(DIR_LaSRC); $(MAKE) install-lasrc);

clean-lasrc:
	echo "make clean in $(DIR_LaSRC)"; \
        (cd $(DIR_LaSRC); $(MAKE) clean-lasrc);

#-----------------------------------------------------------------------------
all-lasrc-aux:
	echo "make all in $(DIR_LaSRC)"; \
        (cd $(DIR_LaSRC); $(MAKE) all-lasrc-aux);

install-lasrc-aux:
	echo "make install in $(DIR_LaSRC)"; \
        (cd $(DIR_LaSRC); $(MAKE) install-lasrc-aux);

clean-lasrc-aux:
	echo "make clean in $(DIR_LaSRC)"; \
        (cd $(DIR_LaSRC); $(MAKE) clean-lasrc-aux);

#-----------------------------------------------------------------------------
all-ledaps:
	echo "make all in $(DIR_LEDAPS)"; \
        (cd $(DIR_LEDAPS); $(MAKE) all-ledaps);

install-ledaps:
	echo "make install in $(DIR_LEDAPS)"; \
        (cd $(DIR_LEDAPS); $(MAKE) install-ledaps);

clean-ledaps:
	echo "make clean in $(DIR_LEDAPS)"; \
        (cd $(DIR_LEDAPS); $(MAKE) clean-ledaps);

#-----------------------------------------------------------------------------
all-ledaps-aux:
	echo "make all in $(DIR_LEDAPS)"; \
        (cd $(DIR_LEDAPS); $(MAKE) all-ledaps-aux);

install-ledaps-aux:
	echo "make install in $(DIR_LEDAPS)"; \
        (cd $(DIR_LEDAPS); $(MAKE) install-ledaps-aux);

clean-ledaps-aux:
	echo "make clean in $(DIR_LEDAPS)"; \
        (cd $(DIR_LEDAPS); $(MAKE) clean-ledaps-aux);

#-----------------------------------------------------------------------------
all-aux: all-ledaps-aux all-lasrc-aux

install-aux: install-ledaps-aux install-lasrc-aux

#-----------------------------------------------------------------------------
check-environment:
ifndef PREFIX
    $(error Environment variable PREFIX is not defined)
endif

