#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing Ledaps.
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ../../..
include $(TOP)/make.config

# Define scripts
SCRIPTS = do_ledaps.py

#-----------------------------------------------------------------------------
all:

install:
	install -d $(link_path)
	install -d $(ledaps_bin_install_path)
	@for script in $(SCRIPTS); do \
		cmd="install -m 755 $$script $(ledaps_bin_install_path)"; \
		echo "$$cmd"; $$cmd || exit 1; \
		cmd="ln -sf $(ledaps_link_source_path)/$$script $(link_path)/$$script"; \
		echo "$$cmd"; $$cmd; \
	done;

clean:

