#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing L8 SR.
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ../../..
include $(TOP)/make.config

# Define scripts
SCRIPTS = updatelads.py generate_monthly_climatology.py updatelads_modis.py

# Installed but only imported by scripts
SCRIPT_IMPORTS = \
    api_interface.py \
    config_utils.py \
    download_lads.py

#-----------------------------------------------------------------------------
all:

install:
	install -d $(link_path)
	install -d $(lasrc_bin_install_path)

	@for script in $(SCRIPTS); do \
		cmd="install -m 755 $$script $(lasrc_bin_install_path)"; \
		echo "$$cmd"; $$cmd || exit 1; \
		cmd="ln -sf $(lasrc_link_source_path)/$$script $(link_path)/$$script"; \
		echo "$$cmd"; $$cmd; \
	done;

	@for script in $(SCRIPT_IMPORTS); do \
		cmd="install -m 644 $$script $(lasrc_bin_install_path)"; \
		echo "$$cmd"; $$cmd || exit 1; \
	done;

clean:

