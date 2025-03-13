## LaSRC Version 3.5.2 Release Notes
Release Date: March, 2025

This release of the LaSRC code has been modified by the NASA-IMPACT team as part of
maintenance for the HLS project.

### Fixed

- Support Sentinel-2C and -2D by updating hard coded platform values
- Include `PROC_ALL_BANDS` definition
- Mask a pixel if _any_ band has invalid data [#17](https://github.com/NASA-IMPACT/espa-surface-reflectance/pull/17)

## LaSRC Version 3.5.1 Release Notes
Release Date: TBD 2023

### Downloads
LaSRC (Landsat Surface Reflectance Code) source code

    git clone https://eroslab.cr.usgs.gov/lsrd/espa-surface-reflectance.git

LaSRC auxiliary files

    http://edclpdsftp.cr.usgs.gov/downloads/auxiliaries/lasrc_auxiliary/lasrc_aux.2013-2017.tar.gz
    http://edclpdsftp.cr.usgs.gov/downloads/auxiliaries/lasrc_auxiliary/MSILUT.tar.gz

See git tag [version_3.5.1]

### Installation
  * Install dependent libraries - ESPA product formatter (https://eroslab.cr.usgs.gov/lsrd/espa-product-formatter.git)
  * Set up environment variables.  Can create an environment shell file or add the following to your bash shell.  For C shell, use 'setenv VAR "directory"'.
```
    export PREFIX="path_to_directory_for_lasrc_build_data"
```

  * Install baseline auxiliary files and set up the environment variables.
```
    tar -xvzf lasrc_auxiliary.2013-2017.tar.gz
    export LASRC_AUX_DIR="directory_saved_auxiliary_files"
    (or in c shell use 
    setenv LASRC_AUX_DIR "directory_saved_auxiliary_files")
```

  * Install the MSILUT files for the Sentinel-2 processing. Untar the MSILUT.tar.gz file into $LASRC\_AUX\_DIR.
```
    tar -xvzf MSILUT.tar.gz
```

  * Download (from Github USGS-EROS surface-reflectance project) and install source files. The following build will create a list of executable files under $PREFIX/bin (tested in Linux with the gcc compiler). It will also copy various scripts from the scripts directory up the the $PREFIX/bin directory.
```
    cd lasrc\c_version\src
    make
    make install
    make clean

    cd lasrc\landsat_aux\src
    make
    make install
    make clean
```

  * Test - Download Landsat Level 1 files.  Run the do\_lasrc Python script in the PREFIX bin directory to run the applications.  Use do\_lasrc.py --help for the usage information.  This script requires that your LaSRC binaries are in your $PATH or that you have a $BIN environment variable set up to point to the PREFIX bin directory.
```
    convert_lpgs_to_espa --mtl <Landsat_MTL_file>
    do_lasrc.py --xml <Landsat_ESPA_XML_file>
```

  * Check output
```
    {scene_name}_toa_*: top-of-atmosphere (TOA) reflectance (or brightness temperature for thermal bands) in internal ESPA file format
    {scene_name}_sr_*: surface reflectance in internal ESPA file format
```

  * Note that the FORTRAN version contains the code as delivered from Eric Vermote and team at NASA Goddard Space Flight Center.  The C version contains the converted FORTRAN code into C to work in the ESPA environment.  It also contains any bug fixes, agreed upon by Eric's team, along with performance enhancements.  The FORTRAN code contains debugging and validation code, which is not needed for production processing.

### Dependencies
  * Python >= 3.6.X
  * ESPA raw binary and ESPA common libraries from ESPA product formatter and associated dependencies
  * XML2 library
  * Auxiliary data products
    1. LAADS Terra and Aqua CMG and CMA data
       OR VIIRS JPSS-1 or SUOMI NPP data
    2. CMGDEM HDF file
    3. Various input files and model information provided with the LaSRC auxiliary .tar.gz file

### Auxiliary Data Updates
The baseline auxiliary files provided don't include the daily climate data.  In order to generate or update the auxiliary files to the most recent day of year (actually the most current auxiliary files available will be 2-3 days prior to the current day of year do to the latency of the underlying LAADS products) the user will want to run the updatelads.py and/or updatelads_modis.py script(s) available in $PREFIX/bin.  This script can be run with the "--help" argument to print the usage information.  In general the --quarterly argument will reprocess/update all the LAADS data backwards in time.  This is good to do every once in a while to make sure any updates to the LAADS data products are captured.  The --today command-line argument will process the LAADS data for the most recent year.  In general, it is suggested to run the script with --quarterly once a quarter.  Then run the script with --today on a nightly basis.  This should provide an up-to-date version of the auxiliary input data for LaSRC.  The easiest way to accomplish this is to set up a nightly and quarterly cron job.

The generate_monthly_climatology script will generate the monthly climatology that is used for filling the gaps in the daily LAADS VIIRS products.  This script should be run nightly starting on the first of every month and going through the fifth of the month.  This will allow the monthly averages for the previous month to be generated and available, even if the previous month isn't complete with all the auxiliary products due to time lags.  By the fifth of the month, all the auxiliary products from the previous month should be availble (under most conditions) and therefore the final monthly average generated should contain data for the entire month.

NASA GSFC has deprecated the use of ftp access to their data.  The updatelads script uses a public https interface.  Scripted downloads need to use LAADS app keys in order to be properly authorized.  Non-ESPA environments will need to obtain an app key for the LAADS DAAC and add that key to the 'TOKEN' variable at the top of the updatelads.py script.  ESPA environments will be able to run the updatelads.py script as-is.

### Data Preprocessing
This version of the LaSRC application requires the input Landsat products to be in the ESPA internal file format.  After compiling the product formatter raw\_binary libraries and tools, the convert\_lpgs\_to\_espa command-line tool can be used to create the ESPA internal file format for input to the LaSRC application.

### Data Postprocessing
After compiling the product-formatter raw\_binary libraries and tools, the convert\_espa\_to\_gtif and convert\_espa\_to\_hdf command-line tools can be used to convert the ESPA internal file format to HDF or GeoTIFF.  Otherwise the data will remain in the ESPA internal file format, which includes each band in the ENVI file format (i.e. raw binary file with associated ENVI header file) and an overall XML metadata file.

### Verification Data

### User Manual

### Product Guide

## Release Notes
1. Fixed a multi-threading bug in the Landsat aerosol processing.
2. Updated the Sentinel-2 print information for the use of the VIIRS default
   auxiliary starting date (to be consistent with Landsat messaging).
3. Added print statements at the beginning of the main applications to indicate
   whether the code was built with multi-threading or not.
