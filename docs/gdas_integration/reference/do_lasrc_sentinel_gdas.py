#! /usr/bin/env python
import sys
import os
import re
import subprocess
import datetime
from optparse import OptionParser
import logging
import glob     # list and manipulate filenames

ERROR = 1
SUCCESS = 0
VIIRS_AUX_STARTING_DATE = "20231001"
GRIB_STARTING_DATE = "20260501"
AUX_MODE_ENV_VAR = "LASRC_AUX_MODE"
DEFAULT_LASRC_BIN = (
    "/usr/local/bin/lasrc"
)


def _parse_timestamp(value):
    value = value.strip().replace("Z", "")
    for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S"):
        try:
            return datetime.datetime.strptime(value, fmt)
        except ValueError:
            pass
    return None


def _extract_scene_datetime(xml_infile, workdir, default_date):
    safe_name = os.path.basename(workdir)
    match = re.search(r"_(\d{8}T\d{6})_", safe_name)
    if match:
        return datetime.datetime.strptime(match.group(1), "%Y%m%dT%H%M%S")

    patterns = [
        r"<PRODUCT_START_TIME>([^<]+)</PRODUCT_START_TIME>",
        r"<SENSING_TIME>([^<]+)</SENSING_TIME>",
        r"<Acquisition_Date>([^<]+)</Acquisition_Date>",
    ]
    try:
        with open(xml_infile, "r") as xml_fd:
            xml_text = xml_fd.read()
        for pattern in patterns:
            match = re.search(pattern, xml_text)
            if match:
                parsed = _parse_timestamp(match.group(1))
                if parsed is not None:
                    return parsed
    except OSError:
        pass

    return datetime.datetime.strptime(default_date, "%Y%m%d")


def _parse_gdas_time(path):
    name = os.path.basename(path)
    match = re.search(r"(\d{6})\.(\d{2})z$", name)
    if match is None:
        return None
    return datetime.datetime.strptime(
        match.group(1) + match.group(2), "%y%m%d%H"
    )


def _choose_grib_pair(scene_dt, full_aux_dir):
    aux_files = sorted(glob.glob("{}/gdas*".format(full_aux_dir)))
    timed_files = []
    for path in aux_files:
        parsed = _parse_gdas_time(path)
        if parsed is not None:
            timed_files.append((parsed, path))

    if len(timed_files) < 2:
        return None, None, aux_files

    timed_files.sort(key=lambda item: item[0])
    before = [item for item in timed_files if item[0] <= scene_dt]
    after = [item for item in timed_files if item[0] >= scene_dt]

    if before and after:
        first = before[-1]
        second = after[0]
        if first[1] != second[1]:
            return os.path.basename(first[1]), os.path.basename(second[1]), aux_files

        first_idx = timed_files.index(first)
        if first_idx > 0:
            return (
                os.path.basename(timed_files[first_idx - 1][1]),
                os.path.basename(first[1]),
                aux_files,
            )
        if first_idx + 1 < len(timed_files):
            return (
                os.path.basename(first[1]),
                os.path.basename(timed_files[first_idx + 1][1]),
                aux_files,
            )

    if len(timed_files) >= 2:
        return (
            os.path.basename(timed_files[-2][1]),
            os.path.basename(timed_files[-1][1]),
            aux_files,
        )

    return None, None, aux_files

#############################################################################
# Created on August 23, 2019 by Gail Schmidt, USGS/EROS
# Created Python script to run the Sentinel surface reflectance code based
# on the inputs specified by the user.  This script will determine the input
# auxiliary file needed for processing, based on the date of the Sentinel
# input file.
#
# Usage: do_lasrc_sentinel.py --help prints the help message
############################################################################
class SurfaceReflectance():

    def __init__(self):
        pass


    ########################################################################
    # Description: runSr will use the parameters passed for the input and
    # output files.  If input/output files are None (i.e. not specified) then
    # the command-line parameters will be parsed for this information.  The
    # surface reflectance application is then executed to generate the desired
    # outputs on the specified input file.  If a log file was specified, then
    # the output from this application will be logged to that file.
    #
    # Inputs:
    #   xml_infile - name of the input XML file
    #   viirs_aux_starting_date - YYYYmmdd string to identify the date when
    #     VIIRS auxiliary products should start being used. Default is to
    #     pull the starting date from the VIIRS_AUX_STARTING_DATE environment
    #     variable. If that isn't defined, then the hardcoded date will be used.
    #
    # Returns:
    #     ERROR - error running the surface reflectance application
    #     SUCCESS - successful processing
    #
    # Notes:
    #   1. The script obtains the path of the XML file and changes
    #      directory to that path for running the surface reflectance
    #      application.  If the XML file directory is not writable, then this
    #      script exits with an error.
    #   2. If the XML file is not specified and the information is
    #      going to be grabbed from the command line, then it's assumed all
    #      the parameters will be pulled from the command line.
    #######################################################################
    def runSr (self, xml_infile=None, viirs_aux_starting_date=None):
        use_grib = False
        # if no parameters were passed then get the info from the
        # command line
        if xml_infile is None:
            # Get version number
            cmdstr = ('lasrc --version')
            (exit_code, self.version) = subprocess.getstatusoutput(cmdstr)

            # get the command line argument for the XML file
            parser = OptionParser(version = self.version)
            parser.add_option ("-i", "--xml", type="string",
                dest="xml",
                help="name of XML file", metavar="FILE")
            parser.add_option ("--viirs_aux_starting_date", type="string",
                dest="viirs_aux_starting_date",
                action='store',
                metavar='YYYYMMDD',
                help="Acquisition date at which to begin using VIIRS "
                     "auxiliary data instead of MODIS data. The default is "
                     "to pull the date from the VIIRS_AUX_STARTING_DATE "
                     "environment variable. If that isn't set, then the "
                     "hardcoded date will be used.")
            (options, args) = parser.parse_args()
    
            # XML input file
            xml_infile = options.xml
            if xml_infile is None:
                parser.error ('missing input XML file command-line argument');
                return ERROR

            # options
            viirs_aux_starting_date = options.viirs_aux_starting_date

        # get the logger
        logger = logging.getLogger(__name__)
        msg = ('Surface reflectance processing of Sentinel-2 file: {}'
               .format(xml_infile))
        logger.info (msg)

        aux_mode = os.environ.get(AUX_MODE_ENV_VAR, "").strip().lower()
        grib_starting_date = os.environ.get(
            "GRIB_STARTING_DATE", GRIB_STARTING_DATE
        )

        # determine the VIIRS starting date for processing. The user-specified
        # value is priority. Then look at the VIIRS_AUX_STARTING_DATE
        # environment variable. Then use the default/hardcoded date.
        if viirs_aux_starting_date is None:
            logger.debug ('User did not specify the VIIRS aux starting date')
            viirs_aux_starting_date = os.environ.get('VIIRS_AUX_STARTING_DATE')
            if viirs_aux_starting_date is None:
                logger.debug ('VIIRS_AUX_STARTING_DATE environment variable is '
                              'not set. Using default VIIRS starting date.')
                viirs_aux_starting_date = VIIRS_AUX_STARTING_DATE

        msg = ('VIIRS auxiliary processing start date: {}'
               .format(viirs_aux_starting_date))
        logger.info (msg)
        logger.info('GRIB auxiliary processing start date: {}'.format(
            grib_starting_date))
        logger.info('Auxiliary mode override: {}'.format(
            aux_mode if aux_mode else 'auto'))

        # make sure the XML file exists
        if not os.path.isfile(xml_infile):
            msg = ('XML file does not exist or is not accessible: {}'
                   .format(xml_infile))
            logger.error (msg)
            return ERROR

        # use the base XML filename and not the full path.
        base_xmlfile = os.path.basename (xml_infile)
        msg = 'Processing XML file: {}'.format(base_xmlfile)
        logger.info (msg)
        
        # get the path of the XML file and change directory to that location
        # for running this script.  save the current working directory for
        # return to upon error or when processing is complete.  Note: use
        # abspath to handle the case when the filepath is just the filename
        # and doesn't really include a file path (i.e. the current working
        # directory).
        mydir = os.getcwd()
        xmldir = os.path.dirname (os.path.abspath (xml_infile))
        if not os.access(xmldir, os.W_OK):
            msg = ('Path of XML file is not writable: {}. Script needs '
                   'write access to the XML directory.'.format(xmldir))
            logger.error (msg)
            return ERROR
        msg = ('Changing directories for surface reflectance processing: {}'
               .format(xmldir))
        logger.info (msg)
        os.chdir (xmldir)

        # determine the auxiliary directory for the data
        auxdir = os.environ.get('LASRC_AUX_DIR')
        if auxdir is None:
            logger.error('LASRC_AUX_DIR environment variable not set. Exiting.')
            return ERROR

        # pull the date from the XML filename to determine which auxiliary
        # file should be used for input.
        # Example: S2A_MSI_L1C_T10TFR_20180816_20180903.xml uses the
        # VJ[12]04ANC.A2018228.001.*.h5 HDF5 file or VNP04ANC.A2018228.001.*.h5.
        s2_prefixes_collection = ['S2A', 'S2B', 'S2C', 'S2D']
        if base_xmlfile[0:3] in s2_prefixes_collection:
            # Collection naming convention. Pull the year, month, day from the
            # XML filename. It should be the 4th group, separated by
            # underscores. Then convert month, day to DOY.
            aux_date = base_xmlfile.split('_')[4]
            aux_year = aux_date[0:4]
            aux_month = aux_date[4:6]
            aux_day = aux_date[6:8]
            myday = datetime.date(int(aux_year), int(aux_month), int(aux_day))
            aux_doy = myday.strftime("%j")
            scene_dt = _extract_scene_datetime(xml_infile, xmldir, aux_date)

            # Select MODIS or VIIRS atmospheric aux data based on date acquired
            if aux_mode == "modis":
                logger.info('Using MODIS auxiliary due to {}={}'.format(
                    AUX_MODE_ENV_VAR, aux_mode))
                aux_file = 'L8ANC{}{}.hdf_fused'.format(aux_year, aux_doy)
            elif aux_mode == "viirs":
                logger.info('Using VIIRS auxiliary due to {}={}'.format(
                    AUX_MODE_ENV_VAR, aux_mode))
                full_aux_dir = ('{}/LADS/{}'.format(auxdir, aux_year))
                aux_files = glob.glob('{}/V*04ANC.A{}{}.*.h5'
                                      .format(full_aux_dir, aux_year, aux_doy))
                if len(aux_files) == 0:
                    logger.error('No auxiliary files were found to match '
                                 'V*04ANC.A{}{}.*.h5'.format(aux_year, aux_doy))
                    return ERROR
                aux_file = os.path.basename(aux_files[0])
            elif aux_mode == "grib":
                use_grib = True
                logger.info('Using GRIB auxiliary due to {}={}'.format(
                    AUX_MODE_ENV_VAR, aux_mode))
                full_aux_dir = ('{}/LADS'.format(auxdir))
                aux_file1, aux_file2, aux_files = _choose_grib_pair(
                    scene_dt, full_aux_dir
                )
                logger.debug('Found {} GRIB candidates: {}'.format(
                    len(aux_files), aux_files))
                if aux_file1 is None or aux_file2 is None:
                    logger.error('Not enough GRIB auxiliary files were found in '
                                 '{}'.format(full_aux_dir))
                    return ERROR
                logger.info('Selected GRIB auxiliary pair: {}, {}'.format(
                    aux_file1, aux_file2))
            elif aux_date < viirs_aux_starting_date:
                # Use MODIS L8ANC files
                logger.debug('Using MODIS auxiliary')
                aux_file = 'L8ANC{}{}.hdf_fused'.format(aux_year, aux_doy)
            elif aux_date < grib_starting_date:
                # Use VIIRS files
                logger.debug('Using VIIRS auxiliary')
                full_aux_dir = ('{}/LADS/{}'.format(auxdir, aux_year))

                # The auxiliary file could be VJ[12]04ANC or VNP04ANC, however
                # there should only be one, with VJ204ANC being the priority.
                # And LaSRC only wants the base filename, not the entire path.
                aux_files = glob.glob('{}/V*04ANC.A{}{}.*.h5'
                                      .format(full_aux_dir, aux_year, aux_doy))
                if len(aux_files) == 0:
                    logger.error('No auxiliary files were found to match '
                                 'V*04ANC.A{}{}.*.h5'.format(aux_year, aux_doy))
                    return ERROR
                aux_file = os.path.basename(aux_files[0])
            else:
                # use GRIB files. Note GRIB files must be located under dir /LADS/
                use_grib = True
                logger.debug('Using GRIB auxiliary')
                full_aux_dir = ('{}/LADS'.format(auxdir))
                aux_file1, aux_file2, aux_files = _choose_grib_pair(
                    scene_dt, full_aux_dir
                )
                logger.debug('Found {} GRIB candidates: {}'.format(
                    len(aux_files), aux_files))
                if aux_file1 is None or aux_file2 is None:
                    logger.error('Not enough GRIB auxiliary files were found in '
                                 '{}'.format(full_aux_dir))
                    return ERROR
                logger.info('Selected GRIB auxiliary pair: {}, {}'.format(
                    aux_file1, aux_file2))
        else:
            msg = ('Base XML filename is not recognized as a valid Sentinel-2 '
                   'scene name'.format(base_xmlfile))
            logger.error (msg)
            os.chdir (mydir)
            return ERROR

        # run surface reflectance algorithm, checking the return status.  exit
        # if any errors occur.
        lasrc_bin = os.environ.get('LASRC_BIN', DEFAULT_LASRC_BIN)
        if use_grib:
             cmdstr = ('{} --xml={} --aux={} --aux2={} --verbose'
                  .format(lasrc_bin, xml_infile, aux_file1, aux_file2))
        else:
                
            cmdstr = ('{} --xml={} --aux={} --verbose'
                  .format(lasrc_bin, xml_infile, aux_file))
        msg = 'Executing lasrc command: {}'.format(cmdstr)
        logger.debug (msg)
        # Stream LaSRC output so long aerosol stages are visible while running.
        process = subprocess.Popen(
            cmdstr, shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, universal_newlines=True, bufsize=1)
        for line in iter(process.stdout.readline, ''):
            logger.info(line.rstrip())
        process.stdout.close()
        exit_code = process.wait()
        if exit_code != 0:
            logger.error('LaSRC process exit code: {}'.format(exit_code))
            msg = 'Error running lasrc.  Processing will terminate.'
            logger.error (msg)
            os.chdir (mydir)
            return ERROR
        
        # successful completion.  return to the original directory.
        os.chdir (mydir)
        msg = 'Completion of surface reflectance.'
        logger.info (msg)
        return SUCCESS

######end of SurfaceReflectance class######

if __name__ == "__main__":
    # determine the logging level. Default is INFO.
    espa_log_level = os.environ.get('ESPA_LOG_LEVEL')
    if espa_log_level == 'DEBUG':
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    # setup the default logger format and level. log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=log_level)
    sys.exit (SurfaceReflectance().runSr())
