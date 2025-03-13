#! /usr/bin/env python

'''
    PURPOSE: Determine which executable to run and then pass all arguments
             through to the appropriate script.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3

    AUTHOR: ngenetzky@usgs.gov

    NOTES:
        This script is intended to be ran in the directory containing the data.
        This script does not have its own help message and will just return the
            help from underlying executables where appropriate.
        If this script has a required argument then only the usage for that
            argument will be shown if that argument is not included.
        All output from the underlying script will be given to the logger as an
            info message.
'''

import os
import sys
import logging
import argparse
import subprocess


class ExecuteError(Exception):
    '''Raised when command in execute_cmd returns with error'''

    def __init__(self, message, *args):
        self.message = message
        Exception.__init__(self, message, *args)


def execute_cmd(cmd_string):
    '''Execute a command line and return the terminal output

    Raises:
        ExecuteError (Stdout/Stderr)

    Returns:
        output:The stdout and/or stderr from the executed command.
    '''

    (status, output) = subprocess.getstatusoutput(cmd_string)

    if status < 0:
        message = ('Application terminated by signal [{0}]'
                   .format(cmd_string))
        if len(output) > 0:
            message = ' Stdout/Stderr is: '.join([message, output])
        raise ExecuteError(message)

    if status != 0:
        message = 'Application failed to execute [{0}]'.format(cmd_string)
        if len(output) > 0:
            message = ' Stdout/Stderr is: '.join([message, output])
        raise ExecuteError(message)

    if os.WEXITSTATUS(status) != 0:
        message = ('Application [{0}] returned error code [{1}]'
                   .format(cmd_string, os.WEXITSTATUS(status)))
        if len(output) > 0:
            message = ' Stdout/Stderr is: '.join([message, output])
        raise ExecuteError(message)

    return output


def parse_cmd_line():
    '''Will only parse --xml XML_FILENAME from cmdline.

    Precondition:
        '--xml FILENAME' exists in command line arguments
    Postcondition:
        returns xml_filename

    Note: Help is not included because the program will return
          the help from the underlying program.
    '''

    # Try to parse out the XML so the application can be determined
    parse_xml = argparse.ArgumentParser(add_help=False)
    parse_xml.add_argument('--xml', action='store',
                           dest='xml_filename', required=True,
                           help='Input XML metadata file',
                           metavar='FILE')
    (temp, extra_args) = parse_xml.parse_known_args()

    return temp.xml_filename


def get_science_application_name(satellite_sensor_code):
    '''Returns name of executable that needs to be called'''

    s2_prefixes_collection = ['S2A_', 'S2B_', 'S2C_', 'S2D_']
    l89_prefixes_collection = ['LC08', 'LO08', 'LT08', 'LC09', 'LO09', 'LT09']
    other_prefixes_collection = ['LT04', 'LT05', 'LE07']

    if satellite_sensor_code in l89_prefixes_collection:
        return 'do_lasrc_landsat.py'
    elif satellite_sensor_code in s2_prefixes_collection:
        return 'do_lasrc_sentinel.py'
    elif satellite_sensor_code in other_prefixes_collection:
        return 'do_ledaps.py'
    else:
        raise Exception('Satellite-Sensor code ({0}) not understood'
                        .format(satellite_sensor_code))


def main():
    '''Determines executable, and calls it with all input arguments '''
    # Determine the logging level. Default is INFO.
    espa_log_level = os.environ.get('ESPA_LOG_LEVEL')
    if espa_log_level == 'DEBUG':
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    # Setup the default logger format and level.  Log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=log_level,
                        stream=sys.stdout)

    # Get the logger
    logger = logging.getLogger(__name__)

    xml_filename = parse_cmd_line()
    satellite_sensor_code = os.path.basename(xml_filename)[0:4]

    # Get the science application
    cmd = [get_science_application_name(satellite_sensor_code)]

    # Pass all arguments through to the science application
    cmd.extend(sys.argv[1:])

    # Convert the list to a string
    cmd_string = ' '.join(cmd)
    try:
        logger.info(' '.join(['EXECUTING SCIENCE APPLICATION:',
                              cmd_string]))
        output = execute_cmd(cmd_string)

        if len(output) > 0:
            logger.info('\n{0}'.format(output))
    except ExecuteError:
        logger.exception('Error running {0}.'
                         'Processing will terminate.'
                         .format(os.path.basename(__file__)))
        raise  # Re-raise so exception message will be shown.

if __name__ == '__main__':
    main()
