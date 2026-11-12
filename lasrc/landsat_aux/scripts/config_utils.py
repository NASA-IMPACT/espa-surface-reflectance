## This file was pulled from the espa-processing project to be used with
## the LAADS processing.

import os
from configparser import ConfigParser

def get_cfg_file_path(filename):
    """Build the full path to the config file
    Args:
        filename (str): The name of the file to append to the full path.
    Raises:
        Exception(message)
    """

    # Use the users home directory as the base source directory for
    # configuration
    if 'HOME' not in os.environ:
        raise Exception('[HOME] not found in environment')
    home_dir = os.environ.get('HOME')

    # Build the full path to the configuration file
    config_path = os.path.join(home_dir, '.usgs', 'espa', filename)

    return config_path


def retrieve_cfg(cfg_filename):
    """Retrieve the configuration for the cron
    Returns:
        cfg (ConfigParser): Configuration for ESPA cron.
    Raises:
        Exception(message)
    """

    # Build the full path to the configuration file
    config_path = get_cfg_file_path(cfg_filename)

    if not os.path.isfile(config_path):
        raise Exception('Missing configuration file [{}]'.format(config_path))

    # Create the object and load the configuration
    cfg = ConfigParser()
    cfg.read(config_path)

    return cfg
