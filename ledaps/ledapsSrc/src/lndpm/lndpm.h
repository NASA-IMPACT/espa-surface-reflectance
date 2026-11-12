/*****************************************************************************
FILE: lndpm.h
  
PURPOSE: Contains defines and prototypes to read the ESPA XML metadata file
and generate the metadata and text files for downstream LEDAPS processing.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/16/2014    Gail Schmidt     Original development

NOTES:
*****************************************************************************/
#ifndef LNDPM_H
#define LNDPM_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <getopt.h>
#include "espa_metadata.h"
#include "parse_metadata.h"
#include "error_handler.h"
#include "envi_header.h"
#include "gctp_defines.h"

/* Defines */
/* LEDAPS VERSION definitions */
#define LEDAPS_VERSION "5.1.1 (Collection 2)"

/* define useful constants */
#define MAX_STRING_LENGTH 1000
#define DIR_BUF_SIZE (MAXNAMLEN*20)
#define LEDAPS_MSS 1
#define LEDAPS_TM 2
#define LEDAPS_ETM 3
#define MAX_N_BANDS 7
#define N_BANDS_MSS 4

/* define the list of output metadata PARAMETER IDs */
#define HEADER_FILE          0
#define FILE_TYPE            1
#define DATA_PROVIDER        2
#define SATELLITE            3
#define INSTRUMENT           4
#define ACQUISITION_DATE     5
#define ACQUISITION_TIME     6
#define PRODUCTION_DATE      7
#define SOLAR_ZENITH         8
#define SOLAR_AZIMUTH        9
#define WRS_SYSTEM          10
#define WRS_PATH            11
#define WRS_ROW             12
#define NBAND               13
#define BANDS               14
#define GAIN_SETTINGS       15
#define NSAMPLE             16
#define NLINE               17
#define FILE_NAMES          18
#define NBAND_TH            19
#define BANDS_TH            20
#define GAIN_SETTINGS_TH    21
#define NSAMPLE_TH          22
#define NLINE_TH            23
#define FILE_NAMES_TH       24
#define PROJECTION_NUMBER   25
#define PIXEL_SIZE          26
#define UPPER_LEFT_CORNER   27
#define PROJECTION_ZONE     28
#define PROJECTION_SPHERE   29
#define GAIN                30
#define BIAS                31
#define GAIN_TH             32
#define BIAS_TH             33
#define VERT_LON_FROM_POLE  34
#define TRUE_SCALE_LAT      35
#define FALSE_EASTING       36
#define FALSE_NORTHING      37
#define PROJECTION_PARAMS   38
#define FILE_END            39
#define MAX_CM_PARAMS       (FILE_END + 1)
#define NMETA               MAX_CM_PARAMS

typedef struct metadata
{
  char name[MAX_STRING_LENGTH];
  char val[MAX_STRING_LENGTH];
} METADATA;

#endif
