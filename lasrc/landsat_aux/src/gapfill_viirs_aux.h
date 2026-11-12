#ifndef _GAPFILL_VIIRS_AUX_H_
#define _GAPFILL_VIIRS_AUX_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <libgen.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include "hdf5.h"
#include "error_handler.h"

#define MAXLENGTH 128
#define MAXLENGTH2 5000

/* Define the input fill value */
#define VIIRS_FILL 0
#define IFILL -1

/* Prototypes */
int get_args
(
    int argc,              /* I: number of cmd-line args */
    char *argv[],          /* I: string of cmd-line args */
    int *month,            /* O: month (1-12) of auxiliary file to be proc'd */
    int *day,              /* O: day (1-31) of auxiliary file to be proc'd */
    int *year,             /* O: year of auxiliary file to be processed */
    char **viirs_aux_file  /* O: address of input VIIRS auxiliary file */
);

void usage();

int open_oz_wv_datasets
(
    char *filename,       /* I: VIIRS file to be read */
    hid_t *file_id,       /* O: VIIRS file id */
    hid_t *ozone_dsid,    /* O: ozone dataset ID */
    hid_t *wv_dsid        /* O: water vapor dataset ID */
);

#endif
