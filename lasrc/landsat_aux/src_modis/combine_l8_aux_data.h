#ifndef _COMBINE_L8_AUX_DATA_H_
#define _COMBINE_L8_AUX_DATA_H_

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <ctype.h>
#include <libgen.h>
#include <math.h>
#include <stdbool.h>
#include "mfhdf.h"
#include "error_handler.h"

/* Defines and typedefs */
enum {UNSET, TERRA, AQUA, BOTH};

#define MAXLENGTH 128
#define MAXLENGTH2 5000

/* SRC_DIRECTORY is the location of the output files to be written */
#define LAADS_FILL 0
#define FFILL -999.0
#define IFILL -1
#define SRC_DIRECTORY  "./"

typedef struct{
   int32 sd_id;
   int32 sds_id;
   int32 data_type;
   int sds_dims[2];
   void *data;
   char sdsname[100];
} io_param;


/* Prototypes */
int get_args
(
    int argc,               /* I: number of cmd-line args */
    char *argv[],           /* I: string of cmd-line args */
    char **terra_cmg_file,  /* O: address of input Terra CMG file */
    char **aqua_cmg_file,   /* O: address of input Aqua CMG file */
    char **terra_cma_file,  /* O: address of input Terra CMA file */
    char **aqua_cma_file,   /* O: address of input Aqua CMA file */
    char **output_dir,      /* O: address of output directory */
    bool *verbose           /* O: verbose flag */
);

void usage();

int parse_sds_info
(
    char *filename,            /* I: Aqua/Terra file to be read */
    io_param Terra_params[],   /* O: array of structs for Terra params */
    io_param Aqua_params[]     /* O: array of structs for Aqua params */
);

void parse_lgid
(
    char lgid[],            /* I: local granule ID */
    char product_type[],    /* O: MODIS product type */
    char yearday[]          /* O: year/day string */
);

int metareader
(
    int sd_id,               /* I: file ID */
    char *type_of_meta,      /* I: which metadata will be read */
    char *metastring,        /* I: metadata variable to be found */
    int32 *count,            /* O: count of characters in output string */
    char *data               /* O: string returned for metastring */
);

void get_a_line
(
    char text[MAXLENGTH2],  /* I: text string to be read from */
    int *start,             /* I/O: location where to start reading the line;
                                    updated for where to start reading the next
                                    line after the current line is read */
    char *line              /* O: line that was read from the text string */
);

void make_outfile_name
(
    char *yearday_str,      /* I: string containing the year and DOY */
    char *output_dir,       /* I: output directory for the auxiliary prods */
    char outfile[STR_SIZE]  /* O: output filename for the auxiliary products */
);

void copy_param
(
    void *dest,        /* O: destination data array */
    void *source,      /* I: source data array */
    int32 data_type,   /* I: data type */
    int32 offset       /* I: pixel in source/dest data arrays to be copied */
);

void interpolate
(
    int32 data_type,     /* I: data type of the data array */
    void *data,          /* I: data array */
    long lineoffset,     /* I: pixel location for the start of this line */
    int left,            /* I: location in the line of the left pixel */
    int right            /* I: location in the line of the right pixel */
);

#endif
