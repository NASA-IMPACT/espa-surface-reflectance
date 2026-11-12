#ifndef OUTPUT_H
#define OUTPUT_H

#include "common.h"
#include "input.h"

/* Define some of the constants to use in the output data products */
#define FILL_VALUE 0
#define CLOUD_FILL_VALUE 0
#define FILL_VALUE_AERO -9999

/* written to XML file for users of the SR data */
/*#define SCALE_FACTOR_REFL 0.0000275*/
#define SCALE_FACTOR_REFL 0.0001
/*#define OFFSET_REFL -0.20*/
#define OFFSET_REFL 0.0
#define SCALE_FACTOR_TH 0.00341802
#define OFFSET_TH 149.0
#define SCALE_FACTOR_AERO 0.001
#define OFFSET_AERO 0.0

/* applied to the SR data before writing to the output file */
#define MULT_FACTOR_REFL (1.0 / SCALE_FACTOR_REFL)
#define BAND_OFFSET_REFL (-OFFSET_REFL)
#define MULT_FACTOR_TH (1.0 / SCALE_FACTOR_TH)
#define BAND_OFFSET_TH (-OFFSET_TH)
#define MULT_FACTOR_AERO (1.0 / SCALE_FACTOR_AERO)
#define BAND_OFFSET_AERO (-OFFSET_AERO)

/* min/max valid values (unscaled) */
#define MIN_VALID_REFL -0.2
#define MAX_VALID_REFL 1.60
#define MIN_VALID_TH 150.0
#define MAX_VALID_TH 350.0
#define MIN_VALID_AERO 0.0
#define MAX_VALID_AERO 5.0

/* Define the output product types */
typedef enum {OUTPUT_TOA=0, OUTPUT_SR=1} Myoutput_t;

/* Structure for the 'output' data type */
typedef struct {
  bool open;            /* Flag to indicate whether output file is open;
                           'true' = open, 'false' = not open */
  Inst_t inst;          /* instrument */
  int nband;            /* Number of output bands */
  int nlines;           /* Number of output lines */
  int nsamps;           /* Number of output samples */
  Espa_internal_meta_t metadata;  /* Metadata container to hold the band
                           metadata for the output bands; global metadata
                           won't be valid */
  FILE *fp_bin[NBAND_TTL_OUT];  /* File pointer for binary files; see common.h
                           for the bands and order of bands in the output */
} Output_t;

/* Prototypes */
Output_t *open_output
(
    Espa_internal_meta_t *in_meta,  /* I: input metadata structure */
    Input_t *input,                 /* I: input band data structure */
    Myoutput_t output_type,         /* I: are we processing TOA, SR outputs? */
    char *aux_name                  /* I: auxiliary filename (for SR) */
);

int close_output
(
    Sat_t sat,              /* I: satellite */
    Output_t *output,       /* I/O: Output data structure to close */
    Myoutput_t output_type  /* I: are we processing TOA, SR outputs? */
);

int free_output
(
    Output_t *output,       /* I/O: Output data structure to free */
    Myoutput_t output_type  /* I: are we processing TOA, SR outputs? */
);

int put_output_lines
(
    Output_t *output,  /* I: Output data structure; buf contains the line to
                             be written */
    void *buf,         /* I: buffer to be written */
    int iband,         /* I: current band to be written (0-based) */
    int iline,         /* I: current line to be written (0-based) */
    int nlines,        /* I: number of lines to be written */
    int nbytes         /* I: number of bytes per pixel in this band */
);

int get_output_lines
(
    Output_t *output, /* I: pointer to output data structure */
    int iband,        /* I: current band to read (0-based) */
    int iline,        /* I: current line to read (0-based) */
    int nlines,       /* I: number of lines to read */
    int nbytes,       /* I: number of bytes per pixel in this band */
    void *buf         /* I: pointer to the buffer to be returned */
);

char *upper_case_str
(
    char *str    /* I: string to be converted to upper case */
);

void convert_output
(
    float **sband,      /* I: unscaled SR or TOA bands */
    int band,           /* I: band number to convert */
    int nlines,         /* I: number of lines */
    int nsamps,         /* I: number of samples */
    bool thermal,       /* I: flag to specify if processing a thermal band,
                              for correct scale/offset */
    uint16 *out_band    /* O: scaled output for the processed band */
);

void convert_aerosol_output
(
    float *aero,        /* I: unscaled aerosol band */
    uint16 *qaband,     /* I: QA band for the input image, nlines x nsamps */
    int nlines,         /* I: number of lines */
    int nsamps,         /* I: number of samples */
    int16 *out_band     /* O: scaled output for the aerosol band */
);

#endif
