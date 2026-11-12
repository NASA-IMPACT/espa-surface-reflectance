#ifndef _AERO_INTERP_H_
#define _AERO_INTERP_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "lasrc.h"

void aerosol_interp_landsat
(
    Espa_internal_meta_t *xml_metadata, /* I: XML metadata information */
    int aero_window,   /* I: size of the aerosol window */
    int half_aero_window, /* I: size of half the aerosol window */
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps.  It is expected that the ipflag
                               values are computed for the center of the
                               aerosol windows. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                          It is expected that the aerosol values are computed
                          for the center of the aerosol windows.  This routine
                          will fill in the pixels for the remaining, non-center
                          pixels of the window. */
    int nlines,        /* I: number of lines in qaband & taero bands */
    int nsamps         /* I: number of samps in qaband & taero bands */
);

int fix_invalid_aerosols_landsat
(
    uint8 *ipflag,     /* I: QA flag to assist with aerosol interpolation,
                             nlines x nsamps.  It is expected that the ipflag
                             values are computed for the center of the aerosol
                             windows. */
    float *taero,      /* I: aerosol values for each pixel, nlines x nsamps
                             It is expected that the aerosol values are computed
                             for the center of the aerosol windows */
    float *teps,       /* I: angstrom coeff for each pixel, nlines x nsamps
                             It is expected that the eps values are computed
                             for the center of the aerosol windows */
    int aero_window,   /* I: size of the aerosol window (NxN) */
    int half_aero_window, /* I: size of half the aerosol window */
    int nlines,        /* I: number of lines in taero band */
    int nsamps         /* I: number of samps in taero band */
);

void aerosol_interp_sentinel
(
    int aero_window,   /* I: size of the aerosol window */
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps.  It is expected that the ipflag
                               values are computed for the UL of the aerosol
                               windows to match taero. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                          It is expected that the aerosol values are computed
                          for the UL of the aerosol windows.  This routine
                          will fill in the pixels for the remaining pixels of
                          the window. */
    int nlines,        /* I: number of lines in ipflag & taero bands */
    int nsamps         /* I: number of samps in ipflag & taero bands */
);

void ipflag_expand_failed_sentinel
(
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps. It is expected that the ipflag
                               values are computed for all pixels. */
    int nlines,        /* I: number of lines in ipflag band */
    int nsamps         /* I: number of samps in ipflag band */
);

int aero_avg_failed_sentinel
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps. It is expected that the ipflag
                               values are computed for all pixels. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                               updated with average values if needed */
    float *teps,       /* I/O: angstrom coeff for each pixel, nlines x nsamps
                               updated with average values if needed */
    int nlines,        /* I: number of lines in ipflag band */
    int nsamps         /* I: number of samps in ipflag band */
);

#endif
