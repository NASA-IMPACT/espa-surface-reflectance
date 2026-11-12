#ifndef UTMTODEG_H
#define UTMTODEG_H

#include "espa_geoloc.h"

/* Prototypes */
void utmtodeg
(
    Space_def_t *space_def,  /* I: space definition structure */
    int line,       /* I: image line (zero-based) */
    int samp,       /* I: image sample (zero-based) */
    float *lat,     /* O: latitude for line, sample in the image */
    float *lon      /* O: longitude for line, sample in the image */
);

#endif
