/*****************************************************************************
FILE: utmtodeg.c
  
PURPOSE: Contains functions for handling the UTM to lat/long degrees conversion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "utmtodeg.h"

/******************************************************************************
MODULE:  utmtodeg

PURPOSE:  Converts a line, sample location in the image to a latitude, longitude
location in geographic lat/long (degrees), using WGS84 spheroid.


RETURN VALUE:
Type  =  N/A

NOTES:
******************************************************************************/
void utmtodeg
(
    Space_def_t *space_def,  /* I: space definition structure */
    int line,       /* I: image line (zero-based) */
    int samp,       /* I: image sample (zero-based) */
    float *lat,     /* O: latitude for line, sample in the image */
    float *lon      /* O: longitude for line, sample in the image */
)
{
    const double sa = 6378137.0;                 /* WGS84 radius of sphere */
    const double inv_flattening = 298.257223563; /* WGS84 inverse flattening */
    const double false_easting = 500000.0;       /* UTM false easting */
    const double false_northing = 10000000.0;    /* UTM false northing for
                                                    southern hemisphere */
    const float scale_fact = 0.9996;       /* UTM scale factor */
    double cos_lat;                        /* cosine of latitude */
    double sqr_cos_lat;                    /* (cosine of latitude) ^ 2 */
    float x, y;              /* x, y proj coordinates for input line, sample */
    float central_meridian;  /* central meridian for the longitude */
    float sb, c;             /* sphere-based variables */
    float e2, e2cuadrada;    /* sphere-based variables */
    float a, a1, a2, v;      /* intermediate variables */
    float j2, j4, j6;        /* intermediate variables */
    float alpha, beta;       /* intermediate variables */
    float gama, delta;       /* intermediate variables */
    float b, bm, nab;        /* intermediate variables */
    float eps, epsi;         /* intermediate variables */
    float ta0;               /* intermediate variables */
    float senoheps;          /* intermediate variables */
    int zone;                /* UTM zone (abs value of) */

    /* Calculate variables from the WGS84 sphere */
    sb = sa - (sa / inv_flattening);
    e2 = sqrt (pow(sa, 2) - pow(sb, 2)) / sb;
    e2cuadrada = pow(e2, 2);
    c = pow(sa, 2) / sb;

    /* Determine the UTM projection x,y location within the image */
    x = space_def->ul_corner.x + (samp * space_def->pixel_size[0]);
    y = space_def->ul_corner.y - (line * space_def->pixel_size[1]);

    /* Subtract the false easting and handle the false northing if this is
       a south zone product */
    x -= false_easting;
    if (space_def->zone < 0)
        y -= false_northing;

    /* Get the central meridian using the absolute value of the zone;
       6 degree zones, 183 degree offset */
    zone = abs(space_def->zone);
    central_meridian = (zone * 6.0) - 183.0;

    /* Determine the latitude and associated variables */
    *lat = y / (6366197.724 * scale_fact);
    cos_lat = cos(*lat);
    sqr_cos_lat = pow(cos_lat, 2);

    /* Compute intermediate variables */
    v = (c / sqrt (1.0 + e2cuadrada * sqr_cos_lat)) * scale_fact;
    a = x / v;
    a1 = sin(2.0 * *lat);
    a2 = a1 * sqr_cos_lat;
    j2 = *lat + (a1 / 2.0);
    j4 = (3.0 * j2 + a2) / 4.0;
    j6 = (5.0 * j4 + (a2 * sqr_cos_lat)) / 3.0;
    alpha = (3.0 / 4.0) * e2cuadrada;
    beta = (5.0 / 3.0) * pow(alpha, 2);
    gama = (35.0 / 27.0) * pow(alpha, 3);
    bm = scale_fact * c * (*lat - alpha*j2 + beta*j4 - gama*j6);
    b = (y - bm) / v;
    epsi = e2cuadrada * pow(a, 2) / 2.0 * sqr_cos_lat;
    eps = a * (1.0 - epsi / 3.0);
    nab = b * (1.0 - epsi) + *lat;
    senoheps = (exp (eps) - exp (-eps)) / 2.0;
    delta = atan(senoheps / cos(nab));
    ta0 = atan(cos(delta) * tan(nab));

    /* Update the latitude and compute the longitude */
    *lon = (delta * (180.0 / M_PI)) + central_meridian;
    *lat = (*lat + (1 + e2cuadrada * sqr_cos_lat - 1.5 * e2cuadrada *
          sin(*lat) * cos_lat * (ta0 - *lat)) * (ta0 - *lat)) * (180.0 / M_PI);

    return;
}
