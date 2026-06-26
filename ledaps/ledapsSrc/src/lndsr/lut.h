/*
!C****************************************************************************

!File: lut.h

!Description: Header file for 'input.c' - see 'input.c' for more information.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

!Team Unique Header:
  This software was developed by the MODIS Land Science Team Support 
  Group for the Laboratory for Terrestrial Physics (Code 922) at the 
  National Aeronautics and Space Administration, Goddard Space Flight 
  Center, under NASA Task 92-012-00.

 ! References and Credits:

  ! MODIS Science Team Member:
      Christopher O. Justice
      MODIS Land Science Team           University of Maryland
      justice@hermes.geog.umd.edu       Dept. of Geography
      phone: 301-405-1600               1113 LeFrak Hall
                                        College Park, MD, 20742

  ! Developers:
      Robert E. Wolfe (Code 922)
      MODIS Land Team Support Group     Raytheon ITSS
      robert.e.wolfe.1@gsfc.nasa.gov    4400 Forbes Blvd.
      phone: 301-614-5508               Lanham, MD 20770  

 ! Updates:
   09/19/2012 Gail Schmidt, USGS EROS
   Removed the cloud fill, water, land, shadow, snow, cloud values since the
   ACCA product is no longer used

 ! Design Notes:
   1. Structure is declared for the 'input' data type.
  
!END****************************************************************************
*/

#ifndef LUT_H
#define LUT_H

#include <espa_geoloc.h>
#include "input.h"
#include <stdbool.h>

#define NBAND_SR_LUT (3)

/* Structure for the 'lut' data type */

typedef struct {
  char *file_name;         /* Lookup table file name */
  int nband;               /* Number of bands */
  int output_fill;         /* Output fill value */
  int output_fill_opacity; /* Output fill value for opacity band */
  int aerosol_fill;        /* Aerosol fill value */
  int in_satu;             /* Input saturated value */
  int output_satu;         /* Output saturation value (Feng, 3/23/09) */
  Img_coord_int_t ar_region_size;  
                           /* Size of the aerosol retreval region */
  Img_coord_int_t ar_size;  
                           /* Size of the aerosol retreval image */
  int min_valid_sr;        /* Minimum valid surface reflectance */
  int max_valid_sr;        /* Maximum valid surface reflectance */
  int min_valid_opacity;   /* Minimum valid surface reflectance for opacity band */
  int max_valid_opacity;   /* Maximum valid surface reflectance for opacity band */
  Input_meta_t meta;       /* Input metadata */
  char* long_name_prefix;  /* long name prefix (append band num) */
  char* units;             /* units */
  double scale_factor;     /* scale factor */
  double mult_factor;      /* multiplication factor */
  double atmos_opacity_scale_factor;  /* Atmospheric opacity scale factor */
  double add_offset;       /* add offset */
  double b6_scale_factor;
  double b6_add_offset;
} Lut_t;

/* Prototypes */

Lut_t *GetLut(int nband, Input_meta_t *input_meta, Input_meta_t *b6_meta, 
    Img_coord_int_t *input_size);
void FreeLut(Lut_t *this);

#endif
