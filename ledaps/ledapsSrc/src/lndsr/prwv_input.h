/*
!C****************************************************************************

!File: input.h

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

 ! Design Notes:
   1. Structure is declared for the 'input' data type.
  
!END****************************************************************************
*/

#ifndef PRWV_INPUT_H
#define PRWV_INPUT_H

#include <stdlib.h>
#include <stdio.h>
#include "lndsr.h"
#include "myhdf.h"
#include "const.h"
#include "date.h"
#include "read_grib_tools.h"

typedef struct {
  int ntime;               /* number of time samples */
  int nlat;                /* number of lat samples */
  int nlon;                /* number of lon samples */
} Grid_coord_int_t;

typedef struct {
  int year;                /* year */
  int doy;                 /* Day of year */
  char platform[256];      /* Platform for the ozone data */
  Sat_t sat;               /* Satellite */
  Inst_t inst;             /* Instrument */
  Date_t acq_date;         /* Acqsition date/time (scene center) */
  Date_t prod_date;        /* Production date (must be available for ETM) */
  float sun_zen;           /* Solar zenith angle (radians; scene center) */
  float sun_az;            /* Solar azimuth angle (radians; scene center) */
  Wrs_t wrs_sys;	   /* WRS system */
  int ipath;               /* WRS path number */
  int irow;                /* WRS row number */
  int fill;                /* Fill value */
  int iband[NBAND_PRWV_MAX];
                           /* Band numbers */
} InputPrwv_meta_t;

/* Structure for the 'input' data type */

typedef struct {
  char *file_name;         /* Input image file name */
  bool open;               /* Open file flag; open = true */
  InputPrwv_meta_t meta;   /* Input metadata */
  int nband;               /* Number of input image files (bands) */
  Grid_coord_int_t size;   /* Input file size */
  int32 sds_file_id;       /* SDS file id */
  Myhdf_sds_t sds[NBAND_PRWV_MAX];
                           /* SDS data structures */
  float *buf[NBAND_PRWV_MAX];
                           /* Input data buffer (one line of data) */
  float scale_factor[NBAND_PRWV_MAX];
  float add_offset[NBAND_PRWV_MAX];
} InputPrwv_t;

typedef struct {
  char *file_name;         /* Input image file name */
  bool open;               /* Open file flag; open = true */
  InputPrwv_meta_t meta;   /* Input metadata */
  int nband;               /* Number of input image files (bands) */
  Grid_coord_int_t size;   /* Input file size */
  int32 sds_file_id;       /* SDS file id */
  Myhdf_sds_t sds[NBAND_PRWV_MAX];
                           /* SDS data structures */
  int16 *buf[NBAND_PRWV_MAX];
                           /* Input data buffer (one line of data) */
  float scale_factor[NBAND_PRWV_MAX];
  float add_offset[NBAND_PRWV_MAX];
} InputOzon_t;

/* Prototypes */

InputPrwv_t *OpenInputPrwv(char *file_name);
bool GetInputPrwv(InputPrwv_t *this, int iband, float *read_buffer);
bool CloseInputPrwv(InputPrwv_t *this);
void FreeInputPrwv(InputPrwv_t *this);
bool GetInputPrwvMeta(InputPrwv_t *this);
int get_prwv_anc(t_ncep_ancillary *anc,InputPrwv_t *this, float *data, int index);

InputOzon_t *OpenInputOzon(char *file_name);
bool GetInputOzon(InputOzon_t *this, int iband, int *read_buffer);
bool CloseInputOzon(InputOzon_t *this);
void FreeInputOzon(InputOzon_t *this);
bool GetInputOzonMeta(InputOzon_t *this);
int get_ozon_anc(t_ncep_ancillary *anc,InputOzon_t *this, int*data, int index);

#endif
