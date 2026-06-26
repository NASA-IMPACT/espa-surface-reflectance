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

#ifndef INPUT_H
#define INPUT_H

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "lndcal.h"
#include "const.h"
#include "date.h"
#include "param.h"

#define ANGLE_FILL (-999.0)
#define WRS_FILL (-1)
#define GAIN_BIAS_FILL (-999.0)


/* Structure for the 'input metadata' data type */
typedef struct {
  Sat_t sat;               /* Satellite */
  Inst_t inst;             /* Instrument */
  Date_t acq_date;         /* Acquisition date/time (scene center) */
  bool time_fill;          /* Acquisition time fill; true = fill value (0h) */
  Date_t prod_date;        /* Production date (must be available for ETM) */
  float sun_zen;           /* Solar zenith angle (radians; scene center) */
  float sun_az;            /* Solar azimuth angle (radians; scene center) */
  double szen_scale;       /* solar zenith angle scale factor */
  double szen_offset;      /* solar zenith angle offset */
  float earth_sun_dist;    /* Earth-sun distance */
  Wrs_t wrs_sys;           /* WRS system */
  int ipath;               /* WRS path number */
  int irow;                /* WRS row number */
  unsigned char fill;      /* Fill value */
  int iband[NBAND_REFL_MAX]; /* Band numbers */
  int iband_th;            /* Thermal Band number= (6) */
  float rad_gain[NBAND_REFL_MAX]; /* TOA radiance band gain */
  float rad_bias[NBAND_REFL_MAX]; /* TOA radiance band bias */
  float rad_gain_th[2];    /* Thermal TOA radiance band gain
                              [0] is TM band 6  OR
                              [0] is ETM+ band 6L, [1] is ETM+ band 6H */
  float rad_bias_th[2];    /* Thermal TOA radiance band bias
                              [0] is TM band 6  OR
                              [0] is ETM+ band 6L, [1] is ETM+ band 6H */
  bool use_toa_refl_consts;    /* Are the TOA reflectance gain/bias and K1/K2
                                  constants available? Same with earth-sun
                                  distance */
  float refl_gain[NBAND_REFL_MAX]; /* TOA reflectance band gain */
  float refl_bias[NBAND_REFL_MAX]; /* TOA reflectance band bias */
  float k1_const[2];       /* K1 thermal constant
                              [0] is TM band 6  OR
                              [0] is ETM+ band 6L, [1] is ETM+ band 6H */
  float k2_const[2];       /* K2 thermal constant
                              [0] is TM band 6  OR
                              [0] is ETM+ band 6L, [1] is ETM+ band 6H */
} Input_meta_t;

/* Structure for the 'input' data type */

typedef struct {
  Input_meta_t meta;       /* Input metadata */
  int nband;               /* Number of input image files (bands) */
  int nband_th;            /* Number of thermal input image files (0, 1, 2) */
  Img_coord_int_t size;    /* Input file size */
  Img_coord_int_t size_th; /* Input (thermal) file size */
  char *file_name[NBAND_REFL_MAX];  
                           /* Name of the input image files */
  char *file_name_th[2];   /* Name of the thermal input image files
                              [0] is TM band 6  OR
                              [0] is ETM+ band 6L, [1] is ETM+ band 6H */
  char *file_name_sun_zen; /* Name of the represetative per-pixel solar zenith
                              file */
  char *file_name_band_qa; /* Name of the L1 QA band*/
  bool open[NBAND_REFL_MAX]; 
                           /* Flag to indicate whether the specific input file 
			      is open for access; 'true' = open, 
			     'false' = not open */
  bool open_th;            /* thermal open flag */
  bool open_sun_zen;       /* solar zenith open flag */
  bool open_band_qa;       /* band QA open flag */
  FILE *fp_bin[NBAND_REFL_MAX];  /* File pointer for binary files */
  FILE *fp_bin_th[2];      /* File pointer for thermal binary file;
                              [0] is TM band 6  OR
                              [0] is ETM+ band 6L, [1] is ETM+ band 6H */
  FILE *fp_bin_sun_zen;    /* File pointer for the representative per-pixel
                              array solar zenith band */
  FILE *fp_bin_band_qa;    /* File pointer for the L1 QA band */
} Input_t;

/* Prototypes */

Input_t *OpenInput(Espa_internal_meta_t *metadata);
bool GetInputLine(Input_t *this, int iband, int iline, unsigned char *line);
bool GetInputLineTh(Input_t *this, int iline, unsigned char *line,
  unsigned char *line_b6h);
bool GetInputLineQA(Input_t *this, int iline, uint16_t*line);
bool GetInputLineSunZen(Input_t *this, int iline, int16 *line);
bool CloseInput(Input_t *this);
bool FreeInput(Input_t *this);
bool GetXMLInput(Input_t *this, Espa_internal_meta_t *metadata);

#endif
