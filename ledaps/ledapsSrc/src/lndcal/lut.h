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

 ! Design Notes:
   1. Structure is declared for the 'input' data type.
  
!END****************************************************************************
*/

#ifndef LUT_H
#define LUT_H

#include <stdbool.h>
#include "lndcal.h"
#include "input.h"
#include "param.h"

#define VALID_MIN_REF        (-0.01) /* Unscaled */
#define VALID_MAX_REF        (1.6)   /* Unscaled */
#define VALID_MIN_TH         (150)   /* Unscaled */
#define VALID_MAX_TH         (350)   /* Unscaled */

/* Structure for the 'lut' data type */

typedef struct {
  char *file_name;             /* Lookup table file name                    */
  int out_fill;                /* Output fill value                         */
  int out_satu;                /* Output saturation value (Feng, 3/23/09)   */
  Input_meta_t meta;           /* Input metadata                            */
  float cos_sun_zen;           /* Cosine of the solar zenith angle          */
  float esun[NBAND_REFL_MAX];  /* Mean solar exoatmospheric irradiances     */ 
  float dsun2;                 /* Earth sun distance squared                */
  float K1[2];                 /* K1 constant
                                  [0] is TM band 6  OR
                                  [0] is ETM+ band 6L, [1] is ETM+ band 6H  */
  float K2[2];                 /* K2 constant
                                  [0] is TM band 6  OR
                                  [0] is ETM+ band 6L, [1] is ETM+ band 6H  */
  char* long_name_prefix_ref;  /* ref long name prefix (append band num)    */
  char* units_ref;             /* ref units                                 */
  int valid_range_ref[2];      /* ref valid range                           */
  double scale_factor_ref;     /* ref scale factor                          */
  double mult_factor_ref;      /* ref multiplication factor                 */
  double add_offset_ref;       /* ref add offset                            */
  char* long_name_prefix_th;   /* thermal long name prefix (append band num)*/ 
  char* units_th;              /* thermal units                             */
  int valid_range_th[2];       /* thermal valid range                       */
  double scale_factor_th;      /* thermal scale factor                      */
  double mult_factor_th;       /* thermal multiplication factor             */
  double add_offset_th;        /* thermal add offset                        */
  double refl_conv[NBAND_REFL_MAX];
} Lut_t;

/* Prototypes */

Lut_t *GetLut(Param_t *param, int nband, Input_t *input);
bool FreeLut(Lut_t *this);

#endif
