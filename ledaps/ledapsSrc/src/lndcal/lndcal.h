#ifndef LNDCAL_H
#define LNDCAL_H

#define bounded(A,B,C) (A>B?(A<C?A:C):B)
#define min(A,B) (A>B?B:A)
#define max(A,B) (A>B?A:B)

#include <stdint.h>
#include "mystring.h"
#include "espa_metadata.h"
#include "parse_metadata.h"
#include "write_metadata.h"
#include "envi_header.h"

#define NBAND_REFL_MAX (6)

typedef signed short int16;
typedef unsigned char uint8;

typedef enum {
  ALG_NASA,
  ALG_NASA_CPF,
  ALG_CCRS 
} Algorithm_t;

/* Satellite type definition */

typedef enum {
  SAT_NULL = -1,
  SAT_LANDSAT_1 = 0, 
  SAT_LANDSAT_2, 
  SAT_LANDSAT_3, 
  SAT_LANDSAT_4, 
  SAT_LANDSAT_5, 
  SAT_LANDSAT_7, 
  SAT_MAX
} Sat_t;

extern const Key_string_t Sat_string[SAT_MAX];

/* Instrument type definition */

typedef enum {
  INST_NULL = -1,
  INST_MSS = 0, 
  INST_TM,
  INST_ETM, 
  INST_MAX
} Inst_t;

extern const Key_string_t Inst_string[INST_MAX];

/* World Reference System (WRS) type definition */

typedef enum {
  WRS_NULL = -1,
  WRS_1 = 0, 
  WRS_2,
  WRS_MAX
} Wrs_t;

extern const Key_string_t Wrs_string[WRS_MAX];

/* Band gain settings (ETM+ only) */

typedef enum {
  GAIN_NULL = -1,
  GAIN_HIGH = 0, 
  GAIN_LOW, 
  GAIN_MAX
} Gain_t;

extern const Key_string_t Gain_string[GAIN_MAX];

#ifndef IMG_COORD_INT_TYPE_DEFINED

#define IMG_COORD_INT_TYPE_DEFINED

/* Integer image coordinates data structure */

typedef struct {
  int l;                /* line number */
  int s;                /* sample number */
} Img_coord_int_t;

#endif

typedef struct {
  float gains[7];
  bool valid_flag;
} Gains_t;

typedef struct {
  double min_lon;  /* Geodetic longitude coordinate (degrees) */ 
  double min_lat;  /* Geodetic latitude coordinate (degrees) */ 
  double max_lon;  /* Geodetic longitude coordinate (degrees) */ 
  double max_lat;  /* Geodetic latitude coordinate (degrees) */ 
  bool is_fill;    /* Flag to indicate whether the point is a fill value; */
} Geo_bounds_t;

#define SCALE_FACTOR_REF     (0.0000275)
#define SCALE_FACTOR_TH      (0.00341802)
#define ADD_OFFSET_REF       (-0.2)
#define ADD_OFFSET_TH        (149.0)

/* data retreval functions*/
double get_scale_refl();    /* scale for reflective bands */
double get_offset_refl();   /* add offset for reflective bands */
double get_scale_therm();   /* scale for thermal bands */
double get_offset_therm();  /* add offset for thermal bands */

#endif
