#ifndef _COMMON_H_
#define _COMMON_H_

#include <stdlib.h>
#include "hdf.h"
#include "mfhdf.h"
typedef char byte;

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

/* Surface reflectance version */
#define SR_VERSION "3.5.2.0 (Collection 2)"

/* Define the default aerosol and EPS value */
#define DEFAULT_AERO 0.05
#define DEFAULT_EPS 1.5

/* Define the size of the aerosol window that will be used when running the
   aerosol inversion.  The aerosols will be inverted for the center of the
   NxN window (with special handling for fill, clouds, water) and then used
   to fill the rest of the window.  Aerosols are fairly homogenious over a
   reasonable area.  Landsat windows represent the center of the window.
   The following are example NxN window setups:
   9x9: AERO_WINDOW 9 and HALF_AERO_WINDOW 4
   7x7: AERO_WINDOW 7 and HALF_AERO_WINDOW 3
   5x5: AERO_WINDOW 5 and HALF_AERO_WINDOW 2
   3x3: AERO_WINDOW 3 and HALF_AERO_WINDOW 1
   1x1: AERO_WINDOW 1 and HALF_AERO_WINDOW 0
*/
/*#define LAERO_WINDOW 3
#define LHALF_AERO_WINDOW 1
#define SAERO_WINDOW 6*/
#define LAERO_WINDOW 1
#define LHALF_AERO_WINDOW 0
#define SAERO_WINDOW 1

/* Define the size of the window used for fixing the invalid aerosols, using
   an average of the representative pixels in this window. Define the minimum
   number of clear/valid pixels needed in the window in order for the average
   to be valid. The MIN_CLEAR_PIX is ~16.7% of the number of representative
   pixels in the FIX_AERO_WINDOW.
   Landsat ==> 15 x 30m pixels = 450m window for average
               There are 25 (5x5) possible representative pixels in a window
               of 15x15 using 3x3 representative pixel windows. 4/24 = .167
   Sentinel ==> 45 x 10m pixels = 450m window for average
               There are 49 (7x7) possible representative pixels in a window
               of 45x45 using 6x6 representative pixel windows. 8/48 = .167
*/
#define LFIX_AERO_WINDOW 15
#define LHALF_FIX_AERO_WINDOW 7
#define LMIN_CLEAR_PIX 4
#define SFIX_AERO_WINDOW 45
#define SHALF_FIX_AERO_WINDOW 22
#define SMIN_CLEAR_PIX 8

/* How many lines of data should be processed at one time */
#define PROC_NLINES 10

/* For angle conversions -
   degrees to radians = PI/180
   radians to degrees = 180/PI */
#define PI 3.14159265
#define PIx2 (PI * 2.0)
#define DEG2RAD 0.017453293
#define RAD2DEG 57.29577951

/* For divisions - to minimize processing time */
#define ATMOS_PRES_0 1013.0 /* mean atmospheric pressure (mbar) at sea level */
#define ONE_DIV_ATMOS_PRES_0 0.000987166
#define ONE_DIV_8500 0.000117647

/* Identify if all the Sentinel-2 bands will be processed. Default is to
   leave off for ESPA processing, therefore bands 9 and 10 will not be
   processed. */
#define PROC_ALL_BANDS 1

/* Number of bands corrected to surface reflectance
   * Landsat (bands 1-7).  The atmospheric correction variables store
   information for 8 bands, so we will go with that for the array size.
   * Sentinel-2 (bands 1-13, possibly skipping band 9 and 10).  The
   atmospheric correction variables store information for 11/13 bands, so we
   will go with that for the array size. */
#define NSRL_BANDS 8
#ifdef PROC_ALL_BANDS
/* Process all bands if turned on */
#define NSRS_BANDS 13
#else
/* Skip bands 9 and 10 as default for ESPA */
#define NSRS_BANDS 11
#endif

/* Get the larger number of bands between the two instruments */
#define NSR_BANDS MAX(NSRL_BANDS, NSRS_BANDS)

/* Landsat Level-1 products have 8 reflectance bands (bands 1-7, and 9),
   2 thermal bands (band 10 and 11), 1 pan band (band 8), and 1 QA band
   (band 12)
   Sentinel-2 Level-1 products have 13 reflectance bands, but we might not
   process bands 9 and 10 */
#define NBANDL_REFL_MAX 8
#define NBANDL_THM_MAX 2
#define NBANDL_PAN_MAX 1
#define NBANDL_QA_MAX 1
#define NBANDL_TTL_MAX (NBANDL_REFL_MAX + NBANDL_THM_MAX + NBANDL_PAN_MAX + NBANDL_QA_MAX)

#ifdef PROC_ALL_BANDS
/* Process all bands if turned on */
#define NBANDS_REFL_MAX 13
#else
/* Skip bands 9 and 10 as default for ESPA */
#define NBANDS_REFL_MAX 11
#endif
#define NBANDS_TTL_MAX NBANDS_REFL_MAX

/* Get the larger number of bands between the two instruments */
#define NBAND_REFL_MAX MAX(NBANDL_REFL_MAX, NBANDS_REFL_MAX)

/* Landsat surface reflectance products have 8 reflectance bands, 2 thermal
   bands, 0 pan bands, 1 aerosol band, and 1 QA band.
   Sentinel-2 surface reflectance products have 13 reflectance bands,
   1 aerosol band, and 1 QA band. Bands 9 (water vapor) and 10 (cirrus) might
   not be processed. */
#define NBANDL_REFL_OUT 8
#define NBANDL_THM_OUT 2
#define NBANDL_PAN_OUT 0
#define NBANDL_AEROSOL_OUT 1
#define NBANDL_QA_OUT 1
#define NBANDL_TTL_OUT (NBANDL_REFL_OUT + NBANDL_THM_OUT + NBANDL_PAN_OUT + NBANDL_AEROSOL_OUT + NBANDL_QA_OUT)

#ifdef PROC_ALL_BANDS
/* Process all bands if turned on */
#define NBANDS_REFL_OUT 13
#else
/* Skip bands 9 and 10 as default for ESPA */
#define NBANDS_REFL_OUT 11
#endif
#define NBANDS_AEROSOL_OUT 1
#define NBANDS_QA_OUT 1
#define NBANDS_TTL_OUT (NBANDS_REFL_OUT + NBANDS_AEROSOL_OUT + NBANDS_QA_OUT)

/* CMG and DEM files are lat/long images where each pixel represents 0.05 deg x
   0.05 deg */
/* DEM information */
#define DEM_NBLAT 3600
#define DEM_NBLON 7200

/* Ratio file information */
#define RATIO_NBLAT 3600
#define RATIO_NBLON 7200

/* Ozone and water vapor information */
/* Expected size of the CMG dataset (VIIRS auxiliary) is
   3600 lines x 7200 samples */
#define CMG_NDIMS 2
#define CMG_NBLAT 3600
#define CMG_NBLON 7200

/* Lookup table index value */
#define NPRES_VALS 7
#define NAOT_VALS 22
#define NSOLAR_VALS 8000
#define NSUNANGLE_VALS 22
#define NVIEW_ZEN_VALS 20
#define NSOLAR_ZEN_VALS 22
#define NAOTxNSOLAR_VALS (NAOT_VALS * NSOLAR_VALS)
#define NAOTxNSUNANGLE_VALS (NAOT_VALS * NSUNANGLE_VALS)

/* Coefficients for determining atmospheric values */
#define NCOEF 4
#define NREFLL_BANDS 7
#ifdef PROC_ALL_BANDS
/* Process all bands if turned on */
#define NREFLS_BANDS 13
#else
/* Skip bands 9 and 10 as default for ESPA */
#define NREFLS_BANDS 11
#endif

/* Get the larger number of bands between the two instruments */
#define NREFL_BANDS MAX(NREFLL_BANDS, NREFLS_BANDS)

/* Define the full list of bands for Sentinel-2, even though we aren't
   processing all these bands. Some look-up tables and input arrays are
   defined for the full list of bands. */
typedef enum {SBAND1=0, SBAND2, SBAND3, SBAND4, SBAND5, SBAND6,
    SBAND7, SBAND8, SBAND8A, SBAND9, SBAND10, SBAND11, SBAND12,
    SENTINEL_TTL} Mysband_t;
extern char SENTINEL_FULL_BANDNAME[SENTINEL_TTL][3]; /* defined in lut_subr.c */

/* Define the input products to be processed */
typedef enum {DNL_BAND1=0, DNL_BAND2, DNL_BAND3, DNL_BAND4,
    DNL_BAND5, DNL_BAND6, DNL_BAND7, DNL_BAND8, DNL_BAND9,
    DNL_BAND10, DNL_BAND11, DNL_QA, DNL_TTL} Mydnl_band_t;
#ifdef PROC_ALL_BANDS
/* Process all bands if turned on */
typedef enum {DNS_BAND1=0, DNS_BAND2, DNS_BAND3, DNS_BAND4, DNS_BAND5,
    DNS_BAND6, DNS_BAND7, DNS_BAND8, DNS_BAND8A, DNS_BAND9, DNS_BAND10,
    DNS_BAND11, DNS_BAND12, DNS_QA, DNS_TTL} Mydns_band_t;
#else
/* Skip bands 9 and 10 as default for ESPA */
typedef enum {DNS_BAND1=0, DNS_BAND2, DNS_BAND3, DNS_BAND4,
    DNS_BAND5, DNS_BAND6, DNS_BAND7, DNS_BAND8, DNS_BAND8A,
    DNS_BAND11, DNS_BAND12, DNS_QA, DNS_TTL} Mydns_band_t;
#endif

/* Define the output products to be processed */
typedef enum {SRL_BAND1=0, SRL_BAND2, SRL_BAND3, SRL_BAND4, SRL_BAND5,
    SRL_BAND6, SRL_BAND7, SRL_BAND9, SRL_BAND10, SRL_BAND11, SRL_AEROSOL,
    SRL_AEROSOL_QA, SRL_TTL} Mysrl_band_t;
#ifdef PROC_ALL_BANDS
/* Process all bands if turned on */
typedef enum {SRS_BAND1=0, SRS_BAND2, SRS_BAND3, SRS_BAND4, SRS_BAND5,
    SRS_BAND6, SRS_BAND7, SRS_BAND8, SRS_BAND8A, SRS_BAND9, SRS_BAND10,
    SRS_BAND11, SRS_BAND12, SRS_AEROSOL, SRS_AEROSOL_QA, SRS_TTL} Mysrs_band_t;
#else
/* Skip bands 9 and 10 as default for ESPA */
typedef enum {SRS_BAND1=0, SRS_BAND2, SRS_BAND3, SRS_BAND4, SRS_BAND5,
    SRS_BAND6, SRS_BAND7, SRS_BAND8, SRS_BAND8A, SRS_BAND11,
    SRS_BAND12, SRS_AEROSOL, SRS_AEROSOL_QA, SRS_TTL} Mysrs_band_t;
#endif
extern char SENTINEL_BANDNAME[NREFLS_BANDS][3];  /* defined in output.c */

/* For the minimal arrays that use this value, we will just use the largest
   of the Landsat and Sentinel number of output bands */
#define NBAND_TTL_OUT MAX((int)SRL_TTL, (int)SRS_TTL)

/* High confidence Level-1 QA values */
#define L1QA_HIGH_CONF 3

/* Low and average aerosol reflectance thresholds */
#define LOW_AERO_THRESH 0.015
#define AVG_AERO_THRESH 0.03

/* Bit values of ipflag (interpolation flag) QA, which includes aerosol
   levels */
typedef enum {
  IPFLAG_FILL=0,        /* fill value */
  IPFLAG_CLEAR=1,       /* valid aerosol retrieval (Landsat: water and land) */
  IPFLAG_VALID=1,       /* valid aerosol retrieval (Sentinel: water and land) */
  IPFLAG_WATER=2,       /* valid water pixel */
  IPFLAG_FAILED=3,      /* failed aerosol retrieval pixel (Sentinel) */
  IPFLAG_FIXED=4,       /* invalid retrieval which was fixed with a local
                           average of valid aerosols (Landsat is internal &
                           Sentinel is external -- slightly different fix/avg
                           algorithms) */
  IPFLAG_FAILED_TMP=5,  /* temp flag for expanding possible failed pixels
                           (Sentinel) */
  IPFLAG_INTERP_WINDOW=5, /* aerosol was interpolated using the center
                             of the NxN window (Landsat) */
  AERO1_QA=6,    /* these two AERO bits mark the amount of aerosols and = 64 */
  AERO2_QA=7     /* reflect the level of atmospheric correction made    = 128 */
} Ipflag_t;

/* Satellite type definitions, mainly to allow future satellites to be
   supported if needed */
typedef enum {
  SAT_NULL = -1,
  SAT_LANDSAT_8 = 0, 
  SAT_LANDSAT_9, 
  SAT_SENTINEL_2,
} Sat_t;

/* Instrument type definition */
typedef enum {
  INST_NULL = -1,
  INST_OLI_TIRS = 0, 
  INST_OLI, 
  INST_MSI,
} Inst_t;

/* World Reference System (WRS) type definition */
typedef enum {
  WRS_NULL = -1,
  WRS_1 = 0, 
  WRS_2,
} Wrs_t;

typedef struct {
  int nlines;
  int nsamps;
  double pixsize[2];
} Img_coord_info_t;

/* Atmospheric Aux Data Source used to index the following arrays */
typedef enum {
   AUXNULL = -1,
   MODIS = 0,
   VIIRS,
   GRIB
} aux_src_t;

/** Defines for MODIS and VIIRS scale and defaults **/
/* Water vapor scale factor */
//#define WV_SCALE_FACTOR (int[]) {200, 100}  /** USGS OPS with MODIS bug **/
//#define WV_SCALE_FACTOR (int[]) {100, 100}  /** Fixed wv scale for MODIS **/
#define WV_SCALE_FACTOR (int[]) { 200, 200, 200}  /** USGS OPS with MODIS bug, and
                                                added bug for VIIRS **/

/* Water Vapor Default (g/cm^2) - used for scene center computations */
#define WV_DEFAULT (float[]) { 0.5, 2.5, 0.5}

/* Water Vapor DN Default - use if WV is fill */
#define WV_DEFAULT_DN (int[]) {0, 250, 100}

/* Water Vapor fill */
#define WV_FILL (int[]) {0, 0, 0}

/* Ozone scale factor */
#define OZ_SCALE_FACTOR (int[]) {400, 400, 400}

/* Ozone Default (cm atm) - used for scene center computations */
#define OZ_DEFAULT (float[]) {0.3, 0.3, 0.3}

/* Ozone DN Default - use if OZ is fill */
#define OZ_DEFAULT_DN (int[]) {120, 110, 125}

/* Ozone fill */
#define OZ_FILL (int[]) {0, 0, 0}
#endif
