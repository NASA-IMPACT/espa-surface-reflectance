#ifndef AR_H
#define AR_H

#include <stdbool.h>
#include "lndsr.h"
#include "lut.h"
#include "sixs_runs.h"

/* QA values for ddv_line */
#define AR_DDV 0x01
#define AR_ADJ_CLOUD 0x04
#define AR_FILL 0x08
#define AR_WATER 0x10
#define AR_CLOUD 0x20
#define AR_CLOUD_OR_ADJ_CLOUD 0x24
#define AR_CLOUD_SHADOW 0x40
#define AR_ADJ_CLOUD_OR_SHADOW 0x44
#define AR_CLOUD_OR_ADJ_CLOUD_OR_SHADOW 0x64
#define AR_SNOW 0x80
#define AR_RESET_ADJ_CLOUD 0xfb
#define AR_RESET_SHADOW 0xbf
#define AR_RESET_CLOUD 0xdf
#define AR_RESET_WATER 0xef
#define AR_CLEAR

/* Band number index for access to line_in */
#define AR_B1 0
#define AR_B2 1
#define AR_B3 2
#define AR_B4 3
#define AR_B5 4
#define AR_B7 5

typedef struct {
  bool first;
  int ar_min;
  int ar_max;
  long nfill;
} Ar_stats_t;

typedef struct
{
    unsigned short b[3];   /* collect_band */
    unsigned short b7;     /* collect_band7 */
} collect_bands_t;

bool Ar(int il_ar, Lut_t *lut, Img_coord_int_t *size_in, uint16_t ***line_in,
        char **ddv_line, atmos_t *atmos_coef_ar, collect_bands_t *cbands,
        int **line_ar, Ar_stats_t *ar_stats, Ar_gridcell_t *ar_gridcell,
        sixs_tables_t *sixs_tables);

int ArInterp(Lut_t *lut, Img_coord_int_t *loc, int ***line_ar, int *inter_aot);
int Fill_Ar_Gaps(Lut_t *lut, int ***line_ar, int ib);

/* Compute the surface reflectance for a given point and atmospheric
   coefficients. */
static inline float compute_rho(float unscaled_val, float tgOG, float tgH2O,
                                float td_ra, float tu_ra, float rho_ra,
                                float S_ra)
{
    float rho = unscaled_val / tgOG - rho_ra;
    rho /= tgH2O*td_ra*tu_ra + S_ra*rho;

    return rho;
}

#endif
