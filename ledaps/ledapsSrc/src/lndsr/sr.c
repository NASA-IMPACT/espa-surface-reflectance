#include <string.h>
#include "sr.h"
#include "ar.h"
#include "const.h"
#include "sixs_runs.h"
#include "read_level1_qa.h"

/* !Revision:
 *
 * revision 1.2.1 3/22/2013  Gail Schmidt, USGS
 * - writing UL and LR corners to the output metadata to be able to detect
 *   ascending scenes or scenes where the image is flipped North to South
 * revision 8/11/2015  Gail Schmidt, USGS
 * - input saturated pixels are flagged as such and output as saturated
 */

void SrInterpAtmCoef(Lut_t *lut, double grid_line, double grid_sample,
                     atmos_t *atmos_coef, atmos_t *interpol_atmos_coef);

bool Sr
(
    Lut_t *lut,           /* I: lookup table information */
    int nsamp,            /* I: number of samples to be processed */
    int il,               /* I: current line being processed */
    atmos_t *atmos_coef,  /* I: atmospheric coefficients */
    atmos_t *interpol_atmos_coef, /* I: storage space for interpolated
                                        atmospheric coefficients */
    uint16_t **line_in,   /* I: array of input lines, one for each band */
    uint16_t *qa_line,    /* I: array of QA data for the current line */
    uint16_t **line_out,  /* O: array of output lines, one for each band */
    Sr_stats_t *sr_stats  /* O: statistics for this line */
)
{
    int is;                   /* current sample in the line */
    int ib;                   /* current band for this pixel */
    double grid_line, grid_sample; /* interpolation grid line and sample */
    float rho;                /* surface reflectance value */

    /* loop through the samples in this line */
    grid_line = (double)il/lut->ar_region_size.l - 0.5;
    double sample_step = 1./lut->ar_region_size.s;
    for (is = 0, grid_sample = -0.5; is < nsamp;
         is++, grid_sample += sample_step) {

        /* Interpolate the atmospheric coefficients for the current line/sample
           location */
        /* NAZMI 6/2/04 : correct even cloudy pixels */
        SrInterpAtmCoef(lut, grid_line, grid_sample, atmos_coef,
                        interpol_atmos_coef);

        /* Loop through each band, correcting the pixel.  Fill and saturated
           pixels are skipped and flagged. */
        for (ib = 0; ib < lut->nband; ib++) {
            if (level1_qa_is_fill(qa_line[is])) {
                /* fill pixel */
                line_out[ib][is] = lut->output_fill;
                sr_stats->nfill[ib]++;
                continue;
            }
            if (line_in[ib][is] == lut->in_satu) {
                /* saturated pixel */
                line_out[ib][is] = lut->output_satu;
                sr_stats->nsatu[ib]++;
                continue;
            }

            rho = compute_rho(line_in[ib][is] * lut->scale_factor 
                              + lut->add_offset,
                              *interpol_atmos_coef->tgOG[ib],
                              *interpol_atmos_coef->tgH2O[ib],
                              *interpol_atmos_coef->td_ra[ib],
                              *interpol_atmos_coef->tu_ra[ib],
                              *interpol_atmos_coef->rho_ra[ib],
                              *interpol_atmos_coef->S_ra[ib]);

            /* Scale the reflectance value for output and store it as 
               an uint16 */
            float temp;
            temp = (rho - lut->add_offset) * lut->mult_factor;
    
            /* Verify the reflectance value is within the valid range */
                if (temp < lut->min_valid_sr) {
                sr_stats->nout_range[ib]++;
                line_out[ib][is] = lut->min_valid_sr;
            }
                else if (temp > lut->max_valid_sr) {
                sr_stats->nout_range[ib]++;
                line_out[ib][is] = lut->max_valid_sr;
            }
                else {
                    line_out[ib][is] = (unsigned short) temp;
            }
    
            /* Keep track of the min/max value for the stats */
            if (sr_stats->first[ib]) {
                sr_stats->sr_min[ib] = sr_stats->sr_max[ib] = line_out[ib][is];
                sr_stats->first[ib] = false;
            }
            else if (line_out[ib][is] < sr_stats->sr_min[ib])
                sr_stats->sr_min[ib] = line_out[ib][is];
            else if (line_out[ib][is] > sr_stats->sr_max[ib])
                sr_stats->sr_max[ib] = line_out[ib][is];
        }  /* end for ib */
    }  /* end for is */

    return true;
}


void SrInterpAtmCoef
(
    Lut_t *lut,                    /* I: lookup table info */
    double grid_line,              /* I: grid line location */
    double grid_sample,            /* I: grid sample location */
    atmos_t *atmos_coef,           /* I: actual atmospheric coefficients */
    atmos_t *interpol_atmos_coef   /* O: interpolated atmospheric coefficients
                                         for the current line/samp */
)
/* 
  Point order:

    0 ---- 1    +--> sample
    |      |    |
    |      |    v
    2 ---- 3   line

NOTE: A handful of the coefficients are never used in the interpolated form.
  Therefore, to save computation time, they will be left out of the
  interpolation.
 */
{
    Img_coord_int_t p0, p1, p2;      /* 3 of the corner points for the
                                        aerosol interpolation */
    int i, n, ipt[4], ib;
    int lindex0, lindex1;
    double dl, ds, w[4];
    double sum_w;              /* sum of the weights */

    /* Set the four corner point indices. */
    p0.l = (int)grid_line;
    p2.l = p0.l + 1;
    lindex0 = p0.l*lut->ar_size.s;
    lindex1 = lindex0 + lut->ar_size.s;
    if (p2.l >= lut->ar_size.l) {
        p2.l = lut->ar_size.l - 1;
        lindex1 = p2.l*lut->ar_size.s;
        if (p0.l > 0)
        {
            p0.l--;
            lindex0 -= lut->ar_size.s;
        }
    }

    p0.s = (int)grid_sample;
    p1.s = p0.s + 1;
    if (p1.s >= lut->ar_size.s) {
        p1.s = lut->ar_size.s - 1;
        if (p0.s > 0)
            p0.s--;
    }

    /* Initialize the variables to 0 */
    n = 0;
    sum_w = 0.0;
    for (ib = 0; ib < 6; ib++) {
        *interpol_atmos_coef->tgOG[ib] = 0;
        *interpol_atmos_coef->tgH2O[ib] = 0;
        *interpol_atmos_coef->td_ra[ib] = 0;
        *interpol_atmos_coef->tu_ra[ib] = 0;
        *interpol_atmos_coef->rho_ra[ib] = 0;
        *interpol_atmos_coef->S_ra[ib] = 0;
    }

    /* Compute the fractional grid cell offset of the sample point from
       the upper-left corner of the cell. */
    dl = grid_line - p0.l;
    ds = grid_sample - p0.s;
    w[0] = (1 - dl)*(1 - ds);
    w[1] = 1 - dl - w[0];
    w[2] = 1 - ds - w[0];
    w[3] = dl - w[2];

    /* Loop through the four points to be used in the interpolation */
    ipt[0] = lindex0 + p0.s;
    ipt[1] = lindex0 + p1.s;
    ipt[2] = lindex1 + p0.s;
    ipt[3] = lindex1 + p1.s;
    for (i = 0; i < 4; i++) {
        if (!atmos_coef->computed[ipt[i]])
            continue;

        /* Increment the count of valid points, and add the current weight
           to the sum of weights. */
        n++;
        sum_w += w[i];

        /* Loop through each band, and add in the coefficient * weight. */
        for (ib = 0; ib < 6; ib++) {
            *interpol_atmos_coef->tgOG[ib] +=
                                          atmos_coef->tgOG[ib][ipt[i]]*w[i];
            *interpol_atmos_coef->tgH2O[ib] +=
                                          atmos_coef->tgH2O[ib][ipt[i]]*w[i];
            *interpol_atmos_coef->td_ra[ib] +=
                                          atmos_coef->td_ra[ib][ipt[i]]*w[i];
            *interpol_atmos_coef->tu_ra[ib] +=
                                          atmos_coef->tu_ra[ib][ipt[i]]*w[i];
            *interpol_atmos_coef->rho_ra[ib] +=
                                          atmos_coef->rho_ra[ib][ipt[i]]*w[i];
            *interpol_atmos_coef->S_ra[ib] +=
                                          atmos_coef->S_ra[ib][ipt[i]]*w[i];
        }
    }

    /* Divide by the sum of the weights if it's not zero or one. */
    if (n > 0 && n < 4) {
        double inv_sum_w = 1/sum_w;
        for (ib = 0; ib < 6; ib++) {
            *interpol_atmos_coef->tgOG[ib] *= inv_sum_w;
            *interpol_atmos_coef->tgH2O[ib] *= inv_sum_w;
            *interpol_atmos_coef->td_ra[ib] *= inv_sum_w;
            *interpol_atmos_coef->tu_ra[ib] *= inv_sum_w;
            *interpol_atmos_coef->rho_ra[ib] *= inv_sum_w;
            *interpol_atmos_coef->S_ra[ib] *= inv_sum_w;
        }
    }
}
