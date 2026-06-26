#include "ar.h"
#include "const.h"
#include "error.h"
#include "sixs_runs.h"
#include "clouds.h"
#include "read_level1_qa.h"

/* #define VRA_THRESHOLD 0.1 */
#define VRA_THRESHOLD 0.08

int allocate_mem_atmos_coeff(int nbpts,atmos_t *atmos_coef);
int free_mem_atmos_coeff(atmos_t *atmos_coef);
void SrInterpAtmCoef(Lut_t *lut, double grid_line, double grid_sample,
                     atmos_t *atmos_coef, atmos_t *interpol_atmos_coef);


bool cloud_detection_pass1
(
    Lut_t *lut,              /* I: lookup table informat */
    int nsamp,               /* I: number of samples to be processed */
    int il,                  /* I: current line being processed */
    uint16_t **line_in,      /* I: array of input lines, one for each band */
    uint16_t *qa_line,       /* I: QA data (bqa_pixel) for the current line */
    uint16_t *qa2_line,      /* I: QA data (bqa_radsat) for the current line */
    uint16_t *b6_line,       /* I: array of thermal data for the current line */
    float *atemp_line,       /* I: auxiliary temperature for the line */
    atmos_t *atmos_coef,     /* I: atmospheric coefficients */
    atmos_t *interpol_atmos_coef, /* I: storage space for interpolated
                                        atmospheric coefficients */
    cld_diags_t *cld_diags   /* I/O: cloud diagnostics (stats are updated) */
)
{
    int is;                   /* current sample in the line */
    float tmpflt;             /* temporary floating point value */
    float rho1, rho3, rho4, rho5, rho7, t6;  /* reflectance and temp values */
    int C1, C2=0, C3, C4=0, C5=0, water;  /* cloud and water
                                                            indicators */
    int cld_row, cld_col;     /* cloud line, sample location */
    float vra, ndvi;          /* NDVI value */
    double grid_line, grid_sample; /* interpolation grid line and sample */

    /* Start the location for the current line and current cloud row. */
    cld_row = il / cld_diags->cellheight;

    /* Calculate a threshold that was a hardcoded magic number before.
       This was set to 5000 with the original scaling factor. */
    int thresh_tm = (0.5 - lut->add_offset) * lut->mult_factor;

    /* Loop through the samples in this line */
    grid_line = (double)il/lut->ar_region_size.l - 0.5;
    double sample_step = 1./lut->ar_region_size.s;
    int next_cld_step = cld_diags->cellwidth;
    for (is = 0, grid_sample = -0.5, cld_col = 0; is < nsamp;
         is++, grid_sample += sample_step) {
        if (is == next_cld_step) {
            cld_col++;
            next_cld_step += cld_diags->cellwidth;
        }

        if (!level1_qa_is_fill(qa_line[is]) &&
            (!level1_qa_is_saturated(qa2_line[is],3) ||
             (lut->meta.inst == INST_TM && line_in[AR_B3][is] < thresh_tm)))
        { /* not fill and no saturation in band 3 */
            /* Interpolate the atmospheric coefficients for the current
               pixel */
            SrInterpAtmCoef(lut, grid_line, grid_sample, atmos_coef,
                            interpol_atmos_coef);

            /* Compute the reflectance for each band using the interpolated
               atmospheric coefficients */
            uint16 *lin[5] = {line_in[AR_B1], line_in[AR_B3], line_in[AR_B4],
                              line_in[AR_B5], line_in[AR_B7]};
            float *rho[5] = {&rho1, &rho3, &rho4, &rho5, &rho7};
            float tgog[5] = {*interpol_atmos_coef->tgOG[0],
                             *interpol_atmos_coef->tgOG[2],
                             *interpol_atmos_coef->tgOG[3],
                             *interpol_atmos_coef->tgOG[4],
                             *interpol_atmos_coef->tgOG[5]};
            float tgh2o[5] = {*interpol_atmos_coef->tgH2O[0],
                              *interpol_atmos_coef->tgH2O[2],
                              *interpol_atmos_coef->tgH2O[3],
                              *interpol_atmos_coef->tgH2O[4],
                              *interpol_atmos_coef->tgH2O[5]};
            float rho_ra[5] = {*interpol_atmos_coef->rho_ra[0],
                               *interpol_atmos_coef->rho_ra[2],
                               *interpol_atmos_coef->rho_ra[3],
                               *interpol_atmos_coef->rho_ra[4],
                               *interpol_atmos_coef->rho_ra[5]};
            float td_ra[5] = {*interpol_atmos_coef->td_ra[0],
                              *interpol_atmos_coef->td_ra[2],
                              *interpol_atmos_coef->td_ra[3],
                              *interpol_atmos_coef->td_ra[4],
                              *interpol_atmos_coef->td_ra[5]};
            float tu_ra[5] = {*interpol_atmos_coef->tu_ra[0],
                              *interpol_atmos_coef->tu_ra[2],
                              *interpol_atmos_coef->tu_ra[3],
                              *interpol_atmos_coef->tu_ra[4],
                              *interpol_atmos_coef->tu_ra[5]};
            float s_ra[5] = {*interpol_atmos_coef->S_ra[0],
                             *interpol_atmos_coef->S_ra[2],
                             *interpol_atmos_coef->S_ra[3],
                             *interpol_atmos_coef->S_ra[4],
                             *interpol_atmos_coef->S_ra[5]};
            int i;
            for (i=0; i<5; i++)
                *rho[i] = compute_rho(lin[i][is] * lut->scale_factor 
                                      + lut->add_offset, tgog[i], tgh2o[i],
                                      td_ra[i], tu_ra[i], rho_ra[i], s_ra[i]);

            /* Get the temperature */
            t6 = b6_line[is] * lut->b6_scale_factor + lut->b6_add_offset;

            /* Compute cloud coefficients */
            vra = rho1 - rho3 * 0.5;
                    
            C1 = (int)(vra > VRA_THRESHOLD);
            C2 = (t6 < atemp_line[is] - 7);
   
            tmpflt = rho4 / rho3;
            C3 = (tmpflt >= 0.9 && tmpflt <= 1.3);

            C4 = (rho7 > 0.03);
            C5 = (rho3 > 0.6 || rho4 > 0.6);
                    
            /**
               Water test :
               ndvi < 0 => water
               ((0<ndvi<0.1) or (b4<5%)) and b5 < 0.01 => turbid water
            **/
            if (rho4 + rho3 != 0)
                ndvi = (rho4 - rho3)/(rho4 + rho3);
            else
                ndvi = 0.01;
            water = ndvi < 0 || (((ndvi > 0 && ndvi < 0.1) || rho4 < 0.05) &&
                                 rho5 < 0.02);
            if (!water &&
                (t6 > atemp_line[is] - 20. && !C5) &&
                !((C1 || C3) && C2 && C4))
            {
                /* not water and clear */
                cld_diags->avg_t6_clear[cld_row][cld_col] += t6;
                cld_diags->std_t6_clear[cld_row][cld_col] += t6*t6;
                cld_diags->avg_b7_clear[cld_row][cld_col] += rho7;
                cld_diags->std_b7_clear[cld_row][cld_col] += rho7*rho7;
                cld_diags->nb_t6_clear[cld_row][cld_col]++;
            }
        }  /* not fill and no saturation in band 3 */
    }  /* end for is */

    return true;
}


bool cloud_detection_pass2
(
    Lut_t *lut,              /* I: lookup table informat */
    int nsamp,               /* I: number of samples to be processed */
    int il,                  /* I: current line being processed */
    uint16_t **line_in,      /* I: array of input lines, one for each band */
    uint16_t *qa_line,       /* I: QA data (bqa_pixel) for the current line */
    uint16_t *qa2_line,      /* I: QA data (bqa_radsat) for the current line */
    uint16_t *b6_line,       /* I: array of thermal data for the current line */
    atmos_t *atmos_coef,     /* I: atmospheric coefficients */
    atmos_t *interpol_atmos_coef, /* I: storage space for interpolated
                                        atmospheric coefficients */
    cld_diags_t *cld_diags,  /* I: cloud diagnostics */
    char *ddv_line           /* O: dark dense vegetation line */
            /**
            use ddv_line to store internal cloud screening info
            bit 2 = adjacent cloud 1=yes 0=no
            bit 3 = fill value 1=fill 0=valid
            bit 4 = land/water mask 1=land 0=water
            bit 5 = cloud 0=clear 1=cloudy
            bit 6 = cloud shadow 
            bit 7 = snow
            **/
)
{
    int is;                   /* current sample in the line */
    int il_ar, is_ar;         /* line/sample in the aerosol region */
    bool thermal_band;        /* is thermal data available */
    float rho1, rho2, rho3, rho4, rho5, rho7;  /* reflectance & temp values */
    float t6 = 0.0;           /* temperature values */
    int C1, C2=0, C4=0, C5=0, water;  /* cloud and water
                                                            indicators */
    int cld_row, cld_col;     /* cloud line, sample location */
    float vra, ndvi, ndsi, temp_snow_thshld; /* NDVI, NDSI, snow threshold */
    float temp_b6_clear,temp_thshld1,temp_thshld2,atemp_ancillary;
    float tmpflt_arr[2];      /* temporary floats */
    double grid_line, grid_sample; /* interpolation grid line and sample */

    /* Initialize the thermal band information and snow threshold */
    thermal_band = true;
    if (b6_line == NULL) 
        thermal_band = false;

    /* This is an unscaled value in Kelvin.  Note: This threshold is outside
       the current maximum valid value for the thermal band. */
    temp_snow_thshld = 380.;  /* now flag snow and possibly salt pan */

    /* Calculate a threshold that was a hardcoded magic number before.
       This was set to 2000 with the original scaling factor. */
    int band5_thresh = (0.2 - lut->add_offset) * lut->mult_factor;

    /* Calculate a threshold that was a hardcoded magic number before.
       This was set to 5000 with the original scaling factor. */
    int thresh_tm = (0.5 - lut->add_offset) * lut->mult_factor;

    /* Start the location for the current line and current cloud row. */
    cld_row = il / cld_diags->cellheight;

    /* Initialize the aerosol line location value */
    il_ar = il / lut->ar_region_size.l;
    if (il_ar >= lut->ar_size.l)
        il_ar = lut->ar_size.l - 1;

    /* Loop through the samples in this line */
    int ar_sample = 0;
    grid_line = (double)il/lut->ar_region_size.l - 0.5;
    double sample_step = 1./lut->ar_region_size.s;
    int next_cld_step = cld_diags->cellwidth;
    for (is = 0, cld_col = 0, is_ar = 0, grid_sample = -0.5; is < nsamp;
         is++, ar_sample++, grid_sample += sample_step) {
        if (is == next_cld_step)
        {
            cld_col++;
            next_cld_step += cld_diags->cellwidth;
        }
        if (ar_sample == lut->ar_region_size.s && is_ar < lut->ar_size.s - 1)
        {
            is_ar++;
            ar_sample = 0;
        }

        /* If fill, continue with the next sample. */
        if (level1_qa_is_fill(qa_line[is])) {
            ddv_line[is] = AR_FILL;
            continue;
        }

        ddv_line[is] &= AR_ADJ_CLOUD_OR_SHADOW; /* reset all bits except cloud
                                                   shadow and adjacent cloud */ 

        if (level1_qa_is_saturated(qa2_line[is],3) ||
            ((lut->meta.inst == INST_TM) &&
            (line_in[AR_B3][is] >= thresh_tm))) {
            if (thermal_band) {
                t6 = b6_line[is]*lut->b6_scale_factor + lut->b6_add_offset;

                /* Interpolate the cloud diagnostics for current pixel */
                interpol_clddiags_1pixel (cld_diags, il, is, tmpflt_arr);
                temp_b6_clear = tmpflt_arr[0];
                atemp_ancillary = tmpflt_arr[1];
                if (temp_b6_clear < 0.) {
                    temp_thshld1 = atemp_ancillary - 20.;
                    temp_thshld2 = atemp_ancillary - 20.;
                }
                else if (cld_diags->std_t6_clear[cld_row][cld_col] > 0.) {
                    temp_thshld1 = temp_b6_clear -
                        (cld_diags->std_t6_clear[cld_row][cld_col] + 4.);
                    temp_thshld2 = temp_b6_clear -
                        cld_diags->std_t6_clear[cld_row][cld_col];
                }
                else {
                    temp_thshld1 = temp_b6_clear - 4.;
                    temp_thshld2 = temp_b6_clear - 2.;
                }

                if ((level1_qa_is_saturated(qa2_line[is],5) ||
                     ((lut->meta.inst == INST_TM) &&
                      (line_in[AR_B5][is] >= thresh_tm))) &&
                      (t6 < temp_thshld1)) {
                    /* saturated band 5 and t6 < threshold => cloudy */
                    ddv_line[is] &= AR_RESET_SHADOW; /* reset shadow bit */
                    ddv_line[is] &= AR_RESET_ADJ_CLOUD; /* reset adjcloud bit */
                    ddv_line[is] |= AR_CLOUD; /* set cloud bit */
                }
                else if ((line_in[AR_B5][is] < band5_thresh) &&
                         (t6 < temp_snow_thshld)) { /* snow */
                    ddv_line[is] |= AR_SNOW;
                }
                else { /* assume cloudy */
                    ddv_line[is] &= AR_RESET_SHADOW; /* reset shadow bit */
                    ddv_line[is] &= AR_RESET_ADJ_CLOUD; /* reset adjcloud bit */
                    ddv_line[is] |= AR_CLOUD; /* set cloud bit */
                }
            }  /* end if thermal */
        }  /* end if band 3 saturated .... */
        else {
            /* Interpolate the atmospheric conditions for current pixel */
            SrInterpAtmCoef(lut, grid_line, grid_sample, atmos_coef,
                            interpol_atmos_coef);

            /* Compute the reflectance for each band using the interpolated
               atmospheric coefficients */
            float *rho[6] = {&rho1, &rho2, &rho3, &rho4, &rho5, &rho7};
            float **tgog = interpol_atmos_coef->tgOG;
            float **tgh2o = interpol_atmos_coef->tgH2O;
            float **rho_ra = interpol_atmos_coef->rho_ra;
            float **td_ra = interpol_atmos_coef->td_ra;
            float **tu_ra = interpol_atmos_coef->tu_ra;
            float **s_ra = interpol_atmos_coef->S_ra;
            int i;
            for (i=0; i<6; i++)
            {
                *rho[i] = compute_rho(line_in[i][is] * lut->scale_factor
                                      + lut->add_offset, *tgog[i], *tgh2o[i],
                                      *td_ra[i], *tu_ra[i], *rho_ra[i],
                                      *s_ra[i]);
            }

            /* Get the temperature */
            if (thermal_band)
                t6 = b6_line[is] * lut->b6_scale_factor + lut->b6_add_offset;

            /* Interpolate the cloud diagnostics for the current pixel */
            interpol_clddiags_1pixel (cld_diags, il, is, tmpflt_arr);
            temp_b6_clear = tmpflt_arr[0];
            atemp_ancillary = tmpflt_arr[1];

            if (temp_b6_clear < 0.) {
                temp_thshld1 = atemp_ancillary - 20.;
                temp_thshld2 = atemp_ancillary - 20.;
            }
            else if (cld_diags->std_t6_clear[cld_row][cld_col] > 0.) {
                temp_thshld1 = temp_b6_clear -
                    (cld_diags->std_t6_clear[cld_row][cld_col] + 4.);
                temp_thshld2 = temp_b6_clear -
                    cld_diags->std_t6_clear[cld_row][cld_col];
            }
            else {
                temp_thshld1 = temp_b6_clear - 4.;
                temp_thshld2 = temp_b6_clear - 2.;
            }

            if (thermal_band) {
                /* Compute cloud coefficients */
                vra = rho1 - rho3 * 0.5;

                C1 = (int)(vra > VRA_THRESHOLD);
                C2 = (t6 < temp_thshld1);
   
                C4 = (rho7 > 0.03);
                C5 = (t6 < temp_thshld2) && C1;
            }

            /**
               Water test :
               ndvi < 0 => water
               ((0<ndvi<0.1) or (b4<5%)) and b5 < 0.01 => turbid water
            **/
            if (rho4 + rho3 != 0)
                ndvi = (rho4 - rho3) / (rho4 + rho3);
            else
                ndvi = 0.01;
            water = ndvi < 0 || (((ndvi > 0 && ndvi < 0.1) || rho4 < 0.05) &&
                                 rho5 < 0.02);



            if (thermal_band) {
                if (!water) { /* if not water */
                    ddv_line[is] |= AR_WATER;
                    if ((C2 || C5) && C4) { /* cloudy */
                        ddv_line[is] &= AR_RESET_SHADOW; /* reset shadow bit */
                        ddv_line[is] &= AR_RESET_ADJ_CLOUD; /* reset adjcloud */
                        ddv_line[is] |= AR_CLOUD; /* set cloud bit */
                    }
                    else { /* clear */
                        ddv_line[is] &= AR_RESET_CLOUD;
                        ndsi = (rho2 - rho5) / (rho2 + rho5);
                        if ((ndsi > 0.3) && (t6 < temp_snow_thshld) &&
                            (rho4 > 0.2))
                            ddv_line[is] |= AR_SNOW;
                    }
                }
                else 
                    ddv_line[is] &= AR_RESET_WATER;
            }
            else { /* no thermal band - cannot run cloud mask */
                ddv_line[is] &= AR_RESET_CLOUD; /* assume clear */
                if (!water) { /* if not water */
                    ddv_line[is] |= AR_WATER;
                }
                else {
                    ddv_line[is] &= AR_RESET_WATER;
                }
            }
        }  /* end else saturated band 3 */
    }  /* end for is */

    return true;
}


bool dilate_cloud_mask
(
    Lut_t *lut,
    int nsamp,
    char ***cloud_buf,
    int dilate_dist
)
{
    int il,is,il_adj,is_adj,buf_ind;
    int k;

    for (il = 0; il < lut->ar_region_size.l; il++) {
        int k_start = il - dilate_dist;
        for (is = 0; is < nsamp; is++) {
            /* If not cloudy, continue with next sample. */
            if (!(cloud_buf[1][il][is] & AR_CLOUD))
                continue;

            /* Cloudy, so dilate. */
            if (k_start <= 0)
            {
                buf_ind = 0;
                il_adj = k_start + lut->ar_region_size.l;
            }
            else
            {
                buf_ind = 1;
                il_adj = k_start;
            }
            for (k = k_start; k < il + dilate_dist; k++, il_adj++) {
                if (k == 0) {
                    buf_ind = 1;
                    il_adj = 0;
                }
                else if (k == lut->ar_region_size.l) {
                    buf_ind = 2;
                    il_adj = 0;
                }

                /* If il_adj is out of range, continue with the next value. */
                if (il_adj < 0 || il_adj >= lut->ar_region_size.l)
                    continue;

                int is_start, is_end;
                if (is > dilate_dist)
                    is_start = is - dilate_dist;
                else
                    is_start = 0;
                if (is < nsamp - dilate_dist)
                    is_end = is + dilate_dist;
                else
                    is_end = nsamp;
                for (is_adj = is_start; is_adj < is_end; is_adj++) {
                    if (!(cloud_buf[buf_ind][il_adj][is_adj] & AR_CLOUD)) {
                        /* not cloudy */
                        /* reset shadow bit */
                        cloud_buf[buf_ind][il_adj][is_adj] &= AR_RESET_SHADOW;
                        /* set adjacent cloud bit */
                        cloud_buf[buf_ind][il_adj][is_adj] |= AR_ADJ_CLOUD;
                    }
                }  /* for is_adj */
            }  /* for k */
        }  /* for is */
    }  /* for il */
    return true;
}


void cast_cloud_shadow
(
    Lut_t *lut,
    int nsamp,
    int il_start,
    uint16_t ***line_in,
    uint16_t **b6_line,
    cld_diags_t *cld_diags,
    char ***cloud_buf,
    Ar_gridcell_t *ar_gridcell,
    float pixel_size,
    float adjust_north
)
{
    int il,is,il_ar,is_ar,shd_buf_ind;
    float t6,temp_b6_clear,atemp_ancillary,tmpflt_arr[2];
    float conv_factor,cld_height,ts,fs,dx,dy;
    int shd_x,shd_y;
    float shadow_factor = 1000/pixel_size;

/***
    Cloud Shadow
***/
    il_ar = il_start / lut->ar_region_size.l;
    if (il_ar >= lut->ar_size.l)
        il_ar = lut->ar_size.l - 1;
    int index = il_ar*lut->ar_size.s;
    for (il = 0; il <lut->ar_region_size.l; il++) {
        int ar_sample = 0;
        for (is = 0, is_ar = 0; is < nsamp; is++, ar_sample++) {
            if (ar_sample == lut->ar_region_size.s &&
                is_ar < lut->ar_size.s - 1)
            {
                is_ar++;
                ar_sample = 0;
            }

            /* Get the thermal info */
            t6 = b6_line[il][is] * lut->b6_scale_factor + lut->b6_add_offset;

            /* Interpolate the cloud diagnostics for this pixel */
            interpol_clddiags_1pixel (cld_diags, il+il_start, is, tmpflt_arr);
            temp_b6_clear = tmpflt_arr[0];
            atemp_ancillary = tmpflt_arr[1];

            /* If not cloudy, continue with next sample. */
            if (!(cloud_buf[1][il][is] & AR_CLOUD))
                continue;

            conv_factor = 6.;
            while (conv_factor <= 6.) {
                /* Determine the cloud height */
                if (temp_b6_clear > 0)
                    cld_height = (temp_b6_clear - t6) / conv_factor;
                else
                    cld_height = (atemp_ancillary - t6) / conv_factor;

                /* If the cloud height is less than or equal to 0, there is
                   no shadow.  So continue with the next value.  */
                conv_factor++;
                if (cld_height <= 0)
                    continue;

                ts = ar_gridcell->sun_zen[index+is_ar]/DEG;
                fs = (ar_gridcell->rel_az[index+is_ar] - adjust_north)/DEG;

                dy = cos(fs) * tan(ts) * cld_height;
                dx = sin(fs) * tan(ts) * cld_height;
                shd_x = is - dx*shadow_factor;
                shd_y = il + dy*shadow_factor;

                if ((shd_x >= 0) && (shd_x < nsamp)) {
                    shd_buf_ind = 1;
                    if (shd_y < 0) {
                        shd_buf_ind--;
                        shd_y += lut->ar_region_size.l;
                    }
                    if (shd_y >= lut->ar_region_size.l) {
                        shd_buf_ind++;
                        shd_y -= lut->ar_region_size.l;
                    }
                    /* Mask as cloud shadow */
                    if (shd_y >= 0 && shd_y < lut->ar_region_size.l) {
                        /* if not cloud, adjacent cloud or cloud
                           shadow */
                        if (!((cloud_buf[shd_buf_ind][shd_y][shd_x] &
                               AR_CLOUD) ||
                              (cloud_buf[shd_buf_ind][shd_y][shd_x] &
                               AR_ADJ_CLOUD) ||
                              (cloud_buf[shd_buf_ind][shd_y][shd_x] &
                               AR_CLOUD_SHADOW)))
                            /* set cloud shadow bit */
                            cloud_buf[shd_buf_ind][shd_y][shd_x] |=
                                AR_CLOUD_SHADOW;
                    }
                }
            } /* while conv_fact <= 6. */
        } /* sample loop */
    } /* line loop */
}

bool dilate_shadow_mask
(
    Lut_t *lut,          /* I: lookup table */
    int nsamp,           /* I: number of samples in the current line */
    char *fill_mask,     /* I: storage for fill mask */
    char ***cloud_buf,   /* I/O: cloud buffer */
    int dilate_dist      /* I: size of dilation window */
)
{
    int il,is,il_adj,is_adj,buf_ind;
    int k;

    /* Initialize the fill mask. */
    memset(fill_mask, 0, lut->ar_region_size.l*nsamp);

    int index = 0;
    for (il = 0; il < lut->ar_region_size.l; il++) {
        int k_start;
        if (il > dilate_dist)
            k_start = il - dilate_dist;
        else
            k_start = 0;

        for (is = 0; is < nsamp; is++, index++) {
            /* If not cloud shadow, continue with the next sample.
               Otherwise, dilate. */
            if (!(cloud_buf[0][il][is] & AR_CLOUD_SHADOW) || fill_mask[index])
                continue;

            for (k = k_start, il_adj = k, buf_ind = 0; k <= il + dilate_dist;
                 k++, il_adj++) {
                if (k == lut->ar_region_size.l) {
                    buf_ind = 1;
                    il_adj = 0;
                }

                int adj_start, adj_end;
                if (is > dilate_dist)
                    adj_start = is - dilate_dist;
                else
                    adj_start = 0;
                if (is < nsamp - dilate_dist - 1)
                    adj_end = is + dilate_dist;
                else
                    adj_end = nsamp - 1;
                int index_adj = il_adj*nsamp + adj_start;
                for (is_adj = adj_start; is_adj <= adj_end;
                     is_adj++, index_adj++) {
                    /* if not cloud, adjacent cloud or cloud shadow */
                    if (!(cloud_buf[buf_ind][il_adj][is_adj] &
                          AR_CLOUD_OR_ADJ_CLOUD_OR_SHADOW)) {
                        /* set adjacent cloud shadow bit */
                        cloud_buf[buf_ind][il_adj][is_adj] |= AR_CLOUD_SHADOW;
                        fill_mask[index_adj] = 1;
                    }
                }
            }
        }  /* for is */
    }  /* for il */

    return true;
}

int allocate_cld_diags
(
    cld_diags_t *cld_diags,   /* O: cloud diagnostics */
    int cell_height,          /* I: cell height of cloud diagnostics */
    int cell_width,           /* I: cell width of cloud diagnostics */
    int scene_height,         /* I: number of lines in the scene */
    int scene_width           /* I: number of samples in the scene */
)
{
    int i;
    
    /* Set the size of the diagnostic cells and the number of rows/cols in
       the diagnostics */
    cld_diags->cellheight = cell_height;
    cld_diags->cellwidth = cell_width;
    cld_diags->nbrows = (scene_height - 1) / cell_height + 1;
    cld_diags->nbcols = (scene_width - 1) / cell_width + 1;

    if ((cld_diags->avg_t6_clear = malloc (cld_diags->nbrows*sizeof(double *)))
        ==NULL)
        return -1;
    for (i = 0; i < cld_diags->nbrows; i++)
        if ((cld_diags->avg_t6_clear[i] = calloc(cld_diags->nbcols,
            sizeof(double)))==NULL)
            return -1;

    if ((cld_diags->std_t6_clear = malloc (cld_diags->nbrows*sizeof(double *)))
        ==NULL)
        return -1;
    for (i = 0; i < cld_diags->nbrows; i++)
        if ((cld_diags->std_t6_clear[i] = calloc(cld_diags->nbcols,
            sizeof(double)))==NULL)
            return -1;

    if ((cld_diags->avg_b7_clear = malloc (cld_diags->nbrows*sizeof(double *)))
        ==NULL)
        return -1;
    for (i = 0; i < cld_diags->nbrows; i++) 
        if ((cld_diags->avg_b7_clear[i] = calloc(cld_diags->nbcols,
            sizeof(double)))==NULL)
            return -1;

    if ((cld_diags->std_b7_clear = malloc (cld_diags->nbrows*sizeof(double *)))
        ==NULL)
        return -1;
    for (i = 0; i < cld_diags->nbrows; i++) 
        if ((cld_diags->std_b7_clear[i] = calloc(cld_diags->nbcols,
            sizeof(double)))==NULL)
            return -1;

    if ((cld_diags->airtemp_2m = malloc (cld_diags->nbrows*sizeof(float *)))
        ==NULL)
        return -1;
    for (i = 0; i < cld_diags->nbrows; i++) 
        if ((cld_diags->airtemp_2m[i] = calloc(cld_diags->nbcols,
            sizeof(float)))==NULL)
            return -1;

    if ((cld_diags->nb_t6_clear = malloc (cld_diags->nbrows*sizeof(int *)))
        ==NULL)
        return -1;
    for (i = 0; i < cld_diags->nbrows; i++) 
        if ((cld_diags->nb_t6_clear[i] = calloc(cld_diags->nbcols,
            sizeof(int)))==NULL)
            return -1;

    return 0;
}

void free_cld_diags(cld_diags_t *cld_diags) {
    
    int i;
    
    for (i = 0; i < cld_diags->nbrows; i++) {
        free(cld_diags->avg_t6_clear[i]);
        free(cld_diags->std_t6_clear[i]);
        free(cld_diags->avg_b7_clear[i]);
        free(cld_diags->std_b7_clear[i]);
        free(cld_diags->airtemp_2m[i]);
        free(cld_diags->nb_t6_clear[i]);
    }

    free(cld_diags->avg_t6_clear);
    free(cld_diags->std_t6_clear);
    free(cld_diags->avg_b7_clear);
    free(cld_diags->std_b7_clear);
    free(cld_diags->airtemp_2m);
    free(cld_diags->nb_t6_clear);
}

void fill_cld_diags(cld_diags_t *cld_diags) {
/*
!Description: fill in missing values in the T6Clear grid based
on existing values (spatial interpolation). Missing values have been previously
set to -9999. A filling can be distincted from a original value by looking at
the standard deviation of the optical depth which is set to -9999 for a filling.
!END****************************************************************************
*/
    int i,j,k,l,pass,count;
    float lastt6=0.0, dist, sumt6=0.0, sumdist=0.0, lastb7=0.0, sumb7=0.0;
    char missing_flag[300][300];
    int min_nb_values,n,max_distance;
   
    count=0;
    for (i=0;i<cld_diags->nbrows;i++) {
        for (j=0;j<cld_diags->nbcols;j++) {
            missing_flag[i][j]=1;
            if (cld_diags->avg_t6_clear[i][j]!=CLOUD_FILL)  {
                count++;
                lastt6=cld_diags->avg_t6_clear[i][j];
                lastb7=cld_diags->avg_b7_clear[i][j];
                missing_flag[i][j]=0;
            }
        }
    }

    if (count==0)
        return;

    else if (count==1) {
        for (i=0;i<cld_diags->nbrows;i++)
            for (j=0;j<cld_diags->nbcols;j++) {
                cld_diags->avg_t6_clear[i][j]=lastt6;
                cld_diags->avg_b7_clear[i][j]=lastb7;
            }

        return;
    }
    
    /* Loop through the lines and samples in the cloud diagnostics */
    for (i=0;i<cld_diags->nbrows;i++) {
        for (j=0;j<cld_diags->nbcols;j++) {
            /**
            Look for at least 3 neighboring valid values within 4 GPs
            **/
            min_nb_values=3;
            max_distance=4;    
            pass=0;
            while ((cld_diags->avg_t6_clear[i][j] == CLOUD_FILL) &&
                   (pass<max_distance)) {
                pass++;
                sumdist=0.;
                sumt6=0.;
                sumb7=0.;
                n=0;
                for (k=i-pass;k<=(i+pass);k++) {
                    if ((k>=0)&&(k<cld_diags->nbrows)) {
                        for (l=j-pass;l<=(j+pass);l++) {
                            if ((l>=0)&&(l<cld_diags->nbcols)) {
                                if (!missing_flag[k][l]) {
                                    dist=sqrt((k-i)*(k-i)+(l-j)*(l-j));
                                    sumdist += dist;
                                    sumt6 += (dist *
                                        cld_diags->avg_t6_clear[k][l]);
                                    sumb7 += (dist *
                                        cld_diags->avg_b7_clear[k][l]);
                                    n++;
                                }
                            }
                        }
                    }
                }

                if ((n>=min_nb_values)&&(sumdist!=0.)) {
                    cld_diags->avg_t6_clear[i][j] = sumt6 / sumdist;
                    cld_diags->avg_b7_clear[i][j] = sumb7 / sumdist;
                }
            }  /* end while */

            /**
            Look for at least 2 neighboring valid values within 6 GPs
            **/
            min_nb_values=2;
            max_distance=6;    
            pass=0;
            while ((cld_diags->avg_t6_clear[i][j] == CLOUD_FILL) &&
                   (pass<max_distance)) {
                pass++;
                sumdist=0.;
                sumt6=0.;
                sumb7=0.;
                n=0;
                for (k=i-pass;k<=(i+pass);k++) {
                    if ((k>=0)&&(k<cld_diags->nbrows)) {
                        for (l=j-pass;l<=(j+pass);l++) {
                            if ((l>=0)&&(l<cld_diags->nbcols)) {
                                if (!missing_flag[k][l]) {
                                    dist=sqrt((k-i)*(k-i)+(l-j)*(l-j));
                                    sumdist += dist;
                                    sumt6 += (dist *
                                        cld_diags->avg_t6_clear[k][l]);
                                    sumb7 += (dist *
                                        cld_diags->avg_b7_clear[k][l]);
                                    n++;
                                }
                            }
                        }
                    }
                }

                if ((n>=min_nb_values)&&(sumdist!=0.)) {
                    cld_diags->avg_t6_clear[i][j]=sumt6/sumdist;
                    cld_diags->avg_b7_clear[i][j]=sumb7/sumdist;
                }
            }  /* end while */

            /**
            Look for at least 1 neighboring valid values within 10 GPs
            **/
            min_nb_values=1;
            max_distance=10;    
            pass=0;
            while ((cld_diags->avg_t6_clear[i][j] == CLOUD_FILL) &&
                   (pass<max_distance)) {
                pass++;
                sumdist=0.;
                sumt6=0.;
                sumb7=0.;
                n=0;
                for (k=i-pass;k<=(i+pass);k++) {
                    if ((k>=0)&&(k<cld_diags->nbrows)) {
                        for (l=j-pass;l<=(j+pass);l++) {
                            if ((l>=0)&&(l<cld_diags->nbcols)) {
                                if (!missing_flag[k][l]) {
                                    dist=sqrt((k-i)*(k-i)+(l-j)*(l-j));
                                    sumdist += dist;
                                    sumt6 += (dist *
                                        cld_diags->avg_t6_clear[k][l]);
                                    sumb7 += (dist *
                                        cld_diags->avg_b7_clear[k][l]);
                                    n++;
                                }
                            }
                        }
                    }
                }

                if ((n>=min_nb_values)&&(sumdist!=0.)) {
                    cld_diags->avg_t6_clear[i][j]=sumt6/sumdist;
                    cld_diags->avg_b7_clear[i][j]=sumb7/sumdist;
                }
            }  /* end while */
        }  /* for j */
    }  /* for i */
}

void interpol_clddiags_1pixel
(
    cld_diags_t *cld_diags,  /* I: cloud diagnostics */
    int img_line,            /* I: current line in image */
    int img_sample,          /* I: current sample in image */
    float *inter_value       /* O: interpolated cloud diagnostic value */
)
/* 
  Point order:

    0 ---- 1    +--> sample
    |      |    |
    |      |    v
    2 ---- 3   line

    inter_value[0] => t6_clear
    inter_value[1] => airtemp_2m

    Updated by Gail Schmidt, USGS EROS, on 10/20/2014
    We want the airtemp_2m to be calculated regardless of whether the band6
        temp is available (i.e. there were clear pixels to compute the average
        thermal temp).  Many of the users of these interpolated values use the
        airtemp_2m as the default if the band6 clear temps are not valid.
 */

{
    typedef struct {
      int l;                /* line number */
      int s;                /* sample number */
    } Img_coord_int_t;
    
    Img_coord_int_t p[FOUR_PTS];
    int i, n, n_anc;
    float dl, ds, w[FOUR_PTS];
    float grid_line, grid_sample;
    float sum_w, sum_anc_w;

    int cell_half_height, cell_half_width;

    inter_value[0] = inter_value[1] = 0;

    cell_half_height = (cld_diags->cellheight + 1) >> 1;  /* divide by 2 */
    cell_half_width = (cld_diags->cellwidth + 1) >> 1;  /* divide by 2 */

    /* Set the four corner point indices. */
    grid_line = (float)(img_line  - cell_half_height)/cld_diags->cellheight;
    grid_sample = (float)(img_sample - cell_half_width)/cld_diags->cellwidth;
    p[0].l = (int)grid_line;
    if (p[0].l < 0)
        p[0].l = 0;
    p[2].l = p[0].l + 1;
    if (p[2].l >= cld_diags->nbrows) {
        p[2].l = cld_diags->nbrows - 1;
            if (p[0].l > 0)
                p[0].l--;
    }
    p[1].l = p[0].l;
    p[3].l = p[2].l;

    p[0].s = (int)grid_sample;
    if (p[0].s < 0)
        p[0].s = 0;
    p[1].s = p[0].s + 1;
    if (p[1].s >= cld_diags->nbcols) {
        p[1].s = cld_diags->nbcols - 1;
        if (p[0].s > 0)
            p[0].s--;
    }
    p[2].s = p[0].s;
    p[3].s = p[1].s;

    /* Initialize the counting and sum variables */
    n = 0;
    n_anc = 0;
    sum_w = 0.0;
    sum_anc_w = 0.0;

    /* Compute the fractional grid cell offset of the sample point from
       the upper-left corner of the cell and the weights for each corner
       sammple. */
    dl = grid_line - p[0].l;
    ds = grid_sample - p[0].s;
    w[0] = (1 - dl)*(1 - ds);
    w[1] = 1 - dl - w[0];
    w[2] = 1 - ds - w[0];
    w[3] = ds - w[1];

    /* Apply the bilinear interpolation for the sample point. */
    for (i = 0; i < 4; i++)
    {
        if (cld_diags->avg_t6_clear[p[i].l][p[i].s] != CLOUD_FILL)
        {
            n++;
            sum_w += w[i];
            inter_value[0] += cld_diags->avg_t6_clear[p[i].l][p[i].s]*w[i];
        }

        if (cld_diags->airtemp_2m[p[i].l][p[i].s] != CLOUD_FILL)
        {
            n_anc++;
            sum_anc_w += w[i];
            inter_value[1] += cld_diags->airtemp_2m[p[i].l][p[i].s]*w[i];
        }
    }

    /* Divide by the sum of the weights if it's not zero or one. */
    if (sum_w > 0 && n < 4)
        inter_value[0] /= sum_w;
    else if (sum_w == 0)
        inter_value[0] = CLOUD_FILL;

    if (sum_anc_w > 0 && n_anc < 4)
        inter_value[1] /= sum_anc_w;
    else if (sum_anc_w == 0)
        inter_value[1] = CLOUD_FILL;

}
