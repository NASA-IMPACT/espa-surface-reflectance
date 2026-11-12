#ifndef CLOUDS_H
#define CLOUDS_H

#include <stdbool.h>
#include "lndsr.h"
#include "lut.h"
#include "sixs_runs.h"

#define CLOUD_FILL -9999
#define FOUR_PTS 4
#define CLDDIAGS_CELLHEIGHT_1KM 40
#define CLDDIAGS_CELLWIDTH_1KM 40
#define CLDDIAGS_CELLHEIGHT_5KM 160
#define CLDDIAGS_CELLWIDTH_5KM 160
#define CLDDIAGS_CELLHEIGHT_10KM 330
#define CLDDIAGS_CELLWIDTH_10KM 330

typedef struct cld_diags_t {
	int nbrows,nbcols,cellheight,cellwidth;
	double **avg_t6_clear,**std_t6_clear,**avg_b7_clear,**std_b7_clear;
	float **airtemp_2m;
	int **nb_t6_clear;
}cld_diags_t;

int allocate_cld_diags(struct cld_diags_t *cld_diags,int cell_height,
                       int cell_width, int scene_height, int scene_width);
void free_cld_diags(struct cld_diags_t *cld_diags);
void fill_cld_diags(cld_diags_t *cld_diags);

void interpol_clddiags_1pixel(cld_diags_t *cld_diags, int img_line,
                              int img_sample,float *inter_value);
bool cloud_detection_pass1(Lut_t *lut, int nsamp, int il, uint16_t **line_in,
                           uint16_t *qa_line, uint16_t *qa2_line,
                           uint16_t *b6_line,float *atemp_line,
                           atmos_t *atmos_coef, atmos_t *interpol_atmos_coef,
                           cld_diags_t *cld_diags);
bool cloud_detection_pass2(Lut_t *lut, int nsamp, int il, uint16_t **line_in,
                           uint16_t *qa_line, uint16_t *qa2_line,
                           uint16_t *b6_line, atmos_t *atmos_coef,
                           atmos_t *interpol_atmos_coef, cld_diags_t *cld_diags,
                           char *ddv_line);
void cast_cloud_shadow(Lut_t *lut, int nsamp, int il_start, uint16_t ***line_in,
                       uint16_t **b6_line, cld_diags_t *cld_diags,
                       char ***cloud_buf, Ar_gridcell_t *ar_gridcell,
                       float pixel_size, float adjust_north);
bool dilate_cloud_mask(Lut_t *lut, int nsamp, char ***cloud_buf, int dilate_dist);
bool dilate_shadow_mask(Lut_t *lut, int nsamp, char *fill_mask,
                        char ***cloud_buf, int dilate_dist);

#endif
