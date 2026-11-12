
#include <stdlib.h>
#include "lndsr.h"
#include "lut.h"
#include "input.h"
#include "mystring.h"
#include "error.h"

#define OUTPUT_FILL (0)
#define OUTPUT_FILL_OPACITY (-9999)
#define OUTPUT_SATU (USHRT_MAX)
#define MIN_VALID_SR (-0.2)         /* unscaled */
#define MAX_VALID_SR (1.6)          /* unscaled */
#define MIN_VALID_OPACITY (-2000)
#define MAX_VALID_OPACITY (16000)
#define AEROSOL_FILL (-9999)
#define AEROSOL_REGION_NLINE (40)
#define AEROSOL_REGION_NSAMP (AEROSOL_REGION_NLINE)
#define LONG_NAME_PREFIX ("band %d reflectance")
#define UNITS            ("reflectance")
#define ATMOS_OPACITY_SCALE_FACTOR (0.001)

Lut_t *GetLut(int nband, Input_meta_t *meta, Input_meta_t *b6_meta, 
              Img_coord_int_t *input_size) {
  Lut_t *this;

  /* Create the lookup table data structure */

  this = (Lut_t *)malloc(sizeof(Lut_t));
  if (this == NULL) 
    RETURN_ERROR("allocating Input data structure", "OpenInput", NULL);

  /* Populate the data structure */
  this->nband = nband;
  this->output_fill = OUTPUT_FILL;
  this->output_fill_opacity = OUTPUT_FILL_OPACITY;
  this->in_satu = meta->saturate_value;
  this->output_satu = OUTPUT_SATU;
  this->aerosol_fill = AEROSOL_FILL;
  this->ar_region_size.l = AEROSOL_REGION_NLINE;
  this->ar_region_size.s = AEROSOL_REGION_NSAMP;
  this->ar_size.l = ((input_size->l - 1) / this->ar_region_size.l) + 1;
  this->ar_size.s = ((input_size->s - 1) / this->ar_region_size.s) + 1;

  this->scale_factor = get_scale_refl();        /* scale factor            */
  this->mult_factor = 1 / this->scale_factor;   /* multiplication factor   */
  this->add_offset = get_offset_refl();         /* add offset              */
  this->min_valid_sr = (MIN_VALID_SR - this->add_offset) * this->mult_factor;
  if (this->min_valid_sr < 0)
    this->min_valid_sr = 0;
  this->max_valid_sr = (MAX_VALID_SR - this->add_offset) * this->mult_factor; 
  if (this->max_valid_sr > USHRT_MAX)
    this->max_valid_sr = USHRT_MAX;

  this->min_valid_opacity = MIN_VALID_OPACITY;
  this->max_valid_opacity = MAX_VALID_OPACITY;
  this->atmos_opacity_scale_factor= ATMOS_OPACITY_SCALE_FACTOR;
  this->b6_scale_factor = b6_meta->scale_factor;
  this->b6_add_offset =  b6_meta->add_offset;

  this->long_name_prefix = DupString(LONG_NAME_PREFIX);
  if (this->long_name_prefix == NULL) {
    free(this);
    RETURN_ERROR("duplicating long name prefix", "GetLut", NULL);
  }

  this->units = DupString(UNITS);
  if (this->units == NULL) {
    free(this);
    RETURN_ERROR("duplicating ref units", "GetLut", NULL);
  }

  if (!InputMetaCopy(meta, nband, &this->meta)) {
    free(this);
    RETURN_ERROR("copying input metadata", "GetLut", NULL);
  }

  return this;
}

void FreeLut(Lut_t *this) {
    free(this);
    this = NULL;
}
