#ifndef CAL_H
#define CAL_H

#include <stdbool.h>
#include "lndcal.h"
#include "lut.h"
#include "input.h"

/* NOTE: Be careful of modifying these high end saturation values. They
   currently jive with the values used to produce the Level-1 RADSAT band,
   but if they are changed they will be out of sync with the Level-1 RADSAT
   band used in lndsr. */
static const int SATU_VAL[7]={255,255,255,255,255,255,255};
static const int SATU_VAL6= 255;
static const int SATU_LOW_VAL6=1;

bool Cal(Param_t *param, Lut_t *lut, int iband, Input_t *input,
         unsigned char *line_in, int16 *line_in_sun_zen, uint16_t *line_out,
         uint16_t *line_qa, int iy);

bool Cal_TM_thermal(Lut_t *lut, Input_t *input, unsigned char *line_in,
          uint16_t *line_out, uint16_t *line_qa, int iy);

bool Cal_ETM_thermal(Lut_t *lut, Input_t *input, unsigned char *line_b6l,
          unsigned char *line_b6h, uint16_t *line_out, uint16_t *line_qa,
          int iy);

#endif
