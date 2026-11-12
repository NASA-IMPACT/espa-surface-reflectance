#include "cal.h"
#include "const.h"
#include "error.h"
#include "read_level1_qa.h"
#include "local_defines.h"

#define nint(A)(A<0?(int)(A-0.5):(int)(A+0.5))

/* Functions */
/* !Revision:
 *
 * NOTES:
 * 1. TOA radiance and reflectance equations for Landsat 7 are available in
 *    http://landsathandbook.gsfc.nasa.gov/data_prod/prog_sect11_3.html
 * 2. The TOA reflectance gain/bias values from the MTL file (stored in the
 *    XML file after converting from LPGS to ESPA) do not account for the
 *    solar angle.  Thus the gain and bias need to be applied and then we
 *    still need to account for the solar angle.
 */

bool Cal(Param_t *param, Lut_t *lut, int iband, Input_t *input,
         unsigned char *line_in, int16 *line_in_sun_zen, uint16_t *line_out,
         uint16_t *line_qa, int iy) {
  int is,val;
  float rad_gain = 0, rad_bias = 0;   /* TOA radiance gain/bias */
  float refl_gain = 0.0,
        refl_bias = 0.0;              /* TOA reflectance gain/bias */
  float rad;                          /* TOA radiance value */
  float ref_conv = 0.0;               /* TOA reflectance conversion value */
  float ref;                          /* TOA reflectance value */
  float fval;                         /* temporary float value */
  float sun_zen;                      /* solar zenith angle for the current
                                         pixel (radians) */
  float temp;

  int nsamp= input->size.s;

  /* Get the TOA reflectance gain/bias if they are available, otherwise use
     the TOA reflectance equation from the Landsat handbook. */
  if (input->meta.use_toa_refl_consts) {
    refl_gain = lut->meta.refl_gain[iband];
    refl_bias = lut->meta.refl_bias[iband];

    if ( iy==0 ) {
      printf("*** band=%1d refl gain=%f refl bias=%f "
            "cos_sun_zen(scene center)=%f\n", lut->meta.iband[iband],
            refl_gain, refl_bias, lut->cos_sun_zen);
      fflush(stdout);
    }
  }
  else {
    /* Get the TOA radiance gain/bias */
      rad_gain = lut->meta.rad_gain[iband];
      rad_bias = lut->meta.rad_bias[iband];

    ref_conv = (PI * lut->dsun2) / (lut->esun[iband] * lut->cos_sun_zen);
  
    if ( iy==0 ) {
      printf("*** band=%1d rad gain=%f rad bias=%f dsun2=%f\n"
             "    ref_conv=%f=(PI*%f)/(%f*%f) ***\n", lut->meta.iband[iband],
             rad_gain, rad_bias, lut->dsun2, ref_conv, lut->dsun2,
             lut->esun[iband], lut->cos_sun_zen);
      fflush(stdout);
    }
  }

  /* Loop through the samples in the line */
  for (is = 0; is < nsamp; is++) {
    if (level1_qa_is_fill(line_qa[is])) {
      line_out[is] = lut->out_fill;
      continue;
    }

    /* flag saturated pixels, added by Feng (3/23/09) */
    val = line_in[is];
    if (val == SATU_VAL[iband]) {
      line_out[is] = lut->out_satu;
      continue;
    }

    fval= (float)val;

    /* If the TOA reflectance gain/bias values are available, then use them.
       Otherwise compute the TOA radiance then reflectance, per the Landsat
       handbook equations. */
    if (input->meta.use_toa_refl_consts) {
      /* use per-pixel angles - convert the degree values to radians and then
         unscale */
      sun_zen = (line_in_sun_zen[is]*lut->meta.szen_scale
                 + lut->meta.szen_offset)*RAD;
      ref = ((refl_gain * fval) + refl_bias) / cos (sun_zen);
    }
    else {
      rad = rad_gain*fval + rad_bias;
      ref = rad * ref_conv;
    }

    /* Apply scaling. Values are set up in lut.c */
    temp = ((ref - lut->add_offset_ref) * lut->mult_factor_ref + 0.5);

    /* Cap the output using the min/max values */
    if (temp < lut->valid_range_ref[0]) {
      line_out[is] = lut->valid_range_ref[0];
    }
    else if (temp > lut->valid_range_ref[1]) {
      line_out[is] = lut->valid_range_ref[1];
    }
    else
      line_out[is] = temp;

  }  /* end for is */

  return true;
}

/* Brightness Temp Calibration for TM.  Band 6 is the only band available for
   computing brightness temperature. */
bool Cal_TM_thermal(Lut_t *lut, Input_t *input, unsigned char *line_in,
          uint16_t *line_out, uint16_t *line_qa, int iy) {
  int is, val;
  float rad_gain, rad_bias, rad, temp, temp2;
  int nsamp= input->size_th.s;

  /* The gain and bias are only for band 6 in this case */
  rad_gain = lut->meta.rad_gain_th[0];
  rad_bias = lut->meta.rad_bias_th[0];
  
  if ( iy==0 ) {
    printf("*** band=%1d gain=%f bias=%f ***\n", 6, rad_gain, rad_bias);
  }

  /* Loop through the pixels in this line */
  for (is = 0; is < nsamp; is++) {
    if (level1_qa_is_fill(line_qa[is])) {
      line_out[is] = lut->out_fill;
      continue;
    }

    /* for saturated pixels */
    val = line_in[is];
    if (val >= SATU_VAL6) {
      line_out[is] = lut->out_satu;
      continue;
    }

    /* compute the TOA brightness temperature in Kelvin and apply scaling.
       Values are set up in lut.c
       as well. */
    rad = (rad_gain * (float)val) + rad_bias;
    temp = lut->K2[0] / log(1.0 + (lut->K1[0]/rad));
    temp2 = (temp - lut->add_offset_th) * lut->mult_factor_th + 0.5;

    /* Cap the output using the min/max values */
    if (temp2 < lut->valid_range_th[0]) {
      line_out[is] = lut->valid_range_th[0];
    }
    else if (temp2 > lut->valid_range_th[1]) {
      line_out[is] = lut->valid_range_th[1];
    }
    else
        line_out[is] = (uint16_t) temp2;

  }  /* end for is */

  return true;
}


/* Brightness Temp Calibration for ETM+.  Band 6H is the primary band used for
   brightness temp.  If Band 6H is saturated, then Band 6L is used.  If it is
   also saturated then the pixel is masked as saturated. */
bool Cal_ETM_thermal(Lut_t *lut, Input_t *input, unsigned char *line_b6l,
          unsigned char *line_b6h, uint16_t *line_out, uint16_t *line_qa,
          int iy) {
  int is;
  int val[2];  /* input thermal value for band 6L, 6H */
  int cband;   /* index for 6L or 6H to be used in the combined band */
  float rad_gain[2], rad_bias[2], rad, temp, temp2;

  int nsamp= input->size_th.s;

  /* Set up the calibration coefficients */
  rad_gain[BAND_6L] = lut->meta.rad_gain_th[BAND_6L];
  rad_bias[BAND_6L] = lut->meta.rad_bias_th[BAND_6L];
  rad_gain[BAND_6H] = lut->meta.rad_gain_th[BAND_6H];
  rad_bias[BAND_6H] = lut->meta.rad_bias_th[BAND_6H];

  if ( iy==0 ) {
    printf("*** band=6L gain=%f bias=%f ***\n", rad_gain[BAND_6L],
      rad_bias[BAND_6L]);
    printf("*** band=6H gain=%f bias=%f ***\n", rad_gain[BAND_6H],
      rad_bias[BAND_6H]);
  }

  /* Loop through the pixels in this line */
  for (is = 0; is < nsamp; is++) {

    if (level1_qa_is_fill(line_qa[is])) {
      line_out[is] = lut->out_fill;
      continue;
    }

    val[BAND_6L] = line_b6l[is];
    val[BAND_6H] = line_b6h[is];

    /* If band 6H is not saturated, it will be used for the brightness temp */
    if((val[BAND_6H] < SATU_VAL6) && (val[BAND_6H] > SATU_LOW_VAL6))
    {
        cband = BAND_6H;
    }
    /* Otherwise if band 6L is not saturated, it will be used for the
       brightness temp */
    else if((val[BAND_6L] < SATU_VAL6) && (val[BAND_6L] > SATU_LOW_VAL6))
    {
        cband = BAND_6L;
    }
    /* If both are saturated, then the output is saturated. */
    else
    {
        line_out[is] = lut->out_satu;
        continue;
    }

    /* Compute the TOA brightness temperature in Kelvin for the selected band
     * and apply scaling.  Valid ranges are set up in lut.c as well. */
    rad = (rad_gain[cband] * (float)val[cband]) + rad_bias[cband];
    temp = lut->K2[cband] / log(1.0 + (lut->K1[cband]/rad));
    temp2 = (temp - lut->add_offset_th) * lut->mult_factor_th + 0.5;

    /* Cap the output using the min/max values */
    if (temp2 < lut->valid_range_th[0]) {
      line_out[is] = lut->valid_range_th[0];
    }
    else if (temp2 > lut->valid_range_th[1]) {
      line_out[is] = lut->valid_range_th[1];
    }
    else
        line_out[is] = (uint16_t) temp2;

  }  /* end for is */

  return true;
}
