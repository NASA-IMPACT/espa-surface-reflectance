#include "stdlib.h"
#include "lut.h"
#include "input.h"
#include "mystring.h"
#include "error.h"
#include "const.h"

#define OUTPUT_FILL (0)
#define OUTPUT_SATU (USHRT_MAX)
#define LONG_NAME_PREFIX_REF ("band %d reflectance")
#define UNITS_REF            ("reflectance")
#define LONG_NAME_PREFIX_TH  ("band %d temperature")
#define UNITS_TH             ("temperature (kelvin)")

/* Relative sun-earth distance variation (Ref. L7 handbook)
 http://lpwww.gsfc.nasa.gov/IAS/handbook/handbook_htmls/chaper11/chapter11.html
   (web site visited on August 22, 2002)
*/

#define NDSUN (25)
typedef struct {
  int doy;         /* Day of year */
  float dsun;      /* Relative sun-earth distance */
} Dsun_table_t;

const Dsun_table_t dsun_table[NDSUN] = {
  {  1, 0.9832}, { 15, 0.9836}, { 32, 0.9853}, { 46, 0.9878}, { 60, 0.9909},
  { 74, 0.9945}, { 91, 0.9993}, {106, 1.0033}, {121, 1.0076}, {135, 1.0109}, 
  {152, 1.0140}, {166, 1.0158}, {182, 1.0167}, {196, 1.0165}, {213, 1.0149},
  {227, 1.0128}, {242, 1.0092}, {258, 1.0057}, {274, 1.0011}, {288, 0.9972},
  {305, 0.9925}, {319, 0.9892}, {335, 0.9860}, {349, 0.9843}, {366, 0.9833}
};


/* Solar Spectral Irradiances 
   Units = watts / (meter squared * micrometers)
   Reference:
     ETM+ values are from the Landsat 7 User Handbook (see above) 
     TM and MSS values from EOSAT Landsat Technical Notes, August 1986 
     
     New version of L4 and L5 given Solar Irradinaces in  (mW/m^2-nm) to Jeff Masek from Brian Markham on 11/3/06 : 
>L7: 1969      1840      1551      1044    225.7      82.07     1368
>L5: 1957      1826      1554      1036    215.0      80.67  -
>L4: 1957      1825      1557      1033    214.9      80.72  -
*/

const float esun_etm[8] = {
  1969.000, 1840.000, 1551.000, 1044.000, 225.700, -1.0, 82.07, 1368.000
};

/*
const float esun_tm_4[7] = {
  1958.0, 1828.0, 1559.0, 1045.0, 219.1, -1.0, 74.57
};
*/
const float esun_tm_4[7] = {
  1957.0, 1825.0, 1557.0, 1033.0, 214.9, -1.0, 80.72
};

/*
const float esun_tm_5[7] = {
  1957.0, 1829.0, 1557.0, 1047.0, 219.3, -1.0, 74.57
};
*/

const float esun_tm_5[7] = {
  1957.0, 1826.0, 1554.0, 1036.0, 215.0, -1.0, 80.67
};

const float esun_mss_1[4] = {
  1852.0, 1584.0, 1276.0, 904.0
};

const float esun_mss_2[4] = {
  1856.0, 1559.0, 1269.0, 906.0
};

const float esun_mss_3[4] = {
  1860.0, 1571.0, 1289.0, 910.0
};

const float esun_mss_4[4] = {
  1851.0, 1593.0, 1260.0, 878.0
};

const float esun_mss_5[4] = {
  1849.0, 1595.0, 1253.0, 870.0
};

const float K1_tm_4=  671.62;
const float K2_tm_4= 1284.30;
const float K1_tm_5=  607.76;
const float K2_tm_5= 1260.56;
const float K1_etm =  666.09;
const float K2_etm = 1282.71;

Lut_t *GetLut(Param_t *param, int nband, Input_t *input) {
  Lut_t *this;
  int ib, iband;
  int jdoy, i;
  float dsun;
  char msgbuf[1024];
  Input_meta_t *input_meta= &(input->meta);

  /* Create the lookup table data structure */
  this = (Lut_t *)calloc(1, sizeof(Lut_t));
  if (this == NULL) 
    RETURN_ERROR("allocating Input data structure", "OpenInput", NULL);

  /* Make a copy of the input metadata for the look-up table metadata.
     Note that if pointers are added to the Input_meta_t structure, their
     values will need to be copied to the new struct here. */
  this->meta = *input_meta;
  this->out_fill = OUTPUT_FILL;
  this->out_satu = OUTPUT_SATU;

  /* Check the dates for MSS
     Notes: values for Landsats 4, 5 are only valid for data acquired after
     4/1/83 and 11/9/84, respectively.  Values for data *processed* during
     1972-1978 may be different for Landsat 1-3. */                            
  if (input_meta->inst == INST_MSS) {
    if ( input_meta->sat == SAT_LANDSAT_1 ||
         input_meta->sat == SAT_LANDSAT_2 ||
         input_meta->sat == SAT_LANDSAT_3 ) {
      if ( input_meta->prod_date.year > 1971 &&
           input_meta->prod_date.year < 1979 ) {
        sprintf(msgbuf,"Landsat_%1.1d production year=%4.4d valid only "
          "between 1972 and 1978", (int)(input_meta->sat+1),
          input_meta->prod_date.year );
        RETURN_ERROR(msgbuf, "GetLut", false);
      }
    }
    if ( input_meta->sat == SAT_LANDSAT_4 ) {
      if ( input_meta->acq_date.year < 1983 || 
         ( input_meta->acq_date.year == 1983 &&
           input_meta->acq_date.doy <= 91 ) ) {
        sprintf(msgbuf,"Landsat_%1.1d acq date (%2.2d/%2.2d/%4.4d) invalid "
          "before 4/1/83", (int)(input_meta->sat+1),
          input_meta->acq_date.month, input_meta->acq_date.day,
          input_meta->acq_date.year);
        RETURN_ERROR(msgbuf, "GetLut", false);
      }
    }
    if ( input_meta->sat == SAT_LANDSAT_5 ) {
      if ( input_meta->acq_date.year < 1984 || 
         ( input_meta->acq_date.year == 1984 &&
           input_meta->acq_date.doy <= 314 ) ) {
        sprintf(msgbuf,"Landsat_%1.1d acq date (%2.2d/%2.2d/%4.4d) invalid "
          "before 11/19/84", (int)(input_meta->sat+1),
          input_meta->acq_date.month, input_meta->acq_date.day,
          input_meta->acq_date.year);
        RETURN_ERROR(msgbuf, "GetLut", false);
      }
    }
  }  /* end if inst == MSS */
  
  /* Compute the coefficients for the reflectance. Use the K1,K2 thermal consts
     from the XML file if they are available. */
  for (ib = 0; ib < nband; ib++) {
    iband = input_meta->iband[ib] - 1;
    switch (input_meta->sat) {
      case SAT_LANDSAT_1:
        this->esun[ib] = esun_mss_1[iband];
        break;

      case SAT_LANDSAT_2:
        this->esun[ib] = esun_mss_2[iband];
        break;

      case SAT_LANDSAT_3:
        this->esun[ib] = esun_mss_3[iband];
        break;

      case SAT_LANDSAT_4:
        if (input_meta->inst == INST_TM) {
          this->esun[ib] = esun_tm_4[iband];
          if ( ib==0 ){
            if (input_meta->use_toa_refl_consts) {
              this->K1[0] = input_meta->k1_const[0];
              this->K2[0] = input_meta->k2_const[0];
            }
            else {
              this->K1[0] = K1_tm_4;
              this->K2[0] = K2_tm_4;
            }
          }
        }
        else
          this->esun[ib] = esun_mss_4[iband];
        break;

      case SAT_LANDSAT_5:
        if (input_meta->inst == INST_TM) {
          this->esun[ib] = esun_tm_5[iband];
          if ( ib==0 ){
            if (input_meta->use_toa_refl_consts) {
              this->K1[0] = input_meta->k1_const[0];
              this->K2[0] = input_meta->k2_const[0];
            }
            else {
              this->K1[0] = K1_tm_5;
              this->K2[0] = K2_tm_5;
            }
          }
        }
        else
          this->esun[ib] = esun_mss_5[iband];
        break;

      case SAT_LANDSAT_7:
        this->esun[ib] = esun_etm[iband];
        if ( ib==0 ){
          if (input_meta->use_toa_refl_consts) {
            this->K1[0] = input_meta->k1_const[0];
            this->K1[1] = input_meta->k1_const[1];
            this->K2[0] = input_meta->k2_const[0];
            this->K2[1] = input_meta->k2_const[1];
          }
          else {
            this->K1[0] = K1_etm;
            this->K1[1] = K1_etm;
            this->K2[0] = K2_etm;
            this->K2[1] = K2_etm;
          }
        }
        break;

      case SAT_NULL: break;
      case SAT_MAX: break;
    }  /* switch satellite */
  }  /* for iband for coefficients */

  this->cos_sun_zen = cos(input_meta->sun_zen);
  jdoy = input_meta->acq_date.doy;

  /* Use the earth sun distance from the XML file if it was available otherwise
     calculate from the DSUN table */
  if (input_meta->use_toa_refl_consts) {
    dsun = input_meta->earth_sun_dist;
  }
  else { /* calculate from the table */
    for (i = 0; i < (NDSUN - 1); i++) {
      if (jdoy >= dsun_table[i].doy  &&  jdoy <= dsun_table[i+1].doy)
        break;
    }

    if (i >= (NDSUN - 1)) {
      free(this);
      RETURN_ERROR("finding sun-earth distance", "GetLut", false);
    }

    dsun = dsun_table[i].dsun + 
           ((jdoy - dsun_table[i].doy) * 
            ((dsun_table[i + 1].dsun - dsun_table[i].dsun) /
         (dsun_table[i + 1].doy - dsun_table[i].doy)));
  }
  this->dsun2 = dsun * dsun;

  if (input_meta->inst == INST_MSS) {
    for (ib = 0; ib < nband; ib++) {
      this->refl_conv[ib] = (PI * this->dsun2) /
        (this->esun[ib] * this->cos_sun_zen);  
    }
  }

  this->long_name_prefix_ref = DupString(LONG_NAME_PREFIX_REF);
  if (this->long_name_prefix_ref == NULL) {
    free(this);
    RETURN_ERROR("duplicating long name prefix ref", "GetLut", NULL);
  }

  this->units_ref = DupString(UNITS_REF);
  if (this->units_ref == NULL) {
    free(this);
    RETURN_ERROR("duplicating ref units", "GetLut", NULL);
  }

  this->scale_factor_ref =      get_scale_refl();
  this->mult_factor_ref =       1 / this->scale_factor_ref;
  this->add_offset_ref =        get_offset_refl();
  this->valid_range_ref[0] =    (VALID_MIN_REF - this->add_offset_ref) *
                                this->mult_factor_ref;
  this->valid_range_ref[1] =    (VALID_MAX_REF - this->add_offset_ref) *
                                this->mult_factor_ref;
  if (this->valid_range_ref[0] < 0)
      this->valid_range_ref[0] = 0;
  if (this->valid_range_ref[1] > USHRT_MAX)
      this->valid_range_ref[1] = USHRT_MAX;

  this->long_name_prefix_th = DupString(LONG_NAME_PREFIX_TH);
  if (this->long_name_prefix_th == NULL) {
    free(this);
    RETURN_ERROR("duplicating long name prefix th", "GetLut", NULL);
  }

  this->units_th = DupString(UNITS_TH);
  if (this->units_th == NULL) {
    free(this);
    RETURN_ERROR("duplicating th units", "GetLut", NULL);
  }

  this->scale_factor_th =      get_scale_therm();
  this->mult_factor_th =       1 / this->scale_factor_th;
  this->add_offset_th =        get_offset_therm();
  this->valid_range_th[0] =    (VALID_MIN_TH - this->add_offset_th) *
                               this->mult_factor_th;
  this->valid_range_th[1] =    (VALID_MAX_TH - this->add_offset_th) *
                               this->mult_factor_th;
  if (this->valid_range_th[0] < 0)
      this->valid_range_th[0] = 0;
  if (this->valid_range_th[1] > USHRT_MAX)
      this->valid_range_th[1] = USHRT_MAX;

  return this;
}

bool FreeLut(Lut_t *this) {
  free(this);
  return true;
}
