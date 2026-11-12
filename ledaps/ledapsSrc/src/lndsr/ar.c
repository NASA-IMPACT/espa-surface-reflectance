/***************************************************************
03-08-2005: add other gases to band 7 correction
08-01-2013: Modified divide by 10000. to multiply by 0.0001 since
            that is faster.  Gail Schmidt, USGS EROS LSRD
***************************************************************/
#include <stdlib.h>
#include "lndsr.h"
#include "ar.h"
#include "const.h"
#include "error.h"
#include "sixs_runs.h"

#define AOT_MIN_NB_SAMPLES 100

#ifndef  HPUX
#define chand chand_
#define csalbr csalbr_
#endif
void chand(float *phi, float *muv, float *mus, float *tau_ray,
           float *actual_rho_ray);
void csalbr(float *tau_ray,float *actual_S_r);
int compute_aot(float toarhoblue, float toarhored, float ts, float tv,
    float phi, float uoz, float uwv, float spres, sixs_tables_t *sixs_tables,
    float *aot);

static int cmp_collects(const void *p1, const void *p2)
{
    const collect_bands_t *a = p1;
    const collect_bands_t *b = p2;
    if (a->b[0] < b->b[0])
        return -1;
    if (a->b[0] > b->b[0])
        return 1;
    return 0;
}


bool Ar(int il_ar, Lut_t *lut, Img_coord_int_t *size_in, uint16_t ***line_in,
        char **ddv_line, atmos_t *atmos_coef_ar, collect_bands_t *cbands,
        int **line_ar, Ar_stats_t *ar_stats, Ar_gridcell_t *ar_gridcell,
        sixs_tables_t *sixs_tables)
{
    /***
      ddv_line contains results of cloud_screening when this routine is called
      bit 2 = adjacent clouds 1=yes 0=no
      bit 3 = fill value 1=fill 0=valid
      bit 4 = land/water mask 1=land 0=water
      bit 5 = cloud 0=clear 1=cloudy
      bit 6 = cloud shadow 
      bit 7 = snow

      The DDV flag in ddv_line (bit 0) is updated in this routine

     ***/
    int is, il,i;
    int is_ar;
    int is_start, is_end;
    int ib;
    int collect_nbsamps;

    int nb_all_pixs,nb_water_pixs,nb_fill_pixs,nb_cld_pixs,nb_cldshadow_pixs,
        nb_snow_pixs;
    float fraction_water,fraction_clouds;
    bool is_fill;
    int n,water;

    double avg_band[3],std_band[3];
    float avg_srefl,std_srefl;
    float fts,ftv,phi;
    float uwv,uoz,spres;
    float avg_aot;
    int start_i;

    float T_h2o_b7,T_g_b7,rho7,rho4,MP;
    float a_h2o_b7=-3.7338, b_h2o_b7=0.76348,c_h2o_b7=-.030233;
    float a_CO2_b7=0.0071958, b_CO2_b7=0.55665;
    float a_NO2_b7=0.0013383, b_NO2_b7=0.95109;
    float a_CH4_b7=0.030172, b_CH4_b7=0.79652;


    float rho;
    int nb_negative_red,nb_red_obs,ipt;

    /* Do for each region along a line */

    int iipt = il_ar*ar_gridcell->nbcols;
    for (is_start = 0, is_ar = 0; 
        is_start < size_in->s; 
        is_start += lut->ar_region_size.s, is_ar++, iipt++) {

        is_end = is_start + lut->ar_region_size.s - 1;
        if (is_end >= size_in->s) is_end = size_in->s - 1;

        n = 0;

        collect_nbsamps=0;

        fts=ar_gridcell->line_sun_zen[is_ar];
        ftv=ar_gridcell->line_view_zen[is_ar];
        phi=ar_gridcell->line_rel_az[is_ar];
        uwv=ar_gridcell->line_wv[is_ar];
        uoz=ar_gridcell->line_ozone[is_ar];
        spres=ar_gridcell->line_spres[is_ar];

        /**
          compute wv transmittance for band 7
         **/
        MP=(1./cos(fts/DEG)+1./cos(ftv/DEG));
        T_h2o_b7=log(MP*uwv);
        T_h2o_b7=a_h2o_b7+b_h2o_b7*T_h2o_b7+c_h2o_b7*T_h2o_b7*T_h2o_b7;
        T_h2o_b7=exp(-exp(T_h2o_b7));
        T_g_b7=T_h2o_b7;
        T_g_b7 /= (1.+a_CO2_b7*pow(MP,b_CO2_b7));
        T_g_b7 /= (1.+a_NO2_b7*pow(MP,b_NO2_b7));
        T_g_b7 /= (1.+a_CH4_b7*pow(MP,b_CH4_b7));


        nb_all_pixs=0;
        nb_water_pixs=0;
        nb_cld_pixs=0;
        nb_cldshadow_pixs=0;
        nb_snow_pixs=0;
        nb_fill_pixs=0;
        for (il = 0; il < lut->ar_region_size.l; il++) {
            for (is = is_start; is < (is_end + 1); is++) {
                nb_all_pixs++;
                if (ddv_line[il][is] & AR_FILL) {  /* fill */
                    nb_fill_pixs++;
                    continue;
                }

                is_fill = false;

                water = ((ddv_line[il][is] & AR_WATER)==0);  /* water */
                if (water) {
                    nb_water_pixs++;
                    is_fill= true; 
                }

                /***
                  Exclude clouds, cloud shadow & snow pixels flagged by
                  the internal cloud mask.
                 ***/
                if ((ddv_line[il][is] & AR_CLOUD_OR_ADJ_CLOUD) != 0)
                { /* clouds or adjacent clouds */
                    is_fill=true;
                    nb_cld_pixs++;
                }
                if ((ddv_line[il][is] & AR_CLOUD_SHADOW) != 0)
                { /* cloud shadow */
                    is_fill=true;
                    nb_cldshadow_pixs++;
                }
                if ((ddv_line[il][is] & AR_SNOW) != 0)
                { /* snow */
                    is_fill=true;
                    nb_snow_pixs++;
                }

                if (is_fill)
                    continue;

                /* band 7 water vapor correction */
                rho7=line_in[il][AR_B7][is] * lut->scale_factor +
                    lut->add_offset;
                rho4=line_in[il][AR_B4][is] * lut->scale_factor +
                    lut->add_offset;
                rho7 /= T_g_b7;  /* correct for water vapor and other gases*/

                /* update sums if dark target, dark target if not water
                   and 0.015 < rho7 < 0.05 */
                ddv_line[il][is] &= 0xfe;  /* unset DDV bit 0 */

                if ((rho7>0.015) && (rho4 >0.10) /* &&(rho7<0.05) */) {
                    n++;
                    for (ib=0;ib<3;ib++) {
                        cbands[collect_nbsamps].b[ib] = line_in[il][ib][is];
                    }

                    /* Allow the result to go outside the expected
                       range since it is only used for statistics, not an
                       output pixel value. */
                    cbands[collect_nbsamps].b7 = (unsigned short)
                        ((rho7 - lut->add_offset) * lut->mult_factor);

                    collect_nbsamps++;
                    if (rho7<0.05)
                        ddv_line[il][is] |= AR_DDV; /* set DDV bit */
                }
             }  /* end for is */
        }  /* end for il */
    
        if (collect_nbsamps == 0) {
            line_ar[0][is_ar] = lut->aerosol_fill;
            line_ar[1][is_ar] = lut->aerosol_fill;
            line_ar[2][is_ar] = lut->aerosol_fill;
            ar_stats->nfill++;
        } else {

            /**
               Sort collected observations
            **/
            qsort(cbands, collect_nbsamps, sizeof(*cbands), cmp_collects);

            if (collect_nbsamps >= 2*AOT_MIN_NB_SAMPLES) {
                start_i=AOT_MIN_NB_SAMPLES;
                collect_nbsamps=AOT_MIN_NB_SAMPLES; /* Take the first
                                                       AOT_MIN_NB_SAMPLES
                                                       samples only */
                double sum_band[3] = {0, 0, 0},
                       sum_band_sq[3] = {0, 0, 0};
                for (ib=0;ib<3;ib++) {
                    for (i=0;i<collect_nbsamps;i++) {
                        double dtemp = cbands[i+start_i].b[ib] 
                            * lut->scale_factor + lut->add_offset;
                        sum_band[ib] += dtemp;
                        sum_band_sq[ib] += (dtemp * dtemp);
                    }
                }
                double sum_srefl = 0, sum_srefl_sq = 0;
                for (i=0;i<collect_nbsamps;i++) {
                    double dtemp = cbands[i+start_i].b7 * lut->scale_factor
                        + lut->add_offset;
                    sum_srefl += dtemp;
                    sum_srefl_sq += (dtemp * dtemp);
                }

                /* update stats line */
      avg_srefl = (sum_srefl) / collect_nbsamps; 
                for (ib=0;ib<3;ib++)
			avg_band[ib]=sum_band[ib]/collect_nbsamps; 
            
                if (collect_nbsamps>1) {
                    std_srefl= (sum_srefl_sq -
                                collect_nbsamps*avg_srefl*avg_srefl)
                             / (collect_nbsamps - 1);
                    if (std_srefl>0)
                        std_srefl=sqrt(std_srefl);
                    else 
                        std_srefl=0;
                    for (ib=0;ib<3;ib++) {
                        std_band[ib] = (sum_band_sq[ib] -
                                        collect_nbsamps*avg_band[ib]*
                                        avg_band[ib])
                                     / (collect_nbsamps - 1);
                        if (std_band[ib]>0)
                            std_band[ib]=sqrt(std_band[ib]);
                        else 
                            std_band[ib]=0;
                    }
                } else {
                    std_srefl=0.;
                    for (ib=0;ib<3;ib++) {
                        std_band[ib]=0;
                    }
                }
                fraction_water=(float)nb_water_pixs/(nb_all_pixs-nb_fill_pixs);
                fraction_clouds=(float)nb_cld_pixs/(nb_all_pixs-nb_fill_pixs);
        
/**
   Compute AOT blue band
***/

                if (std_srefl <= 1.015 && avg_srefl <= 0.15 &&
                    nb_snow_pixs < 5 && fraction_water < 0.3 &&
                    fraction_clouds < 1e-10) {
                
                    compute_aot(avg_band[0], avg_band[2], fts, ftv, phi,
                                uoz, uwv, spres, sixs_tables, &avg_aot);
        
                    line_ar[0][is_ar] = (int)(avg_aot*1000.);

/***
    Filter aot : Correct red band using retreived aot. if over 30% of the
    corrected refelctances are negative, reject aot.
***/
                    update_gridcell_atmos_coefs(iipt, atmos_coef_ar,
                                                ar_gridcell, sixs_tables,
                                                line_ar[0][is_ar],
                                                lut, 6, 0);

                    ib=2; /*  test with red band */
                    nb_red_obs=0;
                    nb_negative_red=0;
                    ipt=il_ar*lut->ar_size.s+is_ar;
                    for (il = 0; il < lut->ar_region_size.l; il++) {
                        for (is = is_start; is < (is_end + 1); is++) {
                            if (!(ddv_line[il][is] & AR_FILL)) {  /* fill */
                                rho7 = (float)line_in[il][AR_B7][is] 
                                    * lut->scale_factor + lut->add_offset;
                                rho7 /= T_g_b7;  /* correct for water vapor
                                                    and other gases*/
                                rho = compute_rho(line_in[il][ib][is] 
                                               * lut->scale_factor 
                                               + lut->add_offset,
                                               atmos_coef_ar->tgOG[ib][ipt],
                                               atmos_coef_ar->tgH2O[ib][ipt],
                                               atmos_coef_ar->td_ra[ib][ipt],
                                               atmos_coef_ar->tu_ra[ib][ipt],
                                               atmos_coef_ar->rho_ra[ib][ipt],
                                               atmos_coef_ar->S_ra[ib][ipt]);
                                nb_red_obs++;
            
                                if (rho < 0. || rho > rho7)
                                    /*eric introduced that to get rid of the
                                      salt pan */
                                    nb_negative_red++;
                            }
                        }
                    }
                    if (((float)nb_negative_red/(float)nb_red_obs) > 0.01) {
                        line_ar[0][is_ar]=lut->aerosol_fill;
                    }

                    if (ar_stats->first) {

                        ar_stats->ar_min = ar_stats->ar_max = line_ar[0][is_ar];
                        ar_stats->first = false;

                    } else {

                        if (line_ar[0][is_ar] < ar_stats->ar_min)
                            ar_stats->ar_min = line_ar[0][is_ar];

                        if (line_ar[0][is_ar] > ar_stats->ar_max)
                            ar_stats->ar_max = line_ar[0][is_ar];
                    }
                } else {
                    line_ar[0][is_ar] = lut->aerosol_fill;
                    line_ar[1][is_ar] = lut->aerosol_fill;
                    line_ar[2][is_ar] = lut->aerosol_fill;
                }
            } else {
                line_ar[0][is_ar] = lut->aerosol_fill;
                line_ar[1][is_ar] = lut->aerosol_fill;
                line_ar[2][is_ar] = lut->aerosol_fill;
                ar_stats->nfill++;
            }
        }
    }

    return true;
}


void aot_correct_band(int band, float toarho, float phi, float mus, float muv,
    float tau_ray, sixs_tables_t *sixs_tables, float surrho[SIXS_NB_AOT])
{
    int i;
    float actual_rho_ray,actual_T_ray,actual_S_r;

    chand(&phi,&muv,&mus,&tau_ray,&actual_rho_ray);

    actual_T_ray = (2./3. + mus + (2./3. - mus)*exp(-tau_ray/mus))
                 / (4./3. + tau_ray); /* downward */
    actual_T_ray *= (2./3. + muv + (2./3. - muv)*exp(-tau_ray/muv))
                  / (4./3. + tau_ray); /* total */

    csalbr(&tau_ray,&actual_S_r);

    for (i=0;i<SIXS_NB_AOT;i++) {
        surrho[i] = toarho/sixs_tables->T_g_og[band];
        surrho[i] -= actual_rho_ray + sixs_tables->rho_ra[band][i]
                      - sixs_tables->rho_r[band];
        surrho[i] *= sixs_tables->T_a[band][i]/actual_T_ray
                     * sixs_tables->T_g_wv[band];
        surrho[i] /= 1 + (actual_S_r + sixs_tables->S_ra[band][i] -
                             sixs_tables->S_r[band])*surrho[i];
    }
}


int compute_aot(float toarhoblue, float toarhored, float ts, float tv,
    float phi, float uoz, float uwv, float spres, sixs_tables_t *sixs_tables,
    float *aot)
{
    int i,iaot;
    float minimum,temp,eratio;
    float slope;
    float surrhoblue[SIXS_NB_AOT],surrhored[SIXS_NB_AOT];
        float temp1,temp2;
    float mus,muv,tau_ray,ratio;
    float tau_ray_sealevel[7] = {0.16511,0.08614,0.04716,0.01835,0.00113,
                                 0.00037}; /* bands 1,2,3,4,5,7 */
    
/* correct the blue band */ 
    mus=cos(ts*RAD);
    muv=cos(tv*RAD);
    ratio=spres/1013.;
    tau_ray=tau_ray_sealevel[AR_B1]*ratio;
    aot_correct_band (AR_B1, toarhoblue, phi, mus, muv, tau_ray, sixs_tables,
        surrhoblue);

    /* correct the red band */  
    tau_ray=tau_ray_sealevel[AR_B3]*ratio;
    aot_correct_band (AR_B3, toarhored, phi, mus, muv, tau_ray, sixs_tables,
        surrhored);
    
/* compute aot based on empirical ratio */  
    minimum=9999999;
    eratio=0.66;
    iaot=-1;
    for (i=0;i<SIXS_NB_AOT;i++) {
        if ( surrhoblue[i] > 0. ) {
            /* Find index of the wavelength with temp nearest 0 */
            temp=surrhoblue[i]-eratio*surrhored[i];
            if (fabs(temp) < minimum ) {
                minimum=temp;
                iaot=i;
        }
         }
    
    }
    if (iaot==-1)  {
        *aot=0.01;
    }
    else
    {
        if (iaot==0) iaot=1;
        temp1=surrhoblue[iaot-1]-eratio*surrhored[iaot-1];
        temp2=surrhoblue[iaot]-eratio*surrhored[iaot];
        /* Calculate slope of line with wavelength as function of temp. */
        slope = (sixs_tables->aot_wavelength[1][iaot] -
                 sixs_tables->aot_wavelength[1][iaot-1])/(temp2 - temp1);

        /* let's hope all the aot table are independent of wavelength but
           that need to be checked */        
        /* Calculate wavelength where temp = 0 */
        *aot = sixs_tables->aot_wavelength[1][iaot-1] - temp1*slope;
    }
    if (*aot < 0.01)
        *aot=0.01;

    return 0;       
}

int ArInterp(Lut_t *lut, Img_coord_int_t *input_loc, int ***line_ar,
             int *inter_aot) 
/* 
  Point order:

    0 ---- 1    +--> sample
    |      |    |
    |      |    v
    2 ---- 3   line

 */

{
    Img_coord_int_t p[4];
    int i, n;
    float dl, ds, w[4];
    float grid_line, grid_sample;
    float sum_w;
    Img_coord_int_t ar_region_half;

    *inter_aot = lut->aerosol_fill;

    /* Note the right shift by 1 is a faster way of divide by 2 */
    ar_region_half.l = (lut->ar_region_size.l + 1) >> 1;
    ar_region_half.s = (lut->ar_region_size.s + 1) >> 1;

    /* Set the four corner point indices. */
    /* AOT values are assumed to be for the center of the
     * ar_region_size.l x ar_region_size.s cell.
     * Find the cells whose centers surround this point */
    grid_line = (float)(input_loc->l - ar_region_half.l)
              / lut->ar_region_size.l;
    grid_sample = (float)(input_loc->s - ar_region_half.s)
                / lut->ar_region_size.s;
    p[0].l = (int)grid_line;
    p[2].l = p[0].l + 1;
    if (p[2].l >= lut->ar_size.l) {
        p[2].l = lut->ar_size.l - 1;
        if (p[0].l > 0)
            p[0].l--;
    }
    p[1].l = p[0].l;
    p[3].l = p[2].l;

    p[0].s = (int)grid_sample;
    p[1].s = p[0].s + 1;
    if (p[1].s >= lut->ar_size.s) {
        p[1].s = lut->ar_size.s - 1;
        if (p[0].s > 0)
            p[0].s--;
    }
    p[2].s = p[0].s;
    p[3].s = p[1].s;

    /* Initialize the counting and sum variables */
    n = 0;
    sum_w = 0;

    /* Compute the fractional grid cell offset of the sample point from
       the upper-left corner of the cell. */
    dl = grid_line - p[0].l;
    /* Lines above the 0th cell are assumed to have distance 0 to it,
     * effectively putting all the weight on w[0] and w[1] */
    if (dl < 0) dl = 0.0;
    ds = grid_sample - p[0].s;
    /* Likewise for samples to the left of the 0th cell */
    if (ds < 0) ds = 0.0;
    w[0] = (1 - dl)*(1 - ds);
    w[1] = 1 - dl - w[0];
    w[2] = 1 - ds - w[0];
    w[3] = ds - w[1];;

    float val = 0;
    for (i = 0; i < 4; i++) {
        if (line_ar[p[i].l][0][p[i].s] != lut->aerosol_fill)
        {
            val += line_ar[p[i].l][0][p[i].s]*w[i];
            sum_w += w[i];
            n++;
        }
    }

    if (n > 0)
    {
        if (n < 4)
            *inter_aot = rint(val/sum_w);
        else
            *inter_aot = rint(val);
    }

    return 0;
}

int Old_Fill_Ar_Gaps(Lut_t *lut, int ***line_ar, int ib) 
/* 
  Point order:

    0 ---- 1    +--> sample
    |      |    |
    |      |    v
    2 ---- 3   line

 */

{
    Img_coord_int_t *gaps_loc;
    int *gaps_neighbors,nb_gaps;
    int i, j, i_aot,j_aot,n,more_gaps;
    float w, sum, sum_w;


    /**
       allocate memory for gaps location and number of neighbors
    **/
    gaps_loc = malloc(lut->ar_size.l*lut->ar_size.s*sizeof(Img_coord_int_t));
    if (gaps_loc == NULL)
        EXIT_ERROR("failed to allocate memory for gaps_loc","Fill_Ar_Gaps");
    gaps_neighbors=(int *)malloc(lut->ar_size.l*lut->ar_size.s*sizeof(int));
    if (gaps_neighbors == NULL)
        EXIT_ERROR("failed to allocate memory for gaps_neighbors",
                   "Fill_Ar_Gaps");

    more_gaps=1;
    while (more_gaps) {
        more_gaps=0;
        /**
           get location of gaps and number of their valid neighbors
        **/
        nb_gaps=0;
        for (i=0;i<lut->ar_size.l;i++)
            for (j=0;j<lut->ar_size.s;j++)
                if (line_ar[i][ib][j] == lut->aerosol_fill) {
                    /* check aot in the red band since it is used to fill the
                       others */
                    gaps_loc[nb_gaps].l=i;
                    gaps_loc[nb_gaps].s=j;
                    gaps_neighbors[nb_gaps]=0;
                    for (i_aot=(i-1);i_aot<=(i+1);i_aot++) 
                        if ((i_aot>=0)&&(i_aot<lut->ar_size.l))
                            for (j_aot=(j-1);j_aot<=(j+1);j_aot++) 
                                if ((j_aot>=0)&&(j_aot<lut->ar_size.s))
                                    if (line_ar[i_aot][ib][j_aot] !=
                                        lut->aerosol_fill)
                                        gaps_neighbors[nb_gaps]++;
                    nb_gaps++;
                }

        /**
           sort list of gaps by decreasing number of neighbors
        **/
        for (i=0;i<(nb_gaps-1);i++) 
            for (j=(i+1);j<nb_gaps;j++)
                if (gaps_neighbors[i]<gaps_neighbors[j]) {
                    n=gaps_neighbors[j];
                    gaps_neighbors[j]=gaps_neighbors[i];
                    gaps_neighbors[i]=n;
                    n=gaps_loc[j].l;
                    gaps_loc[j].l=gaps_loc[i].l;
                    gaps_loc[i].l=n;
                    n=gaps_loc[j].s;
                    gaps_loc[j].s=gaps_loc[i].s;
                    gaps_loc[i].s=n;
                }

        for (i=0;i<nb_gaps;i++) { 
            n=0;
            sum=0.;
            sum_w=0.;
            int i_start, i_end, j_start, j_end;
            if (gaps_loc[i].l > 0)
                i_start = -1;
            else
                i_start = 0;
            if (gaps_loc[i].l < lut->ar_size.l - 1)
                i_end = 1;
            else
                i_end = 0;
            if (gaps_loc[i].s > 0)
                j_start = -1;
            else
                j_start = 0;
            if (gaps_loc[i].s < lut->ar_size.s - 1)
                j_end = 1;
            else
                j_end = 0;
            int lindex = gaps_loc[i].l + i_start;
            int sindex = gaps_loc[i].s + j_start;
            for (i_aot = i_start; i_aot <= i_end; i_aot++, lindex++)
                for (j_aot = j_start; j_aot <= j_end; j_aot++, sindex++)
                    if (line_ar[lindex][ib][sindex]
                        != lut->aerosol_fill) {
                        w = (3 - abs(i_aot))*(3 - abs(j_aot));
                        n++;
                        sum_w += w;
                        sum += line_ar[lindex][ib][sindex]*w;
                    }

            if (n>0) {
                line_ar[gaps_loc[i].l][ib][gaps_loc[i].s] = rint(sum/sum_w);
                if (line_ar[gaps_loc[i].l][ib][gaps_loc[i].s] < 20)
                    line_ar[gaps_loc[i].l][ib][gaps_loc[i].s]=20;
            } else
                more_gaps=1;
        }
    }
  
    free(gaps_loc);
    free(gaps_neighbors);
    return 0;
}


int Fill_Ar_Gaps(Lut_t *lut, int ***line_ar, int ib) {
/*
!Description: fill in missing values in the OZONE grid based
on existing values (spatial interpolation). 

!END****************************************************************************
*/
   int i,j,k,l,count,More_Gaps,nbfills;
   int last_value=-99;
   float dist,sum_value,sum_dist;
   char **missing_flag,*missing_flag_array;
   int min_nb_values,n,max_distance;

    if ((missing_flag_array = calloc(lut->ar_size.l*lut->ar_size.s,
                                     sizeof(char))) == NULL)
        return 0;
    if ((missing_flag=(char **)calloc(lut->ar_size.l,sizeof(char *)))==NULL)
        return 0;
    for (i=0;i<lut->ar_size.l;i++)
        missing_flag[i]=&missing_flag_array[i*lut->ar_size.s];
/**
Start by counting valid values
if nb gaps = 0 do nothing
if nb gaps = 1 duplicate value everywhere
**/
    count=0;
    for (i=0;i<lut->ar_size.l;i++) {
        for (j=0;j<lut->ar_size.s;j++) {
            if (line_ar[i][ib][j] != lut->aerosol_fill)  {
                count++;
                last_value=line_ar[i][ib][j];
            }
        }
    } 
    if (count==0)
      return 0;
    if (count==1) {
        for (i=0;i<lut->ar_size.l;i++)
            for (j=0;j<lut->ar_size.s;j++) {
                line_ar[i][ib][j]=last_value;
            }
        return 0;
    }
    More_Gaps=1;
    nbfills=1;
    while (More_Gaps && nbfills != 0) {
        More_Gaps=0;
        nbfills=0;
        for (i=0;i<lut->ar_size.l;i++) {
            for (j=0;j<lut->ar_size.s;j++) {
                missing_flag[i][j]=0;
                if (line_ar[i][ib][j]==lut->aerosol_fill)  {
                    missing_flag[i][j]=1;
                    More_Gaps=1;
                }
            }
        }
        if (!More_Gaps)
            break;

        for (i=0;i<lut->ar_size.l;i++) {
            for (j=0;j<lut->ar_size.s;j++) {
                /**
                   Look for at least 3 neighboring valid values within
                   4 GPs smooth
                **/
                min_nb_values=3;
                max_distance=3; 
                if (missing_flag[i][j]) {
                    sum_dist=0.;
                    sum_value=0.;
                    n=0;
                    int k_start, k_end;
                    if (i > max_distance)
                        k_start = i - max_distance;
                    else
                        k_start = 0;
                    if (i < lut->ar_size.l - max_distance - 1)
                        k_end = i + max_distance;
                    else
                        k_end = lut->ar_size.l - 1;
                    for (k = k_start; k <= k_end; k++) {
                        int l_start, l_end;
                        if (j > max_distance)
                            l_start = j - max_distance;
                        else
                            l_start = 0;
                        if (j < lut->ar_size.s - max_distance - 1)
                            l_end = j + max_distance;
                        else
                            l_end = lut->ar_size.s - 1;
                        for (l = l_start; l <= l_end; l++) {
                            if (!missing_flag[k][l]) {
                                dist=sqrt((k-i)*(k-i)+(l-j)*(l-j));
                                sum_dist += dist;
                                sum_value += dist*line_ar[k][ib][l];
                                n++;
                            }
                        }
                    }
                    if (n >= min_nb_values && sum_dist != 0) {
                        line_ar[i][ib][j]=sum_value/sum_dist;
                        nbfills++;
                    }
                } /* fill missing */
            } /* sample loop */
        } /* line loop */
    } /* gaps loop */

    if (More_Gaps && nbfills == 0)
        for (i=0;i<lut->ar_size.l;i++)
            for (j=0;j<lut->ar_size.s;j++)
                if (line_ar[i][ib][j] == lut->aerosol_fill)
                    line_ar[i][ib][j]=60;

    return 0;
}
