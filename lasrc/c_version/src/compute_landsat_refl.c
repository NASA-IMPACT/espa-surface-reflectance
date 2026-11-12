/******************************************************************************
FILE: compute_landsat_refl.c

PURPOSE: Contains functions for handling the Landsat 8/9 TOA reflectance and
surface reflectance corrections.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
******************************************************************************/

#include "time.h"
#include "aero_interp.h"
#include "poly_coeff.h"
#include "read_level1_qa.h"
#include "read_level2_qa.h"
//#define WRITE_TAERO 1

/******************************************************************************
MODULE:  compute_landsat_toa_refl

PURPOSE:  Computes the TOA reflectance and TOA brightness temps for all
the Landsat bands except the pan band. Uses a per-pixel solar zenith angle for
the TOA corrections.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error computing the reflectance
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
  1. These TOA and BT algorithms match those as published by the USGS Landsat
     team in http://landsat.usgs.gov/Landsat8_Using_Product.php
******************************************************************************/
int compute_landsat_toa_refl
(
    Input_t *input,     /* I: input structure for the Landsat product */
    Espa_internal_meta_t *xml_metadata,
                        /* I: XML metadata structure */
    uint16 *qaband,     /* I: QA band for the input image, nlines x nsamps */
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    char *instrument,   /* I: instrument to be processed (OLI, TIRS) */
    int16 *sza,         /* I: scaled per-pixel solar zenith angles (degrees),
                              nlines x nsamps */
    float **sband       /* O: output TOA reflectance and brightness temp
                              values (unscaled) */
)
{
    char errmsg[STR_SIZE];                   /* error message */
    char FUNC_NAME[] = "compute_landsat_toa_refl";   /* function name */
    int i;               /* looping variable for pixels */
    int ib;              /* looping variable for input bands */
    int th_indx;         /* index of thermal band and K constants */
    int sband_ib;        /* looping variable for output bands */
    int iband;           /* current band */
    long npixels;        /* number of pixels to process */
    float rotoa;         /* top of atmosphere reflectance */
    float tmpf;          /* temporary floating point value */
    float refl_mult;     /* reflectance multiplier for bands 1-9 */
    float refl_add;      /* reflectance additive for bands 1-9 */
    float xcals;         /* radiance multiplier for bands 10 and 11 */
    float xcalo;         /* radiance additive for bands 10 and 11 */
    float k1;            /* K1 temperature constant for thermal bands */
    float k2;            /* K2 temperature constant for thermal bands */
    float xmus;          /* cosine of solar zenith angle (per-pixel) */
    float sza_mult;      /* sza gain value (for unscaling) */
    float sza_add;       /* sza offset value (for unscaling) */
    uint16 *uband = NULL;  /* array for input image data for a single band,
                              nlines x nsamps */
    time_t mytime;       /* time variable */

    /* Start the processing */
    mytime = time(NULL);
    printf ("Start TOA reflectance corrections: %s", ctime(&mytime));

    /* Allocate memory for band data */
    npixels = nlines * nsamps;
    uband = calloc (npixels, sizeof (uint16));
    if (uband == NULL)
    {
        sprintf (errmsg, "Error allocating memory for uband");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through all the bands (except the pan band) and compute the TOA
       reflectance and TOA brightness temp */
    for (ib = DNL_BAND1; ib <= DNL_BAND11; ib++)
    {
        /* Don't process the pan band */
        if (ib == DNL_BAND8)
            continue;
        printf ("%d ... ", ib+1);

        /* Read the current band and calibrate bands 1-9 (except pan) to
           obtain TOA reflectance. Bands are corrected for the sun angle. */
        if (ib <= DNL_BAND9)
        {
            if (ib <= DNL_BAND7)
            {
                iband = ib;
                sband_ib = ib;
            }
            else
            {  /* don't count the pan band */
                iband = ib - 1;
                sband_ib = ib - 1;
            }

            if (get_input_refl_lines (input, iband, 0, nlines, nsamps, uband)
                != SUCCESS)
            {
                sprintf (errmsg, "Error reading Landsat band %d", ib+1);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Get TOA reflectance coefficients for this reflectance band from
               XML file */
            refl_mult = input->meta.gain[iband];
            refl_add = input->meta.bias[iband];
            sza_mult = input->meta.gain_sza;
            sza_add = input->meta.bias_sza;

#ifdef _OPENMP
            #pragma omp parallel for private (i, xmus, rotoa)
#endif
            for (i = 0; i < npixels; i++)
            {
                /* If this pixel is fill, continue with the next pixel. */
                if (level1_qa_is_fill(qaband[i]))
                {
                    sband[sband_ib][i] = FILL_VALUE;
                    continue;
                }

                /* Compute the TOA reflectance based on the per-pixel
                   sun angle (need to unscale the DN value) */
                xmus = cos((sza[i] * sza_mult + sza_add) * DEG2RAD);
                rotoa = (uband[i] * refl_mult) + refl_add;
                rotoa /= xmus;

                /* Save the scaled TOA reflectance value, but make
                   sure it falls within the defined valid range since it will
                   get used for SR computations */
                if (rotoa < MIN_VALID_REFL)
                    sband[sband_ib][i] = MIN_VALID_REFL;
                else if (rotoa > MAX_VALID_REFL)
                    sband[sband_ib][i] = MAX_VALID_REFL;
                else
                    sband[sband_ib][i] = rotoa;
            }  /* for i */
        }  /* end if band <= band 9 */

        /* Read the current band and calibrate thermal bands.  Not available
           for OLI-only scenes. */
        else if ((ib == DNL_BAND10 || ib == DNL_BAND11) &&
                 strcmp (instrument, "OLI"))
        {
            /* Handle index differences between bands */
            if (ib == DNL_BAND10)
            {
               th_indx = 0;
               sband_ib = SRL_BAND10;
            }
            else
            {  /* if (ib == DNL_BAND11) */
               th_indx = 1;
               sband_ib = SRL_BAND11;
            }

            /* Read the input thermal lines */
            if (get_input_th_lines (input, th_indx, 0, nlines, uband)
                != SUCCESS)
            {
                sprintf (errmsg, "Reading band %d", ib+1);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Get brightness temp coefficients for this band from XML file */
            xcals = input->meta.gain_th[th_indx];
            xcalo = input->meta.bias_th[th_indx];
            k1 = input->meta.k1_const[th_indx];
            k2 = input->meta.k2_const[th_indx];

            /* Compute brightness temp for band 10.  Make sure it falls
               within the min/max range for the thermal bands. */
#ifdef _OPENMP
            #pragma omp parallel for private (i, tmpf)
#endif
            for (i = 0; i < npixels; i++)
            {
                /* If this pixel is fill, continue with the next pixel. */
                if (level1_qa_is_fill(qaband[i]))
                {
                    sband[sband_ib][i] = FILL_VALUE;
                    continue;
                }

                /* Compute the TOA spectral radiance */
                tmpf = xcals * uband[i] + xcalo;

                /* Compute TOA brightness temp (K) */
                tmpf = k2 / log (k1 / tmpf + 1.0);

                /* Make sure the brightness temp falls within the specified
                   range, since it will get used for the SR computations */
                if (tmpf < MIN_VALID_TH)
                    sband[sband_ib][i] = MIN_VALID_TH;
                else if (tmpf > MAX_VALID_TH)
                    sband[sband_ib][i] = MAX_VALID_TH;
                else
                    sband[sband_ib][i] = tmpf;
            }
        }  /* end if band 10 || band 11*/
    }  /* end for ib */
    printf ("\n");

    /* The input data has been read and calibrated. The memory can be freed. */
    free (uband);

    /* Successful completion */
    mytime = time(NULL);
    printf ("End of TOA reflectance computations: %s", ctime(&mytime));
    return (SUCCESS);
}


/******************************************************************************
MODULE:  compute_landsat_sr_refl

PURPOSE:  Computes the surface reflectance for all the Landsat reflectance
bands.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error computing the reflectance
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
1. Initializes the variables and data arrays from the lookup table and
   auxiliary files.
2. The tauray array was originally read in from a static ASCII file, but it is
   now hardcoded to save time from reading the file each time.  This file was
   generated (like many of the other auxiliary input tables) by running 6S and
   storing the coefficients.
4. Aerosols are retrieved for all non-fill pixels.  If the aerosol fails the
   model residual or NDVI test, then the pixel is flagged as water.  All water
   pixels are run through a water-specific aerosol retrieval.  If the model
   residual fails, then that pixel is marked as failed aerosol retrieval.  Any
   pixel that failed retrieval is then interpolated using an average of the
   clear (valid land pixel aerosols) and water (valid water pixel aerosols).
   Those final aerosol values are used for the surface reflectance corrections.
5. Cloud-based QA information is not processed in this algorithm.
******************************************************************************/
int compute_landsat_sr_refl
(
    Input_t *input,     /* I: input structure for the Landsat product */
    Espa_internal_meta_t *xml_metadata,
                        /* I: XML metadata structure */
    char *xml_infile,   /* I: input XML filename */
    uint16 *qaband,     /* I: QA band for the input image, nlines x nsamps */
    uint16 *out_band,   /* I: allocated array for writing scaled output */
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    float pixsize,      /* I: pixel size for the reflectance bands */
    float **sband,      /* I/O: input TOA (unscaled) and output surface
                                reflectance (unscaled) */
    int16 *sza,         /* I: scaled per-pixel solar zenith angles (degrees),
                              nlines x nsamps */
    int16 *saa,         /* I: scaled per-pixel solar azimuth angles (degrees),
                              nlines x nsamps */
    int16 *vza,         /* I: scaled per-pixel view zenith angles (degrees),
                              nlines x nsamps */
    int16 *vaa,         /* I: scaled per-pixel view azimuth angles (degrees),
                              nlines x nsamps */
    float xts_center,   /* I: scene center solar zenith angle (deg) */
    float xmus_center,  /* I: cosine of scene center solar zenith angle */
    bool use_orig_aero, /* I: use the original aerosol handling if specified,
                              o/w use the semi-empirical approach */
    char *anglehdf,     /* I: angle HDF filename */
    char *intrefnm,     /* I: intrinsic reflectance filename */
    char *transmnm,     /* I: transmission filename */
    char *spheranm,     /* I: spherical albedo filename */
    char *cmgdemnm,     /* I: climate modeling grid (CMG) DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm,        /* I: auxiliary filename for ozone and water vapor */
    char *aux2nm,       /* I: 2nd auxiliary filename for ozone and water vapor (23-FEB-26, JPR) */
    aux_src_t aux_src   /* I: identifies the source of atmospheric aux data */
)
{
    char errmsg[STR_SIZE];  /* error message */
    char FUNC_NAME[] = "compute_landsat_sr_refl";   /* function name */
    Sat_t sat = input->meta.sat;  /* satellite */
    int retval;          /* return status */
    int i, j;            /* looping variable for pixels */
    int ib;              /* looping variable for input bands */
    int iband;           /* current band */
    int curr_pix = -99;  /* current pixel in 1D arrays of nlines * nsamps */
    int center_pix;      /* current pixel in 1D arrays of nlines * nsamps for
                            the center of the aerosol window */
    int center_line;     /* line for the center of the aerosol window */
    int center_samp;     /* sample for the center of the aerosol window */
    int nearest_line;    /* line for nearest non-fill/cloud pixel in the
                            aerosol window */
    int nearest_samp;    /* samp for nearest non-fill/cloud pixel in the
                            aerosol window */
    long npixels;        /* number of pixels to process */
    float tmpf;          /* temporary floating point value */
    float rotoa;         /* top of atmosphere reflectance */
    float roslamb;       /* lambertian surface reflectance */
    float tgo;           /* other gaseous transmittance (tgog * tgoz) */
    float roatm;         /* intrinsic atmospheric reflectance */
    float ttatmg;        /* total atmospheric transmission */
    float satm;          /* atmosphere spherical albedo */
    float tgo_x_roatm;   /* variable for tgo * roatm */
    float tgo_x_ttatmg;  /* variable for tgo * ttatmg */
    float xrorayp;       /* reflectance of the atmosphere due to molecular
                            (Rayleigh) scattering */
    float erelc[NSR_BANDS];   /* band ratio variable for refl bands */
    float troatm[NSR_BANDS];  /* atmospheric reflectance table for refl bands */
    float btgo[NSR_BANDS];    /* other gaseous transmittance for refl bands */
    float broatm[NSR_BANDS];  /* atmospheric reflectance for refl bands */
    float bttatmg[NSR_BANDS]; /* ttatmg for refl bands */
    float bsatm[NSR_BANDS];   /* atmosphere spherical albedo for refl bands */

    int iband1;         /* band index (zero-based) */
    float raot;         /* AOT reflectance */
    float sraot1, sraot3;
                        /* raot values for three different eps values */
    float residual;     /* model residual */
    float residual1, residual2, residual3;
                        /* residuals for 3 different eps values */
    float rsurf;        /* surface reflectance */
    float corf;         /* aerosol impact (higher values represent high
                           aerosol) */
    float ros1, ros4, ros5; /* surface reflectance for bands 1, 4, and 5 */
#ifndef _OPENMP
    int tmp_percent;      /* current percentage for printing status */
    int curr_tmp_percent; /* percentage for current line */
#endif

    float lat, lon;       /* pixel lat, long location */
    int lcmg, scmg;       /* line/sample index for the CMG */
    int lcmg1, scmg1;     /* line+1/sample+1 index for the CMG */
    float u, v;           /* line/sample index for the CMG */
    float one_minus_u;    /* 1.0 - u */
    float one_minus_v;    /* 1.0 - v */
    float one_minus_u_x_one_minus_v;  /* (1.0 - u) * (1.0 - v) */
    float one_minus_u_x_v;  /* (1.0 - u) * v */
    float u_x_one_minus_v;  /* u * (1.0 - v) */
    float u_x_v;          /* u * v */
    float ndwi_th1, ndwi_th2; /* values for NDWI calculations */
    float xcmg, ycmg;     /* x/y location for CMG */
    float xndwi;          /* calculated NDWI value */
    int uoz11, uoz21, uoz12, uoz22;  /* ozone at line,samp; line, samp+1;
                                        line+1, samp; and line+1, samp+1 */
    float pres11, pres12, pres21, pres22;  /* pressure at line,samp;
                             line, samp+1; line+1, samp; and line+1, samp+1 */
    float wv11, wv12, wv21, wv22;  /* water vapor at line,samp;
                             line, samp+1; line+1, samp; and line+1, samp+1 */
    uint8 *ipflag = NULL; /* QA flag to assist with aerosol interpolation,
                             nlines x nsamps */
    float *twvi = NULL;   /* interpolated water vapor value,
                             nlines x nsamps */
    float *tozi = NULL;   /* interpolated ozone value, nlines x nsamps */
    float *tp = NULL;     /* interpolated pressure value, nlines x nsamps */
    float *taero = NULL;  /* aerosol values for each pixel, nlines x nsamps */
    float *teps = NULL;   /* angstrom coeff for each pixel, nlines x nsamps */
    float *aerob1 = NULL; /* atmospherically corrected band 1 data
                             (unscaled TOA refl), nlines x nsamps */
    float *aerob2 = NULL; /* atmospherically corrected band 2 data
                             (unscaled TOA refl), nlines x nsamps */
    float *aerob4 = NULL; /* atmospherically corrected band 4 data
                             (unscaled TOA refl), nlines x nsamps */
    float *aerob5 = NULL; /* atmospherically corrected band 5 data
                             (unscaled TOA refl), nlines x nsamps */
    float *aerob7 = NULL; /* atmospherically corrected band 7 data
                             (unscaled TOA refl), nlines x nsamps */

    /* Vars for forward/inverse mapping space */
    Geoloc_t *space = NULL;       /* structure for geolocation information */
    Space_def_t space_def;        /* structure to define the space mapping */
    Img_coord_float_t img;        /* coordinate in line/sample space */
    Geo_coord_t geo;              /* coordinate in lat/long space */

    /* Lookup table variables */
    float eps;           /* angstrom coefficient */
    float eps1, eps2, eps3;  /* eps values for three runs */
    float xts;           /* solar zenith angle (deg) */
    float xmus;          /* cosine of solar zenith angle */
    float xtv_center;    /* scene center observation zenith angle (deg) */
    float xmuv_center;   /* cosine of scene center observation zenith angle */
    float xtv;           /* observation zenith angle (deg) */
    float xmuv;          /* cosine of observation zenith angle */
    float xfi_center;    /* azimuthal difference between the sun and
                            observation angle at scene center (deg) */
    float cosxfi_center; /* cosine of azimuthal difference at scene center */
    float xfi;           /* azimuthal difference between the sun and
                            observation angle (deg) */
    float cosxfi;        /* cosine of azimuthal difference */
    float xtsstep;       /* solar zenith step value */
    float xtsmin;        /* minimum solar zenith value */
    float xtvstep;       /* observation step value */
    float xtvmin;        /* minimum observation value */
    float *rolutt = NULL;  /* intrinsic reflectance table
                          [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float *transt = NULL;  /* transmission table
                       [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUNANGLE_VALS] */
    float *sphalbt = NULL; /* spherical albedo table 
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext = NULL; /* aerosol extinction coefficient at the current
                              wavelength (normalized at 550nm)
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *tsmax = NULL;   /* maximum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *tsmin = NULL;   /* minimum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfi = NULL;    /* number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfic = NULL;   /* communitive number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *ttv = NULL;     /* view angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tts[22];         /* sun angle table */
    int32 indts[22];       /* index for sun angle table */
    int iaots;             /* index for AOTs */

    /* Atmospheric correction coefficient variables (semi-empirical approach) */
    float tgo_arr[NREFL_BANDS];     /* per-band other gaseous transmittance */
    float roatm_arr[NREFL_BANDS][NAOT_VALS];  /* per band AOT vals for roatm */
    float ttatmg_arr[NREFL_BANDS][NAOT_VALS]; /* per band AOT vals for ttatmg */
    float satm_arr[NREFL_BANDS][NAOT_VALS];   /* per band AOT vals for satm */
    float roatm_coef[NREFL_BANDS][NCOEF];  /* per band poly coeffs for roatm */
    float ttatmg_coef[NREFL_BANDS][NCOEF]; /* per band poly coeffs for ttatmg */
    float satm_coef[NREFL_BANDS][NCOEF];   /* per band poly coeffs for satm */
    float normext_p0a3_arr[NREFL_BANDS];   /* per band normext[iband][0][3] */
    int roatm_iaMax[NREFL_BANDS];          /* ??? */
    int ia;                                /* looping variable for AOTs */
    int iaMaxTemp;                         /* max temp for current AOT level */

    /* Auxiliary file variables */
    int16 *dem = NULL;        /* CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 *andwi = NULL;      /* avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *sndwi = NULL;      /* standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob1 = NULL;    /* mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob2 = NULL;    /* mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob7 = NULL;    /* mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob1 = NULL;   /* intercept band1 ratio,
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *intratiob2 = NULL;   /* intercept band2 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *intratiob7 = NULL;   /* intercept band7 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *slpratiob1 = NULL;   /* slope band1 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *slpratiob2 = NULL;   /* slope band2 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *slpratiob7 = NULL;   /* slope band7 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    uint16 *wv = NULL;       /* water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 *oz = NULL;        /* ozone values [CMG_NBLAT x CMG_NBLON] */
    float raot550nm;    /* nearest input value of AOT */
    float uoz;          /* total column ozone */
    float uwv;          /* total column water vapor (precipital water vapor) */
    float pres;         /* surface pressure */
    float rb1;          /* band ratio 1 (unscaled) */
    float rb2;          /* band ratio 2 (unscaled) */
    float slpr11, slpr12, slpr21, slpr22;  /* band ratio slope at line,samp;
                           line, samp+1; line+1, samp; and line+1, samp+1 */
    float intr11, intr12, intr21, intr22;  /* band ratio intercept at line,samp;
                           line, samp+1; line+1, samp; and line+1, samp+1 */
    float slprb1, slprb2, slprb7;  /* interpolated band ratio slope values for
                                      band ratios 1, 2, 7 */
    float intrb1, intrb2, intrb7;  /* interpolated band ratio intercept values
                                      for band ratios 1, 2, 7 */
    int ratio_pix11;  /* pixel location for ratio products [lcmg][scmg] */
    int ratio_pix12;  /* pixel location for ratio products [lcmg][scmg+1] */
    int ratio_pix21;  /* pixel location for ratio products [lcmg+1][scmg] */
    int ratio_pix22;  /* pixel location for ratio products [lcmg+1][scmg+1] */
    int cmg_pix11;    /* pixel location for CMG/DEM products [lcmg][scmg] */
    int cmg_pix12;    /* pixel location for CMG/DEM products [lcmg][scmg+1] */
    int cmg_pix21;    /* pixel location for CMG/DEM products [lcmg+1][scmg] */
    int cmg_pix22;    /* pixel location for CMG/DEM products [lcmg+1][scmg+1] */

    /* Variables for finding the eps that minimizes the residual */
    double xa, xb, xc, xd, xe, xf;  /* coefficients */
    double coefa, coefb;            /* coefficients */
    float epsmin;                   /* eps which minimizes the residual */

    /* Output file info */
    time_t mytime;               /* timing variable */
    Output_t *sr_output = NULL;  /* output structure and metadata for the SR
                                    product */
    Envi_header_t envi_hdr;      /* output ENVI header information */
    char envi_file[STR_SIZE];    /* ENVI filename */
    char *cptr = NULL;           /* pointer to the file extension */

    /* Table constants */
    float aot550nm[NAOT_VALS] =  /* AOT look-up table */
        {0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00, 1.20,
         1.40, 1.60, 1.80, 2.00, 2.30, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00};
    float tpres[NPRES_VALS] =    /* surface pressure table */
        {1050.0, 1013.0, 900.0, 800.0, 700.0, 600.0, 500.0};

    /* Atmospheric correction variables */
    /* Look up table for atmospheric and geometric quantities.  Taurary comes
       from tauray-ldcm/msi.ASC and the oz, wv, og variables come from
       gascoef-modis/msi.ASC. */
    float tauray[NSRL_BANDS] =  /* molecular optical thickness coefficients --
                                   produced by running 6S */
        {0.23638, 0.16933, 0.09070, 0.04827, 0.01563, 0.00129, 0.00037,
         0.07984};
    double oztransa[NSRL_BANDS] =   /* ozone transmission coeff */
        {-0.00255649, -0.0177861, -0.0969872, -0.0611428, 0.0001, 0.0001,
          0.0001, -0.0834061};
    double wvtransa[NSRL_BANDS] =   /* water vapor transmission coeff */
        {2.29849e-27, 2.29849e-27, 0.00194772, 0.00404159, 0.000729136,
         0.00067324, 0.0177533, 0.00279738};
    double wvtransb[NSRL_BANDS] =   /* water vapor transmission coeff */
        {0.999742, 0.999742, 0.775024, 0.774482, 0.893085, 0.939669, 0.65094,
         0.759952};
    double ogtransa1[NSRL_BANDS] =  /* other gases transmission coeff */
        {4.91586e-20, 4.91586e-20, 4.91586e-20, 1.04801e-05, 1.35216e-05,
         0.0205425, 0.0256526, 0.000214329};
    double ogtransb0[NSRL_BANDS] =  /* other gases transmission coeff */
        {0.000197019, 0.000197019, 0.000197019, 0.640215, -0.195998, 0.326577,
         0.243961, 0.396322};
    double ogtransb1[NSRL_BANDS] =  /* other gases transmission coeff */
        {9.57011e-16, 9.57011e-16, 9.57011e-16, -0.348785, 0.275239, 0.0117192,
         0.0616101, 0.04728};

#ifdef WRITE_TAERO
    FILE *aero_fptr=NULL;   /* file pointer for aerosol files */
#endif

    /* Start processing */
    mytime = time(NULL);
    printf ("Start surface reflectance corrections: %s", ctime(&mytime));

    /* Allocate memory for the many arrays needed to do the surface reflectance
       computations */
    npixels = nlines * nsamps;
    retval = landsat_memory_allocation_sr (nlines, nsamps, &aerob1, &aerob2,
        &aerob4, &aerob5, &aerob7, &ipflag, &twvi, &tozi, &tp, &taero, &teps,
        &dem, &andwi, &sndwi, &ratiob1, &ratiob2, &ratiob7, &intratiob1,
        &intratiob2, &intratiob7, &slpratiob1, &slpratiob2, &slpratiob7, &wv,
        &oz, &rolutt, &transt, &sphalbt, &normext, &tsmax, &tsmin, &nbfic,
        &nbfi, &ttv);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error allocating memory for the data arrays needed "
            "for surface reflectance calculations.");
        error_handler (false, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Initialize the geolocation space applications */
    if (!get_geoloc_info (xml_metadata, &space_def))
    {
        sprintf (errmsg, "Getting the space definition from the XML file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    space = setup_mapping (&space_def);
    if (space == NULL)
    {
        sprintf (errmsg, "Setting up the geolocation mapping");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Initialize the look up tables and atmospheric correction variables.
       view zenith initialized to 0.0 (xtv)
       azimuthal difference between sun and obs angle initialize to 0.0 (xfi)
       surface pressure is initialized to the pressure at the center of the
           scene (using the DEM) (pres)
       water vapor is initialized to the value at the center of the scene (uwv)
       ozone is initialized to the value at the center of the scene (uoz) */
    retval = init_sr_refl (nlines, nsamps, input, &space_def, space, anglehdf,
        intrefnm, transmnm, spheranm, cmgdemnm, rationm, auxnm, aux2nm, aux_src, &eps,
        &iaots, &xtv_center, &xmuv_center, &xfi_center, &cosxfi_center,
        &raot550nm, &pres, &uoz, &uwv, &xtsstep, &xtsmin, &xtvstep, &xtvmin,
        tsmax, tsmin, tts, ttv, indts, rolutt, transt, sphalbt, normext, nbfic,
        nbfi, dem, andwi, sndwi, ratiob1, ratiob2, ratiob7, intratiob1,
        intratiob2, intratiob7, slpratiob1, slpratiob2, slpratiob7, wv, oz);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error initializing the lookup tables and "
            "atmospheric correction variables.");
        error_handler (false, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through all the reflectance bands and perform atmospheric
       corrections based on climatology */
    mytime = time(NULL);
    printf ("Performing atmospheric corrections for each Landsat reflectance "
        "band ... %s", ctime(&mytime));

    /* rotoa is not defined for the atmcorlamb2 call, which is ok, but the
       roslamb value is not valid upon output. Just set it to 0.0 to be
       consistent. */
    rotoa = 0.0;
    raot550nm = aot550nm[1];
    eps = 2.5;
    for (ib = 0; ib <= SRL_BAND7; ib++)
    {
        /* Get the parameters for the atmospheric correction */
        retval = atmcorlamb2 (input->meta.sat, xts_center, xtv_center,
            xmus_center, xmuv_center, xfi_center, cosxfi_center, raot550nm, ib,
            pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
            xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
            ttv, uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
            wvtransb, oztransa, rotoa, &roslamb, &tgo, &roatm, &ttatmg, &satm,
            &xrorayp, eps);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Performing lambertian atmospheric correction "
                "type 2.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Save these band-related parameters for later */
        btgo[ib] = tgo;
        broatm[ib] = roatm;
        bttatmg[ib] = ttatmg;
        bsatm[ib] = satm;
        tgo_x_roatm = tgo * roatm;
        tgo_x_ttatmg = tgo * ttatmg;

        /* Perform atmospheric corrections for bands 1-7 */
#ifdef _OPENMP
        #pragma omp parallel for private (i, roslamb)
#endif
        for (i = 0; i < npixels; i++)
        {
            /* Skip fill pixels, which have already been marked in the
               TOA calculations. */
            if (level1_qa_is_fill(qaband[i]))
            {
                if (ib == DNL_BAND1)
                    /* Initialize the fill flag, only need to do for band 1 */
                    ipflag[i] = (1 << IPFLAG_FILL);
                continue;
            }

            /* Store the unscaled TOA reflectance values for later use before
               completing atmospheric corrections */
            if (ib == DNL_BAND1)
                aerob1[i] = sband[ib][i];
            else if (ib == DNL_BAND2)
                aerob2[i] = sband[ib][i];
            else if (ib == DNL_BAND4)
                aerob4[i] = sband[ib][i];
            else if (ib == DNL_BAND5)
                aerob5[i] = sband[ib][i];
            else if (ib == DNL_BAND7)
                aerob7[i] = sband[ib][i];

            /* Apply the atmospheric corrections (ignoring the Rayleigh
               scattering component and water vapor), and store the
               unscaled value for further corrections.  (NOTE: the full
               computations are in atmcorlamb2) */
            roslamb = sband[ib][i] - tgo_x_roatm;
            roslamb /= tgo_x_ttatmg + satm * roslamb;

            /* Save the unscaled surface reflectance value */
            if (roslamb < MIN_VALID_REFL)
                sband[ib][i] = MIN_VALID_REFL;
            else if (roslamb > MAX_VALID_REFL)
                sband[ib][i] = MAX_VALID_REFL;
            else
                sband[ib][i] = roslamb;
        }  /* end for i */
    }  /* for ib */

    /* Start the retrieval of atmospheric correction parameters for each band */
    mytime = time(NULL);
    printf ("Starting retrieval of atmospheric correction parameters ... %s",
        ctime(&mytime));

    /* Get the coefficients for the semi-empirical atmospheric correction */
    if (!use_orig_aero)
    {
        mytime = time(NULL);
        printf ("Obtaining the coefficients for the semi-empirical approach "
            "... %s", ctime(&mytime));
        for (ib = 0; ib <= SRL_BAND7; ib++)
        {
            normext_p0a3_arr[ib] = normext[ib * NPRES_VALS * NAOT_VALS + 0 + 3];
                /* normext[ib][0][3]; */
            rotoa = 0.0;
            eps = 2.5;
            for (ia = 0; ia < NAOT_VALS; ia++)
            {
                raot550nm = aot550nm[ia];
                retval = atmcorlamb2 (input->meta.sat, xts_center, xtv_center,
                    xmus_center, xmuv_center, xfi_center, cosxfi_center,
                    raot550nm, ib, pres, tpres, aot550nm, rolutt, transt,
                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                    ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                    oztransa, rotoa, &roslamb, &tgo, &roatm, &ttatmg, &satm,
                    &xrorayp, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2 for band %d.", ib);
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
    
                /* Store the AOT-related variables for use in the atmospheric
                   corrections */
                roatm_arr[ib][ia] = roatm;
                ttatmg_arr[ib][ia] = ttatmg;
                satm_arr[ib][ia] = satm;
            }
    
            /* Store the band-related variables for use in the atmospheric
               corrections. tgo and xrorayp are the same for each AOT, so just
               save the last set for this band. */
            tgo_arr[ib] = tgo;
        }

        /* Setup the 3rd order polynomial coefficients for the semi-empirical
           approach in the aerosol inversion */
        for (ib = 0; ib <= SRL_BAND7; ib++)
        {
            /* Determine the maximum AOT index */
            iaMaxTemp = 1;
            for (ia = 1; ia < NAOT_VALS; ia++)
            {
                if (ia == NAOT_VALS-1)
                    iaMaxTemp = NAOT_VALS-1;
    
                if ((roatm_arr[ib][ia] - roatm_arr[ib][ia-1]) > ESPA_EPSILON)
                    continue;
                else
                {
                    iaMaxTemp = ia-1;
                    break;
                }
            }
    
            /* Get the polynomial coefficients for roatm */
            roatm_iaMax[ib] = iaMaxTemp;
            get_3rd_order_poly_coeff (aot550nm, roatm_arr[ib], iaMaxTemp,
                roatm_coef[ib]);
    
            /* Get the polynomial coefficients for ttatmg */
            get_3rd_order_poly_coeff (aot550nm, ttatmg_arr[ib], NAOT_VALS,
                ttatmg_coef[ib]);
    
            /* Get the polynomial coefficients for satm */
            get_3rd_order_poly_coeff (aot550nm, satm_arr[ib], NAOT_VALS,
                satm_coef[ib]);
        }
    }  /* if !use_orig_aero */

    /* If using the original aerosol approach we need some auxiliary data to
       be interpolated for every pixel so it's available for the final aerosol
       correction */
    if (use_orig_aero)
    {
        mytime = time(NULL);
        printf ("Interpolating the auxiliary data ... %s", ctime(&mytime));

#ifdef _OPENMP
    #pragma omp parallel for private (i, j, curr_pix, img, geo, lat, lon, xcmg, ycmg, lcmg, scmg, lcmg1, scmg1, u, v, one_minus_u, one_minus_v, one_minus_u_x_one_minus_v, one_minus_u_x_v, u_x_one_minus_v, u_x_v, cmg_pix11, cmg_pix12, cmg_pix21, cmg_pix22, wv11, wv12, wv21, wv22, uoz11, uoz12, uoz21, uoz22, pres11, pres12, pres21, pres22)
#endif
        for (i = 0; i < nlines; i++)
        {
            curr_pix = i * nsamps;
            for (j = 0; j < nsamps; j++, curr_pix++)
            {
                /* If this pixel is fill, do not process */
                if (qaband[curr_pix] == 1)
                {
                    ipflag[curr_pix] |= (1 << IPFLAG_FILL);
                    continue;
                }
    
                /* Get the lat/long for the current pixel */
                img.l = i - 0.5;
                img.s = j + 0.5;
                img.is_fill = false;
                if (!from_space (space, &img, &geo))
                {
                    sprintf (errmsg, "Mapping line/sample (%d, %d) to "
                        "geolocation coords", i, j);
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
                lat = geo.lat * RAD2DEG;
                lon = geo.lon * RAD2DEG;
    
                /* Use that lat/long to determine the line/sample in the
                   CMG-related lookup tables, using the center of the UL
                   pixel. Note, we are basically making sure the line/sample
                   combination falls within -90, 90 and -180, 180 global climate
                   data boundaries.  However, the source code below uses lcmg+1
                   and scmg+1, which for some scenes may wrap around the
                   dateline or the poles.  Thus we need to wrap the CMG data
                   around to the beginning of the array. */
                /* Each CMG pixel is 0.05 x 0.05 degrees.  Use the center of the
                   pixel for each calculation.  Negative latitude values should
                   be the largest line values in the CMG grid.  Negative
                   longitude values should be the smallest sample values in the
                   CMG grid. */
                /* The line/sample calculation from the x/ycmg values are not
                   rounded.  The interpolation of the value using line+1 and
                   sample+1 are based on the truncated numbers, therefore
                   rounding up is not appropriate. */
                ycmg = (89.975 - lat) * 20.0;   /* vs / 0.05 */
                xcmg = (179.975 + lon) * 20.0;  /* vs / 0.05 */
                lcmg = (int) ycmg;
                scmg = (int) xcmg;

                /* Handle the edges of the lat/long values in the CMG grid */
                if (lcmg < 0)
                    lcmg = 0;
                else if (lcmg >= CMG_NBLAT)
                    lcmg = CMG_NBLAT;
    
                if (scmg < 0)
                    scmg = 0;
                else if (scmg >= CMG_NBLON)
                    scmg = CMG_NBLON;
    
                /* If the current CMG pixel is at the edge of the CMG array,
                   then allow the next pixel for interpolation to wrap around
                   the array */
                if (scmg >= CMG_NBLON-1)  /* 180 degrees so wrap around */
                    scmg1 = 0;
                else
                    scmg1 = scmg + 1;
    
                if (lcmg >= CMG_NBLAT-1)  /* -90 degrees so wrap around */
                    lcmg1 = 0;
                else
                    lcmg1 = lcmg + 1;
    
                /* Determine the four CMG pixels to be used for the current
                   Landsat pixel */
                cmg_pix11 = lcmg * CMG_NBLON + scmg;
                cmg_pix12 = lcmg * CMG_NBLON + scmg1;
                cmg_pix21 = lcmg1 * CMG_NBLON + scmg;
                cmg_pix22 = lcmg1 * CMG_NBLON + scmg1;
    
                /* Get the water vapor pixels. If the water vapor value is
                   fill, then use a default value. */
                wv11 = wv[cmg_pix11];
                if (wv11 == WV_FILL[aux_src])
                    wv11 = WV_DEFAULT_DN[aux_src];
                wv12 = wv[cmg_pix12];
                if (wv12 == WV_FILL[aux_src])
                    wv12 = WV_DEFAULT_DN[aux_src];
                wv21 = wv[cmg_pix21];
                if (wv21 == WV_FILL[aux_src])
                    wv21 = WV_DEFAULT_DN[aux_src];
                wv22 = wv[cmg_pix22];
                if (wv22 == WV_FILL[aux_src])
                    wv22 = WV_DEFAULT_DN[aux_src];

                /* Get the ozone pixels. If the ozone value is fill, then
                   use a default value. */
                uoz11 = oz[cmg_pix11];
                if (uoz11 == OZ_FILL[aux_src])
                    uoz11 = OZ_DEFAULT_DN[aux_src];
    
                uoz12 = oz[cmg_pix12];
                if (uoz12 == OZ_FILL[aux_src])
                    uoz12 = OZ_DEFAULT_DN[aux_src];
    
                uoz21 = oz[cmg_pix21];
                if (uoz21 == OZ_FILL[aux_src])
                    uoz21 = OZ_DEFAULT_DN[aux_src];
    
                uoz22 = oz[cmg_pix22];
                if (uoz22 == OZ_FILL[aux_src])
                    uoz22 = OZ_DEFAULT_DN[aux_src];
    
                /* Get the surface pressure from the global DEM.  Set to 1013.0
                   (sea level) if the DEM is fill (= -9999), which is likely
                   ocean. The dimensions on the DEM array is the same as that
                   of the CMG arrays. Use the current pixel locations already
                   calculated. */
                if (dem[cmg_pix11] != -9999)
                    pres11 = 1013.0 * exp (-dem[cmg_pix11] * ONE_DIV_8500);
                else
                    pres11 = 1013.0;
    
                if (dem[cmg_pix12] != -9999)
                    pres12 = 1013.0 * exp (-dem[cmg_pix12] * ONE_DIV_8500);
                else
                    pres12 = 1013.0;
    
                if (dem[cmg_pix21] != -9999)
                    pres21 = 1013.0 * exp (-dem[cmg_pix21] * ONE_DIV_8500);
                else
                    pres21 = 1013.0;
    
                if (dem[cmg_pix22] != -9999)
                    pres22 = 1013.0 * exp (-dem[cmg_pix22] * ONE_DIV_8500);
                else
                    pres22 = 1013.0;
    
                /* Determine the fractional difference between the integer
                   location and floating point pixel location to be used for
                   interpolation */
                u = (ycmg - lcmg);
                v = (xcmg - scmg);
                one_minus_u = 1.0 - u;
                one_minus_v = 1.0 - v;
                one_minus_u_x_one_minus_v = one_minus_u * one_minus_v;
                one_minus_u_x_v = one_minus_u * v;
                u_x_one_minus_v = u * one_minus_v;
                u_x_v = u * v;
    
                /* Interpolate water vapor, and unscale */
                twvi[curr_pix] = wv11 * one_minus_u_x_one_minus_v +
                                 wv12 * one_minus_u_x_v +
                                 wv21 * u_x_one_minus_v +
                                 wv22 * u_x_v;
                twvi[curr_pix] = twvi[curr_pix] * 0.01;   /* vs / 100 */

                /* Interpolate ozone, and unscale */
                tozi[curr_pix] = uoz11 * one_minus_u_x_one_minus_v +
                                 uoz12 * one_minus_u_x_v +
                                 uoz21 * u_x_one_minus_v +
                                 uoz22 * u_x_v;
                tozi[curr_pix] = tozi[curr_pix] * 0.0025;   /* vs / 400 */
    
    
                /* Interpolate surface pressure */
                tp[curr_pix] = pres11 * one_minus_u_x_one_minus_v +
                               pres12 * one_minus_u_x_v +
                               pres21 * u_x_one_minus_v +
                               pres22 * u_x_v;
            }  /* end for j */
        }  /* end for i */
    }  /* if use_orig_aero */

    /* Initialize and compute some EPS values */
    eps1 = LOW_EPS;
    eps2 = MOD_EPS;
    eps3 = HIGH_EPS;
    xa = (eps1 * eps1) - (eps3 * eps3);
    xd = (int) ((eps2 * eps2) - (eps3 * eps3));
    xb = eps1 - eps3;
    xe = eps2 - eps3;

    /* Start the aerosol inversion */
    mytime = time(NULL);
    printf ("Aerosol Inversion using %d x %d aerosol window ... %s",
        LAERO_WINDOW, LAERO_WINDOW, ctime(&mytime));
#ifdef _OPENMP
    #pragma omp parallel for private (i, j, center_line, center_samp, nearest_line, nearest_samp, curr_pix, center_pix, img, geo, lat, lon, xcmg, ycmg, lcmg, scmg, lcmg1, scmg1, u, v, one_minus_u, one_minus_v, one_minus_u_x_one_minus_v, one_minus_u_x_v, u_x_one_minus_v, u_x_v, ratio_pix11, ratio_pix12, ratio_pix21, ratio_pix22, rb1, rb2, slpr11, slpr12, slpr21, slpr22, intr11, intr12, intr21, intr22, slprb1, slprb2, slprb7, intrb1, intrb2, intrb7, xndwi, ndwi_th1, ndwi_th2, ib, xtv, xmuv, xts, xmus, xfi, cosxfi, iband, iband1, iaots, pres, uoz, uwv, retval, eps, residual, residual1, residual2, residual3, raot, sraot1, sraot3, xc, xf, epsmin, corf, rotoa, raot550nm, roslamb, tgo, roatm, ttatmg, satm, xrorayp, ros1, ros5, ros4, erelc, troatm)
#else
    tmp_percent = 0;
#endif
    for (i = LHALF_AERO_WINDOW; i < nlines; i += LAERO_WINDOW)
    {
#ifndef _OPENMP
        /* update status, but not if multi-threaded */
        curr_tmp_percent = 100 * i / nlines;
        if (curr_tmp_percent > tmp_percent)
        {
            tmp_percent = curr_tmp_percent;
            if (tmp_percent % 10 == 0)
            {
                printf ("%d%% ", tmp_percent);
                fflush (stdout);
            }
        }
#endif

        curr_pix = i * nsamps + LHALF_AERO_WINDOW;
        for (j = LHALF_AERO_WINDOW; j < nsamps;
             j += LAERO_WINDOW, curr_pix += LAERO_WINDOW)
        {
            /* Keep track of the center pixel for the current aerosol window;
               may need to return here if this is fill, cloudy or water */
            center_line = i;
            center_samp = j;
            center_pix = curr_pix;

            /* If this pixel is fill */
            if (level1_qa_is_fill (qaband[curr_pix]))
            {
                /* Look for other non-fill pixels in the window */
                if (find_closest_non_fill (qaband, nlines, nsamps, center_line,
                    center_samp, LHALF_AERO_WINDOW, &nearest_line,
                    &nearest_samp))
                {
                    /* Use the line/sample location of the non-fill pixel for
                       further processing of aerosols. However we will still
                       write to the center of the aerosol window for the
                       current window. */
                    i = nearest_line;
                    j = nearest_samp;
                    curr_pix = i * nsamps + j;
                }
                else
                {
                    /* No other non-fill pixels found.  Pixel is already
                       flagged as fill. Move to next aerosol window. */
                    continue;
                }
            }

            /* Get the lat/long for the current pixel (which may not be the
               center of the aerosol window), for the center of that pixel */
            img.l = i - 0.5;
            img.s = j + 0.5;
            img.is_fill = false;
            if (!from_space (space, &img, &geo))
            {
                sprintf (errmsg, "Mapping line/sample (%d, %d) to "
                    "geolocation coords", i, j);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
            lat = geo.lat * RAD2DEG;
            lon = geo.lon * RAD2DEG;

            /* Use that lat/long to determine the line/sample in the
               CMG-related lookup tables, using the center of the UL
               pixel. Note, we are basically making sure the line/sample
               combination falls within -90, 90 and -180, 180 global climate
               data boundaries.  However, the source code below uses lcmg+1
               and scmg+1, which for some scenes may wrap around the
               dateline or the poles.  Thus we need to wrap the CMG data
               around to the beginning of the array. */
            /* Each CMG pixel is 0.05 x 0.05 degrees.  Use the center of the
               pixel for each calculation.  Negative latitude values should
               be the largest line values in the CMG grid.  Negative
               longitude values should be the smallest sample values in the
               CMG grid. */
            /* The line/sample calculation from the x/ycmg values are not
               rounded.  The interpolation of the value using line+1 and
               sample+1 are based on the truncated numbers, therefore
               rounding up is not appropriate. */
            ycmg = (89.975 - lat) * 20.0;   /* vs / 0.05 */
            xcmg = (179.975 + lon) * 20.0;  /* vs / 0.05 */
            lcmg = (int) ycmg;
            scmg = (int) xcmg;

            /* Handle the edges of the lat/long values in the CMG grid */
            if (lcmg < 0)
                lcmg = 0;
            else if (lcmg >= CMG_NBLAT)
                lcmg = CMG_NBLAT - 1;

            if (scmg < 0)
                scmg = 0;
            else if (scmg >= CMG_NBLON)
                scmg = CMG_NBLON - 1;

            /* If the current CMG pixel is at the edge of the CMG array, then
               allow the next pixel for interpolation to wrap around the
               array */
            if (scmg >= CMG_NBLON-1)  /* 180 degrees so wrap around */
                scmg1 = 0;
            else
                scmg1 = scmg + 1;

            if (lcmg >= CMG_NBLAT-1)  /* -90 degrees, so set the next pixel
                                         to also use -90 */
                lcmg1 = lcmg;
            else
                lcmg1 = lcmg + 1;

            /* Determine the fractional difference between the integer location
               and floating point pixel location to be used for interpolation */
            u = (ycmg - lcmg);
            v = (xcmg - scmg);
            one_minus_u = 1.0 - u;
            one_minus_v = 1.0 - v;
            one_minus_u_x_one_minus_v = one_minus_u * one_minus_v;
            one_minus_u_x_v = one_minus_u * v;
            u_x_one_minus_v = u * one_minus_v;
            u_x_v = u * v;

            /* Determine the band ratios and slope/intercept */
            ratio_pix11 = lcmg * RATIO_NBLON + scmg;
            ratio_pix12 = lcmg * RATIO_NBLON + scmg1;
            ratio_pix21 = lcmg1 * RATIO_NBLON + scmg;
            ratio_pix22 = lcmg1 * RATIO_NBLON + scmg1;

            rb1 = ratiob1[ratio_pix11] * 0.001;  /* vs. / 1000. */
            rb2 = ratiob2[ratio_pix11] * 0.001;  /* vs. / 1000. */
            if (rb2 > 1.0 || rb1 > 1.0 || rb2 < 0.1 || rb1 < 0.1)
            {
                slpratiob1[ratio_pix11] = 0;
                slpratiob2[ratio_pix11] = 0;
                slpratiob7[ratio_pix11] = 0;
                intratiob1[ratio_pix11] = 550;
                intratiob2[ratio_pix11] = 600;
                intratiob7[ratio_pix11] = 2000;
            }
            else if (sndwi[ratio_pix11] < 200)
            {
                slpratiob1[ratio_pix11] = 0;
                slpratiob2[ratio_pix11] = 0;
                slpratiob7[ratio_pix11] = 0;
                intratiob1[ratio_pix11] = ratiob1[ratio_pix11];
                intratiob2[ratio_pix11] = ratiob2[ratio_pix11];
                intratiob7[ratio_pix11] = ratiob7[ratio_pix11];
            }

            rb1 = ratiob1[ratio_pix12] * 0.001;  /* vs. / 1000. */
            rb2 = ratiob2[ratio_pix12] * 0.001;  /* vs. / 1000. */
            if (rb2 > 1.0 || rb1 > 1.0 || rb2 < 0.1 || rb1 < 0.1)
            {
                slpratiob1[ratio_pix12] = 0;
                slpratiob2[ratio_pix12] = 0;
                slpratiob7[ratio_pix12] = 0;
                intratiob1[ratio_pix12] = 550;
                intratiob2[ratio_pix12] = 600;
                intratiob7[ratio_pix12] = 2000;
            }
            else if (sndwi[ratio_pix12] < 200)
            {
                slpratiob1[ratio_pix12] = 0;
                slpratiob2[ratio_pix12] = 0;
                slpratiob7[ratio_pix12] = 0;
                intratiob1[ratio_pix12] = ratiob1[ratio_pix12];
                intratiob2[ratio_pix12] = ratiob2[ratio_pix12];
                intratiob7[ratio_pix12] = ratiob7[ratio_pix12];
            }

            rb1 = ratiob1[ratio_pix21] * 0.001;  /* vs. / 1000. */
            rb2 = ratiob2[ratio_pix21] * 0.001;  /* vs. / 1000. */
            if (rb2 > 1.0 || rb1 > 1.0 || rb2 < 0.1 || rb1 < 0.1)
            {
                slpratiob1[ratio_pix21] = 0;
                slpratiob2[ratio_pix21] = 0;
                slpratiob7[ratio_pix21] = 0;
                intratiob1[ratio_pix21] = 550;
                intratiob2[ratio_pix21] = 600;
                intratiob7[ratio_pix21] = 2000;
            }
            else if (sndwi[ratio_pix21] < 200)
            {
                slpratiob1[ratio_pix21] = 0;
                slpratiob2[ratio_pix21] = 0;
                slpratiob7[ratio_pix21] = 0;
                intratiob1[ratio_pix21] = ratiob1[ratio_pix21];
                intratiob2[ratio_pix21] = ratiob2[ratio_pix21];
                intratiob7[ratio_pix21] = ratiob7[ratio_pix21];
            }

            rb1 = ratiob1[ratio_pix22] * 0.001;  /* vs. / 1000. */
            rb2 = ratiob2[ratio_pix22] * 0.001;  /* vs. / 1000. */
            if (rb2 > 1.0 || rb1 > 1.0 || rb2 < 0.1 || rb1 < 0.1)
            {
                slpratiob1[ratio_pix22] = 0;
                slpratiob2[ratio_pix22] = 0;
                slpratiob7[ratio_pix22] = 0;
                intratiob1[ratio_pix22] = 550;
                intratiob2[ratio_pix22] = 600;
                intratiob7[ratio_pix22] = 2000;
            }
            else if (sndwi[ratio_pix22] < 200)
            {
                slpratiob1[ratio_pix22] = 0;
                slpratiob2[ratio_pix22] = 0;
                slpratiob7[ratio_pix22] = 0;
                intratiob1[ratio_pix22] = ratiob1[ratio_pix22];
                intratiob2[ratio_pix22] = ratiob2[ratio_pix22];
                intratiob7[ratio_pix22] = ratiob7[ratio_pix22];
            }

            /* Interpolate the slope/intercept for each band, and unscale */
            slpr11 = slpratiob1[ratio_pix11] * 0.001;  /* vs / 1000 */
            intr11 = intratiob1[ratio_pix11] * 0.001;  /* vs / 1000 */
            slpr12 = slpratiob1[ratio_pix12] * 0.001;  /* vs / 1000 */
            intr12 = intratiob1[ratio_pix12] * 0.001;  /* vs / 1000 */
            slpr21 = slpratiob1[ratio_pix21] * 0.001;  /* vs / 1000 */
            intr21 = intratiob1[ratio_pix21] * 0.001;  /* vs / 1000 */
            slpr22 = slpratiob1[ratio_pix22] * 0.001;  /* vs / 1000 */
            intr22 = intratiob1[ratio_pix22] * 0.001;  /* vs / 1000 */
            slprb1 = slpr11 * one_minus_u_x_one_minus_v +
                     slpr12 * one_minus_u_x_v +
                     slpr21 * u_x_one_minus_v +
                     slpr22 * u_x_v;
            intrb1 = intr11 * one_minus_u_x_one_minus_v +
                     intr12 * one_minus_u_x_v +
                     intr21 * u_x_one_minus_v +
                     intr22 * u_x_v;

            slpr11 = slpratiob2[ratio_pix11] * 0.001;  /* vs / 1000 */
            intr11 = intratiob2[ratio_pix11] * 0.001;  /* vs / 1000 */
            slpr12 = slpratiob2[ratio_pix12] * 0.001;  /* vs / 1000 */
            intr12 = intratiob2[ratio_pix12] * 0.001;  /* vs / 1000 */
            slpr21 = slpratiob2[ratio_pix21] * 0.001;  /* vs / 1000 */
            intr21 = intratiob2[ratio_pix21] * 0.001;  /* vs / 1000 */
            slpr22 = slpratiob2[ratio_pix22] * 0.001;  /* vs / 1000 */
            intr22 = intratiob2[ratio_pix22] * 0.001;  /* vs / 1000 */
            slprb2 = slpr11 * one_minus_u_x_one_minus_v +
                     slpr12 * one_minus_u_x_v +
                     slpr21 * u_x_one_minus_v +
                     slpr22 * u_x_v;
            intrb2 = intr11 * one_minus_u_x_one_minus_v +
                     intr12 * one_minus_u_x_v +
                     intr21 * u_x_one_minus_v +
                     intr22 * u_x_v;

            slpr11 = slpratiob7[ratio_pix11] * 0.001;  /* vs / 1000 */
            intr11 = intratiob7[ratio_pix11] * 0.001;  /* vs / 1000 */
            slpr12 = slpratiob7[ratio_pix12] * 0.001;  /* vs / 1000 */
            intr12 = intratiob7[ratio_pix12] * 0.001;  /* vs / 1000 */
            slpr21 = slpratiob7[ratio_pix21] * 0.001;  /* vs / 1000 */
            intr21 = intratiob7[ratio_pix21] * 0.001;  /* vs / 1000 */
            slpr22 = slpratiob7[ratio_pix22] * 0.001;  /* vs / 1000 */
            intr22 = intratiob7[ratio_pix22] * 0.001;  /* vs / 1000 */
            slprb7 = slpr11 * one_minus_u_x_one_minus_v +
                     slpr12 * one_minus_u_x_v +
                     slpr21 * u_x_one_minus_v +
                     slpr22 * u_x_v;
            intrb7 = intr11 * one_minus_u_x_one_minus_v +
                     intr12 * one_minus_u_x_v +
                     intr21 * u_x_one_minus_v +
                     intr22 * u_x_v;

            /* Calculate NDWI variables for the band ratios */
            xndwi = ((double) sband[SRL_BAND5][curr_pix] -
                     (double) (sband[SRL_BAND7][curr_pix] * 0.5)) /
                    ((double) sband[SRL_BAND5][curr_pix] +
                     (double) (sband[SRL_BAND7][curr_pix] * 0.5));

            ndwi_th1 = (andwi[ratio_pix11] + 2.0 *
                        sndwi[ratio_pix11]) * 0.001;
            ndwi_th2 = (andwi[ratio_pix11] - 2.0 *
                        sndwi[ratio_pix11]) * 0.001;
            if (xndwi > ndwi_th1)
                xndwi = ndwi_th1;
            if (xndwi < ndwi_th2)
                xndwi = ndwi_th2;

            /* Initialize the band ratios */
            for (ib = 0; ib < NSR_BANDS; ib++)
            {
                erelc[ib] = -1.0;
                troatm[ib] = 0.0;
            }

            /* Compute the band ratio - coastal aerosol, blue, red, SWIR */
            erelc[DNL_BAND1] = (xndwi * slprb1 + intrb1);
            erelc[DNL_BAND2] = (xndwi * slprb2 + intrb2);
            erelc[DNL_BAND4] = 1.0;
            erelc[DNL_BAND7] = (xndwi * slprb7 + intrb7);

            /* Retrieve the TOA reflectance values for the current pixel */
            troatm[DNL_BAND1] = aerob1[curr_pix];
            troatm[DNL_BAND2] = aerob2[curr_pix];
            troatm[DNL_BAND4] = aerob4[curr_pix];
            troatm[DNL_BAND7] = aerob7[curr_pix];

            /* Determine the solar and view angles for the current pixel */
            if (use_orig_aero)
            {
                xtv = vza[curr_pix] * 0.01;
                xmuv = cos(xtv * DEG2RAD);
                xts = sza[curr_pix] * 0.01;
                xmus = cos(xts * DEG2RAD);
                xfi = saa[curr_pix] * 0.01 - vaa[curr_pix] * 0.01 ;
                cosxfi = cos(xfi * DEG2RAD);
            }

            /* Retrieve the aerosol information for low eps 1.0 */
            iband1 = DNL_BAND4;   /* red band */
            iaots = 0;
            if (use_orig_aero)
            {
                pres = tp[curr_pix];
                uoz = tozi[curr_pix];
                uwv = twvi[curr_pix];
                retval = subaeroret (input->meta.sat, false, iband1, xts, xtv,
                    xmus, xmuv, xfi, cosxfi, pres, uoz, uwv, erelc, troatm,
                    tpres, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
                    sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
                    ttv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                    wvtransb, oztransa, &raot, &residual, &iaots, eps1);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing aerosol retrieval.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
            }
            else
                subaeroret_new (input->meta.sat, false, iband1, erelc, troatm,
                    tgo_arr, roatm_iaMax, roatm_coef, ttatmg_coef, satm_coef,
                    normext_p0a3_arr, &raot, &residual, &iaots, eps1);

            /* Save the data */
            residual1 = residual;
            sraot1 = raot;

            /* Retrieve the aerosol information for moderate eps 1.75 */
            if (use_orig_aero)
            {
                retval = subaeroret (input->meta.sat, false, iband1, xts, xtv,
                    xmus, xmuv, xfi, cosxfi, pres, uoz, uwv, erelc, troatm,
                    tpres, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
                    sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
                    ttv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                    wvtransb, oztransa, &raot, &residual, &iaots, eps2);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing aerosol retrieval.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
            }
            else
                subaeroret_new (input->meta.sat, false, iband1, erelc, troatm,
                    tgo_arr, roatm_iaMax, roatm_coef, ttatmg_coef, satm_coef,
                    normext_p0a3_arr, &raot, &residual, &iaots, eps2);

            /* Save the data */
            residual2 = residual;

            /* Retrieve the aerosol information for high eps 2.5 */
            if (use_orig_aero)
            {
                retval = subaeroret (input->meta.sat, false, iband1, xts, xtv,
                    xmus, xmuv, xfi, cosxfi, pres, uoz, uwv, erelc, troatm,
                    tpres, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
                    sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
                    ttv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                    wvtransb, oztransa, &raot, &residual, &iaots, eps3);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing aerosol retrieval.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
            }
            else
                subaeroret_new (input->meta.sat, false, iband1, erelc, troatm,
                    tgo_arr, roatm_iaMax, roatm_coef, ttatmg_coef, satm_coef,
                    normext_p0a3_arr, &raot, &residual, &iaots, eps3);

            /* Save the data */
            residual3 = residual;
            sraot3 = raot;

            /* Find the eps (angstrom coefficient for AOT) that minimizes the
               residual */
            xc = residual1 - residual3;
            xf = residual2 - residual3;
            coefa = (xc*xe - xb*xf) / (xa*xe - xb*xd);
            coefb = (xa*xf - xc*xd) / (xa*xe - xb*xd);
            epsmin = -coefb / (2.0 * coefa);

            eps = epsmin;
            if (epsmin >= LOW_EPS && epsmin <= HIGH_EPS)
            {
                if (use_orig_aero)
                {
                    retval = subaeroret (input->meta.sat, false, iband1, xts,
                        xtv, xmus, xmuv, xfi, cosxfi, pres, uoz, uwv, erelc,
                        troatm, tpres, rolutt, transt, xtsstep, xtsmin,
                        xtvstep, xtvmin, sphalbt, normext, tsmax, tsmin, nbfic,
                        nbfi, tts, indts, ttv, tauray, ogtransa1, ogtransb0,
                        ogtransb1, wvtransa, wvtransb, oztransa, &raot,
                        &residual, &iaots, epsmin);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Performing aerosol retrieval.");
                        error_handler (true, FUNC_NAME, errmsg);
                        exit (ERROR);
                    }
                }
                else
                    subaeroret_new (input->meta.sat, false, iband1, erelc,
                        troatm, tgo_arr, roatm_iaMax, roatm_coef, ttatmg_coef,
                        satm_coef, normext_p0a3_arr, &raot, &residual, &iaots,
                        epsmin);
            }
            else if (epsmin <= LOW_EPS)
            {
                eps = eps1;
                residual = residual1;
                raot = sraot1;
            }
            else if (epsmin >= HIGH_EPS)
            {
                eps = eps3;
                residual = residual3;
                raot = sraot3;
            }

            teps[center_pix] = eps;
            taero[center_pix] = raot;
            if (use_orig_aero)
                corf = raot / xmus;
            else
                corf = raot / xmus_center;

            /* Check the model residual.  Corf represents aerosol impact.
               Test the quality of the aerosol inversion. */
            if (residual < (0.015 + 0.005 * corf + 0.10 * troatm[DNL_BAND7]))
            {
                /* Test if NIR band 5 makes sense */
                iband = DNL_BAND5;
                rotoa = aerob5[curr_pix];
                raot550nm = raot;
                if (use_orig_aero)
                {
                    retval = atmcorlamb2 (input->meta.sat, xts, xtv, xmus,
                        xmuv, xfi, cosxfi, raot550nm, iband, pres, tpres,
                        aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
                        xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi,
                        tts, indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                        ogtransb1, wvtransa, wvtransb, oztransa, rotoa,
                        &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, eps);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Performing lambertian "
                            "atmospheric correction type 2.");
                        error_handler (true, FUNC_NAME, errmsg);
                        exit (ERROR);
                    }
                }
                else
                    atmcorlamb2_new (input->meta.sat, tgo_arr[iband],
                        aot550nm[roatm_iaMax[iband]], &roatm_coef[iband][0],
                        &ttatmg_coef[iband][0], &satm_coef[iband][0], raot550nm,
                        iband, normext_p0a3_arr[iband], rotoa, &roslamb, eps);
                ros5 = roslamb;

                /* Test if red band 4 makes sense */
                iband = DNL_BAND4;
                rotoa = aerob4[curr_pix];
                raot550nm = raot;
                if (use_orig_aero)
                {
                    retval = atmcorlamb2 (input->meta.sat, xts, xtv, xmus,
                        xmuv, xfi, cosxfi, raot550nm, iband, pres, tpres,
                        aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
                        xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi,
                        tts, indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                        ogtransb1, wvtransa, wvtransb, oztransa, rotoa,
                        &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, eps);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Performing lambertian "
                            "atmospheric correction type 2.");
                        error_handler (true, FUNC_NAME, errmsg);
                        exit (ERROR);
                    }
                }
                else
                    atmcorlamb2_new (input->meta.sat, tgo_arr[iband],
                        aot550nm[roatm_iaMax[iband]], &roatm_coef[iband][0],
                        &ttatmg_coef[iband][0], &satm_coef[iband][0], raot550nm,
                        iband, normext_p0a3_arr[iband], rotoa, &roslamb, eps);
                ros4 = roslamb;

                /* Use the NDVI to validate the reflectance values or flag
                   as water */
                if ((ros5 > 0.1) && ((ros5 - ros4) / (ros5 + ros4) > 0))
                {
                    /* Clear pixel with valid aerosol retrieval */
                    ipflag[center_pix] |= (1 << IPFLAG_CLEAR);
                }
                else
                {
                    /* Flag as water */
                    ipflag[center_pix] |= (1 << IPFLAG_WATER);
                }
            }
            else
            {
                /* Flag as water */
                ipflag[center_pix] |= (1 << IPFLAG_WATER);
            }

            /* Retest any water pixels to verify they are water and obtain
               their aerosol */
            if (lasrc_qa_is_water(ipflag[center_pix]))
            {
                /* Initialize the band ratios */
                for (ib = 0; ib < NSR_BANDS; ib++)
                    erelc[ib] = -1.0;
                troatm[DNL_BAND1] = aerob1[curr_pix];
                troatm[DNL_BAND4] = aerob4[curr_pix];
                troatm[DNL_BAND5] = aerob5[curr_pix];
                troatm[DNL_BAND7] = aerob7[curr_pix];

                /* Set the band ratio - coastal aerosol, red, NIR, SWIR */
                erelc[DNL_BAND1] = 1.0;
                erelc[DNL_BAND4] = 1.0;
                erelc[DNL_BAND5] = 1.0;
                erelc[DNL_BAND7] = 1.0;

                /* Retrieve the water aerosol information for eps 1.5 */
                eps = WATER_EPS;
                iaots = 0;
                if (use_orig_aero)
                {
                    retval = subaeroret (input->meta.sat, true, iband1, xts,
                        xtv, xmus, xmuv, xfi, cosxfi, pres, uoz, uwv, erelc,
                        troatm, tpres, rolutt, transt, xtsstep, xtsmin,
                        xtvstep, xtvmin, sphalbt, normext, tsmax, tsmin, nbfic,
                        nbfi, tts, indts, ttv, tauray, ogtransa1, ogtransb0,
                        ogtransb1, wvtransa, wvtransb, oztransa, &raot,
                        &residual, &iaots, eps);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Performing aerosol retrieval.");
                        error_handler (true, FUNC_NAME, errmsg);
                        exit (ERROR);
                    }
                }
                else
                    subaeroret_new (input->meta.sat, true, iband1, erelc,
                        troatm, tgo_arr, roatm_iaMax, roatm_coef, ttatmg_coef,
                        satm_coef, normext_p0a3_arr, &raot, &residual, &iaots,
                        eps);
                teps[center_pix] = eps;
                taero[center_pix] = raot;
                if (use_orig_aero)
                    corf = raot / xmus;
                else
                    corf = raot / xmus_center;

                /* Test band 1 reflectance to eliminate negative */
                iband = DNL_BAND1;
                rotoa = aerob1[curr_pix];
                raot550nm = raot;
                if (use_orig_aero)
                {
                    retval = atmcorlamb2 (input->meta.sat, xts, xtv, xmus,
                        xmuv, xfi, cosxfi, raot550nm, iband, pres, tpres,
                        aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
                        xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi,
                        tts, indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                        ogtransb1, wvtransa, wvtransb, oztransa, rotoa,
                        &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, eps);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Performing lambertian "
                            "atmospheric correction type 2.");
                        error_handler (true, FUNC_NAME, errmsg);
                        exit (ERROR);
                    }
                }
                else
                    atmcorlamb2_new (input->meta.sat, tgo_arr[iband],
                        aot550nm[roatm_iaMax[iband]], &roatm_coef[iband][0],
                        &ttatmg_coef[iband][0], &satm_coef[iband][0], raot550nm,
                        iband, normext_p0a3_arr[iband], rotoa, &roslamb, eps);
                ros1 = roslamb;

                if (residual > (0.010 + 0.005 * corf) || ros1 < 0)
                {
                    /* Not a valid water pixel (possibly urban). Clear all
                       the QA bits, and leave the IPFLAG_CLEAR bit off to
                       indicate the aerosol retrieval was not valid. */
                    ipflag[center_pix] = 0;  /* IPFLAG_CLEAR bit is 0 */
                }
                else
                {
                    /* Valid water pixel. Set the clear aerosol retrieval bit
                       and turn on the water bit. */
                    ipflag[center_pix] = (1 << IPFLAG_CLEAR);
                    ipflag[center_pix] |= (1 << IPFLAG_WATER);
                }
            }  /* if water pixel */

            /* Reset the looping variables to the center of the aerosol window
               versus the actual non-fill/non-cloud pixel that was processed
               so that we get the correct center for the next aerosol window */
            i = center_line;
            j = center_samp;
            curr_pix = center_pix;
        }  /* end for j */
    }  /* end for i */

#ifndef _OPENMP
    /* update status */
    printf ("100%%\n");
    fflush (stdout);
#endif

    /* Done with the aerob* arrays */
    free (aerob1);  aerob1 = NULL;
    free (aerob2);  aerob2 = NULL;
    free (aerob4);  aerob4 = NULL;
    free (aerob5);  aerob5 = NULL;
    free (aerob7);  aerob7 = NULL;

    /* Done with the ratiob* arrays */
    free (andwi);  andwi = NULL;
    free (sndwi);  sndwi = NULL;
    free (ratiob1);  ratiob1 = NULL;
    free (ratiob2);  ratiob2 = NULL;
    free (ratiob7);  ratiob7 = NULL;
    free (intratiob1);  intratiob1 = NULL;
    free (intratiob2);  intratiob2 = NULL;
    free (intratiob7);  intratiob7 = NULL;
    free (slpratiob1);  slpratiob1 = NULL;
    free (slpratiob2);  slpratiob2 = NULL;
    free (slpratiob7);  slpratiob7 = NULL;

    /* Done with the DEM, water vapor, and ozone arrays */
    free (dem);  dem = NULL;
    free (wv);  wv = NULL;
    free (oz);  oz = NULL;

#ifdef WRITE_TAERO
    /* Write the ipflag values for comparison with other algorithms */
    aero_fptr = fopen ("ipflag.img", "w");
    fwrite (ipflag, npixels, sizeof (uint8), aero_fptr);
    fclose (aero_fptr);

    /* Write the aerosol values for comparison with other algorithms */
    aero_fptr = fopen ("aerosols.img", "w");
    fwrite (taero, npixels, sizeof (float), aero_fptr);
    fclose (aero_fptr);
#endif

    /* Replace the invalid aerosol retrievals (taero and teps) with a local
       average of those values */
    mytime = time(NULL);
    printf ("Filling invalid aerosol values in the 3x3 windows %s",
        ctime(&mytime));
    retval = fix_invalid_aerosols_landsat (ipflag, taero, teps, LAERO_WINDOW,
        LHALF_AERO_WINDOW, nlines, nsamps);
    if (retval != SUCCESS)
    {
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

#ifdef WRITE_TAERO
    /* Write the ipflag values for comparison with other algorithms */
    aero_fptr = fopen ("ipflag_filled.img", "w");
    fwrite (ipflag, npixels, sizeof (uint8), aero_fptr);
    fclose (aero_fptr);

    /* Write the aerosol values for comparison with other algorithms */
    aero_fptr = fopen ("aerosols_filled.img", "w");
    fwrite (taero, npixels, sizeof (float), aero_fptr);
    fclose (aero_fptr);
#endif

    /* Use the center of the aerosol windows to interpolate the remaining
       pixels in the window for taero */
    mytime = time(NULL);
    printf ("Interpolating the aerosol values in the 3x3 windows %s",
        ctime(&mytime));
    aerosol_interp_landsat (xml_metadata, LAERO_WINDOW, LHALF_AERO_WINDOW,
        qaband, ipflag, taero, nlines, nsamps);

#ifdef WRITE_TAERO
    /* Write the ipflag values for comparison with other algorithms */
    aero_fptr = fopen ("ipflag_final.img", "w");
    fwrite (ipflag, npixels, sizeof (uint8), aero_fptr);
    fclose (aero_fptr);

    /* Write the aerosol values for comparison with other algorithms */
    aero_fptr = fopen ("aerosols_final.img", "w");
    fwrite (taero, npixels, sizeof (float), aero_fptr);
    fclose (aero_fptr);
#endif

    /* Use the center of the aerosol windows to interpolate the teps values
       (angstrom coefficient).  The median value used for filling in clouds and
       water will be the default eps value. */
    mytime = time(NULL);
    printf ("Interpolating the teps values in the 3x3 windows %s",
        ctime(&mytime));
    aerosol_interp_landsat (xml_metadata, LAERO_WINDOW, LHALF_AERO_WINDOW,
        qaband, ipflag, teps, nlines, nsamps);

    /* Perform the second level of atmospheric correction using the aerosols */
    mytime = time(NULL);
    printf ("Performing atmospheric correction ... %s", ctime(&mytime));

    /* 0 .. DNL_BAND7 is the same as 0 .. SRL_BAND7 here, since the pan band
       isn't spanned */
    for (ib = 0; ib <= DNL_BAND7; ib++)
    {
#ifdef _OPENMP
        #pragma omp parallel for private (i, rsurf, rotoa, raot550nm, eps, xtv, xmuv, xts, xfi, cosxfi, pres, uwv, uoz, retval, tmpf, roslamb, tgo, roatm, ttatmg, satm, xrorayp)
#endif
        for (i = 0; i < npixels; i++)
        {
            /* If this pixel is fill, then don't process */
            if (level1_qa_is_fill (qaband[i]))
                continue;

            /* Correct all pixels */
            rsurf = sband[ib][i];
            rotoa = (rsurf * bttatmg[ib] / (1.0 - bsatm[ib] * rsurf) +
                broatm[ib]) * btgo[ib];
            raot550nm = taero[i];
            eps = teps[i];

            if (use_orig_aero)
            {
                /* Determine the solar and view angles for the current pixel */
                xtv = vza[i] * 0.01;
                xmuv = cos(xtv * DEG2RAD);
                xts = sza[i] * 0.01;
                xmus = cos(xts * DEG2RAD);
                xfi = saa[i] * 0.01 - vaa[i] * 0.01;
                cosxfi = cos(xfi * DEG2RAD);

                pres = tp[i];
                uwv = twvi[i];
                uoz = tozi[i];
                retval = atmcorlamb2 (input->meta.sat, xts, xtv, xmus, xmuv,
                    xfi, cosxfi, raot550nm, ib, pres, tpres, aot550nm,
                    rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                    normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
                    uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                    wvtransb, oztransa, rotoa, &roslamb, &tgo, &roatm, &ttatmg,
                    &satm, &xrorayp, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian "
                        "atmospheric correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
            }
            else
                atmcorlamb2_new (input->meta.sat, tgo_arr[ib], 
                    aot550nm[roatm_iaMax[ib]], &roatm_coef[ib][0],
                    &ttatmg_coef[ib][0], &satm_coef[ib][0], raot550nm, ib,
                    normext_p0a3_arr[ib], rotoa, &roslamb, eps);

            /* If this is the coastal aerosol band then set the aerosol
               bits in the QA band */
            if (ib == DNL_BAND1)
            {
                /* Set up aerosol QA bits */
                tmpf = fabs (rsurf - roslamb);
                if (tmpf <= LOW_AERO_THRESH)
                {  /* Set the first aerosol bit (low aerosols) */
                    ipflag[i] |= (1 << AERO1_QA);
                }
                else
                {
                    if (tmpf < AVG_AERO_THRESH)
                    {  /* Set the second aerosol bit (average aerosols) */
                        ipflag[i] |= (1 << AERO2_QA);
                    }
                    else
                    {  /* Set both aerosol bits (high aerosols) */
                        ipflag[i] |= (1 << AERO1_QA);
                        ipflag[i] |= (1 << AERO2_QA);
                    }
                }
            }  /* end if this is the coastal aerosol band */

            /* Save the unscaled surface reflectance value */
            if (roslamb < MIN_VALID_REFL)
                sband[ib][i] = MIN_VALID_REFL;
            else if (roslamb > MAX_VALID_REFL)
                sband[ib][i] = MAX_VALID_REFL;
            else
                sband[ib][i] = roslamb;
        }  /* end for i */
    }  /* end for ib */

    /* Free memory for arrays no longer needed */
    if (use_orig_aero)
    {
        free (twvi);
        free (tozi);
        free (tp);
    }
    free (teps);
 
    /* Write the data to the output file */
    mytime = time(NULL);
    printf ("Writing surface reflectance corrected data to the output "
        "files ... %s", ctime(&mytime));

    /* Open the output file */
    sr_output = open_output (xml_metadata, input, OUTPUT_SR, auxnm);
    if (sr_output == NULL)
    {   /* error message already printed */
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through the reflectance bands and write the data */
    for (ib = 0; ib <= DNL_BAND7; ib++)
    {
        /* Scale the output data from float to uint16 */
        convert_output (sband, ib, nlines, nsamps, false, out_band);

        /* Write the scaled product */
        if (put_output_lines (sr_output, out_band, ib, 0, nlines,
            sizeof (uint16)) != SUCCESS)
        {
            sprintf (errmsg, "Writing output data for band %d", ib);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Create the ENVI header file this band */
        if (create_envi_struct (&sr_output->metadata.band[ib],
            &xml_metadata->global, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Creating ENVI header structure.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Write the ENVI header */
        strcpy (envi_file, sr_output->metadata.band[ib].file_name);
        cptr = strchr (envi_file, '.');
        strcpy (cptr, ".hdr");
        if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Writing ENVI header file.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Append the surface reflectance bands (1-7) to the XML file */
    if (append_metadata (7, sr_output->metadata.band, xml_infile) !=
        SUCCESS)
    {
        sprintf (errmsg, "Appending surface reflectance bands to the "
            "XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /** Write the aerosol band **/
    /* Scale the output aerosol data from float to int16 (just use the uint16
       array since it has the same size memory) */
    convert_aerosol_output (taero, qaband, nlines, nsamps, (int16 *) out_band);

    /* Write the aerosol band */
    printf ("  Aerosol Band %d: %s\n", SRL_AEROSOL+1,
            sr_output->metadata.band[SRL_AEROSOL].file_name);
    if (put_output_lines (sr_output, out_band, SRL_AEROSOL, 0, nlines,
        sizeof (int16)) != SUCCESS)
    {
        sprintf (errmsg, "Writing aerosol data");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Free the memory for the taero aerosol band */
    free (taero);

    /* Create the ENVI header for the aerosol band */
    if (create_envi_struct (&sr_output->metadata.band[SRL_AEROSOL],
        &xml_metadata->global, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Creating ENVI header structure.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Write the ENVI header */
    strcpy (envi_file, sr_output->metadata.band[SRL_AEROSOL].file_name);
    cptr = strchr (envi_file, '.');
    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Writing ENVI header file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Append the aerosol band to the XML file */
    if (append_metadata (1, &sr_output->metadata.band[SRL_AEROSOL],
        xml_infile) != SUCCESS)
    {
        sprintf (errmsg, "Appending aerosol band to XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /** Write the aerosol QA band **/
    printf ("  Aerosol QA Band %d: %s\n", SRL_AEROSOL_QA+1,
            sr_output->metadata.band[SRL_AEROSOL_QA].file_name);
    if (put_output_lines (sr_output, ipflag, SRL_AEROSOL_QA, 0, nlines,
        sizeof (uint8)) != SUCCESS)
    {
        sprintf (errmsg, "Writing aerosol QA output data");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Free memory for ipflag data */
    free (ipflag);

    /* Create the ENVI header for the aerosol QA band */
    if (create_envi_struct (&sr_output->metadata.band[SRL_AEROSOL_QA],
        &xml_metadata->global, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Creating ENVI header structure.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Write the ENVI header */
    strcpy (envi_file, sr_output->metadata.band[SRL_AEROSOL_QA].file_name);
    cptr = strchr (envi_file, '.');
    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Writing ENVI header file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Append the aerosol QA band to the XML file */
    if (append_metadata (1, &sr_output->metadata.band[SRL_AEROSOL_QA],
        xml_infile) != SUCCESS)
    {
        sprintf (errmsg, "Appending aerosol QA band to XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the output surface reflectance products */
    close_output (sat, sr_output, OUTPUT_SR);
    free_output (sr_output, OUTPUT_SR);

    /* Free the spatial mapping pointer */
    free (space);

    /* Free the data arrays */
    free (rolutt);
    free (transt);
    free (sphalbt);
    free (normext);
    free (tsmax);
    free (tsmin);
    free (nbfic);
    free (nbfi);
    free (ttv);

    /* Successful completion */
    mytime = time(NULL);
    printf ("Surface reflectance correction complete ... %s\n", ctime(&mytime));
    return (SUCCESS);
}
