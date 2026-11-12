/*****************************************************************************
FILE: subaeroret.c
  
PURPOSE: Contains functions for handling the atmosperic corrections.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/
#include "lut_subr.h"

/******************************************************************************
MODULE:  subaeroret_new

PURPOSE:  Main driver for the atmospheric correction.  This subroutine uses
atmospheric coefficients to determine the atmospheric variables, then performs
the atmospheric corrections.  Updated algorithm to utilize semi-empirical
approach to accessing the look-up table.


RETURN VALUE:
Type = N/A

NOTES:
******************************************************************************/
void subaeroret_new
(
    Sat_t sat,                             /* I: satellite */
    bool water,                            /* I: water pixel flag */
    int iband1,                            /* I: band 1 index (0-based) */
    float erelc[NSR_BANDS],                /* I: band ratio variable */
    float troatm[NSR_BANDS],               /* I: toa reflectance */
    float tgo_arr[NREFL_BANDS],            /* I: per-band other gaseous
                                                 transmittance */
    int roatm_iaMax[NREFL_BANDS],          /* I: roatm_iaMax */
    float roatm_coef[NREFL_BANDS][NCOEF],  /* I: per band polynomial
                                                 coefficients for roatm */
    float ttatmg_coef[NREFL_BANDS][NCOEF], /* I: per band polynomial
                                                 coefficients for ttatmg */
    float satm_coef[NREFL_BANDS][NCOEF],   /* I: per band polynomial
                                                 coefficients for satm */
    float normext_p0a3_arr[NREFL_BANDS],   /* I: normext[iband][0][3] */
    float *raot,     /* O: AOT reflectance */
    float *residual, /* O: model residual */
    int *iaots,      /* I/O: AOT index that is passed in and out for multiple
                             calls (0-based) */
    float eps        /* I: angstroem coefficient; spectral dependency of AOT */
)
{
    int iaot;               /* aerosol optical thickness (AOT) index */
    int ib;                 /* band index */
    int start_band = 0;     /* starting band index for the loop */
    int end_band = 0;       /* ending band index for the loop */
    float raot550nm=0.0;    /* nearest input value of AOT */
    float roslamb;          /* lambertian surface reflectance */
    double ros1;            /* surface reflectance for bands */
    double raot1, raot2;    /* AOT ratios that bracket the predicted ratio */
    float raotsaved;        /* save the raot value */
    double residual1, residual2;  /* residuals for storing and comparing */
    double residualm;       /* local model residual */
    int nbval;              /* number of values meeting criteria */
    bool testth;            /* surface reflectance test variable */
    double xa, xb;          /* AOT ratio values */
    double raotmin;         /* minimum AOT ratio */
    double point_error;     /* residual differences for each pixel */
    int iaot1, iaot2;       /* AOT indices (0-based) */
    float *tth = NULL;      /* pointer to the Landsat or Sentinel tth array */
    float landsat_tth[NSRL_BANDS] = {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 0.0, 0.0,
                                     1.0e-04, 0.0};
    float landsat_tth_water[NSRL_BANDS] =
                                    {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 1.0e-03,
                                     0.0, 1.0e-04, 0.0};
                            /* constant values for comparing against Landsat
                               surface reflectance */
#ifdef PROC_ALL_BANDS
/* Process all bands if turned on */
    float sentinel_tth[NSRS_BANDS] = {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-04};
    float sentinel_tth_water[NSRS_BANDS] =
                                   {1.0e-03, 0.0, 0.0, 1.0e-03, 0.0, 0.0, 0.0,
                                    0.0, 1.0e-3, 0.0, 0.0, 0.0, 1.0e-4};
/* ESPA - I believe the sentinel tth water should be for band 1, 4, 8a, and 12 vs. what I've commented below ... which is a duplicate of the FORTRAN array coeffs.  I think that's a bug.
                                   {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 1.0e-3,
                                    0.0, 1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
*/
                            /* constant values for comparing against Sentinel
                               surface reflectance (different values for land
                               vs. water) */
#else
/* Skip bands 9 and 10 as default for ESPA */
    float sentinel_tth[NSRS_BANDS] = {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 1.0e-04};
    float sentinel_tth_water[NSRS_BANDS] =
                                   {1.0e-03, 0.0, 0.0, 1.0e-03, 0.0, 0.0, 0.0,
                                    0.0, 1.0e-3, 0.0, 1.0e-4};
                            /* constant values for comparing against Sentinel
                               surface reflectance (removed band 9&10) */
#endif
    float aot550nm[NAOT_VALS] = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6,
                                 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3, 2.6,
                                 3.0, 3.5, 4.0, 4.5, 5.0}; /* AOT values */

    /* Initialize variables based on the satellite type */
    if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
    {
        if (water)
            tth = landsat_tth_water;
        else
            tth = landsat_tth;
        start_band = DNL_BAND1;
        end_band = DNL_BAND7;
    }
    else if (sat == SAT_SENTINEL_2)
    {
        if (water)
            tth = sentinel_tth_water;
        else
            tth = sentinel_tth;
        start_band = DNS_BAND1;
        end_band = DNS_BAND12;
    }

    /* Correct input band with increasing AOT (using pre till ratio is equal to
       erelc[2]) */
    iaot = *iaots;
    residual1 = 2000.0;
    residual2 = 1000.0;
    iaot2 = 0;
    iaot1 = 0;
    raot2 = 1.0e-06;
    raot1 = 0.0001;
    ros1 = 1.0;
    raot550nm = aot550nm[iaot];
    testth = false;
    nbval = 0;
    *residual = 0.0;

    /* Atmospheric correction for band 1 */
    ib = iband1;
    atmcorlamb2_new (sat, tgo_arr[ib], aot550nm[roatm_iaMax[ib]],
        &roatm_coef[ib][0], &ttatmg_coef[ib][0], &satm_coef[ib][0], raot550nm,
        ib, normext_p0a3_arr[ib], troatm[ib], &roslamb, eps);

    if (roslamb - tth[iband1] < 0.0)
        testth = true;
    ros1 = roslamb;

    /* Atmospheric correction for each band; handle residual for water-based
       correction differently */
    if (water)
    {  /* expect that we are processing pixels over water */
        for (ib = start_band; ib <= end_band; ib++)
        {
            /* This will process iband1, different from land */
            if (erelc[ib] > 0.0)
            {
                atmcorlamb2_new (sat, tgo_arr[ib], aot550nm[roatm_iaMax[ib]],
                    &roatm_coef[ib][0], &ttatmg_coef[ib][0], &satm_coef[ib][0],
                    raot550nm, ib, normext_p0a3_arr[ib], troatm[ib], &roslamb,
                    eps);

                if (roslamb - tth[ib] < 0.0)
                    testth = true;

                *residual += roslamb*roslamb;
                nbval++;
            }
        }
    }
    else
    {  /* expect that we are processing pixels over land */
        for (ib = start_band; ib <= end_band; ib++)
        {
            /* Don't reprocess iband1 */
            if (ib != iband1 && erelc[ib] > 0.0)
            {
                atmcorlamb2_new (sat, tgo_arr[ib], aot550nm[roatm_iaMax[ib]],
                    &roatm_coef[ib][0], &ttatmg_coef[ib][0], &satm_coef[ib][0],
                    raot550nm, ib, normext_p0a3_arr[ib], troatm[ib], &roslamb,
                    eps);

                if (roslamb - tth[ib] < 0.0)
                    testth = true;

                point_error = roslamb - erelc[ib] * ros1;
                *residual += point_error * point_error;
                nbval++;
            }
        }
    }
    *residual = sqrt (*residual) / nbval;

    /* Loop until we converge on a solution */
    iaot++;
    while ((iaot < NAOT_VALS) && (*residual < residual1) && (!testth))
    {
        /* Reset variables for this loop */
        residual2 = residual1;
        iaot2 = iaot1;
        raot2 = raot1;
        residual1 = *residual;
        raot1 = raot550nm;
        iaot1 = iaot;
        raot550nm = aot550nm[iaot];
        *residual = 0.0;
        nbval = 0;

        /* Atmospheric correction for band 1 */
        ib = iband1;
        testth = false;
        atmcorlamb2_new (sat, tgo_arr[ib], aot550nm[roatm_iaMax[ib]],
            &roatm_coef[ib][0], &ttatmg_coef[ib][0], &satm_coef[ib][0],
            raot550nm, ib, normext_p0a3_arr[ib], troatm[ib], &roslamb, eps);
        if (roslamb - tth[iband1] < 0.0)
            testth = true;
        ros1 = roslamb;

        /* Atmospheric correction for each band; handle residual for water-based
           correction differently */
        if (water)
        {  /* expect that we are processing pixels over water */
            for (ib = start_band; ib <= end_band; ib++)
            {
                /* This will process iband1, different from land */
                if (erelc[ib] > 0.0)
                {
                    atmcorlamb2_new (sat, tgo_arr[ib],
                        aot550nm[roatm_iaMax[ib]], &roatm_coef[ib][0],
                        &ttatmg_coef[ib][0], &satm_coef[ib][0], raot550nm, ib,
                        normext_p0a3_arr[ib], troatm[ib], &roslamb, eps);
    
                    if (roslamb - tth[ib] < 0.0)
                        testth = true;
    
                    *residual += roslamb*roslamb;
                    nbval++;
                }
            }
        }
        else
        {  /* expect that we are processing pixels over land */
            for (ib = start_band; ib <= end_band; ib++)
            {
                /* Don't reprocess iband1 */
                if (ib != iband1 && erelc[ib] > 0.0)
                {
                    atmcorlamb2_new (sat, tgo_arr[ib],
                        aot550nm[roatm_iaMax[ib]], &roatm_coef[ib][0],
                        &ttatmg_coef[ib][0], &satm_coef[ib][0], raot550nm, ib,
                        normext_p0a3_arr[ib], troatm[ib], &roslamb, eps);
    
                    if (roslamb - tth[ib] < 0.0)
                        testth = true;
    
                    point_error = roslamb - erelc[ib] * ros1;
                    *residual += point_error * point_error;
                    nbval++;
                }
            }
        }
        *residual = sqrt (*residual) / nbval;

        /* Move to the next AOT index */
        iaot++;
    }  /* while aot */

    /* If a minimum local was not reached for raot1, then just use the
       raot550nm value.  Otherwise continue to refine the raot. */
    if (iaot == 1)
    {
        *raot = raot550nm;
    }
    else
    {
        /* Refine the AOT ratio.  This is performed by applying a parabolic
           (quadratic) fit to the three (raot, residual) pairs found above:
                   res = a(raot)^2 + b(raot) + c
               The minimum occurs where the first derivative is zero:
                   res' = 2a(raot) + b = 0
                   raot_min = -b/2a

               The a and b coefficients are solved for in the three
               residual equations by eliminating c:
                   r_1 - r = a(raot_1^2 - raot^2) + b(raot_1 - raot)
                   r_2 - r = a(raot_2^2 - raot^2) + b(raot_2 - raot) */
        *raot = raot550nm;
        raotsaved = *raot;
        xa = (residual1 - *residual)*(raot2 - *raot);
        xb = (residual2 - *residual)*(raot1 - *raot);
        raotmin = 0.5*(xa*(raot2 + *raot) - xb*(raot1 + *raot))/(xa - xb);

        /* Validate the min AOT ratio */
        if (raotmin < 0.01 || raotmin > 4.0)
            raotmin = *raot;

        /* Atmospheric correction for band 1 */
        raot550nm = raotmin;
        ib = iband1;
        testth = false;
        residualm = 0.0;
        nbval = 0;
        atmcorlamb2_new (sat, tgo_arr[ib], aot550nm[roatm_iaMax[ib]],
            &roatm_coef[ib][0], &ttatmg_coef[ib][0], &satm_coef[ib][0],
            raot550nm, ib, normext_p0a3_arr[ib], troatm[ib], &roslamb, eps);

        if (roslamb - tth[iband1] < 0.0)
            testth = true;
        ros1 = roslamb;
        if (water && erelc[ib] > 0.0)
        {
            residualm += (roslamb * roslamb);
            nbval++;
        }

        /* Atmospheric correction for each band */
        for (ib = start_band; ib <= end_band; ib++)
        {
            /* Don't reprocess iband1 */
            if (ib != iband1 && erelc[ib] > 0.0)
            {
                atmcorlamb2_new (sat, tgo_arr[ib], aot550nm[roatm_iaMax[ib]],
                    &roatm_coef[ib][0], &ttatmg_coef[ib][0], &satm_coef[ib][0],
                    raot550nm, ib, normext_p0a3_arr[ib], troatm[ib], &roslamb,
                    eps);

                if (roslamb - tth[ib] < 0.0)
                    testth = true;
                if (water)
                    residualm += (roslamb * roslamb);
                else
                {
                    point_error = roslamb - erelc[ib] * ros1;
                    residualm += point_error * point_error;
                }
                nbval++;
            }
        }

        residualm = sqrt (residualm) / nbval;
        *raot = raot550nm;

        /* Check the residuals and reset the AOT ratio */
        if (residualm > *residual)
        {
            residualm = *residual;
            *raot = raotsaved;
        }
        if (residualm > residual1)
        {
            residualm = residual1;
            *raot = raot1;
        }
        if (residualm > residual2)
        {
            residualm = residual2;
            *raot = raot2;
        }
        *residual = residualm;

        /* Check the iaot values */
        if (water && iaot == 1)
            *iaots = 0;
        else
            *iaots = MAX ((iaot2 - 3), 0);
    }
}


/******************************************************************************
MODULE:  subaeroret

PURPOSE:  Main driver for the atmospheric correction.  This subroutine reads
the lookup table (LUT) and performs the atmospheric corrections.  Traditional
FORTRAN algorithm.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading the LUT or doing the correction
SUCCESS        Successful completion

NOTES:
******************************************************************************/
int subaeroret
(
    Sat_t sat,                   /* I: satellite */
    bool water,                  /* I: water pixel flag */
    int iband1,                  /* I: band 1 index (0-based) */
    float xts,                   /* I: solar zenith angle (deg) */
    float xtv,                   /* I: observation zenith angle (deg) */
    float xmus,                  /* I: cosine of solar zenith angle */
    float xmuv,                  /* I: cosine of observation zenith angle */
    float xfi,                   /* I: azimuthal difference between sun and
                                       observation (deg) */
    float cosxfi,                /* I: cosine of azimuthal difference */
    float pres,                  /* I: surface pressure */
    float uoz,                   /* I: total column ozone */
    float uwv,                   /* I: total column water vapor (precipital
                                       water vapor) */
    float erelc[NSR_BANDS],      /* I: band ratio variable */
    float troatm[NSR_BANDS],     /* I: atmospheric reflectance table */
    float tpres[NPRES_VALS],     /* I: surface pressure table */
    float *rolutt,               /* I: intrinsic reflectance table
                          [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float *transt,               /* I: transmission table
                       [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUNANGLE_VALS] */
    float xtsstep,               /* I: solar zenith step value */
    float xtsmin,                /* I: minimum solar zenith value */
    float xtvstep,               /* I: observation step value */
    float xtvmin,                /* I: minimum observation value */
    float *sphalbt,              /* I: spherical albedo table
                                       [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext,              /* I: aerosol extinction coefficient at the
                                       current wavelength (normalized at 550nm)
                                       [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *tsmax,                /* I: maximum scattering angle table
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *tsmin,                /* I: minimum scattering angle table
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfic,                /* I: communitive number of azimuth angles
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfi,                 /* I: number of azimuth angles
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tts[NSOLAR_ZEN_VALS],  /* I: sun angle table */
    int32 indts[NSUNANGLE_VALS], /* I: index for the sun angle table */
    float *ttv,                  /* I: view angle table
                                       [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tauray[NSR_BANDS],     /* I: molecular optical thickness coeff */
    double ogtransa1[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb0[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb1[NSR_BANDS], /* I: other gases transmission coeff */
    double wvtransa[NSR_BANDS],  /* I: water vapor transmission coeff */
    double wvtransb[NSR_BANDS],  /* I: water vapor transmission coeff */
    double oztransa[NSR_BANDS],  /* I: ozone transmission coeff */
    float *raot,                 /* O: AOT reflectance */
    float *residual,             /* O: model residual */
    int *iaots,                  /* I/O: AOT index that is passed in and out
                                         for multiple calls (0-based) */
    float eps                    /* I: angstroem coefficient; spectral
                                       dependency of the AOT */
)
{
    char FUNC_NAME[] = "subaeroret";   /* function name */
    char errmsg[STR_SIZE];       /* error message */
    int iaot;               /* aerosol optical thickness (AOT) index */
    int retval;             /* function return value */
    int ib;                 /* band index */
    int start_band = 0;     /* starting band index for the loop */
    int end_band = 0;       /* ending band index for the loop */
    float raot550nm=0.0;    /* nearest input value of AOT */
    float roslamb;          /* lambertian surface reflectance */
    double ros1;            /* surface reflectance for bands */
    double raot1, raot2;    /* AOT ratios that bracket the predicted ratio */
    double point_error;     /* residual differences for each pixel */
    float raotsaved;        /* save the raot value */
    float tgo;              /* other gaseous transmittance */
    float roatm;            /* intrinsic atmospheric reflectance */
    float ttatmg;           /* total atmospheric transmission */
    float satm;             /* spherical albedo */
    float xrorayp;          /* reflectance of the atmosphere due to molecular
                               (Rayleigh) scattering */
    double residual1, residual2;  /* residuals for storing and comparing */
    double residualm;       /* local model residual */
    int nbval;              /* number of values meeting criteria */
    bool testth;            /* surface reflectance test variable */
    double xa, xb;          /* AOT ratio values */
    double raotmin;         /* minimum AOT ratio */
    int iaot1, iaot2;       /* AOT indices (0-based) */

    float *tth = NULL;      /* pointer to the Landsat or Sentinel tth array */
    float landsat_tth[NSRL_BANDS] = {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 0.0, 0.0,
                                     1.0e-04, 0.0};
    float landsat_tth_water[NSRL_BANDS] =
                                    {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 1.0e-03,
                                     0.0, 1.0e-04, 0.0};
                            /* constant values for comparing against Landsat
                               surface reflectance */
#ifdef PROC_ALL_BANDS
/* Process all bands if turned on */
    float sentinel_tth[NSRS_BANDS] = {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-04};
    float sentinel_tth_water[NSRS_BANDS] =
                                   {1.0e-03, 0.0, 0.0, 1.0e-03, 0.0, 0.0, 0.0,
                                    0.0, 1.0e-3, 0.0, 0.0, 0.0, 1.0e-4};
/* ESPA - I believe the sentinel tth water should be for band 1, 4, 8a, and 12 vs. what I've commented below ... which is a duplicate of the FORTRAN array coeffs.  I think that's a bug.
                                   {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 1.0e-3,
                                    0.0, 1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
*/
                            /* constant values for comparing against Sentinel
                               surface reflectance (different values for land
                               vs. water) */
#else
/* Skip bands 9 and 10 as default for ESPA */
    float sentinel_tth[NSRS_BANDS] = {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 1.0e-04};
    float sentinel_tth_water[NSRS_BANDS] =
                                   {1.0e-03, 0.0, 0.0, 1.0e-03, 0.0, 0.0, 0.0,
                                    0.0, 1.0e-3, 0.0, 1.0e-4};
                            /* constant values for comparing against Sentinel
                               surface reflectance (removed band 9&10) */
#endif
    float aot550nm[NAOT_VALS] = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6,
                                 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3, 2.6,
                                 3.0, 3.5, 4.0, 4.5, 5.0}; /* AOT values */

    /* Initialize variables based on the satellite type */
    if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
    {
        if (water)
            tth = landsat_tth_water;
        else
            tth = landsat_tth;
        start_band = DNL_BAND1;
        end_band = DNL_BAND7;
    }
    else if (sat == SAT_SENTINEL_2)
    {
        if (water)
            tth = sentinel_tth_water;
        else
            tth = sentinel_tth;
        start_band = DNS_BAND1;
        end_band = DNS_BAND12;
    }

    /* Correct input band with increasing AOT (using pre till ratio is equal to
       erelc[2]) */
    iaot = *iaots;
    residual1 = 2000.0;
    residual2 = 1000.0;
    iaot2 = 0;
    iaot1 = 0;
    raot2 = 1.0e-06;
    raot1 = 0.0001;
    ros1 = 1.0;
    raot550nm = aot550nm[iaot];
    testth = false;
    nbval = 0;
    *residual = 0.0;

    /* Atmospheric correction for band 1 */
    ib = iband1;
    retval = atmcorlamb2 (sat, xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm, ib,
        pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
        sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
        uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
        oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp,
        eps);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Performing lambertian atmospheric correction "
            "type 2.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (roslamb - tth[iband1] < 0.0)
        testth = true;
    ros1 = roslamb;

    /* Atmospheric correction for each band; handle residual for water-based
       correction differently */
    if (water)
    {  /* expect that we are processing pixels over water */
        for (ib = start_band; ib <= end_band; ib++)
        {
            /* This will process iband1, different from land */
            if (erelc[ib] > 0.0)
            {
                retval = atmcorlamb2 (sat, xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, ib, pres, tpres, aot550nm, rolutt, transt,
                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                    ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                    oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg,
                    &satm, &xrorayp, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                if (roslamb - tth[ib] < 0.0)
                    testth = true;

                *residual += roslamb*roslamb;
                nbval++;
            }
        }
    }
    else
    {  /* expect that we are processing pixels over land */
        for (ib = start_band; ib <= end_band; ib++)
        {
            /* Don't reprocess iband1 */
            if (ib != iband1 && erelc[ib] > 0.0)
            {
                retval = atmcorlamb2 (sat, xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, ib, pres, tpres, aot550nm, rolutt, transt,
                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                    ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                    oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg,
                    &satm, &xrorayp, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                if (roslamb - tth[ib] < 0.0)
                    testth = true;

                point_error = roslamb - erelc[ib] * ros1;
                *residual += point_error * point_error;
                nbval++;
            }
        }
    }
    *residual = sqrt (*residual) / nbval;

    /* Loop until we converge on a solution */
    iaot++;
    while ((iaot < NAOT_VALS) && (*residual < residual1) && (!testth))
    {
        /* Reset variables for this loop */
        residual2 = residual1;
        iaot2 = iaot1;
        raot2 = raot1;
        residual1 = *residual;
        raot1 = raot550nm;
        iaot1 = iaot;
        raot550nm = aot550nm[iaot];
        *residual = 0.0;
        nbval = 0;

        /* Atmospheric correction for band 1 */
        ib = iband1;
        testth = false;
        retval = atmcorlamb2 (sat, xts, xtv, xmus, xmuv, xfi, cosxfi,
            raot550nm, ib, pres, tpres, aot550nm, rolutt, transt, xtsstep,
            xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax, tsmin, nbfic,
            nbfi, tts, indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
            ogtransb1, wvtransa, wvtransb, oztransa, troatm[ib], &roslamb,
            &tgo, &roatm, &ttatmg, &satm, &xrorayp, eps);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Performing lambertian atmospheric correction "
                "type 2.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (roslamb - tth[iband1] < 0.0)
            testth = true;
        ros1 = roslamb;

        /* Atmospheric correction for each band; handle residual for water-based
           correction differently */
        if (water)
        {  /* expect that we are processing pixels over water */
            for (ib = start_band; ib <= end_band; ib++)
            {
                /* This will process iband1, different from land */
                if (erelc[ib] > 0.0)
                {
                    retval = atmcorlamb2 (sat, xts, xtv, xmus, xmuv, xfi,
                        cosxfi, raot550nm, ib, pres, tpres, aot550nm, rolutt,
                        transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                        normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv,
                        uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1,
                        wvtransa, wvtransb, oztransa, troatm[ib], &roslamb,
                        &tgo, &roatm, &ttatmg, &satm, &xrorayp, eps);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Performing lambertian atmospheric "
                            "correction type 2.");
                        error_handler (true, FUNC_NAME, errmsg);
                        return (ERROR);
                    }
    
                    if (roslamb - tth[ib] < 0.0)
                        testth = true;
    
                    *residual += roslamb*roslamb;
                    nbval++;
                }
            }
        }
        else
        {  /* expect that we are processing pixels over land */
            for (ib = start_band; ib <= end_band; ib++)
            {
                /* Don't reprocess iband1 */
                if (ib != iband1 && erelc[ib] > 0.0)
                {
                    retval = atmcorlamb2 (sat, xts, xtv, xmus, xmuv, xfi,
                        cosxfi, raot550nm, ib, pres, tpres, aot550nm, rolutt,
                        transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                        normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv,
                        uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1,
                        wvtransa, wvtransb, oztransa, troatm[ib], &roslamb,
                        &tgo, &roatm, &ttatmg, &satm, &xrorayp, eps);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Performing lambertian atmospheric "
                            "correction type 2.");
                        error_handler (true, FUNC_NAME, errmsg);
                        return (ERROR);
                    }
    
                    if (roslamb - tth[ib] < 0.0)
                        testth = true;
    
                    point_error = roslamb - erelc[ib] * ros1;
                    *residual += point_error * point_error;
                    nbval++;
                }
            }
        }
        *residual = sqrt (*residual) / nbval;

        /* Move to the next AOT index */
        iaot++;
    }  /* while aot */

    /* If a minimum local was not reached for raot1, then just use the
       raot550nm value.  Otherwise continue to refine the raot. */
    if (iaot == 1)
    {
        *raot = raot550nm;
    }
    else
    {
        /* Refine the AOT ratio.  This is performed by applying a parabolic
           (quadratic) fit to the three (raot, residual) pairs found above:
                   res = a(raot)^2 + b(raot) + c
               The minimum occurs where the first derivative is zero:
                   res' = 2a(raot) + b = 0
                   raot_min = -b/2a

               The a and b coefficients are solved for in the three
               residual equations by eliminating c:
                   r_1 - r = a(raot_1^2 - raot^2) + b(raot_1 - raot)
                   r_2 - r = a(raot_2^2 - raot^2) + b(raot_2 - raot) */
        *raot = raot550nm;
        raotsaved = *raot;
        xa = (residual1 - *residual)*(raot2 - *raot);
        xb = (residual2 - *residual)*(raot1 - *raot);
        raotmin = 0.5*(xa*(raot2 + *raot) - xb*(raot1 + *raot))/(xa - xb);

        /* Validate the min AOT ratio */
        if (raotmin < 0.01 || raotmin > 4.0)
            raotmin = *raot;

        /* Atmospheric correction for band 1 */
        raot550nm = raotmin;
        ib = iband1;
        testth = false;
        nbval = 0;
        residualm = 0.0;
        retval = atmcorlamb2 (sat, xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm,
            ib, pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
            xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
            ttv, uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
            wvtransb, oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg,
            &satm, &xrorayp, eps);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Performing lambertian atmospheric correction "
                "type 2.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (roslamb - tth[iband1] < 0.0)
            testth = true;
        ros1 = roslamb;
        if (water && erelc[ib] > 0.0)
        {
            residualm += (roslamb * roslamb);
            nbval++;
        }

        /* Atmospheric correction for each band */
        for (ib = start_band; ib <= end_band; ib++)
        {
            /* Don't reprocess iband1 */
            if (ib != iband1 && erelc[ib] > 0.0)
            {
                retval = atmcorlamb2 (sat, xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, ib, pres, tpres, aot550nm, rolutt, transt,
                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                    ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                    oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg,
                    &satm, &xrorayp, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                if (roslamb - tth[ib] < 0.0)
                    testth = true;
                if (water)
                    residualm += (roslamb * roslamb);
                else
                {
                    point_error = roslamb - erelc[ib] * ros1;
                    residualm += point_error * point_error;
                }
                nbval++;
            }
        }

        residualm = sqrt (residualm) / nbval;
        *raot = raot550nm;

        /* Check the residuals and reset the AOT ratio */
        if (residualm > *residual)
        {
            residualm = *residual;
            *raot = raotsaved;
        }
        if (residualm > residual1)
        {
            residualm = residual1;
            *raot = raot1;
        }
        if (residualm > residual2)
        {
            residualm = residual2;
            *raot = raot2;
        }
        *residual = residualm;

        /* Check the iaot values */
        if (water && iaot == 1)
            *iaots = 0;
        else
            *iaots = MAX ((iaot2 - 3), 0);
    }

    /* Successful completion */
    return (SUCCESS);
}

