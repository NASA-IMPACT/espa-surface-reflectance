/*****************************************************************************
FILE: lut_subr.c
  
PURPOSE: Contains functions for reading the look-up tables and doing some
of the coefficient computations for the surface reflectance application.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/
#include "lut_subr.h"
#include "hdf.h"
#include "mfhdf.h"

/* Define the full list of band names for Sentinel-2 for the input files */
char SENTINEL_FULL_BANDNAME[SENTINEL_TTL][3] =
    {"1", "2", "3", "4", "5", "6", "7", "8", "8a", "9", "10", "11", "12"};

/* Removed bands 9 and 10 from the Sentinel array */
float landsat_lambda[NREFLL_BANDS] =
    {0.443, 0.480, 0.585, 0.655, 0.865, 1.61, 2.2};

#ifdef PROC_ALL_BANDS
/* Process all bands if turned on */
float sentinel_lambda[NREFLS_BANDS] =
    {0.443, 0.490, 0.560, 0.665, 0.705, 0.740, 0.783, 0.842, 0.865, 0.945,
     1.375, 1.61, 2.19};
#else
/* Skip bands 9 and 10 as default for ESPA */
float sentinel_lambda[NREFLS_BANDS] =
    {0.443, 0.490, 0.560, 0.665, 0.705, 0.740, 0.783, 0.842, 0.865, 1.61, 2.19};
#endif

#define FOUR_PTS 4

/******************************************************************************
MODULE:  atmcorlamb2_new

PURPOSE:  Lambertian atmospheric correction 2, updated to compute the roatm,
ttatmg, and satm from input coefficients.

RETURN VALUE:
Type = N/A

NOTES:
******************************************************************************/
void atmcorlamb2_new
(
    Sat_t sat,                /* I: satellite */
    float tgo,                /* I: other gaseous transmittance  */
    float roatm_upper,        /* I: roatm upper bound poly_fit, given band */
    float roatm_coef[NCOEF],  /* I: poly_fit coefficients for roatm  */
    float ttatmg_coef[NCOEF], /* I: poly_fit coefficients for ttatmg */
    float satm_coef[NCOEF],   /* I: poly_fit coefficients for satm */
    float raot550nm,          /* I: nearest value of AOT */
    int iband,                /* I: band index (0-based) */
    float normext_ib_0_3,     /* I: normext[iband][0][3] */
    float rotoa,              /* I: top of atmosphere reflectance */
    float *roslamb,           /* O: lambertian surface reflectance */
    float eps                 /* I: angstrom coefficient; spectral dependency
                                    of the AOT */
)
{
    float mraot550nm;      /* nearest value of AOT -- modified local variable */
    float mraot550nm_sq;   /* mraot550nm squared */
    float mraot550nm_cube; /* mraot550nm cubed */
    int max_band_indx = 0; /* maximum band index for Landsat or Sentinel */
    float *lambda = NULL;  /* band wavelength pointer for Landsat or Sentinel */
    static const float lambda_sf = 1/0.55; /* lambda scale factor */
    float roatm;           /* intrinsic atmospheric reflectance */
    float ttatmg;          /* total atmospheric transmission */
    float satm;            /* spherical albedo */

    /* Setup Landsat or Sentinel variables */
    if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
    {
        lambda = landsat_lambda;
        max_band_indx = DNL_BAND7;
    }
    else if (sat == SAT_SENTINEL_2)
    {
        lambda = sentinel_lambda;
        max_band_indx = DNS_BAND12;
    }

    /* Modify the AOT value based on the angstrom coefficient and lambda
       values */
    if  (eps < 0.0 || iband > max_band_indx)
        mraot550nm = raot550nm;
    else
        mraot550nm = (raot550nm / normext_ib_0_3) *
            pow (lambda[iband] * lambda_sf, -eps);

    /* Check the upper limit of the modified AOT value */
    if (mraot550nm >= roatm_upper)
        mraot550nm = roatm_upper;

    /* Store the square and cube of the modified AOT value for multiple use */
    mraot550nm_sq = mraot550nm * mraot550nm;
    mraot550nm_cube = mraot550nm_sq * mraot550nm;

    /* Compute the intrinsic atmospheric reflectance from the coefficients */
    roatm = roatm_coef[3] +
            roatm_coef[2] * mraot550nm +
            roatm_coef[1] * mraot550nm_sq +
            roatm_coef[0] * mraot550nm_cube;

    /* Compute the total atmospheric transmission from the coefficients */
    ttatmg = ttatmg_coef[3] +
             ttatmg_coef[2] * mraot550nm +
             ttatmg_coef[1] * mraot550nm_sq +
             ttatmg_coef[0] * mraot550nm_cube;

    /* Compute the spherical albedo from the coefficients */
    satm = satm_coef[3] +
           satm_coef[2] * mraot550nm +
           satm_coef[1] * mraot550nm_sq +
           satm_coef[0] * mraot550nm_cube;

    /* Perform atmospheric correction */
    *roslamb = rotoa - tgo*roatm;
    *roslamb /= tgo*ttatmg + satm*(*roslamb);
}


/******************************************************************************
MODULE:  atmcorlamb2

PURPOSE:  Lambertian atmospheric correction 2.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred doing the atmospheric corrections.
SUCCESS        Successful completion

NOTES:
    1. Standard sea level pressure is 1013 millibars.
******************************************************************************/
int atmcorlamb2
(
    Sat_t sat,                   /* I: satellite */
    float xts,                   /* I: solar zenith angle (deg) */
    float xtv,                   /* I: observation zenith angle (deg) */
    float xmus,                  /* I: cosine of solar zenith angle */
    float xmuv,                  /* I: cosine of observation zenith angle */
    float xfi,                   /* I: azimuthal difference between sun and
                                       observation (deg) */
    float cosxfi,                /* I: cosine of azimuthal difference */
    float raot550nm,             /* I: nearest value of AOT */
    int iband,                   /* I: band index (0-based) */
    float pres,                  /* I: surface pressure */
    float tpres[NPRES_VALS],     /* I: surface pressure table */
    float aot550nm[NAOT_VALS],   /* I: AOT look-up table */
    float *rolutt,               /* I: intrinsic reflectance table
                                       [NSR_BANDS x NPRES_VALS x NAOT_VALS x
                                        NSOLAR_VALS] */
    float *transt,               /* I: transmission table
                                       [NSR_BANDS x NPRES_VALS x NAOT_VALS x
                                        NSUNANGLE_VALS] */
    float xtsstep,               /* I: solar zenith step value */
    float xtsmin,                /* I: minimum solar zenith value */
    float xtvstep,               /* I: observation step value */
    float xtvmin,                /* I: minimum observation value */
    float *sphalbt,              /* I: spherical albedo table
                                       [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext,              /* I: aerosol extinction coefficient at
                                       the current wavelength (normalized
                                       at 550nm)
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
    float uoz,                   /* I: total column ozone */
    float uwv,                   /* I: total column water vapor (precipital
                                       water vapor) */
    float tauray[NSR_BANDS],     /* I: molecular optical thickness coeff */
    double ogtransa1[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb0[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb1[NSR_BANDS], /* I: other gases transmission coeff */
    double wvtransa[NSR_BANDS],  /* I: water vapor transmission coeff */
    double wvtransb[NSR_BANDS],  /* I: water vapor transmission coeff */
    double oztransa[NSR_BANDS],  /* I: ozone transmission coeff */
    float rotoa,                 /* I: top of atmosphere reflectance */
    float *roslamb,              /* O: lambertian surface reflectance */
    float *tgo,                  /* O: other gaseous transmittance */
    float *roatm,                /* O: intrinsic atmospheric reflectance */
    float *ttatmg,               /* O: total atmospheric transmission */
    float *satm,                 /* O: spherical albedo */
    float *xrorayp,              /* O: reflectance of the atmosphere due to
                                       molecular (Rayleigh) scattering */
    float eps                    /* I: angstrom coefficient; spectral
                                       dependency of the AOT */
)
{
    char FUNC_NAME[] = "atmcorlamb2";   /* function name */
    char errmsg[STR_SIZE];  /* error message */
    float xttv;         /* upward transmittance */
    float xtts;         /* downward transmittance */
    float ttatm;        /* total transmission of the atmosphere */
    float tgog;         /* other gases transmission */
    float tgoz;         /* ozone transmission */
    float tgwv;         /* water vapor transmission */
    float tgwvhalf;     /* water vapor transmission, half content */
    float xtaur;        /* rayleigh optical depth for surface pressure */
    float atm_pres;     /* atmospheric pressure at sea level */
    float mraot550nm;   /* nearest value of AOT -- modified local variable */
    int ip;             /* surface pressure looping variable */
    int ip1, ip2;       /* index variables for the surface pressure */
    int iaot;           /* aerosol optical thickness (AOT) looping variable */
    int iaot1, iaot2;   /* index variables for the AOT and spherical albedo
                           arrays */
    int its;            /* index for the sun angle table */
    int itv;            /* index for the view angle table */
    int indx;           /* index for normext array */
    int max_band_indx = 0; /* maximum band index for Landsat or Sentinel */
    float *lambda = NULL;  /* band wavelength pointer for Landsat or Sentinel */
    static const double lambda_sf = 1/0.55; /* lambda scale factor */

    /* Setup Landsat or Sentinel variables */
    if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
    {
        lambda = landsat_lambda;
        max_band_indx = DNL_BAND7;
    }
    else if (sat == SAT_SENTINEL_2)
    {
        lambda = sentinel_lambda;
        max_band_indx = DNS_BAND12;
    }

    /* Modifiy the AOT value based on the angstrom coefficient and lambda
       values */
    if  (eps < 0.0 || iband > max_band_indx)
        mraot550nm = raot550nm;
    else
    {
        indx = iband * NPRES_VALS * NAOT_VALS + 3;
        mraot550nm = (raot550nm / normext[indx]) *
            (pow ((lambda[iband] * lambda_sf), -eps));
    }

    /* Get the pressure and AOT related values for the current surface pressure
       and AOT.  These indices are passed into several functions. */
    /* Look for the appropriate pressure index in the surface pressure table.
       Stop at the second to last item in the table, so that we have the last
       two elements to use as ip1 and ip2, if needed. */
    ip1 = 0;
    for (ip = 0; ip < NPRES_VALS-1; ip++)
    {
        if (pres < tpres[ip])
            ip1 = ip;
    }
    ip2 = ip1 + 1;

    /* Look for the appropriate AOT index in the AOT table.
       Stop at the second to last item in the table, so that we have the last
       two elements to use as iaot1 and iaot2, if needed. */
    iaot1 = 0;
    for (iaot = 0; iaot < 21; iaot++) /* 22 elements in table, stop one short */
    {
        if (mraot550nm > aot550nm[iaot])
            iaot1 = iaot;
    }
    iaot2 = iaot1 + 1;

    /* Determine the index in the view angle table */
    if (xtv <= xtvmin)
        itv = 0;
    else
        itv = (int) ((xtv - xtvmin) / xtvstep + 1.0);

    /* Determine the index in the sun angle table */
    if (xts <= xtsmin)
        its = 0;
    else
        its = (int) ((xts - xtsmin) / xtsstep);
    if (its > 19)
    {
        sprintf (errmsg, "Solar zenith (xts) is too large: %f", xts);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* This routine returns variables for calculating roslamb */
    comproatm (ip1, ip2, iaot1, iaot2, xts, xtv, xmus, xmuv, cosxfi,
        mraot550nm, iband, pres, tpres, aot550nm, rolutt, tsmax, tsmin, nbfic,
        nbfi, tts, indts, ttv, xtsstep, xtvstep, xtvmin, its, itv, roatm);

    /* Compute the transmission for the solar zenith angle */
    comptrans (ip1, ip2, iaot1, iaot2, xts, mraot550nm, iband, pres, tpres,
        aot550nm, transt, xtsstep, xtsmin, tts, &xtts);

    /* Compute the transmission for the observation zenith angle */
    comptrans (ip1, ip2, iaot1, iaot2, xtv, mraot550nm, iband, pres, tpres,
        aot550nm, transt, xtvstep, xtvmin, tts, &xttv);

    /* Compute total transmission (product downward by upward) */
    ttatm = xtts * xttv;

    /* Compute spherical albedo */
    compsalb (ip1, ip2, iaot1, iaot2, mraot550nm, iband, pres, tpres, aot550nm,
        sphalbt, normext, satm);

    /* Compute transmission of water vapor, ozone, and other gasses */
    atm_pres = pres * ONE_DIV_ATMOS_PRES_0;
    comptg (iband, xmus, xmuv, uoz, uwv, atm_pres, ogtransa1, ogtransb0,
        ogtransb1, wvtransa, wvtransb, oztransa, &tgoz, &tgwv, &tgwvhalf,
        &tgog);

    /* Compute rayleigh component (intrinsic reflectance, at p=pres).
       Pressure in the atmosphere is pres / 1013. */
    xtaur = tauray[iband] * atm_pres;
    local_chand (xfi, xmuv, xmus, xtaur, xrorayp);

    /* Perform atmospheric correction */

    *tgo = tgog * tgoz;
    *roatm = (*roatm - *xrorayp)*tgwvhalf + *xrorayp;
    *ttatmg = ttatm * tgwv;
    *roslamb = rotoa/(*tgo) - *roatm;
    *roslamb /= *ttatmg + *satm*(*roslamb);

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  local_chand

PURPOSE:  Computes the atm/molecular reflectance from 0.0 to 1.0, based on the
sun and observation angles.

RETURN VALUE:
Type = None

NOTES:
 1. Here's how the xfd value was originally calculated. Given that these
    are static values, the xfd itself can really be static.
    xdep = 0.0279  // depolarization factor
    xfd = xdep / (2.0 - xdep)
        = 0.014147355
    xfd = (1.0 - xfd) / (1.0 + 2.0 * xfd)
        = .985852645 / 1.02829471
        = .958725777
******************************************************************************/
void local_chand
(
    float xphi,    /* I: azimuthal difference between sun and observation
                         (deg) */
    float xmuv,    /* I: cosine of observation zenith angle */
    float xmus,    /* I: cosine of solar zenith angle */
    float xtau,    /* I: molecular optical depth */
    float *xrray   /* O: molecular reflectance, 0.0 to 1.0 */
)
{
    int i;                             /* looping variable */
    float pl[10];
    float fs0, fs1, fs2;
    float phios;
    float xcosf2, xcosf3;
    float xph1, xph2, xph3;
    float xitm;
    float xp1, xp2, xp3;
    float cfonc1, cfonc2, cfonc3;
    float xlntau;                      /* log molecular optical depth */
    float xitot1, xitot2, xitot3;
    float xmus2, xmuv2;                /* square of xmus and xmuv */

    /* constant vars */
    const float xfd = 0.958725777;
    const float as0[10] = {
         0.33243832, -6.777104e-02, 0.16285370, 1.577425e-03,
        -0.30924818, -1.240906e-02, -0.10324388, 3.241678e-02, 0.11493334,
        -3.503695e-02};
    const float as1[2] = {0.19666292, -5.439061e-02};
    const float as2[2] = {0.14545937, -2.910845e-02};

    phios = xphi * DEG2RAD;
    xcosf2 = -cos (phios);
    xcosf3 = cos (2.0 * phios);

    /* xmus and xmuv squared is used frequently */
    xmus2 = xmus * xmus;
    xmuv2 = xmuv * xmuv;

    xph1 = 1.0 + (3.0 * xmus2 - 1.0) * (3.0 * xmuv2 - 1.0) * xfd * 0.125;
    xph3 = (1.0 - xmus2) * (1.0 - xmuv2);
    xph2 = -xmus * xmuv * sqrt(xph3);
    xph2 = xph2 * xfd * 0.75;
    xph3 = xph3 * xfd * 0.1875;

    /* The original xitm below has an xmus factor that has been removed from
       the following in order to reduce multiplications and the final
       division by xmus in xrray. */
    xitm = (1.0 - exp(-xtau * (1.0 / xmus + 1.0 / xmuv)))/(4*(xmus + xmuv));
    xp1 = xph1 * xitm;
    xp2 = xph2 * xitm;
    xp3 = xph3 * xitm;

    xitm = (1.0 - exp(-xtau / xmus)) * (1.0 - exp(-xtau / xmuv));
    cfonc1 = xph1 * xitm;
    cfonc2 = xph2 * xitm;
    cfonc3 = xph3 * xitm;

    xlntau = log (xtau);
    pl[0] = 1.0;
    pl[1] = xlntau;
    pl[2] = xmus + xmuv;
    pl[3] = xlntau * pl[2];
    pl[4] = xmus * xmuv;
    pl[5] = xlntau * pl[4];
    pl[6] = xmus2 + xmuv2;
    pl[7] = xlntau * pl[6];
    pl[8] = xmus2 * xmuv2;
    pl[9] = xlntau * pl[8];

    fs0 = 0.0;
    for (i = 0; i < 10; i++)
        fs0 += pl[i] * as0[i];
    fs1 = pl[0] * as1[0] + pl[1] * as1[1];
    fs2 = pl[0] * as2[0] + pl[1] * as2[1];
    xitot1 = xp1 + cfonc1 * fs0;
    xitot2 = xp2 + cfonc2 * fs1;
    xitot3 = xp3 + cfonc3 * fs2;

    *xrray = xitot1 + 2*(xitot2*xcosf2 + xitot3*xcosf3);
}


/******************************************************************************
MODULE:  comptg

PURPOSE:  Computes the transmission of the water vapor, ozone, and other gases.

RETURN VALUE:
Type = N/A

NOTES:
1. Standard sea level pressure is 1013 millibars.
******************************************************************************/
void comptg
(
    int iband,                   /* I: band index (0-based) */
    float xmus,                  /* I: cosine of solar zenith angle */
    float xmuv,                  /* I: cosine of view zenith angle */
    float uoz,                   /* I: total column ozone */
    float uwv,                   /* I: total column water vapor (precipital
                                       water vapor) */
    float atm_pres,              /* I: pressure at sea level */
    double ogtransa1[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb0[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb1[NSR_BANDS], /* I: other gases transmission coeff */
    double wvtransa[NSR_BANDS],  /* I: water vapor transmission coeff */
    double wvtransb[NSR_BANDS],  /* I: water vapor transmission coeff */
    double oztransa[NSR_BANDS],  /* I: ozone transmission coeff */
    float *tgoz,                 /* O: ozone transmission */
    float *tgwv,                 /* O: water vapor transmission */
    float *tgwvhalf,             /* O: water vapor transmission, half content */
    float *tgog                  /* O: other gases transmission */
)
{
    float a, b;  /* water vapor transmission coefficient */
    float m;     /* ozone transmission coefficient */
    float x;     /* water vapor transmission coefficient */

    /* Compute ozone transmission */
    m = 1.0 / xmus + 1.0 / xmuv;
    *tgoz = exp(oztransa[iband] * m * uoz);

    /* Compute water vapor transmission */
    a = wvtransa[iband];
    b = wvtransb[iband];

    x = m * uwv;
    if (x > 1.0E-06)
        *tgwv = exp(-a * pow(x, b));
    else
        *tgwv = 1.0;

    /* Compute water vapor transmission half the content */
    x *= 0.5;
    if (x > 1.0E-06)
        *tgwvhalf = exp(-a * pow(x, b));
    else
        *tgwvhalf = 1.0;

    /* Compute other gases transmission */
    *tgog = -(ogtransa1[iband] * atm_pres) *
        pow(m, exp(-(ogtransb0[iband] + ogtransb1[iband] * atm_pres)));
    *tgog = exp(*tgog);
}


/******************************************************************************
MODULE:  compsalb

PURPOSE:  Computes spherical albedo

RETURN VALUE:
Type = N/A

NOTES:
******************************************************************************/
void compsalb
(
    int ip1,            /* I: index variable for surface pressure */
    int ip2,            /* I: index variable for surface pressure */
    int iaot1,          /* I: index variable for AOT */
    int iaot2,          /* I: index variable for AOT */
    float raot550nm,    /* I: nearest value of AOT */
    int iband,          /* I: band index (0-based) */
    float pres,         /* I: surface pressure */
    float tpres[NPRES_VALS],   /* I: surface pressure table */
    float aot550nm[NAOT_VALS], /* I: AOT look-up table */
    float *sphalbt,     /* I: spherical albedo table
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext,     /* I: aerosol extinction coefficient at the current
                              wavelength (normalized at 550nm) 
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *satm         /* O: spherical albedo */
)
{
    float xtiaot1, xtiaot2;         /* spherical albedo trans value */
    float satm1, satm2;             /* spherical albedo value */
    float dpres;                    /* pressure ratio */
    float deltaaot;                 /* AOT ratio */
    int iband_indx;  /* index of the current iband */
    int ip1_indx;    /* index of the current ip1 (without the band) */
    int ip2_indx;    /* index of the current ip2 (without the band) */
    int iband_ip1_iaot1_indx;  /* index for current band, ip1, iaot1 */
    int iband_ip1_iaot2_indx;  /* index for current band, ip1, iaot2 */
    int iband_ip2_iaot1_indx;  /* index for current band, ip2, iaot1 */
    int iband_ip2_iaot2_indx;  /* index for current band, ip2, iaot2 */

    /* Compute the delta AOT */
    deltaaot = raot550nm - aot550nm[iaot1];
    deltaaot /= aot550nm[iaot2] - aot550nm[iaot1];

    /* Compute the ipX and iaotX indices for transt */
    iband_indx = iband * NPRES_VALS * NAOT_VALS;
    ip1_indx = ip1 * NAOT_VALS;
    ip2_indx = ip2 * NAOT_VALS;
    iband_ip1_iaot1_indx = iband_indx + ip1_indx + iaot1;
    iband_ip1_iaot2_indx = iband_indx + ip1_indx + iaot2;
    iband_ip2_iaot1_indx = iband_indx + ip2_indx + iaot1;
    iband_ip2_iaot2_indx = iband_indx + ip2_indx + iaot2;

    /* Compute the spherical albedo */
    xtiaot1 = sphalbt[iband_ip1_iaot1_indx];
    xtiaot2 = sphalbt[iband_ip1_iaot2_indx];
    satm1 = xtiaot1 + (xtiaot2 - xtiaot1) * deltaaot;

    xtiaot1 = sphalbt[iband_ip2_iaot1_indx];
    xtiaot2 = sphalbt[iband_ip2_iaot2_indx];
    satm2 = xtiaot1 + (xtiaot2 - xtiaot1) * deltaaot;

    dpres = (pres - tpres[ip1]) / (tpres[ip2] - tpres[ip1]);
    *satm = satm1 + (satm2 - satm1) * dpres;
}


/******************************************************************************
MODULE:  comptrans

PURPOSE:  Compute transmission

RETURN VALUE:
Type = none

NOTES:
1. This is called by subaeroret for both the solar zenith angle and the
   observation zenith angle.  Thus, xts is not specific to the solar zenith
   angle in this case.
2. This function is heavily dependent upon the input solar zenith and
   observation zenith angles.  At the current time, these are static values
   in the overall application.  Knowing that, speedup is achievable in this
   routine if these values never change.  However, the long-term goal is to
   change the main function to compute the solar zenith angle on a per-pixel
   basis.  Therefore we will leave this routine as-is.
3. This function is also dependent upon surface pressure and AOT.
******************************************************************************/
void comptrans
(
    int ip1,            /* I: index variable for surface pressure */
    int ip2,            /* I: index variable for surface pressure */
    int iaot1,          /* I: index variable for AOT */
    int iaot2,          /* I: index variable for AOT */
    float xts,          /* I: zenith angle */
    float raot550nm,    /* I: nearest value of AOT */
    int iband,          /* I: band index (0-based) */
    float pres,         /* I: surface pressure */
    float tpres[NPRES_VALS],   /* I: surface pressure table */
    float aot550nm[NAOT_VALS], /* I: AOT look-up table */
    float *transt,      /* I: transmission table
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS x
                               NSUNANGLE_VALS] */
    float xtsstep,      /* I: zenith angle step value */
    float xtsmin,       /* I: minimum zenith angle value */
    float tts[NSOLAR_ZEN_VALS], /* I: sun angle table */
    float *xtts         /* O: downward transmittance */
)
{
    char FUNC_NAME[] = "comptrans"; /* function name */
    char errmsg[STR_SIZE];          /* error message */
    float xtiaot1, xtiaot2;         /* spherical albedo trans value */
    float xtts1, xtts2;
    float xmts, xtranst;
    float dpres;                    /* pressure ratio */
    float deltaaot;                 /* AOT ratio */
    int its;                        /* index for the sun angle table */
    int iband_indx;  /* index of the current iband */
    int ip1_indx;    /* index of the current ip1 (without the band) */
    int ip2_indx;    /* index of the current ip2 (without the band) */
    int iaot1_indx;  /* index of the current iaot1 (without the band & ip) */
    int iaot2_indx;  /* index of the current iaot2 (without the band & ip) */
    int iband_ip1_iaot1_indx;  /* index for current band, ip1, iaot1 */
    int iband_ip1_iaot2_indx;  /* index for current band, ip1, iaot2 */
    int iband_ip2_iaot1_indx;  /* index for current band, ip2, iaot1 */
    int iband_ip2_iaot2_indx;  /* index for current band, ip2, iaot2 */

    /* Determine the index in the sun angle table */
    if (xts <= xtsmin) 
        its = 0;
    else
        its = (int) ((xts - xtsmin) / xtsstep);
    if (its > 19)
    {
        sprintf (errmsg, "Zenith angle (xts) is too large: %f", xts);
        error_handler (true, FUNC_NAME, errmsg);
        return;
    }

    /* Compute the ipX and iaotX indices for transt */
    iband_indx = iband * NPRES_VALS * NAOTxNSUNANGLE_VALS;
    ip1_indx = ip1 * NAOTxNSUNANGLE_VALS;
    ip2_indx = ip2 * NAOTxNSUNANGLE_VALS;
    iaot1_indx = iaot1 * NSUNANGLE_VALS;
    iaot2_indx = iaot2 * NSUNANGLE_VALS;
    iband_ip1_iaot1_indx = iband_indx + ip1_indx + iaot1_indx;
    iband_ip1_iaot2_indx = iband_indx + ip1_indx + iaot2_indx;
    iband_ip2_iaot1_indx = iband_indx + ip2_indx + iaot1_indx;
    iband_ip2_iaot2_indx = iband_indx + ip2_indx + iaot2_indx;

    /* Compute for ip1, iaot1 */
    xmts = (xts - tts[its]) * 0.25;
    xtranst = transt[iband_ip1_iaot1_indx + its];
    xtiaot1 = xtranst + (transt[iband_ip1_iaot1_indx + its + 1] - xtranst) *
        xmts;

    /* Compute for ip1, iaot2 */
    xtranst = transt[iband_ip1_iaot2_indx + its];
    xtiaot2 = xtranst + (transt[iband_ip1_iaot2_indx + its + 1] - xtranst) *
        xmts;

    deltaaot = raot550nm - aot550nm[iaot1];
    deltaaot /= aot550nm[iaot2] - aot550nm[iaot1];
    xtts1 = xtiaot1 + (xtiaot2 - xtiaot1) * deltaaot;

    /* Compute for ip2, iaot1 */
    xtranst = transt[iband_ip2_iaot1_indx + its];
    xtiaot1 = xtranst + (transt[iband_ip2_iaot1_indx + its + 1] - xtranst) *
        xmts;

    /* Compute for ip2, iaot2 */
    xtranst = transt[iband_ip2_iaot2_indx + its];
    xtiaot2 = xtranst + (transt[iband_ip2_iaot2_indx + its + 1] - xtranst) *
        xmts;
    xtts2 = xtiaot1 + (xtiaot2 - xtiaot1) * deltaaot;

    dpres = (pres - tpres[ip1]) / (tpres[ip2] - tpres[ip1]);
    *xtts = xtts1 + (xtts2 - xtts1) * dpres;
}


/******************************************************************************
MODULE:  interp_refl_using_scat_angle

PURPOSE:  Interpolate the reflectance as a function of the scattering angle
given four points in the scattering grid.

RETURN VALUE:
Type = float
Value          Description
-----          -----------
any            interpolated reflectance
******************************************************************************/
static float interp_refl_using_scat_angle
(
    int its,          /* I: index for the sun angle table */
    int itv,          /* I: index for the view angle table */
    float *xtsmax,    /* I: tsmax[itv][its], tsmax[itv][its+1],
                            tsmax[itv+1][its], tsmax[itv+1][its+1] */
    float *xtsmin,    /* I: tsmin[itv][its], tsmin[itv][its+1],
                            tsmin[itv+1][its], tsmin[itv+1][its+1] */
    float scaa,       /* I: scattering angle */
    float *nbfic,     /* I: nbfic values */
    float *nbfi,      /* I: nbfi values */
    int32 *indts,     /* I: index for the sun angle table */
    float *rolutt,    /* I: intrinsic reflectance table
                            [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS]*/
    int rolutt_indx,  /* I: index into rolutt array */
    float t,          /* I: sun angle interpolation parameter */
    float u           /* I: view angle interpolation parameter */
)
{
    int i;           /* loop counter */
    int isca;        /* rolutt array index adjustment */
    float sca1;
    float sca2;
    float roinf;
    float rosup;
    float ro[FOUR_PTS]; /* reflectance at four points */
    float ro_interp;    /* interpolated value */

    for (i = 0; i < FOUR_PTS; i++)
    {
        int is = its + i%2; /* its or its + 1 array index */
        int iv = (i < 2) ? itv : itv + 1; /* itv or itv + 1 array index */
        int j = indts[is] + nbfic[i] - nbfi[i]; /* convenience rolutt array
                                                   index offset */

        if (is != 0 && iv != 0)
        {
            isca = (int) ((xtsmax[i] - scaa) * 0.25 + 1); /* * 0.25 vs / 4.0 */
            if (isca <= 0)
                isca = 1;
            if (isca + 1 < nbfi[i])
            {
                sca1 = xtsmax[i] - (isca - 1) * 4.0;
                sca2 = sca1 - 4;
            }
            else
            {
                isca = nbfi[i] - 1;
                sca1 = xtsmax[i] - (isca - 1) * 4.0;
                sca2 = xtsmin[i];
            }

            roinf = rolutt[rolutt_indx + j + isca - 1];
            rosup = rolutt[rolutt_indx + j + isca];
            ro[i] = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
        }
        else
        {
            ro[i] = rolutt[rolutt_indx + j];
        }
    }

    ro_interp = ro[3]
              + u*(ro[1] - ro[3])
              + t*(ro[2] - ro[3])
              + u*t*(ro[0] - ro[1] - ro[2] + ro[3]);
    return ro_interp;
}


/******************************************************************************
MODULE:  comproatm

PURPOSE:  Computes the atmospheric reflectance

RETURN VALUE:
Type = none

NOTES:
1. This function is heavily dependent upon the solar zenith and observation
   zenith angles.  At the current time, these are static values in the
   overall application.  Knowing that, speedup is achievable in this routine
   if these values never change.  However, the long-term goal is to change
   the main function to compute the solar zenith angle on a per-pixel basis.
   Therefore we will leave this routine as-is.
2. This function is also dependent upon surface pressure and AOT.
******************************************************************************/
void comproatm
(
    int ip1,          /* I: index variable for surface pressure */
    int ip2,          /* I: index variable for surface pressure */
    int iaot1,        /* I: index variable for AOT */
    int iaot2,        /* I: index variable for AOT */
    float xts,        /* I: solar zenith angle (deg) */
    float xtv,        /* I: observation zenith angle (deg) */
    float xmus,       /* I: cosine of solar zenith angle */
    float xmuv,       /* I: cosine of observation zenith angle */
    float cosxfi,     /* I: cosine of azimuthal difference */
    float raot550nm,  /* I: nearest value of AOT */
    int iband,        /* I: band index (0-based) */
    float pres,       /* I: surface pressure */
    float tpres[NPRES_VALS],   /* I: surface pressure table */
    float aot550nm[NAOT_VALS], /* I: AOT look-up table */
    float *rolutt,    /* I: intrinsic reflectance table
                            [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS]*/
    float *tsmax,     /* I: maximum scattering angle table
                            [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *tsmin,     /* I: minimum scattering angle table
                            [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfic,     /* I: communitive number of azimuth angles
                            [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfi,      /* I: number of azimuth angles
                            [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tts[NSOLAR_ZEN_VALS],  /* I: sun angle table */
    int32 indts[NSUNANGLE_VALS], /* I: index for the sun angle table */
    float *ttv,       /* I: view angle table
                            [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float xtsstep,    /* I: solar zenith step value */
    float xtvstep,    /* I: observation step value */
    float xtvmin,     /* I: minimum observation value */
    int its,          /* I: index for the sun angle table */
    int itv,          /* I: index for the view angle table */
    float *roatm      /* O: intrinsic atmospheric reflectance */
)
{
    float nbficl[FOUR_PTS]; /* local array of nbfic values */
    float nbfil[FOUR_PTS];  /* local array of nbfi values */
    float rop1, rop2;       /* reflectance at p1 and p2 */
    float xtsmax[FOUR_PTS]; /* tsmax[itv][its], tsmax[itv][its+1],
                               tsmax[itv+1][its], tsmax[itv+1][its+1] */
    float xtsmin[FOUR_PTS]; /* tsmin[itv][its], tsmin[itv][its+1],
                               tsmin[itv+1][its], tsmin[itv+1][its+1] */
    float cscaa;
    float scaa;                     /* scattering angle */
    float dpres;                    /* pressure ratio */
    float deltaaot;                 /* AOT ratio */
    float roiaot1, roiaot2;
    float t, u;
    float logaot550nm[22] =
        {-4.605170186, -2.995732274, -2.302585093,
         -1.897119985, -1.609437912, -1.203972804,
         -0.916290732, -0.510825624, -0.223143551,
          0.000000000, 0.182321557, 0.336472237,
          0.470003629, 0.587786665, 0.693157181,
          0.832909123, 0.955511445, 1.098612289,
          1.252762969, 1.386294361, 1.504077397,
          1.609437912};
    int iband_indx;  /* index of the current iband */
    int ip1_indx;    /* index of the current ip1 (without the band) */
    int ip2_indx;    /* index of the current ip2 (without the band) */
    int iaot1_indx;  /* index of the current iaot1 (without the band & ip) */
    int iaot2_indx;  /* index of the current iaot2 (without the band & ip) */
    int iband_ip_iaot_indx;  /* index for current band, ip, iaot */
    int itv_its_indx;          /* index for [itv][its] */
    int itv_its_indx_1;        /* itv_its_indx + 1 */
    int itv1_its_indx;         /* index for [itv+1][its] */
    int itv1_its_indx_1;       /* itv1_its_indx + 1 */

    /* Initialize some variables */
    cscaa = -xmus * xmuv - cosxfi * sqrt(1.0 - xmus * xmus) *
        sqrt(1.0 - xmuv * xmuv);
    scaa = acos(cscaa) * RAD2DEG;    /* vs / DEG2RAD */

    itv_its_indx = itv*NSOLAR_ZEN_VALS + its;
    itv_its_indx_1 = itv_its_indx + 1;
    itv1_its_indx = itv_its_indx + NSOLAR_ZEN_VALS;
    itv1_its_indx_1 = itv1_its_indx + 1;

    nbficl[0] = nbfic[itv_its_indx];
    nbfil[0] = nbfi[itv_its_indx];
    nbficl[1] = nbfic[itv_its_indx_1];
    nbfil[1] = nbfi[itv_its_indx_1];
    nbficl[2] = nbfic[itv1_its_indx];
    nbfil[2] = nbfi[itv1_its_indx];
    nbficl[3] = nbfic[itv1_its_indx_1];
    nbfil[3] = nbfi[itv1_its_indx_1];

    xtsmax[0] = tsmax[itv_its_indx];
    xtsmin[0] = tsmin[itv_its_indx];
    xtsmax[1] = tsmax[itv_its_indx_1];
    xtsmin[1] = tsmin[itv_its_indx_1];
    xtsmax[2] = tsmax[itv1_its_indx];
    xtsmin[2] = tsmin[itv1_its_indx];
    xtsmax[3] = tsmax[itv1_its_indx_1];
    xtsmin[3] = tsmin[itv1_its_indx_1];

    iband_indx = iband * NPRES_VALS * NAOTxNSOLAR_VALS;
    ip1_indx = ip1 * NAOTxNSOLAR_VALS;
    ip2_indx = ip2 * NAOTxNSOLAR_VALS;
    iaot1_indx = iaot1 * NSOLAR_VALS;
    iaot2_indx = iaot2 * NSOLAR_VALS;

    t = (tts[its+1] - xts) / (tts[its+1] - tts[its]);
    u = (ttv[itv1_its_indx] - xtv) / (ttv[itv1_its_indx] - ttv[itv_its_indx]);

    /* Interpolate points vs scattering angle for ip1, iaot1. */
    iband_ip_iaot_indx = iband_indx + ip1_indx + iaot1_indx;
    roiaot1 = interp_refl_using_scat_angle(its, itv, xtsmax, xtsmin, scaa,
                                           nbficl, nbfil, indts, rolutt,
                                           iband_ip_iaot_indx, t, u);

    /* Interpolate points vs scattering angle for ip1, iaot2. */
    iband_ip_iaot_indx = iband_indx + ip1_indx + iaot2_indx;
    roiaot2 = interp_refl_using_scat_angle(its, itv, xtsmax, xtsmin, scaa,
                                           nbficl, nbfil, indts, rolutt,
                                           iband_ip_iaot_indx, t, u);

    /* Interpolation as log of tau */
    deltaaot = logaot550nm[iaot2] - logaot550nm[iaot1];
    deltaaot = (log (raot550nm) - logaot550nm[iaot1]) / deltaaot;
    rop1 = roiaot1 + (roiaot2 - roiaot1) * deltaaot;

    /* Interpolate points vs scattering angle for ip2, iaot1. */
    iband_ip_iaot_indx = iband_indx + ip2_indx + iaot1_indx;
    roiaot1 = interp_refl_using_scat_angle(its, itv, xtsmax, xtsmin, scaa,
                                           nbficl, nbfil, indts, rolutt,
                                           iband_ip_iaot_indx, t, u);

    /* Interpolate points vs scattering angle for ip2, iaot2. */
    iband_ip_iaot_indx = iband_indx + ip2_indx + iaot2_indx;
    roiaot2 = interp_refl_using_scat_angle(its, itv, xtsmax, xtsmin, scaa,
                                           nbficl, nbfil, indts, rolutt,
                                           iband_ip_iaot_indx, t, u);

    /* Interpolation as log of tau */
    rop2 = roiaot1 + (roiaot2 - roiaot1) * deltaaot;

    dpres = (pres - tpres[ip1]) / (tpres[ip2] - tpres[ip1]);
    *roatm = rop1 + (rop2 - rop1) * dpres;
}


/******************************************************************************
MODULE:  readluts

PURPOSE:  Reads the look-up tables and input atmospheric files

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading the look-up tables or atmospheric files
SUCCESS        Successful completion

NOTES:
******************************************************************************/
int readluts
(
    Sat_t sat,                  /* I: satellite */
    float *tsmax,               /* O: maximum scattering angle table
                                      [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *tsmin,               /* O: minimum scattering angle table
                                      [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *ttv,                 /* O: view angle table
                                      [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tts[NSOLAR_ZEN_VALS], /* O: sun angle table */
    float *nbfic,               /* O: communitive number of azimuth angles
                                      [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfi,                /* O: number of azimuth angles
                                      [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    int32 indts[NSUNANGLE_VALS],/* O: index for the sun angle table */
    float *rolutt,              /* O: intrinsic reflectance table
                                      [NSR_BANDS x NPRES_VALS x NAOT_VALS x
                                       NSOLAR_VALS] */
    float *transt,              /* O: transmission table
                                      [NSR_BANDS x NPRES_VALS x NAOT_VALS x
                                       NSUNANGLE_VALS] */
    float *sphalbt,             /* O: spherical albedo table
                                      [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext,             /* O: aerosol extinction coefficient at the
                                      current wavelength (normalized at 550nm) 
                                      [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float xtsstep,              /* I: solar zenith step value */
    float xtsmin,               /* I: minimum solar zenith value */
    char anglehdf[STR_SIZE],    /* I: angle HDF filename */
    char intrefnm[STR_SIZE],    /* I: intrinsic reflectance filename */
    char transmnm[STR_SIZE],    /* I: transmission filename */
    char spheranm[STR_SIZE]     /* I: spherical albedo filename */
)
{
    char FUNC_NAME[] = "readluts";   /* function name */
    char errmsg[STR_SIZE];  /* error message */
    char tmpstr[STR_SIZE];  /* temporary string variable, not use */
    int i, j;               /* looping variables */
    int nsr_bands = 0;      /* number of SR bands in the input file;
                               for Sentinel the number of bands in the input
                               file might be different than the number of bands
                               we want to store in the array, since we may be
                               skipping bands 9 and 10 */
    int iband;              /* band looping variable */
    int ibndx;              /* index of the current band for writing to array */
    int iaot;               /* aerosol optical thickness (AOT) index */
    int ipres;              /* looping variable for pressure */
    int itau;               /* looping variable for molecular optical thick */
    int ival;               /* looping variable for LUT */
#ifndef PROC_ALL_BANDS
    int iline;              /* current line for consuming, but ignoring, lines
                               in the input ASCII file for Sentinel-2 */
#endif
    int status;             /* return status of the HDF function */
    int start[3];           /* starting point to read SDS data */
    int edges[3];           /* number of values to read in SDS data */
    char fname[STR_SIZE];   /* filename to be read */
    float *rolut = NULL;    /* intrinsic reflectance read from HDF file
                               [NSOLAR_VALS * NAOT_VALS * NPRES_VALS] */
    float ttsr[22];
    float xx;               /* temporary float values, not used */
    int sd_id;              /* file ID for the HDF file */
    int sds_id;             /* ID for the current SDS */
    int sds_index;          /* index for the current SDS */
    int iband_indx;         /* index of current location in the overall array */
    int ipres_indx;         /* index of current pressure (without the band) */
    int itau_indx;          /* index of current itau (without the band & ip */
    int curr_indx;          /* index of current pixel */
    int tmp_lut_indx;       /* index for the LUT */
    FILE *fp = NULL;        /* file pointer for reading ascii files */

    /* Setup Landsat or Sentinel number of SR bands to be read from the input
       LUT files; number of input bands for Sentinel is possibly more than we
       will actually store since we might be skipping bands 9 and 10 */
    if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
        nsr_bands = NSRL_BANDS;
    else if (sat == SAT_SENTINEL_2)
        nsr_bands = SENTINEL_TTL;

    /* Initialize some variables */
    for (i = 0; i < NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS; i++)
        nbfic[i] = 0.0;
    for (j = 0; j < NSUNANGLE_VALS; j++)
        tts[j] = xtsmin + xtsstep * j;

    /* Open as HDF file for reading */
    sd_id = SDstart (anglehdf, DFACC_RDONLY);
    if (sd_id < 0)
    {
        sprintf (errmsg, "Unable to open %s for reading as SDS", anglehdf);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the 2D bands from the angle HDF file */
    start[0] = 0;   /* lines */
    start[1] = 0;   /* samples */
    edges[0] = NVIEW_ZEN_VALS;   /* number of lines */
    edges[1] = NSOLAR_ZEN_VALS;  /* number of samples */

    /* Find the TSMAX SDS */
    sds_index = SDnametoindex (sd_id, "TSMAX");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find TSMAX in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access TSMAX for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    edges[0] = 1;   /* number of lines */
    for (i = 0; i < NVIEW_ZEN_VALS; i++)
    {
        start[0] = i;   /* lines */
        status = SDreaddata (sds_id, start, NULL, edges,
            &tsmax[i*NSOLAR_ZEN_VALS]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the TSMAX SDS");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to TSMAX SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the TSMIN SDS */
    sds_index = SDnametoindex (sd_id, "TSMIN");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find TSMIN in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access TSMIN for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    edges[0] = 1;   /* number of lines */
    for (i = 0; i < NVIEW_ZEN_VALS; i++)
    {
        start[0] = i;   /* lines */
        status = SDreaddata (sds_id, start, NULL, edges,
            &tsmin[i*NSOLAR_ZEN_VALS]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the TSMIN SDS");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to TSMIN SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the TTV SDS */
    sds_index = SDnametoindex (sd_id, "TTV");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find TTV in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access TTV for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    edges[0] = 1;   /* number of lines */
    for (i = 0; i < NVIEW_ZEN_VALS; i++)
    {
        start[0] = i;   /* lines */
        status = SDreaddata (sds_id, start, NULL, edges,
            &ttv[i*NSOLAR_ZEN_VALS]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the TTV SDS");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to TTV SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the NBFI SDS */
    sds_index = SDnametoindex (sd_id, "NBFI");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find NBFI in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access NBFI for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    edges[0] = 1;   /* number of lines */
    for (i = 0; i < NVIEW_ZEN_VALS; i++)
    {
        start[0] = i;   /* lines */
        status = SDreaddata (sds_id, start, NULL, edges,
            &nbfi[i*NSOLAR_ZEN_VALS]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the NBFI SDS");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to NBFI SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the NBFIC SDS */
    sds_index = SDnametoindex (sd_id, "NBFIC");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find NBFIC in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access NBFIC for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    edges[0] = 1;   /* number of lines */
    for (i = 0; i < NVIEW_ZEN_VALS; i++)
    {
        start[0] = i;   /* lines */
        status = SDreaddata (sds_id, start, NULL, edges,
            &nbfic[i*NSOLAR_ZEN_VALS]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the NBFIC SDS");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to NBFIC SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the INDTS SDS */
    sds_index = SDnametoindex (sd_id, "INDTS");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find INDTS in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access INDTS for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    start[0] = 0;   /* lines */
    start[1] = 0;   /* samples */
    edges[0] = 20;  /* number of lines */
    edges[1] = 22;  /* number of samples */
    status = SDreaddata (sds_id, start, NULL, edges, indts);
    if (status == -1)
    {
        sprintf (errmsg, "Reading data from the INDTS SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to INDTS SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the TTS SDS */
    sds_index = SDnametoindex (sd_id, "TTS");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find TTS in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access TTS for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    start[0] = 0;   /* lines */
    start[1] = 0;   /* samples */
    edges[0] = 20;  /* number of lines */
    edges[1] = 22;  /* number of samples */
    status = SDreaddata (sds_id, start, NULL, edges, tts);
    if (status == -1)
    {
        sprintf (errmsg, "Reading data from the TTS SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to TTS SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the HDF file */
    status = SDend (sd_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to HDF file: %s", anglehdf);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Allocate space for rolut [NSOLAR_VALS * NAOT_VALS * NPRES_VALS] */
    rolut = calloc (NAOTxNSOLAR_VALS * NPRES_VALS, sizeof (float));
    if (rolut == NULL)
    {
        sprintf (errmsg, "Error allocating memory for rolut");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Begin read look up table (intrinsic reflectance) */
    /* Open as HDF file for reading */
    sd_id = SDstart (intrefnm, DFACC_RDONLY);
    if (sd_id < 0)
    {
        sprintf (errmsg, "Unable to open %s for reading as SDS", intrefnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through all the bands in the input HDF file */
    start[0] = 0;  /* left-most dimension */
    start[1] = 0;
    start[2] = 0;  /* right-most dimension */
    edges[0] = NSOLAR_VALS;
    edges[1] = NAOT_VALS;
    edges[2] = NPRES_VALS;
    ibndx = -1;
    for (iband = 0; iband < nsr_bands; iband++)
    {
        /* Get the sds name and band index of this band in the rolut array */
        if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
        {
            sprintf (fname, "NRLUT_BAND_%d", iband+1);
            ibndx = iband;
        }
        else if (sat == SAT_SENTINEL_2)
        {
#ifndef PROC_ALL_BANDS
           /* Sentinel-2 we are skipping the processing of bands 9 and 10,
              but otherwise writing the other bands to the rolutt array in
              the same order. */
            if (iband == SBAND9 || iband == SBAND10)
                continue;
#endif

            /* Read this band */
            ibndx++;
            sprintf (fname, "NRLUT_BAND_%s", SENTINEL_FULL_BANDNAME[iband]);
        }

        /* Find the SDS */
        sds_index = SDnametoindex (sd_id, fname);
        if (sds_index == -1)
        {
            sprintf (errmsg, "Unable to find %s in the %s HDF file", fname,
                intrefnm);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    
        /* Open the current band as an SDS */
        sds_id = SDselect (sd_id, sds_index);
        if (sds_id < 0)
        {
            sprintf (errmsg, "Unable to access %s for reading", fname);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    
        /* Read the whole band, as-is, then rearrange the order later */
        status = SDreaddata (sds_id, start, NULL, edges, rolut);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the %s SDS", fname);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    
        /* Close the HDF SDS */
        status = SDendaccess (sds_id);
        if (status == -1)
        {
            sprintf (errmsg, "Ending access to %s SDS", fname);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Flip the look-up table value to allow the index values to be in
           the fastest-changing dimension and the band values to be in the
           slowest-changing dimension.  This allows the access for the data
           in the application to be most efficient.  rolut is read from the
           HDF file (per band) as 8000 x 22 x 7. */
        iband_indx = ibndx * NPRES_VALS * NAOTxNSOLAR_VALS;
        for (ipres = 0; ipres < NPRES_VALS; ipres++)
        {
            ipres_indx = ipres * NAOTxNSOLAR_VALS;
            for (itau = 0; itau < NAOT_VALS; itau++)
            {
                itau_indx = itau * NSOLAR_VALS;
                curr_indx = iband_indx + ipres_indx + itau_indx;
                tmp_lut_indx = itau*NPRES_VALS + ipres;
                for (ival = 0; ival < NSOLAR_VALS;
                     ival++, curr_indx++, tmp_lut_indx += NAOT_VALS*NPRES_VALS)
                {
                    rolutt[curr_indx] = rolut[tmp_lut_indx];
                }
            }
        }
    }  /* for iband */

    /* Free the temporary rolut array */
    free (rolut);

    /* Close the HDF file */
    status = SDend (sd_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to HDF file: %s", intrefnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Begin read look up table (transmission) */
    fp = fopen (transmnm, "r");
    if (fp == NULL)
    {
        sprintf (errmsg, "Opening transmission coefficient file: %s", transmnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through bands of data */
    ibndx = -1;
    for (iband = 0; iband < nsr_bands; iband++)
    {
        /* Get the band index of this band in the transt array */
        if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
            ibndx = iband;
        else if (sat == SAT_SENTINEL_2)
        {
#ifndef PROC_ALL_BANDS
           /* Sentinel-2 we are skipping the processing of bands 9 and 10,
              but otherwise writing the other bands to the transt array in
              the same order. */
            if (iband == SBAND9 || iband == SBAND10)
            {
                /* Given that this is an ASCII file, we need to consume the
                   data for this band even though we aren't using it. Each
                   pressure level has a line for each sunangle. There is a
                   single description line for each band. */
                for (iline = 0; iline < NPRES_VALS * NSUNANGLE_VALS + 1;
                     iline++)
                {
                    /* Consume this line and ignore the data */
                    if (fgets (tmpstr, sizeof (tmpstr), fp) == NULL)
                    {
                        sprintf (errmsg, "Skipping band %s in transmission "
                            "coefficient file: %s",
                            SENTINEL_FULL_BANDNAME[iband], transmnm);
                        error_handler (true, FUNC_NAME, errmsg);
                        return (ERROR);
                    }
                }

                /* Skip to the next band */
                continue;
            }
#endif

            /* Read this band */
            ibndx++;
        }

        /* This first read contains information about the band and source of
           the data; ignore for now */
        if (fgets (tmpstr, sizeof (tmpstr), fp) == NULL)
        {
            sprintf (errmsg, "Skipping data source in transmission data file");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* 7 pressure levels (1050.0 mb, 1013.0 mb, 900.0 mb, 800.0 mb,
           700.0, 600.0 mb, 500.0 mb) */
        iband_indx = ibndx * NPRES_VALS * NAOTxNSUNANGLE_VALS;
        for (ipres = 0; ipres < NPRES_VALS; ipres++)
        {
            /* This next read contains information about the pressure level of
               the data; ignore for now */
            if (fgets (tmpstr, sizeof (tmpstr), fp) == NULL)
            {
                sprintf (errmsg, "Skipping pressure level in transmission data "
                    "file");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* 21 lines of sun angles per pressure level */
            ipres_indx = ipres * NAOTxNSUNANGLE_VALS;
            for (i = 0; i < NSUNANGLE_VALS-1; i++)
            {
                /* Grab the first value in the line.  Basically this is a
                   repeat of the previous pressure level and band, as all the
                   pressure levels have the same values for the ttsr. */
                if (fscanf (fp, "%f ", &ttsr[i]) != 1)
                {
                    sprintf (errmsg, "Reading first transmission value from "
                        "transmission coefficient file: %s", transmnm);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                if (fabsf (tts[i] - ttsr[i]) > 1.0E-5)
                {
                    sprintf (errmsg, "Problem with transmission LUT: %s",
                        transmnm);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                /* Grab the remaining 22 values in the line.  Store the iaot
                   in the more efficient order for processing, not necessarily
                   for reading. */
                curr_indx = iband_indx + ipres_indx + i;
                for (iaot = 0; iaot < NAOT_VALS;
                     iaot++, curr_indx += NSUNANGLE_VALS)
                {
                    if (fscanf (fp, "%f", &transt[curr_indx]) != 1)
                    {
                        sprintf (errmsg, "Reading transmission values from "
                            "transmission coefficient file: %s", transmnm);
                        error_handler (true, FUNC_NAME, errmsg);
                        return (ERROR);
                    }
                }
            }  /* for i */

            /* Clear out the EOL for the last line */
            if (fgets (tmpstr, sizeof (tmpstr), fp) == NULL)
            {
                sprintf (errmsg, "Skipping EOL in last line");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }  /* for ipres */
    }  /* for iband */

    /* Close transmission file */
    fclose (fp);

    /* Begin read look up table (spherical albedo) */
    fp = fopen (spheranm, "r");
    if (fp == NULL)
    {
        sprintf (errmsg, "Opening spherical albedo coefficient file: %s",
            spheranm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through bands of data */
    ibndx = -1;
    for (iband = 0; iband < nsr_bands; iband++)
    {
        /* Get the band index of this band in the sphalbt/normext array */
        if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
            ibndx = iband;
        else if (sat == SAT_SENTINEL_2)
        {
#ifndef PROC_ALL_BANDS
           /* Sentinel-2 we are skipping the processing of bands 9 and 10,
              but otherwise writing the other bands to the transt array in
              the same order. */
            if (iband == SBAND9 || iband == SBAND10)
            {
               /* Given that this is an ASCII file, we need to consume the
                  data for this band even though we aren't using it. Each
                  pressure level has a line for each AOT. There is a
                  single description line for each band. */
                for (iline = 0; iline < NPRES_VALS * (NAOT_VALS+1) + 1; iline++)
                {
                    /* Consume this line and ignore the data */
                    if (fgets (tmpstr, sizeof (tmpstr), fp) == NULL)
                    {
                        sprintf (errmsg, "Skipping band %s in spherical albedo "
                            "file: %s", SENTINEL_FULL_BANDNAME[iband],
                            spheranm);
                        error_handler (true, FUNC_NAME, errmsg);
                        return (ERROR);
                    }
                }

                /* Skip to the next band */
                continue;
            }
#endif

            /* Read this band */
            ibndx++;
        }

        /* This first read contains information about the source of the data;
           ignore for now */
        if (fgets (tmpstr, sizeof (tmpstr), fp) == NULL)
        {
            sprintf (errmsg, "Skipping data source in spherical albedo data "
                "file");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* 7 pressure levels (1050.0 mb, 1013.0 mb, 900.0 mb, 800.0 mb,
           700.0, 600.0 mb, 500.0 mb) */
        iband_indx = ibndx * NPRES_VALS * NAOT_VALS;
        for (ipres = 0; ipres < NPRES_VALS; ipres++)
        {
            /* This next read contains information about the pressure level of
               the data; ignore for now */
            if (fgets (tmpstr, sizeof (tmpstr), fp) == NULL)
            {
                sprintf (errmsg, "Skipping pressure level in spherical albedo "
                    "data file");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* 22 lines of spherical albedo information */
            ipres_indx = ipres * NAOT_VALS;
            curr_indx = iband_indx + ipres_indx;
            for (iaot = 0; iaot < NAOT_VALS; iaot++, curr_indx++)
            {
                if (fscanf (fp, "%f %f %f\n", &xx, &sphalbt[curr_indx],
                    &normext[curr_indx]) != 3)
                {
                    sprintf (errmsg, "Reading spherical albedo values from "
                        "spherical albedo coefficient file: %s", spheranm);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
            }
        }
    }

    /* Close spherical albedo file */
    fclose (fp);

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  memory_allocation_main

PURPOSE:  Allocates memory for all the various arrays within the Landsat surface
reflectance application for the main application.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred allocating memory
SUCCESS        Successful completion

NOTES:
  1. Memory is allocated for each of the input variables, so it is up to the
     calling routine to free this memory.
  2. Each array passed into this function is passed in as the address to that
     1D, 2D, nD array.
******************************************************************************/
int memory_allocation_main
(
    Sat_t sat,           /* I: satellite */
    int nlines,          /* I: number of lines in the scene */
    int nsamps,          /* I: number of samples in the scene */
    bool use_orig_aero,  /* I: use the original aerosol handling if specified,
                               o/w use the semi-empirical approach */
    int16 **sza,         /* O: solar zenith angle, nlines x nsamps  */
    int16 **saa,         /* O: solar azimuth angle table, nlines x nsamps */
    int16 **vza,         /* O: view zenith angle, nlines x nsamps  */
    int16 **vaa,         /* O: view azimuth angle table, nlines x nsamps */
    uint16 **qaband,     /* O: QA band for the input image, nlines x nsamps */
    uint16 **out_band,   /* O: scaled output, nlines x nsamps */
    float ***sband,      /* O: unscaled surface reflectance and brightness temp
                               bands, nlines x nsamps */
    float ***toaband     /* O: Sentinel unscaled TOA reflectance bands,
                               nlines x nsamps */
)
{
    char FUNC_NAME[] = "memory_allocation_main"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    int i;                   /* looping variables */
    int nband_ttl =  0;      /* total number of output bands */

    /* Solar zenith array is only used for Landsat */
    if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
    {
        *sza = calloc (nlines*nsamps, sizeof (int16));
        if (*sza == NULL)
        {
            sprintf (errmsg, "Error allocating memory for sza");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Allocate the rest of the per-pixel angles only if we are using
           the original aerosol application */
        if (use_orig_aero)
        {
            *saa = calloc (nlines*nsamps, sizeof (int16));
            if (*saa == NULL)
            {
                sprintf (errmsg, "Error allocating memory for saa");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            *vza = calloc (nlines*nsamps, sizeof (int16));
            if (*vza == NULL)
            {
                sprintf (errmsg, "Error allocating memory for vza");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            *vaa = calloc (nlines*nsamps, sizeof (int16));
            if (*vaa == NULL)
            {
                sprintf (errmsg, "Error allocating memory for vaa");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }

        nband_ttl = NBANDL_TTL_OUT;
    }
    else if (sat == SAT_SENTINEL_2)
    {
        nband_ttl = NBANDS_TTL_OUT;
        *toaband = calloc (nband_ttl-1, sizeof (float*));
        if (*toaband == NULL)
        {
            sprintf (errmsg, "Error allocating memory for toaband");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        for (i = 0; i < nband_ttl-1; i++)
        {
            (*toaband)[i] = calloc (nlines*nsamps, sizeof (float));
            if ((*toaband)[i] == NULL)
            {
                sprintf (errmsg, "Error allocating memory for toaband");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }
    }

    *qaband = calloc (nlines*nsamps, sizeof (uint16));
    if (*qaband == NULL)
    {
        sprintf (errmsg, "Error allocating memory for qaband");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *out_band = calloc (nlines*nsamps, sizeof (uint16));
    if (*out_band == NULL)
    {
        sprintf (errmsg, "Error allocating memory for out_band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Given that the QA band is its own separate array of uint16s, we need
       one less band for the output image data */
    *sband = calloc (nband_ttl-1, sizeof (float*));
    if (*sband == NULL)
    {
        sprintf (errmsg, "Error allocating memory for sband");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    for (i = 0; i < nband_ttl-1; i++)
    {
        (*sband)[i] = calloc (nlines*nsamps, sizeof (float));
        if ((*sband)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for sband");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  landsat_memory_allocation_sr

PURPOSE:  Allocates memory for all the various arrays needed specifically for
the Landsat surface reflectance corrections.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred allocating memory
SUCCESS        Successful completion

NOTES:
  1. Memory is allocated for each of the input variables, so it is up to the
     calling routine to free this memory.
  2. Each array passed into this function is passed in as the address to that
     1D, 2D, nD array.
******************************************************************************/
int landsat_memory_allocation_sr
(
    int nlines,          /* I: number of lines in the scene */
    int nsamps,          /* I: number of samples in the scene */
    float **aerob1,      /* O: atmospherically corrected band 1 data
                              (unscaled TOA refl), nlines x nsamps */
    float **aerob2,      /* O: atmospherically corrected band 2 data
                              (unscaled TOA refl), nlines x nsamps */
    float **aerob4,      /* O: atmospherically corrected band 4 data
                              (unscaled TOA refl), nlines x nsamps */
    float **aerob5,      /* O: atmospherically corrected band 5 data
                              (unscaled TOA refl), nlines x nsamps */
    float **aerob7,      /* O: atmospherically corrected band 7 data
                              (unscaled TOA refl), nlines x nsamps */
    uint8 **ipflag,      /* O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps */
    float **twvi,        /* O: interpolated water vapor value,
                               nlines x nsamps */
    float **tozi,        /* O: interpolated ozone value, nlines x nsamps */
    float **tp,          /* O: interpolated pressure value, nlines x nsamps */
    float **taero,       /* O: aerosol values for each pixel, nlines x nsamps */
    float **teps,        /* O: eps (angstrom coefficient) for each pixel,
                               nlines x nsamps*/
    int16 **dem,         /* O: CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 **andwi,       /* O: avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 **sndwi,       /* O: standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 **ratiob1,     /* O: mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **ratiob2,     /* O: mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **ratiob7,     /* O: mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **intratiob1,  /* O: band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **intratiob2,  /* O: band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **intratiob7,  /* O: band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **slpratiob1,  /* O: slope band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **slpratiob2,  /* O: slope band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **slpratiob7,  /* O: slope band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    uint16 **wv,         /* O: water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 **oz,          /* O: ozone values [CMG_NBLAT x CMG_NBLON] */
    float **rolutt,      /* O: intrinsic reflectance table
                         [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float **transt,      /* O: transmission table
                        [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUNANGLE_VALS] */
    float **sphalbt,     /* O: spherical albedo table
                               [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float **normext,     /* O: aerosol extinction coefficient at the current
                               wavelength (normalized at 550nm) 
                               [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float **tsmax,       /* O: maximum scattering angle table
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **tsmin,       /* O: minimum scattering angle table
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **nbfic,       /* O: communitive number of azimuth angles
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **nbfi,        /* O: number of azimuth angles
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **ttv          /* O: view angle table
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
)
{
    char FUNC_NAME[] = "landsat_memory_allocation_sr"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    int nsr_bands = 0;       /* number of SR bands */

    /* Setup Landsat number of SR bands */
    nsr_bands = NSRL_BANDS;

    *aerob1 = calloc (nlines*nsamps, sizeof (float));
    if (*aerob1 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob1");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *aerob2 = calloc (nlines*nsamps, sizeof (float));
    if (*aerob2 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob2");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *aerob4 = calloc (nlines*nsamps, sizeof (float));
    if (*aerob4 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob4");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *aerob5 = calloc (nlines*nsamps, sizeof (float));
    if (*aerob5 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob5");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *aerob7 = calloc (nlines*nsamps, sizeof (float));
    if (*aerob7 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob7");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *twvi = calloc (nlines*nsamps, sizeof (float));
    if (*twvi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for twvi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *tozi = calloc (nlines*nsamps, sizeof (float));
    if (*tozi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tozi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *tp = calloc (nlines*nsamps, sizeof (float));
    if (*tp == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tp");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *taero = calloc (nlines*nsamps, sizeof (float));
    if (*taero == NULL)
    {
        sprintf (errmsg, "Error allocating memory for taero");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *teps = calloc (nlines*nsamps, sizeof (float));
    if (*teps == NULL)
    {
        sprintf (errmsg, "Error allocating memory for teps");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ipflag = calloc (nlines*nsamps, sizeof (uint8));
    if (*ipflag == NULL)
    {
        sprintf (errmsg, "Error allocating memory for ipflag");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Allocate memory for all the climate modeling grid files */
    *dem = calloc (DEM_NBLAT * DEM_NBLON, sizeof (int16*));
    if (*dem == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the DEM");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *andwi = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*andwi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the andwi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *sndwi = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*sndwi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the sndwi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ratiob1 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*ratiob1 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the ratiob1");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ratiob2 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*ratiob2 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the ratiob2");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ratiob7 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*ratiob7 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the ratiob7");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *intratiob1 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*intratiob1 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the intratiob1");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *intratiob2 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*intratiob2 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the intratiob2");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *intratiob7 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*intratiob7 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the intratiob7");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *slpratiob1 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*slpratiob1 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the slpratiob1");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *slpratiob2 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*slpratiob2 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the slpratiob2");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *slpratiob7 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*slpratiob7 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the slpratiob7");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *wv = calloc (CMG_NBLAT * CMG_NBLON, sizeof (int16));
    if (*wv == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the wv");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *oz = calloc (CMG_NBLAT * CMG_NBLON, sizeof (uint8));
    if (*oz == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the oz");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* rolutt, transt, sphalbt, and normext */
    *rolutt = calloc (nsr_bands*NPRES_VALS*NAOT_VALS*NSOLAR_VALS,
        sizeof (float));
    if (*rolutt == NULL)
    {
        sprintf (errmsg, "Error allocating memory for rolutt");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *transt = calloc (nsr_bands*NPRES_VALS*NAOT_VALS*NSUNANGLE_VALS,
        sizeof (float));
    if (*transt == NULL)
    {
        sprintf (errmsg, "Error allocating memory for transt");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *sphalbt = calloc (nsr_bands*NPRES_VALS*NAOT_VALS, sizeof (float));
    if (*sphalbt == NULL)
    {
        sprintf (errmsg, "Error allocating memory for sphalbt");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *normext = calloc (nsr_bands*NPRES_VALS*NAOT_VALS, sizeof (float));
    if (*normext == NULL)
    {
        sprintf (errmsg, "Error allocating memory for normext");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* float tsmax, tsmin, nbfic, nbfi, and ttv */
    *tsmax = calloc (NVIEW_ZEN_VALS*NSOLAR_ZEN_VALS, sizeof (float));
    if (*tsmax == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tsmax");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *tsmin = calloc (NVIEW_ZEN_VALS*NSOLAR_ZEN_VALS, sizeof (float));
    if (*tsmin == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tsmin");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *nbfic = calloc (NVIEW_ZEN_VALS*NSOLAR_ZEN_VALS, sizeof (float));
    if (*nbfic == NULL)
    {
        sprintf (errmsg, "Error allocating memory for nbfic");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *nbfi = calloc (NVIEW_ZEN_VALS*NSOLAR_ZEN_VALS, sizeof (float));
    if (*nbfi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for nbfi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ttv = calloc (NVIEW_ZEN_VALS*NSOLAR_ZEN_VALS, sizeof (float));
    if (*ttv == NULL)
    {
        sprintf (errmsg, "Error allocating memory for ttv");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  sentinel_memory_allocation_sr

PURPOSE:  Allocates memory for all the various arrays needed specifically for
the Sentinel surface reflectance corrections.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred allocating memory
SUCCESS        Successful completion

NOTES:
  1. Memory is allocated for each of the input variables, so it is up to the
     calling routine to free this memory.
  2. Each array passed into this function is passed in as the address to that
     1D, 2D, nD array.
******************************************************************************/
int sentinel_memory_allocation_sr
(
    int nlines,          /* I: number of lines in the scene */
    int nsamps,          /* I: number of samples in the scene */
    uint8 **ipflag,      /* O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps */
    float **twvi,        /* O: interpolated water vapor value,
                               nlines x nsamps */
    float **tozi,        /* O: interpolated ozone value, nlines x nsamps */
    float **tp,          /* O: interpolated pressure value, nlines x nsamps */
    float **taero,       /* O: aerosol values for each pixel, nlines x nsamps */
    float **teps,        /* O: eps (angstrom coefficient) for each pixel,
                               nlines x nsamps*/
    int16 **dem,         /* O: CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 **andwi,       /* O: avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 **sndwi,       /* O: standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 **ratiob1,     /* O: mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **ratiob2,     /* O: mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **ratiob7,     /* O: mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **intratiob1,  /* O: band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **intratiob2,  /* O: band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **intratiob7,  /* O: band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **slpratiob1,  /* O: slope band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **slpratiob2,  /* O: slope band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 **slpratiob7,  /* O: slope band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    uint16 **wv,         /* O: water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 **oz,          /* O: ozone values [CMG_NBLAT x CMG_NBLON] */
    float **rolutt,      /* O: intrinsic reflectance table
                         [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float **transt,      /* O: transmission table
                        [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUNANGLE_VALS] */
    float **sphalbt,     /* O: spherical albedo table
                               [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float **normext,     /* O: aerosol extinction coefficient at the current
                               wavelength (normalized at 550nm) 
                               [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float **tsmax,       /* O: maximum scattering angle table
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **tsmin,       /* O: minimum scattering angle table
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **nbfic,       /* O: communitive number of azimuth angles
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **nbfi,        /* O: number of azimuth angles
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float **ttv          /* O: view angle table
                               [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
)
{
    char FUNC_NAME[] = "sentinel_memory_allocation_sr"; /* function name */
    char errmsg[STR_SIZE];         /* error message */
    int nsr_bands = NSRS_BANDS;  /* number of SR bands */

    /* Setup Sentinel number of SR bands */
    nsr_bands = NSRS_BANDS;

    /* Allocate memory for wv, oz, pres, aero, eps, and ipflag */
    *twvi = calloc (nlines*nsamps, sizeof (float));
    if (*twvi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for twvi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *tozi = calloc (nlines*nsamps, sizeof (float));
    if (*tozi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tozi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *tp = calloc (nlines*nsamps, sizeof (float));
    if (*tp == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tp");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *taero = calloc (nlines*nsamps, sizeof (float));
    if (*taero == NULL)
    {
        sprintf (errmsg, "Error allocating memory for taero");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *teps = calloc (nlines*nsamps, sizeof (float));
    if (*teps == NULL)
    {
        sprintf (errmsg, "Error allocating memory for teps");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ipflag = calloc (nlines*nsamps, sizeof (uint8));
    if (*ipflag == NULL)
    {
        sprintf (errmsg, "Error allocating memory for ipflag");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Allocate memory for all the climate modeling grid files */
    *dem = calloc (DEM_NBLAT * DEM_NBLON, sizeof (int16*));
    if (*dem == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the DEM");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *andwi = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*andwi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the andwi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *sndwi = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*sndwi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the sndwi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ratiob1 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*ratiob1 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the ratiob1");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ratiob2 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*ratiob2 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the ratiob2");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ratiob7 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*ratiob7 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the ratiob7");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *intratiob1 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*intratiob1 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the intratiob1");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *intratiob2 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*intratiob2 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the intratiob2");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *intratiob7 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*intratiob7 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the intratiob7");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *slpratiob1 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*slpratiob1 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the slpratiob1");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *slpratiob2 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*slpratiob2 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the slpratiob2");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *slpratiob7 = calloc (RATIO_NBLAT * RATIO_NBLON, sizeof (int16));
    if (*slpratiob7 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the slpratiob7");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *wv = calloc (CMG_NBLAT * CMG_NBLON, sizeof (int16));
    if (*wv == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the wv");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *oz = calloc (CMG_NBLAT * CMG_NBLON, sizeof (uint8));
    if (*oz == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the oz");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* rolutt, transt, sphalbt, and normext */
    *rolutt = calloc (nsr_bands*NPRES_VALS*NAOT_VALS*NSOLAR_VALS,
        sizeof (float));
    if (*rolutt == NULL)
    {
        sprintf (errmsg, "Error allocating memory for rolutt");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *transt = calloc (nsr_bands*NPRES_VALS*NAOT_VALS*NSUNANGLE_VALS,
        sizeof (float));
    if (*transt == NULL)
    {
        sprintf (errmsg, "Error allocating memory for transt");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *sphalbt = calloc (nsr_bands*NPRES_VALS*NAOT_VALS, sizeof (float));
    if (*sphalbt == NULL)
    {
        sprintf (errmsg, "Error allocating memory for sphalbt");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *normext = calloc (nsr_bands*NPRES_VALS*NAOT_VALS, sizeof (float));
    if (*normext == NULL)
    {
        sprintf (errmsg, "Error allocating memory for normext");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* float tsmax, tsmin, nbfic, nbfi, and ttv */
    *tsmax = calloc (NVIEW_ZEN_VALS*NSOLAR_ZEN_VALS, sizeof (float));
    if (*tsmax == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tsmax");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *tsmin = calloc (NVIEW_ZEN_VALS*NSOLAR_ZEN_VALS, sizeof (float));
    if (*tsmin == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tsmin");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *nbfic = calloc (NVIEW_ZEN_VALS*NSOLAR_ZEN_VALS, sizeof (float));
    if (*nbfic == NULL)
    {
        sprintf (errmsg, "Error allocating memory for nbfic");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *nbfi = calloc (NVIEW_ZEN_VALS*NSOLAR_ZEN_VALS, sizeof (float));
    if (*nbfi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for nbfi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ttv = calloc (NVIEW_ZEN_VALS*NSOLAR_ZEN_VALS, sizeof (float));
    if (*ttv == NULL)
    {
        sprintf (errmsg, "Error allocating memory for ttv");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  read_auxiliary_files

PURPOSE:  Reads the auxiliary files required for this application.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading one of the auxiliary files
SUCCESS        Successful completion

NOTES:
  1. It is assumed that memory has already been allocated for the input data
     arrays.
******************************************************************************/
int read_auxiliary_files
(
    char *cmgdemnm,     /* I: climate modeling grid DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm,        /* I: auxiliary filename for ozone and water vapor */
    int16 *dem,         /* O: CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 *andwi,       /* O: avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *sndwi,       /* O: standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob1,     /* O: mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob2,     /* O: mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob7,     /* O: mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob1,  /* O: band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob2,  /* O: band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob7,  /* O: band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob1,  /* O: slope band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob2,  /* O: slope band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob7,  /* O: slope band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    uint16 *wv,         /* O: water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 *oz           /* O: ozone values [CMG_NBLAT x CMG_NBLON] */
)
{
    char FUNC_NAME[] = "read_auxiliary_files"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char sds_name[STR_SIZE]; /* name of the SDS being read */
    int i;               /* looping variable */
    int status;          /* return status of the HDF function */
    int start[5];        /* starting point to read SDS data; handles up to
                            4D dataset */
    int edges[5];        /* number of values to read in SDS data; handles up to
                            4D dataset */
    int sd_id;           /* file ID for the HDF file */
    int sds_id;          /* ID for the current SDS */
    int sds_index;       /* index for the current SDS */

    /*** Read the DEM ***/
    sd_id = SDstart (cmgdemnm, DFACC_RDONLY);
    if (sd_id < 0)
    {
        sprintf (errmsg, "Unable to open %s for reading as SDS", cmgdemnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name */
    strcpy (sds_name, "averaged elevation");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the DEM file %s", sds_name,
            cmgdemnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < DEM_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = DEM_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, &dem[i * DEM_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the DEM file */
    status = SDend (sd_id);
    if (status != 0)
    {
        sprintf (errmsg, "Closing DEM file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /*** Read the RATIO file ***/
    sd_id = SDstart (rationm, DFACC_RDONLY);
    if (sd_id < 0)
    {
        sprintf (errmsg, "Unable to open %s for reading as SDS", rationm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 6) */
    strcpy (sds_name, "average ndvi");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges,
            &andwi[i * RATIO_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 3) */
    strcpy (sds_name, "average ratio b10");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges,
            &ratiob2[i * RATIO_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 2) */
    strcpy (sds_name, "average ratio b9");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges,
            &ratiob1[i * RATIO_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 4) */
    strcpy (sds_name, "average ratio b7");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges,
            &ratiob7[i * RATIO_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 14) */
    strcpy (sds_name, "standard ndvi");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges,
            &sndwi[i * RATIO_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 21) */
    strcpy (sds_name, "slope ratiob9");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges,
            &slpratiob1[i * RATIO_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 22) */
    strcpy (sds_name, "inter ratiob9");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges,
            &intratiob1[i * RATIO_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 24) */
    strcpy (sds_name, "slope ratiob10");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges,
            &slpratiob2[i * RATIO_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 25) */
    strcpy (sds_name, "inter ratiob10");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges,
            &intratiob2[i * RATIO_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 27) */
    strcpy (sds_name, "slope ratiob7");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges,
            &slpratiob7[i * RATIO_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 28) */
    strcpy (sds_name, "inter ratiob7");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges,
            &intratiob7[i * RATIO_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the RATIO file */
    status = SDend (sd_id);
    if (status != 0)
    {
        sprintf (errmsg, "Closing RATIO file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read ozone and water vapor from the user-specified auxiliary file */
    sd_id = SDstart (auxnm, DFACC_RDONLY);
    if (sd_id < 0)
    {
        sprintf (errmsg, "Unable to open %s for reading as SDS", auxnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name */
    strcpy (sds_name, "Coarse Resolution Ozone");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the AUX file %s", sds_name,
            auxnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < CMG_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = CMG_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, &oz[i * CMG_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name */
    strcpy (sds_name, "Coarse Resolution Water Vapor");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the AUX file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < CMG_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = CMG_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, &wv[i * CMG_NBLON]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the AUX file */
    status = SDend (sd_id);
    if (status != 0)
    {
        sprintf (errmsg, "Closing AUX file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Successful completion */
    return (SUCCESS);
}

