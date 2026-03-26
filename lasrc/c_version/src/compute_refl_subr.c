/******************************************************************************
FILE: compute_refl_subr.c

PURPOSE: Contains functions for handling the Landsat TOA reflectance and
Landsat/Sentinel surface reflectance corrections.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
******************************************************************************/
//#define USE_GCTP 1
/* TODO uncomment to use the GCTP library */

#include "time.h"
#include "aero_interp.h"
#include "poly_coeff.h"
#include "read_level1_qa.h"
#ifndef USE_GCTP
#include "utmtodeg.h"
#endif


/******************************************************************************
MODULE:  init_sr_refl

PURPOSE:  Initialization for the atmospheric corrections.  Initialization for
look up tables, auxiliary data, mapping, and geolocation information is used
for the surface reflectance correction.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error initializing the atmospheric parameters
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
1. The view angle is set to 0.0 and this never changes.
2. The DEM is used to calculate the surface pressure.
******************************************************************************/
int init_sr_refl
(
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    Input_t *input,     /* I: input structure for the Landsat product */
    Space_def_t *space_def,  /* I: space definition structure */
    Geoloc_t *space,    /* I: structure for geolocation information */
    char *anglehdf,     /* I: angle HDF filename */
    char *intrefnm,     /* I: intrinsic reflectance filename */
    char *transmnm,     /* I: transmission filename */
    char *spheranm,     /* I: spherical albedo filename */
    char *cmgdemnm,     /* I: climate modeling grid DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm,        /* I: auxiliary filename for ozone and water vapor */
    char *aux2nm,       /* I: 2nd auxiliary filename for ozone and water vapor (23-FEB-26, JPR) */
    aux_src_t aux_src,  /* I: Identifies the source of atmospheric aux data */
    float *eps,         /* O: angstrom coefficient */
    int *iaots,         /* O: index for AOTs */
    float *xtv,         /* O: observation zenith angle (deg) */
    float *xmuv,        /* O: cosine of observation zenith angle */
    float *xfi,         /* O: azimuthal difference between sun and
                              observation (deg) */
    float *cosxfi,      /* O: cosine of azimuthal difference */
    float *raot550nm,   /* O: nearest value of AOT */
    float *pres,        /* O: surface pressure */
    float *uoz,         /* O: total column ozone */
    float *uwv,         /* O: total column water vapor (precipital water
                              vapor) */
    float *xtsstep,     /* O: solar zenith step value */
    float *xtsmin,      /* O: minimum solar zenith value */
    float *xtvstep,     /* O: observation step value */
    float *xtvmin,      /* O: minimum observation value */
    float *tsmax,       /* O: maximum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *tsmin,       /* O: minimum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tts[22],      /* O: sun angle table */
    float *ttv,         /* O: view angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    int32 indts[22],    /* O: index for the sun angle table */
    float *rolutt,      /* O: intrinsic reflectance table
                          [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float *transt,      /* O: transmission table 
                      [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUN_ANGLE_VALS] */
    float *sphalbt,     /* O: spherical albedo table 
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext,     /* O: aerosol extinction coefficient at the current
                              wavelength (normalized at 550nm)
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *nbfic,       /* O: communitive number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfi,        /* O: number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    int16 *dem,         /* O: CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 *andwi,       /* O: avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *sndwi,       /* O: standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob1,     /* O: mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob2,     /* O: mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob7,     /* O: mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob1,  /* O: integer band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob2,  /* O: integer band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob7,  /* O: integer band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob1,  /* O: slope band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob2,  /* O: slope band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob7,  /* O: slope band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    uint16 *wv,         /* O: water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 *oz           /* O: ozone values [CMG_NBLAT x CMG_NBLON] */
)
{
    char errmsg[STR_SIZE];                   /* error message */
    char FUNC_NAME[] = "init_sr_refl";       /* function name */
    int retval;          /* return status */
    int lcmg, scmg;      /* line/sample index for the CMG */
    int cmg_pix;         /* pixel location in the CMG array for [lcmg][scmg] */
    int dem_pix;         /* pixel location in the DEM array for [lcmg][scmg] */
    float xcmg, ycmg;    /* x/y location for CMG */
    Sat_t sat = input->meta.sat; /* satellite */

    /* Vars for forward/inverse mapping space */
#ifdef USE_GCTP
    Img_coord_float_t img;        /* coordinate in line/sample space */
    Geo_coord_t geo;              /* coordinate in lat/long space */
#endif
    float center_lat, center_lon; /* lat/long for scene center */

    /* Initialize the view azimuth and zenith values */
    if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
    {
        /* Landsat values for view zenith and azimuth are 0.0 */
        *xtv = 0.0;
        *xmuv = cos (*xtv * DEG2RAD);
        *xfi = 0.0;
        *cosxfi = cos (*xfi * DEG2RAD);
    }
    else if (sat == SAT_SENTINEL_2)
    {
        /* Sentinel values for view zenith and azimuth are in the input
           metadata */
        *xtv = input->meta.view_zen;  /* degrees */
        *xmuv = cos (*xtv * DEG2RAD);
        *xfi = acos (cos ((input->meta.sun_az - input->meta.view_az)*DEG2RAD));
        *cosxfi = cos (*xfi); /* xfi is still in radians */

        /* Convert xfi to degrees */
        *xfi = *xfi * RAD2DEG;
		
    }

    /* Initialize the look up tables */
    *eps = 1.0;
    *iaots = 0;
    *xtsmin = 0;
    *xtsstep = 4.0;
    *xtvmin = 2.84090;
    *xtvstep = 6.52107 - *xtvmin;
    retval = readluts (sat, tsmax, tsmin, ttv, tts, nbfic, nbfi, indts, rolutt,
        transt, sphalbt, normext, *xtsstep, *xtsmin, anglehdf, intrefnm,
        transmnm, spheranm);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Reading the LUTs");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (sat == SAT_LANDSAT_8 || sat == SAT_LANDSAT_9)
        printf ("The LUTs for urban clean case v2.0 have been read.  We can "
            "now perform atmospheric correction.\n");
    else if (sat == SAT_SENTINEL_2)
        printf ("The LUTs for urban clean case v3.0 have been read.  We can "
            "now perform atmospheric correction.\n");

    /* Read the auxiliary data files used as input to the reflectance
       calculations */
       

    /*printf("Acquisition date, time: %ld, %lf\n", input->meta.acq_date.jday2000, input->meta.acq_date.hour); */


    retval = read_auxiliary_files (cmgdemnm, rationm, auxnm, aux2nm, aux_src, dem,
        andwi, sndwi, ratiob1, ratiob2, ratiob7, intratiob1, intratiob2,
        intratiob7, slpratiob1, slpratiob2, slpratiob7, wv, oz, input->meta.acq_date);  /* for GRIB inputs */
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Reading the auxiliary files");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Getting parameters for atmospheric correction */
    /* Update to get the parameter of the scene center */
    *pres = ATMOS_PRES_0;
    *uoz = 0.30;
    *uwv = 2.5;

    /* Use scene center (and center of the pixel) to compute atmospheric
       parameters */
#ifdef USE_GCTP
    img.l = (int) (nlines * 0.5);
    img.s = (int) (nsamps * 0.5);
    img.is_fill = false;

    if (!from_space (space, &img, &geo))
    {
        sprintf (errmsg, "Mapping scene center to geolocation coords");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    center_lat = geo.lat * RAD2DEG;
    center_lon = geo.lon * RAD2DEG;
#else
    utmtodeg (space_def, (int) (nlines * 0.5), (int) (nsamps * 0.5),
        &center_lat, &center_lon);
#endif
    printf ("Scene center line/sample: %f, %f\n", nlines * 0.5, nsamps * 0.5);
    printf ("Scene center lat/long: %f, %f\n", center_lat, center_lon);

    /* Use the scene center lat/long to determine the line/sample in the
       CMG-related lookup tables, using the center of the UL pixel.
       Negative latitude values should be the largest line values in the CMG
       grid.  Negative longitude values should be the smallest sample values
       in the CMG grid. */
    ycmg = (89.975 - center_lat) * 20.0;    /* vs / 0.05 */
    xcmg = (179.975 + center_lon) * 20.0;   /* vs / 0.05 */
    lcmg = (int) roundf (ycmg);
    scmg = (int) roundf (xcmg);

    /* Handle the edges of the lat/long values in the CMG grid */
    if (lcmg < 0)
        lcmg = 0;
    else if (lcmg >= CMG_NBLAT)
        lcmg = CMG_NBLAT;

    if (scmg < 0)
        scmg = 0;
    else if (scmg >= CMG_NBLON)
        scmg = CMG_NBLON;

    cmg_pix = lcmg * CMG_NBLON + scmg;

    /* If not FILL unscale the water vapor, otherwise use the default */
    if (wv[cmg_pix] != WV_FILL[aux_src])
        *uwv = (float)wv[cmg_pix] / WV_SCALE_FACTOR[aux_src];
    else
        *uwv = WV_DEFAULT[aux_src];

    /* If not FILL unscale the ozone, otherwise use the default */
    if (oz[cmg_pix] != OZ_FILL[aux_src])
        *uoz = (float)oz[cmg_pix] / OZ_SCALE_FACTOR[aux_src];
    else
        *uoz = OZ_DEFAULT[aux_src];

    /* If not FILL compute the atmospheric pressure for the elevation,
       otherwise use the default */
    dem_pix = lcmg * DEM_NBLON + scmg;
    if (dem[dem_pix] != -9999)
        *pres = ATMOS_PRES_0 * exp (-dem[dem_pix] * ONE_DIV_8500);
    else
        *pres = ATMOS_PRES_0;
    *raot550nm = 0.05;
printf ("DEBUG: CMG pixel line,sample: %d, %d\n", lcmg, scmg);
printf ("   uwv, uoz, pres: %f, %f, %f\n", *uwv, *uoz, *pres);

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  find_closest_non_fill

PURPOSE:  Finds the closest non-fill pixel in the aerosol window

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
false           No pixel found
true            Pixel found

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
bool find_closest_non_fill
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    int nlines,        /* I: number of lines in QA band */
    int nsamps,        /* I: number of samps in QA band */
    int center_line,   /* I: line for the center of the aerosol window */
    int center_samp,   /* I: sample for the center of the aerosol window */
    int half_aero_window, /* I: size of half the aerosol window */
    int *nearest_line, /* O: line for nearest non-fill pix in aerosol window */
    int *nearest_samp  /* O: samp for nearest non-fill pix in aerosol window */
)
{
    int curr_pix;            /* looping variable for pixels */
    int line, samp;          /* looping variables for lines and samples */
    int aero_window;         /* looping variabel for the aerosol window */

    /* Loop around the center pixel, moving outward with each loop, searching
       for a pixel that is not of the QA type specified and is not fill */
    for (aero_window = 1; aero_window <= half_aero_window; aero_window++)
    {
        for (line = center_line - aero_window;
             line <= center_line + aero_window; line++)
        {
            /* Make sure the line is valid */
            if (line < 0 || line >= nlines)
                continue;

            curr_pix = line * nsamps + center_samp - aero_window;
            for (samp = center_samp - aero_window;
                 samp <= center_samp + aero_window; samp++, curr_pix++)
            {
                /* Make sure the sample is valid */
                if (samp < 0 || samp >= nsamps)
                    continue;

                /* If this pixel is not fill, then mark it as the closest
                   non-fill pixel and return */
                if (!level1_qa_is_fill (qaband[curr_pix]))
                {
                    *nearest_line = line;
                    *nearest_samp = samp;
                    return (true);
                }
            }
        }
    }

    /* No pixel was found that met the criteria */
    return (false);
}
