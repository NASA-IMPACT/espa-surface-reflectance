#include "aero_interp.h"
#include "quick_select.h"
#include "read_level1_qa.h"
#include "read_level2_qa.h"

typedef enum {
    FORWARD,
    REVERSE
} Fill_direction_t;

/******************************************************************************
MODULE:  aerosol_interp_landsat

PURPOSE:  Interpolates the Landsat aerosol values throughout the image using the
aerosols that were calculated for each NxN window. Also cleans up the fill
pixels in the ipflag.

RETURN VALUE:
Type = N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
void aerosol_interp_landsat
(
    Espa_internal_meta_t *xml_metadata, /* I: XML metadata information */
    int aero_window,   /* I: size of the aerosol window */
    int half_aero_window, /* I: size of half the aerosol window */
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps.  It is expected that the ipflag
                               values are computed for the center of the
                               aerosol windows. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                          It is expected that the aerosol values are computed
                          for the center of the aerosol windows.  This routine
                          will fill in the pixels for the remaining, non-center
                          pixels of the window. */
    int nlines,        /* I: number of lines in qaband & taero bands */
    int nsamps         /* I: number of samps in qaband & taero bands */
)
{
    int i;                 /* looping variable for the bands */
    int line, samp;        /* looping variable for lines and samples */
    int curr_pix;          /* current pixel in 1D arrays of nlines x nsamps */
    int center_line;       /* line for the center of the aerosol window */
    int center_line1;      /* line+1 for the center of the aerosol window */
    int center_lindex;     /* center line array index */
    int center_lindex1;    /* center line 1 array index */
    int center_samp;       /* sample for the center of the aerosol window */
    int center_samp1;      /* sample+1 for the center of the aerosol window */
    int refl_indx = -99;   /* index of band 1 or first band */
    int aero_pix11;        /* pixel location for aerosol window values
                              [lcmg][scmg] */
    int aero_pix12;        /* pixel location for aerosol window values
                              [lcmg][scmg2] */
    int aero_pix21;        /* pixel location for aerosol window values
                              [lcmg2][scmg] */
    int aero_pix22;        /* pixel location for aerosol window values
                              [lcmg2][scmg2] */
    long npixels;          /* number of pixels to process */
    float xaero, yaero;    /* x/y location for aerosol pixel within the overall
                              larger aerosol window grid */
    float aero11;          /* aerosol value at window line, samp */
    float aero12;          /* aerosol value at window line, samp+1 */
    float aero21;          /* aerosol value at window line+1, samp */
    float aero22;          /* aerosol value at window line+1, samp+1 */
    float u, v;            /* line, sample fractional distance from current
                              pixel (weight applied to furthest line, sample) */
    int aero_window_index_step = aero_window*nsamps; /* aerosol window array
                                                        step size */
    float aero_step = 1.0/aero_window; /* fraction of window size representing
                                          1 pixel */

    /* Use band 1 band-related metadata for the reflectance information for
       Landsat (Level 1 products). If band 1 isn't available then just use the
       first band in the XML file. */
    for (i = 0; i < xml_metadata->nbands; i++)
    {
        if (!strcmp (xml_metadata->band[i].name, "b1") &&
            !strncmp (xml_metadata->band[i].product, "L1", 2))
        {
            /* this is the index we'll use for the band info */
            refl_indx = i;
        }
    }
    if (refl_indx == -99)
        refl_indx = 1;   /* default use second band in XML file */

    /* Interpolate the aerosol data for each pixel location */
    for (line = 0, curr_pix = 0; line < nlines; line++)
    {
        /* Determine the line of the representative center pixel in the
           aerosol NxN window array */
        center_line = (int) (line * aero_step) * aero_window + half_aero_window;
        center_lindex = center_line * nsamps;

        /* Determine fractional location of this line in the aerosol window.
           Negative values are at the top of the window. */
        yaero = (line - center_line) * aero_step;
        u = yaero - (int) yaero;

        /* Determine if this pixel is closest to the line below or the line
           above. If the fractional value is in the top part of the aerosol
           window, then use the line above.  Otherwise use the line below. */
        if (u < 0.0)
        {
            center_line1 = center_line - aero_window;
            center_lindex1 = center_lindex - aero_window_index_step;

            /* If the aerosol window line value is outside the bounds of the
               scene, then just use the same line in the aerosol window */
            if (center_line1 < 0)
            {
                center_line1 = center_line;
                center_lindex1 = center_lindex;
            }
        }
        else
        {
            center_line1 = center_line + aero_window;
            center_lindex1 = center_lindex + aero_window_index_step;

            /* If the aerosol window line value is outside the bounds of the
               scene, then just use the same line in the aerosol window */
            if (center_line1 >= nlines-1)
            {
                center_line1 = center_line;
                center_lindex1 = center_lindex;
            }
        }

        /* Make the fractional line distance positive, regardless of where it
           is in the window. */
        u = fabsf (u);

        for (samp = 0; samp < nsamps; samp++, curr_pix++)
        {
            /* If this pixel is fill, then don't process */
            if (level1_qa_is_fill (qaband[curr_pix]))
                continue;

            /* Determine the sample of the representative center pixel in the
               aerosol NxN window array */
            center_samp = (int) (samp * aero_step) * aero_window +
                half_aero_window;

            /* If the current line, sample are the same as the center line,
               sample, then skip to the next pixel.  We already have the
               aerosol value. */
            if (samp == center_samp && line == center_line)
                continue;

            /* Determine fractional location of this sample in the aerosol
               window. Negative values are at the left of the window. */
            xaero = (samp - center_samp) * aero_step;
            v = xaero - (int) xaero;

            /* Determine if this pixel is closest to the sample to the left or
               the sample to the right.  If the fractional value is on the left
               side of the aerosol window, then use the sample to the left.
               Otherwise use the sample to the right. */
            if (v < 0.0)
            {
                center_samp1 = center_samp - aero_window;
    
                /* If the aerosol window sample value is outside the bounds of
                   the scene, then just use the same sample in the aerosol
                   window */
                if (center_samp1 < 0)
                    center_samp1 = center_samp;
            }
            else
            {
                center_samp1 = center_samp + aero_window;
    
                /* If the aerosol window sample value is outside the bounds of
                   the scene, then just use the same sample in the aerosol
                   window */
                if (center_samp1 >= nsamps-1)
                    center_samp1 = center_samp;
            }

            /* Determine the four aerosol window pixels to be used for
               interpolating the current pixel */
            aero_pix11 = center_lindex + center_samp;
            aero_pix12 = center_lindex + center_samp1;
            aero_pix21 = center_lindex1 + center_samp;
            aero_pix22 = center_lindex1 + center_samp1;

            /* Get the aerosol values */
            aero11 = taero[aero_pix11];
            aero12 = taero[aero_pix12];
            aero21 = taero[aero_pix21];
            aero22 = taero[aero_pix22];

            /* Make the fractional sample distance positive, regardless of
               where it is in the window. */
            v = fabsf(v);

            /* Interpolate the aerosol */
            taero[curr_pix] = aero11  +
                              u * (aero21 - aero11) +
                              v * (aero12 - aero11) +
                              u * v * (aero11 - aero12 - aero21 + aero22);

            /* Set the aerosol to window interpolated. Clear anything else. */
            ipflag[curr_pix] = (1 << IPFLAG_INTERP_WINDOW);

            /* If any of the window pixels used in the interpolation were
               water pixels, then mask this pixel with water (in addition to
               the interpolation bit already set) */
            if (lasrc_qa_is_water (ipflag[aero_pix11]) ||
                lasrc_qa_is_water (ipflag[aero_pix12]) ||
                lasrc_qa_is_water (ipflag[aero_pix21]) ||
                lasrc_qa_is_water (ipflag[aero_pix22]))
                ipflag[curr_pix] |= (1 << IPFLAG_WATER);
        }  /* end for samp */
    }  /* end for line */

    /* Clean up the ipflag in the center of the NxN windows, for the fill
       pixels. If an NxN window is a mixture of fill and non-fill, the center
       of the window can be flagged as fill and some other QA based on the
       other pixels in that window. At the end, we want fill to be fill. */
    npixels = nlines * nsamps;
    for (curr_pix = 0; curr_pix < npixels; curr_pix++)
    {
        if (level1_qa_is_fill (qaband[curr_pix]))
            ipflag[curr_pix] = (1 << IPFLAG_FILL);
    }
}


/******************************************************************************
MODULE:  fill_with_local_average_landsat

PURPOSE:  Replaces the invalid Landsat aerosols with a local average of clear
land or water pixel values in a WxW window, for both the aerosol and eps values.
This function can process the image in a forward direction, from upper left
to lower right, or a reverse direction, from lower right to upper left.
This allows filled pixels to propagate throughout the entire image when
use_filled is true;

RETURN VALUE: None

NOTE: "fill" has two meanings here.
   The first meaning is for pixels, usually around the edge of the image,
   that do not have any valid image data; they are marked as fill in the
   Level 1 QA band. No aerosol retrieval is done for them.
   The second meaning is for pixels that have valid image data, but the aerosol
   retrieval failed. For these pixels, the aerosol value will be set ("filled")
   based on an average of nearby valid aerosol pixels.  This is the purpose
   of this function.
******************************************************************************/
static void fill_with_local_average_landsat
(
    Fill_direction_t direction, /* I: Direction to traverse the image */
    int required_clear,  /* I: Number of required clear pixels to use for the
                            local average */
    bool use_filled,     /* I: flag to use or not use previously-filled
                            pixels (as indicated by smflag) in the local
                            average */
    uint8 *ipflag,       /* I: QA flag to assist with aerosol interpolation,
                            nlines x nsamps.  */
    bool *smflag,        /* I/O: flag to indicate if the pixel was an invalid
                            aerosol and filled (true) or not filled (false),
                            nlines x nsamps */
    float *taero,        /* I/O: aerosol values for each pixel, nlines x nsamps
                            It is expected that the aerosol values are computed
                            for the center of the aerosol windows */
    float *teps,         /* I/O: angstrom coeff for each pixel, nlines x nsamps
                            It is expected that the eps values are computed
                            for the center of the aerosol windows */
    int aero_window,     /* I: size of the aerosol window (NxN) */
    int half_aero_window,/* I: size of half the aerosol window */
    int nlines,          /* I: number of lines in taero band */
    int nsamps,          /* I: number of samps in taero band */
    int *nbpixnf,        /* O: Number of pixels left (not filled) */
    int *nbpixtot        /* O: Number of total non-fill pixels */
)
{
    int start_line, start_samp; /* starting line and sample, depending on
                                   direction the image is traversed */
    int step;             /* step to move thru the image forward or backward */
    int line, samp;       /* looping variable for lines and samples for
                             the center pixel in the NxN window */
    int iline, isamp;     /* looping variable for lines and samples for
                             the WxW window around the current center pixel,
                             in increments of NxN */
    int curr_pix;         /* current pixel in 1D array of center nlines
                             and nsamps */
    int ilpix;            /* start of current line in 1D array of expanded
                             WxW window */
    int ipix;             /* current pixel in 1D array of expanded WxW window */
    int nbclrpix;         /* number of clear aerosol pixels in this window */
    int window_offset = LHALF_FIX_AERO_WINDOW - half_aero_window;
                          /* offset from the invalid pixel to the
                             first (UL) NxN aero cell in the WxW fix window */
    double sum_aero;      /* sum of the taero values in the window */
    double sum_eps;       /* sum of the teps values in the window */

    /* Set the line, sample of where we begin processing, and the step
       direction to move through the image */
    if (direction == FORWARD)
    {
        start_line = half_aero_window;
        start_samp = half_aero_window;
        step = aero_window;
    }
    else
    {
        /* When processing in reverse, find the lower right line & sample that
           are center pixels. This is where we will start processing. Set to
           the midpoint of the last NxN window whose center is in the image. */
        start_line = ((int)round((double)nlines / aero_window) - 1) *
            aero_window + half_aero_window;
        start_samp = ((int)round((double)nsamps / aero_window) - 1) *
            aero_window + half_aero_window;
        step = -aero_window;
    }

    /* Loop through the NxN center window values looking for invalid aerosol
       retrievals */
    *nbpixnf = 0;
    *nbpixtot = 0;
    for (line = start_line; line > 0 && line < nlines; line += step)
    {
        curr_pix = line * nsamps + start_samp;
        for (samp = start_samp; samp > 0 && samp < nsamps;
             samp += step, curr_pix += step)
        {
            /* Skip fill pixels */
            if (lasrc_qa_is_fill (ipflag[curr_pix]))
                continue;

            /* Increment the count of non-fill pixels */
            (*nbpixtot)++;

            /* Only process invalid, not-already-filled aerosols */
            if (lasrc_qa_is_valid_aerosol_retrieval (ipflag[curr_pix]) ||
                smflag[curr_pix])
                continue;

            /* Loop through the WxW pixels around the current center pixel,
               in increments of NxN, to look for valid aerosol retrievals
               amongst the surrounding center pixels */
            sum_aero = 0.0;
            sum_eps = 0.0;
            nbclrpix = 0;
            for (iline = line - window_offset;
                 iline <= line + window_offset;
                 iline += aero_window)
            {
                /* Make sure the window line is valid */
                if (iline < 0 || iline >= nlines)
                    continue;

                /* Determine the location of this line in our 2D array */
                ilpix = iline * nsamps;

                /* Loop through the samples of the window */
                for (isamp = samp - window_offset;
                     isamp <= samp + window_offset;
                     isamp += aero_window)
                {
                    /* Make sure the window samp is valid */
                    if (isamp < 0 || isamp >= nsamps)
                        continue;

                    /* Determine the location of this window pixel in our
                       2D array */
                    ipix = ilpix + isamp;

                    /* Check if this window pixel had valid aerosol retrieval
                       or has already been filled */
                    if (lasrc_qa_is_valid_aerosol_retrieval (ipflag[ipix]) ||
                        (use_filled && smflag[ipix]))
                    {
                        /* Add this aerosol and eps to the average */
                        nbclrpix++;
                        sum_aero += taero[ipix];
                        sum_eps += teps[ipix];
                    }
                }  /* for isamp */
            }  /* for iline */

            /* If there are enough clear pixels for computing this average
               then use the average. Otherwise increment the counter for
               number of pixels not filled within the scene. Leave the
               pixel flagged as failed aerosol thus the averaged values
               won't be used to compute other averages in this loop. */
            if (nbclrpix >= required_clear)
            {
                taero[curr_pix] = sum_aero / nbclrpix;
                teps[curr_pix] = sum_eps / nbclrpix;
                smflag[curr_pix] = true;
            }
            else
            {
                (*nbpixnf)++;
                taero[curr_pix] = DEFAULT_AERO;
                teps[curr_pix] = DEFAULT_EPS;
            }
        }  /* for samp */
    }  /* for line */
}


/******************************************************************************
MODULE:  fix_invalid_aerosols_landsat

PURPOSE:  Fixes the invalid Landsat aerosols using multiple passes through the
image, each one replacing invalid aerosols with some type of local average.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error fixing the invalid aerosols
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
int fix_invalid_aerosols_landsat
(
    uint8 *ipflag,     /* I: QA flag to assist with aerosol interpolation,
                             nlines x nsamps.  It is expected that the ipflag
                             values are computed for the center of the aerosol
                             windows. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                             It is expected that the aerosol values are computed
                             for the center of the aerosol windows */
    float *teps,       /* I/O: angstrom coeff for each pixel, nlines x nsamps
                             It is expected that the eps values are computed
                             for the center of the aerosol windows */
    int aero_window,   /* I: size of the aerosol window (NxN) */
    int half_aero_window, /* I: size of half the aerosol window */
    int nlines,        /* I: number of lines in taero band */
    int nsamps         /* I: number of samps in taero band */
)
{
    char errmsg[STR_SIZE];                        /* error message */
    char FUNC_NAME[] = "fix_invalid_aerosols_landsat"; /* function name */
    bool *smflag = NULL;  /* flag to indicate if the pixel was an invalid
                             aerosol and filled (true) or not filled (false),
                             nlines x nsamps */
    int nbpixnf;          /* number of pixels not filled in this scene */
    int nbpixtot;         /* number of non-fill pixels in this scene */

    /* Allocate memory for filled aerosols flag and initialize to false */
    smflag = calloc (nlines * nsamps, sizeof (bool));
    if (smflag == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the fill flag for "
            "invalid aerosols");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* First pass, require at least MIN_CLEAR_PIX valid surrounding values */
    fill_with_local_average_landsat (FORWARD, LMIN_CLEAR_PIX, false, ipflag,
        smflag, taero, teps, aero_window, half_aero_window, nlines, nsamps,
        &nbpixnf, &nbpixtot);
    printf ("First pass: %d pixels were not filled out of a total %d pixels\n",
        nbpixnf, nbpixtot);

    /* If the entire scene is invalid aerosols that aren't able to be filled,
       then all aerosol and eps have been set to default values and the
       ipflag is invalid.  Skip the remaining passes. */
    if (nbpixnf == nbpixtot)
    {
        printf ("Entire scene has invalid aerosols, using default taero and "
                "teps values\n");
        free (smflag);
        return (SUCCESS);
    }

    /** Second pass, require at least 1 valid surrounding value **/
    if (nbpixnf > 0)
    {
        fill_with_local_average_landsat (FORWARD, 1, true, ipflag, smflag,
            taero, teps, aero_window, half_aero_window, nlines, nsamps,
            &nbpixnf, &nbpixtot);
        printf ("Second pass: %d pixels were not filled out of a total %d "
            "pixels\n", nbpixnf, nbpixtot);
    }

    /** Final reverse pass for any remaining invalid retrievals (mostly the
        UL part of the image) **/
    if (nbpixnf > 0)
    {
        fill_with_local_average_landsat(REVERSE, 1, true, ipflag, smflag,
            taero, teps, aero_window, half_aero_window, nlines, nsamps,
            &nbpixnf, &nbpixtot);
        printf ("Final pass: %d pixels were not filled out of a total %d "
            "pixels\n", nbpixnf, nbpixtot);
    }

    /* Free the allocated memory */
    free (smflag);

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  aerosol_interp_sentinel

PURPOSE:  Interpolates the Sentinel aerosol values throughout the image using
the aerosols that were calculated for the UL pixel of each NxN window. Also
cleans up the fill pixels in the ipflag.

RETURN VALUE:
Type = N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
void aerosol_interp_sentinel
(
    int aero_window,   /* I: size of the aerosol window */
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps.  It is expected that the ipflag
                               values are computed for the UL of the aerosol
                               windows to match taero. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                          It is expected that the aerosol values are computed
                          for the UL of the aerosol windows.  This routine
                          will fill in the pixels for the remaining pixels of
                          the window. */
    int nlines,        /* I: number of lines in ipflag & taero bands */
    int nsamps         /* I: number of samps in ipflag & taero bands */
)
{
    int line, samp;    /* looping variable for lines and samples */
    int iline, isamp;  /* looping variable for lines and samples in the
                          aerosol window */
    int curr_pix;      /* current pixel in 1D arrays of nlines x nsamps */
    int curr_win_pix;  /* current pixel in the nxn window for atm corr */
    int next_samp_pix; /* pixel location of the next sample */
    int next_line_pix; /* pixel location of the next line */
    int next_line_samp_pix; /* pixel location of next line and next sample */
    int awline;        /* line for the next aerosol window */
    int awsamp;        /* sample for the next aerosol window */
    int sq_aero_win;   /* square of the aerosol window */
    int awline_iline;  /* awline - iline */
    int iline_line;    /* iline - line */

    /* Use the UL corner of the aerosol windows to interpolate the remaining
       pixels in the window */
    sq_aero_win = aero_window * aero_window;
    for (line = 0; line < nlines; line++)
    {
        /* Determine the current pixel */
        curr_pix = line * nsamps;

        /* Determine the line for the next aerosol window */
        awline = line + aero_window;

        for (samp = 0; samp < nsamps; samp++, curr_pix++)
        {
            /* Skip fill values */
            if (level1_qa_is_fill (qaband[curr_pix]))
                continue;

            /* Determine the sample for the next aerosol window */
            awsamp = samp + aero_window;

            /* Determine the next line and next sample to be used for
               interpolating */
            next_samp_pix = line * nsamps + awsamp;
            next_line_pix = awline * nsamps + samp;
            next_line_samp_pix = awline * nsamps + awsamp;

            /* Loop through an NxN window with the current pixel being the UL
               corner of the window */
            for (iline = line; iline < awline; iline++)
            {
                /* Skip if this isn't a valid line */
                if (iline >= nlines) continue;

                awline_iline = awline - iline;
                iline_line = iline - line;
                for (isamp = samp; isamp < awsamp; isamp++)
                {
                    /* Skip if this isn't a valid sample */
                    if (isamp >= nsamps) continue;

                    /* Skip fill values */
                    curr_win_pix = iline * nsamps + isamp;
                    if (level1_qa_is_fill (qaband[curr_win_pix])) continue;

                    /* Interpolate using surrounding pixels */
                    taero[curr_win_pix] = taero[curr_pix] *
                        (awline_iline) * (awsamp-isamp);

                    if ((awline < nlines) && (awsamp < nsamps))
                        taero[curr_win_pix] +=
                            (isamp-samp)*(awline_iline) * taero[next_samp_pix] +
                            (awsamp-isamp)*(iline_line) * taero[next_line_pix] +
                            (isamp-samp)*(iline_line) *
                             taero[next_line_samp_pix];

                    else if ((awline >= nlines) && (awsamp < nsamps))
                        taero[curr_win_pix] +=
                            (isamp-samp)*(awline_iline) * taero[next_samp_pix] +
                            (awsamp-isamp)*(iline_line) * taero[curr_pix] +
                            (isamp-samp)*(iline_line) * taero[next_samp_pix];

                    else if ((awline < nlines) && (awsamp >= nsamps))
                        taero[curr_win_pix] +=
                            (isamp-samp)*(awline_iline) * taero[curr_pix] +
                            (awsamp-isamp)*(iline_line) * taero[next_line_pix] +
                            (isamp-samp)*(iline_line) * taero[next_line_pix];

                    else if ((awline >= nlines) && (awsamp >= nsamps))
                        taero[curr_win_pix] +=
                            (isamp-samp)*(awline_iline) * taero[curr_pix] +
                            (awsamp-isamp)*(iline_line) * taero[curr_pix] +
                            (isamp-samp)*(iline_line) * taero[curr_pix];

                    /* Compute the average */
                    taero[curr_win_pix] /= sq_aero_win;
                }
            }
        }  /* for samp */
    }  /* for line */
}


/******************************************************************************
MODULE:  ipflag_expand_failed_sentinel

PURPOSE:  For every failed Sentinel pixel (possible urban, barren,
turbid water, etc.), set the surrounding 25x25 clear or interpolated window
pixels to possible failed as well.

RETURN VALUE:
Type = N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
/* Basically look back 2 window pixels and forward 2 window pixels, since each
   window is 6 pixels in size. */
#define HALF_EXPAND_WIN 12
void ipflag_expand_failed_sentinel
(
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps. It is expected that the ipflag
                               values are computed for all pixels. */
    int nlines,        /* I: number of lines in ipflag band */
    int nsamps         /* I: number of samps in ipflag band */
)
{
    int line, samp;         /* looping variable for lines and samples */
    int iline, isamp;       /* looping variable for window lines and samples */
    int win_line, win_samp; /* current line, sample in the window */
    int curr_pix;           /* current pixel in 1D arrays of nlines * nsamps */
    int curr_win_pix;       /* current window pixel in 1D arrays */
    long npixels;           /* number of pixels to process */

    /* Loop through the lines and samples to expand the window around the
       possible failed pixels with a temp value */
    for (line = 0; line < nlines; line++)
    {
        curr_pix = line * nsamps;
        for (samp = 0; samp < nsamps; samp++, curr_pix++)
        {
            /* If not failed then move on */
            if (!btest (ipflag[curr_pix], IPFLAG_FAILED))
                continue;

            /* If this is a failed pixel then expand the failed pixels around
               this pixel */
            for (iline = -HALF_EXPAND_WIN; iline <= HALF_EXPAND_WIN; iline++)
            {
                win_line = line + iline;
                if (win_line < 0 || win_line >= nlines)
                    continue;
                for (isamp = -HALF_EXPAND_WIN; isamp <= HALF_EXPAND_WIN;
                     isamp++)
                {
                    win_samp = samp + isamp;
                    if (win_samp < 0 || win_samp >= nsamps)
                        continue;

                    /* Set non-fill pixels to failed. Note: this will set
                       non-representative pixels to failed. */
                    /* NOTE: FORTRAN code does not fail water pixels */
/** C3 TODO Remove the water check **/
                    curr_win_pix = win_line * nsamps + win_samp;
                    if (!btest (ipflag[curr_win_pix], IPFLAG_FILL) &&
                        !btest (ipflag[curr_win_pix], IPFLAG_WATER) &&
                        !btest (ipflag[curr_win_pix], IPFLAG_FAILED_TMP))
                        ipflag[curr_win_pix] |= (1 << IPFLAG_FAILED_TMP);
                }
            }
        }
    }

    /* Loop through the pixels and reset any temporary possible failed pixel
       to a possible failed pixel */
    npixels = nlines * nsamps;
    for (curr_pix = 0; curr_pix < npixels; curr_pix++)
        if (btest (ipflag[curr_pix], IPFLAG_FAILED_TMP))
        {
            ipflag[curr_pix] |= (1 << IPFLAG_FAILED);
            ipflag[curr_pix] &= ~(1 << IPFLAG_FAILED_TMP);
        }
}


/******************************************************************************
MODULE:  aero_avg_failed_sentinel

PURPOSE:  For every failed Sentinel pixel, determine the 61x61 average of clear
pixels around that pixel. That average value will be used to replace the failed pixels.

RETURN VALUE:
Type = N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
#define HALF_FAILED_WIN 30
#define MIN_VALID_WINDOW_PIX 20
int aero_avg_failed_sentinel
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps. It is expected that the ipflag
                               values are computed for all pixels. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                               updated with average values if needed */
    float *teps,       /* I/O: angstrom coeff for each pixel, nlines x nsamps
                               updated with average values if needed */
    int nlines,        /* I: number of lines in ipflag band */
    int nsamps         /* I: number of samps in ipflag band */
)
{
    char errmsg[STR_SIZE];                         /* error message */
    char FUNC_NAME[] = "aero_avg_failed_sentinel"; /* function name */
    int line, samp;        /* looping variable for lines and samples */
    int iline, isamp;      /* looping variable for window lines and samples */
    int curr_pix;          /* current pixel in 1D arrays of nlines * nsamps */
    int curr_win_pix;      /* current window pixel in 1D arrays */
    int nbpixnf = 0;       /* number of pixels not filled */
    int nbaeroavg = 0;     /* number of pixels counted in the average */
    int ipass;             /* number of passes required to fill failed pixels */
    long npixels;          /* number of pixels to process */
    float taerosum;        /* sum of aerosols in the NxN window */
    float tepssum;         /* sum of angstrom coeff in the NxN window */
    bool one_filled;       /* flag to indicate at least one pixel was filled
                              in the first/stringent fill loop */
    float *taeros=NULL;    /* average aerosol values for each pixel,
                              nlines x nsamps */
    float *tepss=NULL;     /* average angstrom coeff for each pixel,
                              nlines x nsamps */
    bool *smflag=NULL;     /* flag to indicate the pixel was filled (smoothed),
                              nlines x nsamps */

    /* Allocate memory for the intermediate arrays */
    taeros = calloc (nlines * nsamps, sizeof (float));
    if (taeros == NULL)
    {
        sprintf (errmsg, "Error allocating memory for taeros band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    tepss = calloc (nlines * nsamps, sizeof (float));
    if (tepss == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tepss band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    smflag = calloc (nlines * nsamps, sizeof (bool));
    if (smflag == NULL)
    {
        sprintf (errmsg, "Error allocating memory for smflag band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through the lines and samples to determine an average of the
       surrounding window (with a required minimum count of valid pixels) for
       both teps and taero */
    one_filled = false;
#ifdef _OPENMP
    #pragma omp parallel for private (line, samp, curr_pix, taerosum, tepssum, nbaeroavg, iline, isamp, curr_win_pix) reduction(+:nbpixnf)
#endif
    for (line = 0; line < nlines; line++)
    {
        curr_pix = line * nsamps;
        for (samp = 0; samp < nsamps; samp++, curr_pix++)
        {
            smflag[curr_pix] = false;

            /* If this pixel is fill, then don't process */
            if (level1_qa_is_fill (qaband[curr_pix]))
                continue;

            /* Initialize variables */
            taerosum = 0.0;
            tepssum = 0.0;
            nbaeroavg = 0;

            /* Look at the surrounding window and sum up the values for
               non-fill and non-failed pixels (sums water and land) */
            for (iline = -HALF_FAILED_WIN; iline <= HALF_FAILED_WIN; iline++)
            {
                if ((line + iline < 0) || (line + iline >= nlines))
                    continue;
                for (isamp = -HALF_FAILED_WIN; isamp <= HALF_FAILED_WIN;
                     isamp++)
                {
                    if ((samp + isamp < 0) || (samp + isamp >= nsamps))
                        continue;

                    curr_win_pix = (line + iline) * nsamps + (samp + isamp);
                    if (!btest (ipflag[curr_win_pix], IPFLAG_FILL) &&
                        !btest (ipflag[curr_win_pix], IPFLAG_FAILED))
                    {
                        nbaeroavg++;
                        taerosum += taero[curr_win_pix];
                        tepssum += teps[curr_win_pix];
                    }
                }
            }

            /* If more than the required minimum of the surrounding window
               pixels were used in the average, then we can use the average
               for this pixel */
            if (nbaeroavg > MIN_VALID_WINDOW_PIX)
            {
                /* Use the window average. Set the one_filled flag to indicate
                   at least one pixel was filled. */
                taeros[curr_pix] = taerosum / nbaeroavg;
                tepss[curr_pix] = tepssum / nbaeroavg;
                smflag[curr_pix] = true;
#ifdef _OPENMP
                #pragma omp atomic write
#endif
                one_filled = true;
            }
            else
            {
                nbpixnf++;
            }
        }  /* end for samp */
    }  /* end for line */

    /* If none of the pixels were able to be filled, then use default values */
    if (!one_filled)
    {
        npixels = nlines * nsamps;
        for (curr_pix = 0; curr_pix < npixels; curr_pix++)
        {
            if (!btest (ipflag[curr_pix], IPFLAG_FILL))
            {
                taero[curr_pix] = DEFAULT_AERO;
                teps[curr_pix] = DEFAULT_EPS;
                smflag[curr_pix] = true;
            }
        }
        return (SUCCESS);
    }

    /* Second pass through to handle the non-filled pixels */
    printf ("Second Pass");
    ipass = 2;
    while (nbpixnf != 0)
    {
        printf ("  Pass number: %d, number of non-filled pixels: %d\n", ipass,
            nbpixnf);
        nbpixnf = 0;

        /* Loop through the lines and samples to expand the window around the
           failed pixels with a temp value. smflag is updated in the for loops
           as well as used. Therefore we cannot use multi-threading because
           that would randomly update smflag vs. updating it top to bottom,
           left to right in a conrolled manner. The update and use of smflag
           in this loop matches the FORTRAN algorithm. */
        for (line = 0; line < nlines; line++)
        {
            curr_pix = line * nsamps;
            for (samp = 0; samp < nsamps; samp++, curr_pix++)
            {
                /* If this pixel is fill or the pixel has already been filled,
                   then don't process */
                if (level1_qa_is_fill (qaband[curr_pix]) || smflag[curr_pix])
                    continue;

                /* Initialize variables */
                taerosum = 0.0;
                tepssum = 0.0;
                nbaeroavg = 0;

                /* Look at the surrounding window and sum up the averaged
                   values for already filled pixels. Note fill pixels should
                   never be filled in the above loops, so we don't have to
                   specifically check for fill pixels. */
                for (iline = -HALF_FAILED_WIN; iline <= HALF_FAILED_WIN;
                     iline++)
                {
                    if ((line + iline < 0) || (line + iline >= nlines))
                        continue;
                    for (isamp = -HALF_FAILED_WIN; isamp <= HALF_FAILED_WIN;
                         isamp++)
                    {
                        if ((samp + isamp < 0) || (samp + isamp >= nsamps))
                            continue;

                        curr_win_pix = (line + iline) * nsamps + (samp + isamp);
                        if (smflag[curr_win_pix])
                        {
                            nbaeroavg++;
                            taerosum += taeros[curr_win_pix];
                            tepssum += tepss[curr_win_pix];
                        }
                    }
                }

                /* If any of the window pixels were used in the average, then we
                   can use the average for this pixel */
                if (nbaeroavg > 0)
                {
                    /* Use the window average */
                    taeros[curr_pix] = taerosum / nbaeroavg;
                    tepss[curr_pix] = tepssum / nbaeroavg;
                    smflag[curr_pix] = true;
                }
                else
                {
                    nbpixnf++;
                }
            }  /* end for samp */
        }  /* end for line */

        /* Update the pass number */
        ipass++;
    }  /* end do while */

    /* Loop through the lines and samples and use the average pixels to
       fill failed pixels */
    npixels = nlines * nsamps;
    for (curr_pix = 0; curr_pix < npixels; curr_pix++)
    {
        if (btest (ipflag[curr_pix], IPFLAG_FAILED) &&
           !btest (ipflag[curr_pix], IPFLAG_FILL))
        {
            taero[curr_pix] = taeros[curr_pix];
            teps[curr_pix] = tepss[curr_pix];
            ipflag[curr_pix] |= (1 << IPFLAG_FIXED);
        }
    }

    /* Free the allocated memory */
    free (taeros);
    free (tepss);
    free (smflag);

    /* Successful completion */
    return (SUCCESS);
}
