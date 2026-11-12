/******************************************************************************
FILE: gapfill_viirs_aux
  
PURPOSE: Contains functions for gapfilling the VIIRS auxiliary data.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
******************************************************************************/
#include "gapfill_viirs_aux.h"

/* Program will look for these datasets in the VIIRS inputs */
/* Ozone - uint8
   Water vapor - uint16 */
#define N_DATASETS 2
char dataset_path[50] = "/HDFEOS/GRIDS/VIIRS_CMG/Data Fields/";
char list_of_datasets[N_DATASETS][50] = {
    "Coarse Resolution Ozone",
    "Coarse Resolution Water Vapor"};
#define OZONE 0
#define WV 1

/* Expected size of the CMG dataset is 3600 lines x 7200 samples */
#define CMG_NDIMS 2
#define CMG_NLINES 3600
#define CMG_NSAMPS 7200

/******************************************************************************
MODULE:  cmp_vals

PURPOSE:  Compares the absolute difference of the two floating point values to
determine if they are essentially equivalent (close enough due to rounding
issues).

RETURN VALUE:
Type = bool
Value          Description
-----          -----------
false          Values are not the same
true           Values are the same

NOTES:
******************************************************************************/
bool cmp_vals
(
    float val1,   /* I: float value 1 */
    float val2    /* I: float value 2 */
)
{
    if (fabs (val1 - val2) <= 0.0001)
        return true;
    else
        return false;
}


/******************************************************************************
MODULE:  file_exists

PURPOSE:  Determines if the specified file exists.

RETURN VALUE:
Type = bool
Value          Description
-----          -----------
false          File does not exist
true           File exists

NOTES:
******************************************************************************/
bool file_exists
(
    char *filename   /* I: file to check for existence */
)
{
    if (access (filename, F_OK) == 0)
        return true;
    else
        return false;
}


/* DAYSTEP is 50/15 . Weight the monthly data by 50% at day 1, another 50%
 * in 15 days (from day 1 to day 15) and then down to just 50% at day 30
 */
#define DAYSTEP 3.3333333

/******************************************************************************
MODULE:  determine_weights

PURPOSE:  Based on the day of month, the weighting for filling the gaps using
the monthly averages is set up to use the target/current month, the previous
month, and the next month.

RETURN VALUE:
Type = None

NOTES:
 1. Weight the monthly data by 50% at day 1, another 50% in 15 days (from
    day 1 to day 15) and then down to just 50% at day 30.
******************************************************************************/
void determine_weights
(
    int aux_day,             /* I: day of the auxiliary file (1-31) */
    float *prev_weight,      /* O: weight for the previous month */
    float *target_weight,    /* O: weight for the current/target month */
    float *next_weight       /* O: weight for the next month */
)
{
    int i;    /* looping variable for the day within the month */

    /* Determine the weighting for this day (a very crude sawtooth function).
       Start with a 50/50 weight between the current month and the previous
       month. */
    *target_weight = 50.0;
    *prev_weight = 50.0;
    *next_weight = 0.0;

    /* Weighting for the first half of the month (day = 1-15) */
    if (aux_day <= 15)
    {
        for (i = 1; i <= 15; i++)
        {
            if (i >= aux_day)
                break;

            *target_weight += DAYSTEP;
            *prev_weight -= DAYSTEP;
        }
    }
    else
    {
        /* These are the values that come out of the above for loop for days
           that are above 15 */
        *target_weight = 100.0;
        *prev_weight = 0.0;
    }

    /* Weighting for the second half of the month (day = 16-31) */
    if (aux_day >= 15)
    {
        for (i = 16; i <= 31; i++)
        {
            if (i == aux_day)
                break;

            *target_weight -= DAYSTEP;
            *next_weight += DAYSTEP;
        }
    }

    /* Adjust the weights to zero if they are below 3.0 */
    if (*prev_weight < 3.0)
        *prev_weight = 0.0;
    if (*next_weight < 3.0)
        *next_weight = 0.0;
}


/******************************************************************************
MODULE:  read_monthly_avgs

PURPOSE:  Determines the previous, target, and next monthly averages, and
based on the weights, reads those monthly averages to be used for gapfilling.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading the monthly averages or the monthly
               averages are not available
SUCCESS        Successful completion

NOTES:
******************************************************************************/
int read_monthly_avgs
(
    int aux_month,        /* I: month of the auxiliary file (1-12) */
    int aux_year,         /* I: year of the auxiliary file */
    int n_pixels,         /* I: number of pixels in the monthly avg 1D array */
    float prev_weight,    /* I: weight for the previous month */
    float target_weight,  /* I: weight for the current/target month */
    float next_weight,    /* I: weight for the next month */
    uint8_t* monthly_avg_oz_data[3],  /* O: array of the ozone data for
                                previous month, target month, next month. If
                                the weight is 0, then the data will be NULL. */
    uint16_t* monthly_avg_wv_data[3]  /* O: array of the water vapor data
                                for previous month, target month, next month. If
                                the weight is 0, then the data will be NULL. */
)
{
    char FUNC_NAME[] = "read_monthly_avgs"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char envvar[STR_SIZE];   /* LASRC_AUX_DIR environment variable contents */
    char aux_dir[STR_SIZE];  /* location of monthly averages aux files */
    char oz_monthly_file[3][STR_SIZE]; /* input ozone monthly files (previous
                                month, current month, next month) */
    char wv_monthly_file[3][STR_SIZE]; /* input water vapor monthly files
                                (previous month, current month, next month)*/
    int prev_month;          /* previous month for filling gaps */
    int next_month;          /* next month for filling gaps */
    FILE *mavg_fp = NULL;    /* file pointer for monthly averages */

    /* Determine which months will be used for gapfilling */
    prev_month = aux_month - 1;
    if (prev_month <= 0)
        prev_month = 12;

    next_month = aux_month + 1;
    if (next_month > 12)
        next_month = 1;

    /* Get the LASRC auxiliary environment variable */
    if (!getenv ("LASRC_AUX_DIR"))
    {
        sprintf (errmsg, "LASRC_AUX_DIR environment variable is not defined.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    snprintf (envvar, STR_SIZE, "%s", getenv("LASRC_AUX_DIR"));
    sprintf (aux_dir, "%s/monthly_avgs", envvar);

    /* Determine the name of the monthly average files to be opened - previous
       month, current/target month, next month. For the previous month, start
       with the current year.  If that isn't available, then we will have to
       pull the previous year.  For the target and next month, we will have to
       pull the previous year because the monthly averages for those months
       aren't available in the current year. */
    sprintf (oz_monthly_file[0], "%s/%d/monthly_avg_oz_%4d_%02d.img", aux_dir,
        aux_year, aux_year, prev_month);
    sprintf (wv_monthly_file[0], "%s/%d/monthly_avg_wv_%4d_%02d.img", aux_dir,
        aux_year, aux_year, prev_month);
    sprintf (oz_monthly_file[1], "%s/%d/monthly_avg_oz_%4d_%02d.img", aux_dir,
        aux_year-1, aux_year-1, aux_month);
    sprintf (wv_monthly_file[1], "%s/%d/monthly_avg_wv_%4d_%02d.img", aux_dir,
        aux_year-1, aux_year-1, aux_month);
    sprintf (oz_monthly_file[2], "%s/%d/monthly_avg_oz_%4d_%02d.img", aux_dir,
        aux_year-1, aux_year-1, next_month);
    sprintf (wv_monthly_file[2], "%s/%d/monthly_avg_wv_%4d_%02d.img", aux_dir,
        aux_year-1, aux_year-1, next_month);

    /* Initialize the monthly average ozone and water vapor pointers to NULL
       for the previous, target, and next month datasets */
    monthly_avg_oz_data[0] = NULL;
    monthly_avg_oz_data[1] = NULL;
    monthly_avg_oz_data[2] = NULL;
    monthly_avg_wv_data[0] = NULL;
    monthly_avg_wv_data[1] = NULL;
    monthly_avg_wv_data[2] = NULL;

    /* Check for the existence and open the previous month's averages, but
       only if they will be used */
    if (prev_weight > 0.0)
    {
        /* If the file doesn't exist for the current year, then check last
           year */
        /* ozone */
        if (!file_exists (oz_monthly_file[0]))
        {
            sprintf (oz_monthly_file[0], "%s/%d/monthly_avg_oz_%4d_%02d.img",
                aux_dir, aux_year-1, aux_year-1, prev_month);
            if (!file_exists (oz_monthly_file[0]))
            {
                sprintf (errmsg, "Monthly ozone averages for the previous "
                    "month (%d) do not exist for the current year (%d) or the "
                    "previous year (%d). %s", prev_month, aux_year, aux_year-1,
                    oz_monthly_file[0]);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }

        /* Allocate data for the previous monthly avg array */
        monthly_avg_oz_data[0] = calloc (n_pixels, sizeof(uint8_t));
        if (monthly_avg_oz_data[0] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the previous monthly "
                "ozone average");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Read the data for return */
        printf ("Previous monthly averages OZ file: %s\n", oz_monthly_file[0]); 
        mavg_fp = fopen (oz_monthly_file[0], "r");
        if (!mavg_fp)
        {
            sprintf (errmsg, "Not able to open the monthly ozone average: %s",
                oz_monthly_file[0]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (fread (&monthly_avg_oz_data[0][0], sizeof(uint8_t), n_pixels,
            mavg_fp) != n_pixels)
        {
            sprintf (errmsg, "Error reading the monthly ozone average: %s",
                oz_monthly_file[0]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        fclose (mavg_fp);

        /* water vapor */
        if (!file_exists (wv_monthly_file[0]))
        {
            sprintf (wv_monthly_file[0], "%s/%d/monthly_avg_wv_%4d_%02d.img",
                aux_dir, aux_year-1, aux_year-1, prev_month);
            if (!file_exists (wv_monthly_file[0]))
            {
                sprintf (errmsg, "Monthly water vapor averages for the "
                    "previous month (%d) do not exist for the current year "
                    "(%d) or the previous year (%d). %s", prev_month, aux_year,
                    aux_year-1, wv_monthly_file[0]);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }

        /* Allocate data for the previous monthly avg array */
        monthly_avg_wv_data[0] = calloc (n_pixels, sizeof(uint16_t));
        if (monthly_avg_wv_data[0] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the previous monthly "
                "water vapor average");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Read the data for return */
        printf ("Previous monthly averages WV file: %s\n", wv_monthly_file[0]); 
        mavg_fp = fopen (wv_monthly_file[0], "r");
        if (!mavg_fp)
        {
            sprintf (errmsg, "Not able to open the monthly WV average: %s",
                wv_monthly_file[0]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (fread (&monthly_avg_wv_data[0][0], sizeof(uint16_t), n_pixels,
            mavg_fp) != n_pixels)
        {
            sprintf (errmsg, "Error reading the monthly WV average: %s",
                wv_monthly_file[0]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        fclose (mavg_fp);
    }

    /* Check for the existence and open the target month's averages */
    /* ozone */
    if (!file_exists (oz_monthly_file[1]))
    {
        sprintf (errmsg, "Monthly ozone averages for the target month (%d) "
            "do not exist for the previous year (%d). %s", aux_month,
            aux_year-1, oz_monthly_file[1]);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Allocate data for the target monthly avg array */
    monthly_avg_oz_data[1] = calloc (n_pixels, sizeof(uint8_t));
    if (monthly_avg_oz_data[1] == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the target monthly "
            "ozone average");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data for return */
    printf ("Target monthly averages OZ file: %s\n", oz_monthly_file[1]); 
    mavg_fp = fopen (oz_monthly_file[1], "r");
    if (!mavg_fp)
    {
        sprintf (errmsg, "Not able to open the monthly ozone average: %s",
            oz_monthly_file[1]);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
 
    if (fread (&monthly_avg_oz_data[1][0], sizeof(uint8_t), n_pixels,
        mavg_fp) != n_pixels)
    {
        sprintf (errmsg, "Error reading the monthly ozone average: %s",
            oz_monthly_file[1]);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
 
    fclose (mavg_fp);

    /* water vapor */
    if (!file_exists (wv_monthly_file[1]))
    {
        sprintf (errmsg, "Monthly water vapor averages for the target month "
            "(%d) do not exist for the previous year (%d). %s", aux_month,
            aux_year-1, wv_monthly_file[1]);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Allocate data for the target monthly avg array */
    monthly_avg_wv_data[1] = calloc (n_pixels, sizeof(uint16_t));
    if (monthly_avg_wv_data[1] == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the target monthly "
            "water vapor average");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data for return */
    printf ("Target monthly averages WV file: %s\n", wv_monthly_file[1]); 
    mavg_fp = fopen (wv_monthly_file[1], "r");
    if (!mavg_fp)
    {
        sprintf (errmsg, "Not able to open the monthly WV average: %s",
            wv_monthly_file[1]);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (fread (&monthly_avg_wv_data[1][0], sizeof(uint16_t), n_pixels,
        mavg_fp) != n_pixels)
    {
        sprintf (errmsg, "Error reading the monthly WV average: %s",
            wv_monthly_file[1]);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    fclose (mavg_fp);

    /* Check for the existence and open the next month's averages, but only
       if they will be used */
    if (next_weight > 0.0)
    {
        /* ozone */
        if (!file_exists (oz_monthly_file[2]))
        {
            sprintf (errmsg, "Monthly ozone averages for the next month (%d) "
                "do not exist for the previous year (%d). %s", next_month,
                aux_year-1, oz_monthly_file[2]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Allocate data for the next monthly avg array */
        monthly_avg_oz_data[2] = calloc (n_pixels, sizeof(uint8_t));
        if (monthly_avg_oz_data[2] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the next monthly "
                "ozone average");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Read the data for return */
        printf ("Next monthly averages OZ file: %s\n", oz_monthly_file[2]); 
        mavg_fp = fopen (oz_monthly_file[2], "r");
        if (!mavg_fp)
        {
            sprintf (errmsg, "Not able to open the monthly ozone average: %s",
                oz_monthly_file[2]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (fread (&monthly_avg_oz_data[2][0], sizeof(uint8_t), n_pixels,
            mavg_fp) != n_pixels)
        {
            sprintf (errmsg, "Error reading the monthly ozone average: %s",
                oz_monthly_file[2]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        fclose (mavg_fp);

        /* water vapor */
        if (!file_exists (wv_monthly_file[2]))
        {
            sprintf (errmsg, "Monthly water vapor averages for the next month "
                "(%d) do not exist for the previous year (%d). %s", next_month,
                aux_year-1, wv_monthly_file[2]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Allocate data for the next monthly avg array */
        monthly_avg_wv_data[2] = calloc (n_pixels, sizeof(uint16_t));
        if (monthly_avg_wv_data[2] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the next monthly "
                "water vapor average");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Read the data for return */
        printf ("Next monthly averages WV file: %s\n", wv_monthly_file[2]);
        mavg_fp = fopen (wv_monthly_file[2], "r");
        if (!mavg_fp)
        {
            sprintf (errmsg, "Not able to open the monthly WV average: %s",
                wv_monthly_file[2]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (fread (&monthly_avg_wv_data[2][0], sizeof(uint16_t), n_pixels,
            mavg_fp) != n_pixels)
        {
            sprintf (errmsg, "Error reading the monthly WV average: %s",
                wv_monthly_file[2]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        fclose (mavg_fp);
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_fill_value

PURPOSE:  Determines the weighted average to be used for filling the gaps. If
the monthly average is fill, then it won't be used.

RETURN VALUE:
Type = float
Value          Description
-----          -----------
float          Weighted average

NOTES:
******************************************************************************/
float get_fill_value
(
    float prev_weight,    /* I: weight for the previous month */
    float target_weight,  /* I: weight for the current/target month */
    float next_weight,    /* I: weight for the next month */
    float prev_avg,       /* I: previous month average for this pixel */
    float target_avg,     /* I: current/target month average for this pixel */
    float next_avg        /* I: next month average for this pixel */
)
{
    float wgt_avg;     /* weighted average for this pixel */

    /* If any of the monthly averages are fill, then skip their use but
       transfer (for this pixel only) their weight to the other non-fill weight
       so the weights always total 100. If all averages are fill, then the
       transfer of weight doesn't matter because the result will be zero.
       The goal is to add the fill weight to the highest weight, if there are
       more than one non-zero weights. Otherwise, the previous weight will be
       the highest priority to transfer weights to because it's within the
       current year. The next priority will be the target weight, followed by
       the next month's weight.
       Given that the averages are floating point values, we will round them
       to the nearest integer and compare to the VIIRS_FILL value which is
       zero. */
    if (cmp_vals (prev_avg, VIIRS_FILL))
    {
        /* Add to the highest weight if both are not fill */
        if (!cmp_vals (target_avg, VIIRS_FILL) &&
            !cmp_vals (next_avg, VIIRS_FILL))
        {
            if (target_weight >= next_weight)
                target_weight += prev_weight;
            else
                next_weight += prev_weight;
        }
        else if (!cmp_vals (target_avg, VIIRS_FILL))
            target_weight = 100.0;
        else if (!cmp_vals (next_avg, VIIRS_FILL))
            next_weight = 100.0;
    }

    if (cmp_vals (target_avg, VIIRS_FILL))
    {
        /* Add to the highest weight if both are not fill */
        if (!cmp_vals (prev_avg, VIIRS_FILL) &&
            !cmp_vals (next_avg, VIIRS_FILL))
        {
            if (prev_weight >= next_weight)
                prev_weight += target_weight;
            else
                next_weight += target_weight;
        }
        else if (!cmp_vals (prev_avg, VIIRS_FILL))
            prev_weight = 100.0;
        else if (!cmp_vals (next_avg, VIIRS_FILL))
            next_weight = 100.0;
    }

    if (cmp_vals (next_avg, VIIRS_FILL))
    {
        /* Add to the highest weight if both are not fill */
        if (!cmp_vals (prev_avg, VIIRS_FILL) &&
            !cmp_vals (target_avg, VIIRS_FILL))
        {
            if (prev_weight >= target_weight)
                prev_weight += next_weight;
            else
                target_weight += next_weight;
        }
        else if (!cmp_vals (prev_avg, VIIRS_FILL))
            prev_weight = 100.0;
        else if (!cmp_vals (target_avg, VIIRS_FILL))
            target_weight = 100.0;
    }

    /* Determine the sum of the weighted monthly averages */
    wgt_avg = target_avg * (target_weight * 0.01);
    if (prev_weight > 0.0)
        wgt_avg += prev_avg * (prev_weight * 0.01);
    if (next_weight > 0.0)
        wgt_avg += next_avg * (next_weight * 0.01);

    return (wgt_avg);
}


/******************************************************************************
MODULE:  gapfill_viirs_aux

PURPOSE:  Opens the user-specified VIIRS auxiliary product, reads the
identified datasets, fills the gaps in the data, and writes the dataset back
out to the same file.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading the inputs, gapfilling, or writing the
               data back out
SUCCESS        Successful completion

NOTES:
1. This routine will work on the VIIRS JPSS1 (VJ104ANC) and the NPP (VNP04ANC)
   products.
******************************************************************************/
int main (int argc, char **argv)
{    
    char FUNC_NAME[] = "main"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char *viirs_aux_file = NULL;  /* input VIIRS auxiliary file */
    long pix;                /* current pixel location in the 1D array */
    int line, samp;          /* looping variables for line and sample */
    int n_pixels;            /* number of pixels in this 1D array */
    int retval;              /* return status */
    int aux_month;           /* month of the auxiliary file (1-12) */
    int aux_day;             /* day of the auxiliary file (1-31) */
    int aux_year;            /* year of the auxiliary file */
    float target_weight;     /* weight for the current/target month */
    float prev_weight;       /* weight for the previous month */
    float next_weight;       /* weight for the next month */
    float prev_oz;           /* previous ozone value for the monthly avg */
    float target_oz;         /* target ozone value for the monthly avg */
    float next_oz;           /* next ozone value for the monthly avg */
    float prev_wv;           /* previous water vapor value for monthly avg */
    float target_wv;         /* target water vapor value for the monthly avg */
    float next_wv;           /* next water vapor value for the monthly avg */
    uint8_t **oz_2d = NULL;  /* 2D array for reading and writing from the ozone
                                dataset */
    uint16_t **wv_2d = NULL; /* 2D array for reading and writing from the water
                                vapor dataset */
    uint8_t *monthly_avg_oz_data[3];  /* array of 1D ozone data for previous
                                         month, target month, next month */
    uint16_t *monthly_avg_wv_data[3]; /* array of 1D water vapor data for prev
                                         month, target month, next month */
    hid_t file_id;           /* VIIRS file id */
    hid_t ozone_dsid;        /* ozone dataset ID */
    hid_t wv_dsid;           /* water vapor dataset ID */

    /* Read the command-line arguments */
    retval = get_args (argc, argv, &aux_month, &aux_day, &aux_year,
        &viirs_aux_file);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        exit (ERROR);
    }

    /* Open the input VIIRS file and get the dataset IDs along with the
       dimensions.  This routine already confirms the existence of the
       ozone and water vapor datasets, confirms the data types are as expected,
       and confirms the size of the dataset matches our CMG grid size. */
    retval = open_oz_wv_datasets (viirs_aux_file, &file_id, &ozone_dsid,
        &wv_dsid);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error parsing file: %s", viirs_aux_file);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Allocate memory for the 2D ozone and water vapor arrays, however the
       data needs to be contiguous as a 1D array. */
    oz_2d = (uint8_t **) calloc (CMG_NLINES, sizeof(uint8_t*));
    if (oz_2d == NULL)
    {
        sprintf (errmsg, "Error allocating memory for ozone uint8 pointers");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    wv_2d = (uint16_t **) calloc (CMG_NLINES, sizeof(uint16_t*));
    if (wv_2d == NULL)
    {
        sprintf (errmsg, "Error allocating memory for WV uint16 pointers");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    oz_2d[0] = (uint8_t*) calloc (CMG_NLINES * CMG_NSAMPS, sizeof(uint8_t));
    if (oz_2d[0] == NULL)
    {
        sprintf (errmsg, "Error allocating memory for ozone uint8 pixels");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    wv_2d[0] = (uint16_t*) calloc (CMG_NLINES * CMG_NSAMPS, sizeof(uint16_t));
    if (wv_2d[0] == NULL)
    {
        sprintf (errmsg, "Error allocating memory for WV uint16 pixels");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    for (line = 1; line < CMG_NLINES; line++)
    {
        oz_2d[line] = oz_2d[0] + line * CMG_NSAMPS;
        wv_2d[line] = wv_2d[0] + line * CMG_NSAMPS;
    }

    /* Determine the number of pixels to be processed */
    n_pixels = CMG_NSAMPS * CMG_NLINES;

    /* Read the entire VIIRS ozone dataset */
    retval = H5Dread (ozone_dsid, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &oz_2d[0][0]);
    if (retval < 0)
    {
        sprintf (errmsg, "Error reading ozone dataset from file: %s",
            viirs_aux_file);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the entire VIIRS water vapor dataset */
    retval = H5Dread (wv_dsid, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &wv_2d[0][0]);
    if (retval < 0)
    {
        sprintf (errmsg, "Error reading water vapor dataset from file: %s",
            viirs_aux_file);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Determine the gapfill weighting for this day */
    determine_weights (aux_day, &prev_weight, &target_weight, &next_weight);
    printf ("Gapfill weights: Target: %.02f  Previous: %.02f  Next: %.02f\n",
        target_weight, prev_weight, next_weight);

    /* Read the data for each of the previous, target, and next monthly avgs
       for ozone and water vapor. Allocates memory for the arrays. Returns a
       NULL pointer in the array if the previous or next weight is zero. */
    retval = read_monthly_avgs (aux_month, aux_year, n_pixels, prev_weight,
        target_weight, next_weight, monthly_avg_oz_data, monthly_avg_wv_data);
    if (retval == ERROR)
    {
        sprintf (errmsg, "Unable to read the monthly averages needed for "
            "gapfilling the VIIRS file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Fill in the gaps using a weighted value from the monthly averages.  Use
       the ozone data to identify the fill/gap pixels, which are the same for
       both ozone and water vapor. */
    printf ("Gapfilling VIIRS products for WV and OZ ...\n");
    for (line = 0; line < CMG_NLINES; line++)
    {
        pix = line * CMG_NSAMPS;
        for (samp = 0; samp < CMG_NSAMPS; samp++, pix++)
        {
            /* If the pixel is not fill then continue */
            if (oz_2d[line][samp] != VIIRS_FILL)
                continue;

            /* Fill this pixel */
            target_oz = (float) monthly_avg_oz_data[1][pix];
            target_wv = (float) monthly_avg_wv_data[1][pix];

            prev_oz = 0.0;
            prev_wv = 0.0;
            if (prev_weight > 0.0)
            {
                prev_oz = (float) monthly_avg_oz_data[0][pix];
                prev_wv = (float) monthly_avg_wv_data[0][pix];
            }

            next_oz = 0.0;
            next_wv = 0.0;
            if (next_weight > 0.0)
            {
                next_oz = (float) monthly_avg_oz_data[2][pix];
                next_wv = (float) monthly_avg_wv_data[2][pix];
            }

            oz_2d[line][samp] = (uint8_t) get_fill_value (prev_weight,
                target_weight, next_weight, prev_oz, target_oz, next_oz);
            wv_2d[line][samp] = (uint16_t) get_fill_value (prev_weight,
                target_weight, next_weight, prev_wv, target_wv, next_wv);
        }  /* for line */
    }  /* for samp */

    /* Write the entire VIIRS ozone dataset */
    retval = H5Dwrite (ozone_dsid, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &oz_2d[0][0]);
    if (retval < 0)
    {
        sprintf (errmsg, "Error writing ozone dataset to file: %s",
            viirs_aux_file);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Write the entire VIIRS water vapor dataset */
    retval = H5Dwrite (wv_dsid, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &wv_2d[0][0]);
    if (retval < 0)
    {
        sprintf (errmsg, "Error writing water vapor dataset to file: %s",
            viirs_aux_file);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Close and clean up */
    free (viirs_aux_file);
    H5Fclose (file_id); 
    H5Dclose (wv_dsid);
    H5Dclose (ozone_dsid);

    /* Free the ozone and water vapor monthly average data */
    free (monthly_avg_oz_data[0]);
    free (monthly_avg_oz_data[1]);
    free (monthly_avg_oz_data[2]);
    free (monthly_avg_wv_data[0]);
    free (monthly_avg_wv_data[1]);
    free (monthly_avg_wv_data[2]);

    /* Free the VIIRS data arrays */
    free (oz_2d[0]);
    free (wv_2d[0]);
    free (oz_2d);
    free (wv_2d);

    /* Successful completion */
    exit (SUCCESS);
}


/******************************************************************************
MODULE:  open_oz_wv_datasets

PURPOSE:  Opens the daily, global VIIRS CMG file obtaining the ozone and
water vapor dataset IDs

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred opening the input
SUCCESS        Successful completion

NOTES:
******************************************************************************/
int open_oz_wv_datasets
(
    char *filename,       /* I: VIIRS file to be read */
    hid_t *file_id,       /* O: VIIRS file id */
    hid_t *ozone_dsid,    /* O: ozone dataset ID */
    hid_t *wv_dsid        /* O: water vapor dataset ID */
)
{
    char FUNC_NAME[] = "open_oz_wv_datasets"; /* function name */
    char errmsg[STR_SIZE];  /* error message */
    char dataset_name[STR_SIZE];  /* full path name of the dataset */
    int ndims;              /* number of dimensions in this dataset */
    hid_t datatype;         /* datatype id */
    hid_t dataspace;        /* dataspace id */
    H5T_class_t class;      /* datatype class */
    size_t size;            /* size of the data element stored in the file */
    hsize_t dims[2];        /* dimension of desired 2D datasets */

    /* Open the input file for reading and writing */
    *file_id = H5Fopen (filename, H5F_ACC_RDWR, H5P_DEFAULT);
    if (*file_id < 0)
    {
        sprintf (errmsg, "Error opening file: %s", filename);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the ozone dataset */
    sprintf (dataset_name, "%s%s", dataset_path, list_of_datasets[OZONE]);
    *ozone_dsid = H5Dopen (*file_id, dataset_name, H5P_DEFAULT);
    if (*ozone_dsid < 0)
    {
        sprintf (errmsg, "Error opening the ozone dataset: %s", dataset_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    
    /* Get datatype and dataspace handles and then query dataset class, order,
       size, and dimensions. Data type should be a U8-bit integer. */
    datatype = H5Dget_type (*ozone_dsid);     /* datatype handle */ 
    class = H5Tget_class (datatype);
    if (class != H5T_INTEGER)
    {
        sprintf (errmsg, "Unexpected data type of ozone dataset: %d (should be "
            "H5T_INTEGER)", class);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    size = H5Tget_size (datatype);
    if ((int)size != 1)
    {
        sprintf (errmsg, "Unexpected data type of ozone dataset: %d (should be "
            "1 byte)", (int)size);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    H5Tclose (datatype);

    dataspace = H5Dget_space (*ozone_dsid);    /* dataspace handle */
    ndims = H5Sget_simple_extent_dims (dataspace, dims, NULL);
    if (ndims != CMG_NDIMS)
    {
        sprintf (errmsg, "Unexpected number of dimensions for the ozone "
            "dataset: %d (should be %d)", ndims, CMG_NDIMS);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    H5Sclose (dataspace);

    /* Verify the dimensions are as expected for the CMG grid */
    if (dims[0] != CMG_NLINES || dims[1] != CMG_NSAMPS)
    {
        sprintf (errmsg, "Unexpected size of ozone dataset: %d x %d (should be "
            "%d x %d)", (int)dims[0], (int)dims[1], CMG_NLINES, CMG_NSAMPS);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the water vapor dataset */
    sprintf (dataset_name, "%s%s", dataset_path, list_of_datasets[WV]);
    *wv_dsid = H5Dopen (*file_id, dataset_name, H5P_DEFAULT);
    if (*wv_dsid < 0)
    {
        sprintf (errmsg, "Error opening water vapor dataset: %s", dataset_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Get datatype and dataspace handles and then query dataset class, order,
       size, and dimensions. Date type should be a U16-bit integer. */
    datatype = H5Dget_type (*wv_dsid);     /* datatype handle */ 
    class = H5Tget_class (datatype);
    if (class != H5T_INTEGER)
    {
        sprintf (errmsg, "Unexpected data type of water vapor dataset: %d "
            "(should be H5T_INTEGER)", class);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    size = H5Tget_size (datatype);
    if ((int)size != 2)
    {
        sprintf (errmsg, "Unexpected data type of ozone dataset: %d (should be "
            "2 bytes)", (int)size);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    H5Tclose (datatype);

    dataspace = H5Dget_space (*wv_dsid);    /* dataspace handle */
    ndims = H5Sget_simple_extent_dims (dataspace, dims, NULL);
    if (ndims != CMG_NDIMS)
    {
        sprintf (errmsg, "Unexpected number of dimensions for the water vapor "
            "dataset: %d (should be %d)", ndims, CMG_NDIMS);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    H5Sclose (dataspace);

    /* Verify the dimensions are as expected for the CMG grid */
    if (dims[0] != CMG_NLINES || dims[1] != CMG_NSAMPS)
    {
        sprintf (errmsg, "Unexpected size of water vapor dataset: %d x %d "
            "(should be %d x %d)", (int)dims[0], (int)dims[1], CMG_NLINES,
            CMG_NSAMPS);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

NOTES:
******************************************************************************/
void usage ()
{
    printf ("gapfill_viirs_aux reads the ozone and water vapor datasets from "
            "the VIIRS auxiliary data, fills the gaps using monthly "
            "climatology averages, and writes the new data back out to the "
            "HDF5 file.\n\n");
    printf ("usage: gapfill_viirs_aux --viirs_aux=input_viirs_aux_filename "
            "--month=month_of_aux_file --day=day_of_month_of_aux_file "
            "--year=year_of_aux_file\n");

    printf ("\nwhere the following parameters are required:\n");
    printf ("    -viirs_aux: name of the input VIIRS auxiliary file (VNP04ANC "
            "or VJ104ANC) to be processed. The ozone and water vapor datasets "
            "will be modified with the gapfilled data.\n"
            "    -month: month (1-12) of the auxiliary file\n"
            "    -day: day of month (1-31) of the auxiliary file\n"
            "    -year: year of the auxiliary file\n\n");

    printf ("\ngapfill_viirs_aux --help will print the usage statement\n");
}
