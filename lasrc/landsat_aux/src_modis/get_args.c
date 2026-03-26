#include <getopt.h>
#include "combine_l8_aux_data.h"

/******************************************************************************
MODULE:  get_args

PURPOSE:  Gets the command-line arguments and validates that the required
arguments were specified.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error getting the command-line arguments or a command-line
                argument and associated value were not specified
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date          Programmer       Reason
----------    ---------------  -------------------------------------
8/26/2014     Gail Schmidt     Original Development
9/3/2014      Gail Schmidt     Added an output directory option as a cmd-line
                               option for the user

NOTES:
  1. Memory is allocated for the input files.  This should be character a
     pointer set to NULL on input.  The caller is responsible for freeing the
     allocated memory upon successful return.
******************************************************************************/
int get_args
(
    int argc,               /* I: number of cmd-line args */
    char *argv[],           /* I: string of cmd-line args */
    char **terra_cmg_file,  /* O: address of input Terra CMG file */
    char **aqua_cmg_file,   /* O: address of input Aqua CMG file */
    char **terra_cma_file,  /* O: address of input Terra CMA file */
    char **aqua_cma_file,   /* O: address of input Aqua CMA file */
    char **output_dir,      /* O: address of output directory */
    bool *verbose           /* O: verbose flag */
)
{
    int c;                           /* current argument index */
    int option_index;                /* index for the command-line option */
    static int verbose_flag=0;       /* verbose flag */
    char errmsg[STR_SIZE];           /* error message */
    char FUNC_NAME[] = "get_args";   /* function name */
    static struct option long_options[] =
    {
        {"verbose", no_argument, &verbose_flag, 1},
        {"terra_cmg", required_argument, 0, 'a'},
        {"aqua_cmg", required_argument, 0, 'b'},
        {"terra_cma", required_argument, 0, 'c'},
        {"aqua_cma", required_argument, 0, 'd'},
        {"output_dir", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /* Initialize the flags to false */
    *verbose = false;

    /* Loop through all the cmd-line options */
    opterr = 0;   /* turn off getopt_long error msgs as we'll print our own */
    while (1)
    {
        /* optstring in call to getopt_long is empty since we will only
           support the long options */
        c = getopt_long (argc, argv, "", long_options, &option_index);
        if (c == -1)
        {   /* Out of cmd-line options */
            break;
        }

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
     
            case 'h':  /* help */
                usage ();
                return (ERROR);
                break;

            case 'a':  /* Terra CMG input file */
                *terra_cmg_file = strdup (optarg);
                break;
     
            case 'b':  /* Aqua CMG input file */
                *aqua_cmg_file = strdup (optarg);
                break;
     
            case 'c':  /* Terra CMA input file */
                *terra_cma_file = strdup (optarg);
                break;
     
            case 'd':  /* Aqua CMA input file */
                *aqua_cma_file = strdup (optarg);
                break;
     
            case 'o':  /* Output directory */
                *output_dir = strdup (optarg);
                break;
     
            case '?':
            default:
                sprintf (errmsg, "Unknown option %s", argv[optind-1]);
                error_handler (true, FUNC_NAME, errmsg);
                usage ();
                return (ERROR);
                break;
        }
    }

    /* Make sure the Terra/Aqua CMG/CMA files were specified */
    if (*terra_cmg_file == NULL && *aqua_cmg_file == NULL)
    {
        sprintf (errmsg, "Input Terra CMG or Aqua CMG file is a required "
            "argument. At least one must be valid.");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    if (*terra_cma_file == NULL && *aqua_cma_file == NULL)
    {
        sprintf (errmsg, "Input Terra CMA or Aqua CMA file is a required "
            "argument. At least one must be valid.");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    /* If Terra CMA is available, then Terra CMG should be available */
    if (((*terra_cma_file != NULL) && (*terra_cmg_file == NULL)) ||
        ((*terra_cma_file == NULL) && (*terra_cmg_file != NULL)))
    {
        sprintf (errmsg, "Input Terra CMA and CMG files are required "
            "to be valid. If one is valid the other must also be valid.");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    /* If Aqua CMA is available, then Aqua CMG should be available */
    if (((*aqua_cma_file != NULL) && (*aqua_cmg_file == NULL)) ||
        ((*aqua_cma_file == NULL) && (*aqua_cmg_file != NULL)))
    {
        sprintf (errmsg, "Input Aqua CMA and CMG files are required "
            "to be valid. If one is valid the other must also be valid.");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    /* Make sure the output directory was specified */
    if (*output_dir == NULL)
    {
        sprintf (errmsg, "Output directory is a required argument");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    /* Check the flags */
    if (verbose_flag)
        *verbose = true;

    return (SUCCESS);
}

