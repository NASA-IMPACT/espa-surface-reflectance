/*****************************************************************************
FILE: lndpm.c
  
PURPOSE: Contains functions for reading the input XML metadata file and
creating the parameter and metadata files needed for downstream processing
by LEDAPS applications.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
  1. The XML metadata format written via this library follows the ESPA internal
     metadata format found in ESPA Raw Binary Format v1.0.doc.  The schema for
     the ESPA internal metadata format is available at
     http://espa.cr.usgs.gov/schema/espa_internal_metadata_vX_Y.xsd.
*****************************************************************************/
#include <sys/stat.h>
#include "lndpm.h"

int conv_date (int *mm, int *dd, int yyyy);
int find_file (char *path, char *name);
int get_args (int argc, char *argv[], char **xml_infile, bool *process_sr);
void usage ();


/******************************************************************************
MODULE:  main (lndpm)

PURPOSE: Create the parameter files that are needed by the LEDAPS applications
(lndcal, lndsr)

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error with the parameter file generation
SUCCESS         Successfully generated the parameter files

NOTES:
******************************************************************************/
int main (int argc, char *argv[])
{
    char FUNC_NAME[] = "lndpm";    /* function name */
    char errmsg[STR_SIZE];         /* error message */
    char tmpfile[STR_SIZE];        /* temporary storage of the XML file */
    char scene_name[STR_SIZE];     /* name of the scene */
    char lndcal_name[STR_SIZE];    /* name of the lndcal input file */
    char lndsr_name[STR_SIZE];     /* name of the lndsr input file */
    char dem[STR_SIZE];            /* name of DEM file */
    char ozone[STR_SIZE];          /* name of ozone file */
    char reanalysis[STR_SIZE];     /* name of NCEP file */
    char path_buf[DIR_BUF_SIZE];   /* path to the auxiliary/cal file */
    char *xml_infile = NULL;       /* input XML filename */
    char *aux_path = NULL;         /* path for LEDAPS auxiliary data */
    char *token_ptr = NULL;        /* pointer used for obtaining scene name */
    char *file_ptr = NULL;         /* pointer used for obtaining file name */
    int year, month, day;          /* year, month, day of acquisition date */
    int retval;                    /* return status */
    bool anc_missing = false;      /* is the auxiliary data missing? */
    bool process_sr = true;        /* specifies if surface reflectance
                                      processing will be completed (true) or if
                                      only TOA processing will be run (false) */
    FILE *out = NULL;              /* pointer to the output parameter file */
    Espa_internal_meta_t xml_metadata;  /* XML metadata structure */

    printf ("\nRunning lndpm ...\n");

    /* Check the command-line arguments */
    retval = get_args (argc, argv, &xml_infile, &process_sr);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        exit (ERROR);
    }

    /* Validate the input metadata file */
    if (validate_xml_file (xml_infile) != SUCCESS)
    {  /* Error messages already written */
        return (ERROR);
    }

    /* Initialize the metadata structure */
    init_metadata_struct (&xml_metadata);

    /* Parse the metadata file into our internal metadata structure; also
       allocates space as needed for various pointers in the global and band
       metadata */
    if (parse_metadata (xml_infile, &xml_metadata) != SUCCESS)
    {  /* Error messages already written */
        return (ERROR);
    }

    /* Get the scene name.  Strip off the path and file extension as it's
       assumed the XML filename is the scene name followed by the extension. */
    strcpy (tmpfile, xml_infile);
    file_ptr = strrchr (tmpfile, '/');
    if (file_ptr != NULL)
        file_ptr++;
    else
        file_ptr = tmpfile;
    token_ptr = strtok (file_ptr, ".");
    sprintf (scene_name, "%s", token_ptr);
  
    /* Set up the names of the input files for downstream processing */
    sprintf (lndcal_name, "lndcal.%s.txt", scene_name);
    sprintf (lndsr_name, "lndsr.%s.txt", scene_name);

    /* Open the parameter file for lndcal for writing */
    out = fopen (lndcal_name, "w");
    if (out == NULL)
    {
        sprintf (errmsg, "Opening lndcal parameter file for writing: %s",
            lndcal_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Get the path for the LEDAPS auxiliary products (NCEP, TOMS, DEM, etc.)
       from the LEDAPS_AUX_DIR environment variable.  If it isn't defined, then
       assume the products are in the local directory. */
    aux_path = getenv ("LEDAPS_AUX_DIR");
    if (aux_path == NULL)
    {
        aux_path = ".";
        sprintf (errmsg, "LEDAPS_AUX_DIR environment variable isn't defined. "
            "It is assumed the LEDAPS auxiliary products will be available "
            "from the local directory.");
        error_handler (false, FUNC_NAME, errmsg);
    }

    /* Write the parameter data to the lndcal parameter file */
    fprintf (out, "PARAMETER_FILE\n");
    fprintf (out, "XML_FILE = %s\n", xml_infile);
    fprintf (out, "LEDAPSVersion = \"%s\"\n", LEDAPS_VERSION);
    fprintf (out, "END\n");
    fclose (out);

    /* Get year, month, and day from acquisition date */
    token_ptr = strtok (xml_metadata.global.acquisition_date, "-");
    sscanf (token_ptr, "%d", &year);
    token_ptr = strtok (NULL, "-");
    sscanf (token_ptr, "%d", &month);
    token_ptr = strtok (NULL, "-");
    sscanf (token_ptr, "%d", &day);
  
    /* Convert to day of year */
    if (conv_date (&month, &day, year) != SUCCESS)
    {
        sprintf (errmsg, "Not able to convert the month, day, year from the "
            "acquisition date to Julian DOY.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find and prepare auxiliary files.  First look in the expected location.
       If not there, then traverse the auxiliary directory structure. */
    /* DEM file */
    strcpy (dem, "CMGDEM.hdf");
    char full_path[DIR_BUF_SIZE];
    sprintf(full_path, "%s/%s", aux_path, dem);
    if (find_file(full_path, NULL))
    {
        strcpy (dem, full_path);
        printf ("using DEM : %s\n", dem);
    }
    else
    {
        sprintf (errmsg, "Could not find DEM auxiliary data: %s\n  Check "
            "LEDAPS_AUX_DIR environment variable.", dem);
        error_handler (false, FUNC_NAME, errmsg);
        anc_missing = true;
    }

    /* TOMS ozone file */
    sprintf (ozone, "TOMS_%d%03d.hdf", year, day);
    strcpy (path_buf, aux_path);
    sprintf(full_path, "%s/EP_TOMS/ozone_%d/%s", aux_path, year, ozone);
    if (find_file(full_path, NULL))
    {
        strcpy (ozone, full_path);
        printf ("using TOMS: %s\n", ozone);
    }
    else if (find_file (path_buf, ozone))
    {
        strcpy (ozone, path_buf);
        printf ("using TOMS : %s\n", ozone);
    }
    else
    {
        sprintf (errmsg, "Could not find TOMS auxiliary data: %s\n  Check "
            "LEDAPS_AUX_DIR environment variable.", ozone);
        error_handler (false, FUNC_NAME, errmsg);
        anc_missing = true;
    }
    
    /* NCEP file */
    sprintf (reanalysis, "REANALYSIS_%d%03d.hdf", year, day);
    strcpy (path_buf, aux_path);
    sprintf(full_path, "%s/REANALYSIS/RE_%d/%s", aux_path, year, reanalysis);
    if (find_file(full_path, NULL))
    {
        strcpy (reanalysis, full_path);
        printf ("using REANALYSIS : %s\n", reanalysis);
    }
    else if (find_file (path_buf, reanalysis))
    {
        strcpy (reanalysis, path_buf);
        printf ("using REANALYSIS : %s\n", reanalysis);
    }
    else
    {
        sprintf (errmsg, "Could not find NCEP REANALYSIS auxiliary data: %s\n"
            "  Check LEDAPS_AUX_DIR environment variable.", reanalysis);
        error_handler (false, FUNC_NAME, errmsg);
        anc_missing = true;
    }

    /* If processing SR, check to see if missing auxiliary data.  Auxiliary
       data is not used for TOA-only processing. */
    if (anc_missing && process_sr)
    {
        sprintf (errmsg, "Verify the missing auxiliary data products, then "
            "try reprocessing.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the parameter file for lndsr for writing */
    out = fopen (lndsr_name, "w");
    if (out == NULL)
    {
        sprintf (errmsg, "Opening lndsr parameter file for writing: %s",
            lndsr_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    fprintf (out, "PARAMETER_FILE\n");
    fprintf (out, "XML_FILE = %s\n", xml_infile);
    fprintf (out, "DEM_FILE = %s\n", dem);
    if (fopen (ozone, "r") != NULL)
    {
        /* if ozone file doesn't exist then don't write it to the parameter file
           and instead use climatology estimation */
        fprintf (out, "OZON_FIL = %s\n", ozone);
    }
    fprintf (out, "PRWV_FIL = %s\n", reanalysis);
    fprintf (out, "LEDAPSVersion = \"%s\"\n", LEDAPS_VERSION);
    fprintf (out, "END\n");
    fclose (out);

    /* Free the metadata structure and pointers */
    free_metadata (&xml_metadata);
    free (xml_infile);

    /* Successful completion */
    printf ("lndpm complete.\n");
    return (SUCCESS);
}


/******************************************************************************
MODULE:  conv_date

PURPOSE: Convert year-mm-dd to Julian day of year (DOY) or vice versa.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error with the conversion
SUCCESS         Successfully converted the date

NOTES:
  1. *mm and *dd are input/output varialbes
     Input julian day number:
        input *mm = 0;
        input *dd = julian day
        output in mm-dd-yyyy format
     Input in mm-dd-yyyy format
        output *mm = 0;
        output *dd = julian day number
******************************************************************************/
int conv_date
(
    int *mm,         /* I/O: month value */
    int *dd,         /* I/O: day of month or Julian day value */
    int yyyy         /* I: year value */
)
{
    char FUNC_NAME[] = "conv_date";    /* function name */
    char errmsg[STR_SIZE];             /* error message */
    int nm, im;        /* looping vars */
    int ndays[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
                       /* number of days in each month */

    /* Check for leap year */
    if ((yyyy%400 == 0) || ((yyyy%4 == 0) && (yyyy%100 != 0))) ndays[1] = 29; 

    /* Do the conversion */
    if (*mm == 0)
    { /* Convert from Julian DOY to mm-dd-yyyy */
        if (*dd <= 0) {
            sprintf (errmsg, "Invalid input date: %d %d", *dd, yyyy);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        else
        {
            for (im = 0; ((im < 12) && (*dd > 0)); im++)
                *dd -= ndays[im];
            if ((im > 12) || ((im == 12) && (*dd > 0)))
            {
                sprintf (errmsg, "Invalid input date: %d %d", *dd, yyyy);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            else
            {
                *mm = im;
                *dd += ndays[*mm - 1];
            }
        }
    }
    else
    { /* Convert from mm-dd-yyyy to Julian DOY */
        if ((*mm <= 0) || (*dd <= 0))
        {
            sprintf (errmsg, "Invalid input date: %d-%d-%d", *mm, *dd, yyyy);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        else
        {
            nm = *mm - 1;
            for (im = 0; im < nm; im++)
                *dd += ndays[im];
            *mm = 0;
        }
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  scan_dir

PURPOSE: Search the specified directory, recursively, for the specified
filename.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
non-zero        File is found, and path points to full path
zero            File is not found, and path is unchanged

NOTES:
******************************************************************************/
int scan_dir
(
    char *path,   /* I: directory/path to search; upon successfully finding
                        the file, path contains the location of the file */
    char *name    /* I: filename for which to search */
)
{
    char FUNC_NAME[] = "scan_dir";    /* function name */
    char errmsg[STR_SIZE];             /* error message */
    DIR *fd;                  /* pointer to directory */
    struct dirent *dirent_p;  /* pointer to directory */
    char *nbp;                /* pointer to the end of the path */
    int found = 0;            /* was the file found? */
    
    /* Add directory separator to end of directory name if not already there */
    nbp = path + strlen (path);
    if (*nbp != '/')
        *nbp++ = '/';
      
    if (nbp+MAXNAMLEN+2 >= path + DIR_BUF_SIZE)
    {
        sprintf (errmsg, "Path name too long -- cannot search: %s", path);
        error_handler (true, FUNC_NAME, errmsg);
        return (found);
    }

    /* Check the directory to see if it's accessible */
    if ((fd = opendir (path)) == NULL)
    {
        sprintf (errmsg, "Could not read directory: %s", path);
        error_handler (true, FUNC_NAME, errmsg);
        return (found);
    }

    /* Recursively search for the file while traversing through the directory
       structure */
    while ((dirent_p = readdir(fd)) != NULL)   /* search directory */
    {
        if (dirent_p->d_ino == 0)                /* slot not in use */
            continue;
        if (strcmp (dirent_p->d_name, ".") == 0   /* ignore current ... */
        || strcmp (dirent_p->d_name, "..") == 0)  /* and parent directory */
            continue;
        
        strcat (path, dirent_p->d_name);          /* check this path */
        if ((found = find_file (path, name)) != 0)
            break;                                /* found it */
        else
            *nbp = '\0';                          /* restore directory name */
    }

    /* Close the directory pointer and successful completion */
    closedir (fd);
    return (found);
}


/******************************************************************************
MODULE:  find_file

PURPOSE: Look in the current directory for the specified file.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
non-zero        File is found, and path points to full path
zero            File is not found, and path is unchanged

NOTES:
******************************************************************************/
int find_file
(
    char *path,   /* I: directory/path to search; upon successfully finding
                        the file, path contains the location of the file */
    char *name    /* I: filename for which to search
                        If the value is NULL, then path is assumed to be a
                        directory and filename. */
)
{
    char FUNC_NAME[] = "find_file";    /* function name */
    char errmsg[STR_SIZE];             /* error message */
    struct stat stbuf;                 /* buffer for file/directory stat */
    char pbuf[DIR_BUF_SIZE];           /* path buffer */
    int found = 0;                     /* was the file found? */
    
    /* Make sure the path exists */
    if (stat (path, &stbuf) != 0)
    {
        if (name != NULL)
            sprintf(errmsg, "Can't stat directory: %s", path);
        else
            sprintf(errmsg, "Can't stat file: %s", path);
        error_handler (true, FUNC_NAME, errmsg);
        return 0;
    }

    /* If this is a directory, then search it.  Otherwise it's a file. */
    if ((stbuf.st_mode & S_IFMT) == S_IFDIR)
    {
        if (name == NULL)
        {
            error_handler(true, FUNC_NAME, "Filename not specified.");
            return 0;
        }

        strcpy(pbuf, path);
        found = scan_dir (pbuf, name);

        /* If the file was found return the location. */
        if (found)
            strcpy(path, pbuf);
    }
    else if (name != NULL)
    {
        sprintf(errmsg, "Search filename (%s) specified when path (%s) "
                "includes filename.", name, path);
        error_handler(true, FUNC_NAME, errmsg);
        return 0;
    }
    else
        found = 1;

    return found;
}


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

NOTES:
  1. The input files should be character a pointer set to NULL on input. Memory
     for these pointers is allocated by this routine. The caller is responsible
     for freeing the allocated memory upon successful return.
******************************************************************************/
int get_args
(
    int argc,             /* I: number of cmd-line args */
    char *argv[],         /* I: string of cmd-line args */
    char **xml_infile,    /* O: address of input XML file */
    bool *process_sr      /* O: process the surface reflectance products */
)
{
    int c;                           /* current argument index */
    int option_index;                /* index for the command-line option */
    char errmsg[STR_SIZE];           /* error message */
    char FUNC_NAME[] = "get_args";   /* function name */
    static struct option long_options[] =
    {
        {"xml", required_argument, 0, 'i'},
        {"process_sr", required_argument, 0, 'p'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /* Initialize the flags to false */
    *process_sr = true;    /* default is to process SR products */

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

            case 'i':  /* XML input file */
                *xml_infile = strdup (optarg);
                break;

            case 'p':  /* process SR products */
                if (!strcmp (optarg, "true"))
                    *process_sr = true;
                else if (!strcmp (optarg, "false"))
                    *process_sr = false;
                else
                {
                    sprintf (errmsg, "Unknown value for process_sr: %s",
                        optarg);
                    error_handler (true, FUNC_NAME, errmsg);
                    usage ();
                    return (ERROR);
                }
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

    /* Make sure the XML file was specified */
    if (*xml_infile == NULL)
    {
        sprintf (errmsg, "Input XML file is a required argument");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

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
    printf ("lndpm sets up the parameter file for the LEDAPS processing.\n\n");
    printf ("usage: lndpm "
            "--xml=input_xml_filename --process_sr=true:false\n");

    printf ("\nwhere the following parameters are required:\n");
    printf ("    -xml: name of the input XML file to be processed\n");
    printf ("    -process_sr: the default is to process surface reflectance, "
            "however if this flag is set to false then only the TOA "
            "reflectance processing and brightness temperature will be "
            "done.\n");

    printf ("\nlndpm --help will print the usage statement\n");
    printf ("\nExample: lndpm --xml=LC80410272013181LGN00.xml "
            "--process_sr=true\n");
    printf ("   ==> Sets up the parameter file and checks for the auxliary "
            "data files.\n\n");

    printf ("\nExample: lndpm --xml=LC80410272013181LGN00.xml "
            "--process_sr=false\n");
    printf ("   ==> Sets up the paramter file and doesn't check for the "
            "auxiliary data files since only TOA/BT processing will be "
            "run.\n\n");
}

