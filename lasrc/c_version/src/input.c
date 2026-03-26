/*****************************************************************************
FILE: input.c
  
PURPOSE: Contains functions for handling of the input data files for this
application.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/

#include "input.h"

/******************************************************************************
MODULE:  open_input

PURPOSE:  Sets up the input data structure, opens the input reflectance file
for read access, allocates space, and stores some of the metadata for later
reference.

RETURN VALUE:
Type = Input_t*
Value      Description
-----      -----------
NULL       Error occurred opening or reading the file
non-NULL   Successful completion

NOTES:
  1. This routine opens the input Landsat 8/9 files. It also allocates memory
     for pointers in the input structure.  It is up to the caller to use
     close_input and free_input to close the files and free up the memory when
     done using the input data structure.
******************************************************************************/
Input_t *open_input
(
    Espa_internal_meta_t *metadata,   /* I: input metadata */
    bool use_orig_aero,               /* I: use the original aerosol handling
                                            if specified, o/w use the
                                            semi-empirical approach */
    bool process_sr                   /* I: will SR data be processed? */
)
{
    char FUNC_NAME[] = "open_input";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    Input_t *this = NULL;     /* input data structure to be initialized,
                                 populated, and returned to the caller */
    int ib;                   /* loop counter for bands */

    /* Create the Input data structure */
    this = malloc (sizeof (Input_t));
    if (this == NULL) 
    {
        strcpy (errmsg, "Error allocating memory for Input data structure");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Initialize and get input from metadata file */
    if (get_xml_input (metadata, process_sr, use_orig_aero, this) != SUCCESS)
    {
        strcpy (errmsg, "Error getting input information from the metadata "
            "file.");
        free (this);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Make sure the metadata satellite is either Landsat 8/9 or Sentinel-2 */
    if (this->meta.sat != SAT_LANDSAT_8 && this->meta.sat != SAT_LANDSAT_9 &&
        this->meta.sat != SAT_SENTINEL_2)
    {
        strcpy (errmsg, "Error getting satellite information from the input "
            "metadata file. Only Landsat 8/9 and Sentinel 2 are supported.");
        free (this);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Open files for access */
    for (ib = 0; ib < this->nband; ib++)
    {
        this->fp_bin[ib] = open_raw_binary (this->file_name[ib], "rb");
        if (this->fp_bin[ib] == NULL)
        {
            sprintf (errmsg, "Opening reflectance raw binary file: %s",
                this->file_name[ib]);
            error_handler (true, FUNC_NAME, errmsg);
            free_input (this, use_orig_aero);
            return (NULL);
        }
        this->open[ib] = true;
    }

    for (ib = 0; ib < this->nband_th; ib++)
    {  /* NOTE: nband_th will be 0 for OLI-only scenes */
        this->fp_bin_th[ib] = open_raw_binary (this->file_name_th[ib], "rb");
        if (this->fp_bin_th[ib] == NULL)
        {
            sprintf (errmsg, "Opening thermal raw binary file: %s",
                this->file_name_th[ib]);
            error_handler (true, FUNC_NAME, errmsg);
            free_input (this, use_orig_aero);
            return (NULL);
        }
        this->open_th[ib] = true;
    }

    for (ib = 0; ib < this->nband_pan; ib++)
    {
        this->fp_bin_pan[ib] = open_raw_binary (this->file_name_pan[ib], "rb");
        if (this->fp_bin_pan[ib] == NULL)
        {
            sprintf (errmsg, "Opening pan raw binary file: %s",
                this->file_name_pan[ib]);
            error_handler (true, FUNC_NAME, errmsg);
            free_input (this, use_orig_aero);
            return (NULL);
        }
        this->open_pan[ib] = true;
    }

    for (ib = 0; ib < this->nband_qa; ib++)
    {
        this->fp_bin_qa[ib] = open_raw_binary (this->file_name_qa[ib], "rb");
        if (this->fp_bin_qa[ib] == NULL)
        {
            sprintf (errmsg, "Opening QA raw binary file: %s",
                this->file_name_qa[ib]);
            error_handler (true, FUNC_NAME, errmsg);
            free_input (this, use_orig_aero);
            return (NULL);
        }
        this->open_qa[ib] = true;
    }

    /* Open the per-pixel solar and view angle bands for Landsat */
    if (this->meta.sat == SAT_LANDSAT_8 || this->meta.sat == SAT_LANDSAT_9)
    {
        this->fp_bin_sza = open_raw_binary (this->file_name_sza, "rb");
        if (this->fp_bin_sza == NULL)
        {
            sprintf (errmsg, "Opening solar zenith raw binary file: %s",
                this->file_name_sza);
            error_handler (true, FUNC_NAME, errmsg);
            free_input (this, use_orig_aero);
            return (NULL);
        }
        this->open_ppa = true;

        /* Only open the rest of these bands if we are using the original
           aerosol algorithm */
        if (use_orig_aero)
        {
            this->fp_bin_saa = open_raw_binary (this->file_name_saa, "rb");
            if (this->fp_bin_saa == NULL)
            {
                sprintf (errmsg, "Opening solar zenith raw binary file: %s",
                    this->file_name_saa);
                error_handler (true, FUNC_NAME, errmsg);
                free_input (this, use_orig_aero);
                return (NULL);
            }
        
            this->fp_bin_vza = open_raw_binary (this->file_name_vza, "rb");
            if (this->fp_bin_vza == NULL)
            {
                sprintf (errmsg, "Opening view zenith raw binary file: %s",
                    this->file_name_vza);
                error_handler (true, FUNC_NAME, errmsg);
                free_input (this, use_orig_aero);
                return (NULL);
            }
        
            this->fp_bin_vaa = open_raw_binary (this->file_name_vaa, "rb");
            if (this->fp_bin_vaa == NULL)
            {
                sprintf (errmsg, "Opening view zenith raw binary file: %s",
                    this->file_name_vaa);
                error_handler (true, FUNC_NAME, errmsg);
                free_input (this, use_orig_aero);
                return (NULL);
            }
        }
    }

    /* Do a cursory check to make sure the bands and QA band exist and have
       been opened */
    if (!this->open[0])
    {
        sprintf (errmsg, "Reflectance band 1 is not open.");
        error_handler (true, FUNC_NAME, errmsg);
        free_input (this, use_orig_aero);
        return (NULL);
    }

    if (this->nband_th != 0 && !this->open_th[0])
    {
        sprintf (errmsg, "Thermal band 10 is not open.");
        error_handler (true, FUNC_NAME, errmsg);
        free_input (this, use_orig_aero);
        return (NULL);
    }

    /* Landsat should have a QA input band */
    if ((this->meta.sat == SAT_LANDSAT_8 || this->meta.sat == SAT_LANDSAT_9) &&
         !this->open_qa[0])
    {
        sprintf (errmsg, "Landsat QA band is not open.");
        error_handler (true, FUNC_NAME, errmsg);
        free_input (this, use_orig_aero);
        return (NULL);
    }

    /* Landsat should have a per-pixel angle input bands */
    if ((this->meta.sat == SAT_LANDSAT_8 || this->meta.sat == SAT_LANDSAT_9) &&
         !this->open_ppa)
    {
        sprintf (errmsg, "Landsat per-pixel angle bands are not open.");
        error_handler (true, FUNC_NAME, errmsg);
        free_input (this, use_orig_aero);
        return (NULL);
    }

    return this;
}


/******************************************************************************
MODULE:  close_input

PURPOSE:  Ends SDS access and closes the input file.

RETURN VALUE:
Type = None

NOTES:
******************************************************************************/
void close_input
(
    Input_t *this,       /* I: pointer to input data structure */
    bool use_orig_aero   /* I: use the original aerosol handling if specified,
                               o/w use the semi-empirical approach */
)
{
    int ib;      /* loop counter for bands */
  
    /* Close the reflectance files */
    for (ib = 0; ib < this->nband; ib++)
    {
        if (this->open[ib])
        {
            close_raw_binary (this->fp_bin[ib]);
            this->open[ib] = false;
        }
    }

    /* Landsat has thermal, pan, QA, and per-pixel angle bands to close */
    if (this->meta.sat == SAT_LANDSAT_8 || this->meta.sat == SAT_LANDSAT_9)
    {
        /* Close the thermal files */
        for (ib = 0; ib < this->nband_th; ib++)
        {
            if (this->open_th[ib])
            {
                close_raw_binary (this->fp_bin_th[ib]);
                this->open_th[ib] = false;
            }
        }
    
        /* Close the pan files */
        for (ib = 0; ib < this->nband_pan; ib++)
        {
            if (this->open_pan[ib])
            {
                close_raw_binary (this->fp_bin_pan[ib]);
                this->open_pan[ib] = false;
            }
        }
    
        /* Close the QA files */
        for (ib = 0; ib < this->nband_qa; ib++)
        {
            if (this->open_qa[ib])
            {
                close_raw_binary (this->fp_bin_qa[ib]);
                this->open_qa[ib] = false;
            }
        }
    
        /* Close the per-pixel angle band files */
        if (this->open_ppa)
        {
            close_raw_binary (this->fp_bin_sza);
            if (use_orig_aero)
            {
                close_raw_binary (this->fp_bin_saa);
                close_raw_binary (this->fp_bin_vza);
                close_raw_binary (this->fp_bin_vaa);
            }
            this->open_ppa = false;
        }
    }
}


/******************************************************************************
MODULE:  free_input

PURPOSE:  Frees memory in the input data structure.

RETURN VALUE:
Type = None

NOTES:
******************************************************************************/
void free_input
(
    Input_t *this,      /* I: pointer to input data structure */
    bool use_orig_aero  /* I: use the original aerosol handling if specified,
                              o/w use the semi-empirical approach */
)
{
    char FUNC_NAME[] = "free_input";   /* function name */
    char errmsg[STR_SIZE];             /* error message */
    int ib;                            /* loop counter for bands */
   
    if (this != NULL)
    {
        if (this->open[0] || this->open_th[0] || this->open_pan[0]) 
        {
            strcpy (errmsg, "Freeing input data structure, but reflectance, "
                "thermal, and/or pan file(s) is/are still open. Use "
                "close_input to close the file");
            error_handler (false, FUNC_NAME, errmsg);
        }
  
        /* Free image band pointers */
        for (ib = 0; ib < this->nband; ib++)
            free (this->file_name[ib]);

        /* Landsat has thermal, pan, QA, and per-pixel angle bands to close */
        if (this->meta.sat == SAT_LANDSAT_8 || this->meta.sat == SAT_LANDSAT_9)
        {
            for (ib = 0; ib < this->nband_th; ib++)
                free (this->file_name_th[ib]);
            for (ib = 0; ib < this->nband_pan; ib++)
                free (this->file_name_pan[ib]);
            for (ib = 0; ib < this->nband_qa; ib++)
                free (this->file_name_qa[ib]);
            free(this->file_name_sza);
            if (use_orig_aero)
            {
                free(this->file_name_saa);
                free(this->file_name_vza);
                free(this->file_name_vaa);
            }
        }

        /* Free the data structure */
        free (this);
    } /* end if */
}


/******************************************************************************
MODULE:  get_input_refl_lines

PURPOSE:  Reads the reflectance data for the current refl band and lines, and
populates the output buffer.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading data for this band
SUCCESS    Successful completion

NOTES:
  1. The Input_t data structure needs to be populated and memory allocated
     before calling this routine.  Use open_input to do that.
******************************************************************************/
int get_input_refl_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current refl band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    int nsamps,      /* I: number of samples to read (Sentinel nsamps vary
                           depending on the band); if -99 then use the nsamps
                           in the input structure */
    uint16 *out_arr  /* O: output array to populate */
)
{
    char FUNC_NAME[] = "get_input_refl_line";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    long loc;                 /* current location in the input file */
  
    /* Check the parameters */
    if (this == NULL) 
    {
        strcpy (errmsg, "Input structure has not been opened/initialized");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (!this->open[iband])
    {
        strcpy (errmsg, "Reflectance band has not been opened");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iband < 0 || iband >= this->nband)
    {
        strcpy (errmsg, "Invalid reflectance band number for the input date");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iline < 0 || iline >= this->size.nlines)
    {
        strcpy (errmsg, "Invalid line number for reflectance band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Check the number of lines */
    if (nsamps == -99)
        nsamps = this->size.nsamps;
  
    /* Read the data, but first seek to the correct line */
    loc = (long) iline * nsamps * sizeof (uint16);
    if (fseek (this->fp_bin[iband], loc, SEEK_SET))
    {
        strcpy (errmsg, "Seeking to the current line in the input file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (this->fp_bin[iband], nlines, nsamps, sizeof (uint16),
        out_arr) != SUCCESS)
    {
        sprintf (errmsg, "Reading %d lines from reflectance band %d starting "
            "at line %d", nlines, iband, iline);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_input_th_lines

PURPOSE:  Reads the thermal data for the current thermal band and lines, and
populates the output buffer.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading data for this band
SUCCESS    Successful completion

NOTES:
  1. The Input_t data structure needs to be populated and memory allocated
     before calling this routine.  Use open_input to do that.
******************************************************************************/
int get_input_th_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current thermal band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    uint16 *out_arr  /* O: output array to populate */
)
{
    char FUNC_NAME[] = "get_input_th_line";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    long loc;                 /* current location in the input file */
  
    /* Check the parameters */
    if (this == NULL) 
    {
        strcpy (errmsg, "Input structure has not been opened/initialized");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (!this->open_th[iband])
    {
        strcpy (errmsg, "Thermal band has not been opened");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iband < 0 || iband >= this->nband_th)
    {
        strcpy (errmsg, "Invalid thermal band number for the input data");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iline < 0 || iline >= this->size_th.nlines)
    {
        strcpy (errmsg, "Invalid line number for thermal band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Read the data, but first seek to the correct line */
    loc = (long) iline * this->size_th.nsamps * sizeof (uint16);
    if (fseek (this->fp_bin_th[iband], loc, SEEK_SET))
    {
        strcpy (errmsg, "Seeking to the current line in the input file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (this->fp_bin_th[iband], nlines, this->size_th.nsamps,
        sizeof (uint16), out_arr) != SUCCESS)
    {
        sprintf (errmsg, "Reading %d lines from thermal band %d starting at "
            "line %d", nlines, iband, iline);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_input_pan_lines

PURPOSE:  Reads the pan data for the current pan band and lines, and populates
the output buffer.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading data for this band
SUCCESS    Successful completion

NOTES:
  1. The Input_t data structure needs to be populated and memory allocated
     before calling this routine.  Use open_input to do that.
******************************************************************************/
int get_input_pan_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current pan band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    uint16 *out_arr  /* O: output array to populate */
)
{
    char FUNC_NAME[] = "get_input_pan_line";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    long loc;                 /* current location in the input file */
  
    /* Check the parameters */
    if (this == NULL) 
    {
        strcpy (errmsg, "Input structure has not been opened/initialized");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (!this->open_pan[iband])
    {
        strcpy (errmsg, "Pan band has not been opened");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iband < 0 || iband >= this->nband_pan)
    {
        strcpy (errmsg, "Invalid pan band number for the input data");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iline < 0 || iline >= this->size_pan.nlines)
    {
        strcpy (errmsg, "Invalid line number for pan band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Read the data, but first seek to the correct line */
    loc = (long) iline * this->size_pan.nsamps * sizeof (uint16);
    if (fseek (this->fp_bin_pan[iband], loc, SEEK_SET))
    {
        strcpy (errmsg, "Seeking to the current line in the input file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (this->fp_bin_pan[iband], nlines, this->size_pan.nsamps,
        sizeof (uint16), out_arr) != SUCCESS)
    {
        sprintf (errmsg, "Reading %d lines from pan band %d starting at "
            "line %d", nlines, iband, iline);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_input_qa_lines

PURPOSE:  Reads the QA data for the current QA band and lines, and populates
the output buffer.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading data for this band
SUCCESS    Successful completion

NOTES:
  1. The Input_t data structure needs to be populated and memory allocated
     before calling this routine.  Use open_input to do that.
******************************************************************************/
int get_input_qa_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current QA band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    uint16 *out_arr  /* O: output array to populate */
)
{
    char FUNC_NAME[] = "get_input_qa_line";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    long loc;                 /* current location in the input file */
  
    /* Check the parameters */
    if (this == NULL) 
    {
        strcpy (errmsg, "Input structure has not been opened/initialized");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (!this->open_qa[iband])
    {
        strcpy (errmsg, "QA band has not been opened");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iband < 0 || iband >= this->nband_qa)
    {
        strcpy (errmsg, "Invalid QA band number for the input data");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iline < 0 || iline >= this->size_qa.nlines)
    {
        strcpy (errmsg, "Invalid line number for QA band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Read the data, but first seek to the correct line */
    loc = (long) iline * this->size_qa.nsamps * sizeof (uint16);
    if (fseek (this->fp_bin_qa[iband], loc, SEEK_SET))
    {
        strcpy (errmsg, "Seeking to the current line in the input file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (this->fp_bin_qa[iband], nlines, this->size_qa.nsamps,
        sizeof (uint16), out_arr) != SUCCESS)
    {
        sprintf (errmsg, "Reading %d lines from QA band %d starting at "
            "line %d", nlines, iband, iline);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_input_ppa_lines

PURPOSE:  Reads the per-pixel angle data for the current solar/view angle bands
and populates the output buffer.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading data for these bands
SUCCESS    Successful completion

NOTES:
  1. The Input_t data structure needs to be populated and memory allocated
     before calling this routine.  Use open_input to do that.
******************************************************************************/
int get_input_ppa_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    bool use_orig_aero,  /* I: use the original aerosol handling if specified,
                               o/w use the semi-empirical approach */
    int16 *sza_arr,  /* O: output solar zenith array to populate */
    int16 *saa_arr,  /* O: output solar azimuth array to populate */
    int16 *vza_arr,  /* O: output view zenith array to populate */
    int16 *vaa_arr   /* O: output view azimuth array to populate */
)
{
    char FUNC_NAME[] = "get_input_ppa_lines";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    long loc;                 /* current location in the input file */
  
    /* Check the parameters */
    if (this == NULL) 
    {
        strcpy (errmsg, "Input structure has not been opened/initialized");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (!this->open_ppa)
    {
        strcpy (errmsg, "Per-pixel angle bands have not been opened");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iline < 0 || iline >= this->size_ppa.nlines)
    {
        strcpy (errmsg, "Invalid line number for per-pixel angle bands");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Read the solar zenith data, but first seek to the correct line */
    loc = (long) iline * this->size_ppa.nsamps * sizeof (int16);
    if (fseek (this->fp_bin_sza, loc, SEEK_SET))
    {
        strcpy (errmsg, "Seeking to the current line in the sza input file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (this->fp_bin_sza, nlines, this->size_ppa.nsamps,
        sizeof (int16), sza_arr) != SUCCESS)
    {
        sprintf (errmsg, "Reading %d lines from solar zenith band starting "
            "at line %d", nlines, iline);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Only read these per-pixel angles if the original aerosol algorithm is
       being used */
    if (use_orig_aero)
    {
        /* Read the solar azimuth data, but first seek to the correct line */
        if (fseek (this->fp_bin_saa, loc, SEEK_SET))
        {
            strcpy (errmsg, "Seeking to current line in the saa input file");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (read_raw_binary (this->fp_bin_saa, nlines, this->size_ppa.nsamps,
            sizeof (int16), saa_arr) != SUCCESS)
        {
            sprintf (errmsg, "Reading %d lines from solar azimuth band "
                "starting at line %d", nlines, iline);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Read the view zenith data, but first seek to the correct line */
        if (fseek (this->fp_bin_vza, loc, SEEK_SET))
        {
            strcpy (errmsg, "Seeking to current line in the vza input file");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (read_raw_binary (this->fp_bin_vza, nlines, this->size_ppa.nsamps,
            sizeof (int16), vza_arr) != SUCCESS)
        {
            sprintf (errmsg, "Reading %d lines from view zenith band starting "
                "at line %d", nlines, iline);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Read the view azimuth data, but first seek to the correct line */
        if (fseek (this->fp_bin_vaa, loc, SEEK_SET))
        {
            strcpy (errmsg, "Seeking to current line in the vaa input file");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (read_raw_binary (this->fp_bin_vaa, nlines, this->size_ppa.nsamps,
            sizeof (int16), vaa_arr) != SUCCESS)
        {
            sprintf (errmsg, "Reading %d lines from view azimuth band starting "
                "at line %d", nlines, iline);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    return (SUCCESS);
}


#define DATE_STRING_LEN (50)
#define TIME_STRING_LEN (50)

/******************************************************************************
MODULE:  get_xml_input

PURPOSE:  Pulls input values from the XML structure.

RETURN VALUE:
Type = int
Value        Description
-------      -----------
ERROR        Error occurred opening or reading the file
SUCCESS      Successful completion

NOTES:
******************************************************************************/
int get_xml_input
(
    Espa_internal_meta_t *metadata,  /* I: XML metadata */
    bool process_sr,                 /* I: will SR data be processed? */
    bool use_orig_aero,              /* I: use the original aerosol handling if
                                           specified, o/w use the
                                           semi-empirical approach */
    Input_t *this                    /* O: data structure for the input file */
)
{
    char FUNC_NAME[] = "get_xml_input";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    int ib;                   /* looping variable for bands */
    char acq_date[DATE_STRING_LEN+1];    /* acquisition date */
    char prod_date[DATE_STRING_LEN+1];   /* production date */
    char acq_time[TIME_STRING_LEN+1];    /* acquisition time */
    char temp[STR_SIZE]; /* temporary string */
    int i;               /* looping variable */
    int refl_indx = NA;  /* band index in XML file for the reflectance band */
    int th_indx = NA;    /* band index in XML file for the thermal band */
    int pan_indx = NA;   /* band index in XML file for the pan band */
    int qa_indx = NA;    /* band index in XML file for the QA band */
    int sza_indx = NA;   /* band index in XML file for the solar zenith band */
    int saa_indx = NA;   /* band index in XML file for the solar azimuth band */
    int vza_indx = NA;   /* band index in XML file for the view zenith band */
    int vaa_indx = NA;   /* band index in XML file for the view azimuth band */
    Espa_global_meta_t *gmeta = &metadata->global; /* pointer to global meta */

    /* Initialize the input fields */
    this->meta.sat = SAT_NULL;
    this->meta.inst = INST_NULL;
    this->meta.acq_date.fill = true;
    this->meta.time_fill = true;
    this->meta.prod_date.fill = true;
    this->meta.sun_zen = ANGLE_FILL;
    this->meta.sun_az = ANGLE_FILL;
    this->meta.wrs_sys = (Wrs_t)WRS_NULL;
    this->meta.ipath = -1;
    this->meta.irow = -1;
    this->meta.fill = INPUT_FILL;
    this->size.nsamps = this->size.nlines = -1;
    this->meta.gain_set = false;

    this->nband = 0;
    for (ib = 0; ib < NBAND_REFL_MAX; ib++)
    {
        this->meta.gain[ib] = GAIN_BIAS_FILL;
        this->meta.bias[ib] = GAIN_BIAS_FILL;
        this->file_name[ib] = NULL;
        this->open[ib] = false;
        this->fp_bin[ib] = NULL;
    }

    /* use Landsat thermal band count since Sentinel doesn't have thermal */
    this->nband_th = 0;
    for (ib = 0; ib < NBANDL_THM_MAX; ib++)
    {
        this->meta.gain_th[ib] = GAIN_BIAS_FILL;
        this->meta.bias_th[ib] = GAIN_BIAS_FILL;
        this->file_name_th[ib] = NULL;
        this->open_th[ib] = false;
        this->fp_bin_th[ib] = NULL;
    }

    /* use Landsat pan band count as Sentinel doesn't have pan */
    this->nband_pan = 0;
    for (ib = 0; ib < NBANDL_PAN_MAX; ib++)
    {
        this->meta.gain_pan[ib] = GAIN_BIAS_FILL;
        this->meta.bias_pan[ib] = GAIN_BIAS_FILL;
        this->file_name_pan[ib] = NULL;
        this->open_pan[ib] = false;
        this->fp_bin_pan[ib] = NULL;
    }

    /* use Landsat QA band count as Sentinel doesn't have an input QA band */
    this->nband_qa = 0;
    for (ib = 0; ib < NBANDL_QA_MAX; ib++)
    {
        this->file_name_qa[ib] = NULL;
        this->open_qa[ib] = false;
        this->fp_bin_qa[ib] = NULL;
    }

    this->file_name_sza = NULL;
    this->file_name_saa = NULL;
    this->file_name_vza = NULL;
    this->file_name_vaa = NULL;
    this->open_ppa = NULL;
    this->fp_bin_sza = NULL;
    this->fp_bin_saa = NULL;
    this->fp_bin_vza = NULL;
    this->fp_bin_vaa = NULL;
    this->meta.gain_sza = GAIN_BIAS_FILL;
    this->meta.bias_sza = GAIN_BIAS_FILL;
    this->meta.gain_saa = GAIN_BIAS_FILL;
    this->meta.bias_saa = GAIN_BIAS_FILL;
    this->meta.gain_vza = GAIN_BIAS_FILL;
    this->meta.bias_vza = GAIN_BIAS_FILL;
    this->meta.gain_vaa = GAIN_BIAS_FILL;
    this->meta.bias_vaa = GAIN_BIAS_FILL;

    /* Pull the appropriate data from the XML file */
    acq_date[0] = acq_time[0] = '\0';
    prod_date[0] = '\0';
    if (!strcmp (gmeta->satellite, "LANDSAT_8"))
        this->meta.sat = SAT_LANDSAT_8;
    else if (!strcmp (gmeta->satellite, "LANDSAT_9"))
        this->meta.sat = SAT_LANDSAT_9;
    else if (!strncmp (gmeta->satellite, "Sentinel-2", 10))
        this->meta.sat = SAT_SENTINEL_2;
    else
    {
        sprintf (errmsg, "Unsupported satellite: %s", gmeta->satellite);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    printf ("Metadata satellite: %s\n", gmeta->satellite);
    printf ("         satellite (L8=%d, L9=%d, S2=%d): %d\n", SAT_LANDSAT_8,
        SAT_LANDSAT_9, SAT_SENTINEL_2, this->meta.sat);

    if (!strcmp (gmeta->instrument, "OLI_TIRS"))
        this->meta.inst = INST_OLI_TIRS;
    else if (!strcmp (gmeta->instrument, "OLI"))
        this->meta.inst = INST_OLI;
    else if (!strcmp (gmeta->instrument, "MSI"))
        this->meta.inst = INST_MSI;
    else
    {
        sprintf (errmsg, "Unsupported instrument: %s", gmeta->instrument);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    printf ("         instrument (OLI/TIRS=0, OLI=1, MSI=2): %d\n",
        this->meta.inst);

    strcpy (acq_date, gmeta->acquisition_date);
    strcpy (acq_time, gmeta->scene_center_time);
    this->meta.time_fill = false;
    //printf ("     JPR xxxxxxxx %s\n",  gmeta->scene_center_time);
    
    /* new 26-FEB-26 JPR */
    if (this->meta.sat == SAT_SENTINEL_2) {
        ad_hoc_XML_reader(gmeta->scene_center_time);
        //printf ("     JPR xxxxxxxx %s\n",  gmeta->scene_center_time);
       }

    /* Make sure the acquisition time is not too long (i.e. contains too
       many decimal points for the date/time routines).  The time should be
       hh:mm:ss.ssssssZ (see DATE_FORMAT_DATEA_TIME in date.h) which is 16
       characters long.  If the time is longer than that, just chop it off. */
    if (strlen (acq_time) > 16)
        sprintf (&acq_time[15], "Z");

    this->meta.sun_zen = gmeta->solar_zenith;
    if (this->meta.sun_zen < -180.0 || this->meta.sun_zen > 180.0)
    {
        sprintf (errmsg, "Solar zenith angle is out of range: %f",
            this->meta.sun_zen);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    this->meta.sun_az = gmeta->solar_azimuth;
    if (this->meta.sun_az < -360.0 || this->meta.sun_az > 360.0)
    {
        sprintf (errmsg, "Solar azimuth angle is out of range: %f",
            this->meta.sun_az);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (this->meta.sat == SAT_SENTINEL_2)
    {
        /* Sentinel currently is the only satellite that has the view angles in
           the metadata */
        this->meta.view_zen = gmeta->view_zenith;
        if (this->meta.view_zen < -180.0 || this->meta.view_zen > 180.0)
        {
            sprintf (errmsg, "View zenith angle is out of range: %f",
                this->meta.view_zen);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        this->meta.view_az = gmeta->view_azimuth;
        if (this->meta.view_az < -360.0 || this->meta.view_az > 360.0)
        {
            sprintf (errmsg, "View azimuth angle is out of range: %f",
                this->meta.view_az);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    if (this->meta.sat == SAT_LANDSAT_8 || this->meta.sat == SAT_LANDSAT_9)
    {
        switch (gmeta->wrs_system)
        {
            case 1: this->meta.wrs_sys = WRS_1; break;
            case 2: this->meta.wrs_sys = WRS_2; break;
            default:
                sprintf (errmsg, "Invalid WRS system: %d", gmeta->wrs_system);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
        }

        this->meta.ipath = gmeta->wrs_path;
        this->meta.irow = gmeta->wrs_row;
    }

    if (this->meta.inst == INST_OLI_TIRS)
    {
        this->nband = NBANDL_REFL_MAX;     /* number of reflectance bands */
        this->nband_th = NBANDL_THM_MAX;   /* number of thermal bands */
        this->nband_pan = NBANDL_PAN_MAX;  /* number of pan bands */
        this->nband_qa = NBANDL_QA_MAX;    /* number of QA bands */
    }
    else if (this->meta.inst == INST_OLI)
    {  /* No TIRS band */
        this->nband = NBANDL_REFL_MAX;     /* number of reflectance bands */
        this->nband_th = 0;                /* number of thermal bands */
        this->nband_pan = NBANDL_PAN_MAX;  /* number of pan bands */
        this->nband_qa = NBANDL_QA_MAX;    /* number of QA bands */
    }
    else if (this->meta.inst == INST_MSI)
    {
        this->nband = NBANDS_REFL_MAX;     /* number of reflectance bands */
        this->nband_th = 0;     /* number of thermal bands */
        this->nband_pan = 0;    /* number of pan bands */
        this->nband_qa = 0;     /* number of QA bands */
    }

    /* Band differences between Landsat and Sentinel */
    if (this->meta.sat == SAT_LANDSAT_8 || this->meta.sat == SAT_LANDSAT_9)
    {
        /* Use band 1 for the representative band */
        for (i = 0; i < metadata->nbands; i++)
        {
            if (!strcmp (metadata->band[i].name, "b1"))
            {
                /* this is the index we'll use for reflectance band info */
                refl_indx = i;
                break;
            }
        }

        if (refl_indx == NA)
        {
            sprintf (errmsg, "Band 1 not found in the input Landsat XML "
                "metadata file to be used as the representative band.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Find band 1, band 10, and band 8 in the input XML file to obtain
           band-related information for reflectance, thermal, and pan bands */
        for (i = 0; i < metadata->nbands; i++)
        {
            /* Check the band offset for not valid. In this case set it to 0. */
            if (metadata->band[i].add_offset == NOT_VALID)
                metadata->band[i].add_offset = 0.;

            /* Obtain band-related info */
            if (!strcmp (metadata->band[i].name, "b1"))
            {
                /* get the band1 info */
                this->meta.gain[DNL_BAND1] = metadata->band[i].refl_gain;
                this->meta.bias[DNL_BAND1] = metadata->band[i].refl_bias;
                this->file_name[DNL_BAND1] =
                    strdup (metadata->band[i].file_name);
    
                /* get the production date but only the date portion
                   (yyyy-mm-dd) */
                strncpy (prod_date, metadata->band[i].production_date, 10);
                prod_date[10] = '\0';
            }
            else if (!strcmp (metadata->band[i].name, "b2"))
            {
                /* get the band2 info */
                this->meta.gain[DNL_BAND2] = metadata->band[i].refl_gain;
                this->meta.bias[DNL_BAND2] = metadata->band[i].refl_bias;
                this->file_name[DNL_BAND2] =
                    strdup (metadata->band[i].file_name);
            }
            else if (!strcmp (metadata->band[i].name, "b3"))
            {
                /* get the band3 info */
                this->meta.gain[DNL_BAND3] = metadata->band[i].refl_gain;
                this->meta.bias[DNL_BAND3] = metadata->band[i].refl_bias;
                this->file_name[DNL_BAND3] =
                    strdup (metadata->band[i].file_name);
            }
            else if (!strcmp (metadata->band[i].name, "b4"))
            {
                /* get the band4 info */
                this->meta.gain[DNL_BAND4] = metadata->band[i].refl_gain;
                this->meta.bias[DNL_BAND4] = metadata->band[i].refl_bias;
                this->file_name[DNL_BAND4] =
                    strdup (metadata->band[i].file_name);
            }
            else if (!strcmp (metadata->band[i].name, "b5"))
            {
                /* get the band5 info */
                this->meta.gain[DNL_BAND5] = metadata->band[i].refl_gain;
                this->meta.bias[DNL_BAND5] = metadata->band[i].refl_bias;
                this->file_name[DNL_BAND5] =
                    strdup (metadata->band[i].file_name);
            }
            else if (!strcmp (metadata->band[i].name, "b6"))
            {
                /* get the band6 info */
                this->meta.gain[DNL_BAND6] = metadata->band[i].refl_gain;
                this->meta.bias[DNL_BAND6] = metadata->band[i].refl_bias;
                this->file_name[DNL_BAND6] =
                    strdup (metadata->band[i].file_name);
            }
            else if (!strcmp (metadata->band[i].name, "b7"))
            {
                /* get the band7 info */
                this->meta.gain[DNL_BAND7] = metadata->band[i].refl_gain;
                this->meta.bias[DNL_BAND7] = metadata->band[i].refl_bias;
                this->file_name[DNL_BAND7] =
                    strdup (metadata->band[i].file_name);
            }
            else if (!strcmp (metadata->band[i].name, "b9"))
            {
                /* get the band9 info - store in the reflectance array
                   skipping band 8*/
                this->meta.gain[DNL_BAND9-1] = metadata->band[i].refl_gain;
                this->meta.bias[DNL_BAND9-1] = metadata->band[i].refl_bias;
                this->file_name[DNL_BAND9-1] =
                    strdup (metadata->band[i].file_name);
            }
    
            else if (!strcmp (metadata->band[i].name, "b8"))
            {
                /* this is the index we'll use for pan band info */
                pan_indx = i;
    
                /* get the band8 info */
                this->meta.gain_pan[0] = metadata->band[i].refl_gain;
                this->meta.bias_pan[0] = metadata->band[i].refl_bias;
                this->file_name_pan[0] = strdup (metadata->band[i].file_name);
            }
    
            /* NOTE: band10 and band11 won't exist in the input XML file
               for OLI-only products */
            else if (!strcmp (metadata->band[i].name, "b10"))
            {
                /* this is the index we'll use for thermal band info */
                th_indx = i;
    
                /* get the band10 info */
                this->meta.gain_th[0] = metadata->band[i].rad_gain;
                this->meta.bias_th[0] = metadata->band[i].rad_bias;
                this->meta.k1_const[0] = metadata->band[i].k1_const;
                this->meta.k2_const[0] = metadata->band[i].k2_const;
                this->file_name_th[0] = strdup (metadata->band[i].file_name);
            }
            else if (!strcmp (metadata->band[i].name, "b11"))
            {
                /* get the band11 info */
                this->meta.gain_th[1] = metadata->band[i].rad_gain;
                this->meta.bias_th[1] = metadata->band[i].rad_bias;
                this->meta.k1_const[1] = metadata->band[i].k1_const;
                this->meta.k2_const[1] = metadata->band[i].k2_const;
                this->file_name_th[1] = strdup (metadata->band[i].file_name);
            }
    
            else if (!strcmp (metadata->band[i].name, "qa_pixel"))
            {
                /* this is the index we'll use for pixel qa band info */
                qa_indx = i;
    
                /* get the QA band info */
                this->file_name_qa[0] = strdup (metadata->band[i].file_name);
            }
    
            else if (!strcmp (metadata->band[i].name, "sza"))
            {
                /* this is the index we'll use for sza band info */
                sza_indx = i;
    
                /* get the solar zenith band info */
                this->file_name_sza = strdup (metadata->band[i].file_name);
                this->meta.gain_sza = metadata->band[i].scale_factor;
                this->meta.bias_sza = metadata->band[i].add_offset;
                if (this->meta.bias_sza == NOT_VALID)
                    this->meta.bias_sza = 0;
            }

            else if (!strcmp (metadata->band[i].name, "saa") && use_orig_aero)
            {
                /* this is the index we'll use for saa band info */
                saa_indx = i;

                /* get the solar azimuth band info */
                this->file_name_saa = strdup (metadata->band[i].file_name);
                this->meta.gain_saa = metadata->band[i].scale_factor;
                this->meta.bias_saa = metadata->band[i].add_offset;
                if (this->meta.bias_saa == NOT_VALID)
                    this->meta.bias_saa = 0;
            }

            else if (!strcmp (metadata->band[i].name, "vza") && use_orig_aero)
            {
                /* this is the index we'll use for vza band info */
                vza_indx = i;

                /* get the view zenith band info */
                this->file_name_vza = strdup (metadata->band[i].file_name);
                this->meta.gain_vza = metadata->band[i].scale_factor;
                this->meta.bias_vza = metadata->band[i].add_offset;
                if (this->meta.bias_vza == NOT_VALID)
                    this->meta.bias_vza = 0;
            }

            else if (!strcmp (metadata->band[i].name, "vaa") && use_orig_aero)
            {
                /* this is the index we'll use for vaa band info */
                vaa_indx = i;

                /* get the view azimuth band info */
                this->file_name_vaa = strdup (metadata->band[i].file_name);
                this->meta.gain_vaa = metadata->band[i].scale_factor;
                this->meta.bias_vaa = metadata->band[i].add_offset;
                if (this->meta.bias_vaa == NOT_VALID)
                    this->meta.bias_vaa = 0;
            }
        }  /* for i */

        /* Make sure the expected files were found */
        if (this->meta.inst == INST_OLI_TIRS && th_indx == NA)
        {
            sprintf (errmsg, "Band 10 (b10) was not found in the XML file");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    
        if (pan_indx == NA)
        {
            sprintf (errmsg, "Band 8 (b8) was not found in the XML file");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    
        if (qa_indx == NA)
        {
            sprintf (errmsg, "QA band (bqa) was not found in the XML file");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    
        if (sza_indx == NA)
        {
            sprintf (errmsg, "Solar zenith band not found in the XML file");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (use_orig_aero)
        {
            if (saa_indx == NA)
            {
                sprintf (errmsg, "Solar azimuth band not found in XML file");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            if (vza_indx == NA)
            {
                sprintf (errmsg, "View zenith band was found in the XML file");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            if (vaa_indx == NA)
            {
                sprintf (errmsg, "View azimuth band not found in the XML file");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }
    }
    else if (this->meta.sat == SAT_SENTINEL_2)
    {
        /* Store the band-related information. Use band 2 for the
           representative band. Skip bands 9 and 10 if they won't be
           processed. */
        for (i = 0; i < metadata->nbands; i++)
        {
            /* Check the band offset for not valid. In this case set it to 0. */
            if (metadata->band[i].add_offset == NOT_VALID)
                metadata->band[i].add_offset = 0.;

            /* Obtain band-related info */
            if (!strcmp (metadata->band[i].name, "B01"))
                this->file_name[DNS_BAND1] =
                    strdup (metadata->band[i].file_name);
            else if (!strcmp (metadata->band[i].name, "B02"))
            {
                /* This is the index we'll use for reflectance band info,
                   since it is a 30m band. */
                refl_indx = i;
                this->file_name[DNS_BAND2] =
                    strdup (metadata->band[i].file_name);

                /* Get the production date of the level-1 data */
                strncpy (prod_date, metadata->band[i].production_date, 10);
                prod_date[10] = '\0';
            }
            else if (!strcmp (metadata->band[i].name, "B03"))
                this->file_name[DNS_BAND3] =
                    strdup (metadata->band[i].file_name);
            else if (!strcmp (metadata->band[i].name, "B04"))
                this->file_name[DNS_BAND4] =
                    strdup (metadata->band[i].file_name);
            else if (!strcmp (metadata->band[i].name, "B05"))
                this->file_name[DNS_BAND5] =
                    strdup (metadata->band[i].file_name);
            else if (!strcmp (metadata->band[i].name, "B06"))
                this->file_name[DNS_BAND6] =
                    strdup (metadata->band[i].file_name);
            else if (!strcmp (metadata->band[i].name, "B07"))
                this->file_name[DNS_BAND7] =
                    strdup (metadata->band[i].file_name);
            else if (!strcmp (metadata->band[i].name, "B08"))
                this->file_name[DNS_BAND8] =
                    strdup (metadata->band[i].file_name);
            else if (!strcmp (metadata->band[i].name, "B8A"))
                this->file_name[DNS_BAND8A] =
                    strdup (metadata->band[i].file_name);

#ifdef PROC_ALL_BANDS
/* Process all bands if turned on */
            else if (!strcmp (metadata->band[i].name, "B09"))
                this->file_name[DNS_BAND9] =
                    strdup (metadata->band[i].file_name);
            else if (!strcmp (metadata->band[i].name, "B10"))
                this->file_name[DNS_BAND10] =
                    strdup (metadata->band[i].file_name);
#endif

            else if (!strcmp (metadata->band[i].name, "B11"))
                this->file_name[DNS_BAND11] =
                    strdup (metadata->band[i].file_name);
            else if (!strcmp (metadata->band[i].name, "B12"))
                this->file_name[DNS_BAND12] =
                    strdup (metadata->band[i].file_name);
        }

        if (refl_indx == NA)
        {
            sprintf (errmsg, "Band 2 not found in the input Sentinel XML "
                "metadata file to be used as the representative band.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Get the size of the reflectance, thermal, pan, etc. bands by using
       the representative band in the XML file */
    /* reflectance */
    this->size.nsamps = metadata->band[refl_indx].nsamps;
    this->size.nlines = metadata->band[refl_indx].nlines;
    this->size.pixsize[0] = metadata->band[refl_indx].pixel_size[0];
    this->size.pixsize[1] = metadata->band[refl_indx].pixel_size[1];
    this->scale_factor = metadata->band[refl_indx].scale_factor;

    /* thermal */
    if (this->meta.inst == INST_OLI_TIRS)
    {  /* skip for OLI */
        this->size_th.nsamps = metadata->band[th_indx].nsamps;
        this->size_th.nlines = metadata->band[th_indx].nlines;
        this->size_th.pixsize[0] = metadata->band[th_indx].pixel_size[0];
        this->size_th.pixsize[1] = metadata->band[th_indx].pixel_size[1];
        this->scale_factor_th = metadata->band[th_indx].scale_factor;
    }
    else
    {  /* set to 0s */
        this->size_th.nsamps = 0;
        this->size_th.nlines = 0;
        this->size_th.pixsize[0] = 0;
        this->size_th.pixsize[1] = 0;
        this->scale_factor_th = 0;
    }

    /* Landsat-specific input params -- including pan and QA bands */
    if (this->meta.sat == SAT_LANDSAT_8 || this->meta.sat == SAT_LANDSAT_9)
    {
        this->size_pan.nsamps = metadata->band[pan_indx].nsamps;
        this->size_pan.nlines = metadata->band[pan_indx].nlines;
        this->size_pan.pixsize[0] = metadata->band[pan_indx].pixel_size[0];
        this->size_pan.pixsize[1] = metadata->band[pan_indx].pixel_size[1];
        this->scale_factor_pan = metadata->band[pan_indx].scale_factor;
    
        this->size_qa.nsamps = metadata->band[qa_indx].nsamps;
        this->size_qa.nlines = metadata->band[qa_indx].nlines;
        this->size_qa.pixsize[0] = metadata->band[qa_indx].pixel_size[0];
        this->size_qa.pixsize[1] = metadata->band[qa_indx].pixel_size[1];
    
        /* Assume the per-pixel angle bands all have the same size and
           resolution */
        this->size_ppa.nsamps = metadata->band[sza_indx].nsamps;
        this->size_ppa.nlines = metadata->band[sza_indx].nlines;
        this->size_ppa.pixsize[0] = metadata->band[sza_indx].pixel_size[0];
        this->size_ppa.pixsize[1] = metadata->band[sza_indx].pixel_size[1];
    
        /* Check WRS path/rows */
        if (this->meta.wrs_sys == WRS_1)
        {
            if (this->meta.ipath > WRS1_NPATH)
            {
                sprintf (errmsg, "WRS path number out of range: %d",
                    this->meta.ipath);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            else if (this->meta.irow > WRS1_NROW)
            {
                sprintf (errmsg, "WRS row number out of range: %d",
                    this->meta.irow);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }
        else if (this->meta.wrs_sys == WRS_2)
        {
            if (this->meta.ipath > WRS2_NPATH)
            {
                sprintf (errmsg, "WRS path number out of range: %d",
                    this->meta.ipath);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            else if (this->meta.irow > WRS2_NROW)
            {
                sprintf (errmsg, "WRS row number out of range: %d",
                    this->meta.irow);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }
    }
    else if (this->meta.sat == SAT_SENTINEL_2)
    { /* setup output QA */
        this->size_qa.nsamps = this->size.nsamps;
        this->size_qa.nlines = this->size.nlines;
        this->size_qa.pixsize[0] = this->size.pixsize[0];
        this->size_qa.pixsize[1] = this->size.pixsize[1];
    }

    /* Check satellite/instrument combinations */
    if (this->meta.inst == INST_OLI_TIRS ||
        this->meta.inst == INST_OLI)
    {
        if (this->meta.sat != SAT_LANDSAT_8 && this->meta.sat != SAT_LANDSAT_9)
        {
            sprintf (errmsg, "Invalid instrument/satellite combination");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    if (this->meta.inst == INST_MSI)
    {
        if (this->meta.sat != SAT_SENTINEL_2)
        {
            sprintf (errmsg, "Invalid instrument/satellite combination");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Convert the acquisition date/time values from string to date struct */
    if (this->meta.inst != INST_MSI)
    {
        sprintf (temp, "%sT%s", acq_date, acq_time);
        if (!date_init (&this->meta.acq_date, temp, DATE_FORMAT_DATEA_TIME))
        {
            sprintf (errmsg, "Converting the acquisition date and time");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }
    else
    {
        /*sprintf (temp, "%sT00:00:00.00000Z", acq_date);   26-FEB-26, JPR */
        sprintf (temp, gmeta->scene_center_time);  /* read by ad_hoc_XML_reader() */
        if (!date_init (&this->meta.acq_date, temp, DATE_FORMAT_DATEA_TIME))
        {
            sprintf (errmsg, "Converting the acquisition date and time");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Convert the production date value from string to date struct */
    if (!date_init (&this->meta.prod_date, prod_date, DATE_FORMAT_DATEA))
    {
        sprintf (errmsg, "Converting the production date and time");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  convert_to_10m

PURPOSE:  Converts the input band to a 10m band.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred converting the band
SUCCESS    Successful completion

NOTES:
  1. The assumption is that the input resolution is lower than the output
     resolution and therefore we are converting to a higher resolution.
  2. The resolution of the input band must be an even multiple of the output
     10m band.
  3. It is assumed the output array is of the correct size for this resolution.
******************************************************************************/
int convert_to_10m
(
    int in_nlines,   /* I: number of lines in the input product */
    int in_nsamps,   /* I: number of samples in the input product */
    int out_nlines,  /* I: number of lines in the output 10m product */
    int out_nsamps,  /* I: number of samples in the output 10m product */
    uint16 *in_arr,  /* I: input array */
    uint16 *out_arr  /* O: output 10m array */
)
{
    char FUNC_NAME[] = "convert_to_10m";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    int oline, osamp;  /* line/sample indices for the output array */
    int iline, isamp;  /* line/sample indices for the input array */
    int lmult;         /* line multiplier from input to output resolution */
    int smult;         /* sample multiplier from input to output resolution */
    int in_indx;       /* pixel index of input array */
    int out_indx;      /* pixel index of output array */

    /* Make sure the output number of lines and samples is an exact multiple
       of the input number of lines and samples */
    if (out_nlines % in_nlines != 0)
    {
        sprintf (errmsg, "Output number of lines (%d) is not an exact multiple "
            "of the input number of lines (%d)\n", out_nlines, in_nlines);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (out_nsamps % in_nsamps != 0)
    {
        sprintf (errmsg, "Output number of samples (%d) is not an exact "
            "multiple of the input number of samples (%d)\n", out_nsamps,
            in_nsamps);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Determine the multiplier from the input number of lines/samps to the
       output number of lines/samps */
    lmult = out_nlines / in_nlines;
    smult = out_nsamps / in_nsamps;

    /* Loop through the lines in the output array */
    out_indx = 0;
    for (oline = 0; oline < out_nlines; oline++)
    {
        /* Determine the representative line in the input array */
        iline = oline / lmult;

        /* Loop through the samples in the output array */
        for (osamp = 0; osamp < out_nsamps; osamp++, out_indx++)
        {
            /* Determine the representative sample in the input array */
            isamp = osamp / smult;
            in_indx = iline * in_nsamps + isamp;

            /* Copy the input pixel to the output array */
            out_arr[out_indx] = in_arr[in_indx];
        }
    }

    return (SUCCESS);
}
