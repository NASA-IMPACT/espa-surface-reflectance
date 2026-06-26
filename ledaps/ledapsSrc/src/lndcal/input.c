/*
!C****************************************************************************

!File: input.c
  
!Description: Functions reading data from the input data file.

!Team Unique Header:
  This software was developed by the MODIS Land Science Team Support 
  Group for the Laboratory for Terrestrial Physics (Code 922) at the 
  National Aeronautics and Space Administration, Goddard Space Flight 
  Center, under NASA Task 92-012-00.

 ! References and Credits:
  ! MODIS Science Team Member:
      Christopher O. Justice
      MODIS Land Science Team           University of Maryland
      justice@hermes.geog.umd.edu       Dept. of Geography
      phone: 301-405-1600               1113 LeFrak Hall
                                        College Park, MD, 20742

  ! Developers:
      Robert E. Wolfe (Code 922)
      MODIS Land Team Support Group     Raytheon ITSS
      robert.e.wolfe.1@gsfc.nasa.gov    4400 Forbes Blvd.
      phone: 301-614-5508               Lanham, MD 20770  
  
 ! Design Notes:
   1. The following public functions handle the input data:

	OpenInput - Setup 'input' data structure and open file for access.
	CloseInput - Close the input file.
	FreeInput - Free the 'input' data structure memory.

   2. 'OpenInput' must be called before any of the other routines.  
   3. 'FreeInput' should be used to free the 'input' data structure.

!END****************************************************************************
*/

#include <stdlib.h>
#include "input.h"
#include "error.h"
#include "mystring.h"
#include "const.h"
#include "date.h"
#include "local_defines.h"
#define INPUT_FILL (0)

/* Functions */
Input_t *OpenInput(Espa_internal_meta_t *metadata)
/* 
!C******************************************************************************

!Description: 'OpenInput' sets up the 'input' data structure, opens the
 input raw binary files for read access.
 
!Input Parameters:
 metadata     'Espa_internal_meta_t' data structure with XML info

!Output Parameters:
 (returns)      'input' data structure or NULL when an error occurs

!Team Unique Header:

!END****************************************************************************
*/
{
  Input_t *this = NULL;
  char *error_string = NULL;
  int ib;

  /* Create the Input data structure */
  this = (Input_t *)calloc(1, sizeof(Input_t));
  if (this == NULL) 
    RETURN_ERROR("allocating Input data structure", "OpenInput", NULL);

  /* Initialize and get input from header file */
  if (!GetXMLInput (this, metadata)) {
    free(this);
    this = NULL;
    RETURN_ERROR("getting input from header file", "OpenInput", NULL);
  }

  /* Open files for access */
  for (ib = 0; ib < this->nband; ib++) {
    this->fp_bin[ib] = fopen(this->file_name[ib], "r");
    if (this->fp_bin[ib] == NULL) {
      error_string = "opening binary file";
      break;
    }
    this->open[ib] = true;
  }
  if (!error_string) {
      /* TM band 6 */
      if (this->nband_th == 1) {
        this->fp_bin_th[0] = fopen(this->file_name_th[0], "r");
        if (this->fp_bin_th[0] == NULL)
          error_string = "opening thermal binary file";
        else
          this->open_th = true;
      }
      /* ETM+ band 6L & 6H */
      else if (this->nband_th == 2)
      {
        this->open_th = false;
        /* Open B6L */
        this->fp_bin_th[BAND_6L] = fopen(this->file_name_th[BAND_6L], "r");
        if (this->fp_bin_th[BAND_6L] == NULL) {
          error_string = "opening ETM+ band 6L binary file";
        }
        else
        {
          /* Open B6H */
          this->fp_bin_th[BAND_6H] = fopen(this->file_name_th[BAND_6H], "r");
          if (this->fp_bin_th[BAND_6H] == NULL) {
              error_string = "opening ETM+ band 6H binary file";
          }
          else {
              this->open_th = true;
          }
        }
      }
  }
  if (!error_string) {
      this->fp_bin_sun_zen = fopen(this->file_name_sun_zen, "r");
      if (this->fp_bin_sun_zen == NULL) 
        error_string = "opening solar zenith representative band binary file";
      else
        this->open_sun_zen = true;
  }

  if (!error_string) {
      this->fp_bin_band_qa = fopen(this->file_name_band_qa, "r");
      if (this->fp_bin_band_qa == NULL)
        error_string = "opening QA band binary file";
      else
        this->open_band_qa = true;
  }

  if (error_string != NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      free(this->file_name[ib]);
      this->file_name[ib] = NULL;

      if (this->open[ib]) {
        fclose(this->fp_bin[ib]);
        this->fp_bin[ib] = NULL;
        this->open[ib] = false;
      }
    }
    free(this->file_name_th[BAND_6L]);
    this->file_name_th[BAND_6L] = NULL;
    free(this->file_name_th[BAND_6H]);
    this->file_name_th[BAND_6H] = NULL;

    if (this->fp_bin_th[BAND_6L]) {
      fclose(this->fp_bin_th[BAND_6L]);
      this->fp_bin_th[BAND_6L] = NULL;
    }
    if (this->fp_bin_th[BAND_6H]) {
      fclose(this->fp_bin_th[BAND_6H]);
      this->fp_bin_th[BAND_6H] = NULL;
    }
    this->open_th = false;

    free(this->file_name_sun_zen);
    this->file_name_sun_zen = NULL;
    if (this->fp_bin_sun_zen) {
        fclose(this->fp_bin_sun_zen);
        this->fp_bin_sun_zen = NULL;
    }

    free(this->file_name_band_qa);
    this->file_name_band_qa = NULL;
    if (this->fp_bin_band_qa) {
        fclose(this->fp_bin_band_qa);
        this->fp_bin_band_qa = NULL;
    }

    free(this);
    this = NULL;
    RETURN_ERROR(error_string, "OpenInput", NULL);
  }

  return this;
}

bool GetInputLine(Input_t *this, int iband, int iline, unsigned char *line) 
{
  long loc;

  if (this == NULL) 
    RETURN_ERROR("invalid input structure", "GetInputLine", false);
  if (iband < 0  ||  iband >= this->nband) 
    RETURN_ERROR("band index out of range", "GetInputLine", false);
  if (iline < 0  ||  iline >= this->size.l) 
    RETURN_ERROR("line index out of range", "GetInputLine", false);
  if (!this->open[iband])
    RETURN_ERROR("band not open", "GetInputLine", false);

  loc = (long) (iline * this->size.s * sizeof(uint8));
  if (fseek(this->fp_bin[iband], loc, SEEK_SET))
    RETURN_ERROR("error seeking line (binary)", "GetInputLine", false);
  if (fread(line, sizeof(uint8), (size_t)this->size.s, 
            this->fp_bin[iband]) != (size_t)this->size.s)
    RETURN_ERROR("error reading line (binary)", "GetInputLine", false);

  return true;
}

bool GetInputLineTh(Input_t *this, int iline, unsigned char *line,
    unsigned char *line_b6h)
{
  long loc;

  if (this == NULL) 
    RETURN_ERROR("invalid input structure", "GetInputLineTh", false);
  if ( this->nband_th < 1 ) 
    RETURN_ERROR("no thermal input band", "GetInputLineTh", false);
  if (iline < 0  ||  iline >= this->size_th.l) 
    RETURN_ERROR("line index out of range", "GetInputLineTh", false);
  if (!this->open_th)
    RETURN_ERROR("band not open", "GetInputLineTh", false);

  loc = (long) (iline * this->size_th.s * sizeof(uint8));
  if (this->nband_th >= 1) {
    /* Read the first band which is either TM Band 6 or ETM+ Band 6L */
    if (fseek(this->fp_bin_th[BAND_6L], loc, SEEK_SET))
      RETURN_ERROR("error seeking b6/b6L line (binary)", "GetInputLineTh",
        false);
    if (fread(line, sizeof(uint8), (size_t)this->size_th.s, 
      this->fp_bin_th[BAND_6L]) != (size_t)this->size_th.s)
      RETURN_ERROR("error reading b6/b6L line (binary)", "GetInputLineTh",
        false);
  }

  if (this->nband_th > 1) {
    /* Read the second band which is ETM+ Band 6H */
    if (fseek(this->fp_bin_th[BAND_6H], loc, SEEK_SET))
      RETURN_ERROR("error seeking b6H line (binary)", "GetInputLineTh",
        false);
    if (fread(line_b6h, sizeof(uint8), (size_t)this->size_th.s,
      this->fp_bin_th[BAND_6H]) != (size_t)this->size_th.s)
      RETURN_ERROR("error reading b6H line (binary)", "GetInputLineTh",
        false);
  }

  return true;
}

bool GetInputLineQA(Input_t *this, int iline, uint16_t *line)
{
  long loc;

  if (this == NULL)
    RETURN_ERROR("invalid input structure", "GetInputLineQA", false);
  if (iline < 0  ||  iline >= this->size.l)
    RETURN_ERROR("line index out of range", "GetInputLineQA", false);
  if (!this->open_band_qa)
    RETURN_ERROR("band not open", "GetInputLineQA", false);

  loc = (long) (iline * this->size.s * sizeof(uint16_t));
  if (fseek(this->fp_bin_band_qa, loc, SEEK_SET))
    RETURN_ERROR("error seeking line (binary)", "GetInputLineQA", false);
  if (fread(line, sizeof(uint16_t), (size_t)this->size.s,
            this->fp_bin_band_qa) != (size_t)this->size.s)
    RETURN_ERROR("error reading line (binary)", "GetInputLineQA", false);

  return true;
}

bool GetInputLineSunZen(Input_t *this, int iline, int16 *line)
{
  long loc;

  if (this == NULL)
    RETURN_ERROR("invalid input structure", "GetInputLineSunZen", false);
  if (iline < 0  ||  iline >= this->size.l)
    RETURN_ERROR("line index out of range", "GetInputLineSunZen", false);
  if (!this->open_sun_zen)
    RETURN_ERROR("band not open", "GetInputLineSunZen", false);

  loc = (long) (iline * this->size.s * sizeof(int16));
  if (fseek(this->fp_bin_sun_zen, loc, SEEK_SET))
    RETURN_ERROR("error seeking line (binary)", "GetInputLineSunZen", false);
  if (fread(line, sizeof(int16), (size_t)this->size.s,
            this->fp_bin_sun_zen) != (size_t)this->size.s)
    RETURN_ERROR("error reading line (binary)", "GetInputLineSunZen", false);

  return true;
}


bool CloseInput(Input_t *this)
/* 
!C******************************************************************************

!Description: 'CloseInput' ends SDS access and closes the input file.
 
!Input Parameters:
 this           'input' data structure

!Output Parameters:
 this           'input' data structure; the following fields are modified:
                   open
 (returns)      status:
                  'true' = okay
                  'false' = error return

!Team Unique Header:

!END****************************************************************************
*/
{
  int ib;
  bool none_open;

  if (this == NULL) 
    RETURN_ERROR("invalid input structure", "CloseInput", false);

  none_open = true;
  for (ib = 0; ib < this->nband; ib++) {
    if (this->open[ib]) {
      none_open = false;
      fclose(this->fp_bin[ib]);
      this->fp_bin[ib] = NULL;
      this->open[ib] = false;
    }
  }

  /*** now close the thermal file ***/
  if (this->open_th) {
    if (this->fp_bin_th[BAND_6L]) {
      fclose(this->fp_bin_th[BAND_6L]);
      this->fp_bin_th[BAND_6L] = NULL;
    }
    if (this->fp_bin_th[BAND_6H]) {
      fclose(this->fp_bin_th[BAND_6H]);
      this->fp_bin_th[BAND_6H] = NULL;
    }
    this->open_th = false;
  }

  if (this->open_sun_zen) {
    fclose(this->fp_bin_sun_zen);
    this->fp_bin_sun_zen = NULL;
    this->open_sun_zen = false;
  }
  if (this->open_band_qa) {
    fclose(this->fp_bin_band_qa);
    this->fp_bin_band_qa = NULL;
    this->open_band_qa = false;
  }

  if (none_open)
    RETURN_ERROR("no files open", "CloseInput", false);

  return true;
}


bool FreeInput(Input_t *this)
/* 
!C******************************************************************************

!Description: 'FreeInput' frees the 'input' data structure memory.
 
!Input Parameters:
 this           'input' data structure

!Output Parameters:
 (returns)      status:
                  'true' = okay (always returned)

!Team Unique Header:

!END****************************************************************************
*/
{
  int ib;

  if (this != NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      free(this->file_name[ib]);
      this->file_name[ib] = NULL;
    }

    free(this->file_name_th[0]);
    this->file_name_th[0] = NULL;
    free(this->file_name_th[1]);
    this->file_name_th[1] = NULL;

    free(this->file_name_sun_zen);
    this->file_name_sun_zen = NULL;

    free(this);
    this = NULL;
  }

  return true;
}


#define DATE_STRING_LEN (50)
#define TIME_STRING_LEN (50)

bool GetXMLInput(Input_t *this, Espa_internal_meta_t *metadata)
/* 
!C******************************************************************************

!Description: 'GetXMLInput' pulls input values from the XML structure.
 
!Input Parameters:
 this         'Input_t' data structure to be populated
 metadata     'Espa_internal_meta_t' data structure with XML info

!Output Parameters:
 (returns)      status:
                  'true' = okay (always returned)
                  'false' = error determining if the gains/biases were provided

!Team Unique Header:

! Design Notes:
  1. This replaces the previous GetHeaderInput so the input values are pulled
     from the XML file instead of the header file (*.metadata.txt).
  2. Given that LPGS writes the gain values, the gain settings (HIGH, LOW) are
     no longer needed.

!END****************************************************************************
*/
{
    char *error_string = NULL;
    int ib;
    char acq_date[DATE_STRING_LEN + 1];
    char prod_date[DATE_STRING_LEN + 1];
    char acq_time[TIME_STRING_LEN + 1];
    char temp[MAX_STR_LEN + 1];
    int i;               /* looping variable */
    int refl_indx=-99;   /* band index in XML file for the reflectance band */
    int xml_thermal_index=5;  /* band index in XML file for the thermal band */
    Espa_global_meta_t *gmeta = &metadata->global; /* pointer to global meta */

    /* Initialize the input fields.  Set file type to binary, since that is
       the ESPA internal format for the input Level-1 products. */
    this->meta.sat = SAT_NULL;
    this->meta.inst = INST_NULL;
    this->meta.acq_date.fill = true;
    this->meta.time_fill = true;
    this->meta.prod_date.fill = true;
    this->meta.sun_zen = ANGLE_FILL;
    this->meta.sun_az = ANGLE_FILL;
    this->meta.szen_scale = 1.0;
    this->meta.szen_offset = 0.0;
    this->meta.wrs_sys = (Wrs_t)WRS_FILL;
    this->meta.ipath = -1;
    this->meta.irow = -1;
    this->meta.fill = INPUT_FILL;
    this->nband = 0;
    this->size.s = this->size.l = -1;
    for (ib = 0; ib < NBAND_REFL_MAX; ib++)
    {
        this->meta.iband[ib] = -1;
        this->meta.rad_gain[ib] = GAIN_BIAS_FILL;
        this->meta.rad_bias[ib] = GAIN_BIAS_FILL;
        this->meta.refl_gain[ib] = GAIN_BIAS_FILL;
        this->meta.refl_bias[ib] = GAIN_BIAS_FILL;
        this->file_name[ib] = NULL;
        this->open[ib] = false;
        this->fp_bin[ib] = NULL;
    }
    this->nband_th = 0;
    this->open_th = false;
    this->meta.rad_gain_th[0] = GAIN_BIAS_FILL;
    this->meta.rad_gain_th[1] = GAIN_BIAS_FILL;
    this->meta.rad_bias_th[0] = GAIN_BIAS_FILL;
    this->meta.rad_bias_th[1] = GAIN_BIAS_FILL;
    this->meta.k1_const[0] = GAIN_BIAS_FILL;
    this->meta.k1_const[1] = GAIN_BIAS_FILL;
    this->meta.k2_const[0] = GAIN_BIAS_FILL;
    this->meta.k2_const[1] = GAIN_BIAS_FILL;
    this->file_name_th[0] = NULL;
    this->file_name_th[1] = NULL;
    this->fp_bin_th[0] = NULL;
    this->fp_bin_th[1] = NULL;

    this->open_sun_zen = false;
    this->file_name_sun_zen = NULL;
    this->fp_bin_sun_zen = NULL;

    this->open_band_qa = false;
    this->file_name_band_qa = NULL;
    this->fp_bin_band_qa = NULL;

    /* Pull the appropriate data from the XML file */
    acq_date[0] = acq_time[0] = '\0';
    prod_date[0] = '\0';
    if (!strcmp (gmeta->satellite, "LANDSAT_1"))
        this->meta.sat = SAT_LANDSAT_1;
    else if (!strcmp (gmeta->satellite, "LANDSAT_2"))
        this->meta.sat = SAT_LANDSAT_2;
    else if (!strcmp (gmeta->satellite, "LANDSAT_3"))
        this->meta.sat = SAT_LANDSAT_3;
    else if (!strcmp (gmeta->satellite, "LANDSAT_4"))
        this->meta.sat = SAT_LANDSAT_4;
    else if (!strcmp (gmeta->satellite, "LANDSAT_5"))
        this->meta.sat = SAT_LANDSAT_5;
    else if (!strcmp (gmeta->satellite, "LANDSAT_7"))
        this->meta.sat = SAT_LANDSAT_7;
    else
    {
        sprintf (temp, "invalid satellite; value = %s", gmeta->satellite);
        RETURN_ERROR (temp, "GetXMLInput", false);
    }

    if (!strcmp (gmeta->instrument, "TM"))
        this->meta.inst = INST_TM;
    else if (!strncmp (gmeta->instrument, "ETM", 3))
        this->meta.inst = INST_ETM;
    else
    {
        sprintf (temp, "invalid instrument; value = %s", gmeta->instrument);
        RETURN_ERROR (temp, "GetXMLInput", false);
    }

    strcpy (acq_date, gmeta->acquisition_date);
    strcpy (acq_time, gmeta->scene_center_time);
    this->meta.time_fill = false;

    /* Make sure the acquisition time is not too long (i.e. contains too
       many decimal points for the date/time routines).  The time should be
       hh:mm:ss.ssssssZ (see DATE_FORMAT_DATEA_TIME in date.h) which is 16
       characters long.  If the time is longer than that, just chop it off. */
    if (strlen (acq_time) > 16)
        sprintf (&acq_time[15], "Z");

    this->meta.sun_zen = gmeta->solar_zenith;
    if (this->meta.sun_zen < -90.0 || this->meta.sun_zen > 90.0)
    {
        error_string = "solar zenith angle out of range";
        RETURN_ERROR (error_string, "GetXMLInput", false);
    }
    this->meta.sun_zen *= RAD;   /* convert to radians */

    this->meta.sun_az = gmeta->solar_azimuth;
    if (this->meta.sun_az < -360.0 || this->meta.sun_az > 360.0)
    {
        error_string = "solar azimuth angle out of range";
        RETURN_ERROR (error_string, "GetXMLInput", false);
    }
    this->meta.sun_az *= RAD;    /* convert to radians */

    this->meta.earth_sun_dist = gmeta->earth_sun_dist;

    switch (gmeta->wrs_system)
    {
        case 1: this->meta.wrs_sys = WRS_1; break;
        case 2: this->meta.wrs_sys = WRS_2; break;
        default:
            sprintf (temp, "invalid WRS system; value = %d",
                gmeta->wrs_system);
            RETURN_ERROR (temp, "GetXMLInput", false);
    }
    this->meta.ipath = gmeta->wrs_path;
    this->meta.irow = gmeta->wrs_row;

    if (this->meta.inst == INST_TM || this->meta.inst == INST_ETM)
    {
        this->nband = 6;     /* number of reflectance bands */
        this->meta.iband[0] = 1;
        this->meta.iband[1] = 2;
        this->meta.iband[2] = 3;
        this->meta.iband[3] = 4;
        this->meta.iband[4] = 5;
        this->meta.iband[5] = 7;

        /* number of thermal bands; use combination of 6H and 6L for ETM */
        if (this->meta.inst == INST_TM)
            this->nband_th = 1;
        else
            this->nband_th = 2;
        this->meta.iband_th = 6;
    }

    /* Find band 1 and band 6/61 in the input XML file to obtain band-related
       information */
    this->meta.use_toa_refl_consts = false;
    for (i = 0; i < metadata->nbands; i++)
    {
        if (!strcmp (metadata->band[i].name, "b1") &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* Level-1 */
        {
            /* this is the index we'll use for reflectance band info */
            refl_indx = i;

            /* get the band1 info */
            this->meta.rad_gain[0] = metadata->band[i].rad_gain;
            this->meta.rad_bias[0] = metadata->band[i].rad_bias;
            this->file_name[0] = strdup (metadata->band[i].file_name);
            if (!this->file_name[0]) {
                error_string = "allocating memory for filename";
                RETURN_ERROR (error_string, "GetXMLInput", true);
            }

            /* get the production date but only the date portion (yyyy-mm-dd) */
            strncpy (prod_date, metadata->band[i].production_date, 10);
            prod_date[10] = '\0';

            /* are the TOA reflectance coefficients available in XML file? */
            if (existReflGB (metadata))
            {
                this->meta.use_toa_refl_consts = true;
                this->meta.refl_gain[0] = metadata->band[i].refl_gain;
                this->meta.refl_bias[0] = metadata->band[i].refl_bias;
            }
        }
        else if (!strcmp (metadata->band[i].name, "b2") &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* Level-1 */
        {
            /* get the band2 info */
            this->meta.rad_gain[1] = metadata->band[i].rad_gain;
            this->meta.rad_bias[1] = metadata->band[i].rad_bias;
            this->file_name[1] = strdup (metadata->band[i].file_name);
            if (!this->file_name[1]) {
                error_string = "allocating memory for filename";
                RETURN_ERROR (error_string, "GetXMLInput", true);
            }

            if (this->meta.use_toa_refl_consts)
            {
                this->meta.refl_gain[1] = metadata->band[i].refl_gain;
                this->meta.refl_bias[1] = metadata->band[i].refl_bias;
            }
        }
        else if (!strcmp (metadata->band[i].name, "b3") &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* Level-1 */
        {
            /* get the band3 info */
            this->meta.rad_gain[2] = metadata->band[i].rad_gain;
            this->meta.rad_bias[2] = metadata->band[i].rad_bias;
            this->file_name[2] = strdup (metadata->band[i].file_name);
            if (!this->file_name[2]) {
                error_string = "allocating memory for filename";
                RETURN_ERROR (error_string, "GetXMLInput", true);
            }

            if (this->meta.use_toa_refl_consts)
            {
                this->meta.refl_gain[2] = metadata->band[i].refl_gain;
                this->meta.refl_bias[2] = metadata->band[i].refl_bias;
            }
        }
        else if (!strcmp (metadata->band[i].name, "b4") &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* Level-1 */
        {
            /* get the band4 info */
            this->meta.rad_gain[3] = metadata->band[i].rad_gain;
            this->meta.rad_bias[3] = metadata->band[i].rad_bias;
            this->file_name[3] = strdup (metadata->band[i].file_name);
            if (!this->file_name[3]) {
                error_string = "allocating memory for filename";
                RETURN_ERROR (error_string, "GetXMLInput", true);
            }

            if (this->meta.use_toa_refl_consts)
            {
                this->meta.refl_gain[3] = metadata->band[i].refl_gain;
                this->meta.refl_bias[3] = metadata->band[i].refl_bias;
            }
        }
        else if (!strcmp (metadata->band[i].name, "b5") &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* Level-1 */
        {
            /* get the band5 info */
            this->meta.rad_gain[4] = metadata->band[i].rad_gain;
            this->meta.rad_bias[4] = metadata->band[i].rad_bias;
            this->file_name[4] = strdup (metadata->band[i].file_name);
            if (!this->file_name[4]) {
                error_string = "allocating memory for filename";
                RETURN_ERROR (error_string, "GetXMLInput", true);
            }

            if (this->meta.use_toa_refl_consts)
            {
                this->meta.refl_gain[4] = metadata->band[i].refl_gain;
                this->meta.refl_bias[4] = metadata->band[i].refl_bias;
            }
        }
        else if (!strcmp (metadata->band[i].name, "b7") &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* Level-1 */
        {
            /* get the band7 info */
            this->meta.rad_gain[5] = metadata->band[i].rad_gain;
            this->meta.rad_bias[5] = metadata->band[i].rad_bias;
            this->file_name[5] = strdup (metadata->band[i].file_name);
            if (!this->file_name[5]) {
                error_string = "allocating memory for filename";
                RETURN_ERROR (error_string, "GetXMLInput", true);
            }

            if (this->meta.use_toa_refl_consts)
            {
                this->meta.refl_gain[5] = metadata->band[i].refl_gain;
                this->meta.refl_bias[5] = metadata->band[i].refl_bias;
            }
        }

        if (((!strcmp(metadata->band[i].name, "b6") &&
            this->meta.inst == INST_TM) ||
            ((!strcmp(metadata->band[i].name, "b61") ||
            !strcmp(metadata->band[i].name, "b62")) &&
            this->meta.inst == INST_ETM)) &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* Level-1 */
        {
            /* Local indices used to select correct thermal band for TM/ETM */
            int thermal_index;

            /* Select correct thermal index used Input_t struct */
            if (!strcmp(metadata->band[i].name, "b62"))
                thermal_index = BAND_6H;
            else
                thermal_index = BAND_6L;

            /* this is the index we'll use for thermal band xml metadata */
            xml_thermal_index = i;

            /* get the band6 info */
            this->meta.rad_gain_th[thermal_index] = metadata->band[i].rad_gain;
            this->meta.rad_bias_th[thermal_index] = metadata->band[i].rad_bias;
            this->file_name_th[thermal_index] =
                strdup (metadata->band[i].file_name);
            if (!this->file_name_th[thermal_index]) {
                error_string = "allocating memory for filename";
                RETURN_ERROR (error_string, "GetXMLInput", true);
            }

            if (this->meta.use_toa_refl_consts)
            {
                this->meta.k1_const[thermal_index] = metadata->band[i].k1_const;
                this->meta.k2_const[thermal_index] = metadata->band[i].k2_const;
            }
        }
        else if (!strcmp (metadata->band[i].name, "sza"))
        {
            /* get the solar zenith representative band info */
            this->file_name_sun_zen = strdup (metadata->band[i].file_name);
            if (!this->file_name_sun_zen) {
                error_string = "allocating memory for filename";
                RETURN_ERROR (error_string, "GetXMLInput", true);
            }
            if (metadata->band[i].scale_factor != ESPA_FLOAT_META_FILL)
                this->meta.szen_scale = metadata->band[i].scale_factor;
            if (metadata->band[i].add_offset != ESPA_FLOAT_META_FILL)
                this->meta.szen_offset = metadata->band[i].add_offset;
        }
        else if (!strcmp (metadata->band[i].name, "qa_pixel"))
        {
            /* get the L1 qa band file */
            this->file_name_band_qa = strdup (metadata->band[i].file_name);
            if (!this->file_name_band_qa) {
                error_string = "allocating memory for filename";
                RETURN_ERROR (error_string, "GetXMLInput", true);
            }
        }
    }  /* for i */

    if (refl_indx == -99)
    {
        sprintf (temp, "band 1 (b1) was not found in the XML file");
        RETURN_ERROR (temp, "GetXMLInput", false);
    }

    /* Make sure the solar zenith band was found */
    if (this->file_name_sun_zen == NULL)
    {
        sprintf (temp, "Representative band for the solar zenith data was "
            "not found in the XML file.");
        RETURN_ERROR (temp, "GetXMLInput", false);
    }

    /* Make sure the L1 QA band was found */
    if (this->file_name_band_qa == NULL)
    {
        sprintf (temp, "Representative band for the QA data was "
            "not found in the XML file.");
        RETURN_ERROR (temp, "GetXMLInput", false);
    }

    /* Let the user know if the XML params are being used or the LUT params */
    if (this->meta.use_toa_refl_consts)
    {
        printf ("Using the TOA reflectance coefficients, K1/K2 thermal "
            "constants, and earth-sun distance from the XML file.\n");
    }
    else
    {
        printf ("Using the hard-coded TOA reflectance coefficients, K1/K2 "
            "thermal constants, and earth-sun distance table.\n");
    }

    /* Pull the reflectance info from band1 in the XML file */
    this->size.s = metadata->band[refl_indx].nsamps;
    this->size.l = metadata->band[refl_indx].nlines;
    this->size_th.s = metadata->band[xml_thermal_index].nsamps;
    this->size_th.l = metadata->band[xml_thermal_index].nlines;

    /* Check WRS path/rows */
    if (this->meta.wrs_sys == WRS_1)
    {
        if (this->meta.ipath > 251)
            error_string = "WRS path number out of range";
        else if (this->meta.irow > 248)
            error_string = "WRS row number out of range";
    }
    else if (this->meta.wrs_sys == WRS_2)
    {
        if (this->meta.ipath > 233)
            error_string = "WRS path number out of range";
        else if (this->meta.irow > 248)
            error_string = "WRS row number out of range";
    }
    else
        error_string = "invalid WRS system";

    if (error_string != NULL)
    {
        RETURN_ERROR (error_string, "GetXMLInput", true);
    }

    /* Check satellite/instrument combination */
    if (this->meta.inst == INST_MSS)
    {
        if (this->meta.sat != SAT_LANDSAT_1 &&
            this->meta.sat != SAT_LANDSAT_2 &&
            this->meta.sat != SAT_LANDSAT_3 &&
            this->meta.sat != SAT_LANDSAT_4 &&
            this->meta.sat != SAT_LANDSAT_5)
            error_string = "invalid instrument/satellite combination";
    }
    else if (this->meta.inst == INST_TM)
    {
        if (this->meta.sat != SAT_LANDSAT_4 &&
            this->meta.sat != SAT_LANDSAT_5)
            error_string = "invalid instrument/satellite combination";
    }
    else if (this->meta.inst == INST_ETM)
    {
        if (this->meta.sat != SAT_LANDSAT_7)
            error_string = "invalid instrument/satellite combination";
    }
    else
        error_string = "invalid instrument type";

    if (error_string != NULL)
    {
        RETURN_ERROR (error_string, "GetXMLInput", true);
    }

    /* Convert the acquisition date/time values */
    sprintf (temp, "%sT%s", acq_date, acq_time);
    if (!DateInit (&this->meta.acq_date, temp, DATE_FORMAT_DATEA_TIME))
    {
        error_string = "converting acquisition date/time";
        RETURN_ERROR (error_string, "GetXMLInput", false);
    }

    /* Convert the production date value */
    if (!DateInit (&this->meta.prod_date, prod_date, DATE_FORMAT_DATEA))
    {
        error_string = "converting production date";
        RETURN_ERROR (error_string, "GetXMLInput", false);
    }

    return true;
}

