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
	FreeOutput - Free the 'input' data structure memory.

   2. 'OpenInput' must be called before any of the other routines.  
   3. 'FreeInput' should be used to free the 'input' data structure.
   4. The only input file type supported is HDF.

!END****************************************************************************
*/

#include <stdlib.h>
#include <math.h>
#include "prwv_input.h"
#include "error.h"
#include "mystring.h"
#include "myhdf.h"
#include "const.h"

#define SDS_PREFIX ("band")
#define INPUT_FILL (-9999)

#define INPUT_PROVIDER ("DataProvider")
#define INPUT_SAT ("Satellite")
#define INPUT_INST ("Instrument")
#define INPUT_ACQ_DATE ("AcquisitionDate")
#define INPUT_PROD_DATE ("Level1ProductionDate")
#define INPUT_SUN_ZEN ("SolarZenith")
#define INPUT_SUN_AZ ("SolarAzimuth")
#define INPUT_WRS_SYS ("WRS_System")
#define INPUT_WRS_PATH ("WRS_Path")
#define INPUT_WRS_ROW ("WRS_Row")
#define INPUT_NBAND ("NumberOfBands")
#define INPUT_NBANDS (3)
#define INPUT_NBANDS_OZONE (1)
#define INPUT_BANDS ("BandNumbers")
#define INPUT_DOY ("Day Of Year")
#define INPUT_DATE ("base_date")
#define INPUT_PLATFORM ("Platform")
#define INPUT_SCALEF ("scale_factor")
#define INPUT_ADDOFF ("add_offset")
#define INPUT_NAME1 ("slp")
#define INPUT_NAME2 ("pr_wtr")
#define INPUT_NAME3 ("air")
#define INPUT_NAMEOZ ("ozone")

/* Functions */

InputPrwv_t *OpenInputPrwv(char *file_name)
/* 
!C******************************************************************************

!Description: 'OpenInputPrwv' sets up the 'input' data structure, opens the
 input file for read access.
 
!Input Parameters:
 file_name      input file name
 sds_name       name of sds to be read
 iband          band number for application of band offset
 rank           rank of the input SDS
 dim            dimension flags; the line and sample dimensions are 
                indicated by a -1 and -2 in this array, respectively;
		the index in the other dimensions are indicated by a value 
		of zero or greater
 input_space_type  input space type; either 'SWATH_SPACE' for L1 and L2 input 
                data or 'GRID_SPACE' for L2G, L3 and L4 input data

!Output Parameters:
 dim            dimension flags
 (returns)      'input' data structure or NULL when an error occurs

!Team Unique Header:

 ! Design Notes:
   1. When 'OpenInputPrwv' returns, the file is open for HDF access and the 
      SDS is open for access.
   2. For an input space type of 'GRID_SPACE', a band number of -1 should be 
      given.
   3. The only input HDF data types supported are CHAR8, UNIT8, INT16 and
      UNIT16.
   4. An error status is returned when:
       a. the SDS rank is less than 2 or greater than 'MYHDF_MAX_RANK'
       b. the band number is less than -1 or greater than or equal to
          'NBAND_OFFSET'
       c. either none or more than one dimension is given for the line 
          or sample dimension
       d. an invalid dimension field is given
       e. duplicating strings is not successful
       f. errors occur when opening the input HDF file
       g. errors occur when reading SDS dimensions or attributes
       h. errors occur when opening the SDS for read access
       i. the given SDS rank or dimensions do not match the input file
       j. for a input space type of SWATH_SPACE, the dimensions of a swath 
          are not 1, 2 or 4 times the nominal size of a MODIS swath
       k. for a input space type of SWATH_SPACE, the number of lines is not 
          an integral multiple of the size of the scan at the given resolution
       l. memory allocation is not successful
       m. an invalid input data type is not supported.
   5. Error messages are handled with the 'RETURN_ERROR' macro.
   6. 'FreeInputPrwv' should be called to deallocate memory used by the 
      'input' data structures.
   7. 'CloseInputPrwv' should be called after all of the data is written and 
      before the 'input' data structure memory is released.

!END****************************************************************************
*/
{
  InputPrwv_t *this;
  char *error_string = (char *)NULL;
  char sds_name[40];
  int ir;
  bool sds_open[NBAND_REFL_MAX];
  Myhdf_dim_t *dim[3];
  int ib,osize;
  float *buf = NULL;
  const char* input_names[INPUT_NBANDS]={INPUT_NAME1,INPUT_NAME2,INPUT_NAME3};

  /* Create the Input data structure */

  this = (InputPrwv_t *)malloc(sizeof(InputPrwv_t));
  if (this == (InputPrwv_t *)NULL) 
    RETURN_ERROR("allocating InputPrwv data structure", "OpenInputPrwv", 
                 (InputPrwv_t *)NULL);

  /* Populate the data structure */

  this->file_name = DupString(file_name);
  if (this->file_name == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating file name", "OpenInputPrwv", (InputPrwv_t *)NULL);
  }

  /* Open file for SD access */

  this->sds_file_id = SDstart((char *)file_name, DFACC_RDONLY);
  if (this->sds_file_id == HDF_ERROR) {
    free(this->file_name);
    free(this);  
    RETURN_ERROR("opening input file", "OpenInputPrwv", 
                 (InputPrwv_t *)NULL); 
  }
  this->open = true;

  /* Get the input metadata */

  if (!GetInputPrwvMeta(this)) {
    free(this->file_name);
    free(this);  
    RETURN_ERROR("getting input metadata", "OpenInputPrwv", 
                 (InputPrwv_t *)NULL); 
  }

  /* Get SDS information and start SDS access */

  for (ib = 0; ib < this->nband; ib++) {
    this->sds[ib].name = (char *)NULL;
    this->sds[ib].dim[0].name = (char *)NULL;
    this->sds[ib].dim[1].name = (char *)NULL;
    this->sds[ib].dim[2].name = (char *)NULL;
    sds_open[ib] = false;
    this->buf[ib] = (float *)NULL;
  }

  for (ib = 0; ib < this->nband; ib++) {

    if (sprintf(sds_name, "%s", input_names[ib]) < 0) {
      error_string = "creating SDS name";
      break;
    }

    this->sds[ib].name = DupString(sds_name);
    if (this->sds[ib].name == (char *)NULL) {
      error_string = "setting SDS name";
      break;
    }

    if (!GetSDSInfo(this->sds_file_id, &this->sds[ib])) {
      error_string = "getting sds info";
      break;
    }
    sds_open[ib] = true;

    /* Check rank */

    if (this->sds[ib].rank != 2 && this->sds[ib].rank != 3 ) {
      error_string = "invalid rank";
      break;
    }

    /* Check SDS type */

    if (this->sds[ib].type != DFNT_FLOAT) {
      error_string = "invalid data type -- should be float";
      break;
    }

    /* Get dimensions */

    for (ir = 0; ir < this->sds[ib].rank; ir++) {
      dim[ir] = &this->sds[ib].dim[ir];
      if (!GetSDSDimInfo(this->sds[ib].id, dim[ir], ir)) {
        error_string = "getting dimensions";
        break;
      }
    }
    if (error_string != (char *)NULL) break;

    /* Save and check line and sample dimensions */

    if (ib == 0) {
      this->size.ntime = dim[0]->nval;
      this->size.nlat  = dim[1]->nval;
      this->size.nlon  = dim[2]->nval;
    } else {
      if (this->size.ntime != dim[0]->nval) {
        error_string = "all time dimensions do not match";
        break;
      }
      if (this->size.nlat != dim[1]->nval) {
        error_string = "all lat dimensions do not match";
        break;
      }
      if (this->size.nlon != dim[2]->nval) {
        error_string = "all lat dimensions do not match";
        break;
      }
    } 

    /* Set the scale factor and offset since they don't exist in the global
       attributes */

    this->scale_factor[ib] = 1.0;
    this->add_offset[ib] = 0.0;
  }

  /* Allocate input buffer */

  osize= 3 * ( this->size.ntime * this->size.nlat * this->size.nlon );
  buf = calloc((size_t)(osize), sizeof(float));
  if (buf == NULL)
    error_string = "allocating input buffer";
  else {
    this->buf[0] = buf;
    for (ib = 1; ib < this->nband; ib++)
      this->buf[ib] = this->buf[ib - 1] + 
                    ( this->size.ntime * this->size.nlat * this->size.nlon );
  }

  if (error_string != (char *)NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      for (ir = 0; ir < this->sds[ib].rank; ir++) {
        if (this->sds[ib].dim[ir].name != (char *)NULL)
          free(this->sds[ib].dim[ir].name);
      }
      if (sds_open[ib])
        SDendaccess(this->sds[ib].id);
      if (this->sds[ib].name != (char *)NULL)
        free(this->sds[ib].name);
      if (this->buf[ib] != NULL)
        free(this->buf[ib]);
    }

    SDend(this->sds_file_id);
    free(this->file_name);
    free(this);
    RETURN_ERROR(error_string, "OpenInputPrwv", (InputPrwv_t *)NULL);
  }

  return this;
}


bool CloseInputPrwv(InputPrwv_t *this)
/* 
!C******************************************************************************

!Description: 'CloseInputPrwv' ends SDS access and closes the input file.
 
!InputPrwv Parameters:
 this           'input' data structure; the following fields are input:
                   open, sds.id, sds_file_id

!Output Parameters:
 this           'input' data structure; the following fields are modified:
                   open
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. the file is not open for access
       b. an error occurs when closing access to the SDS.
   2. Error messages are handled with the 'RETURN_ERROR' macro.
   3. 'OpenInputPrwv' must be called before this routine is called.
   4. 'FreeInputPrwv' should be called to deallocate memory used by the 
      'input' data structure.

!END****************************************************************************
*/
{
  int ib;

  if (!this->open)
    RETURN_ERROR("file not open", "CloseInputPrwv", false);

  for (ib = 0; ib < this->nband; ib++) {
    if (SDendaccess(this->sds[ib].id) == HDF_ERROR) 
      RETURN_ERROR("ending sds access", "CloseInputPrwv", false);
  }

  SDend(this->sds_file_id);
  this->open = false;

  return true;
}


void FreeInputPrwv(InputPrwv_t *this)
/* 
!C******************************************************************************

!Description: 'FreeInputPrwv' frees the 'input' data structure memory.
 
!InputPrwv Parameters:
 this           'input' data structure; the following fields are input:
                   sds.rank, sds.dim[*].name, sds.name, file_name

!Team Unique Header:

 ! Design Notes:
   1. 'OpenInputPrwv' and 'CloseInputPrwv' must be called before this routine is called.
   2. An error status is never returned.

!END****************************************************************************
*/
{
  int ib, ir;

  if (this != (InputPrwv_t *)NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      for (ir = 0; ir < this->sds[ib].rank; ir++) {
          free(this->sds[ib].dim[ir].name);
      }
      free(this->sds[ib].name);
    }
    free(this->buf[0]);
    free(this->file_name);
    free(this);
    this = NULL;
  }
}


bool GetInputPrwv(InputPrwv_t *this, int iband, float *read_buffer)
/* 
!C******************************************************************************

!Description: 'WriteOutput' writes a line of data to the output file.
 
!InputPrwv Parameters:
 this           'output' data structure; the following fields are input:
                   open, size, sds.id
 iline          output line number
 read_buffer    buffer of data to be written

!Output Parameters:
 this           'output' data structure; the following fields are modified:
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. the file is not open for access
       b. the line number is invalid (< 0; >= 'this->size.l')
       b. an error occurs when writting to the SDS.
   2. Error messages are handled with the 'RETURN_ERROR' macro.
   3. 'OutputFile' must be called before this routine is called.

!END****************************************************************************
*/
{
  int32 start[MYHDF_MAX_RANK], nval[MYHDF_MAX_RANK];

  /* Check the parameters */

  if (this == (InputPrwv_t *)NULL) 
    RETURN_ERROR("invalid input structure", "GetIntputLine", false);
  if (!this->open)
    RETURN_ERROR("file not open", "GetInputPrwvLine", false);
  if (iband < 0  ||  iband >= this->nband)
    RETURN_ERROR("invalid band number", "GetInputPrwvLine", false);

  /* Read the data */

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  nval[0] = this->size.ntime;
  nval[1] = this->size.nlat;
  nval[2] = this->size.nlon;

  if (SDreaddata(this->sds[iband].id, start, NULL, nval, read_buffer) ==
    HDF_ERROR)
    RETURN_ERROR("reading input", "GetInputPrwvLine", false);

  return true;
}


bool GetInputPrwvMeta(InputPrwv_t *this) 
{
  Myhdf_attr_t attr;
  double dval[NBAND_REFL_MAX];
  InputPrwv_meta_t *meta;

  /* Check the parameters */

  if (!this->open)
    RETURN_ERROR("file not open", "GetInputPrwvMeta", false);

  meta = &this->meta;
  meta->fill = INPUT_FILL;

  /* Read the metadata */

  attr.name = INPUT_DOY;
  if (!GetAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("reading input doy", "GetInputPrwvMeta", false);
  meta->doy = (int)floor(dval[0] + 0.5);

  attr.name = INPUT_DATE;
  if (!GetAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("reading input year", "GetInputPrwvMeta", false);
  meta->year = (int)floor(dval[0] + 0.5);

  this->nband = INPUT_NBANDS;
  return true;
}

/***************************************************************************
 * History:
 *   Updated on 7/11/2013 by Gail Schmidt, USGS EROS LSRD Project
 *   Modified the longitude min/max and delta values to be read from the
 *     lat/long dimension arrays in the HDF files.
 *
 ***************************************************************************/
int get_prwv_anc(t_ncep_ancillary *anc,InputPrwv_t *this, float* data,
  int index)
{
  int i,j,jin,osize;
  float *buffer=NULL;
  int32 sds_idx, sds_id;
  float *dim_buf=NULL;
  int32 start[MYHDF_MAX_RANK], nval[MYHDF_MAX_RANK];

  anc->nblayers = this->size.ntime;
  anc->nbrows   = this->size.nlat;
  anc->nbcols   = this->size.nlon;
  anc->doy      = this->meta.doy;
  anc->year     = this->meta.year;
  anc->timeres  =    6.0; 
  strcpy (anc->source, "NASA GSFC");

/***
  anc->latmin   =  -90.0;
  anc->latmax   =   90.0;
  anc->lonmin   = -180.0;
  anc->lonmax   =  180.0;
  anc->deltalat =    2.5;
  anc->deltalon =    2.5;
***/
	
  /* read the min/max values from the latitude dimension, then calculate the
     delta */
  sds_idx = SDnametoindex (this->sds_file_id, "lat");
  if (sds_idx == FAIL)
    RETURN_ERROR("unable to find lat dimension in prwv file", "get_prwv_anc",
      false);
    
  sds_id = SDselect(this->sds_file_id, sds_idx);
  if (sds_id == FAIL)
    RETURN_ERROR("SDselect failed for lat dimension in prwv file",
      "get_prwv_anc", false);

  dim_buf = (float *) calloc (this->size.nlat, sizeof (float));
  if (dim_buf == NULL) 
    RETURN_ERROR ("allocating prwv latitude buffer", "get_prwv_anc", false);

  start[0] = 0;
  nval[0] = this->size.nlat;
  if (SDreaddata(sds_id, start, NULL, nval, dim_buf) == HDF_ERROR)
    RETURN_ERROR("reading lat data from prwv file", "get_prwv_anc", false);

  anc->latmax = dim_buf[0];
  anc->latmin = dim_buf[this->size.nlat-1];
  anc->deltalat = dim_buf[0] - dim_buf[1];

  if (dim_buf)
    free (dim_buf);

  /* read the min/max values from the longitude dimension, then calculate the
     delta */
  sds_idx = SDnametoindex (this->sds_file_id, "lon");
  if (sds_idx == FAIL)
    RETURN_ERROR("unable to find long dimension in prwv file", "get_prwv_anc",
      false);
    
  sds_id = SDselect(this->sds_file_id, sds_idx);
  if (sds_id == FAIL)
    RETURN_ERROR("SDselect failed for long dimension in prwv file",
      "get_prwv_anc", false);

  dim_buf = (float *) calloc (this->size.nlon, sizeof (float));
  if (dim_buf == NULL) 
    RETURN_ERROR ("allocating prwv longitude buffer", "get_prwv_anc", false);

  start[0] = 0;
  nval[0] = this->size.nlon;
  if (SDreaddata(sds_id, start, NULL, nval, dim_buf) == HDF_ERROR)
    RETURN_ERROR("reading lon data from prwv file", "get_prwv_anc", false);

  anc->lonmin = dim_buf[0];
  anc->lonmax = dim_buf[this->size.nlon-1];
  anc->deltalon = dim_buf[1] - dim_buf[0];

  if (dim_buf)
    free (dim_buf);

  /* If longitude values in the prwv file are from 0 to 360, then make them
     from -180 to 180 */
  if (anc->lonmax > 180)
  {
    anc->lonmin -= 180.0;
    anc->lonmax -= 180.0;
  }

  /* Set up the data buffer */
  osize= anc->nblayers * anc->nbrows * anc->nbcols;

  buffer = (float *)calloc((size_t)(osize),sizeof(float));
  if (buffer == (float *)NULL) 
      RETURN_ERROR("allocating prwv buffer", "get_prwv_anc", false);

  anc->data[0]= buffer;
  for ( i=1; i<anc->nblayers; i++)
    anc->data[i]=  anc->data[i-1] + (anc->nbrows * anc->nbcols);

  jin=0;
  for ( i=0; i<anc->nblayers; i++)
    {
    anc->time[i]= (24.0/(float)anc->nblayers)*(float)i;
    for ( j=0; j<(anc->nbrows * anc->nbcols); j++)
      {
      anc->data[i][j]= 
        ((float)data[jin]*this->scale_factor[index])+this->add_offset[index];
      jin++;
      }
    }

  return true;
}

/*********************************************************************/
/*********************************************************************/
/************************** Ozone functions **************************/
/*********************************************************************/
/*********************************************************************/

InputOzon_t *OpenInputOzon(char *file_name)
/* 
!C******************************************************************************

!Description: 'OpenInputOzon' sets up the 'input' data structure, opens the
 input file for read access.
 
!Input Parameters:
 file_name      input file name
 sds_name       name of sds to be read
 iband          band number for application of band offset
 rank           rank of the input SDS
 dim            dimension flags; the line and sample dimensions are 
                indicated by a -1 and -2 in this array, respectively;
		the index in the other dimensions are indicated by a value 
		of zero or greater
 input_space_type  input space type; either 'SWATH_SPACE' for L1 and L2 input 
                data or 'GRID_SPACE' for L2G, L3 and L4 input data

!Output Parameters:
 dim            dimension flags
 (returns)      'input' data structure or NULL when an error occurs

!Team Unique Header:

 ! Design Notes:
   1. When 'OpenInputOzon' returns, the file is open for HDF access and the 
      SDS is open for access.
   2. For an input space type of 'GRID_SPACE', a band number of -1 should be 
      given.
   3. The only input HDF data types supported are CHAR8, UNIT8, INT16 and
      UNIT16.
   4. An error status is returned when:
       a. the SDS rank is less than 2 or greater than 'MYHDF_MAX_RANK'
       b. the band number is less than -1 or greater than or equal to
          'NBAND_OFFSET'
       c. either none or more than one dimension is given for the line 
          or sample dimension
       d. an invalid dimension field is given
       e. duplicating strings is not successful
       f. errors occur when opening the input HDF file
       g. errors occur when reading SDS dimensions or attributes
       h. errors occur when opening the SDS for read access
       i. the given SDS rank or dimensions do not match the input file
       j. for a input space type of SWATH_SPACE, the dimensions of a swath 
          are not 1, 2 or 4 times the nominal size of a MODIS swath
       k. for a input space type of SWATH_SPACE, the number of lines is not 
          an integral multiple of the size of the scan at the given resolution
       l. memory allocation is not successful
       m. an invalid input data type is not supported.
   5. Error messages are handled with the 'RETURN_ERROR' macro.
   6. 'FreeInputOzon' should be called to deallocate memory used by the 
      'input' data structures.
   7. 'CloseInputOzon' should be called after all of the data is written and 
      before the 'input' data structure memory is released.

!END****************************************************************************
*/
{
  InputOzon_t *this;
  char *error_string = (char *)NULL;
  char sds_name[40];
  int ir;
  bool sds_open[NBAND_REFL_MAX];
  Myhdf_dim_t *dim[3];
  int ib,osize;
  int16 *buf;
  const char* input_name={INPUT_NAMEOZ};
  Myhdf_attr_t attr;
  double dval[NBAND_REFL_MAX];

  /* Create the Input data structure */

  this = (InputOzon_t *)malloc(sizeof(InputOzon_t));
  if (this == (InputOzon_t *)NULL) 
    RETURN_ERROR("allocating InputOzon data structure", "OpenInputOzon", 
                 (InputOzon_t *)NULL);

  /* Populate the data structure */

  this->file_name = DupString(file_name);
  if (this->file_name == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating file name", "OpenInputOzon", (InputOzon_t *)NULL);
  }

  /* Open file for SD access */

  this->sds_file_id = SDstart((char *)file_name, DFACC_RDONLY);
  if (this->sds_file_id == HDF_ERROR) {
    free(this->file_name);
    free(this);  
    RETURN_ERROR("opening input file", "OpenInputOzon", 
                 (InputOzon_t *)NULL); 
  }
  this->open = true;

  /* Get the input metadata */

  if (!GetInputOzonMeta(this)) {
    free(this->file_name);
    free(this);  
    RETURN_ERROR("getting input metadata", "OpenInputOzon", 
                 (InputOzon_t *)NULL); 
  }

  /* Get SDS information and start SDS access */

  for (ib = 0; ib < this->nband; ib++) {
    this->sds[ib].name = (char *)NULL;
    this->sds[ib].dim[0].name = (char *)NULL;
    this->sds[ib].dim[1].name = (char *)NULL;
    this->sds[ib].dim[2].name = (char *)NULL;
    sds_open[ib] = false;
    this->buf[ib] = (int16 *)NULL;
  }

  for (ib = 0; ib < this->nband; ib++) {
    if (sprintf(sds_name, "%s", input_name) < 0) {
      error_string = "creating SDS name";
      break;
    }

    this->sds[ib].name = DupString(sds_name);
    if (this->sds[ib].name == (char *)NULL) {
      error_string = "setting SDS name";
      break;
    }

    if (!GetSDSInfo(this->sds_file_id, &this->sds[ib])) {
      error_string = "getting sds info";
      break;
    }
    sds_open[ib] = true;

    /* Check rank */

    if (this->sds[ib].rank != 2 && this->sds[ib].rank != 3 ) {
      error_string = "invalid rank";
      break;
    }

    /* Check SDS type */

    if (this->sds[ib].type != DFNT_INT16) {
      error_string = "invalid number type";
      break;
    }

    /* Get dimensions */

    for (ir = 0; ir < this->sds[ib].rank; ir++) {
      dim[ir] = &this->sds[ib].dim[ir];
      if (!GetSDSDimInfo(this->sds[ib].id, dim[ir], ir)) {
        error_string = "getting dimensions";
        break;
      }
    }
    if (error_string != (char *)NULL) break;

    /* Save and check line and sample dimensions */

    if (ib == 0) {
      this->size.ntime = 1;
      this->size.nlat  = dim[0]->nval;
      this->size.nlon  = dim[1]->nval;
    } else {
        error_string = "ozone should only have one band";
        break;
    }

    attr.name = INPUT_SCALEF;
    if (!GetAttrDouble(this->sds[ib].id, &attr, dval))
      RETURN_ERROR("reading input scale", "OpenInputOzon", false);

    this->scale_factor[ib]= dval[0];

    attr.name = INPUT_ADDOFF;
    if (!GetAttrDouble(this->sds[ib].id, &attr, dval))
      RETURN_ERROR("reading add offset", "OpenInputOzon", false);

    this->add_offset[ib]= dval[0];
  }

  /* Allocate input buffer */

  osize= ( this->size.ntime * this->size.nlat * this->size.nlon );
             
  if (sizeof(int16) != sizeof(int)) {
    buf = (int16 *)calloc((size_t)(osize), sizeof(int16));
    if (buf == (int16 *)NULL)
      error_string = "allocating input buffer";
    else {
      this->buf[0] = buf;
    }
  }

  if (error_string != (char *)NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      for (ir = 0; ir < this->sds[ib].rank; ir++) {
        if (this->sds[ib].dim[ir].name != (char *)NULL)
          free(this->sds[ib].dim[ir].name);
      }
      if (sds_open[ib])
        SDendaccess(this->sds[ib].id);
      if (this->sds[ib].name != (char *)NULL)
        free(this->sds[ib].name);
      if (this->buf[ib] != (int16 *)NULL)
        free(this->buf[ib]);
    }

    SDend(this->sds_file_id);
    free(this->file_name);
    free(this);
    RETURN_ERROR(error_string, "OpenInputOzon", (InputOzon_t *)NULL);
  }

  return this;
}


bool CloseInputOzon(InputOzon_t *this)
/* 
!C******************************************************************************

!Description: 'CloseInputOzon' ends SDS access and closes the input file.
 
!InputOzon Parameters:
 this           'input' data structure; the following fields are input:
                   open, sds.id, sds_file_id

!Output Parameters:
 this           'input' data structure; the following fields are modified:
                   open
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. the file is not open for access
       b. an error occurs when closing access to the SDS.
   2. Error messages are handled with the 'RETURN_ERROR' macro.
   3. 'OpenInputOzon' must be called before this routine is called.
   4. 'FreeInputOzon' should be called to deallocate memory used by the 
      'input' data structure.

!END****************************************************************************
*/
{
  int ib;

  if (!this->open)
    RETURN_ERROR("file not open", "CloseInputOzon", false);

  for (ib = 0; ib < this->nband; ib++) {
    if (SDendaccess(this->sds[ib].id) == HDF_ERROR) 
      RETURN_ERROR("ending sds access", "CloseInputOzon", false);
  }

  SDend(this->sds_file_id);
  this->open = false;

  return true;
}


void FreeInputOzon(InputOzon_t *this)
/* 
!C******************************************************************************

!Description: 'FreeInputOzon' frees the 'input' data structure memory.
 
!InputOzon Parameters:
 this           'input' data structure; the following fields are input:
                   sds.rank, sds.dim[*].name, sds.name, file_name

!Team Unique Header:

 ! Design Notes:
   1. 'OpenInputOzon' and 'CloseInputOzon' must be called before this routine is called.
   2. An error status is never returned.

!END****************************************************************************
*/
{
  int ib, ir;

  if (this != (InputOzon_t *)NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      for (ir = 0; ir < this->sds[ib].rank; ir++) {
          free(this->sds[ib].dim[ir].name);
      }
      free(this->sds[ib].name);
    }
    free(this->buf[0]);
    free(this->file_name);
    free(this);
    this = NULL;
  }
}


bool GetInputOzon(InputOzon_t *this, int iband, int *read_buffer)
/* 
!C******************************************************************************

!Description: 'WriteOutput' writes a line of data to the output file.
 
!InputOzon Parameters:
 this           'output' data structure; the following fields are input:
                   open, size, sds.id
 iline          output line number
 buf            buffer of data to be written

!Output Parameters:
 this           'output' data structure; the following fields are modified:
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. the file is not open for access
       b. the line number is invalid (< 0; >= 'this->size.l')
       b. an error occurs when writting to the SDS.
   2. Error messages are handled with the 'RETURN_ERROR' macro.
   3. 'OutputFile' must be called before this routine is called.

!END****************************************************************************
*/
{
  int32 start[MYHDF_MAX_RANK], nval[MYHDF_MAX_RANK];
  void *buf;
  int is,read_size;
  int32 sds_idx, sds_id;

  /* Check the parameters */
  if (this == (InputOzon_t *)NULL) 
    RETURN_ERROR("invalid input structure", "GetIntputLine", false);
  if (!this->open)
    RETURN_ERROR("file not open", "GetInputOzonLine", false);
  if (iband < 0  ||  iband >= this->nband)
    RETURN_ERROR("invalid band number", "GetInputOzonLine", false);

  if (sizeof(int16) == sizeof(int))
    buf = (void *)read_buffer;
  else
    buf = (void *)this->buf[iband];

  /* Verify the time layer is only one in size, otherwise this code is
     set up incorrectly */
  if (this->size.ntime != 1)
    RETURN_ERROR("invalid number of time layers (expecting only 1)",
      "GetInputOzonLine", false);

  /* Find the ozone band */
  sds_idx = SDnametoindex (this->sds_file_id, "ozone");
  if (sds_idx == FAIL)
    RETURN_ERROR("unable to find ozone band in ozone file", "GetInputOzonLine",
      false);
    
  sds_id = SDselect(this->sds_file_id, sds_idx);
  if (sds_id == FAIL)
    RETURN_ERROR("SDselect failed for ozone band in ozone file",
      "GetInputOzonLine", false);

  /* Read the ozone data */
  start[0] = 0;
  start[1] = 0;
  nval[0] = this->size.nlat;
  nval[1] = this->size.nlon;
  read_size= this->size.nlat * this->size.nlon;

  if (SDreaddata(sds_id, start, NULL, nval, buf) == HDF_ERROR)
    RETURN_ERROR("reading input", "GetInputOzonLine", false);

  if (sizeof(int16) != sizeof(int)) {
    for (is = 0; is < read_size; is++) 
      read_buffer[is] = (int)this->buf[iband][is];
  }

  return true;
}


bool GetInputOzonMeta(InputOzon_t *this)
{
  Myhdf_attr_t attr;
  double dval[NBAND_REFL_MAX];
  InputPrwv_meta_t *meta;

  /* Check the parameters */

  if (!this->open)
    RETURN_ERROR("file not open", "GetInputOzonMeta", false);

  meta = &this->meta;
  meta->fill = INPUT_FILL;

  /* Read the metadata */

  attr.name = INPUT_DOY;
  if (!GetAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("reading input doy", "GetInputOzonMeta", false);
  meta->doy = (int)floor(dval[0] + 0.5);

  attr.name = INPUT_DATE;
  if (!GetAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("reading input year", "GetInputOzonMeta", false);
  meta->year = (int)floor(dval[0] + 0.5);

  attr.name = INPUT_PLATFORM;
  if (!GetAttrString(this->sds_file_id, &attr, meta->platform))
    RETURN_ERROR("reading input platform", "GetInputOzonMeta", false);

  this->nband = INPUT_NBANDS_OZONE;
  return true;
}

/***************************************************************************
 * History:
 *   Updated on 7/10/2013 by Gail Schmidt, USGS EROS LSRD Project
 *   Modified the longitude min/max and delta values to be read from the
 *     lat/long dimension arrays in the HDF files.  These values do change
 *     for OMI vs. pre-OMI platforms.
 *
 ***************************************************************************/
int get_ozon_anc(t_ncep_ancillary *anc,InputOzon_t *this, int* data, int index)
{
  int i,j,jin,osize;
  float *buffer=NULL;
  int32 sds_idx, sds_id;
  float *dim_buf=NULL;
  int32 start[MYHDF_MAX_RANK], nval[MYHDF_MAX_RANK];

  anc->nblayers = this->size.ntime;
  anc->nbrows   = this->size.nlat;
  anc->nbcols   = this->size.nlon;
  anc->doy      = this->meta.doy;
  anc->year     = this->meta.year;
  anc->timeres  =   24.0; 
  strcpy (anc->source, this->meta.platform);

  /* read the min/max values from the latitude dimension, then calculate the
     delta */
  sds_idx = SDnametoindex (this->sds_file_id, "lat");
  if (sds_idx == FAIL)
    RETURN_ERROR("unable to find lat dimension in ozone file", "get_ozon_anc",
      false);
    
  sds_id = SDselect(this->sds_file_id, sds_idx);
  if (sds_id == FAIL)
    RETURN_ERROR("SDselect failed for lat dimension in ozone file",
      "get_ozon_anc", false);

  dim_buf = (float *) calloc (this->size.nlat, sizeof (float));
  if (dim_buf == NULL) 
    RETURN_ERROR ("allocating ozon latitude buffer", "get_ozon_anc", false);

  start[0] = 0;
  nval[0] = this->size.nlat;
  if (SDreaddata(sds_id, start, NULL, nval, dim_buf) == HDF_ERROR)
    RETURN_ERROR("reading lat data from ozone file", "get_ozon_anc", false);

  /* There are some ozone HDF files that were interpolated and not reprocessed,
     thus they will have the latitude values incorrectly represented (upside
     down).  Those latitude values will need to be flipped, as the ozone data
     really goes from north to south. */
  if (dim_buf[0] > 0)
  {
    anc->latmax = dim_buf[0];
    anc->latmin = dim_buf[this->size.nlat-1];
    anc->deltalat = dim_buf[0] - dim_buf[1];
  }
  else
  {
    anc->latmin = dim_buf[0];
    anc->latmax = dim_buf[this->size.nlat-1];
    anc->deltalat = dim_buf[1] - dim_buf[0];
  }

  if (dim_buf)
    free (dim_buf);

  /* read the min/max values from the longitude dimension, then calculate the
     delta */
  sds_idx = SDnametoindex (this->sds_file_id, "lon");
  if (sds_idx == FAIL)
    RETURN_ERROR("unable to find long dimension in ozone file", "get_ozon_anc",
      false);
    
  sds_id = SDselect(this->sds_file_id, sds_idx);
  if (sds_id == FAIL)
    RETURN_ERROR("SDselect failed for long dimension in ozone file",
      "get_ozon_anc", false);

  dim_buf = (float *) calloc (this->size.nlon, sizeof (float));
  if (dim_buf == NULL) 
    RETURN_ERROR ("allocating ozon longitude buffer", "get_ozon_anc", false);

  start[0] = 0;
  nval[0] = this->size.nlon;
  if (SDreaddata(sds_id, start, NULL, nval, dim_buf) == HDF_ERROR)
    RETURN_ERROR("reading lon data from ozone file", "get_ozon_anc", false);

  anc->lonmin = dim_buf[0];
  anc->lonmax = dim_buf[this->size.nlon-1];
  anc->deltalon = dim_buf[1] - dim_buf[0];

  if (dim_buf)
    free (dim_buf);

  /* handle the ozone data differently for the pre-OMI and OMI platforms */
/***
  anc->latmin   =  -89.5;
  anc->latmax   =   89.5;
  anc->deltalat =    1.0;
  if (!strcmp (this->meta.platform, "OMI"))
  { // OMI
    anc->lonmin   = -179.5;
    anc->lonmax   =  179.5;
    anc->deltalon =    1.0;
  }
  else
  { // pre-OMI
    anc->lonmin   = -179.375;
    anc->lonmax   =  179.375;
    anc->deltalon =    1.25;
  }
***/
	
  osize= anc->nblayers * anc->nbrows * anc->nbcols;
  buffer = (float *)calloc((size_t)(osize),sizeof(float));
  if (buffer == (float *)NULL) 
    RETURN_ERROR("allocating ozon buffer", "get_ozon_anc", false);

  anc->data[0]= buffer;
  for ( i=1; i<anc->nblayers; i++)
    anc->data[i]=  anc->data[i-1] + (anc->nbrows * anc->nbcols);

  jin=0;
  for ( i=0; i<anc->nblayers; i++)
    {
    anc->time[i]= (24.0/(float)anc->nblayers)*(float)i;
    for ( j=0; j<(anc->nbrows * anc->nbcols); j++)
      {
      anc->data[i][j]= 
        ((float)data[jin]*this->scale_factor[index])+this->add_offset[index];
      jin++;
      }
    }

  return true;
}
