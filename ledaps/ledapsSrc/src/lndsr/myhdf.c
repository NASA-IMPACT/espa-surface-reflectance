/*
!C****************************************************************************

!File: myhdf.c
  
!Description: Functions for handling HDF files.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

!Team Unique Header:
  This software was developed by the MODIS Land Science Team Support 
  Group for the Labatory for Terrestrial Physics (Code 922) at the 
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
   1. The following functions handle Science Data Sets (SDSs) and attributes
      in HDF files:

       GetSDSInfo - Read SDS information.
       GetSDSDimInfo - Read SDS dimension information.
       PutSDSInfo - Create an SDS and write information.
       PutSDSDimInfo - Write SDS dimension information.
       GetAttrDouble - Get an HDF attribute's value.

!END****************************************************************************
*/

#include <stdlib.h>
#include "myhdf.h"
#include "error.h"
#include "hdf.h"
#include "mfhdf.h"
#include "mystring.h"

/* Constants */

#define DIM_MAX_NCHAR (80)  /* Maximum size of a dimension name */

/* Possible ranges for data types */

#define MYHDF_CHAR8H     (        255  )
#define MYHDF_CHAR8L     (          0  )
#define MYHDF_INT8H      (        127  )
#define MYHDF_INT8L      (       -128  )
#define MYHDF_UINT8H     (        255  )
#define MYHDF_UINT8L     (          0  )
#define MYHDF_INT16H     (      32767  )
#define MYHDF_INT16L     (     -32768  )
#define MYHDF_UINT16H    (      65535u )
#define MYHDF_UINT16L    (          0u )
#define MYHDF_INT32H     ( 2147483647l )
#define MYHDF_INT32L     (-2147483647l )
#define MYHDF_UINT32H    ( 4294967295ul)
#define MYHDF_UINT32L    (          0ul)
#define MYHDF_FLOAT32H   (3.4028234e+38f)
#define MYHDF_FLOAT32L   (1.1754943e-38f)
#define MYHDF_FLOAT64H   (1.797693134862316e+308)
#define MYHDF_FLOAT64L   (2.225073858507201e-308)

/* Functions */

bool GetSDSInfo(int32 sds_file_id, Myhdf_sds_t *sds)
/* 
!C******************************************************************************

!Description: 'GetSDSInfo' reads information for a specific SDS.
 
!Input Parameters:
 sds_file_id    SDS file id

!Output Parameters:
 sds            SDS data structure; the following fields are updated:
                  index, name, id, rank, type, nattr
 (returns)      Status:
                  'true' = okay
		  'false' = error reading the SDS information

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned if the SDS rank is greater than 
      'MYHDF_MAX_RANK'.
   2. On normal returns the SDS is selected for access.
   3. The HDF file is assumed to be open for SD (Science Data) access.
   4. Error messages are handled with the 'RETURN_ERROR' macro.

!END****************************************************************************
*/
{
  int32 dims[MYHDF_MAX_RANK];

/*  printf("DEBUG: SDS name %s\n",sds->name); */
  sds->index = SDnametoindex(sds_file_id, sds->name);
  if (sds->index == HDF_ERROR)
    RETURN_ERROR("getting sds index", "GetSDSInfo", false);

  sds->id = SDselect(sds_file_id, sds->index);
  if (sds->id == HDF_ERROR)
    RETURN_ERROR("getting sds id", "GetSDSInfo", false);

  if (SDgetinfo(sds->id, sds->name, &sds->rank, dims, 
                &sds->type, &sds->nattr) == HDF_ERROR) {
    SDendaccess(sds->id);
    RETURN_ERROR("getting sds information", "GetSDSInfo", false);
  }
  if (sds->rank > MYHDF_MAX_RANK) {
    SDendaccess(sds->id);
    RETURN_ERROR("sds rank too large", "GetSDSInfo", false);
  }
  return true;
}


bool GetSDSDimInfo(int32 sds_id, Myhdf_dim_t *dim, int irank)
/* 
!C******************************************************************************

!Description: 'GetSDSDimInfo' reads information for a specific SDS dimension.
 
!Input Parameters:
 sds_id         SDS id

!Output Parameters:
 dim            Dimension data structure; the following fields are updated:
                   id, nval, type, nattr, name
 (returns)      Status:
                  'true' = okay
		  'false' = error reading the dimension information

!Team Unique Header:

 ! Design Notes:
   1. The HDF file is assumed to be open for SD (Science Data) access.
   2. An dimension name of less than 'DIM_MAX_NCHAR' is expected.
   3. Error messages are handled with the 'RETURN_ERROR' macro.

!END****************************************************************************
*/
{
  char dim_name[DIM_MAX_NCHAR];

  dim->id = SDgetdimid(sds_id, irank);
  if (dim->id == HDF_ERROR) 
    RETURN_ERROR("getting dimension id", "GetSDSDimInfo", false);

  if (SDdiminfo(dim->id, dim_name,
                &dim->nval, &dim->type, 
	        &dim->nattr) == HDF_ERROR)
      RETURN_ERROR("getting dimension information", "GetSDSDimInfo", false);

  dim->name = DupString(dim_name);
  if (dim->name == (char *)NULL)
    RETURN_ERROR("copying dimension name", "GetSDSDimInfo", false);

  return true;
}

bool PutSDSInfo(int32 sds_file_id, Myhdf_sds_t *sds)
/* 
!C******************************************************************************

!Description: 'PutSDSInfo' creates a SDS and writes SDS information.
 
!Input Parameters:
 sds_file_id    SDS file id
 sds            SDS data structure; the following are used:
                   rank, name, type, dim[*].nval

!Output Parameters:
 sds            SDS data structure; the following are updated:
                   id, index
 (returns)      Status:
                  'true' = okay
		  'false' = error writing the SDS information

!Team Unique Header:

 ! Design Notes:
   1. A maximum of 'MYHDF_MAX_RANK' dimensions are expected.
   2. On normal returns the SDS is selected for access.
   3. The HDF file is assumed to be open for SD (Science Data) access.
   4. Error messages are handled with the 'RETURN_ERROR' macro.

!END****************************************************************************
*/
{
  int irank;
  int32 dims[MYHDF_MAX_RANK];

  for (irank = 0; irank < sds->rank; irank++)
    dims[irank] = sds->dim[irank].nval;

  /* Create the SDS */

  sds->id = SDcreate(sds_file_id, sds->name, sds->type, 
                     sds->rank, dims);
  if (sds->id == HDF_ERROR)
    RETURN_ERROR("Creating sds", "PutSDSInfo", false);

  sds->index = SDnametoindex(sds_file_id, sds->name);
  if (sds->index == HDF_ERROR)
    RETURN_ERROR("Getting sds index", "PutSDSInfo", false);

  return true;
}


bool PutSDSDimInfo(int32 sds_id, Myhdf_dim_t *dim, int irank)
/* 
!C******************************************************************************

!Description: 'PutSDSDimInfo' writes information for a SDS dimension.
 
!Input Parameters:
 sds_id         SDS id
 dim            Dimension data structure; the following field is used:
                   name
 irank          Dimension rank

!Output Parameters:
 dim            Dimension data structure; the following field is updated:
                   id
 (returns)      Status:
                  'true' = okay
		  'false' = error writing the dimension information

!Team Unique Header:

 ! Design Notes:
   1. The HDF file is assumed to be open for SD (Science Data) access.
   2. Error messages are handled with the 'RETURN_ERROR' macro.
   3. Setting the type of the dimension is not currently implemented.

!END****************************************************************************
*/
{

  dim->id = SDgetdimid(sds_id, irank);
  if (dim->id == HDF_ERROR) 
    RETURN_ERROR("getting dimension id", "PutSDSDimInfo", false);

  if (SDsetdimname(dim->id, dim->name) == HDF_ERROR)
    RETURN_ERROR("setting dimension name", "PutSDSDimInfo", false);

  /* Set dimension type */

    /* !! do it !! */
 
  return true;
}


bool GetAttrDouble(int32 sds_id, Myhdf_attr_t *attr, double *val)
/* 
!C******************************************************************************

!Description: 'GetAttrDouble' reads an attribute into a parameter of type
 'double'.
 
!Input Parameters:
 sds_id         SDS id
 attr           Attribute data structure; the following field is used:
                   name

!Output Parameters:
 attr           Attribute data structure; the following field is updated:
                   id, type, nval
 val            An array of values from the HDF attribute (converted from the
                  native type to type 'double'.
 (returns)      Status:
                  'true' = okay
		  'false' = error reading the attribute information

!Team Unique Header:

 ! Design Notes:
   1. The values in the attribute are converted from the stored type to 
      'double' type.
   2. The HDF file is assumed to be open for SD (Science Data) access.
   3. If the attribute has more than 'MYHDF_MAX_NATTR_VAL' values, an error
      status is returned.
   4. Error messages are handled with the 'RETURN_ERROR' macro.
!END****************************************************************************
*/
{
  char8 val_char8[MYHDF_MAX_NATTR_VAL];
  uint8 val_int8[MYHDF_MAX_NATTR_VAL];
  uint8 val_uint8[MYHDF_MAX_NATTR_VAL];
  int16 val_int16[MYHDF_MAX_NATTR_VAL];
  uint16 val_uint16[MYHDF_MAX_NATTR_VAL];
  int32 val_int32[MYHDF_MAX_NATTR_VAL];
  uint32 val_uint32[MYHDF_MAX_NATTR_VAL];
  float32 val_float32[MYHDF_MAX_NATTR_VAL];
  float64 val_float64[MYHDF_MAX_NATTR_VAL];
  int i;
  char z_name[80];
  
  if ((attr->id = SDfindattr(sds_id, attr->name)) == HDF_ERROR)
    RETURN_ERROR("getting attribute id", "GetAttrDouble", false);
  if (SDattrinfo(sds_id, attr->id, z_name, &attr->type, &attr->nval) == 
      HDF_ERROR)
    RETURN_ERROR("getting attribute info", "GetAttrDouble", false);
  /* printf("attr name: %s\n", z_name); */

  if (attr->nval < 1)
    RETURN_ERROR("no attribute value", "GetAttrDouble", false);
  if (attr->nval > MYHDF_MAX_NATTR_VAL) 
    RETURN_ERROR("too many attribute values", "GetAttrDouble", false);

  switch (attr->type) {
  case DFNT_CHAR8:
    if (SDreadattr(sds_id, attr->id, val_char8) == HDF_ERROR) 
      RETURN_ERROR("reading attribute (char8)", "GetAttrDouble", false);
    for (i = 0; i < attr->nval; i++) 
      val[i] = (double)val_char8[i];
    break;
  case DFNT_INT8:
    if (SDreadattr(sds_id, attr->id, val_int8) == HDF_ERROR) 
      RETURN_ERROR("reading attribute (int8)", "GetAttrDouble", false);
    for (i = 0; i < attr->nval; i++) 
      val[i] = (double)val_int8[i];
    break;
  case DFNT_UINT8:
    if (SDreadattr(sds_id, attr->id, val_uint8) == HDF_ERROR) 
      RETURN_ERROR("reading attribute (uint8)", "GetAttrDouble", false);
    for (i = 0; i < attr->nval; i++) 
      val[i] = (double)val_uint8[i];
    break;
  case DFNT_INT16:
    if (SDreadattr(sds_id, attr->id, val_int16) == HDF_ERROR) 
      RETURN_ERROR("reading attribute (int16)", "GetAttrDouble", false);
    for (i = 0; i < attr->nval; i++) 
      val[i] = (double)val_int16[i];
    break;
  case DFNT_UINT16:
    if (SDreadattr(sds_id, attr->id, val_uint16) == HDF_ERROR) 
      RETURN_ERROR("reading attribute (uint16)", "GetAttrDouble", false);
    for (i = 0; i < attr->nval; i++) 
      val[i] = (double)val_uint16[i];
    break;
  case DFNT_INT32:
    if (SDreadattr(sds_id, attr->id, val_int32) == HDF_ERROR) 
      RETURN_ERROR("reading attribute (int32)", "GetAttrDouble", false);
    for (i = 0; i < attr->nval; i++) 
      val[i] = (double)val_int32[i];
    break;
  case DFNT_UINT32:
    if (SDreadattr(sds_id, attr->id, val_uint32) == HDF_ERROR) 
      RETURN_ERROR("reading attribute (uint32)", "GetAttrDouble", false);
    for (i = 0; i < attr->nval; i++) 
      val[i] = (double)val_uint32[i];
    break;
  case DFNT_FLOAT32:
    if (SDreadattr(sds_id, attr->id, val_float32) == HDF_ERROR) 
      RETURN_ERROR("reading attribute (float32)", "GetAttrDouble", false);
    for (i = 0; i < attr->nval; i++) 
      val[i] = (double)val_float32[i];
    break;
  case DFNT_FLOAT64:
    if (SDreadattr(sds_id, attr->id, val_float64) == HDF_ERROR) 
      RETURN_ERROR("reading attribute (float64)", "GetAttrDouble", false);
    for (i = 0; i < attr->nval; i++) 
      val[i] = (double)val_float64[i];
    break;
  default:
    RETURN_ERROR("unknown attribute type", "GetAttrDouble", false);
  }

  return true;
}


bool PutAttrDouble(int32 sds_id, Myhdf_attr_t *attr, double *val)
/* 
!C******************************************************************************

!Description: 'PutAttrDouble' writes an attribute from a parameter of type
 'double' to a HDF file.
 
!Input Parameters:
 sds_id         SDS id
 attr           Attribute data structure; the following fields are used:
                   name, type, nval

!Output Parameters:
 val            An array of values from the HDF attribute (converted from
                  type 'double' to the native type
 (returns)      Status:
                  'true' = okay
		      'false' = error writing the attribute information

!History
 Revision 1.1 2013/1/8
 Gail Schmidt - Modified to not add 0.5 to the floating point attribute values
 before writing to the HDF file.  That is only appropriate when converting
 float to int and not float to float.

!Team Unique Header:

 ! Design Notes:
   1. The values in the attribute are converted from the stored type to 
      'double' type.
   2. The HDF file is assumed to be open for SD (Science Data) access.
   3. If the attribute has more than 'MYHDF_MAX_NATTR_VAL' values, an error
      status is returned.
   4. Error messages are handled with the 'RETURN_ERROR' macro.
!END****************************************************************************
*/
{
  char8 val_char8[MYHDF_MAX_NATTR_VAL];
  int8 val_int8[MYHDF_MAX_NATTR_VAL];
  uint8 val_uint8[MYHDF_MAX_NATTR_VAL];
  int16 val_int16[MYHDF_MAX_NATTR_VAL];
  uint16 val_uint16[MYHDF_MAX_NATTR_VAL];
  int32 val_int32[MYHDF_MAX_NATTR_VAL];
  uint32 val_uint32[MYHDF_MAX_NATTR_VAL];
  float32 val_float32[MYHDF_MAX_NATTR_VAL];
  float64 val_float64[MYHDF_MAX_NATTR_VAL];
  int i;
  void *buf;

  if (attr->nval <= 0  ||  attr->nval > MYHDF_MAX_NATTR_VAL) 
    RETURN_ERROR("invalid number of values", "PutAttrDouble", false);

  switch (attr->type) {
    case DFNT_CHAR8:
      for (i = 0; i < attr->nval; i++) {
        if (     val[i] >= ((double)MYHDF_CHAR8H)) val_char8[i] = MYHDF_CHAR8H;
        else if (val[i] <= ((double)MYHDF_CHAR8L)) val_char8[i] = MYHDF_CHAR8L;
        else if (val[i] >= 0.0) val_char8[i] = (char8)(val[i] + 0.5);
        else                    val_char8[i] = -((char8)(-val[i] + 0.5));
      }
      buf = (void *)val_char8;
      break;

    case DFNT_INT8:
      for (i = 0; i < attr->nval; i++) {
        if (     val[i] >= ((double)MYHDF_INT8H)) val_int8[i] = MYHDF_INT8H;
        else if (val[i] <= ((double)MYHDF_INT8L)) val_int8[i] = MYHDF_INT8L;
        else if (val[i] >= 0.0) val_int8[i] =   (int8)( val[i] + 0.5);
        else                    val_int8[i] = -((int8)(-val[i] + 0.5));
      }
      buf = (void *)val_int8;
      break;

    case DFNT_UINT8:
      for (i = 0; i < attr->nval; i++) {
        if (     val[i] >= ((double)MYHDF_UINT8H)) val_uint8[i] = MYHDF_UINT8H;
        else if (val[i] <= ((double)MYHDF_UINT8L)) val_uint8[i] = MYHDF_UINT8L;
        else if (val[i] >= 0.0) val_uint8[i] =   (uint8)( val[i] + 0.5);
        else                    val_uint8[i] = -((uint8)(-val[i] + 0.5));
      }
      buf = (void *)val_uint8;
      break;

    case DFNT_INT16:
      for (i = 0; i < attr->nval; i++) {
        if (     val[i] >= ((double)MYHDF_INT16H)) val_int16[i] = MYHDF_INT16H;
        else if (val[i] <= ((double)MYHDF_INT16L)) val_int16[i] = MYHDF_INT16L;
        else if (val[i] >= 0.0) val_int16[i] =   (int16)( val[i] + 0.5);
        else                    val_int16[i] = -((int16)(-val[i] + 0.5));
      }
      buf = (void *)val_int16;
      break;

    case DFNT_UINT16:
      for (i = 0; i < attr->nval; i++) {
        if (     val[i] >= ((double)MYHDF_UINT16H)) val_uint16[i] = MYHDF_UINT16H;
        else if (val[i] <= ((double)MYHDF_UINT16L)) val_uint16[i] = MYHDF_UINT16L;
        else if (val[i] >= 0.0) val_uint16[i] =   (uint16)( val[i] + 0.5);
        else                    val_uint16[i] = -((uint16)(-val[i] + 0.5));
      }
      buf = (void *)val_uint16;
      break;

    case DFNT_INT32:
      for (i = 0; i < attr->nval; i++) {
        if (     val[i] >= ((double)MYHDF_INT32H)) val_int32[i] = MYHDF_INT32H;
        else if (val[i] <= ((double)MYHDF_INT32L)) val_int32[i] = MYHDF_INT32L;
        else if (val[i] >= 0.0) val_int32[i] =   (int32)( val[i] + 0.5);
        else                    val_int32[i] = -((int32)(-val[i] + 0.5));
      }
      buf = (void *)val_int32;
      break;

    case DFNT_UINT32:
      for (i = 0; i < attr->nval; i++) {
        if (     val[i] >= ((double)MYHDF_UINT32H)) val_uint32[i] = MYHDF_UINT32H;
        else if (val[i] <= ((double)MYHDF_UINT32L)) val_uint32[i] = MYHDF_UINT32L;
        else if (val[i] >= 0.0) val_uint32[i] =   (uint32)( val[i] + 0.5);
        else                    val_uint32[i] = -((uint32)(-val[i] + 0.5));
      }
      buf = (void *)val_uint32;
      break;

    case DFNT_FLOAT32:
      for (i = 0; i < attr->nval; i++) {
        if (     val[i] >= ((double)MYHDF_FLOAT32H)) val_float32[i] = MYHDF_FLOAT32H;
        else if (val[i] <= ((double)MYHDF_FLOAT32L)) val_float32[i] = MYHDF_FLOAT32L;
        else val_float32[i] = (float32) val[i];
      }
      buf = (void *)val_float32;
      break;

    case DFNT_FLOAT64:
      if (sizeof(float64) == sizeof(double))
        buf = (void *)val;
      else {
        for (i = 0; i < attr->nval; i++)
          val_float64[i] = val[i];
        buf = (void *)val_float64;
      }
      break;

    default: 
      RETURN_ERROR("unimplmented type", "PutAttrDouble", false);
  }

  if (SDsetattr(sds_id, attr->name, attr->type, attr->nval, buf) == HDF_ERROR)
    RETURN_ERROR("setting attribute", "PutAttrDouble", false);

  return true;
}


bool GetAttrString(int32 sds_id, Myhdf_attr_t *attr, char *string)
/* 
!C******************************************************************************

!Description: 'GetAttrString' reads an string (char8) attribute.

!Input Parameters:
 sds_id         SDS id
 attr           Attribute data structure; the following field is used:
                   name

!Output Parameters:
 attr           Attribute data structure; the following field is updated:
                   id, type, nval
 val            An array of values from the HDF attribute (converted from the
                  native type to type 'double'.
 (returns)      Status:
                  'true' = okay
		  'false' = error reading the attribute information

!Team Unique Header:

 ! Design Notes:
   1. The values in the attribute are converted from the stored type to 
      'double' type.
   2. The HDF file is assumed to be open for SD (Science Data) access.
   3. If the attribute has more than 'MYHDF_MAX_NATTR_VAL' values, an error
      status is returned.
   4. Error messages are handled with the 'RETURN_ERROR' macro.
!END****************************************************************************
*/
{
  char8 val_char8[MYHDF_MAX_NATTR_VAL];
  int i,i_length;
  char z_name[80];
  void *buf;
  
  if ((attr->id = SDfindattr(sds_id, attr->name)) == HDF_ERROR)
    RETURN_ERROR("getting attribute id", "GetAttrString", false);
  if (SDattrinfo(sds_id, attr->id, z_name, &attr->type, &attr->nval) == 
      HDF_ERROR)
    RETURN_ERROR("getting attribute info", "GetAttrString", false);
  /* printf("attr name: %s\n", z_name); */

  if (attr->nval < 1)
    RETURN_ERROR("no attribute value", "GetAttrString", false);
  if (attr->nval > MYHDF_MAX_NATTR_VAL) 
    RETURN_ERROR("too many attribute values", "GetAttrString", false);
  if (attr->type != DFNT_CHAR8) 
    RETURN_ERROR("invalid type - not string (char8)", "GetAttrString", false);

  if (sizeof(char8) == sizeof(char))
    buf = (void *)string;
  else
    buf = (void *)val_char8;

  if (SDreadattr(sds_id, attr->id, buf) == HDF_ERROR) 
    RETURN_ERROR("reading attribute", "GetAttrString", false);

  if (sizeof(char8) != sizeof(char)) {
    for (i = 0; i < attr->nval; i++) 
      string[i] = (char)val_char8[i];
  }

  i_length= (int)attr->nval;
  string[i_length]= '\0';

  return true;
}


bool PutAttrString(int32 sds_id, Myhdf_attr_t *attr, char *string)
/* 
!C******************************************************************************

!Description: 'PutAttrString' writes an attribute from a parameter of type
 'double' to a HDF file.
 
!Input Parameters:
 sds_id         SDS id
 attr           Attribute data structure; the following fields are used:
                   name, type, nval

!Output Parameters:
 val            An array of values from the HDF attribute (converted from
                  type 'double' to the native type
 (returns)      Status:
                  'true' = okay
		      'false' = error writing the attribute information

!Team Unique Header:

 ! Design Notes:
   1. The values in the attribute are converted from the stored type to 
      'double' type.
   2. The HDF file is assumed to be open for SD (Science Data) access.
   3. If the attribute has more than 'MYHDF_MAX_NATTR_VAL' values, an error
      status is returned.
   4. Error messages are handled with the 'RETURN_ERROR' macro.
!END****************************************************************************
*/
{
  char8 val_char8[MYHDF_MAX_NATTR_VAL];
  int i;
  void *buf;

  if (attr->nval <= 0  ||  attr->nval > MYHDF_MAX_NATTR_VAL) 
    RETURN_ERROR("invalid number of values", "PutAttrString", false);

  if (attr->type != DFNT_CHAR8) 
    RETURN_ERROR("invalid type -- not string (char8)", "PutAttrString", false);

  if (sizeof(char8) == sizeof(char))
    buf = (void *)string;
  else {
    for (i = 0; i < attr->nval; i++) 
      val_char8[i] = (char8)string[i];
    buf = (void *)val_char8;
  }

  if (SDsetattr(sds_id, attr->name, attr->type, attr->nval, buf) == HDF_ERROR)
    RETURN_ERROR("setting attribute", "PutAttrString", false);

  return true;
}
