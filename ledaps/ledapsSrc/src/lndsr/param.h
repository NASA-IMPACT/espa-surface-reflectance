/*
!C****************************************************************************

!File: param.h

!Description: Header file for 'param.c' - see 'param.c' for more information.

!Revision History:
 Revision 1.0 2000/11/07
 Robert Wolfe
 Original Version.

 Revision 1.1 2000/12/13
 Sadashiva Devadiga
 Modified to accept parameters from command line or file.

 Revision 1.2 2001/05/08
 Sadashiva Devadiga
 Cleanup.

 Revision 1.3 2002/05/10
 Robert Wolfe
 Added separate output SDS name.

 Revision 2.0 2014/02/03
 Gail Schmidt
 Modified to support the ESPA internal raw binary format

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
      phone: 301-614-5508               Lanham, MD 20706  

      Sadashiva Devadiga (Code 922)
      MODIS Land Team Support Group     SSAI
      devadiga@ltpmail.gsfc.nasa.gov    5900 Princess Garden Pkway, #300
      phone: 301-614-5549               Lanham, MD 20706

 ! Design Notes:
   1. Structures are declared for the 'input_space_type' and 'param' data types.
   2. The acronym SDS stands for Science Data Set.
  
!END****************************************************************************
*/

#ifndef PARAM_H
#define PARAM_H

#include <stdbool.h>
#include "lndpm.h"  /* For version number */

/* Parameter data structure type definition */

typedef struct {
  char *param_file_name;      /* Parameter file name */
  char *input_xml_file_name;  /* Input XML metadata file name */
  bool thermal_band;          /* True if thermal band file available */
  char *ncep_file_name[4];    /* Bracketing NCEP file names */
  char *prwv_file_name;       /* Bracketing NCEP hdf file names */
  char *ozon_file_name;       /* Ozone hdf file names */
  char *LEDAPSVersion;        /* LEDAPS Version */
  int  num_ncep_files;        /* number of NCEP files */
  int  num_prwv_files;        /* number of PRWV hdf files */
  int  num_ozon_files;        /* number of Ozone hdf files */
  char *dem_file;             /* DEM file name */
  bool dem_flag;              /* false if not present use default */
} Param_t;

/* Prototypes */

Param_t *GetParam(int argc, char *argv[]);
void FreeParam(Param_t *this);

#endif
