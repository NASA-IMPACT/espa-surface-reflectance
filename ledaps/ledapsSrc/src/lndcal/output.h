/*
!C****************************************************************************

!File: output.h

!Description: Header file for output.c - see output.c for more information.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

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
   1. Structure is declared for the 'output' data type.
  
!END****************************************************************************
*/

#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "lndcal.h"
#include "input.h"
#include "lut.h"
#include "param.h"
#include "espa_metadata.h"
#include "raw_binary_io.h"

/* Structure for the 'output' data type */

typedef struct {
  bool open;            /* Flag to indicate whether output files are open 
                           for access; 'true' = open, 'false' = not open */
  int nband;            /* Number of bands */
  Img_coord_int_t size; /* Output image size */
  Espa_internal_meta_t metadata;  /* metadata container to hold the band
                           metadata for the output bands; global metadata
                           won't be valid */
  FILE *fp_bin[NBAND_REFL_MAX];  /* File pointer for binary files */
} Output_t;

/* Prototypes */

Output_t *OpenOutput(Espa_internal_meta_t *metadata, Input_t *input,
  Param_t *param, Lut_t *lut, bool thermal);
bool PutOutputLine(Output_t *this, int iband, int iline, void *line);
bool CloseOutput(Output_t *this);
bool FreeOutput(Output_t *this);

#endif
