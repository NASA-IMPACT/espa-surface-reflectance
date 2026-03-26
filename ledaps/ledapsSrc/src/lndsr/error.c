/*
!C****************************************************************************

!File: error.c
  
!Description: Function for handling errors.

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
   1. See 'error.h' for information on the 'EXIT_ERROR' and 'ERROR_RETURN'
      macros that automatically populate the source code file name, line number
      and exit flag.  

!END****************************************************************************
*/

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "error.h"

void Error(const char *message, const char *module, 
           const char *source, long line, bool done)
/* 
!C******************************************************************************

!Description: 'Error' writes an error message to 'stderr' and optionally 
  exit's the program with a 'EXIT_FAILURE' status.
 
!Input Parameters:
 message        error message
 module         calling module name
 source         source code file name containing the calling module
 line           line number in the source code file
 done           flag indicating if program is to exit with a failure status;
                'true' = exit, 'false' = return

!Output Parameters: (none)
 (returns)      (void)

!Team Unique Header:

 ! Design Notes:
   1. If the 'errno' flag is set, the 'perror' function is first called to 
      print any i/o related errors.
  
   2. The error message is written to 'stdout'.

   3. The module name, source code name and line number are included in the 
      error message.

!END****************************************************************************
*/
{
  if (errno) perror(" i/o error ");
  fprintf(stderr, " error [%s, %s:%ld] : %s\n", module, source, line, message);
  if (done) exit(EXIT_FAILURE);
  else return;
}
