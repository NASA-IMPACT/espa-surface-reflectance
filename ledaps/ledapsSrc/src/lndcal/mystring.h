/*
!C****************************************************************************

!File: mystring.h

!Description: Header file for mystring.c - see mystring.c for more information.

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

 ! Design Notes: (none)
  
!END****************************************************************************
*/

#ifndef MYSTRING_H
#define MYSTRING_H

#include <stdio.h>
#include <stdbool.h>

#define MAX_STR_LEN (256)
#define MAX_NUM_VALUE (20)

/* Key string type definition */

typedef struct {
  int key;
  char *string;
} Key_string_t;

/* Key type definition */

typedef struct {
  char *key;               /* Key string */
  size_t len_key;          /* Length of key */
  int nval;                /* Number of values */
  char *value[MAX_NUM_VALUE];  /* Value strings */
  size_t len_value[MAX_NUM_VALUE];  /* Length of value strings */
} Key_t;


char *DupString(char *string);
int GetLine(FILE *fp, char *s);
bool StringParse(char *s, Key_t *key);
int KeyString(char *key, int len, const Key_string_t *key_string, int null_key, 
              int nkey);

#endif

