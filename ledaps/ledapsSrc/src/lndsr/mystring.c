/*
!C****************************************************************************

!File: mystring.c
  
!Description: Functions for handling strings.

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
   1. The following functions handle strings:

       DupString - Duplicate a string.

!END****************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "mystring.h"
#include "error.h"

char *DupString(char *string) 
/* 
!C******************************************************************************

!Description: 'DupString' duplicates a string.
 
!Input Parameters:
 string         string to duplicate

!Output Parameters:
 (returns)      duplicated string or 
                NULL when an error occurs

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. a null string pointer is input
       b. the input string length is invalid (< 0)
       c. there is an erorr allocating memory for the string
       d. there is an error copying the string.
   2. Memory is allocated for the duplicated string.
   3. Error messages are handled with the 'RETURN_ERROR' macro.

!END****************************************************************************
*/
{
  int len;
  char *s;

  if (string == (char *)NULL) return ((char *)NULL);
  len = strlen(string);
  if (len < 0) 
    RETURN_ERROR("invalid string length", "DupString", (char *)NULL);

  s = (char *)calloc((len + 1), sizeof(char));
  if (s == (char *)NULL) 
    RETURN_ERROR("allocating memory", "DupString", (char *)NULL);
  if (strcpy(s, string) != s) {
    free(s);
    RETURN_ERROR("copying string", "DupString", (char *)NULL);
  }

  return (s);
}


int GetLine(FILE *fp, char *s)
/* 
!C******************************************************************************

!Description: 'GetLine' reads a line from a file. 
 
!Input Parameters:
 fp             file handle

!Output Parameters:
 s              the line read from the file
 (returns)      length of the line or zero for end-of-file

!Team Unique Header:

 ! Design Notes:
   1. A maximum of 'MAX_STR_LEN' characters is returned.

!END****************************************************************************
*/
{
  int i, c;

  i = 0;
  while ((c = fgetc(fp)) != EOF)
  {
    s[i++] = c;
    if (c == '\n') break;
    if (i >= MAX_STR_LEN) break;
  }
  s[i-1] = '\0';
  return i;
}                 

bool StringParse(char *s, Key_t *key)
/* 
!C******************************************************************************

!Description: 'StringParse' parsing an input "key = value" string with multiple 
              values seperated by commas or spaces. 
	      Quoted values are allowed.  Keyword only strings are allowed 
	      (zero values).
 
!Input Parameters:
 s              input string
 len            length of string
 c              character to fine

!Output Parameters:
 key            key and values
 (returns)      status; true = okay; false = error
                    

!Team Unique Header:

 ! Design Notes:
   1. The key values must be after the equals ('=') character.
   2. A maximum of 'MAX_NUM_VALUE' values are returned.
   3. An error is returned if more than 'MAX_NUM_VALUE' values are 
      given.
   4. The option values can either be separated by a comma (',') or by 
      spaces (' ').
   5. Separator characters are removed.
   6. Each option value (including extra spaces) can be at most 
      'MAX_STR_LEN - 1' characters.
   7. Cmma (','), qoutes ('"') and/or spaces (' ') are not allowed in the key.

!END****************************************************************************
*/
{
  typedef enum {
    S0,                      /* start state */
    C0,                      /* comment state */
    K0, K1, K2,              /* key states */
    V0, V1, V2, V3, V4, V5,  /* value states */
    Q0, Q1, Q2, Q3, Q4,      /* quoted string states */
    E0, E1, E2, E3, E4, E5,  /* error states */
    X0                       /* end state */
  } state_t;

  /* State table */
  state_t state_table[16][7] = {
  /* '#'  ' '  ','  '"'  "="  end  other     state */
  /* ---  ---  ---  ---  ---  ---  ----- */
    {C0,  S0,  E0,  E0,  E1,  X0,  K0},  /* S0 - start */
    {C0,  C0,  C0,  C0,  C0,  X0,  C0},  /* C0 - comment */
    {C0,  K2,  E0,  E0,  V0,  X0,  K1},  /* K0 - begin key */
    {C0,  K2,  E0,  E0,  V0,  X0,  K1},  /* K1 - non-blank in key */
    {C0,  K2,  E0,  E0,  V0,  X0,  E2},  /* K2 - blank at end of key */
    {C0,  V3,  V0,  Q0,  E3,  X0,  V1},  /* V0 - begin value */
    {C0,  V4,  V0,  E3,  E3,  X0,  V2},  /* V1 - begin non-blank value (a) */
    {C0,  V4,  V0,  E3,  E3,  X0,  V2},  /* V2 - non-blank in value */
    {C0,  V3,  V0,  Q0,  E3,  X0,  V1},  /* V3 - blank before value */
    {C0,  V4,  V0,  Q1,  E3,  X0,  V5},  /* V4 - blank after value */
    {C0,  V4,  V0,  E3,  E3,  X0,  V2},  /* V5 - begin non-blank value (b) */
    {Q2,  Q2,  Q2,  Q4,  Q2,  E4,  Q2},  /* Q0 - begin quoted string (a) */
    {Q2,  Q2,  Q2,  Q4,  Q2,  E4,  Q2},  /* Q1 - begin quoted string (b) */
    {Q3,  Q3,  Q3,  Q4,  Q3,  E4,  Q3},  /* Q2 - continue quoted string (a) */
    {Q3,  Q3,  Q3,  Q4,  Q3,  E4,  Q3},  /* Q3 - continue quoted string (b) */
    {C0,  V4,  V0,  E5,  E5,  X0,  E5},  /* Q4 - end quoted string (b) */
  };

  int nc, i, n1;
  int is;
  state_t state;

  key->key = (char *)NULL;
  key->len_key = 0;
  key->nval = 0;
  for (i = 0; i < MAX_NUM_VALUE; i++) {
    key->value[i] = (char *)NULL;
    key->len_value[i] = 0;
  }

  nc = (int)strlen(s);
  if (s[nc - 1] == (char)'\n') nc--;
  if (nc < 1) return true;

  key->nval = 0;
  n1 = -1;
  state = S0;
  for (i = 0; i < (nc + 1); i++) {
    if (i == nc) is = 5;
    else {
      switch ((int)*s) {
        case ('#'): is = 0; break;
        case (' '): is = 1; break;
        case (','): is = 2; break;
        case ('"'): is = 3; break;
        case ('='): is = 4; break;
        default: is = 6; break;
      }
    }
    state = state_table[state][is];
    if (state == X0) break;

    switch (state) {
      case S0: break;
      case C0: break;
      case K0: /* begin key */
	key->key = s;
	key->len_key = 1;
	break;
      case K1: /* non-blank in key */
        key->len_key++;
        break;
      case K2: break;
      case Q0: break;
      case V0: /* begin string */
      case Q1: /* begin quoted string (b) */
	if (key->nval >= MAX_NUM_VALUE) 
	  RETURN_ERROR("too many values", "StringParse", false);
        key->nval++;
	n1 = key->nval - 1;
        break;
      case V1: /* begin non-blank value (a) */
      case Q2: /* continue quoted string (b) */
	key->value[n1] = s;
	key->len_value[n1] = 1;
	break;
      case V2: /* non-blank in value */
      case Q3: /* continue quoted string (b) */
	key->len_value[n1]++;
	break;
      case V3: break;
      case Q4: break;
      case V4: break;
      case V5: /* begin non-blank value (b) */
	if (key->nval >= MAX_NUM_VALUE) 
	  RETURN_ERROR("too many values", "StringParse", false);
        key->nval++;
	n1 = key->nval - 1;
	key->value[n1] = s;
	key->len_value[n1] = 1;
	break;
      case E0: 
        RETURN_ERROR("invalid character in key", "StringParse", false);
      case E1: 
        RETURN_ERROR("no key", "StringParse", false);
      case E2: 
        RETURN_ERROR("blank in key", "StringParse", false);
      case E3: 
        RETURN_ERROR("invalid character in value", "StringParse", false);
      case E4: 
        RETURN_ERROR("no end-quote", "StringParse", false);
      case E5: 
        RETURN_ERROR("invalid character after quote", "StringParse", false);
      case X0: break;
      default:
        RETURN_ERROR("invalid state", "StringParse", false);
    }

    s++;
  }
  return true;
}

int KeyString(char *key, int len, const Key_string_t *key_string, int null_key, 
              int nkey) {
  int len_key, len_test;
  int ikey, i;

  len_key = strlen(key);
  if (len < len_key) len_key = len;

  if (nkey < 1) 
    RETURN_ERROR("invalid list", "KeyString", null_key);
  for (ikey = 0; ikey < nkey; ikey++) {
    len_test = strlen(key_string[ikey].string);
    if (len_test != len_key) continue;
    for (i = 0; i < len_key; i++) 
      if (toupper(key[i]) != toupper(key_string[ikey].string[i])) break;
    if (i >= len_key) break;
  }
  if (ikey >= nkey) 
    RETURN_ERROR("string not found in list", "KeyString", null_key);

  return (key_string[ikey].key);
}
