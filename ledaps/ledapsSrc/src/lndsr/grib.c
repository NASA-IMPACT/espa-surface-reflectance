/*
!C******************************************************************

!Description:   grb_routines.c

	Routines for reading data from files in the GRIB format.  The
	primary routine, read_grib_array(), was adapted from the main()
	of a standalone program, wgrib, written by Wesley Ebisuzaki.  
	All routines are from wgrib.c, and most have been edited to
	comment out printf()s or fprintf()s.


!Input Parameters:
!Output Parameters: 
!Revision History:
!Team-unique Header:
!References and Credits 
!Design Notes

 This is a copy of the original file history comments.
 
 *
 * wgrib.c extract/inventory grib records
 *
 *                              Wesley Ebisuzaki
 *
 * 11/94 - v1.0
 * 11/94 - v1.1: arbitary size grids, -i option
 * 11/94 - v1.2: bug fixes, ieee option, more info
 * 1/95  - v1.2.4: fix headers for SUN acc
 * 2/95  - v1.2.5: add num_ave in -s listing
 * 2/95  - v1.2.6: change %d to %ld
 * 2/95  - v1.2.7: more output, added some polar stereographic support
 * 2/95  - v1.2.8: max min format changed %f to %g, tidying up more info
 * 3/95  - v1.3.0: fix bug with bitmap, allow numbers > UNDEFINED
 * 3/95  - v1.3.1: print number of missing points (verbose)
 * 3/95  - v1.3.2: -append option added
 * 4/95  - v1.3.2a,b: more output, polar stereo support (-V option)
 * 4/95  - v1.3.3: added ECMWF parameter table (prelim)
 * 6/95  - v1.3.4: nxny from BDS rather than gds?
 * 9/95  - v1.3.4d: speedup in grib write
 * 11/95 - v1.3.4f: new ECMWF parameter table (from Mike Fiorino), EC logic
 * 2/96  - v1.3.4g-h: prelim fix for GDS-less grib files
 * 2/96  - v1.3.4i: faster missing(), -V: "pos n" -> "n" (field 2)
 * 3/96  - v1.4: fix return code (!inventory), and short records near EOF
 * 6/96  - v1.4.1a: faster grib->binary decode, updated ncep parameter table, mod. in clim. desc
 * 7/96  - v1.5.0: parameter-table aware, -v option changed, added "comments"
 *                 increased NTRY to 100 in seek_grib
 * 11/96 - v1.5.0b: added ECMWF parameter table 128
 * 1/97 - v1.5.0b2: if nxny != nx*ny { nx = nxny; ny = 1 }
 * 3/97 - v1.5.0b5: added: -PDS -GDS, Lambert Conformal
 * 3/97 - v1.5.0b6: added: -verf
 * 4/97 - v1.5.0b7: added -PDS10, -GDS10 and enhanced -PDS -GDS
 * 4/97 - v1.5.0b8: "bitmap missing x" -> "bitmap: x undef"
 * 5/97 - v1.5.0b9: thinned grids meta data
 * 5/97 - v1.5.0b10: changed 0hr fcst to anal for TR=10 and P1=P2=0
 * 5/97 - v1.5.0b10: added -H option
 * 6/97 - v1.5.0b12: thinned lat-long grids -V option
 * 6/97 - v1.5.0b13: -4yr
 * 6/97 - v1.5.0b14: fix century mark Y=100 not 0
 * 7/97 - v1.5.0b15: add ncep opn grib table
 * 12/97 - v1.6.1.a: made ncep_opn the default table
 * 12/97 - v1.6.1.b: changed 03TOT to O3TOT in operational ncep table
 * 1/98  - v1.6.2: added Arakawa E grid meta-data
 * 1/98  - v1.6.2.1: added some mode data, Scan -> scan
 * 4/98  - v1.6.2.4: reanalysis id code: subcenter==0 && process==180
 * 5/98  - v1.6.2.5: fix -H code to write all of GDS
 * 7/98  - v1.7: fix decoding bug for bitmap and no. bits > 24 (theoretical bug)
 * 7/98  - v1.7.0.b1: add km to Mercator meta-data
 *
 *


DESCRIPTION


!END*************************************************************************/

/* NODAAC may need to be defined when running this code
in any non-DAAC processing mode.  This feature was added
by MODIS SDST to ensure compliance with ESDIS standards upon
delivering this software to the DAAC. */

#ifdef NODAAC

#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#include "grib.h"

int read_grib_array(FILE *input, char *what, char *where, int *nx, int *ny, float **narray);
int read_grib_date(FILE *input, char *what, char *where, char *date);
unsigned char *seek_grib(FILE *file, long *pos, long *len_grib, 
                         unsigned char *buffer, unsigned int buf_len);
int read_grib(FILE *file, long pos, long len_grib, unsigned char *buffer);
double ibm2flt(unsigned char *ibm);
void BDS_unpack(float *flt, unsigned char *bits, unsigned char *bitmap,
                int n_bits, int n, double ref, double scale);
double int_power(double x, int y);
char *levels(int, int);
int missing_points(unsigned char *bitmap, int n);
int GDS_grid(unsigned char *gds, int *nx, int *ny, long int *nxny);
int ASCII_TCA_PDS_date(unsigned char *pds, int v_time, char date[]);
int add_time(int *year, int *month, int *day, int *hour, int dtime, int unit);
int verf_time(unsigned char *pds, int *year, int *month, int *day, int *hour);


int read_grib_array(FILE *input, char *what, char *where, 
                    int *n_rows, int *n_cols, float **narray) 
/*
!C*****************************************************************************
!Description:  Routine reads data from a GRIB-format file.  The data to be read
               has the character strings 'what' and 'where' in its internal 
	       GRIB-format label.

!Input Parameters:
input:    pointer to GRIB-format file, returned from fopen()
what:     pointer to character string of what is to be looked for
where:    pointer to character string of the height; e. g., "atmos col" for whole-column 
          data, "850 mb" for data at 850 mb., etc.

          Stuff typically used:
	  
	     variable               what     where                       units
	     ~~~~~~~~               ~~~~     ~~~~~                       ~~~~~

            precipitable water     "PWAT"    "atmos col" (whole-column)  kg/m^2
	    pressure               "PRES"    "sfc"       (surface)       Pa
	    ozone                  "TOZNE"   "atmos col" (whole-column)  Dobsons

          Note: the exact syntax used for 'what' and 'where' depends upon the source
	  of the GRIB file.  See the parameter tables below (lines 1453 - 2771) for 
	  the syntax for 'what', and routine levels() (starting at line 1073) for the 
	  syntax for 'where'.
 

!Output Parameters:
n_rows:   the number of rows and
n_cols:       columns of data retrieved
narray:   array that contains the data

!Input/Output Parameters:
narray -- should be null when this routine is called, should be non-null and filled 
          when routine returns.

!Returns:
	-1	memory error
	-2      cannot find correct data; equivalent to 
	        "missing GRIB record(s)" in original
	-3      file error -- "missing end section" 
	 0	success

!Revision History:
        Revision 1.1, 27-MAY-99, Jim Ray (SSAI)
	1.  Edited "main()" of file wgrib.c into this routine.

	Original version:
	  v1.7.0b1 (8-24-98) Wesley Ebisuzaki (NCEP/NCAR Reanalysis Project).

!Team-unique Header:

   This software developed by the MODIS Land Science Team for 
   the National Aeronautics and Space Administration, Goddard 
   Space Flight Center, under contract NAS5-96062

!References and Credits

		Jim Ray
		Science Systems and Applications Inc.
		NASA's Goddard Space Flight Center Code 923
		Greenbelt MD, 20771
 		jim@cratmos.gsfc.nasa.gov, (301) 614-6684

!Design Notes
!END***************************************************************************
*/
{

    unsigned char *buffer;
    double temp;
    int i, nx, ny;
    long int len_grib, pos = 0, nxny, buffer_size, count = 1;
    unsigned char *msg, *pds, *gds, *bms, *bds, *pointer;
    char level[130];
    float *temp_array;
    int j, jj;
    long lj, ljj;

/* Lots removed!  This was the "main()" routine from file wgrib.c. */
/* Error codes:  -2, file error ("what" is not found) "missing GRIB record(s)"
                 -3, file error "missing end section"
                 -1, memory error */

    if ((buffer = (unsigned char *) malloc(BUFF_ALLOC0)) == NULL) {
        return(-1);
    }
    buffer_size = BUFF_ALLOC0;

    for (;;) {
	msg = seek_grib(input, &pos, &len_grib, buffer, MSEEK);
	if (msg == NULL) {
	    return(-2);
	}

        /* read all whole grib record */
        if (len_grib + msg - buffer > buffer_size) {
            buffer_size = len_grib + msg - buffer + 1000;
            buffer = (unsigned char *) realloc((void *) buffer, buffer_size);
            if (buffer == NULL) {
                return(-1);
            }
        }
        read_grib(input, pos, len_grib, buffer);

	/* parse grib message */

	msg = buffer;
        pds = (msg + 8);
        pointer = pds + PDS_LEN(pds);

        if (PDS_HAS_GDS(pds)) {
            gds = pointer;
            pointer += GDS_LEN(gds);
        }
        else {
            gds = NULL;
        }

        if (PDS_HAS_BMS(pds)) {
            bms = pointer;
            pointer += BMS_LEN(bms);
        }
        else {
            bms = NULL;
        }

        bds = pointer;
        pointer += BDS_LEN(bds);

        /* end section - "7777" in ascii */
        if (pointer[0] != 0x37 || pointer[1] != 0x37 ||
            pointer[2] != 0x37 || pointer[3] != 0x37) {
            return(-3);
        }

	/* figure out size of array */
	if (gds != NULL) {
	    /* this doesn't work for spherical harmonics */
	    GDS_grid(gds, &nx, &ny, &nxny);
	}
	else if (bms != NULL) {
	    nxny = nx = BMS_nxny(bms);
	    ny = 1;
	}
	else {
	    if (BDS_NumBits(bds) == 0) {
                nxny = nx = 1;
                /*fprintf(stderr,"Missing GDS, constant record .. cannot "
                    "determine number of data points\n");*/
	    }
	    else {
	        nxny = nx = BDS_NValues(bds);
	    }
	    ny = 1;
	}

#ifdef CHECK_GRIB
        if (BDS_NumBits(bds) != 0) {
	    i = BDS_NValues(bds);
	    if (bms != NULL) {
		i += missing_points(BMS_bitmap(bms),nx*ny);
	    }
	    if (i != nxny) {
	        /*fprintf(stderr,"grib header at record %ld: two values of nxny %ld %d\n",
			count,nxny,i);
	        fprintf(stderr,"   LEN %d DataStart %d UnusedBits %d #Bits %d nx %d ny %d\n",
		    BDS_LEN(bds), BDS_DataStart(bds),BDS_UnusedBits(bds),
		    BDS_NumBits(bds), nx, ny);*/
		nxny = nx = i;
		ny = 1;
	    }
	}
#endif
 

        strcpy(level, levels(PDS_KPDS6(pds), PDS_KPDS7(pds)));
printf ("GRIB.C: levels is %s\n", level);
        if ((!strcmp(k5toa(pds),what))&&
	    (!strcmp(level,where))) {
            if ((*narray = (float *) malloc(sizeof(float) * nxny)) == NULL) {
                return(-1);
            }
            if ((temp_array = (float *) malloc(sizeof(float) * nxny)) == NULL) {
                return(-1);
            }

	    temp = int_power(10.0, - PDS_DecimalScale(pds));

            BDS_unpack(temp_array, bds + 11, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
	        temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));

            free(buffer);
	    
            /* 15-JUN-99: Reformat the data! Each line now starts at longitude = 0; 
	       reformat to start at longitude = -180 (or so). */
            lj = 0L;
            for (i=0;i<ny;i++)
              for (j=0;j<nx;j++) {
                 jj = j + nx/2;
                 if (jj >= nx) jj -= nx;
                 ljj = i*nx + jj;
                 (*narray)[ljj] = temp_array[lj];
                 lj++;
                                 }
	    
	    free(temp_array);
	    *n_cols = nx;
	    *n_rows = ny;
	    return(0);
	  }
	  
        pos += len_grib;
        count++;
    }
}


int read_grib_date(FILE *input, char *what, char *where, char *date) 
/*
!C*****************************************************************************
!Description:  Routine reads date of data in a GRIB-format file.  The data whose
               date is to be read has the character strings 'what' and 'where' 
	       in its internal GRIB-format label.

!Input Parameters:
input:    pointer to GRIB-format file, returned from fopen()
what:     pointer to character string of what is to be looked for
where:    pointer to character string of the height; e. g., "atmos col" for whole-column 
          data, "850 mb" for data at 850 mb., etc.

          For more information about 'what' amd 'where' syntax, see notes for routine
	  read_grib_array(), above.
 

!Output Parameters:
date      the date of data at 'what' and 'where'.  'date' will be unchanged if data at
          'what' and 'where' is not found; 'date' is in ASCII time code A (28 characters).
	  
!Input/Output Parameters:
none

!Returns:
	-1	memory error
	-2      cannot find correct data; equivalent to 
	        "missing GRIB record(s)" in original
	-3      file error -- "missing end section" 
	 0	success

!Revision History:
        Revision 1.1, 15-JUN-99, Jim Ray (SSAI)
	1.  Edited "read_grib_array()" of file grib.c into this routine.

	Original version:
	  v1.7.0b1 (8-24-98) Wesley Ebisuzaki (NCEP/NCAR Reanalysis Project).

!Team-unique Header:

   This software developed by the MODIS Land Science Team for 
   the National Aeronautics and Space Administration, Goddard 
   Space Flight Center, under contract NAS5-96062

!References and Credits

		Jim Ray
		Science Systems and Applications Inc.
		NASA's Goddard Space Flight Center Code 923
		Greenbelt MD, 20771
 		jim@cratmos.gsfc.nasa.gov, (301) 614-6684

!Design Notes
!END***************************************************************************
*/
{

    unsigned char *buffer;
    int i, nx, ny;
    long int len_grib, pos = 0, nxny, buffer_size, count = 1;
    unsigned char *msg, *pds, *gds, *bms, *bds, *pointer;
    int v_time = 0;
    char level[130];

/* Lots removed!  This was the "main()" routine from file wgrib.c. */
/* Error codes:  -2, file error ("what" is not found) "missing GRIB record(s)"
                 -3, file error "missing end section"
                 -1, memory error */

    if ((buffer = (unsigned char *) malloc(BUFF_ALLOC0)) == NULL) {
        return(-1);
    }
    buffer_size = BUFF_ALLOC0;

    for (;;) {
	msg = seek_grib(input, &pos, &len_grib, buffer, MSEEK);
	if (msg == NULL) {
            free(buffer);
	    return(-2);
	}

        /* read all whole grib record */
        if (len_grib + msg - buffer > buffer_size) {
            buffer_size = len_grib + msg - buffer + 1000;
            buffer = (unsigned char *) realloc((void *) buffer, buffer_size);
            if (buffer == NULL) {
                return(-1);
            }
        }
        read_grib(input, pos, len_grib, buffer);

	/* parse grib message */

	msg = buffer;
        pds = (msg + 8);
        pointer = pds + PDS_LEN(pds);

        if (PDS_HAS_GDS(pds)) {
            gds = pointer;
            pointer += GDS_LEN(gds);
        }
        else {
            gds = NULL;
        }

        if (PDS_HAS_BMS(pds)) {
            bms = pointer;
            pointer += BMS_LEN(bms);
        }
        else {
            bms = NULL;
        }

        bds = pointer;
        pointer += BDS_LEN(bds);

        /* end section - "7777" in ascii */
        if (pointer[0] != 0x37 || pointer[1] != 0x37 ||
            pointer[2] != 0x37 || pointer[3] != 0x37) {
	    free(buffer);
            return(-3);
        }

	/* figure out size of array */
	if (gds != NULL) {
	    /* this doesn't work for spherical harmonics */
	    GDS_grid(gds, &nx, &ny, &nxny);
	}
	else if (bms != NULL) {
	    nxny = nx = BMS_nxny(bms);
	    ny = 1;
	}
	else {
	    if (BDS_NumBits(bds) == 0) {
                nxny = nx = 1;
                /*fprintf(stderr,"Missing GDS, constant record .. cannot "
                    "determine number of data points\n");*/
	    }
	    else {
	        nxny = nx = BDS_NValues(bds);
	    }
	    ny = 1;
	}

#ifdef CHECK_GRIB
        if (BDS_NumBits(bds) != 0) {
	    i = BDS_NValues(bds);
	    if (bms != NULL) {
		i += missing_points(BMS_bitmap(bms),nx*ny);
	    }
	    if (i != nxny) {
	        /*fprintf(stderr,"grib header at record %ld: two values of nxny %ld %d\n",
			count,nxny,i);
	        fprintf(stderr,"   LEN %d DataStart %d UnusedBits %d #Bits %d nx %d ny %d\n",
		    BDS_LEN(bds), BDS_DataStart(bds),BDS_UnusedBits(bds),
		    BDS_NumBits(bds), nx, ny);*/
		nxny = nx = i;
		ny = 1;
	    }
	}
#endif
 

        strcpy(level, levels(PDS_KPDS6(pds), PDS_KPDS7(pds)));
printf ("GRIB.C: levels is %s\n", level);
        if ((!strcmp(k5toa(pds),what))&&
	    (!strcmp(level,where))) {

	    ASCII_TCA_PDS_date(pds, v_time, date);
	    free(buffer);
	    return(0);
	  }
	  
        pos += len_grib;
        count++;
    }
}

#ifndef min
   #define min(a,b)  ((a) < (b) ? (a) : (b))
#endif

#define NTRY 100
/* #define LEN_HEADER_PDS (28+42+100) */
#define LEN_HEADER_PDS (28+8)

unsigned char *seek_grib(FILE *file, long *pos, long *len_grib, 
                         unsigned char *buffer, unsigned int buf_len) 

/*
!C**********************************************************************
!Description:
  finds next grib header.
  
!Input Parameters:
  FILE *file            -- input grib file pointer
  long *pos             -- Position from file start to start reading from 
                           ( = 0 for 1st call).
  unsigned int buf_len  -- Size of buffer
  unsigned char *buffer -- Buffer used to read grib record into.
  
!Output Parameters:
  long *pos             -- Position of 'GRIB' in file, or pos + buf_len -
                             LEN_HEADER_PDS if 'GRIB' not found in NTRY tries.
  long *len_grib        -- Length of grib record as recorded after 'GRIB'.
  unsigned char *buffer -- Buffer used to read grib record into.
  
!Return Value:
  (char *) to start of GRIB header+PDS, NULL if not found


!Revision History:
 * adapted from SKGB (Mark Iredell)
 *
 * v1.1 9/94 Wesley Ebisuzaki
 * v1.2 3/96 Wesley Ebisuzaki handles short records at end of file
 * v1.3 8/96 Wesley Ebisuzaki increase NTRY from 3 to 100 for the folks
 *      at Automation decided a 21 byte WMO bulletin header wasn't long
 *      enough and decided to go to an 8K header.

  Derived from wgrib, a portable program to read grib files that were 
  created by the NCEP/NCAR Reanalysis Project.  

!Team-unique Header:

!References and credits:

!Design Notes:

!END*************************************************************
*/
{

    int i, j, len;

    for (j = 0; j < NTRY; j++) {

        if (fseek(file, *pos, SEEK_SET) == -1) {
            *len_grib = 0;
            return (unsigned char *) NULL;
        }
     
        i = fread(buffer, sizeof (unsigned char), buf_len, file);
     
        len = i - LEN_HEADER_PDS;
     
        for (i = 0; i < len; i++) {
            if (buffer[i] == 'G' && buffer[i+1] == 'R' && buffer[i+2] == 'I'
                && buffer[i+3] == 'B' && buffer[i+7] == 1) {
                    *pos = i + *pos;
                    *len_grib = (buffer[i+4] << 16) + (buffer[i+5] << 8) +
                            buffer[i+6];
                    return (buffer+i);
            }
        }
	*pos = *pos + (buf_len - LEN_HEADER_PDS);
    }

    *len_grib = 0;
    return (unsigned char *) NULL;
}


double ibm2flt(unsigned char *ibm) 
/*
!C**********************************************************************
!Description:
  gets floating point value from four-byte character number representation
  
!Input Parameters:
  unsigned char *ibm  -- the number
  
!Output Parameters:
  none
  
!Return Value:
  double-precision value of number, or 0.0 if error

!Revision History:
 Original version, Wesley Ebisuzaki;
 *
 * v1.1 .. faster
 * v1.1 .. if mant == 0 -> quick return
 *
  Derived from wgrib, a portable program to read grib files that were 
  created by the NCEP/NCAR Reanalysis Project.  

!Team-unique Header:

!References and credits:

!Design Notes:

!END*************************************************************
*/
{
	int positive, power;
	unsigned int abspower;
	long int mant;
	double value, exp;

	mant = (ibm[1] << 16) + (ibm[2] << 8) + ibm[3];
        if (mant == 0) return 0.0;

	positive = (ibm[0] & 0x80) == 0;
	power = (int) (ibm[0] & 0x7f) - 64;
	abspower = power > 0 ? power : -power;


	/* calc exp */
	exp = 16.0;
	value = 1.0;
	while (abspower) {
		if (abspower & 1) {
			value *= exp;
		}
		exp = exp * exp;
		abspower >>= 1;
	}

	if (power < 0) value = 1.0 / value;
	value = value * mant / 16777216.0;
	if (positive == 0) value = -value;
	return value;
}

int read_grib(FILE *file, long pos, long len_grib, unsigned char *buffer) 
/*
!C**********************************************************************
!Description:
  Read full GRIB record of length len_grib into buffer.

!Input Parameters:
  FILE          *file
  long          pos      -- File position from start to begin reading from.
  long          len_grib -- Number of chars to read from file.

!Output Parameters:
  unsigned char *buffer  -- Buffer to read grib record into

!Return Value:
  (int) (i == len_grib) -- if false, error has occured.

!Revision History:
 * 17-JUN-99 Jim Ray, SSAI -- Added this prolog.
 * v1.0 9/94 Wesley Ebisuzaki

!Team-unique Header:

!References and credits:
    Wesley Ebisuzaki

!Design Notes:

!END*************************************************************
*/
{
    int i;

    if (fseek(file, pos, SEEK_SET) == -1) {
	    return 0;
    }

    i = fread(buffer, sizeof (unsigned char), len_grib, file);
    return (i == len_grib);
}

double int_power(double x, int y) 
/*
!C**********************************************************************
!Description:
  returns x**y
  
!Input Parameters:
  double x  --  the number to raise to...
  int y     --    the power
  
!Output Parameters:
  none
  
!Return Value:
  double-precision value of x**y

!Revision History:
 Original version, Wesley Ebisuzaki;
 *
 *  return x**y
 *
 *
 *  input: double x
 *	   int y
 *
  Derived from wgrib, a portable program to read grib files that were 
  created by the NCEP/NCAR Reanalysis Project.  

!Team-unique Header:

!References and credits:

!Design Notes:

!END*************************************************************
*/
{

	double value;

	if (y < 0) {
		y = -y;
		x = 1.0 / x;
	}
	value = 1.0;

	while (y) {
		if (y & 1) {
			value *= x;
		}
		x = x * x;
		y >>= 1;
	}
	return value;
}

struct ParmTable parm_table_ncep_opn[256];
struct ParmTable parm_table_ncep_reanal[256];
struct ParmTable parm_table_omb[256];
struct ParmTable parm_table_ecmwf_128[256];
struct ParmTable parm_table_ecmwf_160[256];
struct ParmTable parm_table_user[256];
enum Def_NCEP_Table def_ncep_table = DEF_T62_NCEP_TABLE;

static struct ParmTable *Parm_Table(unsigned char *pds) 
/*
!C**********************************************************************
!Description:
  returns pointer to the parameter table
  
!Input Parameters:
  unsigned char *pds  -- 
  
!Output Parameters:
  none
  
!Return Value:
  a pointer to the parameter table

!Revision History:
  Original version, Wesley Ebisuzaki;

  Derived from wgrib, a portable program to read grib files that were 
  created by the NCEP/NCAR Reanalysis Project.  

!Team-unique Header:

!References and credits:

!Design Notes:

!END*************************************************************
*/
{

    int i, center, subcenter, ptable, process;
    static int missing_count = 0, reanal_opn_count = 0;

    center = PDS_Center(pds);
    subcenter = PDS_Subcenter(pds);
    ptable = PDS_Vsn(pds);

#ifdef P_TABLE_FIRST
    i = setup_user_table(center, subcenter, ptable);
    if (i == 1) return &parm_table_user[0];
#endif
    /* NCEP opn and reanal */
    if (center == NMC && ptable <= 3) {
	if (subcenter == 1) return &parm_table_ncep_reanal[0];
        process = PDS_Model(pds);
	if (subcenter != 0 || (process != 80 && process != 180) || 
		(ptable != 1 && ptable != 2)) 
            return &parm_table_ncep_opn[0];

	/* at this point could be either the opn or reanalysis table */
	if (def_ncep_table == opn_nowarn) return &parm_table_ncep_opn[0];
	if (def_ncep_table == rean_nowarn) return &parm_table_ncep_reanal[0];
        if (reanal_opn_count++ == 0) {
	    /*fprintf(stderr, "Using NCEP %s table, see -ncep_opn, -ncep_rean options\n",
               (def_ncep_table == opn) ?  "opn" : "reanalysis");*/
	}
        return (def_ncep_table == opn) ?  &parm_table_ncep_opn[0] 
		: &parm_table_ncep_reanal[0];
    }
    if (center == NMC   && ptable == 128) return &parm_table_omb[0];
    if (center == ECMWF && ptable == 128) return &parm_table_ecmwf_128[0];
    if (center == ECMWF && ptable == 160) return &parm_table_ecmwf_160[0];


#ifndef P_TABLE_FIRST
    i = setup_user_table(center, subcenter, ptable);
    if (i == 1) return &parm_table_user[0];
#endif

    if ((ptable > 3 || (PDS_PARAM(pds)) > 127) && missing_count++ == 0) {
	/*fprintf(stderr,
            "\nUndefined parameter table (center %d-%d table %d), using NCEP-opn\n",
            center, subcenter, ptable);*/
    }
    return &parm_table_ncep_opn[0];
}

char *k5toa(unsigned char *pds) 
/*
!C**********************************************************************
!Description:
  returns the field name of the parameter table at pds
  
!Input Parameters:
  unsigned char *pds  -- 
  
!Output Parameters:
  none
  
!Return Value:
  field name of the parameter table at pds

!Revision History:
  Original version, Wesley Ebisuzaki;
  Derived from wgrib, a portable program to read grib files that were 
  created by the NCEP/NCAR Reanalysis Project.  

!Team-unique Header:

!References and credits:

!Design Notes:

!END*************************************************************
*/
{

    return (Parm_Table(pds) + PDS_PARAM(pds))->name;
}


char *k5_comments(unsigned char *pds) 
/*
!C**********************************************************************
!Description:
  returns the comment field of the parameter table at pds
  
!Input Parameters:
  unsigned char *pds  -- 
  
!Output Parameters:
  none
  
!Return Value:
  comment field of the parameter table at pds

!Revision History:
  Original version, Wesley Ebisuzaki;
  Derived from wgrib, a portable program to read grib files that were 
  created by the NCEP/NCAR Reanalysis Project.  

!Team-unique Header:

!References and credits:

!Design Notes:

!END*************************************************************
*/
{

    return (Parm_Table(pds) + PDS_PARAM(pds))->comment;
}

static unsigned int mask[] = {0,1,3,7,15,31,63,127,255};
static unsigned int map_masks[8] = {128, 64, 32, 16, 8, 4, 2, 1};

void BDS_unpack(float *flt, unsigned char *bits, unsigned char *bitmap,
	int n_bits, int n, double ref, double scale) 
/*
!C**********************************************************************
!Description:
  Converts packed grid-data in grib BinaryDataSection to array of float.

!Input Parameters:
  unsigned char *bits   -- BDS octet 12.
  unsigned char *bitmap -- Bit Map (BMS octet 7), or NULL if BMS not defined.
  int           n_bits  -- Nbits/datum, BDS octet 11.
  int           n       -- Number of data to write to flt.
  double        ref     -- (10^-D)*R
  double        scale   -- (10^-D)*2^E
    R = reference value     , BDS octet 7-10
    E = binary scale factor , BDS octet 5-6
    D = decimal scale factor, PDS octet 27-28
    All octet numbers above are 1 based, not 0

!Output Parameters:
  float         *flt    -- decoded grib data.

!Return Value: N/A

!Revision History:
  1996     wesley ebisuzaki -- Original
  1996/04 -- v1.1 faster
  17-JUN-99 Jim Ray, SSAI -- Added this prolog.

!Team-unique Header:

!References and credits:

!Design Notes:

!END*************************************************************
*/
{

    int i, mask_idx, t_bits, c_bits, j_bits;
    unsigned int j, map_mask, tbits, jmask, bbits;
    long int jj;

    tbits = bbits = 0;

    /* assume integer has 32+ bits */
    if (n_bits <= 25) {
        jmask = (1 << n_bits) - 1;
        t_bits = 0;

        if (bitmap) {
	    for (i = 0; i < n; i++) {
		/* check bitmap */
		mask_idx = i & 7;
		if (mask_idx == 0) bbits = *bitmap++;
	        if ((bbits & map_masks[mask_idx]) == 0) {
		    *flt++ = UNDEFINED;
		    continue;
	        }

	        while (t_bits < n_bits) {
	            tbits = (tbits * 256) + *bits++;
	            t_bits += 8;
	        }
	        t_bits -= n_bits;
	        j = (tbits >> t_bits) & jmask;
	        *flt++ = ref + scale*j;
            }
        }
        else {
	    for (i = 0; i < n; i++) {
                while (t_bits < n_bits) {
                    tbits = (tbits * 256) + *bits++;
                    t_bits += 8;
                }
                t_bits -= n_bits;
                flt[i] = (tbits >> t_bits) & jmask;
            }
	    /* at least this vectorizes :) */
	    for (i = 0; i < n; i++) {
		flt[i] = ref + scale*flt[i];
	    }
        }
    }
    else {
	/* older unoptimized code, not often used */
        c_bits = 8;
        map_mask = 128;
        while (n-- > 0) {
	    if (bitmap) {
	        j = (*bitmap & map_mask);
	        if ((map_mask >>= 1) == 0) {
		    map_mask = 128;
		    bitmap++;
	        }
	        if (j == 0) {
		    *flt++ = UNDEFINED;
		    continue;
	        }
	    }

	    jj = 0;
	    j_bits = n_bits;
	    while (c_bits <= j_bits) {
	        if (c_bits == 8) {
		    jj = (jj << 8) + *bits++;
		    j_bits -= 8;
	        }
	        else {
		    jj = (jj << c_bits) + (*bits & mask[c_bits]);
		    bits++;
		    j_bits -= c_bits;
		    c_bits = 8;
	        }
	    }
	    if (j_bits) {
	        c_bits -= j_bits;
	        jj = (jj << j_bits) + ((*bits >> c_bits) & mask[j_bits]);
	    }
	    *flt++ = ref + scale*jj;
        }
    }
    return;
}

char *levels(int kpds6, int kpds7) 
/*
!C**********************************************************************
!Description:
  prints out a simple description of kpds6, kpds7 (level/layer data)
  
!Input Parameters:
  int kpds6  -- octet 10 of the PDS
  int kpds7  -- octet 11 and 12 of the PDS
  
!Output Parameters:
  none
  
!Return Value:
  brief description of level/layer data

!Revision History:
 Original version, Wesley Ebisuzaki;
 *
 * levels.c
 *
 * prints out a simple description of kpds6, kpds7
 *    (level/layer data)
 *  kpds6 = octet 10 of the PDS
 *  kpds7 = octet 11 and 12 of the PDS
 *    (kpds values are from NMC's grib routines)
 *
 * the description of the levels is 
 *   (1) incomplete
 *   (2) include some NMC-only values (>= 200?)
 
  Derived from wgrib, a portable program to read grib files that were 
  created by the NCEP/NCAR Reanalysis Project.  

!Team-unique Header:

!References and credits:

!Design Notes:

!END*************************************************************
*/
{
    int o11, o12;
    static char x[128];
	
	/* octets 11 and 12 */
	o11 = kpds7 / 256;
	o12 = kpds7 % 256;


	switch (kpds6) {

	case 1: sprintf(x, "sfc");
		break;
	case 2: sprintf(x, "cld base");
		break;
	case 3: sprintf(x, "cld top");
		break;
	case 4: sprintf(x, "0C isotherm");
		break;
	case 5: sprintf(x, "cond lev");
		break;
	case 6: sprintf(x, "max wind lev");
		break;
	case 7: sprintf(x, "tropopause");
		break;
	case 8: sprintf(x, "nom. top");
		break;
	case 9: sprintf(x, "sea bottom");
		break;
	case 200:
	case 10: sprintf(x, "atmos col");
		break;

	case 12:
	case 212: sprintf(x, "low cld bot");
		break;
	case 13:
	case 213: sprintf(x, "low cld top");
		break;
	case 14:
	case 214: sprintf(x, "low cld lay");
		break;
	case 22:
	case 222: sprintf(x, "mid cld bot");
		break;
	case 23:
	case 223: sprintf(x, "mid cld top");
		break;
	case 24:
	case 224: sprintf(x, "mid cld lay");
		break;
	case 32:
	case 232: sprintf(x, "high cld bot");
		break;
	case 33:
	case 233: sprintf(x, "high cld top");
		break;
	case 34:
	case 234: sprintf(x, "high cld lay");
		break;


	case 100: sprintf(x, "%d mb",kpds7);
	 	break;
	case 101: sprintf(x, "%d-%d mb",o11*10,o12*10);
	 	break;
	case 102: sprintf(x, "MSL");
	 	break;
	case 103: sprintf(x, "%d m above MSL",kpds7);
	 	break;
	case 104: sprintf(x, "%d-%d m above msl",o11*100,o12*100);
	 	break;
	case 105: sprintf(x, "%d m above gnd",kpds7);
	 	break;
	case 106: sprintf(x, "%d-%d m above gnd",o11*100,o12*100);
	 	break;
	case 107: sprintf(x, "sigma=%.4f",kpds7/10000.0);
	 	break;
	case 108: sprintf(x, "sigma %.2f-%.2f",o11/100.0,o12/100.0);
	 	break;
	case 109: sprintf(x, "hybrid lev %d",kpds7);
	 	break;
	case 110: sprintf(x, "hybrid %d-%d",o11,o12);
	 	break;
	case 111: sprintf(x, "%d cm down",kpds7);
	 	break;
	case 112: sprintf(x, "%d-%d cm down",o11,o12);
	 	break;
	case 113: sprintf(x, "%dK",kpds7);
	 	break;
	case 114: sprintf(x, "%d-%dK",475-o11,475-o12);
	 	break;
	case 115: sprintf(x, "%d mb above gnd",kpds7);
	 	break;
	case 116: sprintf(x, "%d-%d mb above gnd",o11,o12);
	 	break;
	case 121: sprintf(x, "%d-%d mb",1100-o11,1100-o12);
	 	break;
	default:
	 	break;
	}
printf ("GRIB.C: levels = %s\n", x);

    return(x);
}

/*
 *  number of missing data points w. ebisuzaki
 *
 *  v1.1: just faster my dear
 *  v1.2: just faster my dear
 *
 */

static int bitsum[256] = {
    8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 4, 
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 
    5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 
    5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 
    5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 
    5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 
    4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0};


int missing_points(unsigned char *bitmap, int n) 
/*
!C**********************************************************************
!Description:
  
  
!Input Parameters:
  unsigned char *bitmap -- 
  int n                 -- 
  
!Output Parameters:
  none
  
!Return Value:
  

!Revision History:
 Original version, Wesley Ebisuzaki;
 *  number of missing data points w. ebisuzaki
 *
 *  v1.1: just faster my dear
 *  v1.2: just faster my dear
 * 
  Derived from wgrib, a portable program to read grib files that were 
  created by the NCEP/NCAR Reanalysis Project.  

!Team-unique Header:

!References and credits:

!Design Notes:

!END*************************************************************
*/
{

    int count;
    unsigned int tmp;
    if (bitmap == NULL) return 0;

    count = 0;
    while (n >= 8) {
	tmp = *bitmap++;
	n -= 8;
        count += bitsum[tmp];
    }
    tmp = *bitmap | ((1 << (8 - n)) - 1);
    count += bitsum[tmp];

    return count;
}

/*
 * parameter table for NCEP (operations)
 * center = 7, subcenter != 2 parameter table = 1, 2, 3 etc
 * note: see reanalysis parameter table for problems
 */

struct ParmTable parm_table_ncep_opn[256] = {
   /* 0 */ {"var0", "undefined"},
   /* 1 */ {"PRES", "Pressure [Pa]"},
   /* 2 */ {"PRMSL", "Pressure reduced to MSL [Pa]"},
   /* 3 */ {"PTEND", "Pressure tendency [Pa/s]"},
   /* 4 */ {"var4", "undefined"},
   /* 5 */ {"ICAHT", "ICAO Standard Atmosphere Reference Height [M]"},
   /* 6 */ {"GP", "Geopotential [m^2/s^2]"},
   /* 7 */ {"HGT", "Geopotential height [gpm]"},
   /* 8 */ {"DIST", "Geometric height [m]"},
   /* 9 */ {"HSTDV", "Std dev of height [m]"},
   /* 10 */ {"TOZNE", "Total ozone [Dobson]"},
   /* 11 */ {"TMP", "Temp. [K]"},
   /* 12 */ {"VTMP", "Virtual temp. [K]"},
   /* 13 */ {"POT", "Potential temp. [K]"},
   /* 14 */ {"EPOT", "Pseudo-adiabatic pot. temp. [K]"},
   /* 15 */ {"TMAX", "Max. temp. [K]"},
   /* 16 */ {"TMIN", "Min. temp. [K]"},
   /* 17 */ {"DPT", "Dew point temp. [K]"},
   /* 18 */ {"DEPR", "Dew point depression [K]"},
   /* 19 */ {"LAPR", "Lapse rate [K/m]"},
   /* 20 */ {"VIS", "Visibility [m]"},
   /* 21 */ {"RDSP1", "Radar spectra (1) [non-dim]"},
   /* 22 */ {"RDSP2", "Radar spectra (2) [non-dim]"},
   /* 23 */ {"RDSP3", "Radar spectra (3) [non-dim]"},
   /* 24 */ {"PLI", "Parcel lifted index (to 500 hPa) [K]"},
   /* 25 */ {"TMPA", "Temp. anomaly [K]"},
   /* 26 */ {"PRESA", "Pressure anomaly [Pa]"},
   /* 27 */ {"GPA", "Geopotential height anomaly [gpm]"},
   /* 28 */ {"WVSP1", "Wave spectra (1) [non-dim]"},
   /* 29 */ {"WVSP2", "Wave spectra (2) [non-dim]"},
   /* 30 */ {"WVSP3", "Wave spectra (3) [non-dim]"},
   /* 31 */ {"WDIR", "Wind direction [deg]"},
   /* 32 */ {"WIND", "Wind speed [m/s]"},
   /* 33 */ {"UGRD", "u wind [m/s]"},
   /* 34 */ {"VGRD", "v wind [m/s]"},
   /* 35 */ {"STRM", "Stream function [m^2/s]"},
   /* 36 */ {"VPOT", "Velocity potential [m^2/s]"},
   /* 37 */ {"MNTSF", "Montgomery stream function [m^2/s^2]"},
   /* 38 */ {"SGCVV", "Sigma coord. vertical velocity [/s]"},
   /* 39 */ {"VVEL", "Pressure vertical velocity [Pa/s]"},
   /* 40 */ {"DZDT", "Geometric vertical velocity [m/s]"},
   /* 41 */ {"ABSV", "Absolute vorticity [/s]"},
   /* 42 */ {"ABSD", "Absolute divergence [/s]"},
   /* 43 */ {"RELV", "Relative vorticity [/s]"},
   /* 44 */ {"RELD", "Relative divergence [/s]"},
   /* 45 */ {"VUCSH", "Vertical u shear [/s]"},
   /* 46 */ {"VVCSH", "Vertical v shear [/s]"},
   /* 47 */ {"DIRC", "Direction of current [deg]"},
   /* 48 */ {"SPC", "Speed of current [m/s]"},
   /* 49 */ {"UOGRD", "u of current [m/s]"},
   /* 50 */ {"VOGRD", "v of current [m/s]"},
   /* 51 */ {"SPFH", "Specific humidity [kg/kg]"},
   /* 52 */ {"RH", "Relative humidity [%]"},
   /* 53 */ {"MIXR", "Humidity mixing ratio [kg/kg]"},
   /* 54 */ {"PWAT", "Precipitable water [kg/m^2]"},
   /* 55 */ {"VAPP", "Vapor pressure [Pa]"},
   /* 56 */ {"SATD", "Saturation deficit [Pa]"},
   /* 57 */ {"EVP", "Evaporation [kg/m^2]"},
   /* 58 */ {"CICE", "Cloud Ice [kg/m^2]"},
   /* 59 */ {"PRATE", "Precipitation rate [kg/m^2/s]"},
   /* 60 */ {"TSTM", "Thunderstorm probability [%]"},
   /* 61 */ {"APCP", "Total precipitation [kg/m^2]"},
   /* 62 */ {"NCPCP", "Large scale precipitation [kg/m^2]"},
   /* 63 */ {"ACPCP", "Convective precipitation [kg/m^2]"},
   /* 64 */ {"SRWEQ", "Snowfall rate water equiv. [kg/m^2/s]"},
   /* 65 */ {"WEASD", "Accum. snow [kg/m^2]"},
   /* 66 */ {"SNOD", "Snow depth [m]"},
   /* 67 */ {"MIXHT", "Mixed layer depth [m]"},
   /* 68 */ {"TTHDP", "Transient thermocline depth [m]"},
   /* 69 */ {"MTHD", "Main thermocline depth [m]"},
   /* 70 */ {"MTHA", "Main thermocline anomaly [m]"},
   /* 71 */ {"TCDC", "Total cloud cover [%]"},
   /* 72 */ {"CDCON", "Convective cloud cover [%]"},
   /* 73 */ {"LCDC", "Low level cloud cover [%]"},
   /* 74 */ {"MCDC", "Mid level cloud cover [%]"},
   /* 75 */ {"HCDC", "High level cloud cover [%]"},
   /* 76 */ {"CWAT", "Cloud water [kg/m^2]"},
   /* 77 */ {"BLI", "Best lifted index (to 500 hPa) [K]"},
   /* 78 */ {"SNOC", "Convective snow [kg/m^2]"},
   /* 79 */ {"SNOL", "Large scale snow [kg/m^2]"},
   /* 80 */ {"WTMP", "Water temp. [K]"},
   /* 81 */ {"LAND", "Land-sea mask (land=1;sea=0) [fraction]"},
   /* 82 */ {"DSLM", "Deviation of sea level from mean [m]"},
   /* 83 */ {"SFCR", "Surface roughness [m]"},
   /* 84 */ {"ALBDO", "Albedo [%]"},
   /* 85 */ {"TSOIL", "Soil temp. [K]"},
   /* 86 */ {"SOILM", "Soil moisture content [kg/m^2]"},
   /* 87 */ {"VEG", "Vegetation [%]"},
   /* 88 */ {"SALTY", "Salinity [kg/kg]"},
   /* 89 */ {"DEN", "Density [kg/m^3]"},
   /* 90 */ {"WATR", "Water runoff [kg/m^2]"},
   /* 91 */ {"ICEC", "Ice concentration (ice=1;no ice=0) [fraction]"},
   /* 92 */ {"ICETK", "Ice thickness [m]"},
   /* 93 */ {"DICED", "Direction of ice drift [deg]"},
   /* 94 */ {"SICED", "Speed of ice drift [m/s]"},
   /* 95 */ {"UICE", "u of ice drift [m/s]"},
   /* 96 */ {"VICE", "v of ice drift [m/s]"},
   /* 97 */ {"ICEG", "Ice growth rate [m/s]"},
   /* 98 */ {"ICED", "Ice divergence [/s]"},
   /* 99 */ {"SNOM", "Snow melt [kg/m^2]"},
   /* 100 */ {"HTSGW", "Sig height of wind waves and swell [m]"},
   /* 101 */ {"WVDIR", "Direction of wind waves [deg]"},
   /* 102 */ {"WVHGT", "Sig height of wind waves [m]"},
   /* 103 */ {"WVPER", "Mean period of wind waves [s]"},
   /* 104 */ {"SWDIR", "Direction of swell waves [deg]"},
   /* 105 */ {"SWELL", "Sig height of swell waves [m]"},
   /* 106 */ {"SWPER", "Mean period of swell waves [s]"},
   /* 107 */ {"DIRPW", "Primary wave direction [deg]"},
   /* 108 */ {"PERPW", "Primary wave mean period [s]"},
   /* 109 */ {"DIRSW", "Secondary wave direction [deg]"},
   /* 110 */ {"PERSW", "Secondary wave mean period [s]"},
   /* 111 */ {"NSWRS", "Net short wave (surface) [W/m^2]"},
   /* 112 */ {"NLWRS", "Net long wave (surface) [W/m^2]"},
   /* 113 */ {"NSWRT", "Net short wave (top) [W/m^2]"},
   /* 114 */ {"NLWRT", "Net long wave (top) [W/m^2]"},
   /* 115 */ {"LWAVR", "Long wave [W/m^2]"},
   /* 116 */ {"SWAVR", "Short wave [W/m^2]"},
   /* 117 */ {"GRAD", "Global radiation [W/m^2]"},
   /* 118 */ {"BRTMP", "Brightness temperature [K]"},
   /* 119 */ {"LWRAD", "Long wave radiation [W/srm^2]"},
   /* 120 */ {"SWRAD", "Short wave radiation [W/srm2]"},
   /* 121 */ {"LHTFL", "Latent heat flux [W/m^2]"},
   /* 122 */ {"SHTFL", "Sensible heat flux [W/m^2]"},
   /* 123 */ {"BLYDP", "Boundary layer dissipation [W/m^2]"},
   /* 124 */ {"UFLX", "Zonal momentum flux [N/m^2]"},
   /* 125 */ {"VFLX", "Meridional momentum flux [N/m^2]"},
   /* 126 */ {"WMIXE", "Wind mixing energy [J]"},
   /* 127 */ {"IMGD", "Image data []"},
   /* 128 */ {"MSLSA", "Mean sea level pressure (Std Atm) [Pa]"},
   /* 129 */ {"MSLMA", "Mean sea level pressure (MAPS) [Pa]"},
   /* 130 */ {"MSLET", "Mean sea level pressure (ETA model) [Pa]"},
   /* 131 */ {"LFTX", "Surface lifted index [K]"},
   /* 132 */ {"4LFTX", "Best (4-layer) lifted index [K]"},
   /* 133 */ {"KX", "K index [K]"},
   /* 134 */ {"SX", "Sweat index [K]"},
   /* 135 */ {"MCONV", "Horizontal moisture divergence [kg/kg/s]"},
   /* 136 */ {"VWSH", "Vertical speed shear [1/s]"},
   /* 137 */ {"TSLSA", "3-hr pressure tendency (Std Atmos Red) [Pa/s]"},
   /* 138 */ {"BVF2", "Brunt-Vaisala frequency^2 [1/s^2]"},
   /* 139 */ {"PVMW", "Potential vorticity (mass-weighted) [1/s/m]"},
   /* 140 */ {"CRAIN", "Categorical rain [yes=1;no=0]"},
   /* 141 */ {"CFRZR", "Categorical freezing rain [yes=1;no=0]"},
   /* 142 */ {"CICEP", "Categorical ice pellets [yes=1;no=0]"},
   /* 143 */ {"CSNOW", "Categorical snow [yes=1;no=0]"},
   /* 144 */ {"SOILW", "Volumetric soil moisture [fraction]"},
   /* 145 */ {"PEVPR", "Potential evaporation rate [W/m^2]"},
   /* 146 */ {"CWORK", "Cloud work function [J/kg]"},
   /* 147 */ {"U-GWD", "Zonal gravity wave stress [N/m^2]"},
   /* 148 */ {"V-GWD", "Meridional gravity wave stress [N/m^2]"},
   /* 149 */ {"PV", "Potential vorticity [m^2/s/kg]"},
   /* 150 */ {"COVMZ", "Covariance between u and v [m^2/s^2]"},
   /* 151 */ {"COVTZ", "Covariance between u and T [K*m/s]"},
   /* 152 */ {"COVTM", "Covariance between v and T [K*m/s]"},
   /* 153 */ {"CLWMR", "Cloud water [kg/kg]"},
   /* 154 */ {"O3MR", "Ozone mixing ratio [kg/kg]"},
   /* 155 */ {"GFLUX", "Ground heat flux [W/m^2]"},
   /* 156 */ {"CIN", "Convective inhibition [J/kg]"},
   /* 157 */ {"CAPE", "Convective Avail. Pot. Energy [J/kg]"},
   /* 158 */ {"TKE", "Turbulent kinetic energy [J/kg]"},
   /* 159 */ {"CONDP", "Lifted parcel condensation pressure [Pa]"},
   /* 160 */ {"CSUSF", "Clear sky upward solar flux [W/m^2]"},
   /* 161 */ {"CSDSF", "Clear sky downward solar flux [W/m^2]"},
   /* 162 */ {"CSULF", "Clear sky upward long wave flux [W/m^2]"},
   /* 163 */ {"CSDLF", "Clear sky downward long wave flux [W/m^2]"},
   /* 164 */ {"CFNSF", "Cloud forcing net solar flux [W/m^2]"},
   /* 165 */ {"CFNLF", "Cloud forcing net long wave flux [W/m^2]"},
   /* 166 */ {"VBDSF", "Visible beam downward solar flux [W/m^2]"},
   /* 167 */ {"VDDSF", "Visible diffuse downward solar flux [W/m^2]"},
   /* 168 */ {"NBDSF", "Near IR beam downward solar flux [W/m^2]"},
   /* 169 */ {"NDDSF", "Near IR diffuse downward solar flux [W/m^2]"},
   /* 170 */ {"var170", "undefined"},
   /* 171 */ {"var171", "undefined"},
   /* 172 */ {"MFLX", "Momentum flux [N/m^2]"},
   /* 173 */ {"LMH", "Mass point model surface [non-dim]"},
   /* 174 */ {"LMV", "Velocity point model surface [non-dim]"},
   /* 175 */ {"MLYNO", "Model layer number (from bottom up) [non-dim]"},
   /* 176 */ {"NLAT", "Latitude (-90 to +90) [deg]"},
   /* 177 */ {"ELON", "East longitude (0-360) [deg]"},
   /* 178 */ {"var178", "undefined"},
   /* 179 */ {"var179", "undefined"},
   /* 180 */ {"var180", "undefined"},
   /* 181 */ {"LPSX", "x-gradient of log pressure [1/m]"},
   /* 182 */ {"LPSY", "y-gradient of log pressure [1/m]"},
   /* 183 */ {"HGTX", "x-gradient of height [m/m]"},
   /* 184 */ {"HGTY", "y-gradient of height [m/m]"},
   /* 185 */ {"var185", "undefined"},
   /* 186 */ {"var186", "undefined"},
   /* 187 */ {"var187", "undefined"},
   /* 188 */ {"var188", "undefined"},
   /* 189 */ {"VPTMP", "Virtual pot. temp. [K]"},
   /* 190 */ {"HLCY", "Storm relative helicity [m^2/s^2]"},
   /* 191 */ {"PROB", "Prob. from ensemble [non-dim]"},
   /* 192 */ {"PROBN", "Prob. from ensemble norm. to clim. expect. [non-dim]"},
   /* 193 */ {"POP", "Prob. of precipitation [%]"},
   /* 194 */ {"CPOFP", "Prob. of frozen precipitation [%]"},
   /* 195 */ {"CPOZP", "Prob. of freezing precipitation [%]"},
   /* 196 */ {"USTM", "u-component of storm motion [m/s]"},
   /* 197 */ {"VSTM", "v-component of storm motion [m/s]"},
   /* 198 */ {"var198", "undefined"},
   /* 199 */ {"var199", "undefined"},
   /* 200 */ {"var200", "undefined"},
   /* 201 */ {"ICWAT", "Ice-free water surface [%]"},
   /* 202 */ {"var202", "undefined"},
   /* 203 */ {"var203", "undefined"},
   /* 204 */ {"DSWRF", "Downward short wave flux [W/m^2]"},
   /* 205 */ {"DLWRF", "Downward long wave flux [W/m^2]"},
   /* 206 */ {"UVI", "Ultra violet index (1 hour centered at solar noon) [J/m^2]"},
   /* 207 */ {"MSTAV", "Moisture availability [%]"},
   /* 208 */ {"SFEXC", "Exchange coefficient [(kg/m^3)(m/s)]"},
   /* 209 */ {"MIXLY", "No. of mixed layers next to surface [integer]"},
   /* 210 */ {"var210", "undefined"},
   /* 211 */ {"USWRF", "Upward short wave flux [W/m^2]"},
   /* 212 */ {"ULWRF", "Upward long wave flux [W/m^2]"},
   /* 213 */ {"CDLYR", "Non-convective cloud [%]"},
   /* 214 */ {"CPRAT", "Convective precip. rate [kg/m^2/s]"},
   /* 215 */ {"TTDIA", "Temp. tendency by all physics [K/s]"},
   /* 216 */ {"TTRAD", "Temp. tendency by all radiation [K/s]"},
   /* 217 */ {"TTPHY", "Temp. tendency by non-radiation physics [K/s]"},
   /* 218 */ {"PREIX", "Precip index (0.0-1.00) [fraction]"},
   /* 219 */ {"TSD1D", "Std. dev. of IR T over 1x1 deg area [K]"},
   /* 220 */ {"NLGSP", "Natural log of surface pressure [ln(kPa)]"},
   /* 221 */ {"HPBL", "Planetary boundary layer height [m]"},
   /* 222 */ {"5WAVH", "5-wave geopotential height [gpm]"},
   /* 223 */ {"CNWAT", "Plant canopy surface water [kg/m^2]"},
   /* 224 */ {"var224", "undefined"},
   /* 225 */ {"var225", "undefined"},
   /* 226 */ {"BMIXL", "Blackadar's mixing length scale [m]"},
   /* 227 */ {"AMIXL", "Asymptotic mixing length scale [m]"},
   /* 228 */ {"PEVAP", "Pot. evaporation [kg/m^2]"},
   /* 229 */ {"SNOHF", "Snow phase-change heat flux [W/m^2]"},
   /* 230 */ {"var230", "undefined"},
   /* 231 */ {"MFLUX", "Convective cloud mass flux [Pa/s]"},
   /* 232 */ {"DTRF", "Downward total radiation flux [W/m^2]"},
   /* 233 */ {"UTRF", "Upward total radiation flux [W/m^2]"},
   /* 234 */ {"BGRUN", "Baseflow-groundwater runoff [kg/m^2]"},
   /* 235 */ {"SSRUN", "Storm surface runoff [kg/m^2]"},
   /* 236 */ {"var236", "undefined"},
   /* 237 */ {"O3TOT", "Total ozone [kg/m^2]"},
   /* 238 */ {"SNOWC", "Snow cover [%]"},
   /* 239 */ {"SNOT", "Snow temp. [K]"},
   /* 240 */ {"var240", "undefined"},
   /* 241 */ {"LRGHR", "Large scale condensation heating [K/s]"},
   /* 242 */ {"CNVHR", "Deep convective heating [K/s]"},
   /* 243 */ {"CNVMR", "Deep convective moistening [kg/kg/s]"},
   /* 244 */ {"SHAHR", "Shallow convective heating [K/s]"},
   /* 245 */ {"SHAMR", "Shallow convective moistening [kg/kg/s]"},
   /* 246 */ {"VDFHR", "Vertical diffusion heating [K/s]"},
   /* 247 */ {"VDFUA", "Vertical diffusion zonal accel [m/s^2]"},
   /* 248 */ {"VDFVA", "Vertical diffusion meridional accel [m/s^2]"},
   /* 249 */ {"VDFMR", "Vertical diffusion moistening [kg/kg/s]"},
   /* 250 */ {"SWHR", "Solar radiative heating [K/s]"},
   /* 251 */ {"LWHR", "Longwave radiative heating [K/s]"},
   /* 252 */ {"CD", "Drag coefficient [non-dim]"},
   /* 253 */ {"FRICV", "Friction velocity [m/s]"},
   /* 254 */ {"RI", "Richardson number [non-dim]"},
   /* 255 */ {"var255", "undefined"},
};

/*
 * parameter table for the NCEP/NCAR Reanalysis Project
 * center = 7, subcenter = 0/2, parameter table = 1/2
 * in a SNAFU the operational and reanalysis tables diverged
 * and both retained the same parameter table numbers (1,2)
 *
 * some of the Reanalysis files have subcenter=2 while others
 * use subcenter=0  (subcenter field is not standard (7/97))
 *
 * Some ways to tell Reanalysis files from OPN files
 *  Reanalysis: always generated by process 80 - T62 28 level model
 * Original subcenter=0 Reanalysis files had 
 *  2.5x2.5 (144x73) lat-long grid or 192x94 Gaussian grid (PDS grid=255?)
 */

struct ParmTable parm_table_ncep_reanal[256] = {
   /* 0 */ {"var0", "undefined"},
   /* 1 */ {"PRES", "Pressure [Pa]"},
   /* 2 */ {"PRMSL", "Pressure reduced to MSL [Pa]"},
   /* 3 */ {"PTEND", "Pressure tendency [Pa/s]"},
   /* 4 */ {"var4", "undefined"},
   /* 5 */ {"var5", "undefined"},
   /* 6 */ {"GP", "Geopotential [m^2/s^2]"},
   /* 7 */ {"HGT", "Geopotential height [gpm]"},
   /* 8 */ {"DIST", "Geometric height [m]"},
   /* 9 */ {"HSTDV", "Std dev of height [m]"},
   /* 10 */ {"HVAR", "Variance of height [m^2]"},
   /* 11 */ {"TMP", "Temp. [K]"},
   /* 12 */ {"VTMP", "Virtual temp. [K]"},
   /* 13 */ {"POT", "Potential temp. [K]"},
   /* 14 */ {"EPOT", "Pseudo-adiabatic pot. temp. [K]"},
   /* 15 */ {"TMAX", "Max. temp. [K]"},
   /* 16 */ {"TMIN", "Min. temp. [K]"},
   /* 17 */ {"DPT", "Dew point temp. [K]"},
   /* 18 */ {"DEPR", "Dew point depression [K]"},
   /* 19 */ {"LAPR", "Lapse rate [K/m]"},
   /* 20 */ {"VISIB", "Visibility [m]"},
   /* 21 */ {"RDSP1", "Radar spectra (1) [non-dim]"},
   /* 22 */ {"RDSP2", "Radar spectra (2) [non-dim]"},
   /* 23 */ {"RDSP3", "Radar spectra (3) [non-dim]"},
   /* 24 */ {"var24", "undefined"},
   /* 25 */ {"TMPA", "Temp. anomaly [K]"},
   /* 26 */ {"PRESA", "Pressure anomaly [Pa]"},
   /* 27 */ {"GPA", "Geopotential height anomaly [gpm]"},
   /* 28 */ {"WVSP1", "Wave spectra (1) [non-dim]"},
   /* 29 */ {"WVSP2", "Wave spectra (2) [non-dim]"},
   /* 30 */ {"WVSP3", "Wave spectra (3) [non-dim]"},
   /* 31 */ {"WDIR", "Wind direction [deg]"},
   /* 32 */ {"WIND", "Wind speed [m/s]"},
   /* 33 */ {"UGRD", "u wind [m/s]"},
   /* 34 */ {"VGRD", "v wind [m/s]"},
   /* 35 */ {"STRM", "Stream function [m^2/s]"},
   /* 36 */ {"VPOT", "Velocity potential [m^2/s]"},
   /* 37 */ {"MNTSF", "Montgomery stream function [m^2/s^2]"},
   /* 38 */ {"SGCVV", "Sigma coord. vertical velocity [/s]"},
   /* 39 */ {"VVEL", "Pressure vertical velocity [Pa/s]"},
   /* 40 */ {"DZDT", "Geometric vertical velocity [m/s]"},
   /* 41 */ {"ABSV", "Absolute vorticity [/s]"},
   /* 42 */ {"ABSD", "Absolute divergence [/s]"},
   /* 43 */ {"RELV", "Relative vorticity [/s]"},
   /* 44 */ {"RELD", "Relative divergence [/s]"},
   /* 45 */ {"VUCSH", "Vertical u shear [/s]"},
   /* 46 */ {"VVCSH", "Vertical v shear [/s]"},
   /* 47 */ {"DIRC", "Direction of current [deg]"},
   /* 48 */ {"SPC", "Speed of current [m/s]"},
   /* 49 */ {"UOGRD", "u of current [m/s]"},
   /* 50 */ {"VOGRD", "v of current [m/s]"},
   /* 51 */ {"SPFH", "Specific humidity [kg/kg]"},
   /* 52 */ {"RH", "Relative humidity [%]"},
   /* 53 */ {"MIXR", "Humidity mixing ratio [kg/kg]"},
   /* 54 */ {"PWAT", "Precipitable water [kg/m^2]"},
   /* 55 */ {"VAPP", "Vapor pressure [Pa]"},
   /* 56 */ {"SATD", "Saturation deficit [Pa]"},
   /* 57 */ {"EVP", "Evaporation [kg/m^2]"},
   /* 58 */ {"CICE", "Cloud Ice [kg/m^2]"},
   /* 59 */ {"PRATE", "Precipitation rate [kg/m^2/s]"},
   /* 60 */ {"TSTM", "Thunderstorm probability [%]"},
   /* 61 */ {"APCP", "Total precipitation [kg/m^2]"},
   /* 62 */ {"NCPCP", "Large scale precipitation [kg/m^2]"},
   /* 63 */ {"ACPCP", "Convective precipitation [kg/m^2]"},
   /* 64 */ {"SRWEQ", "Snowfall rate water equiv. [kg/m^2/s]"},
   /* 65 */ {"WEASD", "Accum. snow [kg/m^2]"},
   /* 66 */ {"SNOD", "Snow depth [m]"},
   /* 67 */ {"MIXHT", "Mixed layer depth [m]"},
   /* 68 */ {"TTHDP", "Transient thermocline depth [m]"},
   /* 69 */ {"MTHD", "Main thermocline depth [m]"},
   /* 70 */ {"MTHA", "Main thermocline anomaly [m]"},
   /* 71 */ {"TCDC", "Total cloud cover [%]"},
   /* 72 */ {"CDCON", "Convective cloud cover [%]"},
   /* 73 */ {"LCDC", "Low level cloud cover [%]"},
   /* 74 */ {"MCDC", "Mid level cloud cover [%]"},
   /* 75 */ {"HCDC", "High level cloud cover [%]"},
   /* 76 */ {"CWAT", "Cloud water [kg/m^2]"},
   /* 77 */ {"var77", "undefined"},
   /* 78 */ {"SNOC", "Convective snow [kg/m^2]"},
   /* 79 */ {"SNOL", "Large scale snow [kg/m^2]"},
   /* 80 */ {"WTMP", "Water temp. [K]"},
   /* 81 */ {"LAND", "Land-sea mask [1=land; 0=sea]"},
   /* 82 */ {"DSLM", "Deviation of sea level from mean [m]"},
   /* 83 */ {"SFCR", "Surface roughness [m]"},
   /* 84 */ {"ALBDO", "Albedo [%]"},
   /* 85 */ {"TSOIL", "Soil temp. [K]"},
   /* 86 */ {"SOILM", "Soil moisture content [kg/m^2]"},
   /* 87 */ {"VEG", "Vegetation [%]"},
   /* 88 */ {"SALTY", "Salinity [kg/kg]"},
   /* 89 */ {"DEN", "Density [kg/m^3]"},
   /* 90 */ {"RUNOF", "Runoff [kg/m^2]"},
   /* 91 */ {"ICEC", "Ice concentration [ice=1;no ice=0]"},
   /* 92 */ {"ICETK", "Ice thickness [m]"},
   /* 93 */ {"DICED", "Direction of ice drift [deg]"},
   /* 94 */ {"SICED", "Speed of ice drift [m/s]"},
   /* 95 */ {"UICE", "u of ice drift [m/s]"},
   /* 96 */ {"VICE", "v of ice drift [m/s]"},
   /* 97 */ {"ICEG", "Ice growth rate [m/s]"},
   /* 98 */ {"ICED", "Ice divergence [/s]"},
   /* 99 */ {"SNOM", "Snow melt [kg/m^2]"},
   /* 100 */ {"HTSGW", "Sig height of wind waves and swell [m]"},
   /* 101 */ {"WVDIR", "Direction of wind waves [deg]"},
   /* 102 */ {"WVHGT", "Sig height of wind waves [m]"},
   /* 103 */ {"WVPER", "Mean period of wind waves [s]"},
   /* 104 */ {"SWDIR", "Direction of swell waves [deg]"},
   /* 105 */ {"SWELL", "Sig height of swell waves [m]"},
   /* 106 */ {"SWPER", "Mean period of swell waves [s]"},
   /* 107 */ {"DIRPW", "Primary wave direction [deg]"},
   /* 108 */ {"PERPW", "Primary wave mean period [s]"},
   /* 109 */ {"DIRSW", "Secondary wave direction [deg]"},
   /* 110 */ {"PERSW", "Secondary wave mean period [s]"},
   /* 111 */ {"NSWRS", "Net short wave (surface) [W/m^2]"},
   /* 112 */ {"NLWRS", "Net long wave (surface) [W/m^2]"},
   /* 113 */ {"NSWRT", "Net short wave (top) [W/m^2]"},
   /* 114 */ {"NLWRT", "Net long wave (top) [W/m^2]"},
   /* 115 */ {"LWAVR", "Long wave [W/m^2]"},
   /* 116 */ {"SWAVR", "Short wave [W/m^2]"},
   /* 117 */ {"GRAD", "Global radiation [W/m^2]"},
   /* 118 */ {"var118", "undefined"},
   /* 119 */ {"var119", "undefined"},
   /* 120 */ {"var120", "undefined"},
   /* 121 */ {"LHTFL", "Latent heat flux [W/m^2]"},
   /* 122 */ {"SHTFL", "Sensible heat flux [W/m^2]"},
   /* 123 */ {"BLYDP", "Boundary layer dissipation [W/m^2]"},
   /* 124 */ {"UFLX", "Zonal momentum flux [N/m^2]"},
   /* 125 */ {"VFLX", "Meridional momentum flux [N/m^2]"},
   /* 126 */ {"WMIXE", "Wind mixing energy [J]"},
   /* 127 */ {"IMGD", "Image data [integer]"},
   /* 128 */ {"MSLSA", "Mean sea level pressure (Std Atm) [Pa]"},
   /* 129 */ {"MSLMA", "Mean sea level pressure (MAPS) [Pa]"},
   /* 130 */ {"MSLET", "Mean sea level pressure (ETA model) [Pa]"},
   /* 131 */ {"LFTX", "Surface lifted index [K]"},
   /* 132 */ {"4LFTX", "Best (4-layer) lifted index [K]"},
   /* 133 */ {"KX", "K index [K]"},
   /* 134 */ {"SX", "Sweat index [K]"},
   /* 135 */ {"MCONV", "Horizontal moisture divergence [kg/kg/s]"},
   /* 136 */ {"VSSH", "Vertical speed shear [1/s]"},
   /* 137 */ {"TSLSA", "3-hr pressure tendency [Pa/s]"},
   /* 138 */ {"BVF2", "Brunt-Vaisala frequency^2 [1/s^2]"},
   /* 139 */ {"PVMW", "Potential vorticity (mass-weighted) [1/s/m]"},
   /* 140 */ {"CRAIN", "Categorical rain [yes=1;no=0]"},
   /* 141 */ {"CFRZR", "Categorical freezing rain [yes=1;no=0]"},
   /* 142 */ {"CICEP", "Categorical ice pellets [yes=1;no=0]"},
   /* 143 */ {"CSNOW", "Categorical snow [yes=1;no=0]"},
   /* 144 */ {"SOILW", "Volumetric soil moisture [fraction]"},
   /* 145 */ {"PEVPR", "Potential evaporation rate [W/m^2]"},
   /* 146 */ {"CWORK", "Cloud work function [J/kg]"},
   /* 147 */ {"U-GWD", "Zonal gravity wave stress [N/m^2]"},
   /* 148 */ {"V-GWD", "Meridional gravity wave stress [N/m^2]"},
   /* 149 */ {"PV___", "Potential vorticity [m^2/s/kg]"},
   /* 150 */ {"var150", "undefined"},
   /* 151 */ {"var151", "undefined"},
   /* 152 */ {"var152", "undefined"},
   /* 153 */ {"MFXDV", "Moisture flux divergence [gr/gr*m/s/m]"},
   /* 154 */ {"var154", "undefined"},
   /* 155 */ {"GFLUX", "Ground heat flux [W/m^2]"},
   /* 156 */ {"CIN", "Convective inhibition [J/kg]"},
   /* 157 */ {"CAPE", "Convective Avail. Pot. Energy [J/kg]"},
   /* 158 */ {"TKE", "Turbulent kinetic energy [J/kg]"},
   /* 159 */ {"CONDP", "Lifted parcel condensation pressure [Pa]"},
   /* 160 */ {"CSUSF", "Clear sky upward solar flux [W/m^2]"},
   /* 161 */ {"CSDSF", "Clear sky downward solar flux [W/m^2]"},
   /* 162 */ {"CSULF", "Clear sky upward long wave flux [W/m^2]"},
   /* 163 */ {"CSDLF", "Clear sky downward long wave flux [W/m^2]"},
   /* 164 */ {"CFNSF", "Cloud forcing net solar flux [W/m^2]"},
   /* 165 */ {"CFNLF", "Cloud forcing net long wave flux [W/m^2]"},
   /* 166 */ {"VBDSF", "Visible beam downward solar flux [W/m^2]"},
   /* 167 */ {"VDDSF", "Visible diffuse downward solar flux [W/m^2]"},
   /* 168 */ {"NBDSF", "Near IR beam downward solar flux [W/m^2]"},
   /* 169 */ {"NDDSF", "Near IR diffuse downward solar flux [W/m^2]"},
   /* 170 */ {"USTR", "U wind stress [N/m^2]"},
   /* 171 */ {"VSTR", "V wind stress [N/m^2]"},
   /* 172 */ {"MFLX", "Momentum flux [N/m^2]"},
   /* 173 */ {"LMH", "Mass point model surface [integer]"},
   /* 174 */ {"LMV", "Velocity point model surface [integer]"},
   /* 175 */ {"SGLYR", "Nearby model level [integer]"},
   /* 176 */ {"NLAT", "Latitude [deg]"},
   /* 177 */ {"NLON", "Longitude [deg]"},
   /* 178 */ {"UMAS", "Mass weighted u [gm/m*K*s]"},
   /* 179 */ {"VMAS", "Mass weighted v [gm/m*K*s]"},
   /* 180 */ {"XPRATE", "corrected precip [kg/m^2/s]"},
   /* 181 */ {"LPSX", "x-gradient of log pressure [1/m]"},
   /* 182 */ {"LPSY", "y-gradient of log pressure [1/m]"},
   /* 183 */ {"HGTX", "x-gradient of height [m/m]"},
   /* 184 */ {"HGTY", "y-gradient of height [m/m]"},
   /* 185 */ {"STDZ", "Std dev of Geop. hgt. [m]"},
   /* 186 */ {"STDU", "Std dev of zonal wind [m/s]"},
   /* 187 */ {"STDV", "Std dev of meridional wind [m/s]"},
   /* 188 */ {"STDQ", "Std dev of spec. hum. [gm/gm]"},
   /* 189 */ {"STDT", "Std dev of temp. [K]"},
   /* 190 */ {"CBUW", "Covar. u and omega [m/s*Pa/s]"},
   /* 191 */ {"CBVW", "Covar. v and omega [m/s*Pa/s]"},
   /* 192 */ {"CBUQ", "Covar. u and specific hum [m/s*gm/gm]"},
   /* 193 */ {"CBVQ", "Covar. v and specific hum [m/s*gm/gm]"},
   /* 194 */ {"CBTW", "Covar. T and omega [K*Pa/s]"},
   /* 195 */ {"CBQW", "Covar. spec. hum and omega [gm/gm*Pa/s]"},
   /* 196 */ {"CBMZW", "Covar. v and u [m^2/s^2]"},
   /* 197 */ {"CBTZW", "Covar. u and T [K*m/s]"},
   /* 198 */ {"CBTMW", "Covar. v and T [K*m/s]"},
   /* 199 */ {"STDRH", "Std dev of Rel. Hum. [%]"},
   /* 200 */ {"SDTZ", "Std dev of time tend of geop. hgt [m]"},
   /* 201 */ {"ICWAT", "Ice-free water surface [%]"},
   /* 202 */ {"SDTU", "Std dev of time tend of zonal wind [m/s]"},
   /* 203 */ {"SDTV", "Std dev of time tend of merid wind [m/s]"},
   /* 204 */ {"DSWRF", "Downward solar radiation flux [W/m^2]"},
   /* 205 */ {"DLWRF", "Downward long wave flux [W/m^2]"},
   /* 206 */ {"SDTQ", "Std dev of time tend of spec. hum [gm/gm]"},
   /* 207 */ {"MSTAV", "Moisture availability [%]"},
   /* 208 */ {"SFEXC", "Exchange coefficient [kg*m/m^3/s]"},
   /* 209 */ {"MIXLY", "No. of mixed layers next to sfc [integer]"},
   /* 210 */ {"SDTT", "Std dev of time tend of temp. [K]"},
   /* 211 */ {"USWRF", "Upward solar radiation flux [W/m^2]"},
   /* 212 */ {"ULWRF", "Upward long wave flux [W/m^2]"},
   /* 213 */ {"CDLYR", "Non-convective cloud [%]"},
   /* 214 */ {"CPRAT", "Convective precip. rate [kg/m^2/s]"},
   /* 215 */ {"TTDIA", "Temp. tendency by all physics [K/s]"},
   /* 216 */ {"TTRAD", "Temp. tendency by all radiation [K/s]"},
   /* 217 */ {"TTPHY", "Temp. tendency by nonrad physics [K/s]"},
   /* 218 */ {"PREIX", "Precipitation index [fraction]"},
   /* 219 */ {"TSD1D", "Std dev of IR T over 1x1 deg area [K]"},
   /* 220 */ {"NLSGP", "Natural log of surface pressure [ln(kPa)]"},
   /* 221 */ {"SDTRH", "Std dev of time tend of rel hum [%]"},
   /* 222 */ {"5WAVH", "5-wave geopotential height [gpm]"},
   /* 223 */ {"CNWAT", "Plant canopy surface water [kg/m^2]"},
   /* 224 */ {"PLTRS", "Max. stomato plant resistance [s/m]"},
   /* 225 */ {"RHCLD", "RH-type cloud cover [%]"},
   /* 226 */ {"BMIXL", "Blackadar's mixing length scale [m]"},
   /* 227 */ {"AMIXL", "Asymptotic mixing length scale [m]"},
   /* 228 */ {"PEVAP", "Pot. evaporation [kg/m^2]"},
   /* 229 */ {"SNOHF", "Snow melt heat flux [W/m^2]"},
   /* 230 */ {"SNOEV", "Snow sublimation heat flux [W/m^2]"},
   /* 231 */ {"MFLUX", "Convective cloud mass flux [Pa/s]"},
   /* 232 */ {"DTRF", "Downward total radiation flux [W/m^2]"},
   /* 233 */ {"UTRF", "Upward total radiation flux [W/m^2]"},
   /* 234 */ {"BGRUN", "Baseflow-groundwater runoff [kg/m^2]"},
   /* 235 */ {"SSRUN", "Storm surface runoff [kg/m^2]"},
   /* 236 */ {"var236", "undefined"},
   /* 237 */ {"OZONE", "Total column ozone [Dobson]"},
   /* 238 */ {"SNOWC", "Snow cover [%]"},
   /* 239 */ {"SNOT", "Snow temp. [K]"},
   /* 240 */ {"GLCR", "Permanent snow points [mask]"},
   /* 241 */ {"LRGHR", "Large scale condensation heating [K/s]"},
   /* 242 */ {"CNVHR", "Deep convective heating [K/s]"},
   /* 243 */ {"CNVMR", "Deep convective moistening [kg/kg/s]"},
   /* 244 */ {"SHAHR", "Shallow convective heating [K/s]"},
   /* 245 */ {"SHAMR", "Shallow convective moistening [kg/kg/s]"},
   /* 246 */ {"VDFHR", "Vertical diffusion heating [K/s]"},
   /* 247 */ {"VDFUA", "Vertical diffusion zonal accel [m/s^2]"},
   /* 248 */ {"VDFVA", "Vertical diffusion meridional accel [m/s^2]"},
   /* 249 */ {"VDFMR", "Vertical diffusion moistening [kg/kg/s]"},
   /* 250 */ {"SWHR", "Solar radiative heating [K/s]"},
   /* 251 */ {"LWHR", "Longwave radiative heating [K/s]"},
   /* 252 */ {"CD", "Drag coefficient [non-dim]"},
   /* 253 */ {"FRICV", "Friction velocity [m/s]"},
   /* 254 */ {"RI", "Richardson number [non-dim]"},
   /* 255 */ {"var255", "undefined"},
};

struct ParmTable parm_table_ecmwf_128[256] = {
  {"var0", "undefined"}, /* 0 */
  {"PRES", "Pressure [Pa]"},     /* 1 */
  {"PRMSL", "Pressure reduced to MSL [Pa]"},     /* 2 */
  {"PTEND", "Pressure tendency [Pa/s]"},     /* 3 */
  {"var4", "undefined"}, /* 4 */
  {"var5", "undefined"}, /* 5 */
  {"GP", "Geopotential [m**2/s**2]"},     /* 6 */
  {"HGT", "Geopotential height [gpm]"},     /* 7 */
  {"DIST", "Geometric height [m]"},     /* 8 */
  {"HSTDV", "Standard deviation of height [m]"},     /* 9 */
  {"HVAR", "Variance of height [m**2]"},     /* 10 */
  {"TMP", "Temperature [K]"},     /* 11 */
  {"VTMP", "Virtual temperature [K]"},     /* 12 */
  {"POT", "Potential temperature [K]"},     /* 13 */
  {"EPOT", "Pseudo-adiabatic potential temperature [K]"},     /* 14 */
  {"TMAX", "Maximum temperature [K]"},     /* 15 */
  {"TMIN", "Minimum temperature [K]"},     /* 16 */
  {"DPT", "Dew point temperature [K]"},     /* 17 */
  {"DEPR", "Dew point depression [K]"},     /* 18 */
  {"LAPR", "Lapse rate [K/m]"},     /* 19 */
  {"VISIB", "Visibility [m]"},     /* 20 */
  {"RDSP1", "Radar spectra (1) [dimensionless]"},     /* 21 */
  {"RDSP2", "Radar spectra (2) [dimensionless]"},     /* 22 */
  {"RDSP3", "Radar spectra (3) [dimensionless]"},     /* 23 */
  {"var24", "undefined"}, /* 24 */
  {"TMPA", "Temperature anomaly [K]"},     /* 25 */
  {"PRESA", "Pressure anomaly [Pa]"},     /* 26 */
  {"GPA", "Geopotential height anomaly [gpm]"},     /* 27 */
  {"WVSP1", "Wave spectra (1) [dimensionless]"},     /* 28 */
  {"WVSP2", "Wave spectra (2) [dimensionless]"},     /* 29 */
  {"WVSP3", "Wave spectra (3) [dimensionless]"},     /* 30 */
  {"WDIR", "Wind direction [degree]"},     /* 31 */
  {"WIND", "Wind speed [m/s]"},     /* 32 */
  {"UGRD", "u wind [m/s]"},     /* 33 */
  {"VGRD", "v wind [m/s]"},     /* 34 */
  {"STRM", "Stream function [m**2/s]"},     /* 35 */
  {"VPOT", "Velocity potential [m**2/s]"},     /* 36 */
  {"MNTSF", "Montgomery stream function [m**2/s**2]"},     /* 37 */
  {"SGCVV", "Sigma coord. vertical velocity [/s]"},     /* 38 */
  {"VVEL", "Pressure vertical velocity [Pa/s]"},     /* 39 */
  {"DZDT", "Geometric vertical velocity [m/s]"},     /* 40 */
  {"ABSV", "Absolute vorticity [/s]"},     /* 41 */
  {"ABSD", "Absolute divergence [/s]"},     /* 42 */
  {"RELV", "Relative vorticity [/s]"},     /* 43 */
  {"RELD", "Relative divergence [/s]"},     /* 44 */
  {"VUCSH", "Vertical u shear [/s]"},     /* 45 */
  {"VVCSH", "Vertical v shear [/s]"},     /* 46 */
  {"DIRC", "Direction of current [degree]"},     /* 47 */
  {"SPC", "Speed of current [m/s]"},     /* 48 */
  {"UOGRD", "u of current [m/s]"},     /* 49 */
  {"VOGRD", "v of current [m/s]"},     /* 50 */
  {"SPFH", "Specific humidity [kg/kg]"},     /* 51 */
  {"RH", "Relative humidity [percent]"},     /* 52 */
  {"MIXR", "Humidity mixing ratio [kg/kg]"},     /* 53 */
  {"PWAT", "Precipitable water [kg/m**2]"},     /* 54 */
  {"VAPP", "Vapor pressure [Pa]"},     /* 55 */
  {"SATD", "Saturation deficit [Pa]"},     /* 56 */
  {"EVP", "Evaporation [kg/m**2]"},     /* 57 */
  {"CICE", "Cloud Ice [kg/m**2]"},     /* 58 */
  {"PRATE", "Precipitation rate [kg/m**2/s]"},     /* 59 */
  {"TSTM", "Thunderstorm probability [percent]"},     /* 60 */
  {"APCP", "Total precipitation [kg/m**2]"},     /* 61 */
  {"NCPCP", "Large scale precipitation [kg/m**2]"},     /* 62 */
  {"ACPCP", "Convective precipitation [kg/m**2]"},     /* 63 */
  {"SRWEQ", "Snowfall rate water equivalent [kg/m**2/s]"},     /* 64 */
  {"WEASD", "Water equiv. of accum. snow depth [kg/m**2]"},     /* 65 */
  {"SNOD", "Snow depth [m]"},     /* 66 */
  {"MIXHT", "Mixed layer depth [m]"},     /* 67 */
  {"TTHDP", "Transient thermocline depth [m]"},     /* 68 */
  {"MTHD", "Main thermocline depth [m]"},     /* 69 */
  {"MTHA", "Main thermocline anomaly [m]"},     /* 70 */
  {"TCDC", "Total cloud cover [percent]"},     /* 71 */
  {"CDCON", "Convective cloud cover [percent]"},     /* 72 */
  {"LCDC", "Low level cloud cover [percent]"},     /* 73 */
  {"MCDC", "Mid level cloud cover [percent]"},     /* 74 */
  {"HCDC", "High level cloud cover [percent]"},     /* 75 */
  {"CWAT", "Cloud water [kg/m**2]"},     /* 76 */
  {"var77", "undefined"}, /* 77 */
  {"SNOC", "Convective snow [kg/m**2]"},     /* 78 */
  {"SNOL", "Large scale snow [kg/m**2]"},     /* 79 */
  {"WTMP", "Water temperature [K]"},     /* 80 */
  {"LAND", "Land-sea mask (1=land; 0=sea) [integer]"},     /* 81 */
  {"DSLM", "Deviation of sea level from mean [m]"},     /* 82 */
  {"SFCR", "Surface roughness [m]"},     /* 83 */
  {"ALBDO", "Albedo [percent]"},     /* 84 */
  {"TSOIL", "Soil temperature [K]"},     /* 85 */
  {"SOILM", "Soil moisture content [kg/m**2]"},     /* 86 */
  {"VEG", "Vegetation [percent]"},     /* 87 */
  {"SALTY", "Salinity [kg/kg]"},     /* 88 */
  {"DEN", "Density [kg/m**2]"},     /* 89 */
  {"RUNOF", "Runoff [kg/m**2]"},     /* 90 */
  {"ICEC", "Ice concentration (ice=1; no ice=0) [1/0]"},     /* 91 */
  {"ICETK", "Ice thickness [m]"},     /* 92 */
  {"DICED", "Direction of ice drift [degree]"},     /* 93 */
  {"SICED", "Speed of ice drift [m/s]"},     /* 94 */
  {"UICE", "u of ice drift [m/s]"},     /* 95 */
  {"VICE", "v of ice drift [m/s]"},     /* 96 */
  {"ICEG", "Ice growth [m]"},     /* 97 */
  {"ICED", "Ice divergence [/s]"},     /* 98 */
  {"SNOM", "Snow melt [kg/m**2]"},     /* 99 */
  {"HTSGW", "Sig height of wind waves and swell [m]"},     /* 100 */
  {"WVDIR", "Direction of wind waves [degree]"},     /* 101 */
  {"WVHGT", "Significant height of wind waves [m]"},     /* 102 */
  {"WVPER", "Mean period of wind waves [s]"},     /* 103 */
  {"SWDIR", "Direction of swell waves [degree]"},     /* 104 */
  {"SWELL", "Significant height of swell waves [m]"},     /* 105 */
  {"SWPER", "Mean period of swell waves [s]"},     /* 106 */
  {"DIRPW", "Primary wave direction [degree]"},     /* 107 */
  {"PERPW", "Primary wave mean period [s]"},     /* 108 */
  {"DIRSW", "Secondary wave direction [degree]"},     /* 109 */
  {"PERSW", "Secondary wave mean period [s]"},     /* 110 */
  {"NSWRS", "Net short wave radiation (surface) [W/m**2]"},     /* 111 */
  {"NLWRS", "Net long wave radiation (surface) [W/m**2]"},     /* 112 */
  {"NSWRT", "Net short wave radiation (top) [W/m**2]"},     /* 113 */
  {"NLWRT", "Net long wave radiation (top) [W/m**2]"},     /* 114 */
  {"LWAVR", "Long wave radiation [W/m**2]"},     /* 115 */
  {"SWAVR", "Short wave radiation [W/m**2]"},     /* 116 */
  {"GRAD", "Global radiation [W/m**2]"},     /* 117 */
  {"var118", "undefined"}, /* 118 */
  {"var119", "undefined"}, /* 119 */
  {"var120", "undefined"}, /* 120 */
  {"LHTFL", "Latent heat flux [W/m**2]"},     /* 121 */
  {"SHTFL", "Sensible heat flux [W/m**2]"},     /* 122 */
  {"BLYDP", "Boundary layer dissipation [W/m**2]"},     /* 123 */
  {"UFLX", "Zonal component of momentum flux [N/m**2]"},     /* 124 */
  {"VFLX", "Meridional component of momentum flux [N/m**2]"},     /* 125 */
  {"WMIXE", "Wind mixing energy [J]"},     /* 126 */

  /* ECMWF local table */

  {"AT",      /* 127 */ "Atmospheric Tide"},
  {"BV",      /* 128 */ "Budget Values"},
  {"Z",       /* 129 */ "Geopotential [m2 s-2]"},
  {"T",       /* 130 */ "Temperature [K]"},
  {"U",       /* 131 */ "U-component of Wind [ms-1]"},
  {"V",       /* 132 */ "V-component of Wind [ms-1]"},
  {"Q",       /* 133 */ "Specific Humidity [kg/kg]"},
  {"SP",      /* 134 */ "Surface Pressure [Pa]"},
  {"W",       /* 135 */ "Vertical Velocity [Pa s-1]"},
  {"TCW",     /* 136 */ "Total column water (vapor+drops+ice) [kg/m2]"},
  {"TCWV",    /* 137 */ "Total column water vapor [kg/m2]"},
  {"VO",      /* 138 */ "relative vorticity [s-1]"},
  {"STL1",    /* 139 */ "soil temperature level 1 [K]"},
  {"SWL1",    /* 140 */ "soil moisture level 1 [m (H20)]"},
  {"SD",      /* 141 */ "Snow Depth [m]"},
  {"LSP",     /* 142 */ "Large Scale Precipitation [m]"},
  {"CP",      /* 143 */ "Convective Precipitation [m]"},
  {"SF",      /* 144 */ "Snow Fall [m]"},
  {"BLD",     /* 145 */ "Boundary Layer Dissipation [Wm-2]"},
  {"SSHF",    /* 146 */ "Surface Flux of Sensible Heat [Wm-2]"},
  {"SLHF",    /* 147 */ "Surface Flux of Latent Heat [Wm-2]"},
  {"v148",    /* 148 */ "undefined"},
  {"v149",    /* 149 */ "undefined"},
  {"v150",    /* 150 */ "undefined"},
  {"MSL",     /* 151 */ "Mean Sea Level (MSL) Pressure [Pa]    Pa"},
  {"LNSP",    /* 152 */ "Log Surface Pressure"},
  {"v153",    /* 153 */ "undefined"},
  {"v154",    /* 154 */ "undefined"},
  {"D",       /* 155 */ "Divergence [s-1]"},
  {"GH",      /* 156 */ "Height (Geopotential) [m]"},
  {"R",       /* 157 */ "Relative Humidity [%]"},
  {"TSP",     /* 158 */ "Tendency of Surface Pressure [Pa s-1]"},
  {"v159",    /* 159 */ "undefined"},
  {"SDOR",    /* 160 */ "Standard deviation of orography"},
  {"ISOR",    /* 161 */ "Anisotropy of subgrid scale orography"},
  {"ANOR",    /* 162 */ "Angle of subgrid scale orography"},
  {"SLOR",    /* 163 */ "Slope of subgrid scale orography"},
  {"TCC",     /* 164 */ "cloud cover total [0-1]"},
  {"10U",     /* 165 */ "U-wind at 10 [ms-1]"},
  {"10V",     /* 166 */ "V-wind at 10 [ms-1]"},
  {"2T",      /* 167 */ "Temperature at 2 m [K]"},
  {"2D",      /* 168 */ "Dewpoint at 2 m [K]"},
  {"v169",    /* 169 */ "undefined"},
  {"STL2",    /* 170 */ "soil temperature level 2 [K]"},
  {"SWL2",    /* 171 */ "soil wetness level 2 [m (H20)]"},
  {"LSM",     /* 172 */ "land-sea mask [(0,1)]"},
  {"SR",      /* 173 */ "sfc roughness [m]"},
  {"AL",      /* 174 */ "Albedo [0-1]"},
  {"v175",    /* 175 */ "undefined"},
  {"SSR",     /* 176 */ "Net Shortwave Radiation (surface) [Wm-2]"},
  {"STR",     /* 177 */ "Net Longwave Radiation (surface) [Wm-2]"},
  {"TSR",     /* 178 */ "Net Shortwave Radiation (toa) [Wm-2"},
  {"TTR",     /* 179 */ "Net Longwave Radiation (toa) [Wm-2]"},
  {"EWSS",    /* 180 */ "U-component of Surface Wind Stress [Nm-2]"},
  {"NSSS",    /* 181 */ "V-component of Surface Wind Stress [Nm-2]"},
  {"E",       /* 182 */ "Evaporation [m (H2O)]"},
  {"STL3",    /* 183 */ "soil temp level 3 [m (H2O)]"},
  {"SWL3",    /* 184 */ "soil moisture level 3 [K]"},
  {"CCC",     /* 185 */ "cloud convective [0-1]"},
  {"LCC",     /* 186 */ "cloud low [0-1]"},
  {"MCC",     /* 187 */ "cloud mid [0-1]"},
  {"HCC",     /* 188 */ "cloud high [0-1]"},
  {"v189",    /* 189 */ "undefined"},
  {"EWOV",    /* 190 */ "orographic variance e-w [m2]"},
  {"NSOV",    /* 191 */ "orographic variance n-s [m2]"},
  {"NWOV",    /* 192 */ "orographic variance nw-se [m2]"},
  {"NEOV",    /* 193 */ "orographic variance ne-sw [m2]"},
  {"v194",    /* 194 */ "undefined"},
  {"LGWS",    /* 195 */ "gravity wave stress n-s [n/m2s]"},
  {"MGWS",    /* 196 */ "gravity wave stress e-w [n/m2s]"},
  {"GWD",     /* 197 */ "gravity wave diss [w/m2s]"},
  {"SRC",     /* 198 */ "skin resevoir content [m]"},
  {"VEG",     /* 199 */ "sfc vegetation cover [%]"},
  {"VSO",     /* 200 */ "variance of subgrid scale orgography [m2]"},
  {"MX2T",    /* 201 */ "max 2m temp [K]"},
  {"MN2T",    /* 202 */ "min 2m temp [K]"},
  {"v203",    /* 203 */ "undefined"},
  {"PAW",     /* 204 */ "precip analysis weights"},
  {"RO",      /* 205 */ "runoff [m]"},
  {"v206",    /* 206 */ "undefined"},
  {"v207",    /* 207 */ "undefined"},
  {"v208",    /* 208 */ "undefined"},
  {"v209",    /* 209 */ "undefined"},
  {"v210",    /* 210 */ "undefined"},
  {"v211",    /* 211 */ "undefined"},
  {"v212",    /* 212 */ "undefined"},
  {"v213",    /* 213 */ "undefined"},
  {"v214",    /* 214 */ "undefined"},
  {"v215",    /* 215 */ "undefined"},
  {"v216",    /* 216 */ "undefined"},
  {"v217",    /* 217 */ "undefined"},
  {"v218",    /* 218 */ "undefined"},
  {"v219",    /* 219 */ "undefined"},
  {"v220",    /* 220 */ "undefined"},
  {"v221",    /* 221 */ "undefined"},
  {"v222",    /* 222 */ "undefined"},
  {"v223",    /* 223 */ "undefined"},
  {"v224",    /* 224 */ "undefined"},
  {"v225",    /* 225 */ "undefined"},
  {"v226",    /* 226 */ "undefined"},
  {"v227",    /* 227 */ "undefined"},
  {"TP",      /* 228 */ "total precip [m]"},
  {"IEWS",   /* 229 */ "instanteous sfc stress u [Nm-2]"},
  {"INSS",   /* 230 */ "instanteous sfc stress v [Nm-2]"},
  {"ISHF",   /* 231 */ "instanteous sfc sensible heat flux [Wm-2]"},
  {"IE",   /* 232 */ "instanteous sfc latent heat flux [kg/m2s]"},
  {"ASQ",    /* 233 */ "apparent sfc humidity [kg/kg]"},
  {"LSRH",  /* 234 */ "log sfc roughness"},
  {"SKT",     /* 235 */ "skin temperature [K]"},
  {"STL4",    /* 236 */ "soil temperature level 4 [K]"},
  {"SWL4",   /* 237 */ "soil wetness level 4 [m (H2O)]"},
  {"TSN",     /* 238 */ "t of snow layer [K]"},
  {"CSF",    /* 239 */ "convective snow [m]"},
  {"LSF",    /* 240 */ "large scale snow [m]"},
  {"v241",    /* 241 */ "undefined"},
  {"v242",    /* 242 */ "undefined"},
  {"FAL",   /* 243 */ "forecast albedo"},
  {"FSR",    /* 244 */ "forecast sfc roughness [m]"},
  {"FLSR",/* 245 */ "log of forecast sfc roughness"},
  {"CLWC",   /* 246 */ "Cloud liquid water content [kg/kg]"},
  {"CIWC",    /* 247 */ "Cloud ice water content [kg/kg]"},
  {"CC",    /* 248 */ "Cloud cover [0-1]"},
  {"v249",    /* 249 */ "undefined"},
  {"v250",    /* 250 */ "Ice age (0 first year, 1 multi year) [0,1]"},
  {"v251",    /* 251 */ "undefined"},
  {"v252",    /* 252 */ "undefined"},
  {"v253",    /* 253 */ "undefined"},
  {"v254",    /* 254 */ "undefined"},
  {"v255",    /* 255 */ "undefined"}
};


struct ParmTable parm_table_ecmwf_160[256] = {
  {"var0", "undefined"}, /* 0 */
  {"PRES", "Pressure [Pa]"},     /* 1 */
  {"PRMSL", "Pressure reduced to MSL [Pa]"},     /* 2 */
  {"PTEND", "Pressure tendency [Pa/s]"},     /* 3 */
  {"var4", "undefined"}, /* 4 */
  {"var5", "undefined"}, /* 5 */
  {"GP", "Geopotential [m**2/s**2]"},     /* 6 */
  {"HGT", "Geopotential height [gpm]"},     /* 7 */
  {"DIST", "Geometric height [m]"},     /* 8 */
  {"HSTDV", "Standard deviation of height [m]"},     /* 9 */
  {"HVAR", "Variance of height [m**2]"},     /* 10 */
  {"TMP", "Temperature [K]"},     /* 11 */
  {"VTMP", "Virtual temperature [K]"},     /* 12 */
  {"POT", "Potential temperature [K]"},     /* 13 */
  {"EPOT", "Pseudo-adiabatic potential temperature [K]"},     /* 14 */
  {"TMAX", "Maximum temperature [K]"},     /* 15 */
  {"TMIN", "Minimum temperature [K]"},     /* 16 */
  {"DPT", "Dew point temperature [K]"},     /* 17 */
  {"DEPR", "Dew point depression [K]"},     /* 18 */
  {"LAPR", "Lapse rate [K/m]"},     /* 19 */
  {"VISIB", "Visibility [m]"},     /* 20 */
  {"RDSP1", "Radar spectra (1) [dimensionless]"},     /* 21 */
  {"RDSP2", "Radar spectra (2) [dimensionless]"},     /* 22 */
  {"RDSP3", "Radar spectra (3) [dimensionless]"},     /* 23 */
  {"var24", "undefined"}, /* 24 */
  {"TMPA", "Temperature anomaly [K]"},     /* 25 */
  {"PRESA", "Pressure anomaly [Pa]"},     /* 26 */
  {"GPA", "Geopotential height anomaly [gpm]"},     /* 27 */
  {"WVSP1", "Wave spectra (1) [dimensionless]"},     /* 28 */
  {"WVSP2", "Wave spectra (2) [dimensionless]"},     /* 29 */
  {"WVSP3", "Wave spectra (3) [dimensionless]"},     /* 30 */
  {"WDIR", "Wind direction [degree]"},     /* 31 */
  {"WIND", "Wind speed [m/s]"},     /* 32 */
  {"UGRD", "u wind [m/s]"},     /* 33 */
  {"VGRD", "v wind [m/s]"},     /* 34 */
  {"STRM", "Stream function [m**2/s]"},     /* 35 */
  {"VPOT", "Velocity potential [m**2/s]"},     /* 36 */
  {"MNTSF", "Montgomery stream function [m**2/s**2]"},     /* 37 */
  {"SGCVV", "Sigma coord. vertical velocity [/s]"},     /* 38 */
  {"VVEL", "Pressure vertical velocity [Pa/s]"},     /* 39 */
  {"DZDT", "Geometric vertical velocity [m/s]"},     /* 40 */
  {"ABSV", "Absolute vorticity [/s]"},     /* 41 */
  {"ABSD", "Absolute divergence [/s]"},     /* 42 */
  {"RELV", "Relative vorticity [/s]"},     /* 43 */
  {"RELD", "Relative divergence [/s]"},     /* 44 */
  {"VUCSH", "Vertical u shear [/s]"},     /* 45 */
  {"VVCSH", "Vertical v shear [/s]"},     /* 46 */
  {"DIRC", "Direction of current [degree]"},     /* 47 */
  {"SPC", "Speed of current [m/s]"},     /* 48 */
  {"UOGRD", "u of current [m/s]"},     /* 49 */
  {"VOGRD", "v of current [m/s]"},     /* 50 */
  {"SPFH", "Specific humidity [kg/kg]"},     /* 51 */
  {"RH", "Relative humidity [percent]"},     /* 52 */
  {"MIXR", "Humidity mixing ratio [kg/kg]"},     /* 53 */
  {"PWAT", "Precipitable water [kg/m**2]"},     /* 54 */
  {"VAPP", "Vapor pressure [Pa]"},     /* 55 */
  {"SATD", "Saturation deficit [Pa]"},     /* 56 */
  {"EVP", "Evaporation [kg/m**2]"},     /* 57 */
  {"CICE", "Cloud Ice [kg/m**2]"},     /* 58 */
  {"PRATE", "Precipitation rate [kg/m**2/s]"},     /* 59 */
  {"TSTM", "Thunderstorm probability [percent]"},     /* 60 */
  {"APCP", "Total precipitation [kg/m**2]"},     /* 61 */
  {"NCPCP", "Large scale precipitation [kg/m**2]"},     /* 62 */
  {"ACPCP", "Convective precipitation [kg/m**2]"},     /* 63 */
  {"SRWEQ", "Snowfall rate water equivalent [kg/m**2/s]"},     /* 64 */
  {"WEASD", "Water equiv. of accum. snow depth [kg/m**2]"},     /* 65 */
  {"SNOD", "Snow depth [m]"},     /* 66 */
  {"MIXHT", "Mixed layer depth [m]"},     /* 67 */
  {"TTHDP", "Transient thermocline depth [m]"},     /* 68 */
  {"MTHD", "Main thermocline depth [m]"},     /* 69 */
  {"MTHA", "Main thermocline anomaly [m]"},     /* 70 */
  {"TCDC", "Total cloud cover [percent]"},     /* 71 */
  {"CDCON", "Convective cloud cover [percent]"},     /* 72 */
  {"LCDC", "Low level cloud cover [percent]"},     /* 73 */
  {"MCDC", "Mid level cloud cover [percent]"},     /* 74 */
  {"HCDC", "High level cloud cover [percent]"},     /* 75 */
  {"CWAT", "Cloud water [kg/m**2]"},     /* 76 */
  {"var77", "undefined"}, /* 77 */
  {"SNOC", "Convective snow [kg/m**2]"},     /* 78 */
  {"SNOL", "Large scale snow [kg/m**2]"},     /* 79 */
  {"WTMP", "Water temperature [K]"},     /* 80 */
  {"LAND", "Land-sea mask (1=land; 0=sea) [integer]"},     /* 81 */
  {"DSLM", "Deviation of sea level from mean [m]"},     /* 82 */
  {"SFCR", "Surface roughness [m]"},     /* 83 */
  {"ALBDO", "Albedo [percent]"},     /* 84 */
  {"TSOIL", "Soil temperature [K]"},     /* 85 */
  {"SOILM", "Soil moisture content [kg/m**2]"},     /* 86 */
  {"VEG", "Vegetation [percent]"},     /* 87 */
  {"SALTY", "Salinity [kg/kg]"},     /* 88 */
  {"DEN", "Density [kg/m**2]"},     /* 89 */
  {"RUNOF", "Runoff [kg/m**2]"},     /* 90 */
  {"ICEC", "Ice concentration (ice=1; no ice=0) [1/0]"},     /* 91 */
  {"ICETK", "Ice thickness [m]"},     /* 92 */
  {"DICED", "Direction of ice drift [degree]"},     /* 93 */
  {"SICED", "Speed of ice drift [m/s]"},     /* 94 */
  {"UICE", "u of ice drift [m/s]"},     /* 95 */
  {"VICE", "v of ice drift [m/s]"},     /* 96 */
  {"ICEG", "Ice growth [m]"},     /* 97 */
  {"ICED", "Ice divergence [/s]"},     /* 98 */
  {"SNOM", "Snow melt [kg/m**2]"},     /* 99 */
  {"HTSGW", "Sig height of wind waves and swell [m]"},     /* 100 */
  {"WVDIR", "Direction of wind waves [degree]"},     /* 101 */
  {"WVHGT", "Significant height of wind waves [m]"},     /* 102 */
  {"WVPER", "Mean period of wind waves [s]"},     /* 103 */
  {"SWDIR", "Direction of swell waves [degree]"},     /* 104 */
  {"SWELL", "Significant height of swell waves [m]"},     /* 105 */
  {"SWPER", "Mean period of swell waves [s]"},     /* 106 */
  {"DIRPW", "Primary wave direction [degree]"},     /* 107 */
  {"PERPW", "Primary wave mean period [s]"},     /* 108 */
  {"DIRSW", "Secondary wave direction [degree]"},     /* 109 */
  {"PERSW", "Secondary wave mean period [s]"},     /* 110 */
  {"NSWRS", "Net short wave radiation (surface) [W/m**2]"},     /* 111 */
  {"NLWRS", "Net long wave radiation (surface) [W/m**2]"},     /* 112 */
  {"NSWRT", "Net short wave radiation (top) [W/m**2]"},     /* 113 */
  {"NLWRT", "Net long wave radiation (top) [W/m**2]"},     /* 114 */
  {"LWAVR", "Long wave radiation [W/m**2]"},     /* 115 */
  {"SWAVR", "Short wave radiation [W/m**2]"},     /* 116 */
  {"GRAD", "Global radiation [W/m**2]"},     /* 117 */
  {"var118", "undefined"}, /* 118 */
  {"var119", "undefined"}, /* 119 */
  {"var120", "undefined"}, /* 120 */
  {"LHTFL", "Latent heat flux [W/m**2]"},     /* 121 */
  {"SHTFL", "Sensible heat flux [W/m**2]"},     /* 122 */
  {"BLYDP", "Boundary layer dissipation [W/m**2]"},     /* 123 */
  {"UFLX", "Zonal component of momentum flux [N/m**2]"},     /* 124 */
  {"VFLX", "Meridional component of momentum flux [N/m**2]"},     /* 125 */
  {"WMIXE", "Wind mixing energy [J]"},     /* 126 */

  /* ECMWF local table */

  {"at",      /* 127 */ "Atmospheric Tide"},
  {"bdv",     /* 128 */ "Budget Values"},
  {"zg",      /* 129 */ "Geopotential [m2 s-2]"},
  {"ta",      /* 130 */ "Temperature [K]"},
  {"ua",      /* 131 */ "U-component of Wind [ms-1]"},
  {"va",      /* 132 */ "V-component of Wind [ms-1]"},
  {"hus",     /* 133 */ "Specific Humidity [kg/kg]"},
  {"pss",     /* 134 */ "Surface Pressure [Pa]"},
  {"wa",      /* 135 */ "Vertical Velocity [Pa s-1]"},
  {"prwa",    /* 136 */ "preciptable water (vapor+drops+ice) [m]"},
  {"prw",     /* 137 */ "pecipitable water [m]"},
  {"rvort",   /* 138 */ "vorticity [s-1]"},
  {"tso1",    /* 139 */ "soil moisture level 1 [m (H20)]"},
  {"mrso1",   /* 140 */ "soil temperature level 1 [K]"},
  {"snd",     /* 141 */ "Snow Depth [m]"},
  {"prl",     /* 142 */ "Large Scale Precipitation [m]"},
  {"prc",     /* 143 */ "Convective Precipitation [m]"},
  {"prs",     /* 144 */ "Snow Fall"},
  {"bld",     /* 145 */ "Boundary Layer Dissipation [Wm-2]"},
  {"hfss",    /* 146 */ "Surface Flux of Sensible Heat [Wm-2]"},
  {"hfls",    /* 147 */ "Surface Flux of Latent Heat [Wm-2]"},
  {"v148",    /* 148 */ "undefined"},
  {"v149",    /* 149 */ "undefined"},
  {"v150",    /* 150 */ "undefined"},
  {"psl",     /* 151 */ "Mean Sea Level (MSL) Pressure [Pa]    Pa"},
  {"logpsl",  /* 152 */ "Log Surface Pressure"},
  {"v153",    /* 153 */ "undefined"},
  {"v154",    /* 154 */ "undefined"},
  {"div",     /* 155 */ "Divergence [s-1]"},
  {"zg",      /* 156 */ "Height (Geopotential) [m]"},
  {"hur",     /* 157 */ "Relative Humidity [%]"},
  {"pstn",    /* 158 */ "Tendency of Surface Pressure [Pa s-1]"},
  {"v159",    /* 159 */ "undefined"},
  {"v160",    /* 150 */ "undefined"},
  {"v161",    /* 161 */ "undefined"},
  {"v162",    /* 162 */ "undefined"},
  {"v163",    /* 163 */ "undefined"},
  {"clt",     /* 164 */ "cloud cover total [0-1]"},
  {"uas",     /* 165 */ "U-wind at 10 [ms-1]"},
  {"vas",     /* 166 */ "V-wind at 10 [ms-1]"},
  {"tas",     /* 167 */ "Temperature at 2 m [K]"},
  {"tds",     /* 168 */ "Dewpoint at 2 m [K]"},
  {"rsds",    /* 169 */ "Downward SW (sfc) [Wm-2]"},
  {"tso2",    /* 170 */ "soil temperature [K]"},
  {"mrso2",   /* 171 */ "soil wetness level 2 [m (H20)]"},
  {"lsm",     /* 172 */ "land-sea mask [(0,1)]"},
  {"sfr",     /* 173 */ "sfc roughness [m]"},
  {"albs",    /* 174 */ "Albedo [0-1]"},
  {"rlds",    /* 175 */ "Downard LW (sfc) [Wm-2]"},
  {"rss",     /* 176 */ "Net Shortwave Radiation (surface) [Wm-2]"},
  {"rls",     /* 177 */ "Net Longwave Radiation (surface) [Wm-2]"},
  {"rst",     /* 178 */ "Net Shortwave Radiation (toa) [Wm-2"},
  {"rlt",     /* 179 */ "Net Longwave Radiation (toa) [Wm-2]"},
  {"tauu",    /* 180 */ "U-component of Surface Wind Stress [Nm-2]"},
  {"tauv",    /* 181 */ "V-component of Surface Wind Stress [Nm-2]"},
  {"evs",     /* 182 */ "Evaporation [m (H2O)]"},
  {"tso3",    /* 183 */ "soil temp level 3 [m (H2O)]"},
  {"mrso3",   /* 184 */ "soil moisture level 3 [K]"},
  {"clcc",    /* 185 */ "cloud convective [0-1]"},
  {"cll",     /* 186 */ "cloud low [0-1]"},
  {"clm",     /* 187 */ "cloud mid [0-1]"},
  {"clh",     /* 188 */ "cloud high [0-1]"},
  {"v189",    /* 189 */ "undefined"},
  {"orgv",    /* 190 */ "orographic variance"},
  {"orgvew",  /* 191 */ "orographic variance e-w"},
  {"orgvns",  /* 192 */ "orographic variance n-s"},
  {"orgvnwse",/* 193 */ "orographic variance nw-se"},
  {"orgvnesw",/* 194 */ "orographic variance ne-sw"},
  {"gwsv",    /* 195 */ "gravity wave stress n-s"},
  {"gwsu",    /* 196 */ "gravity wave stress e-w"},
  {"gwd",     /* 197 */ "gravity wave diss"},
  {"src",     /* 198 */ "skin resevoir content"},
  {"sfvc",    /* 199 */ "sfc vegetation cover"},
  {"orgvsg",  /* 200 */ "orgographic variance subgrid"},
  {"tasmx",   /* 201 */ "max sfc temp"},
  {"tasmn",   /* 202 */ "min sfc temp"},
  {"v203",    /* 203 */ "undefined"},
  {"praw",    /* 204 */ "precip analysis weights"},
  {"mrro",    /* 205 */ "runoff"},
  {"cvzz",    /* 206 */ "zz variance"},
  {"cvtz",    /* 207 */ "tz covariance"},
  {"cvtt",    /* 208 */ "tt variance"},
  {"cvqz",    /* 209 */ "qz covariance"},
  {"cvqt",    /* 210 */ "qt covariance"},
  {"cvqq",    /* 211 */ "qq variance"},
  {"cvuz",    /* 212 */ "uz covariance"},
  {"cvut",    /* 213 */ "ut covariance"},
  {"cvuq",    /* 214 */ "uq covariance"},
  {"cvuu",    /* 215 */ "uu variance"},
  {"cvvz",    /* 216 */ "vz covariance"},
  {"cvvt",    /* 217 */ "vt covariance"},
  {"cvvq",    /* 218 */ "vq covariance"},
  {"cvvu",    /* 219 */ "vu covariance"},
  {"cvvv",    /* 220 */ "vv variance"},
  {"cvwz",    /* 221 */ "wz covariance"},
  {"cvwt",    /* 222 */ "wt covariance"},
  {"cvwq",    /* 223 */ "wq covariance"},
  {"cvwu",    /* 224 */ "wu covariance"},
  {"cvwv",    /* 225 */ "wv covariance"},
  {"cvww",    /* 226 */ "ww variance"},
  {"cvrr",    /* 227 */ "rh variance"},
  {"pr",      /* 228 */ "total precip"},
  {"tauui",   /* 229 */ "instanteous sfc stress u [Nm-2]"},
  {"tauvi",   /* 230 */ "instanteous sfc stress v [Nm-2]"},
  {"hfssi",   /* 231 */ "instanteous sfc sensible heat flux [Wm-2]"},
  {"hflsi",   /* 232 */ "instanteous sfc latent heat flux"},
  {"husa",    /* 233 */ "apparent sfc humidity"},
  {"logsfr",  /* 234 */ "log sfc roughness"},
  {"tgs",     /* 235 */ "skin temperature [K]"},
  {"tso4",    /* 236 */ "soil temperature level 4 [K]"},
  {"mrso4",   /* 237 */ "soil wetness level 4 [m (H2O)]"},
  {"tgs",     /* 238 */ "t of snow layer [K]"},
  {"prsc",    /* 239 */ "convective snow [m]"},
  {"prsl",    /* 240 */ "large scale snow [m]"},
  {"cllw",    /* 241 */ "cloud liquid water"},
  {"clct",    /* 242 */ "total cloud cover"},
  {"albsf",   /* 243 */ "forecast albedo"},
  {"sfrf",    /* 244 */ "forecast sfc roughness"},
  {"logsfcrf",/* 245 */ "log offorecast sfc roughness"},
  {"wspds",   /* 246 */ "10 m wind speed"},
  {"taum",    /* 247 */ "magnitude of momentum flux"},
  {"v248",    /* 248 */ "undefined"},
  {"gwmf",    /* 249 */ "gravity wave drag momentum flux [Nm-2]"},
  {"v250",    /* 250 */ "undefined"},
  {"v251",    /* 251 */ "undefined"},
  {"v252",    /* 252 */ "undefined"},
  {"v253",    /* 253 */ "undefined"},
  {"v254",    /* 254 */ "undefined"},
  {"v255",    /* 255 */ "undefined"}
};

/* parameter table for ocean modeling branch (OMB) of NCEP */
/* center = 7, subcenter = EMC, parameter table = 128 */

struct ParmTable parm_table_omb[256] = {
 {"var0", "Reserved"},
 {"var1", "Reserved"},
 {"GHz6", "6.6 GHz - K"},
 {"GHz10", "10.7 GHz - K"},
 {"GHz18", "18.0 GHz - K"},
 {"GHz19V", "SSMI 19 GHz, Vertical Polarization - K"},
 {"GHz19H", "SSMI 19 GHz, Horizontal Polarization - K"},
 {"GHz21", "21.0 GHz - K"},
 {"GHz22V", "SSMI 22 GHz, Vertical Polarization - K"},
 {"GHz37V", "SSMI 37 GHz, Vertical Polarization - K"},
 {"GHz37H", "SSMI 37 GHz, Horizontal Polarization - K"},
 {"MSU1", "MSU Ch 1 - 50.30 GHz - K"},
 {"MSU2", "MSU Ch 2 - 53.74 GHz - K"},
 {"MSU3", "MSU Ch 3 - 54.96 GHz - K"},
 {"MSU4", "MSU Ch 4 - 57.95 GHz - K"},
 {"GHz85V", "SSMI 85 GHz, Vertical Polarization - K"},
 {"GHz85H", "SSMI 85 GHz, Horizontal Polarization - K"},
 {"GHz91", "91.65 GHz - K"},
 {"GHz150", "150 GHz - K"},
 {"GHz183pm7", "183 +- 7 GHz - K"},
 {"GHz183pm3", "183 +- 3 GHz - K"},
 {"GHz183pm1", "183 +- 1 GHz - K"},
 {"SSMT1C1", "SSM/T1 - ch 1 - K"},
 {"SSMT1C2", "SSM/T1 - ch 2 - K"},
 {"SSMT1C3", "SSM/T1 - ch 3 - K"},
 {"SSMT1C4", "SSM/T1 - ch 4 - K"},
 {"SSMT1C5", "SSM/T1 - ch 5 - K"},
 {"SSMT1C6", "SSM/T1 - ch 6 - K"},
 {"SSMT1C7", "SSM/T1 - ch 7 - K"},
 {"var29", "Reserved"},
 {"var30", "Reserved"},
 {"var31", "Reserved"},
 {"var32", "Reserved"},
 {"var33", "Reserved"},
 {"var34", "Reserved"},
 {"var35", "Reserved"},
 {"var36", "Reserved"},
 {"var37", "Reserved"},
 {"var38", "Reserved"},
 {"var39", "Reserved"},
 {"var40", "Reserved"},
 {"var41", "Reserved"},
 {"var42", "Reserved"},
 {"var43", "Reserved"},
 {"var44", "Reserved"},
 {"var45", "Reserved"},
 {"var46", "Reserved"},
 {"var47", "Reserved"},
 {"var48", "Reserved"},
 {"var49", "Reserved"},
 {"var50", "Reserved"},
 {"var51", "Reserved"},
 {"var52", "Reserved"},
 {"var53", "Reserved"},
 {"var54", "Reserved"},
 {"var55", "Reserved"},
 {"var56", "Reserved"},
 {"var57", "Reserved"},
 {"var58", "Reserved"},
 {"var59", "Reserved"},
 {"MI14.95", "HIRS/2 ch 1 - 14.95 micron - K"},
 {"MI14.71", "HIRS/2, GOES 14.71 micron - K"},
 {"MI14.49", "HIRS/2 ch 3 - 14.49 micron - K"},
 {"MI14.37", "GOES I-M - 14.37 micron - K"},
 {"MI14.22", "HIRS/2 ch 4 - 14.22 micron - K"},
 {"MI14.06", "GOES I-M - 14.06 micron - K"},
 {"MI13.97", "HIRS/2 ch 5 - 13.97 micron - K"},
 {"MI13.64", "HIRS/2, GOES 13.64 micron - K"},
 {"MI13.37", "GOES I-M - 13.37 micron - K"},
 {"MI13.35", "HIRS/2 ch 7 - 13.35 micron - K"},
 {"MI12.66", "GOES I-M - 12.66 micron - K"},
 {"MI12.02", "GOES I-M - 12.02 micron - K"},
 {"MI12.00", "AVHRR ch 5 - 12.0 micron - K"},
 {"MI11.11", "HIRS/2 ch 8 - 11.11 micron - K"},
 {"MI11.03", "GOES I-M - 11.03 micron - K"},
 {"MI10.80", "AVHRR ch 4 - 10.8 micron - K"},
 {"MI9.71", "HIRS/2, GOES - 9.71 micron - K"},
 {"var77", "Reserved"},
 {"var78", "Reserved"},
 {"var79", "Reserved"},
 {"MI8.16", "HIRS/2 ch 10 - 8.16 micron - K"},
 {"MI7.43", "GOES I-M - 7.43 micron - K"},
 {"MI7.33", "HIRS/2 ch 11 - 7.33 micron - K"},
 {"MI7.02", "GOES I-M - 7.02 micron - K"},
 {"MI6.72", "HIRS/2 ch 12 - 6.72 micron - K"},
 {"MI6.51", "GOES I-M - 6.51 micron - K"},
 {"MI4.57", "HIRS/2, GOES - 4.57 micron - K"},
 {"MI4.52", "HIRS/2, GOES - 4.52 micron - K"},
 {"MI4.46", "HIRS/2 ch 15 - 4.46 micron - K"},
 {"MI4.45", "GOES I-M - 4.45 micron - K"},
 {"MI4.40", "HIRS/2 ch 16 - 4.40 micron - K"},
 {"MI4.24", "HIRS/2 ch 17 - 4.24 micron - K"},
 {"MI4.13", "GOES I-M - 4.13 micron - K"},
 {"MI4.00", "HIRS/2 ch 18 - 4.00 micron - K"},
 {"MI8.16", "GOES I-M - 3.98 micron - K"},
 {"MI8.16", "HIRS/2 Window - 3.76 micron - K"},
 {"MI8.16", "AVHRR, GOES - 3.74 micron - K"},
 {"var97", "Reserved"},
 {"var98", "Reserved"},
 {"var99", "Reserved"},
 {"MI0.91", "AVHRR ch 2 - 0.91 micron - K"},
 {"MI0.696", "GOES I-M - 0.696 micron - K"},
 {"MI0.69", "HIRS/2 Vis - 0.69 micron - K"},
 {"MI0.63", "AVHRR ch 1 - 0.63 micron - K"},
 {"var104", "Reserved"},
 {"var105", "Reserved"},
 {"var106", "Reserved"},
 {"var107", "Reserved"},
 {"var108", "Reserved"},
 {"var109", "Reserved"},
 {"var110", "Reserved"},
 {"var111", "Reserved"},
 {"var112", "Reserved"},
 {"var113", "Reserved"},
 {"var114", "Reserved"},
 {"var115", "Reserved"},
 {"var116", "Reserved"},
 {"var117", "Reserved"},
 {"var118", "Reserved"},
 {"var119", "Reserved"},
 {"var120", "Reserved"},
 {"var121", "Reserved"},
 {"var122", "Reserved"},
 {"var123", "Reserved"},
 {"var124", "Reserved"},
 {"var125", "Reserved"},
 {"var126", "Reserved"},
 {"var127", "Reserved"},
 {"AVDEPTH", "Ocean depth - mean - m"},
 {"DEPTH", "Ocean depth - instantaneous - m"},
 {"ELEV", "Ocean surface elevation relative to geoid - m"},
 {"MXEL24", "Max ocean surface elevation in last 24 hours - m"},
 {"MNEL24", "Min ocean surface elevation in last 24 hours - m"},
 {"var133", "Reserved"},
 {"var134", "Reserved"},
 {"O2", "Oxygen -Mol/kg"},
 {"PO4", "PO4 - Mol/kg"},
 {"NO3", "NO3 - Mol/kg"},
 {"SiO4", "SiO4 - Mol/kg"},
 {"CO2aq", "CO2 (aq) - Mol/kg"},
 {"HCO3", "HCO3 - - Mol/kg"},
 {"CO3", "CO3 -- - Mol/kg"},
 {"TCO2", "TCO2 - Mol/kg"},
 {"TALK", "TALK - Mol/kg"},
 {"var144", "Reserved"},
 {"var145", "Reserved"},
 {"S11", "S11 - 1,1 component of ice stress tensor"},
 {"S12", "S12 - 1,2 component of ice stress tensor"},
 {"S22", "S22 - 2,2 component of ice stress tensor"},
 {"INV1", "T1 - First invariant of stress tensor"},
 {"INV2", "T2 - Second invariant of stress tensor"},
 {"var151", "Reserved"},
 {"var152", "Reserved"},
 {"var153", "Reserved"},
 {"var154", "Reserved"},
 {"WVRGH", "Wave Roughness"},
 {"WVSTRS", "Wave Stresses"},
 {"WHITE", "Whitecap coverage"},
 {"SWDIRWID", "Swell direction width"},
 {"SWFREWID", "Swell frequency width"},
 {"WVAGE", "Wave age"},
 {"PWVAGE", "Physical Wave age"},
 {"var162", "Reserved"},
 {"var163", "Reserved"},
 {"var164", "Reserved"},
 {"LTURB", "Master length scale (turbulence) - m"},
 {"var166", "Reserved"},
 {"var167", "Reserved"},
 {"var168", "Reserved"},
 {"var169", "Reserved"},
 {"AIHFLX", "Net Air-Ice heat flux - W/m^2"},
 {"AOHFLX", "Net Air-Ocean heat flux - W/m^2"},
 {"IOHFLX", "Net Ice-Ocean heat flux - W/m^2"},
 {"IOSFLX", "Net Ice-Ocean salt flux - kg/s"},
 {"var174", "Reserved"},
 {"OMLT", "Ocean Mixed Layer Temperature - K"},
 {"OMLS", "Ocean Mixed Layer Salinity - kg/kg"},
 {"var177", "Reserved"},
 {"var178", "Reserved"},
 {"var179", "Reserved"},
 {"var180", "Reserved"},
 {"var181", "Reserved"},
 {"var182", "Reserved"},
 {"var183", "Reserved"},
 {"var184", "Reserved"},
 {"var185", "Reserved"},
 {"var186", "Reserved"},
 {"var187", "Reserved"},
 {"var188", "Reserved"},
 {"var189", "Reserved"},
 {"var190", "Reserved"},
 {"var191", "Reserved"},
 {"var192", "Reserved"},
 {"var193", "Reserved"},
 {"var194", "Reserved"},
 {"var195", "Reserved"},
 {"var196", "Reserved"},
 {"var197", "Reserved"},
 {"var198", "Reserved"},
 {"var199", "Reserved"},
 {"var200", "Reserved"},
 {"var201", "Reserved"},
 {"var202", "Reserved"},
 {"var203", "Reserved"},
 {"var204", "Reserved"},
 {"var205", "Reserved"},
 {"var206", "Reserved"},
 {"var207", "Reserved"},
 {"var208", "Reserved"},
 {"var209", "Reserved"},
 {"var210", "Reserved"},
 {"var211", "Reserved"},
 {"var212", "Reserved"},
 {"var213", "Reserved"},
 {"var214", "Reserved"},
 {"var215", "Reserved"},
 {"var216", "Reserved"},
 {"var217", "Reserved"},
 {"var218", "Reserved"},
 {"var219", "Reserved"},
 {"var220", "Reserved"},
 {"var221", "Reserved"},
 {"var222", "Reserved"},
 {"var223", "Reserved"},
 {"var224", "Reserved"},
 {"var225", "Reserved"},
 {"var226", "Reserved"},
 {"var227", "Reserved"},
 {"var228", "Reserved"},
 {"var229", "Reserved"},
 {"var230", "Reserved"},
 {"var231", "Reserved"},
 {"var232", "Reserved"},
 {"var233", "Reserved"},
 {"var234", "Reserved"},
 {"var235", "Reserved"},
 {"var236", "Reserved"},
 {"var237", "Reserved"},
 {"var238", "Reserved"},
 {"var239", "Reserved"},
 {"var240", "Reserved"},
 {"var241", "Reserved"},
 {"var242", "Reserved"},
 {"var243", "Reserved"},
 {"var244", "Reserved"},
 {"var245", "Reserved"},
 {"var246", "Reserved"},
 {"var247", "Reserved"},
 {"var248", "Reserved"},
 {"var249", "Reserved"},
 {"var250", "Reserved"},
 {"var251", "Reserved"},
 {"var252", "Reserved"},
 {"var253", "Reserved"},
 {"var254", "Reserved"},
 {"var255", "Reserved"}
};


int GDS_grid(unsigned char *gds, int *nx, int *ny, long int *nxny) 

/*
!C**********************************************************************

!Description:
    Determine number of x & y grid points, and product, from grib GDS.

!Input Parameters:
  unsigned char *gds -- pointer to GDS in grib record.

!Output Parameters:
  int      *nx       -- number of x dim grid points
  int      *ny       -- number of x dim grid points
  long int *nxny     -- nx * ny

!Revision History:
  17-JUN-99 Jim Ray, SSAI -- Added this prolog.
  11/94 - v1.0  Wesley Ebisuzaki

!Team-unique Header:

!References and credits:

!Design Notes:

!END*************************************************************
*/
{

    int i, ix, iy, pl;
    long int isum;

    *nx = ix = GDS_LatLon_nx(gds);
    *ny = iy = GDS_LatLon_ny(gds);
    *nxny = ix * iy;

    /* thin grid */

    if (GDS_Gaussian(gds) || GDS_LatLon(gds)) {
	if (ix == 65535) {
	    *nx = -1;
	    /* reduced grid */
	    isum = 0;
	    pl = GDS_PL(gds);
	    for (i = 0; i < iy; i++) {
		isum += gds[pl+i*2]*256 + gds[pl+i*2+1];
	    }
	    *nxny = isum;
	}
	return 0;
    }
    return 0;
}

#define START -1

static int user_center = 0, user_subcenter = 0, user_ptable = 0;
static enum {filled, not_found, not_checked, no_file, init} status = init;

struct ParmTable parm_table_user[256];

/*
 * sets up user parameter table
 */

int setup_user_table(int center, int subcenter, int ptable) 
/*
!C**********************************************************************
!Description:
  sets up user parameter table
  
!Input Parameters:
  none

!Output Parameters:
  int center    -- 
  int subcenter -- 
  int ptable    -- 
  
!Return Value:
  0 (failure, no file to read from...)
  1 (success)

!Revision History:
  Original version, Wesley Ebisuzaki.
  Derived from wgrib, a portable program to read grib files that were 
  created by the NCEP/NCAR Reanalysis Project.  

!Team-unique Header:

!References and credits:

!Design Notes:

!END*************************************************************
*/
{

    int i, j, c0, c1, c2;
    FILE *input;
    char *filename, line[300];

    if (status == init) {
	for (i = 0; i < 256; i++) {
	    parm_table_user[i].name = parm_table_user[i].comment = NULL;
	}
	status = not_checked;
    }

    if (status == no_file) return 0;

    if ((user_center == -1 || center == user_center) &&
	    (user_subcenter == -1 || subcenter == user_subcenter) &&
	    (user_ptable == -1 || ptable == user_ptable)) {

	if (status == filled) return 1;
	if (status == not_found) return 0;
    }

    /* open gribtab file */

    filename = getenv("GRIBTAB");
    if (filename == NULL) filename = getenv("gribtab");
    if (filename == NULL) filename = "gribtab";

    if ((input = fopen(filename,"r")) == NULL) {
        status = no_file;
        return 0;
    }

    user_center = center;
    user_subcenter = subcenter;
    user_ptable = ptable;

    /* scan for center & subcenter and ptable */
    for (;;) {
        if (fgets(line, 299, input) == NULL) {
	    status = not_found;
            return 0;
        }
	if (atoi(line) != START) continue;
	i = sscanf(line,"%d:%d:%d:%d", &j, &center, &subcenter, &ptable);
        if (i != 4) {
	    /*fprintf(stderr,"illegal gribcap center/subcenter/ptable line: %s\n", line);*/
            continue;
        }
	if ((center == -1 || center == user_center) &&
	    (subcenter == -1 || subcenter == user_subcenter) &&
	    (ptable == -1 || ptable == user_ptable)) break;
    }

    user_center = center;
    user_subcenter = subcenter;
    user_ptable = ptable;

    /* free any used memory */
    if (parm_table_user[i].name != NULL) {
        for (i = 0; i < 256; i++) {
	    free(parm_table_user[i].name);
            free(parm_table_user[i].comment);
        }
    }

    /* read definitions */

    for (;;) {
        if (fgets(line, 299, input) == NULL) break;
	if ((i = atoi(line)) == START) break;
	line[299] = 0;

	/* find the colons and end-of-line */
	for (c0 = 0; line[c0] != ':' && line[c0] != 0; c0++) ;
        if (line[c0] == 0) {
	    /*fprintf(stderr,"illegal gribcap line:%s\n", line);*/
	    continue;
	}
	for (c1 = c0 + 1; line[c1] != ':' && line[c1] != 0; c1++) ;
	c2 = strlen(line);
        if (line[c2-1] == '\n') line[--c2] = '\0';
        if (c2 <= c1) {
	    /*fprintf(stderr,"illegal gribcap line:%s\n", line);*/
	    continue;
	}
	line[c0] = 0;
	line[c1] = 0;

	parm_table_user[i].name = (char *) malloc(c1 - c0);
	parm_table_user[i].comment = (char *) malloc(c2 - c1);
	strcpy(parm_table_user[i].name, line+c0+1);
	strcpy(parm_table_user[i].comment, line+c1+1);
    }

    /* now to fill in undefined blanks */
    for (i = 0; i < 255; i++) {
	if (parm_table_user[i].name == NULL) {
	    parm_table_user[i].name = (char *) malloc(7);
	    sprintf(parm_table_user[i].name, "var%d", i);
	    parm_table_user[i].comment = (char *) malloc(strlen("undefined")+1);
	    strcpy(parm_table_user[i].comment, "undefined");
        }
    }
    status = filled;
    return 1;
}

/* static int msg_count = 0; */

int ASCII_TCA_PDS_date(unsigned char *pds, int v_time, char date[]) 
/*
!C**********************************************************************
!Description:
  returns a string with a date code
  
!Input Parameters:
  unsigned char *pds -- 
  int option         -- flag indicating 2 or 4 digit year
  int v_time         -- flag indicating initial or verification time

!Output Parameters:
  char date[]        -- ASCII, TCA-format date/time
  
!Return Value:
  -1 (failure)
   0 (success)

!Revision History:
 Original version, Wesley Ebisuzaki;
 *
 * prints a string with a date code
 *
 * PDS_date(pds,option, v_time)
 *   options=0  .. 2 digit year
 *   options=1  .. 4 digit year
 *
 *   v_time=0   .. initial time
 *   v_time=1   .. verification time
 *
 * assumption: P1 and P2 are unsigned integers (not clear from doc)
 *
 * v1.2 years that are multiple of 400 are leap years, not 500
 * v1.2.1  make the change to the source code for v1.2
 *
  Derived from wgrib, a portable program to read grib files that were 
  created by the NCEP/NCAR Reanalysis Project.  

!Team-unique Header:

!References and credits:

!Design Notes:
   Adapted from PDS_date() (14-JUL-99, Jim Ray, SSAI).

!END*************************************************************
*/
{

    int year, month, day, hour;

    if (v_time == 0) {
        year = PDS_Year4(pds);
        month = PDS_Month(pds);
        day  = PDS_Day(pds);
        hour = PDS_Hour(pds);
    }
    else {
        if (verf_time(pds, &year, &month, &day, &hour) != 0) {
	    /*if (msg_count++ < 5) fprintf(stderr, "PDS_date: problem\n");*/
	 }
    }

    sprintf(date, "%4.4d-%2.2d-%2.2dT%2.2d:00:00.000000Z", year, month, day, hour);
    
    return 0;
}

#define MINUTE	0
#define HOUR	1
#define DAY	2
#define MONTH	3
#define YEAR	4
#define DECADE	5
#define NORMAL	6
#define CENTURY	7
#define SECOND  254

#define  FEB29   (31+29)
static int monthjday[12] = {
        0,31,59,90,120,151,181,212,243,273,304,334};

static int leap(int year) 
/*
!C**********************************************************************
!Description:
  returns whether year is leap or not
  
!Input Parameters:
  int year         -- four-digit year

!Output Parameters:
  none
  
!Return Value:
   1 (is leap)
   0 (is not leap)

!Revision History:
  Original version, Wesley Ebisuzaki.
  Derived from wgrib, a portable program to read grib files that were 
  created by the NCEP/NCAR Reanalysis Project.  

!Team-unique Header:

!References and credits:

!Design Notes:
 
!END*************************************************************
*/
{
	if (year % 4 != 0) return 0;
	if (year % 100 != 0) return 1;
	return (year % 400 == 0);
}


int add_time(int *year, int *month, int *day, int *hour, int dtime, int unit) 
/*
!C**********************************************************************
!Description:
  adds dtime to current time
  
!Input Parameters:
  int *year         -- four-digit year
  int *month        -- two-digit month
  int *day          -- 
  int *hour         -- 
  int dtime         --
  int unit          --
  
!Output Parameters:
  int *year         -- updated
  int *month        -- updated
  int *day          -- updated
  int *hour         -- updated
  
!Return Value:
   1 (error -- time unit is undefined)
   0 (success)

!Revision History:
  Original version, Wesley Ebisuzaki.
  Derived from wgrib, a portable program to read grib files that were 
  created by the NCEP/NCAR Reanalysis Project.  

!Team-unique Header:

!References and credits:

!Design Notes:
 
!END*************************************************************
*/
{
    int y, m, d, h, jday, i;

    y = *year;
    m = *month;
    d = *day;
    h = *hour;

    if (unit == YEAR) {
	*year = y + dtime;
	return 0;
    }
    if (unit == DECADE) {
	*year =  y + (10 * dtime);
	return 0;
    }
    if (unit == CENTURY) {
	*year = y + (100 * dtime);
	return 0;
    }
    if (unit == NORMAL) {
	*year = y + (30 * dtime);
	return 0;
    }
    if (unit == MONTH) {
	dtime += (m - 1);
	*year = y + (dtime / 12);
	*month = 1 + (dtime % 12);
	return 0;
    }
    if (unit == SECOND) {
	dtime /= 60;
	unit = MINUTE;
    }
    if (unit == MINUTE) {
	dtime /= 60;
	unit = HOUR;
    }
    if (unit == HOUR) {
	dtime += h;
	*hour = dtime % 24;
	dtime = dtime / 24;
        unit = DAY;
    }
    /* this is the hard part */

    if (unit == DAY) {
	/* set m and day to Jan 0, and readjust dtime */
	jday = d + monthjday[m-1];
	if (leap(y) && m > 2) jday++;
        dtime += jday;

        /* 4 year chuncks */
	i = dtime / (4 * 365 + 1);
	if (i) {
	    /* assume century years are leap */
	    y = y + i*4;
	    dtime -= i*(4 * 365 + 1);
            /* see if we have gone past feb 28, 1900, 2000, etc */
	    if ((y - 1) / 100 != (*year-1) / 100) {
	        /* crossed the feb 28, xx00 */
	        /* correct for only one century mark */
		if ((y / 100) % 4 != 0) dtime++;
	    }
	}

	/* one year chunks */
	while (dtime > 365 + leap(y)) {
	    dtime -= (365 + leap(y));
	    y++;
	}

	/* calculate the month and day */

	if (leap(y) && dtime == FEB29) {
	    m = 2;
	    d = 29;
	}
	else {
	    if (leap(y) && dtime > FEB29) dtime--;
	    for (i = 11; monthjday[i] >= dtime; --i);
	    m = i + 1;
	    d = dtime - monthjday[i];
	}
	*year = y;
	*month = m;
	*day = d;
	return 0;
   }
   /*fprintf(stderr,"add_time: undefined time unit %d\n", unit);*/
   return 1;
}


int verf_time(unsigned char *pds, int *year, int *month, int *day, int *hour) 
/*
!C**********************************************************************
!Description:
  gets verification time from pds
  
!Input Parameters:
  unsigned char *pds  -- four-digit year
  
!Output Parameters:
  int *year           -- updated
  int *month          -- updated
  int *day            -- updated
  int *hour           -- updated
  
!Return Value:
   1 (error -- time unit is undefined)
   0 (success)

!Revision History:
 Original version, Wesley Ebisuzaki;
 * verf_time:
 *
 * this routine returns the "verification" time
 * should have behavior similar to gribmap
 *
  Derived from wgrib, a portable program to read grib files that were 
  created by the NCEP/NCAR Reanalysis Project.  

!Team-unique Header:

!References and credits:

!Design Notes:
 
!END*************************************************************
*/
{
    int tr, dtime, unit;

    *year = PDS_Year4(pds);
    *month = PDS_Month(pds);
    *day  = PDS_Day(pds);
    *hour = PDS_Hour(pds);

    /* find time increment */

    dtime = PDS_P1(pds);
    tr = PDS_TimeRange(pds);
    unit = PDS_ForecastTimeUnit(pds);

    if (tr == 10) dtime = PDS_P1(pds) * 256 + PDS_P2(pds);
    if (tr > 1 && tr < 6 ) dtime = PDS_P2(pds);

    if (dtime == 0) return 0;

    return add_time(year, month, day, hour, dtime, unit);
}
