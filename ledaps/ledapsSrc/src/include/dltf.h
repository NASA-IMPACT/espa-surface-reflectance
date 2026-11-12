#ifndef DLTF_H
#define DLTF_H

#include "desc.h"

/*****************************************************************************
 *				DLTF Include File			     *
 *									     *
 *  This C include file defines the support and application structures       *
 *  necessary for use of the Display and Look-up Table File.		     *
 *									     *
 *  By D. Akkerman, USGS/EROS Data Center, 8/88				     *
 *									     *
 *****************************************************************************/

/*  The ASDLTF structure is used by the application programmer.		*/

struct ASDLTF
	{
		long 	numbnd,
			band[256],
			window[4],
			izoom[2],
			gzoom[2],
			ishift[2],
			gshift[2];
	};

/*  The CDLTF structure is used by DLTF support routines and forms the
    character portion of the label services file.			*/

struct CDLTF
	{
		char	window[4][7],
			izoom[2][8],
			gzoom[2][8],
			ishift[2][8],
			gshift[2][8],
			descr[81],
			numbnd[4],
			bands[256][4];
	};

/* The SZDLTF define is the total number of bytes in the CDLTF structure 
   excluding the bands field which is variable.  This is used in the 
   calculation of the total number of bytes in the character portion of 
   a DLTF.  The MXDLTF define is the maximum possible size of the CDLTF
   structure.  This is used to call Label Services routines.		*/

#define SZDLTF 177
#define MXDLTF 1201

/* The RECSIZ2 define is the total number of bytes in the character portion
   of the second type of DLTF record.  The number is derived from the space
   used by the fields of the record.  See the Display and Lookup Table File
   format for details.							*/

#define RECSIZ2 38

/* The DESLEN define is the length of the description field of the DLTF
   record.								*/

#define DESLEN 81

/* The KEYLEN define is the length of the key field of the DLTF record.
   The value of this field corresponds to the value of (DMSLEN + 1) and
   should not be changed because of the possible collision it could cause
   when a DLTF is used by EDC/DMS software.			*/

#define KEYLEN 13

/* The KEYTL define is the true length of the key field used by Label
   Services.  It is needed by DLTF routines in instances where the
   LUT number is appended to the user specified entry name to create
   the name of the second type of DLTF record (LUT name).  The KEYTL
   length is the same as the KEYL length found in the vll.h include
   file and should not be changed unless Label Services is also 
   changed.								*/

#define KEYTL 16

/* The DTYPL define is the length of the data type field of the DLTF
   record.  It is the same length as the TYPL define found in the
   v11.h include file and should not be changed unless Label Services
   is also changed.							*/

#define DTYPL 3


int ckdltf
(
    struct DESC *lasfil,
    struct DESC *key,
    long *fd
);

int cldltf
(
    long *fd,
    struct DESC *lasfil
);

int r1dltf
(
    long *fd,
    struct DESC *host,
    struct DESC *key,
    struct DESC *descr,
    long *asdltf,
    long *dtlen,
    struct DESC *dtype
);

int r2dltf
(
    long *fd,
    struct DESC *key,
    long *band,
    long *startpt,
    long *endpt,
    long *dlen,
    unsigned char *lutdata,
    double *slope,
    double *offset
);

int w1dltf
(
    long *fd,
    struct DESC *key,
    struct DESC *descr,
    struct DESC *type,
    long *asdltf
);

int w2dltf
(
    long *fd,
    struct DESC *type,
    struct DESC *key,
    double *slope,
    double *offset,
    long *startpt,
    long *endpt,
    long *lutnum,
    unsigned char *lutdata
);

int c_ckdltf
(
    char *host,         /* Host file name of the DLTF.  The return
                           string has been "squeezed" for FORTRAN.
                           (I/O)                                */
    char *key,          /* DLTF entry name. The return string
                           contains the key in upper case, and
                           has been "squeezed" for FORTRAN.(I/O)*/
    long *fd            /* File descriptor of the DLTF. Assigned
                           only if the duplicate entry check is
                           unsuccessful. (O)                    */
);

int c_cldltf
(
    long *fd,           /* File descriptor for the DLTF. (I)            */
    char *host          /* Host file name of the DLTF. (I)              */
);

int c_r1dltf
(
    long *fd,           /* File descriptor for the DLTF. (O)            */
    char *host,         /* Host file name of the DLTF. (I)              */
    char *key,          /* DLTF entry name; converted to upper case. 
                           (I/O)                                        */
    char *descr,        /* String containing a user specified description
                           of the DLTF entry. (O)                       */
    struct ASDLTF *asdltf, /* Structure containing DLTF field information.
                              (O)                                       */
    long *dtlen,        /* Total number of bytes required for the data 
                           portion of the second DLTF record type. (O)  */
    char *dtype         /* Data type of the DLTF. (O)                   */
);

int c_r2dltf
(
    long *fd,           /* File descriptor for the DLTF. (I)            */
    char *key,          /* DLTF entry name. (I)                         */
    long *band,         /* Band number of member to be read. (I)        */
    long *startpt,      /* Starting point of the LUT. (O)               */
    long *endpt,        /* Ending point of the LUT. (O)                 */
    long *dlen,         /* Length of the data portion of the second type
                           of DLTF record. (I)                          */
    unsigned char *lutdata, /* Buffer for the data portion of the DLTF
                               entry. (O)                               */
    double *slope,      /* Slope of the LUT. (O)                        */
    double *offset      /* Offset of the LUT. (O)                       */
);

int c_w1dltf
(
    long *fd,           /* File descriptor for the DLTF. (I)            */
    char *key,          /* DLTF entry name. (I)                         */
    char *descr,        /* String containing a user specified description
                           of the DLTF entry. (I)                       */
    char *type,         /* String containing the data type of the LUT(s).
                           (I)                                          */
    struct ASDLTF *asdltf /* Structure containing DLTF field information.
                           (I)                                          */
);

int c_w2dltf
(
    long *fd,           /* File descriptor for the DLTF. (I)            */
    char *type,         /* String containing the data type of the LUT(s);
                           either BYTE or I2. (I)                       */
    char *key,          /* DLTF entry name. (I)                         */
    double *slope,      /* The slope of the image band that corresponds to
                           the LUT band. (I)                            */
    double *offset,     /* The offset of the image band that corresponds to
                           the LUT band. (I)                            */
    long *startpt,      /* The value of the input image that corresponds to
                           the first entry in the LUT data. (I)         */
    long *endpt,        /* The value of the input image that corresponds to
                           the last entry in the LUT data. (I)          */
    long *lutnum,       /* LUT number that corresponds to the image band. 
                           (I)                                          */
    unsigned char *lutdata  /* The LUT to be written to the DLTF. (I)       */
);
#endif
