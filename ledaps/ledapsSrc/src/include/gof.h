#ifndef GOF_H
#define GOF_H

#include "vll.h"

#ifndef FUNCTION
#define FUNCTION
#endif

#define DELREC "DELETED" 		/* Label of a deleted record.	      */
#define GOF_SEP '_'                     /* Used in creation of GOF host file
					   names.                             */
#define LENTYP	(LENL+TYPL+KEYL)	/* Length of type, length, and label. */

#define HDRLEN  538                     /* Header length.                     */
#define MAXATT   35                     /* Maximum number of attributes.      */
#define ATHAT    "@^"   		/* Wild card indicator for GOF
                        		   attribute search                   */
#define ATTIN    22     		/* TAE parameter attribute input
                       			   string length.                     */
#define ATTNAM   9			/* Attribute name text length.	      */
#define ATTYP    4			/* Attribute type text length.	      */
#define ATTVAL   13			/* Attribute value text length.	      */
#define ANTLEN   509			/* Annotation text length.	      */
#define CHARLEN  455			/* Length of character portion of 
					   GOF records.			      */
#define MAXGOF   4   			/* Maximum number of GOF file types   */
#define XYLEN	 8			/* Point pair length.		      */
#define LABLEN   16     		/* Label text length                  */

/* 
    Define constants for read/write types.
 */
#define HDR		0		/* Header record		      */
#define PNT		1		/* Point record			      */
#define HPOLIN		2		/* Line or polygon header record      */
#define FPOLIN		3		/* Full polygon or line record	      */
#define XYPAIRS		4		/* XY point pairs of a polygon or line*/
#define ANOT		5		/* Full annotation record	      */

#define HDRKEY "GOFHDR"                 /* Key for the GOF header record      */

/* 
    Define structures for reads and writes.
*/
struct	HEADER				/* Header record structure 	      */
    {
    unsigned char	mask[MAXATT];
    char		attrb[MAXATT][ATTNAM];
    char		attyp[MAXATT][ATTYP];
    int			nrec;
    };

struct	POINT				/* Point record structure 	      */
    {
    float	x;
    float	y;
    char	label[KEYL];
    char 	attrb[MAXATT][ATTVAL];
    };

struct	POLIN				/* Polygon and line record structure  */
    {
    float	*x;
    float	*y;
    int		numpts;
    char	label[KEYL];
    char	attrb[MAXATT][ATTVAL];
    };

struct	ANNOT				/* Annotation record structure 	      */
    {
    float	x;
    float	y;
    float	angle;
    int		size;
    char	label[KEYL];
    char 	attrb[MAXATT][ATTVAL];
    char	antext[ANTLEN];
    };

union GOF				/* Union of all types of record 
					   structures 			      */
    {
    struct HEADER hdrec;
    struct POINT  ptrec;
    struct POLIN  plrec;
    struct ANNOT  antrec;
    };

struct RWHDR				/* Read/Write header record structure */
    {
    char 		length[LENL];
    char 		type[TYPL];
    char		key[KEYL];
    char		rnum[16];
    char		attrb[MAXATT][ATTNAM];
    char 		attyp[MAXATT][ATTYP];
    unsigned char 	mask[MAXATT];
    };

struct PPLGOF				/* Read/Write structure for lines,    */
   {					/* points, and polygon records.       */
   char 	length[LENL];
   char 	type[TYPL];
   char		label[KEYL];
   char		attrb[MAXATT][ATTVAL];
   float	x;
   float	y;
   float	*xpoints;
   float	*ypoints;
   };

struct ANOT1				/* Read/Write structure for character */
   {					/* portion of the annotation record   */
   char 	length[LENL];
   char 	type[TYPL];
   char		label[KEYL];
   char		attrb[MAXATT][ATTVAL];
   char		anotstr[ANTLEN];
   };

struct ANOT2				/* Read/Write structure for data      */
   {					/* portion of the annotation record   */
   float 	x;
   float 	y;
   float	size;
   float	angle;
   };

int c_attchk
(
    struct HEADER *hdr,         /* Header information from the GOF. (I)       */
    char inatt[][ATTIN],        /* Input attribute string. NAME:VALUE
                                   NOTE: Any of the MAXATT input attributes
                                   that do not have a specified attribute
                                   must be assigned a null string. (I)        */
    char outatt[][ATTVAL],      /* Output attribute VALUES in the order
                                   defined by the header record from the 
                                   GOF. (O)                                   */
    int index[MAXATT],          /* Indices of attributes for searching. (O)   */
    int *search                 /* Flag indicating checking by a search func.
                                   (I)                                        */
);

int c_gofcls
(
    int *fd             /* File descriptor of the file to be closed. (I) */
);

int c_gofwr
(
    int *fd,            /* File descriptor of the GOF to read (I).      */
    int *nbyte,         /* Number of bytes to next record (O).          */
    int *type,          /* Type of record being read (I).               */
    union GOF *gofptr,  /* Pointer to the appropriate structure (I).    */
    int *addrec         /* TRUE - Adding records to the GOF     
                           FALSE - Updating existing records (I).       */
);

int c_gofopn 
(
    int  *fd,           /* File descriptor for the opened file. (O) */
    char *hostname,     /* Host file name of the image whose GOF is to
                           be opened. (I)  On output, the host name of
                           the GOF. (O)*/ 
    char *goflasnm,     /* LAS name of GOF opened. (O)*/
    char *type,         /* Type of GOF being opened. (I)*/
    int  *access,       /* Access mode (IREAD, IWRITE, IUPDATE). (I)*/
    int  *cflag         /* Create flag.  A value of FALSE means that the
                           routine did not have to create the GOF.  A 
                           value of TRUE means that the routine created a 
                           GOF. (O)*/
);

int c_gofrd 
(
    int *fd,                /* File descriptor of the GOF to read (I). */
    int *obyte,             /* Number of bytes read (O).*/
    int *nbyte,             /* Number of bytes to next record (O).*/
    int *type,              /* Type of record being read (I).*/
    union GOF *gofptr       /* Pointer to the proper structure (O). */
);

#endif
