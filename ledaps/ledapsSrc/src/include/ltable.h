#ifndef LTABLE_H
#define LTABLE_H

#define FATAL	 -1		/* fatal error encountered */
#define NONFATAL  1		/* nonfatal error encountered */

#define TAB_LT_OPEN		-1
#define TAB_LT_ACCESS		-2
#define TAB_LT_READ		-3
#define TAB_LT_WRITE		-4
#define TAB_LT_CLOSE		-5
#define TAB_LT_MATRIX_SIZE	-6
#define TAB_LT_FORMAT		-7
#define TAB_LT_DTYPE		-8
#define TAB_LT_EOR		-9
#define TAB_LT_NOEOR		-10
#define TAB_FORTRAN_MAXFLD	-11
#define TAB_ALLOC		-12

#define TAB_LT_SEEK		-13
#define TAB_LT_FTYPE		-14
#define TAB_LT_SUBFIL		-15
#define TAB_LT_PUTVEC		-16
#define TAB_LT_GETDATE		-17
#define TAB_LT_GETTIME		-18

#ifndef MIN
#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#endif

#ifndef FALSE
#define FALSE     0		/* boolean value for FALSE */
#endif
#ifndef TRUE
#define TRUE	  1		/* boolean value for TRUE  */
#endif

#define MAX_SUB  20		/* maximum number of associated LT files    */
#define TAE_VAL  50		/* max number of TAE values for a parameter */

#define FIELD_SEP1 ','		/* logical record field separator 1 */
#define FIELD_SEP2 '\n'		/* logical record field separator 2 */
#define RECORD_SEP ';'		/* logical record separator */
#define CHAR_DEL   '\"'		/* character delimiter */

#define FSEP1	0		/* field separator 1 for FORTRAN programs */
#define FSEP2	1		/* field separator 2 for FORTRAN programs */
#define RSEP	2		/* record separator for FORTRAN programs  */

#define LIN_SIZ         80      /* number characters without a new line */
#define MAX_FIELD_LEN 1024      /* max number of characters in a record field */
#define TAB_SIZ       8192	/* labeled table physical record size */

#define C	0		/* character data type */
#define I1	1		/* 8-bit integer data type */
#define I2	2		/* 16-bit integer data type */
#define I4	3		/* 32-bit integer data type */
#define R4	4		/* 32-bit real data type */
#define R8	5		/* 64-bit real data type */
#define LTDATE  6		/* Date stored as char string mm/dd/yyyy */
#define LTTIME	7		/* Time stored as char string hh:mm:ss */
#define NDTYPE  8		/* total number of data types */

#ifndef OPEN_TAB_C_
#ifdef sun
/* our current sun platform doesn't like the "const" keyword here for some
   reason - merge back with the rest when the sun is updated */
extern long LT_DTYPE_SIZE[NDTYPE];
extern char DATA_TYPE[NDTYPE][7];
#else
extern const long LT_DTYPE_SIZE[NDTYPE];
extern const char DATA_TYPE[NDTYPE][7];
#endif
#endif

#include "desc.h"

struct LT_DATE
   {
   long month;
   long day;
   long year;
   };

struct LT_TIME
   {
   long hour;
   long minute;
   long second;
   };

struct MATRIX			/* Structrue defining the number of rows    */
   {				/*   and columns of a matrix                */
   long nrow;
   long ncol;
   };

struct VECTOR			/* Recursive structure defining the         */
   {				/*   label vectors of the Labeled Table.    */
   char *label;			/* label name                               */
   long dtype;			/* data type of the vector                  */
   char *desc;			/* label description                        */
   union			/* define a union which depends on the data */
      {				/*   type as to which pointer is currently  */
      char   *c;		/*   valid.  There is a separate pointer    */
      char   *i1;		/*   for each of the valid data types.      */
      short  *i2;
      long   *i4;
      float  *r4;
      double *r8;
      } ptr;
   struct MATRIX lt_size;	/* number of rows & columns in labeled table*/
   struct MATRIX size;		/* actual number of rows & columns in data  */
   long null;			/* flag indicating null fields within data  */
   char sep;			/* type of separator                        */
   long flag;			/* boolean flag used for APPEND access      */
   struct VECTOR *next;		/* pointer to next vector                   */
   };

struct TAB_DEF
   {
   char hname[CMLEN];		/* host name of the Labeled Table          */
   char *subfile[MAX_SUB];	/* list of files associated to this LT      */
   long nsubs;			/* number of files associated to this LT    */
   long fd;			/* file descriptor to label table           */
   long access;			/* access the Labeled Table was opened for */
   char *ftype;			/* file type of the labeled table           */
   char *fdesc;			/* file description of the labeled table    */
   struct DESC buf;		/* data buffer                              */
   long ncol;			/* number of columns in the LT              */
   struct VECTOR *vector;	/* column vector                            */
   long nbytes;			/* number of bytes in data buffer           */
   long cur;			/* current byte offset into the data buffer */
   long currec;			/* current logical record being processed   */
   long label_flag;		/* flag indicating data has been written    */
   };


void close_tab
(
    struct TAB_DEF *tab /* Structure defining the Labeled Table files */
);

int cltab
(
    long *fd            /* File descriptor (created by "optab") (input)      */
);

void free_vec
(
    struct VECTOR *vec  /* Labeled Table vector                       */
);

void tab_error
(
    struct TAB_DEF *tab,
    long err_code,
    long severity,
    char msg[],
    char key[]
);

void get_lt_date
(
    struct VECTOR *vec,  /* Defines the LT record element              */
    struct LT_DATE *date /* Date read from the LT I/O buffer           */
);

void get_lt_time
(
    struct VECTOR *vec, /* Structure defining the LT record element   */
    struct LT_TIME *tim /* Time read from the LT I/O buffer           */
);

void get_matrix
(
    struct VECTOR *vec,   /* Structure defining the record to copy      */
    long odtype,          /* Data type of the output matrix             */
    unsigned char *matrix /* Output matrix (any valid matrix data type) */
);

int get_record
(
    struct TAB_DEF *tab /* Labeled Table to be processed              */
);

void get_i1_field
(
    struct TAB_DEF *tab,  /* Labeled Table file to be processed         */
    char *val,            /* Byte value read from Labeled Table I/O buff*/
    long *null,           /* Null flag:                                 */
    long *eor             /* End of record flag:
                                TRUE  - field terminated by record sep.
                                FALSE - field terminated by field sep.  */
);

void get_i1_matrix
(
    struct TAB_DEF *tab,    /* Labeled Table file to be processed         */
    struct MATRIX *lt_size, /* number of rows & columns of LT matrix      */
    struct MATRIX *out,     /* actual num of rows & columns of matrix read*/
    char matrix[],          /* byte matrix read from LT I/O buffer        */
    long *null,             /* Null flag:                                 */
    long *eor               /* End of record flag:                        */
                            /*  TRUE  - matrix terminated by record sep   */
                            /*  FALSE - matrix terminated by field sep    */
);

void get_i2_field
(
    struct TAB_DEF *tab,    /* Labeled table to be processed              */
    short *val,             /* 16-bit value read from LT I/O buffer       */
    long *null,             /* Null flag:                                 */
    long *eor               /* End of record flag:                        */
                            /*   TRUE  - field terminated by record sep   */
                            /*   FALSE - field terminated by field sep    */
);

void get_i2_matrix
(
    struct TAB_DEF *tab,    /* Labeled Table file to be processed         */
    struct MATRIX *lt_size, /* number of rows & columns of LT matrix      */
    struct MATRIX *out,     /* actual num of rows & columns of matrix read*/
    short matrix[],         /* 16-bit integer matrix read from LT I/O buff*/
    long *null,             /* Null flag:                                 */
    long *eor               /* End of record flag:                        */
);

void get_i4_field
(
    struct TAB_DEF *tab,    /* Labeled Table to be processed.             */
    long *val,              /* integer value read from LT I/O buffer      */
    long *null,             /* Null flag:                                 */
    long *eor               /* End of record flag:                        */
                            /*   TRUE  - field terminated by record sep   */
                            /*   FALSE - field terminated by field sep    */
);

void get_i4_matrix
(
    struct TAB_DEF *tab,    /* Labeled Table file to be processed         */
    struct MATRIX *lt_size, /* Number of rows & columns of LT matrix      */
    struct MATRIX *out,     /* Actual num of rows & columns of matrix read*/
    long matrix[],          /* integer matrix read from LT I/O buffer     */
    long *null,             /* Null flag:                                 */
    long *eor               /* End of record flag:                        */
                            /*   TRUE  - matrix terminated by record sep. */
                            /*   FALSE - matrix terminated by field sep.  */
);

void get_r4_field
(
    struct TAB_DEF *tab,    /* Labeled Table file to be processed         */
    float *val,             /* Real value read from LT I/O buffer         */
    long *null,             /* Null flag:                                 */
    long *eor               /* End of record flag:                        */
                            /*   TRUE  - field terminated by record sep   */
                            /*   FALSE - field terminated by field sep    */
);

void get_r4_matrix
(
    struct TAB_DEF *tab,    /* Labeled Table file to be processed.        */
    struct MATRIX *lt_size, /* number of rows & columns of LT matrix      */
    struct MATRIX *out,     /* actual num of rows & columns of matrix read*/
    float matrix[],         /* Real matrix read from LT I/O buffer        */
    long *null,             /* Null flag:                                 */
    long *eor               /* End of record flag:                        */
                            /*   TRUE  - matrix terminated by record sep. */
                            /*   FALSE - matrix terminated by field sep.  */
);

void get_r8_field
(
    struct TAB_DEF *tab,    /* Labeled Table to be processed.             */
    double *val,            /* integer value read from LT I/O buffer      */
    long *null,             /* Null flag:                                 */
    long *eor               /* End of record flag:                        */
                            /*   TRUE  - field terminated by record sep   */
                            /*   FALSE - field terminated by field sep    */
);

void get_r8_matrix
(
    struct TAB_DEF *tab,    /* Labeled Table file to be processed         */
    struct MATRIX *lt_size, /* number of rows & columns of LT matrix      */
    struct MATRIX *out,     /* actual num of rows & columns of matrix read*/
    double matrix[],        /* Double matrix read from LT I/O buffer      */
    long *null,             /* Null flag:                                 */
    long *eor               /* End of record flag:                        */
                            /*   TRUE  - matrix terminated by record sep. */
                            /*   FALSE - matrix terminated by field sep.  */
);

void get_vector
(
    struct TAB_DEF *tab /* struct of LT file to be processed          */
);

struct VECTOR *get_field_ptr
(
    struct TAB_DEF *tab,
    char *name
);

char *get_string_field
(
    struct TAB_DEF *tab,
    long *eor
);

char *make_lt_name
(
    char innam[],
    long nmr
);

void open_tab
(
    struct TAB_DEF *tab,    /* struct of Labeled Table     */
    char *fname,            /* Host name of the LT file to be opened      */
    long access             /* type of access to open the file with       */
                            /*   0 - read    */
                            /*   1 - write     */
                            /*   2 - append     */
);

void put_lt_date
(
    struct VECTOR *vec,
    struct LT_DATE *date
);

void put_lt_time
(
    struct VECTOR *vec,
    struct LT_TIME *tim
);

void put_matrix
(
    struct VECTOR *vec,         /*  VECTOR structuref                     */
    long idtype,                /*  input data type                       */
    struct MATRIX *size,        /*  Matrix structure                      */
    unsigned char *matrix       /*  Output buffer                         */
);

void put_record
(
    struct TAB_DEF *tab
);

void put_i1_matrix
(
    struct TAB_DEF *tab,        /* structure of the LT file being processed   */
    struct MATRIX *in,          /* Number of rows & columns of the matrix     */
    struct MATRIX *lt_size,     /* actual                                     */
    char matrix[],              /* 8-bit integer matrix to be written         */
    int sep                     /* field or reord separator to be written to the
                                   I/O buffer after the real matrix           */
);

void put_i2_matrix
(
    struct TAB_DEF *tab,        /*  structure defining the LT file to be 
                                    processed*/
    struct MATRIX *in,          /*  number of rows & columns of input matrix  */
    struct MATRIX *lt_size,     /*  number of rows & columns of LT matrix     */
    short matrix[],             /*8-bit integer matrix to be written          */
    int sep                     /* field or record separator to be written to */
                                /*the I/O buffer after the real matrix        */
);

void put_i4_field
(
    struct TAB_DEF *tab,        /*  structure defining the LT file to be      */
    int val,
    int sep                     /*  field or record separator to be written   */
                                /*the I/O buffer after the real field         */
);

void put_i4_matrix
(
    struct TAB_DEF *tab,
    struct MATRIX *in,
    struct MATRIX *lt_size,     /*  number of rows & columns of LT matrix    */
    long matrix[],              /*  32-bit integer matrix to be written      */
    int sep                     /*  field or record separator to be written  */
                                /*the I/O buffer after the real matrix        */
);

void put_r4_field
(
    struct TAB_DEF *tab,        /*structure defining LT file to be processed. */
    double val,                 /*real value to be written into the I/O buffer*/
    int sep                     /*field or record separator to be written to  */
                                /*the I/O buffer after the real field         */
);

void put_r4_matrix
(
    struct TAB_DEF *tab,        /*structure defining LT file to be processed  */
    struct MATRIX *in,          /*number of rows & columns of input matrix    */
    struct MATRIX *lt_size,     /* number of rows & columns of LT matrix      */
    float matrix[],             /*real matrix to be written to the I/O buffer */
    int sep                     /*field or record separator to be written to  */
                                /*the I/O buffer after the real matrix        */
);

void put_r8_field
(
    struct TAB_DEF *tab,        /*structure defining LT file to be processed  */
    double val,                 /*double value to be written to the I/O buffer*/
    int sep                     /*field or record separator to be written to  */
                                /*the I/O buffer after the real field         */
);

void put_r8_matrix
(
    struct TAB_DEF *tab,        /*structure defining LT file to be processed. */
    struct MATRIX *in,          /*number of rows & columns of input matrix    */
    struct MATRIX *lt_size,     /*number of rows & columns of LT matrix       */
    double matrix[],            /*double matrix written into the I/O buffer   */
    int sep                     /*field or record separator to be written to  */
                                /*the I/O buffer after the real matrix        */
);

void put_vector
(
    struct TAB_DEF *tab,        /*structure defining LT file to be processed. */
    char label[],               /*label name of the vector                    */
    long dtype,                 /*data type of the label vector               */
    char desc[],                /*description of the label vector             */
    long nrow,                  /*number of rows                              */
    long ncol,                  /*number of columns                           */
    int sep                     /*separator (FIELD_SEP1,FIELD_SEP2,RECORD_SEP)*/
);

void put_null_field
(
    struct TAB_DEF *tab,
    int sep
);

void put_null_matrix
(
    struct TAB_DEF *tab,
    struct MATRIX *size,
    int sep
);

void put_string_field
(
    struct TAB_DEF *tab,        /*structure defining LT file to be processed. */
    char *val,                  /*character string  written into I/O buffer   */
    int sep,                    /*field or record separator to be written to  */
    long flag                   /*flag specifying whether or not the string is*/
);

void read_tab
(
    struct TAB_DEF *tab /*structure defining LT file to be processed */
);

int rdtab
(
    long *fd,                   /* File descriptor array (Created by "optab") */
    long *nbytes,               /* Number of bytes to read      (input)       */
    struct DESC *buffer         /* Buffer to read data into     (output)      */
);

void upcase
(
    char *buf   /*character buffer to be converted to upper case      */
);

void write_tab
(
    struct TAB_DEF *tab /*structure defining LT file to be processed. */
);

int wrtab
(
    long *fd,                   /* File descriptor array        (input)      */
    long *nbytes,               /* Number of bytes to write     (input)      */
    struct DESC *buffer         /* Buffer to write data from    (input)      */
);
#endif
