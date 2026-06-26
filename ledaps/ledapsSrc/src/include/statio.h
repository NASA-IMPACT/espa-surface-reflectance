#ifndef STATIO_H
#define STATIO_H

#include <stdio.h>

/*
                          STATIO include file

This is the include file needed by application programs

*/

/*--------------------------- #defines -------------------------------------*/

#define   MAXCLS           256
#define   MAXSIT           256
#define   STNAME           100		/* length of class and site names     */
#define   ECOMP            7		/* complex data type                  */
#define   EDCOMP           8		/* double complex data type           */

#define   IUNKNOWN         3		/* Access on stopen                   */

/*-------------------------- Status return codes -----------------------------*/

/*
          E_SUCC           0		   successful operation
          E_FAIL          -1		   unsuccessful operation
*/
#define   E_GNLIO          1		/* general I/O error                  */
#define   E_NOCL           2		/* class ic not found                 */
#define   E_NOSI           3		/* site is not found                  */
#define   E_NODT           4		/* data-type not found                */
#define   E_NOCLN          5		/* class name not found               */
#define   E_NOSIN          6		/* site name not found                */
#define   E_OPOLD          7		/* old file opened                    */
#define   E_OPNEW          8		/* new file created - no old file     */
#define   E_DTOLD          9		/* old data replaced                  */
#define   E_DTNEW          10		/* new data entered                   */
#define   E_INVCL          11		/* class ic not found/invalid         */
#define   E_INVSI          12		/* site is not found/invalid          */
#define   E_INVDT          13		/* error in data_type.                */
#define   E_INVNM          14		/* class or site name too long        */


/*--------------------------- typedefs -------------------------------------*/

typedef   long    SHAPE_VEC[9];

/*--------------------------- prototypes -----------------------------------*/

void clprompt
(
    FILE **fd,               /* i   stats file descriptor                  */
    long nclass,             /* o   # classes in stats file                */
    long *nuse,              /* i/o # classes selected                     */
    long clind[]             /* o   class indices                          */
);

float determ
(
    double *matrix,
    long norder
);

void egvcvl
(
    long nm,           /* dimension of z                             */
    long n,            /* order of z                                 */
    float d[],         /* diagonal elements of tridiagonal matrix    */
    float e[],         /* subdiagonal elements of tridiagonal matrix */
    float *z[],        /* input tridiagonal matrix                   */
    long ierr          /* return status                              */
                       /* output eigen vector matrix                 */
);

void egvorder
(
    float d[],
    float *z[],
    long m,
    long n
);

void get_panel
(
    FILE **fd,
    long panel,
    long *ic,
    char scnams[MAXCLS][17],    /* Names to write on screeen          */
    long loc[MAXCLS][2],        /* Screen location of output                  */
    long nclass
);

void matdiv
(
    float *a[],
    float b[],
    long m,
    long n
);

void matmul
(
    float *a[],
    float *b[],
    float *c[],
    long m,
    long n
);

void matrans
(
    float *a[],
    float *b[],
    long m,
    long n
);

void ret_wait(void);

long c_stopen 
(
    FILE **fd,          /* file descriptor for stats file     */
    char hos[],         /* host name of stats file            */
    long *access,       /* type of open access (read, write,..*/
    long *nclass        /* number of classes in file          */
);

int c_stclse 
(
    FILE **fd       /* file descriptor for stats file     */
);

int c_stget 
(
    FILE **fd,          /* file descriptor for stats file     */
    char data_type[],   /* name of data to be read            */
    long *ic,           /* class index of node to read from   */
    long *is,           /* site index of node to read from    */
    SHAPE_VEC shape,    /* shape of data read in              */
    unsigned char *buf  /* the actual data read in            */
);

int c_stput 
(
    FILE **fd,          /* file descriptor for stats file     */
    char data_type[],   /* name of data to write              */
    long *ic,           /* class index to write data to       */
    long *is,           /* site index to write data to        */
    SHAPE_VEC shape,    /* shape of data to be written        */
    unsigned char *buf  /* the actual data to be written      */
);

int c_stfdcl 
(
    FILE **fd,          /* file descriptor for stats file     */
    char class_name[],  /* class name to find                 */
    long *ic            /* class index of class               */
);

int c_stfdsi 
(
    FILE **fd,          /* file descriptor for stats file     */
    char class_name[],  /* class name to be searched          */
    char site_name[],   /* site name to be searched           */
    long *ic,           /* class index found                  */
    long *is            /* site index found                   */
);

int c_stdldt 
(
    FILE **fd,          /* file descriptor for stats file     */
    char data_type[],   /* name of data to be deleted         */
    long *ic,           /* class index                        */
    long *is            /* site index                         */
);

int c_stdlnd 
(
    FILE **fd,      /* file descriptor for stats file     */
    long *ic,       /* class index                        */
    long *is        /* site index                         */
);

int c_stnode 
(
    FILE **fd,          /* file descriptor for stats file     */
    long *ic,           /* class index                        */
    long *is,           /* site index                         */
    long *id,           /* data item index within node        */
    char data_type[],   /* name of the next data              */
    SHAPE_VEC shape     /* shape of data                      */
);

void tridia
(
    long nm,                    /* row dimension of matrices a and z  */
    long n,                     /* order of a                         */
    float *a[],                 /* input matrix to be transformed     */
    float d[],                  /* diagonal elements of tridiagonal   */
    float e[],                  /* subdiagonal elements of tridiognal */
    float *z[]                  /* orthognal transformation matrix    */
);

#endif
