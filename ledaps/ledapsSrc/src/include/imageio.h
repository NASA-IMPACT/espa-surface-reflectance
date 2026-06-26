#ifndef IMAGEIO_H
#define IMAGEIO_H

#include <sys/types.h>  /* include for the system specific seek types */
#include <unistd.h>     /* include for the system specific seek types */

#define MAXBND 256   /* Maximum number of bands to be processed            */
#define MAXNS  50000 /* Maximum number of samples per line to be processed */
#define MAXNL  50000 /* Maximum number of lines to be processed	           */

/* Arrays used for specifying windows contain four fields, sl, ss, nl, ns.
   Here are the definitions of those array indices.
---------------------------------------------------------------------------*/
#define SL 0
#define SS 1
#define NL 2
#define NS 3
#define WINDOW 4

#ifdef aiws
#define MAXIMG 20    /* Maximum number of images to be processed           */
#endif
#ifdef gould
#define MAXIMG 54    /* Maximum number of images to be processed           */
#endif
#ifdef _IBMR2
#define MAXIMG 54    /* Maximum number of images to be processed           */
#endif
#ifdef vms
#define MAXIMG 56    /* Maximum number of images to be processed           */
#endif
#ifdef sun
#define MAXIMG 54    /* Maximum number of images to be processed           */
#endif
#ifdef sgi
#define MAXIMG 54    /* Maximum number of images to be processed           */
#endif
#ifdef DGUX
#define MAXIMG 54    /* Maximum number of images to be processed           */
#endif
#ifdef linux
#define MAXIMG 54    /* Maximum number of images to be processed           */
#endif

/*  Maximum buffer size to be used for most memory allocations
--------------------------------------------------------------*/
#ifdef vms
#define MAXBUFS 65535  	 	  /* 65535 is the largest i/o vms can do */
#else
#define MAXBUFS 1024 * 1024 * 2.0 /* 2.0 megabyte of memory -estimated guess
  			             of the best amount of memory to allocate */
#endif

/* method types 
---------------*/
#define STEP 1
#define RANDOM 2

#define NLABELS 0
#define LABSIZ 512
#define ORG I_CS

/* fdesc flags 
--------------*/
#define CONVERSION 1
#define EDONE 2
#define MOVEDATA 4


/* BLOCK IMAGE PROCESSING CONSTANTS
-----------------------------------*/

/* boundary fill techniques
---------------------------*/
#define NOFIL 0
#define MIRR  1

/* end of image window return status
------------------------------------*/
#define E_EOIW   1

/* no data conversion required
------------------------------*/
#define SAME  0

/* 'next flag' options
----------------------*/
#define CREAT  -1
#define INIT    0
#define CURNT   1
#define NEXT    2

/* Compression type
-------------------*/
#define COMPRES 5
#define NO_COMPRES 6

/* Select the version of fseek and lseek to use for the current operating
   system.  This may need to be adjusted for different versions of 
   operating systems. */
#undef LARGE_FILES
#ifdef sgi
/* SGI supports the unix 98 standard */
#define FSEEK fseeko
typedef off_t FSEEK_OFF_T;
#define LSEEK lseek
typedef off_t LSEEK_OFF_T;
#define FTELL ftello
#define LARGE_FILES 1
#elif defined sun
/* the current EDC suns only support small files, it is suspected that
   the newer sun platforms support the unix 98 standard, but none are
   available to test on currently */
#define FSEEK fseek
typedef long FSEEK_OFF_T;
#define LSEEK lseek
typedef off_t LSEEK_OFF_T;
#define FTELL ftell
#elif defined vax
/* vax supports only small files */
#define FSEEK fseek
typedef long FSEEK_OFF_T;
#define LSEEK lseek
typedef long LSEEK_OFF_T;
#define FTELL ftell
#else
/* for all other platforms, assume no large file support */
#define FSEEK fseek
typedef long FSEEK_OFF_T;
#define LSEEK lseek
typedef off_t LSEEK_OFF_T;
#define FTELL ftell
#endif

/* forward reference to the FDESC and GDESC structures to allow
   the compile to reference them without knowing the real definition */
struct FDESC;
struct GDESC;

char *c_ealloc
(
    struct FDESC *fdarray[],    /* I: file descriptor array        */
    long         *buffersz      /* O: size of user buffer pool */
);

int c_eopenr
(
    struct FDESC **fdesc,   /* O: file descrptor array */
    char *fname,            /* I: host file name */
    long *nbands,           /* I/O: number of bands */
    long *sl,               /* I: start line of window specification */
    long *ss,               /* I: start sample of window specification */
    long *nl,               /* I/O: number of lines in window spec */
    long *ns,               /* I/O: number of samples in window spec */
    long *acc,              /* I: type of access */
    long *dtype,            /* I/O: data type */
    long opt[],             /* I: array of options */
    long *maxline           /* I: maximum number of lines to be read */
);

int c_eopens
(
    struct FDESC **fdesc,  /* I/O: file descriptor array */
    char *fname,           /* I: host file name */
    long bands[],          /* I/O: band array */
    long *sl,              /* I: start line of window specification */
    long *ss,              /* I: start sample of window specification */
    long *nl,              /* I/O: number of lines in window spec */
    long *ns,              /* I/O: number of samples in window spec */
    long offset[],         /* I: band offset array */
    long *acc,             /* I: access mode */
    long *dtype,           /* I/O: data type */
    long *nbuf,            /* I: number of previously read lines */
    long opt[]             /* I: option array */
);

long c_eread
(
    struct FDESC **fdesc,    /* I: file descriptor array */
    long *eband,             /* I: band number to be read */
    long *eline,             /* I: line number to be read */
    unsigned char *euserbuf, /* I/O: buffer */
    long *nlines             /* I: number of lines to read */
);

long c_ewrite
(
    struct FDESC **fdesc,    /* I: file descriptor */
    long *eband,             /* I: band number to write */
    long *eline,             /* I: line number to write */
    unsigned char *euserbuf, /* I: buffer */
    long *nlines             /* I: number of lines to write */
);

int c_eclose
(
    struct FDESC **fdesc  /* I: file descriptor array       */
);

int c_egroup
(
    struct FDESC  *fdarray[],  /* I: array of file descriptors */
    struct GDESC  **gdesc,     /* O: group descriptor */
    unsigned char *buffer,     /* I: address of user buffer pool */
    long          *buffersz    /* I: size of user buffer pool */
);

int c_estep
(
    struct GDESC **gdesc  /*  group descriptor */
);

int c_egclse
(
    struct GDESC **gdesc /* I: group desc from egroup */
);

int get_buf
(
    unsigned char   *buffer,
    long    *nlbuf,
    long    *nsbuf,
    long    *nlav,
    long    *nsav,
    long    *sli,
    long    *ssi,
    long    *nli,
    long    *nsi,
    long    *nlb,
    long    *nsb,
    long    *slcb,
    long    *sscb,
    long    *nextflg,
    long    *band,
    long    *fill,
    struct FDESC **fdesc,
    char    *host,
    long    *dtype
);

int put_buf
(
    unsigned char *buffer,
    long    *nlbuf,
    long    *nsbuf,
    long    *nlav,
    long    *nsav,
    long    *slbuf,
    long    *ssbuf,
    long    *nlti,
    long    *nsti,
    long    *nextflg,
    long    *band,
    long    *bands,
    struct FDESC **fdesc,
    char    *host,
    long    *dtype,
    long    *wtype
);

#endif
