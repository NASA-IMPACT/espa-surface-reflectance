#ifndef DCW_H
#define DCW_H

#include <stdio.h>

#define DCW_NUMCD	4       /* number of DCW CD's			      */
#define DCW_NUMCOVER    6       /* number of DCW coverages supported          */
#define DCW_VDTFILE     "int.vdt"
                                /* value definition filename                  */
#define DCW_VDTRECSIZE  80      /* number of CD containing DCW data           */
#define LINEBLOCK    1000       /* blocks size for line allocation            */

#define PO		0
#define DN		1
#define RD		2
#define RR		3
#define HY		4
#define HS		5

/*---------------------- DCW Structures ------------------------------------*/
struct DCW_CD
   {
   char         *label;         /* label of CD                               */
   char         *lib;           /* corresponding name                        */
   float        minlon,
                minlat,
                maxlon,
                maxlat;         /* bounding rectangle                       */
   };
 
/*
 *  *** Vector Map structure.
 */
struct DCW_LINE
 {
 long           attrib;         /* feature code or identification           */
 long           npts;           /* number of points in the line             */
 double         *pts;           /* list of node locations                   */
 };

struct DCW_VECTOR
 {
 struct
 DCW_LINE       *line;          /* list of lines                            */
 int            maxline;        /* number of lines currently allocated      */
 int            numlines;       /* number of lines defined                  */
 };

struct DCW_LFT
   {
   int          id;             /* line id                                   */
   int          linetype;       /* feature type of line                      */
   int          linestatus;     /* statue of line                            */
   short        tile;           /* ID tile line resides in                   */
   int          edge;           /* line ID number                            */
   double       mbr[4];         /* minimum bounding rectangle of line        */
   };

struct DCW_VDT
   {
   int          id;             /* row identifier                          */
   int          attribute;      /* attribute value                         */
   int          attindex;       /* attribute index			   */
   char         desc[50];       /* attribute value description             */
   unsigned
   char         selected;       /* flag denoting this line type was selected */
   unsigned
   char         loaded;         /* flag denoting this line type was selected */
   };
 

struct DCW_COVERAGE
   {
   char         *desc;          /* description of feature                   */
   char         *name;          /* name/directory                           */
   unsigned
   int          numvdt;         /* number of line types in list             */
   int          numstatus;      /* number of line status features           */
   int          level;          /* topology level                           */
   int          getstatus;      /* flag for line status definition to be used */
   struct
   DCW_VDT      *vdt;           /* list of line type for this feature       */
   struct
   DCW_VDT      *status;        /* list of line status features             */
   };



int dcw_checkstatus
(
    struct DCW_COVERAGE *coverage, /* DCW coverage information     */
    int status                     /* DCW line status              */
);

int dcw_checktype
(
    struct DCW_COVERAGE *coverage, /* DCW coverage information     */
    int type                       /* DCW line type                */
);

int dcw_definecover
(
    char *cdpath,                  /* location of CD drive         */
    char *library,                 /* DCW library (ex. noamer)     */
    struct DCW_COVERAGE *coverage  /* DCW coverage                 */
);

int dcw_freecover
(
    struct DCW_COVERAGE *dcw_coverage  /* DCW coverage structure to be freed */
);

int dcw_freelines
(
    struct DCW_VECTOR *dcw_vector  /* DCW coverage structure to be freed */
);

int dcw_getcd
(
    char *cdpath,               /* location of cd drive         */
    int *curlib,                /* current library              */
    struct DCW_CD *dcw_cd       /* DCW CD information structure */
);

int dcw_getlat
(
    char *tile                  /* pathname of tile of interest         */
);

int dcw_getlon
(
    char *tile                  /* pathname of tile of interest         */
);

char dcw_getrowcol
(
    double coord,               /* lat/lon coordinate of interest       */
    int minflag                 /* flag denoting min value requested    */
);

int dcw_gettileid
(
    char *cdpath,               /* location of cd drive                 */
    char *lib,                  /* DCW library (ex. eurnasia)           */
    char *tilename              /* DCW library (ex. eurnasia)           */
);

int dcw_getoffset
(
    char *cdpath,               /* location of cd drive                 */
    char *lib,                  /* current DCW library (i.e. CD         */
    char *cover,                /* coverage of interest                */
    int tileid,                 /* id of tile of interest               */
    long *recnum,               /* byte recnum file (output)            */
    long *numelem               /* number of elements in lft (output)   */
);

char dcw_lattile
(
    double lat,                 /* latgitude to find DCW directory for  */
    long minflag                /* flag denoting minimum boundary       */
);

char dcw_lontile
(
    double lon,                 /* longitude to find DCW directory for  */
    long minflag                /* denoting minimum boundary            */
);

int dcw_linequery
(
    char *cdpath,               /* location of CD drive                  */
    char *lib,                  /* current DCW library (i.e. CD          */
    struct DCW_COVERAGE *coverage, /* coverage of interest               */
    int **line,                 /* list of request line ID's             */
    int **type,                 /* list of request line types            */
    int tileid,                 /* only for debuggin                     */
    int startrec,               /* first record to read                  */
    int numrec                  /* number of records to read             */
);

int dcw_query
(
    char *cdpath,               /* location of CD drive                 */
    char *library,              /* DCW library name                     */
    struct DCW_COVERAGE *dcw_coverage,/* defines DCW categories information  */
    struct DCW_VECTOR *dcw_vector, /* resulting vector information           */
    double lllat,               /* bounding rectangle of area of interest */
    double lllon,
    double urlat,
    double urlon
);

int dcw_readlines
(
    char *cdpath,               /* location of CD drive                       */
    char *lib,                  /* DCW library (ex. noamer)                   */
    struct DCW_COVERAGE *coverage,/* DCW coverage                             */
    struct DCW_LINE **dcw_line, /* line data structure                        */
    int *dcw_numlines,          /* number of lines loaded into memory         */
    int *dcw_maxline,           /* amount of space allocated for lines        */
    char *tile,                 /* pathname of tile                           */
    int *line,                  /* list of line ID fitting query              */
    int *type,                  /* corresponding type of lines                */
    int numline,                /* number of lines in query                   */
    double lllat,               /* lower left latitude                        */
    double lllon,               /* lower left longitude                       */
    double urlat,               /* upper right latitude                       */
    double urlon,               /* upper right longitude                      */
    int clip                    /* flag denoting clipping to be performed     */
);

int dcw_read
(
    unsigned char *ptr,         /* pointer to block of data	*/
    int		size,               /* size of item in block	*/
    int		num,                /* number of items in block	*/
    FILE*	fd                  /* file descriptor of CD	*/
);

void dcw_writelt
(
    char *hdrfile,              /* name of header labeled table file          */
    struct DCW_LINE *dcw_line,  /* current set of lines                       */
    int numlines,               /* number of lines                            */
    struct DCW_COVERAGE *dcw_coverage, /* DCW coverages                       */
    double lllat,               /* lower left latitude                        */
    double lllon,               /* lower left longitude                       */
    double urlat,               /* upper right latitude                       */
    double urlon                /* upper right longitude                      */
);

#endif
