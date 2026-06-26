#ifndef LASTM_H
#define LASTM_H

#define MAX_QUADS 4	/* maximum number of quads in a TM-P scene	*/
#define MAX_BANDS 7     /* maximum number of bands in a TM-P scene	*/

#define TMLINE 5965 	/* Number of lines in a TM scene		*/
#define TMSAMP 6967	/* Number of sampels in a TM scene		*/

#define TIPS_MOUNT  -1	/* Mount tape error				*/
#define TIPS_READ   -2	/* Tape read error				*/
#define TIPS_DENS   -3	/* Invalid tape density				*/
#define TIPS_TYPE   -4	/* Invalid type					*/
#define TIPS_FORMAT -5	/* Must be BSQ or BIL tape			*/
#define TIPS_WRONG  -6	/* Wrong tape error				*/
#define TIPS_OPEN   -7	/* Tape open error			 	*/
#define TIPS_WRITE  -8	/* Tape write error				*/
#define TIPS_ALLOC  -9	/* Could not allocate buffer space		*/
#define TIPS_GROUP  -10	/* Could not group the images			*/
#define TIPS_PRT    -11	/* Could not open print file			*/

/* organizations
----------------*/
#define BIL 1
#define BSQ 2

/*
** FUNCTION PROTOTYPES
*/

int ibmvx4
(
    int *inval,
    double *outval,
    long *swap
);

int ibmvx8
(
    int *inval1,
    int *inval2,
    double *outval,
    long *swap
);

void c_pronum
(
    char *string,
    long *proj_code
);

int get_annot
(
    long *unit,         /* tape unit number                           */
    long *flag,         /* flag indicating whether to save annotation */
    long bands[],       /* band specification                         */
    char *date,         /* acquisition date                           */
    char *lat_long,     /* lat long at scene center                   */
    char *path_row,     /* WRS path row                               */
    char *scene_id,     /* scene id                                   */
    char *source,       /* data source                                */
    char *resamp,       /* resampling technique                       */
    char *quality       /* overall geometric quality code             */
);

int get_calib
(
    long *unit,         /* tape unit number                             */
    long bands[],       /* Band specification                           */
    long *format        /* image file format (BSQ or BIL)               */
);

int get_head
(
    long *unit,         /* tape unit number                             */
    long *flag,         /* flag indicating whether to save header info  */
    long *format,       /* image file format (BSQ or BIL)               */
    long bands[],       /* Band specification                           */
    double fcenll[],    /* Full scene center latitude/longitude degrees */
    double fcenls[],    /* Full scene center line/sample                */
    double fcproj[],    /* Full scene center in projection coordinates  */
    double *angle,      /* Angle of rotation converted into radians     */
    long *path,         /* WRS path number                              */
    long *utmzone,      /* UTM zone number                              */
    char *direction,    /* Orbital direction                            */
    char *maprj,        /* Map projection identifier                    */
    long *line_center,  /* Line number at scene center                  */
    long *samp_center   /* sample number at scene center                */
);

int get_imgfd
(
    long *unit,         /* pointer to tape unit number                       */
    long *nl,           /* pointer to number of output lines                 */
    long *ns            /* pointer to number of output samples per line      */
);

int get_trail
(
    long *unit,         /* tape unit number                           */
    long *flag          /* flag indicating whether to save annotation */
);

int get_voldir
(
    long *unit,         /* tape unit number                             */
    long *cur_quad,     /* pointer to current quad number               */
    char *sceneid,      /* scene id                                     */
    long *format        /* image file format (BSQ or BIL)               */
);

int homprm
(
    long *path,         /* WRS path designator                          */
    long *zone,         /* HOM zone code                                */
    double parms[]      /* HOM projection parameters                    */
);

void homzon
(
    long *row,          /* WRS row                      */
    long *zone          /* HOM zone                     */
);

void ibmvx_4
(
    int *inval,
    float *outval,
    long *swap
);

void ibmvx_8
(
    int *inval1,
    int *inval2,
    double *outval,
    long *swap
);

void ls2proj
(
    double *angle,      /* Display rotation angle                             */
    long *proj_code,    /* Map projection type (SOM or UTM)                   */
    double fcproj[],    /* Full scene center projection coords (y,x)          */
    double fcenls[],    /* Full scene center line/sample                      */
    double proj[],      /* Projection coordinates calculated (y,x)            */
    long linsamp[]      /* Line/sample values calculated                      */
);

int imgproj
(
    double pwin[],      /* Input image corner projection coordinates          */
    long iwin[],        /* Input image corner line/sample coordinates         */
    long linsamp[],     /* Line/sample values calculated                      */
    double proj[]       /* Projection coordinates calculated (y,x)            */
);

int somprm
(
    long *path,         /* WRS path designator                          */
    double parms[]      /* SOM projection parameters                    */
);

int utmprm
(
    long *zone,         /* UTM zone number                              */
    double *lat,        /* Scene center latitude                        */
    double parms[]      /* UTM parameters                               */
);

void utmzon
(
    double *wrslon,     /* WRS center longitude                 */
    long *zone          /* UTM zone                             */
);

#endif
