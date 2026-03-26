#ifndef DDR_H
#define DDR_H

#include "tae.h"
#include "imageio.h"

/*This is the 'C' include file for the DDR routines			     */


struct DDR                /* ddr for integer data                    */
    {
    long nl;              /* number of lines in image                */
    long ns;              /* number of samples in image              */
    long nbands;          /* number of bands in image                */
    long dtype;           /* data type of pixels:                    */
                          /*   =1:  unsigned char                    */
			  /*   =2:  short                            */
                          /*   =3:  long                             */
                          /*   =4:  float                            */
    long master_line;     /* line relative to master image	     */
    long master_sample;   /* sample relative to master image         */
    long valid[8];	  /* valid flags:                            */
                          /*   =0:  not valid                        */
                          /*   =1:  valid                            */
 		      	  /*   =2:  unknown		             */
    long proj_code;       /* GCTP projection code                    */
                          /*  refer to GCTP document for valid codes */
    long zone_code;       /* UTM or State Plane zone                 */
                          /*  refer to GCTP document for valid codes */
    long datum_code;      /* GCTP datum code                         */
                          /*  refer to GCTP document for valid codes */
    long spare;		  /* spare integer value; for future use or
			     expansion				     */

    char system[12];       /* computer system data is on (with NULL)  */
    char proj_units[12];   /* Projection units (GCTP units+other) 
			     (with NULL)    */
    char last_used_date[12]; /* last access date	(with NULL)	     */
			    /* NOTE: All ddr dates are stored as     */
			    /* "dd-mmm-yy".  For example, a date     */
			    /* of December 31, 1986 is stored as     */
			    /* "31-dec-86" with the month in lower   */
			    /* case.				     */
    char last_used_time[12]; /* last access time	(with NULL)          */
			    /* NOTE: All ddr times are stored        */
			    /* using a twenty-four hour clock.	     */
			    /* Seconds are speperated by a colon.    */
			    /* Example: 1:05:55 pm is stored as      */
			    /* 1305:55				     */

    double proj_coef[15]; /* GCTP projection parameters              */
                          /*  refer to GCTP document for field definitions*/
    double upleft[2];     /* Corner coordinates entered              */
    double loleft[2];     /*   in latitude/longitude or              */
    double upright[2];    /*   northing/easting order or             */
    double loright[2];    /*   projection coordinates (y/x)          */
    double pdist_y;       /* projection distance/pixel (y-direction) */
    double pdist_x;       /* projection distance/pixel (x-direction) */
    double line_inc;      /* line increment for sampling             */
    double sample_inc;    /* sample increment for sampling           */
                          /* NOTE:  The line/sample increments are   */
                          /*   the values applied to the original    */
                          /*   image to obtain this image.  If re-   */
                          /*   or sub-sampling was not applied, the  */
                          /*   value contained in these fields will  */
                          /*   be 1.0                                */
    };

struct BDDR
    {
    long bandno;	  /* band no. of the record 		     */
    long valid;		  /* min/max validity flag of the band	     */
                          /*   =0:  not valid                        */
                          /*   =1:  valid                            */
 		      	  /*   =2:  bounded		             */

    double minval;	  /* minimum value of the band		     */
    double maxval;	  /* maximum value of the band		     */

    char source[32];	  /* source of data (with NULL)		     */
    char instrument[32];  /* type of sensor (with NULL)		     */
    char direction[64];	  /* direction of capture process(with NULL) */
    char date[10];	  /* capture date (with NULL)		     */
			  /* NOTE: All ddr dates are stored as       */
			  /* "dd-mmm-yy".  For example, a date       */
			  /* of December 31, 1986 is stored as	     */
			  /* "31-dec-86" with the month in lower     */
			  /* case.				     */
    char time[8];	  /* capture time (with NULL)		     */
			  /* NOTE: All ddr times are stored          */
			  /* using a twenty-four hour clock.	     */
			  /* Seconds are speperated by a colon.	     */
			  /* Example: 1:05:55 pm is stored as        */
			  /* 1305:55				     */
    };

/*  Below are the constants to be used as keys within GETDDF and PUTDDF
-----------------------------------------------------------------------*/
#define DDNL 0		  /* number of lines			       */
#define DDNS 1		  /* number of samples			       */
#define DDNB 2		  /* number of bands			       */
#define DDDTYP 3	  /* Data type               	 	       */
#define DDML 4		  /* Master line    			       */
#define DDMS 5		  /* Master sample 		               */
#define DDVFLG 6	  /* Validity flag array     		       */
#define DDPCOD 7	  /* Projection code			       */
#define DDZCOD 8	  /* Zone code      			       */
#define DDDCOD 9	  /* Datum code  		     	       */
#define DDSYS 10	  /* System      		               */
#define DDUNIT 11	  /* Projection Unit		               */
#define DDLDAT 12	  /* Last use dated   		               */
#define DDLTIM 13	  /* Last used time		               */
#define DDPCF 14	  /* Projection coeficients		       */
#define DDUL 15		  /* Upper left     		               */
#define DDLL 16		  /* Lower left       		               */
#define DDUR 17		  /* Upper right       		               */
#define DDLR 18		  /* Lower right            		       */
#define DDPDY 19	  /* Projection distance-y  		       */
#define DDPDX 20	  /* Projection distance-x   		       */
#define DDLINC 21	  /* Line increment          		       */
#define DDSINC 22	  /* Sample increment        		       */

/*  Below are the constants to be used as keys within GETBDF and PUTBDF
-----------------------------------------------------------------------*/
#define DDBAND 100	  /* band number    			       */
#define DDMMV 101	  /* Min/max validity flag 		       */
#define DDMIN 102	  /* Minimum value    	      		       */
#define DDMAX 103	  /* Maximum value    	      		       */
#define DDSRC 104	  /* Source           	      		       */
#define DDINST 105	  /* Instrument      		     	       */
#define DDDIR 106	  /* Direction       		      	       */
#define DDCDAT 107	  /* Capture date    		      	       */
#define DDCTIM 108	  /* Capture time    			       */

/* Below are the constants to retrieve the correct validity flag from the
   validity flag array						         
-----------------------------------------------------------------------*/
#define DDPCV 0	  	/* projection code validity flag	       */
#define DDZCV 1	  	/* zone code validity flag		       */
#define DDDCV 2	  	/* datum code validity flag		       */
#define DDPPV 3	  	/* projection parameters validity flag         */
#define DDPUV 4	  	/* ground units validity flag		       */
#define DDPDV 5	  	/* ground distance validity flag	       */
#define DDCCV 6	  	/* corner coordinates validity flag	       */
#define DDINCV  7	/* line/sample increments validity flag        */

/*  Below of Misc. DDR constants
------------------------------------------------------------------------*/
#define E_PROT -2         /* return status of write protected           */
#define INVAL 0	 	  /* invalid validity flag status		*/
#define VALID 1		  /* valid validity flag			*/
#define UNKNOW 2	  /* unknown and equal validity flag		*/
#define BOUND 2 	  /* bounded min/max validity flag              */
#define COMB 1		  /* images are being combined into one band    */
#define NOCOMB 0 	  /* images are not being combined into one band*/

/*  Below are constants used only within the DDR support routines
------------------------------------------------------------------------*/
#define USAME 2		  /* unknown and equal validity flag	 	*/
#define UDIFF 3		  /* unknown and not equal validity flag	*/

#define DISIZE 18	  /* number of integers in record 1 of DDR      */
#define DDSIZE 27	  /* number of doubles in record 2 of DDR       */
#define DBSIZE 2	  /* number of doubles in band record(s) of DDR */

#define DDTSLN 152	  /* length of all strings in band record of DDR*/
#define DDSRLN 32	  /* length of the SOURCE character string      */
#define DDINLN 32	  /* length of the INSTRUMENT character string  */
#define DDDRLN 64	  /* length of the DIRECTION character string   */
#define DDCDLN 10	  /* length of the capture DATE character string*/
#define DDCTLN 8	  /* length of the capture TIME character string*/
#define DDBNLN 4	  /* length of the band no. character string    */
#define DDVLLN 2	  /* length of the valid flag character string  */

#define DDSYLN 12	  /* length of the SYSTEM character string      */
#define DDPULN 12	  /* length of the PROJ_UNITS character string  */
#define DDLDLN 12	  /* length of the LAST_USED_DATE char. string  */
#define DDLTLN 12	  /* length of the LAST_USED_TIME char. string  */
#define DDSTCT 4	  /* No. of string in record 1 of DDR file      */
#define DDNVAL 8	  /* No. of valid flags stored in the DDR file  */

/*  Below are constants used in the process of finding file corners
    when image corners are known
------------------------------------------------------------------------*/
#define UPLEFT_Y 0	/* upper left y/line				*/
#define UPLEFT_X 1	/* upper left x/sample				*/
#define UPRIGHT_Y 2	/* upper right y/line				*/
#define UPRIGHT_X 3	/* upper right x/sample				*/
#define LOLEFT_Y 4	/* lower left y/line				*/
#define LOLEFT_X 5	/* lower left x/sample				*/
#define LORIGHT_Y 6	/* lower right y/line				*/
#define LORIGHT_X 7	/* lower right x/sample				*/

/*
** FUNCTION PROTOTYPES
*/

/* c_calcor.c
-----------------------------------------------------------------*/
void c_calcor
(
   double *sl,      /*  I:  starting line                               */
   double *ss,      /*  I:  starting sample                             */
   double *nl,      /*  I:  number of lines                             */
   double *ns,      /*  I:  number of samples                           */
   struct DDR *ddr, /*  I:  input DDR structure                         */
   double *upleft,  /*  O:  Upper left corner                           */
   double *upright, /*  O:  Upper right corner                          */
   double *loleft,  /*  O:  Lower left corner                           */
   double *loright  /*  O:  Lower right corner                          */
);


/* c_comsys.c
-----------------------------------------------------------------*/
int c_comsys
(
   char *hostname,  /*  I:  host file name                              */
   long *flag       /*  O:  output flag                                 */
);


/* c_getbdf.c
-----------------------------------------------------------------*/
int c_getbdf
(
   struct BDDR *bddr,    /* I:  input band depenent record               */
   long *key,            /* I:  specified key of field needed            */
   unsigned char *value  /* O:  output value of specified field          */
);


/* c_getbdr.c
-----------------------------------------------------------------*/
int c_getbdr
(
   char *hname,        /* I:  input host image name                       */
   struct BDDR *bddr,  /* O:  output band depenent record                 */
   long *band          /* I:  specified band number                       */
);


/* c_getddf.c
-----------------------------------------------------------------*/
int c_getddf
(
   struct DDR *ddr,      /* I:  input band depenent record                */
   long *key,            /* I:  specified key of field needed             */
   unsigned char *value  /* O:  output value of specified field           */
);


/* c_getddr.c
-----------------------------------------------------------------*/
int c_getddr
(
   char *hname,     /* I:   Host image file name                          */
   struct DDR *ddr  /* O:  DDR structure                                  */
);


/* c_intbdr.c
-----------------------------------------------------------------*/
void c_intbdr
(
   struct BDDR *bddr  /* I/O:  bddr data                                  */
);


/* c_intddr.c
-----------------------------------------------------------------*/
void c_intddr
(
   struct DDR *ddr    /* I/O:  ddr data in record 1 and 2                 */
);


/* c_lused.c
-----------------------------------------------------------------*/
int c_lused
(
   char *hostname     /* I/O:  host name of associated DDR                */
);


/* c_minmax.c
-----------------------------------------------------------------*/
int c_minmax
(
   char    hname[][CMLEN],      /* I:  string descriptor of host name     */
   long    *nimg,               /* I:  total number of images             */
   long    window[][4],         /* I:  starting line of image             */
   long    bands[][MAXBND + 1], /* I:  array of band numbers needed       */
   long    *nbands,             /* I:  number of bands for each image     */
   double  *minval,             /* O:  minimum value array                */
   double  *maxval              /* O:  maximum value array                */
);


int ddrconv
(
   char *hostname   /* I:  Host name of the image file.                   */
);


/* c_putbdf.c
-----------------------------------------------------------------*/
int c_putbdf
(
   struct BDDR *bddr,   /* I/O:  input band depenent record               */
   long *key,           /* I:  specified key of field needed              */
   unsigned char *value /* I:  output value of specified field            */
);


/* c_putbdr.c
-----------------------------------------------------------------*/
int c_putbdr
(
   char *hname,       /* I:  input host image name                       */
   struct BDDR *bddr  /* I:  band data for the DDR                       */
);


/* c_putddf.c
-----------------------------------------------------------------*/
int c_putddf
(
   struct DDR *ddr,     /* I/O:  input band depenent record               */
   long *key,           /* I:  specified key of field needed              */
   unsigned char *value /* I:  output value of specified field            */
);


/* c_putddr.c
-----------------------------------------------------------------*/
int c_putddr
(
   char *hname,      /* I:  host image file name                         */
   struct DDR *ddr   /* I:  ddr data in record 1 and 2                   */
);


/* c_rotcor.c
-----------------------------------------------------------------*/
int c_rotcor
(
      /* The following parameters with 8 elements are
         ordered as follows:
           [0] = upper left y
           [1] = upper left x
           [2] = lower left y
           [3] = lower left x
           [4] = upper right y
           [5] = upper right x
           [6] = lower right y
           [7] = lower right x     */

   double proj[8],     /* I:  Input image corner projection coordinates  */
   long image[8],      /* I:  Input image corner line,sample coordinates */
   double window[4],   /* I:  Window applied to image                    */
   double proj_wind[8] /* O:  Output image corner projection coordinates */
);


/* c_rtcoef.c
-----------------------------------------------------------------*/
int c_rtcoef
(
      /* The following parameters with 8 elements are
         ordered as follows:
           [0] = upper left y
           [1] = upper left x
           [2] = lower left y
           [3] = lower left x
           [4] = upper right y
           [5] = upper right x
           [6] = lower right y
           [7] = lower right x     */

   double proj[8],  /* I:  Input image corner projection coordinates    */
   long image[8],   /* I:  Input image corner line,sample coordinates   */
   double coef[6],  /* O:  Transformation coefficients                  */
   long *flag       /* I:  Flag for error messages                      */
);


/* c_rtiddr.c
-----------------------------------------------------------------*/
int c_rtiddr
(
   struct DDR *inddr,  /* I/O:  Input DDR structure to be updated.      */
   long *rot           /* I:  Rotation amount in 90 degree increments:
                                        1 =  90 degrees
                                        2 = 180 degrees
                                        3 = 270 degrees
                                        4 = Horizontal flip
                                        5 = Vertical flip               */
);


/* c_svmnmx.c
-----------------------------------------------------------------*/
int c_svmnmx
(
   char *hname,     /* I:  input host image name                        */
   long *totbnd,    /* I:  number of bands to be processed              */
   double minval[], /* I:  minimum intensity values                     */
   double maxval[], /* I:  maximum intensity values                     */
   char *m_flag     /* I:  min/max validity flag ("YES"/"NO"/"BOUND")   */
                       /*  from $MINMAX global                          */
);


/* c_upddr.c
-----------------------------------------------------------------*/
int c_upddr
(
   struct PARBLK *parblk,   /* I:  TAE parameter block                  */
   char host_in[][CMLEN],   /* I:  array of input host names            */
   char *host_out,          /* I:  Output host name                     */
   long *nmrimg,            /* I:  Total number of input images         */
   long window[][4],        /* I:  Window specifications for each image */
   long bands[][MAXBND +1], /* I:  band numbers for each input image    */
   long nbands[],           /* I:  number of bands for each input image */
   double *line_inc,        /* I:  increment in the line direction      */
   double *samp_inc,        /* I:  increment in the sample direction    */
   long *out_nbands,        /* O:  number of output bands               */
   long *cflag              /* I:  combination of bands flag            */
);


/* ckmnmx.c
-----------------------------------------------------------------*/
int ckmnmx
(
   char *hname,     /* host image name                                  */
   long *bands,     /* array containing bands that were processed       */
   long nbands,     /* no. of bands processed                           */
   double *minval,  /* pointer to array of minimum pixel values         */
   double *maxval,  /* pointer to array of maximum pixel values         */
   long *valid      /* array containing a flag stating if valid         */
);


/* clmnmx.c
-----------------------------------------------------------------*/
int clmnmx
(
   char *hostname,  /* host name of input image                         */
   long *window,    /* window specified for the input image             */
   long *cl_bands,  /* bands that the min/max values are needed         */
   long *mmbands,   /* the correct place in the min/max arrays that
                                   the min/max values should be placed  */
   double *minval,  /* minimum value array                              */
   double *maxval   /* maximum value array                              */
);


/* com_ppar.c
-----------------------------------------------------------------*/
void com_ppar
(
   struct DDR ddr1,    /*  First input DDR                              */
   struct DDR ddr2,    /*  First input DDR                              */
   long *ppar_flag,    /*  Processing projection parameter flag         */
   long default_flag,  /*  default flag -- UNKNOW or VALID              */
   long image1,        /*  number of first input image                  */
   long image2         /*  number of second input image                 */
);

void comp_prec
(
   double ddr1_value,    /*  first ddr value                             */
   double ddr2_value,    /*  second ddr values                           */
   long   *ppar_flag,    /*  Processing projection parameter flag        */
   long   prec,          /*  Precision value                             */
   long   default_flag,  /*  default flag -- UNKNOW or VALID             */
   long   image1,        /*  number of first input image                 */
   long   image2         /*  number of second input image                */
);

void comp_sm_esq
(
   struct DDR ddr1,    /*  First input ddr                               */
   struct DDR ddr2,    /*  Second input ddr                              */
   long *ppar_flag,    /*  Processing projection paramter flag           */
   long default_flag,  /*  Default flag UNKNOW or VALID                  */
   long image1,        /*  number of first input image                   */
   long image2         /*  number of second input image                  */
);

double ecc_sq
(
   double semi_major,   /*  semi major axis                              */
   double semi_minor    /*  semi minor axis                              */
);

double smajor_tbl
(
   long datum    /*  datum code                                           */
);

double sminor_tbl
(
   long datum     /*  datum code                                        */
);


/* convert_ddr.c
-----------------------------------------------------------------*/
int convert_ddr
(
   struct DDR *ddr,  /* input ddr                                      */
   char *units       /* units the ddr fields are to be converted to    */
);


/* set_capture.c
-----------------------------------------------------------------*/
void set_capture
(
   struct BDDR *in_bddr,   /*  input band depentent DDR                   */
   struct BDDR *out_bddr,  /*  output band depentent DDR                  */
   long nimg               /*  number of images                           */
);


/* set_corners.c
-----------------------------------------------------------------*/
void set_corners
(
   long corner_flag,    /*  Processing corner validity flag               */
   long proj_flag,      /*  Processing projection validity flag           */
   long pdist_flag,     /*  Processing proj. dist. validity flag          */
   long zone_flag,      /*  Processing zone validity flag                 */
   long ppar_flag,      /*  Processing proj. parameters validity flag     */
   struct DDR *in_ddr,  /*  Input DDR structure                           */
   struct DDR *out_ddr, /*  Output DDR structure                          */
   double line_inc,     /*  Line increment                                */
   double samp_inc,     /*  Sample increment                              */
   long window[]        /*  Window specifications                         */
);


/* set_datum.c
-----------------------------------------------------------------*/
void set_datum
(
   struct DDR in_ddr[],  /*  input DDRs                                   */
   struct DDR *out_ddr,  /*  Output DDR                                   */
   long nimg             /*  number of input images                       */
);


/* set_inc_flag.c
-----------------------------------------------------------------*/
void set_inc_flag
(
   struct DDR ddr[], /*  Input DDRs                                       */
   long *inc_flag,   /*  Processing increment flag                        */
   long *inc_index,  /*  Index to first DDR with a valid inc flg          */
   long nimg         /*  Total number of input DDRs                       */
);


/* set_master.c
-----------------------------------------------------------------*/
void set_master
(
   struct DDR *in_ddr,  /* input DDR structure                            */
   struct DDR *out_ddr, /* output DDR structure                           */
   long inc_flag,       /* increment flag                                 */
   long window[]        /* window specification of input image            */
);


/* set_pdist.c
-----------------------------------------------------------------*/
void set_pdist
(
   long pdist_flag,     /*  input projection distance flag                 */
   struct DDR *in_ddr,  /*  input DDR structure                            */
   struct DDR *out_ddr, /*  output DDR structure                           */
   double line_inc,     /*  input line increment value                     */
   double samp_inc      /*  input sample increment value                   */
);


/* set_ppar.c
-----------------------------------------------------------------*/
void set_ppar
(
   struct DDR in_ddr[],  /*  input DDRs                                   */
   struct DDR *out_ddr,  /*  output DDR                                   */
   long *ppar_flag,      /*  Processing projection pararmeter flag        */
   long *proj_flag,      /*  Projection flag                              */
   long nimg             /*  Total number of input DDRs                   */
);


/* set_proj.c
-----------------------------------------------------------------*/
void set_proj
(
   struct DDR in_ddr[],  /* input DDRs                                    */
   struct DDR *out_ddr,  /* output DDR                                    */
   long *proj_flag,      /* Processing projection flag                    */
   long nimg             /* number of input DDRs                          */
);


/* set_punit.c
-----------------------------------------------------------------*/
void set_punit
(
   struct DDR in_ddr[],  /*  input DDRs                                   */
   struct DDR *out_ddr,  /*  output DDR                                   */
   long nimg             /*  number of input DDRs                         */
);


/* set_zone.c
-----------------------------------------------------------------*/
void set_zone
(
   struct DDR in_ddr[],  /*  input DDRs                                   */
   struct DDR *out_ddr,  /*  output DDR                                   */
   long *zone_flag,      /*  Processing zone validity flag                */
   long nimg             /*  number of input DDRs                         */
);


/* setcornflag.c
-----------------------------------------------------------------*/
void set_corn_flag
(
   struct DDR ddr[],  /* input integer DDR fields                          */
   long *corner_flag, /* processing flag for the corner field              */
   long (*window)[4], /* input window                                      */
   long *copy_index,  /* index of the first valid corner field             */
   long pdist_flag,   /* processing flag for the proj. dist. field         */
   double tgbl,       /* tolerance global                                  */
   long pdist_index,  /* index of the first valid proj. dist. field        */
   long gbl_flag,     /* TAE global validity flag                          */
   long nimg          /* total number of input images                      */
);

void compare
(
   struct DDR *ddr1,  /*  first input ddr                                  */
   struct DDR *ddr2,  /*  second input ddr                                 */
   long window1[],    /*  first input window                               */
   long window2[],    /*  second input window                              */
   double pdist_y,    /*  projection distance in Y direction               */
   double pdist_x,    /*  projection distance in X direction               */
   double tgbl,       /*  TAE tolerance global                             */
   long *corner_flag, /*  Processing corner flag                           */
   long gbl_flag,     /*  TAE global validity flag                         */
   long image1,       /*  current image                                    */
   long image2,       /*  comparison image                                 */
   int line_inv,
   int samp_inv
);


/* setincrement.c
-----------------------------------------------------------------*/
int set_increment
(
   long inc_flag,       /*  processing increment flag                      */
   struct DDR *in_ddr,  /*  input DDR                                      */
   struct DDR *out_ddr, /*  output DDR                                     */
   double in_line_inc,  /*  line increment                                 */
   double in_samp_inc   /*  sample increment                               */
);


/* setpdisflag.c
-----------------------------------------------------------------*/
void set_pdis_flag
(
   struct DDR ddr[],   /*  input DDRs                                       */
   long *pdist_flag,   /*  Processing projection distance flag              */
   long *pdist_index,  /*  Index of first DDR with a valid pdist flag       */
   long nimg           /*  Number of input DDRs                             */
);


/* upmnmx.c
-----------------------------------------------------------------*/
int upmnmx
(
   char *hname,      /* input host image name                              */
   long *cl_bands,   /* bands of image that were processed                 */
   long *mmbands,    /* bands of min/max arrays to be used                 */
   double minval[],  /* minimum intensity values                           */
   double maxval[]   /* maximum intensity values                           */
);


#endif
