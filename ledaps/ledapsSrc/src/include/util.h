#ifndef UTIL_H
#define UTIL_H

/*===================================================
        This header file contains prototypes to both lasutil and
lasutil2 directories.
===================================================*/

#include <stdio.h>
#include "imageio.h"   /* MAXBND definition */
#include "cm.h"        /* PRS_LST definition */


/* c_4digitdt.c 
===================================================*/
long c_4digitdt
(
   long date       /* I: Date in format YYJJJ to be converted          */
);


/* c_4digityr.c
===================================================*/
long c_4digityr
(
   long year       /* I: year value to be converted                    */
);


/* c_abort.c
===================================================*/
void c_abort
(
   struct PARBLK *block,  /* I:  TAE parameter block                   */
   char message[],        /* I:  error message                         */
   char key[],            /* I:  error message key                     */
   char prog_name[],      /* I:  program name                          */
   char files[][CMLEN],   /* I:  files to be deleted                   */
   long *num_files        /* I:  number of files to be deleted         */
);


/* c_asc2eb.c
===================================================*/
void c_asc2eb
(
   char *buf,      /* I/O: string to converte, ASCII to EBCDIC         */
   long *size      /* I:   number of bytes in string                   */
);


/* c_bitops.c
===================================================*/
int c_bitand
(
   int *in1,        /* I:   The first integer value.                   */
   int *in2         /* I:   The second integer value.                  */
);

int c_bitor
(
   int *in1,        /* I:   The first integer value.                   */
   int *in2         /* I:   The second integer value.                  */
);

int c_bitxor
(
   int *in1,        /* I:   The first integer value.                   */
   int *in2         /* I:   The second integer value.                  */
);


/* c_bndchk.c
===================================================*/
int c_bndchk
(
   char *host,
   long *inbnds,
   long *vlbnds,
   long *nbands
);


/* c_clean.c
===================================================*/
void c_clean
(
   struct PARBLK *vblock, /* I:  TAE parameter block                   */
   char *pname,           /* I:  program name                          */
   char hnames[][CMLEN],  /* I:  output file names                     */
   long *count            /* I:  number of output files                */
);


/* c_complt.c
===================================================*/
void c_complt
(
   char  *pname,          /* I:  program name                          */
   char  tnames[][CMLEN], /* I:  output file names                     */
   long  *count           /* I:  number of output files                */
);


/* c_cpasso.c
===================================================*/
void c_cpasso
(
   struct PARBLK *vblock, /* I:  pointer to TAE parameter block        */
   char *hosout,          /* I:  output image name                     */
   char (*hosin)[CMLEN],  /* I:  input image names                     */
   long *numin            /* I:  number of images specified in "hosin" */
);


/* c_creat_name.c
===================================================*/
void c_creat_name
(
   char *temp_file,
   char *ext 
);


/* c_delfil.c
===================================================*/
void c_delfil
(
   char *host        /* I: host file name to be deleted                */
);


/* c_display.c
===================================================*/
void c_inittm
(
   char  *subject,    /* I: subject of message (e.g. process name)     */
   char  *verb,       /* I: verb of message (e.g. 'Processing')        */
   char  *entity,     /* I: entity of display message (e.g. 'line')    */
   char  *cat,        /* I: category of message (e.g. 'band')          */
   char  *text,       /* I: appended text string of processing msg.    */
   long  *max_entity, /* I: maximum number of 'entity' processed       */
   long  *max_cat,    /* I: maximum number of bands processed          */
   long  *seconds,    /* I: number of seconds between messages         */
   long  *flag        /* I: specifies format of display message.
                            = ALLCAT -- processing all bands at once
                            = ONECAT -- processing one band at a time
                            = GEN1   -- general display message one
                            = GEN2   -- general display message two    */
);

void c_startm(void);

void c_stoptm(void);

/* c_eb2asc.c
===================================================*/
void c_eb2asc
(
   char *buf,         /* I/O: string to converte, EBCDIC to ASCII      */
   long *size         /* I:   number of bytes in string                */
);


/* c_getcnt.c
===================================================*/
int c_getcnt
(
   struct PARBLK *block, /* I:  TAE parameter block                    */
   char *param_buf,      /* I:  Input image parameter name             */
   long *count           /* O:  Element array number of param.         */
);


/* c_getgbl.c
===================================================*/
void c_getgbl
(
   struct PARBLK *vblock,
   char histry[],
   char minmax[],
   long *msgtim
);


/* c_getin.c
===================================================*/
int c_getin
(
   struct PARBLK *block,    /* I:  TAE parameter block                 */
   char *param_buf,         /* I:  Input image parameter name          */
   long *element,           /* I:  element array number of param.      */
   char tnames[][CMLEN],    /* O:  TAE names of images specified       */
   char hnames[][CMLEN],    /* O:  host names of images specified      */
   long *images,            /* O:  Number of images                    */
   long (*window)[4],       /* O:  Image windows (sl ss nl ns)         */
   long (*bands)[MAXBND+1], /* O:  Image bands                         */
   long *nbands,            /* O:  Number of bands per image           */
   long *totbnd,            /* O:  total number of bands specified     */
   long *imglimit,          /* I:  Maximum images caller can handle    */
   long *same_window,       /* O:  1 if all windows the same           */
   long (*qflag)[2]         /* O:  1 user defined window & bands else: */
                                /* 0 if default window and bands       */
);

int c_getin2
(
   char inparm[],           /* I:  Input image string with             */
   char tnames[][CMLEN],    /* O:  TAE names of images specified       */
   char hnames[][CMLEN],    /* O:  host names of images specified      */
   long *images,            /* O:  Number of images                    */
   long (*window)[4],       /* O:  Image windows (sl ss nl ns)         */
   long (*bands)[MAXBND+1], /* O:  Image bands                         */
   long *nbands,            /* O:  Number of bands per image           */
   long *totbnd,            /* O:  total number of bands specified     */
   long *imglimit,          /* I:  Maximum images caller can handle    */
   long *same_window,       /* O:  1 if all windows the same           */
   long (*qflag)[2]         /* O:  1 user defined window & bands else: */
                                /* 0 if default window and bands       */
);


/* c_getinimg.c
===================================================*/
int c_getinimg
(
   struct PARBLK *block,    /* I:  TAE parameter block                 */
   char *param_buf,         /* I:  Input image parameter name          */
   long *element,           /* I:  element array number of param.      */
   char tnames[][CMLEN],    /* O:  TAE names of images specified       */
   char hnames[][CMLEN],    /* O:  host names of images specified      */
   long *images,            /* O:  Number of images                    */
   long (*window)[4],       /* O:  Image windows (sl ss nl ns)         */
   long (*bands)[MAXBND+1], /* O:  Image bands                         */
   long *nbands,            /* O:  Number of bands per image           */
   long *totbnd,            /* O:  total number of bands specified     */
   long *imglimit,          /* I:  Maximum images caller can handle    */
   long *same_window,       /* O:  1 if all windows the same           */
   long (*qflag)[2]         /* O:  1 user defined window & bands else: */
                                /* 0 if default window and bands       */
);


/* c_getnam.c
===================================================*/
int c_getnam
(
   struct PARBLK *vblock, /* I:  TAE parameter block                   */
   char *param,           /* I:  TAE parameter name                    */
   long *element,         /* I:  element of parameter                  */
   char *fname,           /* O:  Host file name to be returned         */
   char *ext,             /* I:  default extension                     */
   long *access           /* I:  Type of I/O access                    */
);


/* c_getnum.c
===================================================*/
int c_getnum
(
   struct PARBLK *vblock,  /* TAE parameter block                      */
   char parname[],         /* The parameter name                       */
   long *count,            /* Number of entries to read in             */
   long *dtype,            /* NEWLAS data type                         */
   unsigned char buf[]     /* Buffer to read into                      */
);


/* c_getout.c
===================================================*/
int c_getout
(
   struct PARBLK *vblock, /* I:  TAE parameter block                   */
   char *param,           /* I:  TAE parameter name                    */
   char taeout[][CMLEN],  /* O:  TAE image name(s)                     */
   char hosout[][CMLEN]   /* O:  Host image file name(s)               */
);


/* c_getparm.c
===================================================*/
int c_getparm
(
   char *filename,             /* Name of the project parameter file       */
   char *parm_name,            /* Name of the variable to be retrieved     */
   char string_value[][CMLEN], /* String values of variable being returned */
   char *type,                 /* Type of the parameter being returned     */
   double double_value[],      /* Double values of the variable returned   */
   long long_value[],          /* variable being returned                  */
   long position               /* The nth occurrence to retrieve           */
);


/* c_getstr.c
===================================================*/
int c_getstr
(
   struct PARBLK *vblock, /* TAE parameter block                           */
   char parname[],        /* The parameter name                            */
   long *count,           /* Number of entries to read in                  */
   long *slen,            /* Dimensioned length of string(s)               */
   char buf[]             /* Buffer to read into                           */
);


/* c_getsys.c
===================================================*/
int c_getsys
(
   char *system         /*  system of CPU                               */
);


/* c_gettyp.c
===================================================*/
int c_gettyp
(
   char hosin[][CMLEN], /* I:  input host file names                      */
   long *nmrimg,        /* I:  number of input images                     */
   long *odtype         /* O:  output data type -- data type of input img */
);


/* c_hsclos.c
===================================================*/
int c_hsclos
(
   long *fd,            /* I:  file descriptor of input images's history   */
   char hosin[]         /* I:  Input host name                             */
);


/* c_hsetpr.c
===================================================*/
int c_hsetpr
(
   struct PARBLK *vblock,       /* TAE parameter block                     */
   long *term,                  /* Boolean for terminal output             */
   long *file,                  /* Boolean for file output                 */
   long  *lp,                   /* Boolean for line printer output         */
   char  *prtfile,              /* file name of file or lp working file    */
   char  *filename              /* file name of file or lp output file     */
);


/* c_hsopen.c
===================================================*/
long c_hsopen
(
   long *fd,         /* I:  file descriptor of input images's history      */
   char hosin[],     /* I:  input image name                               */
   long *access      /* I:  type of access to input image history file     */
);


/* c_hsread.c
===================================================*/
int c_hsread
(
   long *fd,         /* I:  file descriptor of input images's history      */
   char key[],       /* I/O:  key buffer                                   */
   char **buf,       /* O:  history buffer                                 */
   long *len         /* I:  char length of label service record            */
);


/* c_hstout.c
===================================================*/
int c_hstout
(
   struct PARBLK *vblock, /* I:  pointer to TAE parameter block            */
   char *hosout,          /* I:  image whose history is to be updated      */
   char (*hosin)[CMLEN],  /* I:  image names whose history is to be copied */
   long *numin,           /* I:  number of images specified in "hosin"     */
   char *histext          /* I:  text to be included in history record     */
);

int cphist
(
   long fd,               /* I:  file descriptor to opened history file    */
   char (*hosin)[CMLEN],  /* I:  input image names                         */
   long *numin            /* I:  number of input images                    */
);

int formhist
(
   long fd,               /* I:  file descriptor to opened history file    */
   struct PARBLK *vblock, /* I:  pointer to TAE parameter block            */
   char *text             /* I:  history text supplied by calling program  */
);


/* c_imgprs.c
===================================================*/
int c_imgprs
(
   char *i_specs,
   struct PRS_LST *list,
   long *images
);

void compress
(
   char *filled,
   char *small
);

int image_info
(
   char *spec,
   struct PRS_LST *i_record
);

int make_window
(
   char *numerals,
   double window[]
);

int make_bands
(
   char *numerals,
   long bands[]
);

int make_proj
(
   char *numerals,
   long proj[]
);


/* c_index.c
===================================================*/
void c_index
(
   char string1[], /* I: String for comparison.                            */
   char string2[], /* I: Possible substring.                               */
   int  *pos       /* O: Pos. of string2 in string1 (if not a substring)   */
);


/* c_inlas.c
===================================================*/
void c_inlas
(
   struct PARBLK *vblock,       /* O:  TAE parameter block                 */
   long *sout                   /* O:  unit number for opening stdout      */
);


/* c_inmnmx.c
===================================================*/
int c_inmnmx
(
   long *dtype,       /* I:  data type                                 */
   long *totbnd,      /* I:  total number of bands                     */
   double *minval,    /* O:  minimum intensity values array            */
   double *maxval     /* O:  maximum intensity values array            */
);


/* c_isleap.c
===================================================*/
int c_isleap
(
   long *year             /* I:  Year to test                           */
);


/* c_julian.c
===================================================*/
void c_julgre
(
   long *year,          /* I:  Year                                     */
   long *jdate,         /* I:  Julian date                              */
   long *month,         /* O:  Gregorian month                          */
   long *day            /* O:  Gregorian day                            */
);

void c_grejul
(
   long *year,          /* I:  Year                                     */
   long *month,         /* I:  Gregorian month                          */
   long *day,           /* I:  Gregorian day                            */
   long *jdate          /* O:  Julian date                              */
);


/* c_low2up.c
===================================================*/
void c_low2up
(
   char *buf,  /* I/O: buffer to convert                       */
   long *size  /* I:  size of string in buffer                 */
);


/* c_mkimg.c
===================================================*/
int c_mkimg
(
   char *in_string,
   char *ext,
   char *tname,
   char *hname,
   long *flag
);


/* c_mknam.c
===================================================*/
int c_mknam
(
   char *filename,      /* I:  Input file name                         */
   char *outfile,       /* O:  Host file name to be returned           */
   char *ext,           /* I:  default extension                       */
   long *access         /* I:  Type of I/O access                      */
);


/* c_pterm.c
===================================================*/
void c_pterm
(
   char *buf,    /* O: character buffer to output                       */
   long *flag    /* I: flag to query user for command                   */
);


/* c_putparm.c
===================================================*/
int c_putparm
(
   char *filename,         /* I: Hostname of the project parameter file */
   char *parmname,         /* I: Parameter name                         */
   char string_value[][CMLEN], /* I: String values of the parameter     */
   char *type,             /* I: Type of parameter to write (I, R, S)   */
   double double_value[],  /* I: The double values of the parameter     */
   long dimension,         /* I: Dimension of the parameter             */
   long long_value[]       /* I: The long values of the parameter       */
);


/* c_setprt.c
===================================================*/
int c_setprt
(
   struct PARBLK *vblock, /* I:  TAE parameter block                     */
   long *term,            /* O:  Boolean for terminal output             */
   long *file,            /* O:  Boolean for file output                 */
   long *lp,              /* O:  Boolean for line printer output         */
   char *prtfile,         /* I/O:  file name of file or lp working file  */
   char *filename         /* O:  file name of file or lp output file     */
);


/* c_sigblk.c
===================================================*/
void c_sigblk
(
   int *mask
);

void c_sigunblk
(
   int *mask
);


/* c_splitname.c
===================================================*/
void c_splitname
(
   char filespec[CMLEN],  /* I:  File specification                    */
   char directory[CMLEN], /* O:  Directory portion of filespec         */
   char root[CMLEN],      /* O:  Root portion of the file name         */
   char extension[CMLEN]  /* O:  Extension portion of the file name    */
);


/* c_sysset.c
===================================================*/
int c_sysset
(
   char *csys,       /* I/O: character coded system (e.g. "gould-utx")  */
   long *nsys,       /* I/O: system numerical code (e.g. SYS_GOULD_UTX) */
   long *sysflg      /* I:  System conversion flag:                     */
                          /*  SETN = 0 - convert system string to value */
                          /*  SETS = 1 - convert system value to string */
);


/* c_txtout.c
===================================================*/
void c_txtout
(
   char  buf[],   /* I:  Text to be written out                     */
   long  print[], /* I:  Print options (FILE, LP, TERM)             */
   FILE  *fp      /* I:  File descriptor for print FILE or LP       */
                  /*     if fp == 0 then no file open               */
);


/* c_up2low.c
===================================================*/
void c_up2low
(
   char *buf,  /* I/O: Buffer converted from upper - lower */
   long *size  /* I:  munber of bytes to convert in buffer */
);


/* c_wndchk.c
===================================================*/
int c_wndchk
(
   char *host,      /* I:  Host name                                 */
   double inwind[], /* I:  Input window to be checked                */
   long vlwind[]    /* O:  Valid window                              */
);


/* c_wndconv.c
===================================================*/
int c_wndconv
(
   char   hname[],   /* I:  host names of images specified              */
   double window[],  /* I:  Image windows (sl ss nl ns)                 */
   long   units[]    /* I:  winow units                                 */
);


/* getusern.c
===================================================*/
int getusern
(
   char *username
);


/* logit.c
===================================================*/
int logit
(
   char *progname
);


/* terminal.c
===================================================*/
void init_term
(
   long *term_type, /* O:  terminal type                               */
   long *totlin,    /* O:  total number of lines                       */
   long *totcol     /* O:  total number of columns                     */
);

void write_term
(
   char *buf,
   long cc
);

void read_term
(
   char *buf,
   long buf_len,
   long *term
);

void position_term
(
   long line,
   long column
);

void clear_term(void);

void bell_term(void);

void out_term
(
   long line,
   long col,
   char *buf
);

int in_term
(
   long line,
   long col,
   char *buf,
   long *len,
   long *term
);


#endif
