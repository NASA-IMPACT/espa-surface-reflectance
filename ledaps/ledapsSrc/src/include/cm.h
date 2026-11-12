#ifndef CM_H
#define CM_H

#include "tae.h"
#include "imageio.h" /* MAXBND definition */
#include "desc.h"


/*****************************************************************************
 *                                                                           *
 *      This is the include file for the 'cm' routines.                      *
 *                                                                           *
 *      Modified by D. Akkerman, USGS/EROS Data Center, 5/90.                *
 *       - Added the TAE_UPATH define for changes to the cmhost and          *
 *         cmtae routines.                                                   *
 *                                                                           *
 *      Modified by B. Davis, USGS/EROS Data Center, 9/90.                   *
 *       - added the cm tape structures and defines for cmstor/cmrestor.     *
 *                                                                           *
 *      Modified by D. Etrheim, USGS/EROS Data Center, 11/90.                *
 *       - added online/offline/both defines                                 *
 *****************************************************************************/

#define TAE_EXT  ';'
#define HOST_EXT '.'
#define TAE_SEP  '.'
#define HOST_SEP '_'

#define TAE_DIR_BEG '['
#define TAE_DIR_SEP '.'
#define TAE_DIR_END ']'
#define TAE_PARENT  '-'
#define TAE_UPATH   '%'
#ifdef vms
#define HOST_DIR_BEG '['
#define HOST_DIR_SEP '.'
#define HOST_DIR_END ']'
#define HOST_PARENT  '-'
#else
#define HOST_DIR_BEG '/'
#define HOST_DIR_SEP '/'
#define HOST_DIR_END '/'
#define HOST_PARENT  '.'
#endif

#define ONLINE 1
#define OFFLINE 2
#define BOTH 3

/*
   for Catalog Manager tapes
*/
struct TAPELIST
        {
        char    filestatus;             /* status of the file to stor   */
                                        /* blank - stored on tape       */
                                        /* D     - marked for Deletion  */
                                        /* H     - Hidden               */
        char    lassys[11];             /* system name from syschar     */
        char    data_type[4];           /* code for image data type     */
        char    lines[5];               /* number of image lines        */
        char    samples[5];             /* number of image samples      */
        char    bands[3];               /* number of image bands        */
        char    tape_id[6];             /* name of stor tape            */
        char    timestamp[11];          /* file create date YYDDDHHMMSS */
        char    feet[9];                /* approx feet written by file  */
        char    totft[9];               /* total feet utilized this tape*/
        char    bpi[5];                 /* tape density                 */
        char    storetime[11];          /* date and time the file       */
                                        /* described by this record was */
                                        /* stored - YYDDDHHMMSS format  */
        char    file_no[4];             /* file no. in this tapes sequence*/
        char    tapeloc;                /* physical tape location code  */
        char    end_space;              /* for null-terminated raw i/o  */
        } ;

struct TAPENAME
        {
        char    filestatus;             /* T - denotes this record is a */
                                        /*     tape name                */
        char    loc_chain[2];           /* location in tape chain       */
                                        /*   "HD" - head tape           */
                                        /*   "CO" - continuation tape   */
        char    tape_no[2];             /* number of tape in the chain  */
        char    back_chain[6];          /* tape chain back pointer:     */
                                        /* previous tape label          */
        char    blk_count[6];           /* number of blocks written     */
        char    space1[12];             /* blank filled                 */
        char    tape_id[6];             /* name of stor tape            */
        char    space2[11];             /* blank filled                 */
        char    feet[9];                /* approx feet written on tape  */
        char    totft[9];               /* total feet on this tape      */
        char    bpi[5];                 /* tape density                 */
        char    timestamp[11];          /* tape init date YYDDDHHMMSS   */
        char    file_no[4];             /* files on this tape           */
        char    tapeloc;                /* physical tape location code  */
        char    end_space;              /* for null-terminated raw i/o  */
        } ;


#define COMPUTER_ROOM '0'               /* single digit character code  */
                                        /* for physical location of tape*/
#define TAPE_LIBRARY '1'                /* single digit character code  */
                                        /* for physical location of tape*/


#define TPLSLEN 86                      /* each field of TAPELIST is one
                                          greater in length than required
                                          for the data it holds in order
                                          to allow it to contain a null-
                                          terminated string             */

#define MAXFILES 100                    /* maximum number of files to be
                                           associated with and image,
                                           in order to help limit the
                                           size of array declarations.  */

#define DUPE 2                          /* return status of a tape dir-
                                           ectory call, meaning the file
                                           already exists on tape, or in
                                           the catalog.  Sometimes this
                                           is desired, sometimes not    */

#define MAXTPWT 32768                   /* maximum number of bytes that
                                           can be streamed together and
                                           read/written to/from a single
                                           buffer by a single read/write.
                                           This value is the smallest one
                                           of the current systems,
                                           sun-bsd.  If LAS is ported to
                                           a new different machine with a
                                           smaller MAX, this value will
                                           need to be adjusted.         */

#define TPIDLEN 7                       /* to help limit the size of 2d   
                                           arrays by not using CMLEN.
                                           corresopnds to TAPELIST.tape_id,
                                           with an extra char for
                                           null-terminated strings.     
                                           there became occasions when a
                                           tapeid variable was needed in
                                           a function in which there were
                                           no structures declared that 
                                           knew the size of tape_id.    */

/* Structure for parse list. */
struct  PRS_LST
    {
    char    name[CMLEN+1];          /* image name */
    double  window[4];              /* image window */
    long    bands[MAXBND+1];        /* image bands */
    long    proj[MAXBND+1];         /* image bands */
    long    qflag[2];               /* window and band flags */
    };


int antoi
(
    register char *srce,        /* source string to convert     */
    register long *sz,          /* number of chars to convert   */
    register long *value        /* decimal value to return      */
);

void apendmsg
(
    char *msg,          /* current message containing mount prompt      */
    int loc             /* one char code from tapeloc field of user catalog */
);

int basenm
(
    char *hnam          /* pointer to host file name                  */
);

void bldname
(
    char *fname,        /* current host filename (i)    */
    char *newspec,      /* new host name specified (i)  */
    char *newname       /* new filename (no wildcard)(o)*/
);

int bldtpls
(
    char *catfile       /* name of user's tapelist file */
);

int c_cmhost
(
    char *tae,          /* Input TAE file name                       */
    char *hname         /* Output host file name                     */
);

int c_cmhostcs
(
    char *tae,          /* Input TAE file name                       */
    char *hname         /* Output host file name                     */
);

int c_cmtae
(
    char *hostname,     /* input host file name                      */
    char *tname         /* output TAE file name                      */
);

int c_getfiles
(
    char *filespec,          /* Host file specification for files to find  */
    char file_names[][CMLEN],/* Array of file names to be output          */
    long *nfiles,            /* Number of files found matching filespec    */
    long maxfiles,           /* Max number of files to return at one time  */
    long first_call,         /* Flag for first call to this function (0,1) */
    long ck_lock             /* 1 = return list of unlocked files only,    */
);

int c_getfiles
(
    char *filespec,          /* Specification of desired files to find     */
    char file_names[][CMLEN],/* Names of all the files found              */
    long *nfiles,            /* Number of files found                      */
    long maxfiles,           /* Max number of files to report              */
    long startf,             /* Flag for first call or not (0,1)           */
    long ck_lock             /* 1 = return list of unlocked files only,    */
);

int cmlock
(
    char *hname,             /* host name of file to open    */
    int *fd,                 /* file descriptor              */
    int flags,               /* file access                  */
    int maxtime              /* max time to wait (seconds)   */
);

int cmtaedir
(
    char *hostdir,           /* host directory specification */
    char *taedir             /* TAE directory specification  */
);

int copyfile
(
    char *hname,             /* host name of existing file   */
    char *hnewname,          /* host name of new file        */
    long report              /* report file?                 */
);

int deletefile
( 
    char *hname,             /* host name of file    */
    char *type,              /* "file" or "dir"      */
    long verify,             /* verify? (yes=1)      */
    long report,             /* report it? (yes=1)   */
    long chkflg              /* Check flag for on/off line*/
);

int dirlist
(
    char *nameoffile,        /* pointer to temp file of file names */
    char *filespec,          /* input spec                         */
    char file_ary[][CMLEN],  /* array of file names to be output   */
    long *nfiles,            /* number of files generated in list  */
    long maxfiles,           /* max number to list at a time       */
    long startf              /* flag for first call                */
);

int dirnm
(
    char *hnam               /* pointer to host path name                  */
);

int findfile
(
    char *filename,          /* filename to search for       */
    char *filestat			/* file status (single character) */
                            /* blank - stored on tape       */
                            /* D     - marked for Deletion  */
                            /* H     - Hidden               */
);

int get_assoc
(
    char *tmpnam,               /*ptr to name for temp file passed to fileutil*/
    char *host,                 /* pointer to host image name                 */
    char (*assoc_array)[CMLEN], /*pointer to array for files assoc with img   */
    long *nassoc,               /* pointer to number of associated files      */
    long wild,                  /* flag for wild carded file name             */
    long diskfiles              /* flag for displaying disk or tape files     */
);

int getcredate
(
    char *hostname,           /* Host file name               */
    char *credate             /* Date host file created       */
);

int getcredate
(
    char *hostname,
    char *credate
);

int getdat
(
    char *hostname,           /* Host file name               */
    char *strdate,            /* Date host file created       */
    char *timflg              /* flag for which time format, c or l */
);

int get_def_dir
(
    char *usercode,
    char *dir
);

int get_def_dir
(
    char *usercode,
    char *dir
);

int getdirs
(
    char *filespec,           /* host description of whole path and file name */
    char directories[][CMLEN],/* array of directory names that are found      */
    int *ndirs,               /* number of directories found                  */
    int maxdirs,              /* maximum number of directories to look for    */
    int startf                /* number of directories already found          */
);

int getmydef
(
    char *dir
);

int getmydef
(
    char *dir
);

int getprnm
(
    char *hostnm,             /* ptr to host file name        */
    long *nmlen               /* pointer length of file name allowed  */
);

char *gettoken
(
    char *string,
   char *delim
);

int gettplsnm
(
    char *hostn,              /* pointer to host file name              */
    char (*file_array)[CMLEN],/* pointer to array of input file list    */
    long *nfiles,
    long wild
);

int getver
(
    char *log_name,           /* logical name to be searched  */
    char *version             /* user name to be returned     */
);

int getver
(
    char *env_name,           /* environment variable to be searched*/
    char *version             /* user name to be returned     */
);

int isassoc
(
    char *hnam,               /* pointer to host file name                  */
    char *ext                 /* pointer to desired assoc file extension    */
);

int markstat
(
    char *filename,           /* filename to search for       */
    char *newstatus           /* new status of file in list   */
);

int marktpls
(
    char *hostn,              /* pointer to catalog entry name*/
    char *filestatus,         /* pointer to code for filestatus*/
    char *timestamp           /* pointer to file create date/time*/
);

int mkdr
(
    char *hostnm              /* pointer to directory path    */
);

int ops_bldtp
(
    char *catfile,            /* name of user's tapelist file */
    char *user                /* user name                    */
);

void processfile
(
    struct PARBLK *vblock,    /* ptr to TAE param block*/
    char *func,               /* function to perform  */
    char *name,               /* current filename     */
    char *newname,            /* new filename         */
    long imgset,              /* img & assoc? (1=yes) */
    long imageflg,            /* Process image file   */
    char aftype[][CMLEN],     /* Assoc file types     */
    long *acount,             /* # assoc file types   */
    long verify,              /* verify? (1=yes)      */
    long *nproc,              /* # files processed    */
    long chkflg,              /* Check flag for on/off line*/
    long copyflg              /* flag if image data is copied or renamed    */
);

int pttpls
(
    char *hosnm,              /* host name of file to store   */
    char *tplsrec,            /* pointer to buf for file info rec*/
    long *wflag               /* flag to write to file -  1   */
                              /*      or find file name - 0   */
);

int renamefile
(
    char *hname,              /* current host name    */
    char *hnewname,           /* new host name        */
    long verify,              /* verify? (1=yes)      */
    long report,              /* report file? (1=yes) */
    long copyflg              /* flag if image data is copied or renamed    */
);

int srchtpls
(
    char *hostroot,           /* pointer to host img name root(no extension)*/
    char file_array[][CMLEN], /* array for assoc files found in catalog     */
    long file_no_list[],      /*array of file_no fields of assoc_array files*/
    long *nfiles,             /* pointer to number of files to loop through */
    char tapeid[][TPIDLEN],   /* ptr to the tape on which these files exist */
    char tploc[]              /* pointer to array of tape locations   */
);

int unlock
(
    int *fd,                  /* file descriptor              */
    int maxtime               /* max time to wait (seconds)   */
);

int uphist
(
    struct PARBLK *vblock,    /* tae parameter block                  */
    char *hosthis             /* newly ouput his file of same name    */
);

int uptpls
(
    char *catf,               /* user catalog file name       */
    char *tapeid,             /* tape name-tapelist rec name  */
    long *dens,               /* tape density                 */
    float *feet,              /* pointer to inches of tape used */
    float *totft,             /* pointer to inches of tape available*/
    long *blk_count,          /* total blocks written to tape */
    long *file_no             /* sequence of this file on tape*/
);

int findtpls 
(
    FILE        *catf,        /* tapelist file  (input)               */
    char        *hname,       /* host image name  (input)             */
    char        *tplsnam,     /* taplist record  (output)             */
    char        *tplsrec,     /* taplist record  (output)             */
    long        *nbytes       /* num bytes read to reach name (output)*/
);

BOOL isdir
(
   char *path                 /* host name            */
);

int isimage
(
    char *taen               /* pointer to tae name                        */
);

int ops_tpop
(
    FILE        **catfd,      /* file descriptor for tapelist file    */
    long        nflag,        /* numeric flag for cmlock call         */
    char        *cflag,       /* pointer to string flag for cmlock call*/
    char        *user         /* user name                            */
);

int rdtplsrc
(
    FILE        **catfd,      /* file descriptor for tapelist file    */
    char        *buf,         /* pointer to char buf to place record in*/
    long        len           /* length of record to read             */
);

int tplsclos
(
    FILE        **catfd       /* file descriptor for tapelist file    */
);

int tplsopen
(
    FILE        **catfd,      /* file descriptor for tapelist file    */
    long        nflag,        /* numeric flag for cmlock call         */
    char        *cflag        /* pointer to string flag for cmlock call*/
);

#endif
