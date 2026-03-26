#ifndef TAPEIO_H
#define TAPEIO_H

/*
*       Include File History:
*       version         Date            Contrib         Change
*       -------         ----            -------         ------
*        1.0                                            Original development
*        ??             6/89            K.Gacke         Port to SUN 4
*        ??             1/90            B.Davis         mc68010 for sun2 only,
*                                                       mc68000 for 2 and 3
*         2             7/90            K.Gacke         Primos Defintions
*	  3	       11/90		K.Gacke		Tape config file
*	  4		3/91		K.Gacke		Add Silicon Graphics
*	  5		4/91		K.Gacke		Add data general
*         6             7/91		K.Gacke		Add network, add drive
*							type field in 
*							tape.config file
*	  7		8/92		K.Gacke		Modified TAPEDESC
*							to add tape_len,
*							tape_bpi, and tape_irg
*							fields.
*	  8		8/95		K.Gacke		Added CDR and IRIX53_
*							EXABYTE as device types
*	  9		8/95		K.Gacke		Added DG54_GEN as a
*							device type
*        10            11/95            K.Gacke         Modified for IRIX 5.3 
*							(Fugitsu 9-track)
*        11		1/96		K.Gacke		Solaris SCSI supported
*        12		9/96		K.Gacke		DLT supported
*	 13		2/97		K.Gacke		CD-R flags added
*	 14		9/97		K.Gacke		Increased MAXDRIV from
*							16 to 32 (sgs6 has more
*							than 16 devices)
*/
#ifndef FUNCTION
#define FUNCTION
#endif
#ifndef TRUE
#define TRUE -1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef NULLVAL
#define NULLVAL -1
#endif

#ifdef sun
#ifdef sparc
#define sun4
#endif
#ifdef mc68010
#define sun2
#endif
#ifdef mc68000
#define sun2
#define sun3
#endif
#endif

#define MTGRP "rpops"		/* group name to send mount messages on unix */

#define TAPE_SUCC       0
#define TAPE_FAIL       -1
#define TAPE_EOV        -2
#define TAPE_EOT        -3
#define TPOPEN_SKIP     -4
#define TPFSF_FSF       -5
#define TPOPEN_REW      -6
#define TPREW_REW       -7
#define TP_SKIP         -8
#define TPWRIT_WRITE    -9
#define TPREAD_READ     -10
#define TPWEOF_WEOF     -11
#define TPCLOS_EOF      -12
#define TPOPEN_OPEN     -13
#define TPREW_OPEN      -14
#define TPREW_OFFL      -15
#define TPCLOS_OFFL     -16
#define TPCLOS_CLOSE    -17
#define TPREW_CLOSE     -18
#define TPWRIT_SIZE     -19
#define TPBSR_BSR       -20
#define TPFSR_FSR       -21
#define TPBSF_BSF       -22
#define TPOPEN_NOAVAIL  -23
#define TPOPEN_MNT      -24
#define TPOPEN_MNTSZ    -25
#define TPOPEN_ALL      -26
#define TPOPEN_WLOCK    -27
#define TPOPEN_OFFLINE  -28
#define TPOPEN_OPS      -29
#define TPOPEN_OPSWR    -30
#define TPOPEN_OPSRD    -31
#define TPOPEN_SOCKET   -32
#define TPOPEN_ASN      -33
#define TPOPEN_DENS     -34
#define TPOPEN_DRNUM    -35
#define TPOPEN_BADENS   -36
#define TPSK_CLOSE      -37
#define TPSK_OPEN       -38
#define TAPE_SMAL       -39
#define TPOPEN_LHOST    -40
#define TPBSF_NO        -41		/* BSF not supported */
#define TPBSR_NO	-42		/* BSR not supported */
#define TPFSF_NO	-43		/* FSF not supported */
#define TPFSR_NO	-44		/* FSR not supported */
#define TPEOT_NO	-45		/* EOT may not be detected properly */
#define TPAPP_NO	-46		/* unable to append data */
#define TP_NONET	-47		/* tape i/o not built for network */
#define TAPE_BOF        -97
#define TAPE_BOT        -98
#define TAPE_EOF        -99

#define NOTOKT   44              /* abort status for VMS */
#define NOREW   0
#define REWIND  1
#define REWOFFL 2
#define NO_HANDLE 0
#define HANDLE 1
#define TREAD 0
#define TWRITE 1
#define TUPDATE 2

#define TOPEN   1
#define TREOPEN 2

/*
   Device types
*/
#define GOULD_TU77  	0	/* pn tape drive */
#define PRIME_GEN   	1	/* prime tape drive */
#define VAX_GEN     	2	/* vax tape drive  (use FSR to skip files) */
#define VAX_FSF        	3	/* vax tape drive  (use FSF to skip files) */
#define SUN_EXABYTE 	4	/* sun exabyte */
#define SUN_XY	    	5	/* sun Xylogics 472 */
#define SUN_TM	    	6	/* sun TapeMaster */
#define SUN_ST	    	7	/* sun 1/4" streamer */
#define SUN_SCSA_EX	8	/* sun SCSA exabyte */
#define IBMRS_EXABYTE   9	/* ibm rs6000 exabyte */
#define SGI_EXABYTE    10	/* SGI exabyte */
#define SGI_XY	       11	/* SGI Xylogics 772 1/2 inch */
#define SGI_SCSI       12	/* SGI generic SCSI 1/2 inch */
#define IRIX40_EXABYTE 13	/* SGI IRIX 4.0 exabyte */
#define IRIX40_SCSI    14	/* SGI IRIX 4.0 SCSI drives */
#define IRIX40_JAG     15	/* SGI IRIX 4.0 JAG (diff scsi) drives */
#define IRIX53_EXABYTE 16	/* SGI IRIX 5.3 exabyte */
#define IRIX53_SCSI    17	/* SGI IRIX 5.3 generic SCSI */
#define DG_GEN	       18	/* data general */
#define DG54_GEN       19	/* data general DGUX >= 5.4R3.10 */
#define SOLARIS_SCSI   20	/* Solaris SCSI */
#define CD	       21       /* CD reader */
#define CDR	       22       /* CD writer */
#define NMR_TAPE_DEV   23

/* the following constants are highly device dependent */
/* see selio/mt.h  ; the lower 16 bits of status are in dp->mt_dsreg*/
/* the high 16 bits are in dp->mt_erreg */

#define EOF_GOULD_TU77 0x20	/* pn */
#define EOT_GOULD_TU77 0x02
#define BOT_GOULD_TU77 0x04
#define MIN_GOULD_TU77 1
#define MAX_GOULD_TU77 64511

#define EOF_SUN_EXABYTE 0x0	/* sun exabyte */
#define EOT_SUN_EXABYTE 0x0
#define BOT_SUN_EXABYTE 0x0
#define MIN_SUN_EXABYTE 1
#define MAX_SUN_EXABYTE 65535

#define EOF_SUN_SCSA_EX 0x0	/* sun exabyte */
#define EOT_SUN_SCSA_EX 0x0
#define BOT_SUN_SCSA_EX 0x0
#define MIN_SUN_SCSA_EX 1
#define MAX_SUN_SCSA_EX 64512

#define EOF_SUN_XY 0x01		/* sun Xylogics */
#define EOT_SUN_XY 0x80
#define BOT_SUN_XY 0x40
#define MIN_SUN_XY 1
#define MAX_SUN_XY 32768

#define EOF_SUN_TM 0x80		/* sun TapeMaster */
#define EOT_SUN_TM 0x10
#define BOT_SUN_TM 0x20
#define MIN_SUN_TM 1
#define MAX_SUN_TM 32768

#define EOF_SUN_ST 0x80		/* sun 1/4" streamer */
#define EOT_SUN_ST 0x10
#define BOT_SUN_ST 0x20
#define MIN_SUN_ST 1
#define MAX_SUN_ST 8192

#define EOF_IBMRS_EXABYTE 0x0	/* ibm rs6000 exabyte */
#define EOT_IBMRS_EXABYTE 0x0
#define BOT_IBMRS_EXABYTE 0x0
#define MIN_IBMRS_EXABYTE 1
#define MAX_IBMRS_EXABYTE 65535

#define EOF_IRIX40_EXABYTE 0x8000 /* silicon graphics IRIX 4.0 exabyte */
#define EOT_IRIX40_EXABYTE 0x0
#define BOT_IRIX40_EXABYTE 0x0002
#define MIN_IRIX40_EXABYTE 1
#define MAX_IRIX40_EXABYTE 245760 

#define EOF_IRIX40_SCSI 0x8000	/* silicon graphics IRIX 4.0 SCSI drives */
#define EOT_IRIX40_SCSI	0x08
#define BOT_IRIX40_SCSI 0x02
#define MIN_IRIX40_SCSI 1
#define MAX_IRIX40_SCSI 204800

#define EOF_IRIX40_JAG 0x21	/* silicon graphics IRIX 4.0 JAG drives */
#define EOT_IRIX40_JAG	0x08
#define BOT_IRIX40_JAG 0x02
#define MIN_IRIX40_JAG 1
#define MAX_IRIX40_JAG 204800

#define EOF_SGI_EXABYTE 0x8000	/* silicon graphics IRIX 3.xx exabyte */
#define EOT_SGI_EXABYTE 0x0
#define BOT_SGI_EXABYTE 0x0002
#define MIN_SGI_EXABYTE 1
#define MAX_SGI_EXABYTE 131072

#define EOF_SGI_XY      0x0	/* silicon graphics Xylogics 772 1/2 inch */
#define EOT_SGI_XY	0x02
#define BOT_SGI_XY 	0x0
#define MIN_SGI_XY 	1
#define MAX_SGI_XY 	65535

#define EOF_SGI_SCSI    0x8000	/* silicon graphics IRIX 4.0 SCSI drives */
#define EOT_SGI_SCSI	0x08
#define BOT_SGI_SCSI 	0x02
#define MIN_SGI_SCSI 	1
#define MAX_SGI_SCSI 	511987

#define EOF_DG_GEN 0x08		/* data general tapes */
#define EOT_DG_GEN 0x10
#define BOT_DG_GEN 0x04
#define MIN_DG_GEN 1
#define MAX_DG_GEN 65535

#define EOF_SOLARIS_SCSI 0x0	/* sun solais */
#define EOT_SOLARIS_SCSI 0x0
#define BOT_SOLARIS_SCSI 0x0
#define MIN_SOLARIS_SCSI 1
#define MAX_SOLARIS_SCSI 245000

#define EOF_VAX_GEN 0x0		/* vax */
#define EOT_VAX_GEN 0x0
#define BOT_VAX_GEN 0x0
#define MIN_VAX_GEN 14
#define MAX_VAX_GEN 65535

#define EOF_PRIME_GEN 0x0	/* prime */
#define EOT_PRIME_GEN 0x0
#define BOT_PRIME_GEN 0x0
#define MIN_PRIME_GEN 1
#define MAX_PRIME_GEN 12286  /* 6144 half words minus 2 */

#define MINTAPEWT 1
#define MAXTAPEWT 65535 

#define MNTSZ 10
#define MAXDRIV 32

#define DEN_DEF  0
#define DEN_800  1
#define DEN_1600 2
#define DEN_6250 3
#define DEN_LO   4
#define DEN_MED  5
#define DEN_HI   6
#define DENSCT   7              /*  number of different tape densities */

/*
*       Define default tape lengths for various type tapes.  Measured in feet.
*/
#define TPLEN_9TRACK    2400
#define TPLEN_3480       600
#define TPLEN_3490      1100
#define TPLEN_8MM        365
#define TPLEN_4MM        270
#define TPLEN_DLT       1200
#define TPLEN_D3	1100


struct TAPEDESC
	{
	int  bot;		/* beginning of tape status bits */
	int  dev_type;		/* magtape device type */
	int  eof_flg;		/* boolean flag, TRUE=>write EOF on close */
	int  eof;		/* end of file status bits */
	int  eot;		/* end of tape status bits */
	int  err;		/* system error code */
        int  fd;		/* tape descriptor */
	int  flags;		/* tape drive flags */
        int  inout;		/* i/o access */
	int  maxbyte;		/* max # bytes per physical record */
	int  minbyte;		/* min # bytes per physical record */
	int  nd;		/* network descriptor */
	int  netflg;		/* network flag */
	float tape_bpi;		/* tape bpi */
	float tape_irg;		/* tape inter record gap */
	float tape_len;		/* length of tape available (inches) */
	char density[8];	/* tape density */
        char dtype[CMLEN];	/* drive type (9-track, 8mm, etc) */
	char lhost[CMLEN];	/* local host name */
	char lname[CMLEN];	/* logical magtape device name */
        char name[CMLEN];	/* physical magtape device name */
	char thost[CMLEN];	/* tape drive host name */
	};

#define E_EOL	-101		/* end of line for tape configuration file */

struct TPREC 
	{
	int  dev_type;
	int  flags;
	char devnam[CMLEN];
	char drive[CMLEN];
	char drvnmr[DENSCT][CMLEN];
	char host[CMLEN];
        char dtype[CMLEN];	/* drive type (9-track, 8mm, etc) */
	};

struct TMNT
    {
    char lname[CMLEN];		/* logical magtape device name */
    char host[CMLEN];		/* computer system name */
    char user[CMLEN];		/* user name */
    char name[CMLEN];		/* physical magtape device name */
    char dtype[CMLEN];		/* drive type (9-track, 8mm, etc) */
    };

#ifdef prime

fortran t$mt();

struct TPREQ
    {
    long userid;
    long tapnam;
    };

#define TP_WAIT       -2          /* wait for operator to reply to mnt req */
#define TPST_PARITY  0x8000       /* bit 1 */
#define TPST_RUNAWAY 0x4000       /* bit 2 */
#define TPST_CRC     0x2000       /* bit 3 */
#define TPST_LRC     0x1000       /* bit 4 */
#define TPST_SMALL   0x0800       /* bit 5 */
#define TPST_ERROR   0x0400       /* bit 6 */
#define TPST_WRITE   0x0200       /* bit 7 */
#define TPST_EOF     0x0100       /* bit 8 */
#define TPST_READY   0x0080       /* bit 9 */
#define TPST_ONLINE  0x0040       /* bit 10 */
#define TPST_EOT     0x0020       /* bit 11 */
#define TPST_REWIND  0x0010       /* bit 12 */
#define TPST_BOT     0x0008       /* bit 13 */
#define TPST_WRPROT  0x0004       /* bit 14 */
#define TPST_OVERRUN 0x0002       /* bit 15 */
#define TPST_RWDONE  0x0001       /* bit 16 */


#define TP_BSF     0x2440         /* back space one file */
#define TP_BSR     0x6440         /* back space one record */
#define TP_DEN1600 0x8040         /* set density to 1600 bpi */
#define TP_DEN6250 0x8050         /* set density to 6250 bpi */
#define TP_FSF     0x2480         /* forward skip one file */
#define TP_FSR     0x6480         /* forward skip one record */
#define TP_OFFLIN  0x8020         /* rewind tape and unload */
#define TP_READ    0x4180         /* read tape record */
#define TP_REWIND  0x0020         /* rewind tape to BOT */
#define TP_STATUS  0x8000         /* check drive status */
#define TP_WEOF    0x2490         /* write a file mark */
#define TP_WRITE   0x4190         /* write tape record */
#endif 

#ifndef MTNBSF
#define MTNBSF          11      /* backward space file to BOF */
#endif

#define NOBSF	0x01		/* tape drive does not handle BSF opcode */
#define NOBSR   0x02		/* tape drive does not handle BSR opcode */
#define NOFSF	0x04		/* tape drive does not handle FSF opcode */
#define NOFSR	0x08		/* tape drive does not handle FSR opcode */
#define NOEOT	0x10		/* tape drive does not handle EOT properly */
#define NOAPP	0x20		/* tape drive is not able to append data */
#define CDRPR   0x01		/* CD-R has a label printer attached */
#define CDR2X   0x02		/* CD-R is a 2x speed CD writer */
#define CDR6X   0x04		/* CD-R is a 6x speed CD writer */


int c_tapelabel
(
    char *thebuf[11],  /* IN - buffer of information for LIBRARY FORMAT   */
                       /*      [00] code for the tape                     */
                       /*      [01] density of the tape                   */
                       /*      [02] description of the tape               */
                       /*      [03] drive used                            */
                       /*      [04] project the tape was created for      */
                       /*      [05] requestor.  Initials of tape creator  */
                       /*      [06] sequence number of the tape           */
                       /*      [07] tapeid                                */
                       /*      [08] error flag                            */
                       /*      [09] unused                                */
                       /*      [10] unused                                */
                       /* IN - buffer of information for FREE FORMAT      */
                       /*      [0-10] text in free format                 */

    long retain,   /* IN - Months to retain (used in library format only) */
                   /*      0    = permanent storage in archive            */
                   /*      1-60 = number of months to retain in arhcive   */
                   /*      > 60 = reset to 60                             */
    long printflg, /* IN - Print flag                                     */
                   /*      FALSE = do not print                           */
                   /*      TRUE  = print                                  */
    long formflg,  /* IN - File format flag (0 = FREE, 1 = STANDARD)      */

    char tae_name[] /* OUT - TAE name of file containing the tape label    */
);

int c_tpansi
(
    long *unit,         /* (I/O) pointer to tape descriptor */
    long *mode,         /* (I)   error mode */
    char *dens,         /* (I)   tape density */
    long *inout,        /* (I)   tape access (TWRITE,TUPDATE,TREAD) */
    char *message,      /* (I/O) message to the operator */
    char *tapeid        /* (I/O) tape ID */
);

int c_tpbsf
(
    long *tapefd,       /*  tape file descriptor                */
    long *count,        /*  number of files to backspace        */
    long *mode          /*  Error mode                          */
);

int c_tpbsr
(
    long *tapefd,       /*  tape file descriptor                */
    long *count,        /*  number of records to backspace      */
    long *mode          /*  error mode                          */
);

int c_tpclos
(
    long *tapefd,  /* tape drive unit number */
    long *mode,    /* error message option:
                          = handle   = 1 = report errors to user
                          = nohandle = 0 = return error status to
                                           calling program */
    long *opt1,    /* rewind option:
                          = norew   = 0 = do not rewind tape 
                          = rewind  = 1 = rewind tape after close
                          = rewoffl = 2 = rewind tape drive and 
                                          take drive offline after close */
    long *opt2    /* end of file option:
                          = 1, 2 or 3 tape eof marks to write
                            at tape write device close. */
);

int c_tpfsf
(
    long *tapefd,       /*  tape file descriptor                */
    long *count,        /*  number of files to forward skip     */
    long *mode          /*  error mode                          */
);

int c_tpfsr
(
    long *tapefd,       /*  tape file descriptor                */
    long *count,        /*  number of records to forward skip   */
    long *mode          /*  error mode                          */
);

int c_tplblpr
(
    char *printer       /* output printer name                          */
);

long c_tplsrd
(
    long *unit,         /* tape unit number             (input)         */
    char *key,          /* input key of record to be read (in/out)      */
    long *clen,         /* number of char bytes to read (in/out)        */
    long *dlen,         /* number of data bytes to read (in/out)        */
    char *cbuf,         /* buffer to read char data into  (output)      */
    unsigned char *dbuf,/* data buffer to read data into                */
    char *dtype         /* datatype of data record (I2,B,R4,...) (out)  */
);

int c_tplswr
(
    long *unit,         /* tape unit number                     (input)       */
    char *key,          /* key of record to be written          (input)       */
    long *clen,         /* number of char bytes to be written   (input)       */
    long *dlen,         /* number of data bytes to be written   (input)       */
    char *cbuf,         /* character buffer to be written       (input)       */
    char *dbuf,         /* data buffer to be written            (input)       */
    char *dtype         /* datatype of data record (I2,B,R4,...)(input)       */
);

int c_tpopen
(
    long *unit,         /*  (O)   tape unit */
    long *skip,         /*  (I)   files to be skipped */
    long *mode,         /*  (I)   error mode */
    long *opt1,         /*  (I)   rewind flag */
    char *dens,         /*  (I)   tape density */
    long *inout,        /*  (I)   type of access */
    char *message,      /*  (I/O) message to the operator */
    char *tapeid        /*  (O)   tape id */
);

int c_tpread
(
    long *tapefd,          /*  tape unit number            */
    unsigned char *buffer, /*  buffer                      */
    long *bufsz,           /*  buffer size                 */
    long *mode             /*  error mode                  */
);

int c_tpreop
(
    long *unit,            /*  tape unit number                    */
    long *skip,            /*  number of files to skip             */
    long *mode,            /*  error mode                          */
    long *opt1,            /*  rewind flag                         */
    char *dens,            /*  tape density                        */
    long *inout,           /*  tape access                         */
    char *tapename         /*  tape drive name                     */
);

int c_tprew
(
    long *tapefd,          /*  tape unit number                    */
    long *mode,            /*  error mode                          */
    long *opt1             /*  rewind offline flag                 */
);

int c_tpweof
(
    long *tapefd,          /*  tape unit number                    */
    long *count,           /*  number of end-of-file marks to be written */
    long *mode             /*  error mode                          */
);

int c_tpwrit
(
    long *tapefd,           /*  Tape unit number                    */
    unsigned char *buffer,  /*  buffer                              */
    long *bufsz,            /*  buffer size -- number of bytes      */
    long *mode              /*  error mode                          */
);


#endif
