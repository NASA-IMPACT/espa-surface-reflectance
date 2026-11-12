#ifndef DITTO_H
#define DITTO_H

struct TAPE_FD			/* TAPE DESCRIPTOR FOR DITTO DEVICE           */
 	{
	char	drive[20];	/*	tape drive name			      */
	long	unit;		/*	unit number tape is opened to         */
	long	access;		/*      type of access (read, write, update)  */
	long	blksiz;		/*      block size (physical record size)     */
	long	recsiz;		/*      record size (logical record size)     */
	long    rectyp;		/*      record type (fixed or variable)       */
	long    cc;		/*      type carriage control                 */
	long	cur_rec;	/*      current record tape is positioned at  */
  	long	cur_file;	/*      current file tape is positioned at    */
	};

struct DISK_FD			/* DISK DESCRIPTOR FOR DITTO DEVICE           */
	{
	char	fname[CMLEN];	/*	disk file name                        */
	FILE	*fp;		/*      file pointer to opened disk file      */
	long	fd;		/* 	file descriptor to opened disk file   */
	long	recsiz;		/*	record size (logical record size)     */
	long	blksiz;		/*	block size (physical record size)     */
	long 	rectyp;		/* 	record type (fixed or variable)       */
	long 	cc;		/*	type of carriage control              */
	long	cur_rec;	/*	current record file is positioned at  */
	long	cur_file;       /*      current file disk is positioned at    */
				/*	= 0  ---> fname                       */
				/*      > 0  ---> fname_x                     */
	};

struct DITTO_DEV		/* DITTO descriptor of input and output device*/
	{
	struct TAPE_FD *in_tape;
	struct TAPE_FD *out_tape;
	struct DISK_FD *in_disk;
	struct DISK_FD *out_disk;
	char           *input_dev;
	char	       *output_dev;
	};
	

#define DCREATE 0		/* create new disk file       (write)         */
#define DREAD	1		/* read from disk file        (read)          */
#define DWRITE	2		/* write to disk file         (write)         */
#define DAPPEND 3		/* append to end of disk file (write)         */
#define DUPDATE 4		/* update disk file           (read & write)  */

#define NONE	 0		/* no carriage control                        */
#define PRINT	 1		/* print carriage control                     */
#define FORTRAN  2		/* fortran carriage control                   */
#define CARRIAGE 3		/* carriage control                           */

#define NOVAL   -1		/* no value specified by user                 */
#define VAR	 0		/* variable length records                    */
#define FIX      1		/* fixed length records                       */

#define ERROR_SEV	-3	/* severe error encountered                   */
#define ERROR_VERIFY    -2	/* error encountered - user is queried        */
#define ERROR_WARN	-1	/* warning encountered                        */


/* 
*   CONSTANTS USED FOR THE DITTO ERROR ROUTINE
*/
#define UPDATE_FD	0	/* update the DITTO I/O descriptor            */
#define NO_UPDATE_FD    1	/* don't update the DITTO I/O descriptor      */

#define DITTO_DOUBLE	100	/* device has already been allocated          */
#define DITTO_DSTAT	101	/* error retrieving disk file attributes      */
#define DITTO_DFIX	102	/* disk file does not have fixed length rec   */
#define DITTO_DVAR	103	/* disk file does not have variable length rec*/
#define DITTO_TPOPEN	104	/* error mounting tape                        */
#define DITTO_TPCLOS	105	/* error closing tape                         */
#define DITTO_TPREOP	106	/* error reopening tape                       */
#define DITTO_NALLOC	107	/* device was not allocated                   */
#define DITTO_ALLOC	108	/* device is not allocated                    */
#define DITTO_DOPEN	109	/* error opening disk file                    */
#define DITTO_CALLOC	110	/* error allocating dynamic buffer space      */
#define DITTO_DREAD	111	/* error reading from disk file               */
#define DITTO_DCLOSE	112	/* error closing disk file                    */
#define DITTO_TPREAD	113	/* error reading from tape                    */
#define DITTO_TALLOC	114	/* device is not allocated to tape            */
#define DITTO_BSF	115	/* error back spacing files                   */
#define DITTO_BSR	116	/* error back spacing records                 */
#define DITTO_FSF	117	/* error forward spacing files                */
#define DITTO_FSR	118	/* error forward spacing records              */
#define DITTO_TPREW	119	/* error rewinding device                     */
#define DITTO_TPWRIT	120	/* error writing to tape                      */
#define DITTO_DWRITE	121	/* error writing to disk file                 */
#define DITTO_RECTYP	123	/* record types are not compatiable           */
#define DITTO_RECSIZ	124	/* logical record sizes are not compatiable   */
#define DITTO_BLKSIZ    125	/* physical record sizes are not compatiable  */
#define DITTO_DALLOC	126	/* device is not allocated to a disk file     */
#define DITTO_VERIFY    127	/* verify error -- files are not the same     */
#define DITTO_SWAP      128	/* byte swap error encountered                */
#define DITTO_ABORT     129	/* user abort                                 */
#define DITTO_EOF	130	/* EOF detected prior to processing records   */
#define DITTO_HEXDIG    131     /* too many hex digits for byte value         */
#define DITTO_HEXINV    132	/* invalid hex digit specified                */
#define DITTO_PADCHR    133	/* more than one pad character specified      */
#define DITTO_INEOF     134	/* input EOF detected but not output          */
#define DITTO_OUTEOF	135	/* output EOF detected but not input          */
#define DITTO_UNEQ	136	/* unequal record sizes read                  */
#define DITTO_MODLEN    137	/* requested to modify beyond end of record   */
#define DITTO_INVMOD	138	/* user specified to modify an invalid record */
#define DITTO_MODREC	139	/* copy-modrec command run in interactive mode*/
#define DITTO_WEOF	140	/* error writing EOF mark		      */


int dopen
(
    struct DISK_FD *fd,         /* DITTO disk file descriptor        */
    long  amode                 /* type of access mode               */
);

int dclose
(
    struct DISK_FD *fd          /* DITTO disk file descriptor        */
);

int get_disk_fd
(
    char *which,                /* DITTO input/output device         */
    struct DISK_FD *disk        /* disk descriptor                   */
);

void put_disk_fd
(
    char *which,                /* DITTO input/output device         */
    struct DISK_FD *disk        /* disk descriptor                   */
);

int dread
(
    struct DISK_FD *fd,         /* DITTO disk file descriptor        */
    unsigned char *buf,         /* buffer to read data into          */
    long nbytes                 /* number of bytes to be read        */
);

int dskip
(
    struct DISK_FD *fd,         /* DITTO disk file descriptor        */
    long *count                 /* number of records to skip         */
);

int dwrite
(
    struct DISK_FD *fd,         /* DITTO disk file descriptor        */
    unsigned char *buf,         /* buffer to write from              */
    long nbytes                 /* number of bytes to be written     */
);

int dstat
(
	struct DISK_FD *fd, /* DITTO disk file descriptor        */
	long  mode          /* access mode                       */
);

void outmsg
(
    char *msg,				/* message to be printed             */
    long *print,			/* destination of the message        */
    FILE *lp_fp 			/* file pointer to line printer      */
);

void ditto_error
(
    char *which,                /* which device (input or output)    */
    long err_code,              /* error code                        */
    long severity,              /* error severity                    */
    long update                 /* attribute update flag             */
);

void hex2byt
(
    char *string,               /* character string of hex digits    */
    unsigned char *byte_val     /* hex string converted to byte value*/
);

int reopen_tape
(
    char *which,                /* DITTO input/output device         */
    long rewopt,                /* rewind option on tape open        */
    long access,                /* I/O access to the tape            */
    struct TAPE_FD *tape        /* tape descriptor                   */
);

void close_tape
(
    char *which,                /* DITTO input/output device         */
    struct TAPE_FD *tape,       /* tape descriptor                   */
    long rewopt,                /* tape rewind option                */
    long eofopt                 /* tape end of file option for close */
);

#endif
