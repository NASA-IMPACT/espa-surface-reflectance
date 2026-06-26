#ifndef DISKIO_H
#define DISKIO_H

/*
 * Return status values from disk handling routines
 */
#define	DISK_EXISTS	(-2L)			/* file exists		  */
#define	DISK_GONE	(-3L)			/* file does not exist	  */
#define	DISK_NOSPACE	(-4L)			/* no space for file	  */
#define MAXPATH		50			/* maximum number of paths*/

#ifdef gould
#define UNIX_DISK				/* define unix disk	  */
#endif
#ifdef sparc
#define UNIX_DISK				/* define unix disk	  */
#endif
#ifdef sgi
#define UNIX_DISK				/* define unix disk	  */
#endif
#ifdef DGUX
#define UNIX_DISK				/* define unix disk	  */
#endif

/* structure for determining what disk is presently being reserved	  */

struct DISKRESERVE
{
	char filename[CMLEN];		/* name of file being reserved	  */
	char active;			/* "A" active "D" deleted	  */
	int required_space;		/* amount of disk space reseved	  */
	int ppid;			/* pid of reserving process	  */
#if defined (DGUX) || defined (sgi)
	int dev;			/* device id of the filesystem */
#else
	short dev;			/* device id of the filesystem */
#endif
};

int c_dkbuf
(
    long *bufsize       /* size of buffer */
);

long c_dkcopy
(
    char *from,         /* current file name */
    char *to,           /* new file name */
    char directories[][CMLEN]  /* list of destination directories */
);

long c_dkcre
(
    char directories[][CMLEN],  /* array of directories to search for space */
    char *filename,             /* name of file to be created */
    long filesize,              /* size of file to be created */
    long wait                   /* wait > 0, wait for space */
);

long c_dkdel
(
    char *filename      /* name of file to be deleted */
);

long c_dkfree
(
    char *path          /* pathname to unreserve space for*/
);

long c_dkmove
(
    char *from,         /* name of file to be moved */
    char *to,           /* name of output file */
    char directories[][CMLEN] /* list of destination directories */
);

long c_dkresv
(
    char *path,         /* pathname to reserve space for */
    int size            /* Size of the file to reserve space for*/
);

double c_dkrspac
(
    char *path          /* pathname to determine free space          */
);

long c_dkusage
(
    char *path,         /* pathname to check usage thereof */
    long *usage_counter /* Utilization counter for path */
);

#endif
