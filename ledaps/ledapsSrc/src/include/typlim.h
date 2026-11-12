/*  Data type limits 
--------------------*/
#define	MAXUTINY	255		/* 8 bit: 2^8 - 1                     */
#define MINUTINY	0		/* 8 bit: 0                           */

#define MAXWORD		32767		/* 16 bit: 2^15 - 1                   */
#define MINWORD		-32768		/* 16 bit: -(2^15)                    */

#define MAXUWORD	65535		/* 16 bit: 2^16 - 1                   */
#define MINUWORD	0		/* 16 bit: 0                          */

#define MAXULONG	(unsigned long)4294967295	/* 32 bit: 2^32 - 1   */
#define MINULONG	0				/* 32 bit: 0          */

/* Modified MINLONG to be the same on all systems.  -2147483648 will produce
   an overflow on most systems, if not all.  (-2147483647 - 1) will work
   on any system.
------------------------------------------------------------------------------*/
#define MAXLONG    2147483647        /* 32 bit: 2 to the 31st -1              */
#define MINLONG    (-2147483647 - 1) /* VAX compiler wouldn't except          */

#ifdef gould
#define MAXREAL	   7.2370e75         /* 32 bit floating point (gould)         */
#define MINREAL	  -7.2370e75         /* 32 bit floating point (gould)         */
#endif
#ifdef sun
#define MAXREAL	   3.402823e38       /* 32 bit floating point (sun)           */
#define MINREAL	  -3.402823e38       /* 32 bit floating point (sun)           */
#endif
#ifdef vms
#define MAXREAL	   1.7e38            /* 32 bit floating point (vax)           */
#define MINREAL	  -1.7e38            /* 32 bit floating point (vax)           */
#endif
#ifdef _IBMR2
#define MAXREAL	   3.402823e38       /* 32 bit floating point (ibm aix)       */
#define MINREAL	  -1.401298e45       /* 32 bit floating point (ibm aix)       */
#endif
#ifdef sgi
#define MAXREAL	   3.4028e38         /* 32 bit floating point (sgi)           */
#define MINREAL	  -3.4028e38         /* 32 bit floating point (sgi)           */
#endif
#ifdef DGUX
#define MAXREAL	   3.4028e38         /* 32 bit floating point (dgux)          */
#define MINREAL	  -1.1755e38         /* 32 bit floating point (dgux)          */
#endif
#ifdef aiws
#define MAXREAL	   3.402823e38       /* 32 bit floating point (ibm rt)        */
#define MINREAL	  -3.402823e38       /* 32 bit floating point (ibm rt)        */
#endif
#ifdef linux
#define MAXREAL    3.402823e38       /* 32 bit floating point (linux)         */
#define MINREAL   -3.402823e38       /* 32 bit floating point (linux)         */
#endif
