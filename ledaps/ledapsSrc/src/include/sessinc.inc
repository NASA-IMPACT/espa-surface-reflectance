/****************************************************************************
 *      Copyright (c) 1993, 1994
 *      Century Computing, Inc.
 *      ALL RIGHTS RESERVED
 *
 *      The software (programs, data bases and/or documentation) on or in
 *      any media can not be reproduced, disclosed, or used except under
 *      the terms of the TAE Plus Software License Agreement.
 *
 ***************************************************************************/



/*
 * Session log specific definitions & structures.
 *
 * CHANGE LOG:
 *
 */

#ifndef	I_SESSINC		/* prevent double include	*/
#define I_SESSINC	0

/*
 * define the struct SESSLOG for session logging .
 */
    struct	  SESSLOG
	{
 	BOOL		exist;			/* file existence flag */
	BOOL		enable;			/* logging enabled flag */
	struct  SFILE   sfile;			/* SFILE for session log file */
	};

#endif
