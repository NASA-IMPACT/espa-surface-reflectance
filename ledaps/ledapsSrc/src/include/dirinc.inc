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



/* TDM CHECKOUT FILE_TIME=31-OCT-1983 15:57 DUA1:[TAEV1.INC]DIRINC.INC;1 */
/* TDM CHECKOUT FILE_TIME=11-MAY-1983 10:13 DUA0:[TAEV1.INC]DIRINC.INC;11 */
/*TDM         CHECKOUT FILE_TIME=28-APR-1983 16:06 DMA1:[TAEV1.INC]DIRINC.INC;6 */
/*TDM         CHECKOUT FILE_TIME=30-JAN-1983 13:22 DMA1:[TAEV1.INC]DIRINC.INC;1 */

/*
 * INCLUDE FILE FOR THE DIRECTIVE PACKAGE 
 *
 *
 *	CHANGE LOG:
 *
 *	28-apr-83	additional fields for .INCLUDE implementation...dm
 *	11-may-83	Delete .saveline and .saveflag from DIRBLK...dm
 *	27-oct-83	implement conditional directives...dm
 *	22-apr-85	Add field .remctx to access remote include files...dm
 *	05-nov-86	Eliminated #ifdef TAE_RCJM...nhe
 *
 */

#ifndef  I_DIRBLK			/* prevent double include */
#define  I_DIRBLK	0

#include	"taeconf.inp"
#include	"tmhost.inp"
#include	"tminc.inc"


    struct  DIRBLK			/* directive block structure */
	{
	struct  SFILE	*sfileptr;	/* pointer to SFILE of input file   */
 	TEXT	libr[FLIBRSIZ+1];	/* current  library 		    */
	struct  POSCTX  incpos;		/* posctx of current include file   */
	struct  SFILE	incfile;	/* SFILE for current include file   */
	TEXT		**sys_char;	/* system characteristics	    */
	COUNT		char_cnt;	/* number of char. strings	    */
	BOOL		inblock;	/* true if within if/else block	    */
	BOOL		c_valid;	/* true if condition already valid  */
	struct  CONTXT	*remctx;	/* context block for remote include */
	};


/*
 *...........    STATUS CODES FOR DIRECTIVE PACKAGE   ............ */
#define		DIRINC_BASE	2000
  
#define		D_EOT	 2001		/* end of text */
#define		D_PAGE	 2002		/* .PAGE within text */  
#define		D_EOF 	 2003		/* end of file */
#define 	D_ERRINC 2004		/* error in include file open */
#define		D_ERROR  2005		/* all other errors */

#endif
