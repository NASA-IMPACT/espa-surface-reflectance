/*
   include file for ansi tape information needed by Catalog Manager
*/
#ifdef vms
#include		<file.h>
#include		<time.h>
#else
#include		<setjmp.h>
#include		<stdio.h>
#include		<sys/file.h>
#include		<sys/time.h>
#ifdef _IBMR2
#include                <time.h>
#endif
#endif

#include		"header.h"		/* header label layout	*/
#include		"trailer.h"		/* trailer label layout	*/
#include		"volume.h"		/* volume label layout	*/

#define	IGNORE		void			/* ignore return value	*/
#define LABELSZ		80			/* size of all labels	*/
#define	MIN(a,b)	((a) < (b) ? (a) : (b))	/* macro for minimum	*/
#define	OK		(0)			/* good exit status	*/
#define	VOIDE					/* function that does	*/
						/* not return a value	*/

#define	BOF		1			/* processing headers	*/
#define	DATA		2			/* processing data	*/
#define	SWITCH		1			/* time to reelswitch	*/
#define	BUFSZ		65535			/* tape buffer size	*/

#ifdef	GLOBALS
#define	DECLARE				/* global definition		*/
#else
#define	DECLARE		extern		/* global declaration		*/
#endif
#undef DECLARE
