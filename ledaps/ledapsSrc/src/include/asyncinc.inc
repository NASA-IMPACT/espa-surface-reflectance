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



/* >>VAX/VMS<< */
/********  INCLUDE FILE FOR ASYNCHRONOUS JOB PROCESSING    *************
 *
 *	CHANGE LOG:
 *
 *	30-jan-85	Initial release (extracted from TMINC)...dm
 *	22-feb-85	Add fields .rpcb and .killed to ACB ...dm
 *	04-mar-85	Add MO_TUTHLP message code...nhe
 *	11-mar-85	Add .m_parblk as union to MON_MSG struct...dm
 *	13-mar-85	Changed sfi and skey to m_sfi, m_skey etc...dm
 *	16-apr-85	Add 'path_num' field  to rpcb...dm
 *	19-apr-85	Change acb.proc to acb.pdf to save pdf spec(needed
 *			for dynamic tutoring of an async process)...dm
 *	26-apr-85	Add .remlibr to acb, change MO_CONT to MO_SYNC...dm
 *	14-jun-85	Update for UNIX implementation...dm
 *	26-aug-85	Add MO_DYNEXIT for dynamic tutor exit msg code...dab
 *	15-jul-86	Add ASYNC_PROCESS defn...palm
 *	16-jul-86	Add .aprocess, .child_job, and .parent_job to ACB...palm
 *	10-dec-86	Correct UNIX-only define to "...child_pid"...peb
 *	23-may-90	Remove all RCJM stuff...ljn
 *
 */

#ifndef  I_ASYNC
#define  I_ASYNC	0			/* prevent double definition */

#include	"taeconf.inp"			/* TAE configuration defs    */
#include	"tminc.inc"				/* TAE standard definitions  */
#include	"parblk.inc"			/* parameter block defs	     */
#include	"fileinc.inp"			/* file package definitions  */
#include	"eventinc.inp"			/* event control definitions */

/*  MON_MSG
 *
 *  Used for monitor-to-monitor messages (async).  Note that all message types
 *  use the same structure and there are no unions; we may resort to unions in 
 *  the future to save space.
 */
    struct 	MON_MSG			/* message from one monitor to another */
    	{
    	CODE	type;			/* see below			*/

	union
	    {
	    struct
		{
#ifdef UNIX
		TAEINT  child_pid;		/* child's process id	*/
#endif
    		TEXT	filespec[FSPECSIZ+1];	/* name of save file 	*/
		CODE	parmode;		/* Parent's run_type	*/
		TAEINT	sfi;			/* $SFI value 		*/
		TEXT	skey[STRINGSIZ+1];	/* $SKEY 		*/
		TEXT	parent_job[JOBNAMESIZ+1];  /* used at init	*/
		TEXT	child_job[JOBNAMESIZ+1]; 
		} ma_msg;
	    } mu_asy;
    	};

/* redefine MON_MSG members with unique names
*/
 
#define  m_parblk	mu_asy.parblk
#define  m_childpid	mu_asy.ma_msg.child_pid 
#define  m_filespec	mu_asy.ma_msg.filespec    
#define  m_parmode 	mu_asy.ma_msg.parmode  
#define  m_sfi		mu_asy.ma_msg.sfi	
#define  m_skey		mu_asy.ma_msg.skey
#define  m_rempdf	mu_asy.ma_msg.rempdf
#define  m_parent_job   mu_asy.ma_msg.parent_job
#define  m_child_job    mu_asy.ma_msg.child_job

/* define async monitor message types
 */
#define MO_INIT		1		/* initiation message (to spawned monitor) */
#define MO_REQUEST	2		/* request message (to parent monitor)  */
#define MO_REPLY	3		/* reply to request (to spawned)	*/
#define MO_TERM		4		/* job termination message (either way)	*/
#define MO_STATUS	5		/* status message (to parent)		*/
#define MO_SYNC		6		/* handshake (sync) msg (either way)	*/
#define MO_GETPAR	7		/* getpar request msg (to parent mon)   */
#define MO_EXIT		8		/* parent agent terminate (either way) */
#define MO_ABORT	9		/* job abort interrupt msg (to agent)		*/
#define MO_START	10		/* req to establish a remote path (to agent)	*/
#define MO_PROCEXEC	11		/* request for proc execution (to agent)	*/
#define MO_GETFILE	12		/* request for a remote file (to agent)	*/
#define MO_DYNEXIT	13		/* exit from dyn tut (to spawned monitor)	*/

#define ASYNC_PROCESS	(ASYNC+100)	/* for job_type arg to create_async   */

/* NOTE: MO_ABORT sent from parent tm to agent over the interrupt channel only 
 */

/* 	Async Job Control Block Structure					    */

    struct	ACB		
	{
	struct 	ACB	*link;		/* forward link			    */
	TEXT	name[JOBNAMESIZ+1];	/* TAE job name			    */
	struct  FSBLOCK	 pdf;		/* pdf associated with the proc	    */
	BOOL	active;			/* true if this job is active	    */
	BOOL	waiting;		/* true if waiting for REPLY	    */
    	BOOL	done;			/* true if monitor has exited	    */
	BOOL	notify;			/* true if user to be told on jobend*/
	TAEINT	sfi;			/* final $SFI value		    */
	TEXT	skey[STRINGSIZ+1];	/* final $SKEY			    */
	COUNT 	seq;			/* sequence no. on last SHOW rep't  */
    	CODE	act_mode;		/* AC_TUTOR or AC_CMD		    */
    	TEXT	requ_file[HNAMESIZ+1];	/* name of par file for a dyn request */
    	TEXT	preface[MAXSTRSIZ];	/* preface for tutor on dyn request */
    				/* TBD: make these pointers to dyn storage  */
	CODE	wait_num;		/* ticket number for request waits  */
	unsigned remote;		/* TRUE if remote job		    */
	struct  ATCB *atcb;		/* ptr to corresponding ATCB	    */
	BOOL	aprocess;		/* true: ASYNC-PROCESS type job	    */
	};

#define  ACBSIZ  sizeof(struct ACB)


/* 	Async Task Control Block Structure				    */

    struct	ATCB		
	{
    	TEXT	session_id[PNAMESIZ+1];	/* session id of new session	    */
    	GENPTR  inmsg;			/* ptr to msg from spawned monitor  */
    	struct  TCB	*tcb;		/* points to tcb		    */
	int	pid;			/* VMS process ID		    */
	int	send_chan;		/* send mailbox channel		    */
	};

#define ATCBSIZ sizeof(struct ATCB)

/* 	path control block	*/

        struct  RPCB				/* remote path control block */
	{
	struct  RPCB	*link;			/* link to next entry	     */
	TEXT  	node[NODENAMESIZ+1];		/* remote node name	     */
	TEXT	user[FINAMESIZ+1];		/* remote login user name    */
	TEXT	path_id[4];			/* unique path id as 'pnn'   */
	BOOL	active;				/* TRUE if a job active	     */
	struct  ATCB	*atcb;			/* pointer to atcb for agent */
	};

#endif
