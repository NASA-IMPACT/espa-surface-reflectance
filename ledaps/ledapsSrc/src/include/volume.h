/*
 *	ANSI X3.27-1978  American National Standard Magnetic
 *		Tape Labels and File Structure for Information Interchange
 *			pp 13-22
 */

/*
 *	Info for User-Volume-Label #1
 */
typedef struct
	{
	char	host_id[16];		/* host system identifier 	*/
	char	cm_version[30];		/* CM version number  		*/
	char	tae_version[30];	/* TAE verstion number 		*/
	} UVL1INFO;

/*
 *	Info for User-Volume-Label #1
 */
typedef struct
	{
	char	init_dtime[11];		/* tape initialization date and time */
					/* YYDDDHHMMSS,e. g. 89181135216*/
	char	density[4];		/* recording density  		*/
	char	loc_chain[2];		/* location in tape chain 	*/
					/*   "HD" - head tape	  	*/
					/*   "CO" - continuation tape	*/
	char	tape_no[2];		/* number of this tape in the chain */
	char	back_chain[6];		/* tape chain back pointer:	*/
					/* previous tape label		*/
	char	storetime[11];		/* date/time this tape was last */
					/* written to			*/
					/*  YYDDD, e. g. 89181		*/
	char	blk_count[6];		/* number of blocks written	*/
	char	feet[9];		/* number of feet of tape utilized */
	char	totft[9];		/* total number of feet available */
	char	uvl2_spaces[33];	/* reserved - blanks		*/
	} UVL2INFO;

/*
 *	Volume-Header #1 Layout
 */
struct VOL1
	{
	char	label_id[3];		/* "VOL" 			*/
	char	label_no;		/* "1" 				*/
	char	volume_ident[6];	/* volume identifier 		*/
	char	accessibility;		/* accessability 		*/
	char	vol1_spaces1[26];	/* spaces 			*/
	char	owner_ident[14];	/* owner identifier 		*/
	char	vol1_spaces2[28];	/* spaces 			*/
	char	label_version;		/* label-standard version - "3" */
	} ;

/*
 *	User Volume #1 Layout
 */
struct UVL1
	{
	char	 label_id[3];		/* "UVL" 			*/
	char	 label_no;		/* "1", "2", ... "9" 		*/
	UVL1INFO info;
	} ;

/*
 *	User Volume #2 Layout
 */
struct UVL2
	{
	char	 label_id[3];		/* "UVL" 			*/
	char	 label_no;		/* "1", "2", ... "9" 		*/
	UVL2INFO info;
	} ;

/*
 *	End-of-Volume #1 Layout
 */
struct EOV1
	{
	char	 label_id[3];		/* "EOV" */
	char	 label_no;		/* "1" */
	UVL1INFO info;
	} ;
/*
 *	End-of-Volume #2 Layout
 */
struct EOV2
	{
	char	 label_id[3];		/* "EOV" */
	char	 label_no;		/* "2" */
	UVL2INFO info;
	} ;
/*
 *	End-of-Volume #n Layout
 */
struct EOVN
	{
	char	label_id[3];		/* "EOV" */
	char	label_no;		/* "3", "4", ... "9" */
	char	eovn_data[76];		/* reserved for system use */
	} ;

#define	LABEL_VERSION	'3'
