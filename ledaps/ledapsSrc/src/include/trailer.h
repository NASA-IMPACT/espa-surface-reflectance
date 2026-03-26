/*
 *	ANSI X3.27-1978  American National Standard Magnetic
 *		Tape Labels and File Structure for Information Interchange
 *			pp 13-22
 */

/*
 *	End-of-File #1 Layout
 */
struct EOF1
	{
	char	label_id[3];		/* "EOF" */
	char	label_no;		/* "1" */
	INFO1	info;
	} ;
/*
 *	End-of-File #2 Layout
 */
struct EOF2
	{
	char	label_id[3];		/* "EOF" */
	char	label_no;		/* "2" */
	INFO2	info;
	} ;
/*
 *	End-of-File #n Layout
 */
struct EOFN
	{
	char	label_id[3];		/* "EOF" */
	char	label_no;		/* "3", "4", ... "9" */
	char	eofn_data[76];		/* reserved for system use */
	} ;
/*
 *	User Trailer Labels
 */
struct UTLA
	{
	char	label_id[3];		/* "UTL" */
	char	label_no;		/* "a" character */
	char	utla_data[76];		/* reserved for user application */
	} ;
