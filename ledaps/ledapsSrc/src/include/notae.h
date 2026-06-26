#ifndef NOTAE_H
#define NOTAE_H

#include "getargv.h" /* struct PARAMETER definition */
#include "worgen.h"  /* CMLEN definition */


void c_getfval
(
    char *valuestring,	/* space delimited string of floats	*/
    long max,		/* maximum number of values to extract	*/
    long *numspec,	/* actual number ov values extraced	*/
    float *output	/* array of floats extracted		*/
);

void c_getival
(
    char *valuestring,	/* space delimited string of integers	*/
    long max,		/* maximum number of values to extract	*/
    long *numspec,	/* actual number ov values extraced	*/
    long *output	/* array of integers extracted		*/
);

int c_getsval
(
    char *name,		/* parameter name to search for		*/
    char list[][CMLEN],	/* list of parameter specifications	*/
    char *value		/* parameter value(s)			*/
);

int c_parsepar
(
    long argc,			/* IN: argument count			*/
    char *argv[],		/* IN: argument list			*/
    struct PARAMETER *parms,
    char *scmdlist[],
    char list[][CMLEN]		/* OUT: list of input parsed parameters	*/
);

int check_name
(
    char *str,		/* IN: string to test for match		*/
    char *list[],	/* IN: list of strings for comparison	*/
    long numlist,	/* IN: number of strings in list	*/
    char *name		/* OUT: full matching name from list	*/
);

int check_type
(
    struct PARAMETER *parms,	/* IN: parameter structure		*/
    char list[][CMLEN],		/* IN: list of parsed parameters	*/
    long numparm		/* number of parameters			*/
);

int valdens
(
    long density
);

#endif
