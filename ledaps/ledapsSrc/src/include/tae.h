#ifndef LASTAE_H
#define LASTAE_H

/**********************************************************************
 *
 *  This is an aggregate include file for TAE applications.
 *   
 ********************************************************************** */

#include "taeconf.inp"
#include "parblk.inc"
#include "terminc.inc"

void c_errmsg
(
    char *message,	/* Message associated with the error	      */
    char *key,		/* LAS error key, null string, "UNIX" or      */
    long *vrity		/* Severity:				*/
    				/*    <= 0   fatal			*/
	    			/*     > 0   non-fatal			*/
);


/* NOTE: the string returned by squeeze should be free'd by the calling
         routine to prevent a memory leak */
char *squeeze
(
   char *str,
   long len
);
 
/* TEMPORARY TAE function prototypes
--------------------------*/
long t_init
(
  long *,
  long *,
  long *
);

void t_pos
(
  long,
  long
);

void t_output
(
  long,
  long,
  char *
);

int t_read
(
    char *,
    int *
);

int t_write
(
  char *,
  long
);

int t_bell(void);

void t_clear
(
  void
);

void t_input
(
  long,
  long,
  char *,
  long *
);

struct VARIABLE *p_find
(
  struct PARBLK *,
  char *
);

void p_inifor
(
  struct PARBLK *
);


int p_forp
(
   struct PARBLK *,
   char *,
   int
);


int q_intg(struct PARBLK *p, TEXT name[], FUNINT count, 
           TAEINT intg[], FUNINT mode);

int q_string(struct PARBLK *p, TEXT name[], FUNINT count, 
             TEXT *vector[], FUNINT mode);

int q_out(struct PARBLK *p);

#endif
