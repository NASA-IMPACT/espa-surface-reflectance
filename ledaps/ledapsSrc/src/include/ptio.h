#ifndef PTIO_H
#define PTIO_H

struct PTDESC
   {
   long fptr;		/* File pointer */
   long acc;		/* Access mode:  (IREAD, IUPDATE, or IWRITE) */
   long (*routine)();	/* Pointer to read/write routine */
   char name[CMLEN];	/* Name of tie point file */
   };

#endif
