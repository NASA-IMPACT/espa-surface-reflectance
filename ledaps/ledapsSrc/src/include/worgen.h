/*  This is a general WORLD include file.  The constants within this include 
    file are used by many if not all the LAS support and application routines
------------------------------------------------------------------------------*/
/* Boolean true and false 
-------------------------*/
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* define NULL 
--------------*/
#ifndef NULL
#define NULL 0
#endif

/* SUCCESS or FAIL return status,
   End of File, & Buffer to small
---------------------------------*/
#define E_SUCC 0
#define E_FAIL -1
#define E_EOF  -2
#define E_SMAL -3

/* character string lengths
---------------------------*/
#define ERRLEN 256
#define CMLEN 256

/* valid data types 
-------------------*/
#define EBYTE 1
#define EWORD 2
#define ELONG 3
#define EREAL 4
#define EDOUBLE 5
#define ECHAR 6
#define EUWORD 7
#define EULONG 8
#define NDTYPES 8

/* sizes of data types
----------------------*/
#define EBYTESIZE (sizeof(unsigned char))
#define EWORDSIZE (sizeof(short))
#define ELONGSIZE (sizeof(long))
#define EREALSIZE (sizeof(float))
#define EDOUBLESIZE (sizeof(double))
#define ECHARSIZE (sizeof(char))
#define EUWORDSIZE (sizeof(unsigned short))
#define EULONGSIZE (sizeof(unsigned long))

/* file access types 
--------------------*/
#define IREAD 0
#define IWRITE 1
#define IUPDATE 2

/*  define the key word FUNCTION  -- This is placed at the beginning of each 
    function declared in 'C' in order to designate the functions from each other
------------------------------------------------------------------------------*/
#ifndef FUNCTION		/*  if not defined			      */
#define FUNCTION		/*  define FUNCTION			      */
#endif

/*  Define unix for IBM
-----------------------*/
#ifdef _IBMR2
#define unix
#endif

