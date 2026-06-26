/*
 * Warning, this file was automatically created by the TIFF configure script
 * VERSION:	 v3.5.2
 * DATE:	 Sat Feb 24 23:43:49 2001
 * TARGET:	 i586-unknown-Cygwin1.1.8(0.34/3/2)
 * CCOMPILER:	 /usr/bin/gcc-2.95.2-6 19991024 (cygwin experimental)
 */
#ifndef _PORT_
#define _PORT_ 1
#ifdef __cplusplus
extern "C" {
#endif
#include <sys/types.h>
#define HOST_FILLORDER FILLORDER_LSB2MSB
#define HOST_BIGENDIAN	0
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
typedef double dblparam_t;
#ifdef __STRICT_ANSI__
#define	INLINE	__inline__
#else
#define	INLINE	inline
#endif
#define GLOBALDATA(TYPE,NAME)	extern TYPE NAME
#ifdef __cplusplus
}
#endif
#endif
