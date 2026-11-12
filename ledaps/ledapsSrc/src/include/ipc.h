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



/*****************************************************************************
 *
 *  ipc.h >>> Header for ipc.c
 *
 *
 *	CHANGE LOG
 *
 *  10-mar-93	Copied from palm's directory...swd
 *  23-apr-93	Added extern "C" block, include <kr_ansi.h>...swd
 *  03-jun-93	Changed kr_ansi.h to TAEProto.h, ARGS to _TAE_PROTO...kbs
 *
 *****************************************************************************/
#ifndef  Ipc__Hdr
#define  Ipc__Hdr 0

#include	<taeproto.h>

#ifdef __cplusplus
extern "C" {		// Tell C++ not to overload these functions.
#endif

/*  server side (event driven)  */

int Ipc_NewService _TAE_PROTO ((char *serviceName)); /* returns fd of service */

int Ipc_Accept _TAE_PROTO ((int serviceFd));	/* returns fd of new client */

int Ipc_Receive _TAE_PROTO ((int fd, void *buffer, int bufferSize));

void Ipc_DeleteService _TAE_PROTO ((int fd));


/*  client side (sequential)  	*/

int Ipc_Connect _TAE_PROTO ((char *serviceName));	/* returns zero if ok */

int Ipc_Send _TAE_PROTO ((int fd, void *buffer, int buffersize));

int Ipc_Close _TAE_PROTO ((int fd));
#ifdef __cplusplus
}			// end the extern "C" block
#endif



#endif
