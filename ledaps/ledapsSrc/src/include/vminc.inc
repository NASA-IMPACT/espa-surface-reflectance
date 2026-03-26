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



/*
 *	include file for vm_ functions
 *
 *	CHANGE LOG:
 *                23-APR-87 Created...tpl
 *
 */

#include	"taeconf.inp"	/* TAE configuration (REQUIRED)		*/
#include	"parblk.inc"	/* parameter block definitions		*/

#ifndef       I_VMINC
#define       I_VMINC   0

    struct    VM_STRUCT {
                        struct NP_PARBLK    npblk;
                        };
#endif
