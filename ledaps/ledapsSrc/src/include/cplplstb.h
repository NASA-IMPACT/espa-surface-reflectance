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


// These avoid C++ runtime library references.
// CHANGE LOG
// 01-oct-92	Initial...ljn/rt
// 19-oct-92	link1 should be linkl...rt
// 04-apr-94	Changed filename to cplpl_stubs.h so it's <=14 characters...dag
extern "C"
{
void 	_entry() {}			// dummy C++ entry point (OASYS)
int 	__main() {return 1;}		// dummy C++ entry point (GNU)
void 	_main() {}			// dummy C++ entry point (cfront)
}
struct __linkl *__head;			// dummy C++ global (cfront)
