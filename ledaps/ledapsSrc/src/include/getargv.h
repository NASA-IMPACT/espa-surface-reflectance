#ifndef LASGETARGV_H
#define LASGETARGV_H

/*
  These definitions will be global and included in the local function
  getargv() and the support functions such as c_parsepar() and c_getsval().

  The enumerated types are supposed to represent the allowed TAE types.
  The structure defined here will be initialized in each local getargv().
--------------------------------------------------------------------------*/
enum PTYPE { NULLP, INTEGER, FLOAT, STRING };


struct PARAMETER {
    char **plist;
    long *pcount;
    enum PTYPE *ptype;
} PARAMETER;

#endif
