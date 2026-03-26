#ifndef _QUICK_SELECT_H_
#define _QUICK_SELECT_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

void swap
(
    float *x,   /* I/O: first value to swap */
    float *y    /* I/O: second value to swap */
);

float quick_select
(
    float *arr,  /* I/O: input floating point array to be sorted; output
                         sorted floating point array */
    int n        /* I: number of elements in the array */
);

#endif
