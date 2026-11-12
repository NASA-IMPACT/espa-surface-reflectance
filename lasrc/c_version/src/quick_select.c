/*****************************************************************************
FILE: quick_select.c

PURPOSE: Contains functions for handling the quick_select algorithm for finding
the median of a floating point array.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
This Quickselect routine is based on the algorithm described in
  "Numerical recipes in C", Second Edition,
  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
  This code by Nicolas Devillard - 1998. Public domain.
*****************************************************************************/
#include <quick_select.h>

/******************************************************************************
MODULE:  swap

PURPOSE:  Swaps two items in the the floating point array.

RETURN VALUE:
Type = N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
void swap
(
    float *x,   /* I/O: first value to swap */
    float *y    /* I/O: second value to swap */
)
{
    float temp;    /* temporary storage location */

    temp = *x;
    *x = *y;
    *y = temp;
}


/******************************************************************************
MODULE:  quick_select

PURPOSE:  Routine to return the median of the input array.

RETURN VALUE:
Type = N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
  1. The input array is modified/sorted during this process.
******************************************************************************/
float quick_select
(
    float *arr,  /* I/O: input floating point array to be sorted; output
                         sorted floating point array */
    int n        /* I: number of elements in the array */
)
{
    int low, high;   /* starting and ending elements of array to be sorted */
    int median;      /* median index of the input array */
    int middle;      /* middle of the current indices */
    int ll, hh;      /* low and high indices of the current array */

    /* Set up the array low, high, and median indices */
    low = 0;
    high = n - 1;
    median = (low + high) / 2;

    /* Loop until the array is sorted and median is found */
    for (;;)
    {
        /* One element array */
        if (high <= low)
            return arr[median];

        /* Two element array */
        if (high == low + 1)
        {
            if (arr[low] > arr[high])
                swap (&arr[low], &arr[high]);
            return arr[median];
        }

        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if (arr[middle] > arr[high])
            swap (&arr[middle], &arr[high]);
        if (arr[low] > arr[high])
            swap (&arr[low], &arr[high]);
        if (arr[middle] > arr[low])
            swap (&arr[middle], &arr[low]);

        /* Swap low item (now in position middle) into position (low+1) */
        swap (&arr[middle], &arr[low+1]);

        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        for (;;)
        {
            do ll++; while (arr[low] > arr[ll]);
            do hh--; while (arr[hh]  > arr[low]);

            if (hh < ll)
                break;

            swap (&arr[ll], &arr[hh]);
        }

        /* Swap middle item (in position low) back into correct position */
        swap (&arr[low], &arr[hh]);

        /* Re-set active partition */
        if (hh <= median)
            low = ll;
        if (hh >= median)
            high = hh - 1;
    }
}
