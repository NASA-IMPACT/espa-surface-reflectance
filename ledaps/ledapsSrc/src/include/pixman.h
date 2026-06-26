#ifndef PIXMAN_H
#define PIXMAN_H

#define TOOBIG		77		/* Output value to large for data type*/
#define TOOSMALL	88		/* Output value to small for data type*/
#define BIGSMALL	95		/* Both TOOBIG and TOOSMALL errors
					   were encountered in array	      */
#define TYPERR		66		/* not a correct data type	      */

#include "typlim.h"			/* include the data type limits	      */


/* c_pxadd.c
---------------------------------------------------------------------*/
int c_pxadd
(
   long *ns,      /* I:  number of samples                              */
   long *in1,     /* I:  input array one                                */
   long *in2,     /* I:  input array two                                */
   long *out,     /* O:  output array                                   */
   float *scale1, /* I:  scale applied to in1 before addition           */
   float *scale2, /* I:  scale applied to in2 before addition           */
   float *offset, /* I:  offset added in after sum                      */
   long *dtype    /* I:  data type (EBYTE,EWORD,ELONG,EREAL)            */
);


/* c_pxbtan.c
---------------------------------------------------------------------*/
int c_pxbtan
(
   long *in,   /* I:  input array                                       */
   long *map,  /* I:  input mapping array                               */
   long *out,  /* O:  output "anded" array                              */
   long *ns,   /* I:  number of samples                                 */
   long *dtype /* I:  data type (EBYTE,EWORD,ELONG,EREAL)               */
);


/* c_pxbtor.c
---------------------------------------------------------------------*/
int c_pxbtor
(
   long *in,   /* I:  input array                                        */
   long *map,  /* I:  input mapping array                                */
   long *out,  /* O:  output "ored" array                                */
   long *ns,   /* I:  number of samples                                  */
   long *dtype /* I:  data type (EBYTE,EWORD,ELONG,EREAL)                */
);


/* c_pxbxor.c
---------------------------------------------------------------------*/
int c_pxbxor
(
   long *in,   /* I:  input array                                        */
   long *map,  /* I:  input mapping array                                */
   long *out,  /* O:  output "xored" array                               */
   long *ns,   /* I:  number of samples                                  */
   long *dtype /* I:  data type (EBYTE,EWORD,ELONG,EREAL)                */
);


/* c_pxcomp.c
---------------------------------------------------------------------*/
int c_pxcomp
(
   long *in,   /* I:  input array                                        */
   long *out,  /* O:  output array                                       */
   long *ns,   /* I:  number of samples                                  */
   long *dtype /* I:  data type (EBYTE,EWORD,ELONG,EREAL)                */
);


/* c_pxconv.c
---------------------------------------------------------------------*/
int c_pxconv
(
    long *ftyp,             /* I: from type                              */
    long *totyp,            /* I: to type                                */
    unsigned char *frombuf, /* I: from buffer                            */
    unsigned char *tobuf,   /* O: to buffer                              */
    long *size              /* I: number of samples in the input buffer  */
);


/* c_pxcopy.c
---------------------------------------------------------------------*/
void c_pxcopy
(
    unsigned char *frombuf, /* I: source buffer of the data to be copied */
    unsigned char *tobuf,   /* O: destination buffer of the data         */
    long *dtype,            /* I: data type of the buffers               */
    long *ns                /* I: number of elements or pixels           */
);


/* c_pxdiv.c
---------------------------------------------------------------------*/
int c_pxdiv
(
   long *ns,      /* I: number of samples                                 */
   long *in1,     /* I: input array one                                   */
   long *in2,     /* I: input array two                                   */
   long *out,     /* O: output array                                      */
   float *scale,  /* I: scale applied to input arrays after mult.         */
   float *offset, /* I: offset added in after sum                         */
   long *dtype,   /* I: data type (EBYTE,EWORD,ELONG,EREAL)               */
   long *errval   /* I: value to place in array where div zero occurs     */
);


/* c_pxe10.c
---------------------------------------------------------------------*/
int c_pxe10
(
   long *in,   /* I: Input array                                          */
   long *out,  /* O: output array                                         */
   long *ns,   /* I: Number of elements in the input array                */
   long *dtype /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxe255.c
---------------------------------------------------------------------*/
int c_pxe255
(
   long *in,   /* I: Input array                                          */
   long *out,  /* O: output array                                         */
   long *ns,   /* I: Number of elements in the input array                */
   long *dtype /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxe32k.c
---------------------------------------------------------------------*/
int c_pxe32k
(
   long *in,   /* I: Input array                                          */
   long *out,  /* O: output array                                         */
   long *ns,   /* I: Number of elements in the input array                */
   long *dtype /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxfill.c
---------------------------------------------------------------------*/
void c_pxfill
(
   unsigned char *tobuf,  /* O: the buffer to be filled with fillval      */
   long *dtype,           /* I: data type (EBYTE,EWORD,ELONG,EREAL)       */
   long *ns,              /* I: number of elements or pixels to be filled */
   unsigned char *fillval /* I: fill value for the buffer                 */
);


/* c_pxhiss.c
---------------------------------------------------------------------*/
int c_pxhiss
(
   long *ns,      /* I: pointer to number of samples in buffer             */
   long *bufptr,  /* I: 4 byte buffer pointer                              */
   long *datatyp, /* I: pointer to the type of buffer                      */
   long *histo,   /* O: pointer to a 4 byte histogram array                */
   long *count,   /* O: total number of values the histogramed             */
   long *skip     /* I: sampling step                                      */
);


/* c_pxhist.c
---------------------------------------------------------------------*/
int c_pxhist
(
   long *ns,      /* I: pointer to number of samples in buffer             */
   long *bufptr,  /* I: 4 byte buffer pointer                              */
   long *datatyp, /* I: data type (EBYTE,EWORD,ELONG,EREAL)                */
   long *histo,   /* O: pointer to a 4 byte histogram array                */
   long *count    /* O: umber of values the histogram array contains       */
);


/* c_pxindr.c
---------------------------------------------------------------------*/
int c_pxindr
(
   long *in,    /* I: input array                                          */
   long *val,   /* I: value to be searched for                             */
   long *ns,    /* I: number of samples                                    */
   long *index, /* O: index to be returned                                 */
   long *dtype  /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxinx.c
---------------------------------------------------------------------*/
int c_pxinx
(
   long *in,    /* I: input array                                          */
   long *val,   /* I: value to be searched for                             */
   long *ns,    /* I: number of samples                                    */
   long *index, /* O: index to be returned                                 */
   long *dtype  /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxinxd.c
---------------------------------------------------------------------*/
int c_pxinxd
(
   long *in,    /* I: input array                                          */
   long *val,   /* I: value to be searched for                             */
   long *ns,    /* I: number of samples                                    */
   long *index, /* O: index to be returned                                 */
   long *dtype  /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxinxr.c
---------------------------------------------------------------------*/
int c_pxinxr
(
   long *in,    /* I: input array                                          */
   long *val,   /* I: value to be searched for                             */
   long *ns,    /* I: number of samples                                    */
   long *index, /* O: index to be returned                                 */
   long *dtype  /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxl10.c
---------------------------------------------------------------------*/
int c_pxl10
(
   long *in,   /* I: input array                                          */
   long *out,  /* O: output array                                         */
   long *ns,   /* I: number of samples                                    */
   long *dtype /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxl255.c
---------------------------------------------------------------------*/
int c_pxl255
(
   long *in,   /* I: input array                                          */
   long *out,  /* O: output array                                         */
   long *ns,   /* I: number of samples                                    */
   long *dtype /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxland.c
---------------------------------------------------------------------*/
int c_pxland
(
   long *in,   /* I: input array                                          */
   long *map,  /* I: input mapping array                                  */
   long *out,  /* O: output "anded" array                                 */
   long *ns,   /* I: number of samples                                    */
   long *dtype /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxline.c
---------------------------------------------------------------------*/
int c_pxline
(
   long *in,     /* I: input array                                          */
   long *out,    /* O: output array                                         */
   long *ns,     /* I: number of samples                                    */
   long *dtype,  /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
   float *scale, /* I: scale factor                                         */
   float *offset /* I: offset to be added                                   */
);


/* c_pxlkup.c
---------------------------------------------------------------------*/
int c_pxlkup
(
   long *in,     /* I: input array                                          */
   long *out,    /* O: output array                                         */
   long *ns,     /* I: number of samples                                    */
   long *dtype,  /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
   long *lut,    /* I: lookup table array                                   */
   long *luttype /* I: data type of lookup table                            */
);


/* c_pxlor.c
---------------------------------------------------------------------*/
int c_pxlor
(
   long *in,   /* I: input array                                          */
   long *map,  /* I: input mapping array                                  */
   long *out,  /* O: output "anded" array                                 */
   long *ns,   /* I: number of samples                                    */
   long *dtype /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxlxor.c
---------------------------------------------------------------------*/
int c_pxlxor
(
   long *in,   /* I: input array                                          */
   long *map,  /* I: input mapping array                                  */
   long *out,  /* O: output "xored" array                                 */
   long *ns,   /* I: number of samples                                    */
   long *dtype /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxmax.c
---------------------------------------------------------------------*/
int c_pxmax
(
   long *in,   /* I: input array                                          */
   long *map,  /* I: input mapping array                                  */
   long *out,  /* O: output  array                                        */
   long *ns,   /* I: number of samples                                    */
   long *dtype /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxmin.c
---------------------------------------------------------------------*/
int c_pxmin
(
   long *in,   /* I: input array                                          */
   long *map,  /* I: input mapping array                                  */
   long *out,  /* O: output  array                                        */
   long *ns,   /* I: number of samples                                    */
   long *dtype /* I: data type (EBYTE,EWORD,ELONG,EREAL)                  */
);


/* c_pxmnmx.c
---------------------------------------------------------------------*/
int c_pxmnmx
(
   long *in,     /* I: input array                                         */
   long *ns,     /* I: number of samples                                   */
   double *dmin, /* I/O: minimum pixel value                               */
   double *dmax, /* I/O: maximum pixel value                               */
   long *dtype   /* I: data type (EBYTE,EWORD,ELONG,EREAL)                 */
);


/* c_pxmul.c
---------------------------------------------------------------------*/
int c_pxmul
(
   long *ns,      /* I: number of samples                                  */
   long *in1,     /* I: input array one                                    */
   long *in2,     /* I: input array two                                    */
   long *out,     /* O: output array                                       */
   float *scale,  /* I: scale applied to input arrays after multiplication */
   float *offset, /* I: offset added in after sum                          */
   long *dtype    /* I: data type (EBYTE,EWORD,ELONG,EREAL)                */
);


/* c_pxsamp.c
---------------------------------------------------------------------*/
void c_pxsamp
(
   long *inns,        /* I: number of samples in input buffer              */
   long *dtype,       /* I: data type (EBYTE,EWORD,ELONG,EREAL)            */
   long *inc,         /* I: sample increment to subsample by               */
   unsigned char *buf /* I/O: buffer to be subsampled                      */
);


/* c_pxsub.c
---------------------------------------------------------------------*/
int c_pxsub
(
   long *ns,      /* I: number of samples                                   */
   long *in1,     /* I: input array one                                     */
   long *in2,     /* I: input array two                                     */
   long *out,     /* O: output array                                        */
   float *scale1, /* I: scale applied to input array one before subtraction */
   float *scale2, /* I: scale applied to input array two before subtraction */
   float *offset, /* I: offset added in after difference                    */
   long *dtype    /* I: data type (EBYTE,EWORD,ELONG,EREAL)                 */
);


/* c_pxswap.c
---------------------------------------------------------------------*/
int c_pxswap
(
   unsigned char *buf, /* I/O: data buffer which is to be byte swapped     */
   long *ns,           /* I: number of entities in the buffer              */
   long *size          /* I: number of bytes per entity                    */
);


/* c_pxsys.c
---------------------------------------------------------------------*/
int c_pxsys
(
   long *insys,        /* I: input computer system                         */
   long *outsys,       /* I: output computer system                        */
   unsigned char *buf, /* I/O: array of real numbers to be converted       */
   long *dtype,        /* I: data type (EBYTE,EWORD,ELONG,EREAL)           */
   long *size,         /* I: number of elements in array                   */
   void *under,        /* I: value assigned if underflow detected          */
   void *over,         /* I: value assigned if overflow detected           */
   long *flag          /* I: flag used for error detection                 */
);


#endif
