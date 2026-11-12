#ifndef LASMATH_H
#define LASMATH_H

double c_randnm
(
    long *seed                /* seed value to be reused */
);

void ggud
(
    long *seed,
    long k,
    long nr,
    long ir[]
);

void linv1f
(
    float a[],
    long *n,
    long *ia,
    float ainv[]
);

void linv2p
(
    double *a,
    long *n,
    double *ainv,
    float det[2]
);

void c_qrdecomp
(
    double *A,      /* Matrix A (stored by columns)    */
    long m,         /* Number of rows in matrix A      */
    long n,         /* Number of columns in matrix A   */
    double *v,      /* Vector v                        */
    long flag       /* if True, order matrix A by columns */
);

void c_qrsolve
(
    double *A,      /* Matrix A (stored by columns)    */
    long m,         /* Number of rows in matrix A      */
    long n,         /* Number of columns in matrix A   */
    double *v,      /* Vector v                        */
    double *b,      /* Vector b                        */
    long iflag      /* flag (see comments)             */
);

void house
(
    long n,         /* Input--dimension of vectors x and v                  */
    double *x,      /* Input--vector x                                      */
    double *v,      /* Output--vector v                                     */
    double *vtv,    /* Output--Product of vector v transpose and v          */
    double *Px      /* Output--Product of Householder matrix P and vector x */
);

void row_house
(
    double *A,      /* Matrix A (stored by columns)                         */
    long m,         /* Number of rows in matrix A                           */
    long n,         /* Number of columns in matrix A                        */
    long mdim,      /* Row dimension                                        */
    double *v,      /* Vector v                                             */
    double vtv      /* Product of vector v transpose and vector v           */
);

double norm_2
(
    double *x,
    long n
);

void mdsti
(
    float *q,
    long *f,
    float *x
);

float ran2
(
    long *idum
);

void sort
(
    int n,
    long ra[]
);

void vabmxf
(
    float v[],      /* input,  vector                          */
    long *length,   /* input, length of the subset vector      */
    long *inc,      /* input, increment between successive     */
    long *index,    /* output, index of the element having the */
    float *vmax     /* output,maximum abs value                */
);

void vabsmf
(
    float v[],      /* input,  vector                          */
    long *length,   /* input, length of the subset vector      */
    long *inc,      /* input, increment between successive     */
    float *vsum     /* output, sum of abs values               */
);

void viprff
(
    float x[],      /* vector a                        */
    float y[],      /* vector y                        */
    long *n,        /* no. of elements in x or y       */
    long *ix,       /* increment between elements of x */
    long *iy,       /* increment between elements of y */
    float *c        /* inner product c                 */
);

void vmulff
(
    float a[],      /* matrix a                     */
    float b[],      /* matrix b                     */
    long *l,        /* no. rows in a                */
    long *m,        /* no. columns in a             */
    long *n,        /* no. columns in b             */
    long *ja,       /* row dimension of a           */
    long *jb,       /* row dimension of b           */
    float c[],      /* product matrix c             */
    long *jc        /* row dimension of c           */
);

void vmulfp
(
    float a[],      /* matrix a                     */
    float b[],      /* matrix b                     */
    long *l,        /* no. rows in a                */
    long *m,        /* no. columns in a             */
    long *n,        /* no. columns in b             */
    long *ja,       /* row dimension of a           */
    long *jb,       /* row dimension of b           */
    float c[],      /* product matrix c             */
    long *jc        /* row dimension of c           */
);

void vmulfs
(
    float a[],      /* matrix a                     */
    float b[],      /* symmetric matrix b stored in */
    long *l,        /* no. rows in a                */
    long *m,        /* no. columns in a             */
    long *ja,       /* row dimension of a           */
    float c[],      /* product matrix c             */
    long *jc        /* row dimension of c           */
);

void vmulsf
(
    float a[],      /* symmetric matrix a stored in */
    long *n,        /* order of a and no. rows in b */
    float b[],      /* matrix b                     */
    long *m,        /* no. columns in b             */
    long *jb,       /* row dimension of b           */
    float c[],      /* product matrix c             */
    long *jc        /* row dimension of c           */
);

void vmulss
(
    float a[],      /* matrix a in symmetric storage mode */
    float b[],      /* matrix b in symmetric storage mode */
    long *n,        /* order of a and b             */
    float c[],      /* product matrix c             */
    long *jc        /* row dimension of c           */
);

#endif
