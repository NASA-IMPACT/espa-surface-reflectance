#include "poly_coeff.h"

/* LU decomposition functions, taken from the Wikipedia page for
   LU decomposition */

/* INPUT: A - array of pointers to rows of a square matrix having dimension N
          Tol - small tolerance number to detect failure when the matrix is
                near degenerate
   OUTPUT: Matrix A is changed, it contains both matrices L-E and U as
           A=(L-E)+U such that P*A=L*U.
           The permutation matrix is not stored as a matrix, but in an
           integer vector P of size N containing column indexes where the
           permutation matrix has "1". */
static int LUPDecompose(float **a, int n, float tol, int *p)
{
    int i, j, k, imax;
    float maxa, *ptr, absa;

    for (i = 0; i < n; i++)
        p[i] = i; // Unit permutation matrix, p[n] initialized with n

    for (i = 0; i < n; i++)
    {
        maxa = 0;
        imax = i;

        for (k = i; k < n; k++)
            if ((absa = fabs(a[k][i])) > maxa)
            {
                maxa = absa;
                imax = k;
            }

        if (maxa < tol)
            return 0;   //failure, matrix is degenerate

        if (imax != i)
        {
            //pivoting p
            j = p[i];
            p[i] = p[imax];
            p[imax] = j;

            //pivoting rows of a
            ptr = a[i];
            a[i] = a[imax];
            a[imax] = ptr;
        }

        for (j = i + 1; j < n; j++)
        {
            a[j][i] /= a[i][i];

            for (k = i + 1; k < n; k++)
                a[j][k] -= a[j][i]*a[i][k];
        }
    }

    return 1;  //decomposition done
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
static void LUPSolve(float **a, int *p, float *b, int n, float *x)
{
    int i, j;

    for (i = 0; i < n; i++)
    {
        x[i] = b[p[i]];

        for (j = 0; j < i; j++)
            x[i] -= a[i][j]*x[j];
    }

    for (i = n - 1; i >= 0; i--)
    {
        for (j = i + 1; j < n; j++)
            x[i] -= a[i][j]*x[j];

        x[i] = x[i]/a[i][i];
    }
}


/*****************************************************************
MODULE:   get_3rd_order_poly_coeff
*****************************************************************/
void get_3rd_order_poly_coeff
(
    float aot[], /* I: array of AOT values [NAOT_VALS] */
    float *atm,  /* I: atmospheric variable, corresponding to aerosol vector */
    int nAtm,    /* I: dimension of atm array, max NAOT_VALS  */
    float *coeff /* O: regression coefficient, quadratic [NCOEF] */
)
{
    float y1[NCOEF] = {0.0, 0.0, 0.0, 0.0};
    float x[NAOT_VALS][NCOEF];
    float z[NCOEF][NCOEF];
    float *zptr[NCOEF];
    int i,j,k;
    int p[4];  /* permutation vector */

    for (i=0; i<nAtm; i++)
    {
        x[i][0] = aot[i] * aot[i] * aot[i];
        x[i][1] = aot[i] * aot[i];
        x[i][2] = aot[i];
        x[i][3] = 1.0;
    }
  
    for (i=0; i<NCOEF; i++)
    {
        for (j=0; j<NCOEF; j++)
        {
           z[i][j] = 0.0;
           for (k=0; k<nAtm; k++)
               z[i][j] += x[k][i] * x[k][j];
        }

        zptr[i] = z[i];
    }

    for (i=0; i<NCOEF; i++)
        for (j=0; j<nAtm; j++)
            y1[i] += x[j][i] * atm[j];

    LUPDecompose(zptr, NCOEF, 1e-10, p);
    LUPSolve(zptr, p, y1, NCOEF, coeff);
}
