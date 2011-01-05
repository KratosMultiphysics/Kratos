
#include <stdio.h>
#include <stdlib.h>
#include "slu_zdefs.h"
#include "slu_util.h"


void
zreadtriple(int *m, int *n, int *nonz,
	    doublecomplex **nzval, int **rowind, int **colptr)
{
/*
 * Output parameters
 * =================
 *   (a,asub,xa): asub[*] contains the row subscripts of nonzeros
 *	in columns of matrix A; a[*] the numerical values;
 *	row i of A is given by a[k],k=xa[i],...,xa[i+1]-1.
 *
 */
    int    i, j, k, jsize, nz, lasta;
    doublecomplex *a, *val;
    int    *asub, *xa, *row, *col;
    
    /* 	Matrix format:
     *    First line:  #rows, #cols, #non-zero
     *    Triplet in the rest of lines:
     *                 row, col, value
     */

    scanf("%d%d%d", m, n, nonz);
#ifdef DEBUG
    printf("zreadtriple(): *m %d, *n %d, *nonz, %d\n", *m, *n, *nonz);
#endif
    zallocateA(*n, *nonz, nzval, rowind, colptr); /* Allocate storage */
    a    = *nzval;
    asub = *rowind;
    xa   = *colptr;

    val = (doublecomplex *) SUPERLU_MALLOC(*nonz * sizeof(doublecomplex));
    row = (int *) SUPERLU_MALLOC(*nonz * sizeof(int));
    col = (int *) SUPERLU_MALLOC(*nonz * sizeof(int));

    /* Read into the triplet array from a file */
    for (i = 0; i < *n+1; ++i) xa[i] = 0;
    for (nz = 0; nz < *nonz; ++nz) {
	scanf("%d%d%lf%lf\n", &row[nz], &col[nz], &val[nz].r, &val[nz].i);
	if (row[nz] < 0 || row[nz] >= *m || col[nz] < 0 || col[nz] >= *n) {
	    fprintf(stderr, "(%d, %d) out of bound!\n", row[nz], col[nz]);
	    exit (-1);
	}
	++xa[col[nz]]; /* Count number of nonzeros in each column */
    }

    /* Initialize the array of column pointers */
    k = 0;
    jsize = xa[0];
    xa[0] = 0;
    for (j = 1; j < *n; ++j) {
	k += jsize;
	jsize = xa[j];
	xa[j] = k;
    }
    
    /* Copy the triplets into the column oriented storage */
    for (nz = 0; nz < *nonz; ++nz) {
	j = col[nz];
	k = xa[j];
	asub[k] = row[nz];
	a[k] = val[nz];
	++xa[j];
    }

    /* Reset the column pointers to the beginning of each column */
    for (j = *n; j > 0; --j)
	xa[j] = xa[j-1];
    xa[0] = 0;

    SUPERLU_FREE(val);
    SUPERLU_FREE(row);
    SUPERLU_FREE(col);

#ifdef CHK_INPUT
    for (i = 0; i < *n; i++) {
	printf("Col %d, xa %d\n", i, xa[i]);
	for (k = xa[i]; k < xa[i+1]; k++)
	    printf("%d\t%16.10f\n", asub[k], a[k]);
    }
#endif

}
