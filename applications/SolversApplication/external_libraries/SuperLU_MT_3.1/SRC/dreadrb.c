/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file dreadrb.c
 * \brief Read a matrix stored in Rutherford-Boeing format
 *
 * <pre>
 * -- SuperLU routine (version 3.0) --
 * Lawrence Berkeley National Laboratory.
 * June 30, 2009
 * </pre>
 *
 * Purpose
 * =======
 *
 * Read a DOUBLE PRECISION matrix stored in Rutherford-Boeing format 
 * as described below.
 *
 * Line 1 (A72, A8)
 *      Col. 1 - 72   Title (TITLE)
 *      Col. 73 - 80  Matrix name / identifier (MTRXID)
 *
 * Line 2 (I14, 3(1X, I13))
 *      Col. 1 - 14   Total number of lines excluding header (TOTCRD)
 *      Col. 16 - 28  Number of lines for pointers (PTRCRD)
 *      Col. 30 - 42  Number of lines for row (or variable) indices (INDCRD)
 *      Col. 44 - 56  Number of lines for numerical values (VALCRD)
 *
 * Line 3 (A3, 11X, 4(1X, I13))
 *      Col. 1 - 3    Matrix type (see below) (MXTYPE)
 *      Col. 15 - 28  Compressed Column: Number of rows (NROW)
 *                    Elemental: Largest integer used to index variable (MVAR)
 *      Col. 30 - 42  Compressed Column: Number of columns (NCOL)
 *                    Elemental: Number of element matrices (NELT)
 *      Col. 44 - 56  Compressed Column: Number of entries (NNZERO)
 *                    Elemental: Number of variable indeces (NVARIX)
 *      Col. 58 - 70  Compressed Column: Unused, explicitly zero
 *                    Elemental: Number of elemental matrix entries (NELTVL)
 *
 * Line 4 (2A16, A20)
 *      Col. 1 - 16   Fortran format for pointers (PTRFMT)
 *      Col. 17 - 32  Fortran format for row (or variable) indices (INDFMT)
 *      Col. 33 - 52  Fortran format for numerical values of coefficient matrix
 *                    (VALFMT)
 *                    (blank in the case of matrix patterns)
 *
 * The three character type field on line 3 describes the matrix type.
 * The following table lists the permitted values for each of the three
 * characters. As an example of the type field, RSA denotes that the matrix
 * is real, symmetric, and assembled.
 *
 * First Character:
 *      R Real matrix
 *      C Complex matrix
 *      I integer matrix
 *      P Pattern only (no numerical values supplied)
 *      Q Pattern only (numerical values supplied in associated auxiliary value
 *        file)
 *
 * Second Character:
 *      S Symmetric
 *      U Unsymmetric
 *      H Hermitian
 *      Z Skew symmetric
 *      R Rectangular
 *
 * Third Character:
 *      A Compressed column form
 *      E Elemental form
 *
 * </pre>
 */
#include <stdio.h>
#include "slu_mt_ddefs.h"

/*! \brief Eat up the rest of the current line */
static int_t dDumpLine(FILE *fp)
{
    register int_t c;
    while ((c = fgetc(fp)) != '\n') ;
    return 0;
}

static int_t dParseIntFormat(char *buf, int_t *num, int_t *size)
{
    char *tmp;

    tmp = buf;
    while (*tmp++ != '(') ;
    *num = atoi(tmp);
    while (*tmp != 'I' && *tmp != 'i') ++tmp;
    ++tmp;
    *size = atoi(tmp);
    return 0;
}

static int_t dParseFloatFormat(char *buf, int_t *num, int_t *size)
{
    char *tmp, *period;

    tmp = buf;
    while (*tmp++ != '(') ;
    *num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
    while (*tmp != 'E' && *tmp != 'e' && *tmp != 'D' && *tmp != 'd'
           && *tmp != 'F' && *tmp != 'f') {
        /* May find kP before nE/nD/nF, like (1P6F13.6). In this case the
           num picked up refers to P, which should be skipped. */
        if (*tmp=='p' || *tmp=='P') {
           ++tmp;
           *num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
        } else {
           ++tmp;
        }
    }
    ++tmp;
    period = tmp;
    while (*period != '.' && *period != ')') ++period ;
    *period = '\0';
    *size = atoi(tmp); /*sscanf(tmp, "%2d", size);*/

    return 0;
}

static int_t ReadVector(FILE *fp, int_t n, int_t *where, int_t perline, int_t persize)
{
    register int_t i, j, item;
    char tmp, buf[100];

    i = 0;
    while (i < n) {
        fgets(buf, 100, fp);    /* read a line at a time */
        for (j=0; j<perline && i<n; j++) {
            tmp = buf[(j+1)*persize];     /* save the char at that place */
            buf[(j+1)*persize] = 0;       /* null terminate */
            item = atoi(&buf[j*persize]); 
            buf[(j+1)*persize] = tmp;     /* recover the char at that place */
            where[i++] = item - 1;
        }
    }

    return 0;
}

static int_t dReadValues(FILE *fp, int_t n, double *destination, int_t perline,
        int_t persize)
{
    register int_t i, j, k, s;
    char tmp, buf[100];

    i = 0;
    while (i < n) {
        fgets(buf, 100, fp);    /* read a line at a time */
        for (j=0; j<perline && i<n; j++) {
            tmp = buf[(j+1)*persize];     /* save the char at that place */
            buf[(j+1)*persize] = 0;       /* null terminate */
            s = j*persize;
            for (k = 0; k < persize; ++k) /* No D_ format in C */
                if ( buf[s+k] == 'D' || buf[s+k] == 'd' ) buf[s+k] = 'E';
            destination[i++] = atof(&buf[s]);
            buf[(j+1)*persize] = tmp;     /* recover the char at that place */
        }
    }

    return 0;
}



void
dreadrb(int_t *nrow, int_t *ncol, int_t *nonz,
        double **nzval, int_t **rowind, int_t **colptr)
{

    register int_t i, j, numer_lines = 0;
    int_t tmp, colnum, colsize, rownum, rowsize, valnum, valsize;
    char buf[100], type[4];
    FILE *fp;

    fp = stdin;

    /* Line 1 */
    fgets(buf, 100, fp);
    fputs(buf, stdout);

    /* Line 2 */
    for (i=0; i<4; i++) {
        j = fscanf(fp, "%14c", buf); buf[14] = 0;
        tmp = atoi(buf); /*sscanf(buf, "%d", &tmp);*/
        if (i == 3) numer_lines = tmp;
    }
    dDumpLine(fp);

    /* Line 3 */
    j = fscanf(fp, "%3c", type);
    j = fscanf(fp, "%11c", buf); /* pad */
    type[3] = 0;
#ifdef DEBUG
    printf("Matrix type %s\n", type);
#endif

    fscanf(fp, "%14c", buf); *nrow = atoi(buf); 
    fscanf(fp, "%14c", buf); *ncol = atoi(buf); 
    fscanf(fp, "%14c", buf); *nonz = atoi(buf); 
    fscanf(fp, "%14c", buf); tmp = atoi(buf);   

    if (tmp != 0)
        printf("This is not an assembled matrix!\n");
    if (*nrow != *ncol)
        printf("Matrix is not square.\n");
    dDumpLine(fp);

    /* Allocate storage for the three arrays ( nzval, rowind, colptr ) */
    dallocateA(*ncol, *nonz, nzval, rowind, colptr);

    /* Line 4: format statement */
    j = fscanf(fp, "%16c", buf);
    dParseIntFormat(buf, &colnum, &colsize);
    j = fscanf(fp, "%16c", buf);
    dParseIntFormat(buf, &rownum, &rowsize);
    j = fscanf(fp, "%20c", buf);
    dParseFloatFormat(buf, &valnum, &valsize);
    dDumpLine(fp);

#ifdef DEBUG
    printf("%d rows, %d nonzeros\n", *nrow, *nonz);
    printf("colnum %d, colsize %d\n", colnum, colsize);
    printf("rownum %d, rowsize %d\n", rownum, rowsize);
    printf("valnum %d, valsize %d\n", valnum, valsize);
#endif

    ReadVector(fp, *ncol+1, *colptr, colnum, colsize);
    ReadVector(fp, *nonz, *rowind, rownum, rowsize);
    if ( numer_lines ) {
        dReadValues(fp, *nonz, *nzval, valnum, valsize);
    }

    fclose(fp);
}
