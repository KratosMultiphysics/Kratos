/**
 *
 * @file p_spm_tests.c
 *
 * @copyright 2011-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Tests and validate the spm_convert routines.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <spm_tests.h>
#include "cblas.h"
#include "lapacke.h"
#include <p_spm.h>

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
void
p_spm_print_check( char *filename, const spmatrix_t *spm )
{
    char *file;
    FILE *f;
    int rc;

    rc = asprintf( &file, "expand_%s_sparse_cp.dat", filename );
    if ( (f = fopen( file, "w" )) == NULL ) {
        perror("p_spm_print_check:sparse_cp");
        return;
    }
    p_spmPrint( f, spm );
    fclose(f);
    free(file);

    if ( spm->dof != 1 ) {
        spmatrix_t *espm = p_spmExpand( spm );

        rc = asprintf( &file, "expand_%s_sparse_ucp.dat", filename );
        if ( (f = fopen( file, "w" )) == NULL ) {
            perror("p_spm_print_check:sparse_ucp");
            return;
        }
        p_spmPrint( f, espm );
        fclose(f);
        free(file);

        spmExit( espm );
        free(espm);
    }

    (void)rc;
    return;
}
