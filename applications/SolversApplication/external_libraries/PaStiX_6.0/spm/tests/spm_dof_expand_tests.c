/**
 *
 * @file spm_dof_expand_tests.c
 *
 * Tests and validate the spmNorm routines when the spm_tests.hold constant and/or variadic dofs.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
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
#include <time.h>
#include "spm_tests.h"

int z_spm_norm_check( const spmatrix_t *spm );
int c_spm_norm_check( const spmatrix_t *spm );
int d_spm_norm_check( const spmatrix_t *spm );
int s_spm_norm_check( const spmatrix_t *spm );

void z_spm_print_check( char *filename, const spmatrix_t *spm );
void c_spm_print_check( char *filename, const spmatrix_t *spm );
void d_spm_print_check( char *filename, const spmatrix_t *spm );
void s_spm_print_check( char *filename, const spmatrix_t *spm );
void p_spm_print_check( char *filename, const spmatrix_t *spm );

#define PRINT_RES(_ret_)                        \
    if(_ret_) {                                 \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

char* fltnames[] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
char* fmtnames[] = { "CSC", "CSR", "IJV" };
char* mtxnames[] = { "General", "Symmetric", "Hermitian" };

int main (int argc, char **argv)
{
    spmatrix_t    original, *spm;
    spm_driver_t driver;
    char *filename;
    spm_mtxtype_t spmtype;
    spm_mtxtype_t mtxtype;
    spm_fmttype_t fmttype;
    int baseval, i, rc, dofmax = 3;

    /**
     * Get options from command line
     */
    spmGetOptions( argc, argv,
                   &driver, &filename );

    rc = spmReadDriver( driver, filename, &original );
    free(filename);

    if ( rc != SPM_SUCCESS ) {
        fprintf(stderr, "ERROR: Could not read the file, stop the test !!!\n");
        return EXIT_FAILURE;
    }

    spmtype = original.mtxtype;
    printf(" -- SPM Dof Expand Test --\n");

    for( i=0; i<2; i++ )
    {
        for( mtxtype=SpmGeneral; mtxtype<=SpmHermitian; mtxtype++ )
        {
            if ( (mtxtype == SpmHermitian) &&
                 ( ((original.flttype != SpmComplex64) && (original.flttype != SpmComplex32)) ||
                   (spmtype != SpmHermitian) ) )
            {
                continue;
            }
            if ( (mtxtype != SpmGeneral) &&
                 (spmtype == SpmGeneral) )
            {
                continue;
            }
            original.mtxtype = mtxtype;

            for( baseval=0; baseval<2; baseval++ )
            {
                spmBase( &original, baseval );

                for( fmttype=SpmCSC; fmttype<=SpmIJV; fmttype++ )
                {
                    spmConvert( fmttype, &original );
                    spm = spmDofExtend( &original, i, dofmax );

                    rc = asprintf( &filename, "%d_%s_%d_%s_%s",
                                   i, fmtnames[fmttype], baseval,
                                   mtxnames[mtxtype - SpmGeneral],
                                   fltnames[spm->flttype] );

                    printf( "-- %s --\n", filename );
                    switch( spm->flttype ){
                    case SpmComplex64:
                        z_spm_print_check( filename, spm );
                        break;

                    case SpmComplex32:
                        c_spm_print_check( filename, spm );
                        break;

                    case SpmFloat:
                        s_spm_print_check( filename, spm );
                        break;

                    case SpmPattern:
                        p_spm_print_check( filename, spm );
                        break;

                    case SpmDouble:
                    default:
                        d_spm_print_check( filename, spm );
                    }
                    free(filename);

                    spmExit( spm );
                    free(spm);
                    spm = NULL;
                }
            }
        }
    }
    spmExit( &original );

    (void)rc;
    return EXIT_SUCCESS;
}
