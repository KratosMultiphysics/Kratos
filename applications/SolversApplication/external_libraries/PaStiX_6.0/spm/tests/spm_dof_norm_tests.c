/**
 *
 * @file spm_dof_norm_tests.c
 *
 * Tests and validate the spm_norm routines when the spm_tests.hold constant and/or variadic dofs.
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
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <spm_tests.h>

int z_spm_norm_check( const spmatrix_t *spm );
int c_spm_norm_check( const spmatrix_t *spm );
int d_spm_norm_check( const spmatrix_t *spm );
int s_spm_norm_check( const spmatrix_t *spm );

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
    spm_driver_t  driver;
    char         *filename;
    spm_mtxtype_t spmtype, mtxtype;
    spm_fmttype_t fmttype;
    int baseval;
    int rc = SPM_SUCCESS;
    int err = 0;
    int i, dofmax = 4;

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

    if ( original.flttype == SpmPattern ) {
        spmGenFakeValues( &original );
    }

    spmtype = original.mtxtype;
    printf(" -- SPM Norms Dof Test --\n");

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

                    printf(" Case: %d / %s / %s / %d / %s\n",
                           i, fltnames[spm->flttype],
                           fmtnames[spm->fmttype], baseval,
                           mtxnames[mtxtype - SpmGeneral] );

                    switch( spm->flttype ){
                    case SpmComplex64:
                        rc = z_spm_norm_check( spm );
                        break;

                    case SpmComplex32:
                        rc = c_spm_norm_check( spm );
                        break;

                    case SpmFloat:
                        rc = s_spm_norm_check( spm );
                        break;

                    case SpmDouble:
                    default:
                        rc = d_spm_norm_check( spm );
                    }
                    PRINT_RES(rc);

                    spmExit( spm );
                    free(spm);
                    spm = NULL;
                }
            }
        }
    }
    spmExit( &original );

    if( err == 0 ) {
        printf(" -- All tests PASSED --\n");
        return EXIT_SUCCESS;
    }
    else
    {
        printf(" -- %d tests FAILED --\n", err);
        return EXIT_FAILURE;
    }
}
