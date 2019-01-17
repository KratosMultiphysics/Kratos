/**
 *
 * @file spm_norm_tests.c
 *
 * Tests and validate the spm_norm routines.
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
    spmatrix_t    spm;
    spm_driver_t driver;
    char *filename;
    spm_mtxtype_t spmtype, mtxtype;
    spm_fmttype_t fmttype;
    int baseval;
    int rc = SPM_SUCCESS;
    int err = 0;

    /**
     * Get options from command line
     */
    spmGetOptions( argc, argv,
                   &driver, &filename );

    rc = spmReadDriver( driver, filename, &spm );
    free(filename);

    if ( rc != SPM_SUCCESS ) {
        fprintf(stderr, "ERROR: Could not read the file, stop the test !!!\n");
        return EXIT_FAILURE;
    }

    if ( spm.flttype == SpmPattern ) {
        spmGenFakeValues( &spm );
    }

    spmtype = spm.mtxtype;
    printf(" -- SPM Norms Test --\n");

    for( fmttype=SpmCSC; fmttype<=SpmIJV; fmttype++ ) {

        spmConvert( fmttype, &spm );

        for( baseval=0; baseval<2; baseval++ )
        {
            spmBase( &spm, baseval );

            for( mtxtype=SpmGeneral; mtxtype<=SpmHermitian; mtxtype++ )
            {
                if ( (mtxtype == SpmHermitian) &&
                     ( ((spm.flttype != SpmComplex64) && (spm.flttype != SpmComplex32)) ||
                       (spmtype != SpmHermitian) ) )
                {
                    continue;
                }
                if ( (mtxtype != SpmGeneral) &&
                     (spmtype == SpmGeneral) )
                {
                    continue;
                }
                spm.mtxtype = mtxtype;

                printf(" Case: %s / %s / %d / %s\n",
                       fltnames[spm.flttype],
                       fmtnames[spm.fmttype], baseval,
                       mtxnames[mtxtype - SpmGeneral] );

                switch( spm.flttype ){
                case SpmComplex64:
                    rc = z_spm_norm_check( &spm );
                    break;

                case SpmComplex32:
                    rc = c_spm_norm_check( &spm );
                    break;

                case SpmFloat:
                    rc = s_spm_norm_check( &spm );
                    break;

                case SpmDouble:
                default:
                    rc = d_spm_norm_check( &spm );
                }
                PRINT_RES(rc);
            }
        }
    }
    spmExit( &spm  );

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
