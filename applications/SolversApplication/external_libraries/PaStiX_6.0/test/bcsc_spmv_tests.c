/**
 *
 * @file bcsc_spmv_tests.c
 *
 * @copyright 2011-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Tests and validate the bcsc_spmv routines.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2018-07-16
 *
 **/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <pastix.h>
#include "common.h"
#include <spm.h>
#include <bcsc.h>
#include "sopalin_data.h"

int z_bcsc_spmv_check( spm_trans_t trans, const spmatrix_t *spm, const pastix_data_t *pastix_data );
int c_bcsc_spmv_check( spm_trans_t trans, const spmatrix_t *spm, const pastix_data_t *pastix_data );
int d_bcsc_spmv_check( spm_trans_t trans, const spmatrix_t *spm, const pastix_data_t *pastix_data );
int s_bcsc_spmv_check( spm_trans_t trans, const spmatrix_t *spm, const pastix_data_t *pastix_data );

#define PRINT_RES(_ret_)                        \
    if(_ret_) {                                 \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

char* fltnames[] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
char* transnames[] = { "NoTrans", "Trans", "ConjTrans" };
char* mtxnames[] = { "General", "Symmetric", "Hermitian" };

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix  */
    pastix_int_t    iparm[IPARM_SIZE];  /* integer parameters for pastix                    */
    double          dparm[DPARM_SIZE];  /* floating parameters for pastix                   */
    spm_driver_t    driver;             /* Matrix driver(s) requested by user               */
    spmatrix_t     *spm, spm2;
    char *filename;                     /* Filename(s) given by user                        */
    int t;
    int ret = PASTIX_SUCCESS;
    int err = 0;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Get options from command line
     */
    pastixGetOptions( argc, argv,
                      NULL, NULL,
                      NULL, &driver, &filename );

    spm = malloc( sizeof( spmatrix_t ) );
    spmReadDriver( driver, filename, spm );
    free(filename);

    ret = spmCheckAndCorrect( spm, &spm2 );
    if ( ret != 0 ) {
        spmExit( spm );
        *spm = spm2;
    }

    spmBase( spm, 0 );

    if ( spm->flttype == SpmPattern ) {
        spmGenFakeValues( spm );
    }

    /**
     * Run preprocessing steps required to generate the blocked csc
     */
    pastix_task_analyze( pastix_data, spm );

    /**
     * Generate the blocked csc
     */
    pastix_data->bcsc = malloc( sizeof(pastix_bcsc_t) );
    bcscInit( spm,
              pastix_data->ordemesh,
              pastix_data->solvmatr,
              spm->mtxtype == SpmGeneral,
              pastix_data->bcsc );

    printf(" -- BCSC MatVec Test --\n");
    for( t=PastixNoTrans; t<=PastixConjTrans; t++ )
    {
        if ( (t == PastixConjTrans) &&
             ((spm->flttype != SpmComplex64) && (spm->flttype != SpmComplex32)) )
        {
            continue;
        }
        if ( (spm->mtxtype != SpmGeneral) && (t != PastixNoTrans) )
        {
            continue;
        }
        printf("   Case %s - %s - %s:\n",
               fltnames[spm->flttype],
               mtxnames[spm->mtxtype - SpmGeneral],
               transnames[t - PastixNoTrans] );

        switch( spm->flttype ){
        case SpmComplex64:
            ret = z_bcsc_spmv_check( t, spm, pastix_data );
            break;

        case SpmComplex32:
            ret = c_bcsc_spmv_check( t, spm, pastix_data );
            break;

        case SpmFloat:
            ret = s_bcsc_spmv_check( t, spm, pastix_data );
            break;

        case SpmDouble:
        default:
            ret = d_bcsc_spmv_check( t, spm, pastix_data );
        }
        PRINT_RES(ret);
    }

    spmExit( spm );
    free( spm );

    pastixFinalize( &pastix_data );

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
