/**
 *
 * @file bcsc_norm_tests.c
 *
 * @copyright 2011-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Tests and validate the bcsc_norm routines.
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

int z_bcsc_norm_check( const spmatrix_t   *spm, const pastix_bcsc_t *bcsc );
int c_bcsc_norm_check( const spmatrix_t   *spm, const pastix_bcsc_t *bcsc );
int d_bcsc_norm_check( const spmatrix_t   *spm, const pastix_bcsc_t *bcsc );
int s_bcsc_norm_check( const spmatrix_t   *spm, const pastix_bcsc_t *bcsc );

#define PRINT_RES(_ret_)                        \
    if(_ret_) {                                 \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

char* fltnames[] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
char* mtxnames[] = { "General", "Symmetric", "Hermitian" };

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix  */
    pastix_int_t    iparm[IPARM_SIZE];  /* integer parameters for pastix                    */
    double          dparm[DPARM_SIZE];  /* floating parameters for pastix                   */
    spm_driver_t    driver;             /* Matrix driver(s) requested by user               */
    spmatrix_t     *spm, spm2;
    pastix_bcsc_t   bcsc;
    char *filename;                     /* Filename(s) given by user                        */
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
    bcscInit( spm,
              pastix_data->ordemesh,
              pastix_data->solvmatr,
              spm->mtxtype == SpmGeneral, &bcsc );

    printf(" -- BCSC Norms Test --\n");
    printf(" Datatype: %s\n", fltnames[spm->flttype] );
    spmBase( spm, 0 );

    printf("   Matrix type : %s\n", mtxnames[spm->mtxtype - SpmGeneral] );

    switch( spm->flttype ){
    case SpmComplex64:
        ret = z_bcsc_norm_check( spm, &bcsc );
        break;

    case SpmComplex32:
        ret = c_bcsc_norm_check( spm, &bcsc );
        break;

    case SpmFloat:
        ret = s_bcsc_norm_check( spm, &bcsc );
        break;

    case SpmDouble:
    default:
        ret = d_bcsc_norm_check( spm, &bcsc );
    }
    PRINT_RES(ret);

    spmExit( spm );
    free( spm );
    bcscExit( &bcsc );
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
