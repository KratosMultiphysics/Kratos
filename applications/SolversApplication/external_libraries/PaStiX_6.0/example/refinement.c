/**
 * @file refinement.c
 *
 * @brief A refinement example that runs iterative methods without applying
 *        the preconditioner (factorization and solve steps are not called).
 *        Based on the step-by-step example.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 * @ingroup pastix_examples
 * @code
 *
 */
#include <pastix.h>
#include <spm.h>

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t    iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    double          dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    spm_driver_t driver;
    char           *filename;
    pastix_spm_t   *spm, spm2;
    void           *x, *b, *x0 = NULL;
    size_t          size;
    int             check = 1;
    int             nrhs  = 1;
    int             rc    = 0;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );

    /**
     * Get options from command line
     */
    pastixGetOptions( argc, argv,
                      iparm, dparm,
                      &check, &driver, &filename );

    /**
     * Read the sparse matrix with the driver
     */
    spm = malloc( sizeof( pastix_spm_t ) );
    spmReadDriver( driver, filename, spm );
    free( filename );

    spmPrintInfo( spm, stdout );

    rc = spmCheckAndCorrect( spm, &spm2 );
    if ( rc != 0 ) {
        spmExit( spm );
        *spm = spm2;
    }

    /**
     * Generate a Fake values array if needed for the numerical part
     */
    if ( spm->flttype == PastixPattern ) {
        spmGenFakeValues( spm );
    }

    /**
     * Startup PaStiX
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    pastix_subtask_order( pastix_data, spm, NULL );
    pastix_subtask_symbfact( pastix_data );
    pastix_subtask_reordering( pastix_data );
    pastix_subtask_blend( pastix_data );

    /**
     * Normalize A matrix (optional, but recommended for low-rank functionality)
     */
    double normA = spmNorm( SpmFrobeniusNorm, spm );
    spmScalMatrix( 1./normA, spm );

    size = pastix_size_of( spm->flttype ) * spm->n * nrhs;
    x = malloc( size );
    b = malloc( size );
    if ( check > 1 ) {
        x0 = malloc( size );
    }

    pastix_subtask_spm2bcsc( pastix_data, spm );
    //pastix_subtask_bcsc2ctab( pastix_data );
    //pastix_subtask_sopalin( pastix_data );

    /**
     * Generates the b and x vector such that A * x = b
     * Compute the norms of the initial vectors if checking purpose.
     */
    if ( check )
      {
        spmGenRHS( SpmRhsRndX, nrhs, spm, x0, spm->n, b, spm->n );
        memcpy( x, b, size );
      }
    else
      {
        spmGenRHS( SpmRhsRndB, nrhs, spm, NULL, spm->n, x, spm->n );

        /* Apply also normalization to b vectors */
        spmScalVector( spm->flttype, 1./normA, spm->n * nrhs, b, 1 );

        /* Save b for refinement */
        memcpy( b, x, size );
      }

    /**
     * Solve the linear system
     */
    //pastix_task_solve( pastix_data, nrhs, x, spm->n );
    pastix_task_refine( pastix_data, spm->n, nrhs, b, spm->n, x, spm->n );

    if ( check ) {
      rc |= spmCheckAxb( dparm[DPARM_EPSILON_REFINEMENT], nrhs, spm, x0, spm->n, b, spm->n, x, spm->n );
    }
    if (iparm[IPARM_NBITER]>=iparm[IPARM_ITERMAX]) {
      rc |= 1;
    }

    spmExit( spm );
    free( spm );
    free( b );
    free( x );
    if ( x0 ) {
        free( x0 );
    }
    pastixFinalize( &pastix_data );

    return rc;
}

/**
 * @endcode
 */
