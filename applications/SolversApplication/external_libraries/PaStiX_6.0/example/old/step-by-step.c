/**
 *
 * @file old/step-by-step.c
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * This an example calling PaStiX in step-by-step mode.
 * If runs one full analyze (ordering, symbolic factorization, analyze), then
 * it loops over 2 factorizations that are both used for 2 solves each.
 *
 * @version 6.0.1
 * @author Hastaran Matias
 * @date 2018-07-16
 *
 **/
#include <pastix.h>
#include <pastix/old_api.h>
#include <spm.h>

int main (int argc, char **argv)
{

    pastix_data_t  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix           */
    pastix_float_t *b           = NULL; /* right hand side                                           */
    pastix_int_t    iparm[IPARM_SIZE]; /* integer parameters for pastix                             */
    double          dparm[DPARM_SIZE]; /* floating parameters for pastix                            */
    spm_driver_t    driver;    /* Matrix driver(s) requested by user                        */
    char           *filename;  /* Filename(s) given by user                                 */
    long            i;
    int             j;
    int             nfact       = 2;
    int             nsolv       = 2;
    int             nrhs        = 1;
    spmatrix_t     *spm;
    spmatrix_t      spm2;
    void           *x           = NULL;
    void           *x0          = NULL;
    size_t          size;
    int             check       = 1;
    int             rc          = 0;

    /*
     * Initialize parameters to default values
     */
    iparm[IPARM_MODIFY_PARAMETER] = API_NO;
    pastix( &pastix_data, MPI_COMM_WORLD, -1, NULL, NULL, NULL,
            NULL, NULL, NULL, 1, iparm, dparm );

    /*
     * Update options from command line, and get the matrix filename
     */
    pastixGetOptions( argc, argv,
                      iparm, dparm,
                      &check, &driver, &filename );

    /*
     * Read Matrice
     */
    spm = malloc( sizeof( spmatrix_t ) );
    spmReadDriver( driver, filename, spm );
    free( filename );

    spmPrintInfo( spm, stdout );

    /*
     * Check Matrix format
     */
    rc = spmCheckAndCorrect( spm, &spm2 );
    if ( rc != 0 ) {
        spmExit( spm );
        *spm = spm2;
    }

    /*
     * Generate a Fake values array if needed for the numerical part
     */
    if ( spm->flttype == SpmPattern ) {
        spmGenFakeValues( spm );
    }

    iparm[IPARM_FLOAT]    = spm->flttype;
    iparm[IPARM_MTX_TYPE] = spm->mtxtype;
    iparm[IPARM_DOF_NBR]  = spm->dof;

    /*
     * Step 0 - Initialize pastix
     */
    iparm[IPARM_START_TASK] = API_TASK_INIT;
    iparm[IPARM_END_TASK]   = API_TASK_INIT;
    pastix( &pastix_data, MPI_COMM_WORLD,
            -1, NULL, NULL, NULL,
            NULL, NULL, NULL, 1, iparm, dparm );

    /*
     * Step 1 - Ordering / Scotch
     * Perform it only when the pattern of matrix change.
     * eg: mesh refinement
     * In many cases users can simply go from API_TASK_ORDERING to API_TASK_ANALYSE
     * in one call.
     */
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK]   = API_TASK_ORDERING;
    pastix(&pastix_data, MPI_COMM_WORLD,
           spm->n, spm->colptr, spm->rowptr, spm->values,
           NULL, NULL, NULL, nrhs, iparm, dparm );


    /*
     * Step 2 - Symbolic factorization
     * Perform it only when the pattern of matrix change.
     */
    iparm[IPARM_START_TASK] = API_TASK_SYMBFACT;
    iparm[IPARM_END_TASK]   = API_TASK_SYMBFACT;

    pastix(&pastix_data, MPI_COMM_WORLD,
           spm->n, spm->colptr, spm->rowptr, spm->values,
           NULL, NULL, NULL, nrhs, iparm, dparm );

    /*
     * Step 3 - Mapping and Compute scheduling
     * Perform it only when the pattern of matrix change.
     */
    iparm[IPARM_START_TASK] = API_TASK_ANALYSE;
    iparm[IPARM_END_TASK]   = API_TASK_ANALYSE;
    pastix(&pastix_data, MPI_COMM_WORLD, spm->n, spm->colptr, spm->rowptr, spm->values,
           NULL, NULL, NULL, nrhs, iparm, dparm );

    size = pastix_size_of( spm->flttype ) * spm->n;
    x = malloc( size );
    b = malloc( size );
    if ( check > 1 ) {
        x0 = malloc( size );
    }

    /**
     * Normalize A matrix (optional, but recommended for low-rank functionality)
     */
    double normA = spmNorm( SpmFrobeniusNorm, spm );
    spmScalMatrix( 1./normA, spm );

    /* Do nfact factorization */
    for (i = 0; i < nfact; i++)
    {
        /*
         * Step 4 - Numerical Factorisation
         * Perform it each time the values of the
         * matrix changed.
         */
        fprintf(stdout, "\t> Factorisation number %ld <\n", (long)(i+1));
        iparm[IPARM_START_TASK] = API_TASK_NUMFACT;
        iparm[IPARM_END_TASK]   = API_TASK_NUMFACT;
        pastix(&pastix_data, MPI_COMM_WORLD, spm->n, spm->colptr, spm->rowptr, spm->values,
               NULL, NULL, NULL, nrhs, iparm, dparm );

        /* Do two solve */
        for (j = 0; j < nsolv; j++)
        {

            /*
             * Generates the b and x vector such that A * x = b
             * Compute the norms of the initial vectors if checking purpose.
             */
            if ( check )
            {
                spmGenRHS( SpmRhsRndX, nrhs, spm, x0, spm->n, b, spm->n );
                memcpy( x, b, size );
            }
            else {
                spmGenRHS( SpmRhsRndB, nrhs, spm, NULL, spm->n, x, spm->n );

                /* Apply also normalization to b vector */
                spmScalVector( spm->flttype, 1./normA, spm->n * nrhs, b, 1 );

                /* Save b for refinement: TODO: make 2 examples w/ or w/o refinement */
                memcpy( b, x, size );
            }

            /*
             * Step 5.1 - Solve
             * If you don't need iterative refinement
             * x contains the RHS b as input
             * x returns the solution as output
             */
            iparm[IPARM_START_TASK] = API_TASK_SOLVE;
            iparm[IPARM_END_TASK]   = API_TASK_SOLVE;

            fprintf(stdout, "\t>> Solve step number %ld  <<\n", (long)(j+1));
            pastix(&pastix_data, MPI_COMM_WORLD,
                   spm->n, spm->colptr, spm->rowptr, spm->values,
                   NULL, NULL, x, nrhs, iparm, dparm );

            /*
             * Step 5.2 - Refinnement
             * b contains the RHS b as input
             * b returns the soluton as output
             */
            iparm[IPARM_START_TASK] = API_TASK_REFINE;
            iparm[IPARM_END_TASK]   = API_TASK_REFINE;

            fprintf(stdout, "\t>> Refine step number %ld  <<\n", (long)(j+1));
            pastix(&pastix_data, MPI_COMM_WORLD,
                   spm->n, spm->colptr, spm->rowptr, spm->values,
                   NULL, NULL, b, nrhs, iparm, dparm );
            if (check) {
                rc = spmCheckAxb( dparm[DPARM_EPSILON_REFINEMENT], nrhs, spm, x0, spm->n, b, spm->n, x, spm->n );
            }
        }
    }

    /*
     * Step 6 - Clean structures
     * When you don't need PaStiX anymore
     */
    iparm[IPARM_START_TASK] = API_TASK_CLEAN;
    iparm[IPARM_END_TASK]   = API_TASK_CLEAN;

    pastix(&pastix_data, MPI_COMM_WORLD,
           spm->n, spm->colptr, spm->rowptr, spm->values,
           NULL, NULL, x, nrhs, iparm, dparm );

    spmExit( spm );
    free( spm );
    free( b );
    free( x );
    if ( x0 ) {
        free( x0 );
    }
    return rc;
}
