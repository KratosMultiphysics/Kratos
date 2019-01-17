/**
 *
 * @file bvec_tests.c
 *
 * Tests for preformance of parallel subroutines.
 * We test each subroutines of vectorial functions (copy, axpy, scale, ...)
 * and we compute the average time over 50 calls to each functions.
 * Size of the matrices/vectors are given by the first two dimensions of the
 * Laplacian driver.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Vincent Bridonneau
 * @date 2018-07-16
 *

 **/
#include <pastix.h>
#include <stdlib.h>
#include <string.h>
#include "refinement/z_refine_functions.h"
#include "bcsc/bcsc.h"
#include "bcsc/bcsc_z.h"
#include <time.h>
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */

char* fltnames[] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
char* schednames[] = { "Sequential", "Static" };

int
z_bvec_check( pastix_data_t *pastix_data,
                        int m );
int
c_bvec_check( pastix_data_t *pastix_data,
                        int m );
int
d_bvec_check( pastix_data_t *pastix_data,
                        int m );
int
s_bvec_check( pastix_data_t *pastix_data,
                        int m );

int main ( int argc, char **argv )
{
    pastix_data_t      *pastix_data = NULL;
    pastix_int_t        iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    pastix_fixdbl_t     dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    spm_driver_t        driver = SpmDriverRSA;
    char               *filename = NULL;
    int                 check = 1;
    pastix_int_t        m, n, s, t;
    int rc = 0;
    pastix_int_t    dim3;
    spm_coeftype_t  flttype;
    pastix_fixdbl_t alpha, beta;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );
    iparm[IPARM_VERBOSE] = 0;

    /**
     * Get options from command line
     * Prevent from failing if no arguments is given
     */
    if ( argc > 1 ) {
        pastixGetOptions( argc, argv,
                          iparm, dparm,
                          &check, &driver, &filename );
    }

    if ( driver != SpmDriverLaplacian ) {
        flttype = SpmDouble;
        m = 10000;
        n = iparm[IPARM_GMRES_IM];
    }
    else {
        spmParseLaplacianInfo( filename, &flttype, &m, &n, &dim3, &alpha, &beta );
    }

    /**
     * Startup pastix to start the scheduler
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    if ( check ) {
        for( s=PastixSchedSequential; s<=PastixSchedStatic; s++ )
        {
            pastix_data->iparm[IPARM_SCHEDULER] = s;
            for( t=SpmFloat; t<=SpmComplex64; t++ )
            {
                printf("  Case %s - %s:\n ",
                       fltnames[t], schednames[s] );
                switch( flttype ){
                case SpmComplex64:
                    rc = z_bvec_check( pastix_data, m );
                    break;

                case SpmComplex32:
                    rc = c_bvec_check( pastix_data, m );
                    break;

                case SpmFloat:
                    rc = s_bvec_check( pastix_data, m );
                    break;

                case SpmDouble:
                default:
                    rc = d_bvec_check( pastix_data, m );
                }
            }
        }
    }
    else {
        switch( flttype ){
        case SpmComplex64:
            rc = z_bvec_check( pastix_data, m );
            break;

        case SpmComplex32:
            rc = c_bvec_check( pastix_data, m );
            break;

        case SpmFloat:
            rc = s_bvec_check( pastix_data, m );
            break;

        case SpmDouble:
        default:
            rc = d_bvec_check( pastix_data, m );
        }
    }

    pastixFinalize( &pastix_data );

    return rc;
}
