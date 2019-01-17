/**
 *
 * @file pastix.c
 *
 * PaStiX main interface for compatibility with former releases
 *
 * @copyright 2011-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Mathias Hastaran
 * @date 2018-07-16
 *
 **/
#include "common.h"
#include "spm.h"
#include "pastix/order.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 * @brief Main function for compatibility with former releases
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data_ptr
 *          The pastix data structure of the solver to store the state of the
 *          solver at every call.
 *
 * @param[in] pastix_comm
 *          The MPI communicator to use for the distributed solver.
 *
 * @param[in] n
 *          The size of the sparse matrix system to solve.
 *
 * @param[inout] colptr
 *          The pointer to the column index of the sparse matrix in the CSC
 *          format.
 *          On exit, the base value of the array might have changed, and/or the
 *          pointer might have been freed if IPARM_FREE_CSCUSER is set, and the
 *          factorization step is performed.
 *
 * @param[inout] row
 *          The pointer to the row array of the sparse matrix in the CSC
 *          format.
 *          On exit, the base value of the array might have changed, and/or the
 *          pointer might have been freed if IPARM_FREE_CSCUSER is set, and the
 *          factorization step is performed.
 *
 * @param[inout] avals
 *          The pointer to the values array of the sparse matrix in the CSC
 *          format.
 *          On exit, the pointer might have been freed if IPARM_FREE_CSCUSER is
 *          set, and the factorization step is performed.
 *
 * @param[inout] perm
 *          The pointer to the permutation array.
 *          On entry: the pointer might be allocated to store the generated
 *          permutation on exit, or to provide the user permutation.
 *          On exit, the permutation used by the solver is returned if perm is
 *          not NULL.
 *
 * @param[inout] invp
 *          The pointer to the inverse permutation array.
 *          On entry: the pointer might be allocated to store the generated
 *          inverse permutation on exit, or to provide the user permutation.
 *          On exit, the inverse permutation used by the solver is returned if
 *          invp is not NULL.
 *
 * @param[inout] b
 *          Array of size n -by- nrhs
 *          On entry, contains the nrhs vectors of the problem.
 *          On exit, contains the nrhs solution vectors of the problem.
 *
 * @param[in] nrhs
 *          The number of right hand side in the problem.
 *
 * @param[inout] iparm
 *          Array of size IPARM_SIZE
 *          On entry, contains all the integer parameters of the solver.
 *          On exit, the aray is updated with integer outputs of the solver.
 *
 * @param[inout] dparm
 *          Array of size DPARM_SIZE
 *          On entry, contains all the double parameters of the solver.
 *          On exit, the aray is updated with double outputs of the solver.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on succesful exit,
 * @retval PASTIX_ERR_BADPARAMETER on incorrect input parameter,
 * @retval PASTIX_ERR_NOTIMPLEMENTED on variadic dofs,
 * @retval PASTIX_ERR_UNKNOWN on undefined behaviors.
 *
 *******************************************************************************
 */
int
pastix( pastix_data_t **pastix_data_ptr,
        MPI_Comm        pastix_comm,
        pastix_int_t    n,
        pastix_int_t   *colptr,
        pastix_int_t   *row,
        void           *avals,
        pastix_int_t   *perm,
        pastix_int_t   *invp,
        void           *b,
        pastix_int_t    nrhs,
        pastix_int_t   *iparm,
        double         *dparm )
{
    pastix_data_t *pastix_data;
    spmatrix_t   *spm = NULL;
    int ret;
    size_t size;

    /*
     * Initialize iparm/dparm to default values and exit when called with
     * IPARM_MODIFY_PARAMETER set to 0
     */
    if (!iparm[IPARM_MODIFY_PARAMETER])
    {
        pastixInitParam(iparm, dparm);
        iparm[IPARM_MODIFY_PARAMETER] = 1;
        return PASTIX_SUCCESS;
    }

    /*
     * Initialization step
     * Create the pastix_data structure and initialize the runtimes
     */
    if (iparm[IPARM_END_TASK] < PastixTaskInit) {
        return PASTIX_SUCCESS;
    }

    if (iparm[IPARM_START_TASK] == PastixTaskInit) {
        if (*pastix_data_ptr != NULL)
        {
            /*
             * Let's consider the user want to restart pastix with different
             * parameters
             */
            if (iparm[IPARM_VERBOSE] > PastixVerboseNo) {
                pastix_print( 0, 0, "WARNING: PaStiX schedulers restarted\n" );
            }
            pastixFinalize( pastix_data_ptr );
        }

        /*
         * Initialize pastix_data structure, and start scheduler(s)
         */
        pastixInit( pastix_data_ptr, pastix_comm, iparm, dparm );

        iparm[IPARM_START_TASK]++;
    }

    /*
     * Return now if only initialization is required
     */
    if (iparm[IPARM_END_TASK] < PastixTaskOrdering) {
        return PASTIX_SUCCESS;
    }

    pastix_data = *pastix_data_ptr;

    /*
     * Initialize the internal spm structure.
     * We perform if only if starting step is lower than numerical
     * factorization, because further steps are using the internal bcsc for
     * computations with A.
     */
    if (iparm[IPARM_START_TASK] <= PastixTaskNumfact) {
        if ( (pastix_data->csc != NULL) &&
             ((pastix_data->csc->n      != n)                       ||
              (pastix_data->csc->nnz    != (colptr[n] - colptr[0])) ||
              (pastix_data->csc->colptr != colptr)                  ||
              (pastix_data->csc->rowptr != row)) )
        {
            /*
             * This is a new csc, we need to delete the old one stored in pastix_data, and create a new one
             * We do not use spmExit, because the user allocated the fields.
             */
            memFree_null( pastix_data->csc );
        }

        if ( pastix_data->csc == NULL )
        {
            /*
             * Check and set the matrix type
             */
            if (iparm[IPARM_FLOAT] == -1) {
                printf("Pastix old interface: you have to set iparm[IPARM_FLOAT]\n");
                return PASTIX_ERR_BADPARAMETER;
            }
            if (iparm[IPARM_MTX_TYPE] == -1) {
                printf("Pastix old interface: you have to set iparm[IPARM_MTX_TYPE]\n");
                return PASTIX_ERR_BADPARAMETER;
            }
            if (iparm[IPARM_DOF_NBR] < 1) {
                fprintf(stderr,
                        "pastix: Variadic dofs are not supported in old pastix interface.\n"
                        "        Please switch to the new interface to use this feature, \n"
                        "        or set to a positive value\n");
                return PASTIX_ERR_BADPARAMETER;
            }

            spm = malloc(sizeof( spmatrix_t ));
            spmInit( spm );

            spm->mtxtype = iparm[IPARM_MTX_TYPE];
            spm->flttype = iparm[IPARM_FLOAT];
            spm->fmttype = SpmCSC;

            spm->n    = n;
            spm->nnz  = colptr[n] - colptr[0];
            spm->dof  = iparm[IPARM_DOF_NBR];

            spm->colptr = colptr;
            spm->rowptr = row;
            spm->values = avals;

            spmUpdateComputedFields( spm );

            pastix_data->csc = spm;
        }
        else {
            /* Cast to overwrite the const */
            spm = (spmatrix_t*)(pastix_data->csc);
        }

        /*
         * Update value field if given only at numerical steps
         */
        if ( spm->values == NULL ) {
            spm->values = avals;
        }
    }
    else {
        /* Cast to overwrite the const */
        spm = (spmatrix_t*)(pastix_data->csc);
    }

    /*
     * Ordering
     */
    if (iparm[IPARM_START_TASK] == PastixTaskOrdering)
    {
        pastix_order_t *o = NULL;

        if ( (perm != NULL) || (invp != NULL) ) {
            o = malloc( sizeof(pastix_order_t) );
            ret = pastixOrderAlloc( o, 0, 0 );
            if (PASTIX_SUCCESS != ret) {
                return ret;
            }

            if ( perm != NULL ) {
                MALLOC_INTERN(o->permtab, n, pastix_int_t);
                memcpy( o->permtab, perm, n * sizeof(pastix_int_t) );
                o->vertnbr = n;
            }
            if ( invp != NULL ) {
                MALLOC_INTERN(o->peritab, n, pastix_int_t);
                memcpy( o->peritab, invp, n * sizeof(pastix_int_t) );
                o->vertnbr = n;
            }
        }

        ret = pastix_subtask_order( pastix_data, spm, o );
        if (PASTIX_SUCCESS != ret) {
            return ret;
        }

        if ( o != NULL ) {
            if ( perm != NULL ) {
                assert( o->permtab != NULL );
                assert( o->vertnbr == n );
                memcpy( perm, o->permtab, n * sizeof(pastix_int_t));
            }
            if ( invp != NULL ) {
                assert( o->peritab != NULL );
                assert( o->vertnbr == n );
                memcpy( invp, o->peritab, n * sizeof(pastix_int_t));
            }
            pastixOrderExit(o);
            free(o);
        }
        iparm[IPARM_START_TASK]++;
    }

    /*
     * Symbolic factorization
     */
    if (iparm[IPARM_END_TASK] < PastixTaskSymbfact) {
        return PASTIX_SUCCESS;
    }

    if (iparm[IPARM_START_TASK] == PastixTaskSymbfact)
    {
        ret = pastix_subtask_symbfact( pastix_data );
        if (PASTIX_SUCCESS != ret)
        {
            return ret;
        }
        iparm[IPARM_START_TASK]++;
    }

    /*
     * Analyze step
     */
    if (iparm[IPARM_END_TASK] < PastixTaskAnalyze) {
        return PASTIX_SUCCESS;
    }

    if (iparm[IPARM_START_TASK] == PastixTaskAnalyze)
    {
        ret = pastix_subtask_blend( pastix_data );
        if (PASTIX_SUCCESS != ret)
        {
            return ret;
        }
        iparm[IPARM_START_TASK]++;
    }

    /*
     * Numerical factorisation
     */
    if (iparm[IPARM_END_TASK] < PastixTaskNumfact) {
        return PASTIX_SUCCESS;
    }

    if (iparm[IPARM_START_TASK] == PastixTaskNumfact)
    {
        ret = pastix_task_numfact( pastix_data, spm );
        if (PASTIX_SUCCESS != ret) {
            return ret;
        }
        iparm[IPARM_START_TASK]++;
    }

    /*
     * Solve
     */
    if (iparm[IPARM_END_TASK] < PastixTaskSolve) {
        return PASTIX_SUCCESS;
    }

    if (iparm[IPARM_START_TASK] == PastixTaskSolve) {
        size = pastix_size_of( spm->flttype ) * spm->n;
        if ( pastix_data->x0 ) {
            free(pastix_data->x0);
            pastix_data->x0 = NULL;
        }
        if ( pastix_data->b ) {
            free(pastix_data->b);
            pastix_data->b = NULL;
        }

        /*
         * Backup the initial b if we need to perform an iterative
         * refinement after the solve step
         */
        if (iparm[IPARM_END_TASK] > PastixTaskSolve) {
            pastix_data->b = malloc(size);
            memcpy(pastix_data->b, b, size);
        }
        pastix_task_solve( pastix_data, nrhs, b, spm->n );
        iparm[IPARM_START_TASK]++;

        /*
         * Backup the first solution x0 if the user wants to come back later for
         * iterative refinement
         */
        if (iparm[IPARM_END_TASK] == PastixTaskSolve) {
            pastix_data->x0 = malloc(size);
            memcpy(pastix_data->x0, b, size);
        }
    }

    /*
     * Refinement
     */
    if (iparm[IPARM_END_TASK] < PastixTaskRefine) {
        return PASTIX_SUCCESS;
    }

    if (iparm[IPARM_START_TASK] == PastixTaskRefine) {
        void *refineB  = pastix_data->b;
        void *refineX0 = pastix_data->x0;
        size = pastix_size_of( spm->flttype ) * spm->n;
        if ( !refineB ) {
            if ( !refineX0 ) {
                /*
                 * Neither b or x0 have been saved.
                 * Then, we need to start with x0 as a null vector. For that, we
                 * backup the original b, and we use the given b as x in the
                 * refinement step to store the solution in place.
                 */
                /* refineB = malloc(size); */
                /* memcpy(refineB, b, size); */

                /* refineX0 = b; */
                /* memset(refineX0, 0, size); */
                /* exit(0);  */
                /*
                 * Neither b and x0 have been saved, this should never happen.
                 */
                fprintf(stderr, "Neither b and x0 have been saved, this should never happen\n");
                return PASTIX_ERR_UNKNOWN;
            }
            else {
                /*
                 * x0 is saved, but not b. It means that we exit the pastix
                 * function call between the solve and refinemnet
                 * step. Therefor, b holds the original b.
                 */
                refineB = b;
            }
        }
        else {
            if ( !refineX0 ) {
                /*
                 * b is saved, but not x0. It means that we did not exit the
                 * pastix function call between solve and refinement steps.
                 * Therefor, b holds the initial solution x0 from the solve step.
                 */
                refineX0 = b;
            }
            else {
                /*
                 * Both x0 and b are saved. This should never happen.
                 */
                fprintf(stderr, "Both b and x0 are defined, this should never happen\n");
                return PASTIX_ERR_UNKNOWN;
            }
        }
        pastix_task_refine( pastix_data, spm->n, nrhs,
                            refineB, spm->n, refineX0, spm->n );
        iparm[IPARM_START_TASK]++;

        /*
         * Let's return the solution to the user
         */
        if ( b != refineX0 ) {
            memcpy(b, refineB, size);
        }
    }

    if ( pastix_data->x0 != NULL ) {
        free( pastix_data->x0 );
        pastix_data->x0 = NULL;
    }
    if ( pastix_data->b != NULL ) {
        free( pastix_data->b );
        pastix_data->b = NULL;
    }

    /*
     * Cleaning
     */
    if (iparm[IPARM_END_TASK] < PastixTaskClean) {
        return PASTIX_SUCCESS;
    }

    if (iparm[IPARM_START_TASK] == PastixTaskClean) {
        if ( pastix_data->csc != NULL ) {
            spmatrix_t *spm = (spmatrix_t*)(pastix_data->csc);
            free( spm );
        }
        pastixFinalize( pastix_data_ptr );
        iparm[IPARM_START_TASK]++;
    }

    return PASTIX_SUCCESS;
}
