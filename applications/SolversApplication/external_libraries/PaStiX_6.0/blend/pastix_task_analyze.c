/**
 *
 * @file pastix_task_analyze.c
 *
 * PaStiX analyse task function
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_users
 *
 * @brief Perform all the preprocessing steps: ordering, symbolic factorization,
 * reordering, proportionnal mapping, ...
 *
 * This function combines the calls to all the preprocessing steps. See
 * pastix_subtask_order() for further details on the ordering step,
 * pastix_subtask_symbfact() for further details on the symbolic factotrization
 * step, pastix_subtask_reordering() for further details on the reordering process,
 * and pastix_subtask_blend() for the final analyse step that performs
 * proportionnal mapin and generates the numerical factorization data
 * structures.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[in]    spm
 *          The sparse matrix descriptor that describes problem instance.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastix_task_analyze( pastix_data_t *pastix_data,
                     spmatrix_t    *spm )
{
    int rc;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_task_analyze: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_INIT) ) {
        errorPrint("pastix_task_analyze: pastixInit() has to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    /*
     * Ordering step
     */
    rc = pastix_subtask_order( pastix_data, spm, NULL );
    if (rc != PASTIX_SUCCESS) {
        return rc;
    }

    /*
     * Symbolic factorization step
     */
    rc = pastix_subtask_symbfact( pastix_data );
    if (rc != PASTIX_SUCCESS) {
        return rc;
    }

    /*
     * Reordering step
     */
    rc = pastix_subtask_reordering( pastix_data );
    if (rc != PASTIX_SUCCESS) {
        return rc;
    }

    /*
     * Blend step (Proportional mapping)
     */
    rc = pastix_subtask_blend( pastix_data );
    if (rc != PASTIX_SUCCESS) {
        return rc;
    }

    return PASTIX_SUCCESS;
}
