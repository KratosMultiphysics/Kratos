/**
 *
 * @file sopalin/schur.c
 *
 * PaStiX schur interface functions
 *
 * @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Mathias Hastaran
 * @date 2018-07-16
 *
 * @addtogroup pastix_schur
 * @{
 *
 **/
#include "common.h"
#include "spm.h"
#include "blend/solver.h"
#include "sopalin/coeftab_z.h"
#include "sopalin/coeftab_c.h"
#include "sopalin/coeftab_d.h"
#include "sopalin/coeftab_s.h"

/**
 *******************************************************************************
 *
 * @brief Set the list of unknowns that belongs to the schur complement.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix data structure of the solver to store the list of Schur
 *          unknowns.
 *
 * @param[in] n
 *          The number of unknowns in the Schur complement.
 *
 * @param[in] list
 *          Array of integer of size n.
 *          The list of unknowns belonging to the Schur complement with the same
 *          baseval as the associated spm.
 *
 *******************************************************************************/
void
pastixSetSchurUnknownList( pastix_data_t      *pastix_data,
                           pastix_int_t        n,
                           const pastix_int_t *list)
{
    if ( n > 0 ) {
        pastix_data->schur_n    = n;
        pastix_data->schur_list = (pastix_int_t*)malloc(n * sizeof(pastix_int_t));
        memcpy( pastix_data->schur_list, list, n * sizeof(pastix_int_t) );
    }
}

/**
 *******************************************************************************
 *
 * @brief Return the Schur complement.
 *
 * The Schur complement is returned in the column major layout used by the
 * classic linear algebra libraries such as Blas or Lapack.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix data structure of the problem solved.
 *
 * @param[inout] S
 *          Array of size spm->n -by- lds of arithmetic spm->flttype, where spm
 *          is the spm of the original problem.
 *          On exit, the array contains the Schur complement of the factorized
 *          matrix.
 *
 * @param[in] lds
 *          The leading dimension of the S array.
 *
 ********************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastixGetSchur( const pastix_data_t *pastix_data,
                void                *S,
                pastix_int_t         lds )
{
    pastix_int_t *iparm;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_getSchur: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (S == NULL) {
        errorPrint("pastix_getSchur: S parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (lds <= 0) {
        errorPrint("pastix_getSchur: lds parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        errorPrint("pastix_getSchur: All steps from pastix_task_init() to pastix_task_numfact() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    iparm = pastix_data->iparm;
    switch(iparm[IPARM_FLOAT])
    {
    case PastixPattern:
        break;
    case PastixFloat:
        coeftab_sgetschur( pastix_data->solvmatr, S, lds );
        break;
    case PastixComplex32:
        coeftab_cgetschur( pastix_data->solvmatr, S, lds );
        break;
    case PastixComplex64:
        coeftab_zgetschur( pastix_data->solvmatr, S, lds );
        break;
    case PastixDouble:
    default:
        coeftab_dgetschur( pastix_data->solvmatr, S, lds );
    }
    return PASTIX_SUCCESS;
}

/**
 * @}
 */

