/**
 *
 * @file sopalin/diag.c
 *
 * PaStiX diagonal interface functions
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
 * @addtogroup pastix
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
 * @brief Return the diagonal of the matrix.
 *
 * The diagonal is returned in the D vector such that D[incD *i] = A(i, i)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix data structure of the problem solved.
 *
 * @param[inout] D
 *          The pointer to the allocated vector array that will store the diagonal.
 *          D must be of size spm->nexp * incD.
 *
 * @param[in] incD
 *          The leading dimension of the D array.
 *
 ********************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastixGetDiag( const pastix_data_t *pastix_data,
               void                *D,
               pastix_int_t         incD )
{
    pastix_int_t *iparm;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_getDiag: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (D == NULL) {
        errorPrint("pastix_getDiag: D parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (incD <= 0) {
        errorPrint("pastix_getDiag: incD parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        errorPrint("pastix_getDiag: All steps from pastix_task_init() to pastix_task_numfact() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    iparm = pastix_data->iparm;
    switch(iparm[IPARM_FLOAT])
    {
    case PastixPattern:
        break;
    case PastixFloat:
        coeftab_sgetdiag( pastix_data->solvmatr, D, incD );
        break;
    case PastixComplex32:
        coeftab_cgetdiag( pastix_data->solvmatr, D, incD );
        break;
    case PastixComplex64:
        coeftab_zgetdiag( pastix_data->solvmatr, D, incD );
        break;
    case PastixDouble:
    default:
        coeftab_dgetdiag( pastix_data->solvmatr, D, incD );
    }
    return PASTIX_SUCCESS;
}

/**
 * @}
 */

