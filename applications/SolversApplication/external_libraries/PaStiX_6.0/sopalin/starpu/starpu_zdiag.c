/**
 *
 * @file starpu_zdiag.c
 *
 * PaStiX zdiag StarPU wrapper.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Vincent Bridonneau
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-16
 * @precisions normal z -> s d c
 *
 * @addtogroup starpu_diag_solve
 * @{
 *
 **/

#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"
#include "pastix_starpu.h"
#include "pastix_zstarpu.h"

/**
 *******************************************************************************
 *
 * @brief Apply the diagonal solve on a problem of 1 dimension.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          Provide informations about starpu and the schur solving mode.
 *
 * @param[in] sopalin_data
 *          Provide descriptor of b and the sparse matrix and the SolverMatrix
 *          structure.
 *
 *******************************************************************************/
void
starpu_zdiag_sp1dplus( pastix_data_t               *pastix_data,
                       sopalin_data_t              *sopalin_data )
{
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    pastix_int_t        k, cblknbr;
    pastix_solv_mode_t  mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];

    cblk = datacode->cblktab;
    cblknbr = (mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;
    for (k=0; k<cblknbr; k++, cblk++) {
        starpu_stask_cblk_zdiag( sopalin_data, cblk, cblknbr - k);
    }
}

/**
 *******************************************************************************
 *
 * @brief Apply the diagonal solve (StarPU version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          Provide informations about starpu and the schur solving mode.
 *
 * @param[in] sopalin_data
 *          Provide descriptor of b and the sparse matrix and the SolverMatrix
 *          structure.
 *
 * @param[in] nrhs
 *          The number of right hand side
 *
 * @param[inout] b
 *          The pointer to vectors of the right hand side
 *
 * @param[in] ldb
 *          The leading dimension of b
 *
 *******************************************************************************/
void
starpu_zdiag( pastix_data_t      *pastix_data,
              sopalin_data_t     *sopalin_data,
              int                 nrhs,
              pastix_complex64_t *b,
              int                 ldb )
{
    starpu_sparse_matrix_desc_t *sdesc = sopalin_data->solvmtx->starpu_desc;
    starpu_dense_matrix_desc_t  *ddesc;

    /*
     * Start StarPU if not already started
     */
    if (pastix_data->starpu == NULL) {
        int argc = 0;
        pastix_starpu_init( pastix_data, &argc, NULL, NULL );
    }

    if ( sdesc == NULL ) {
        /* Create the sparse matrix descriptor */
        starpu_sparse_matrix_init( sopalin_data->solvmtx,
                                   sizeof( pastix_complex64_t ), PastixHermitian,
                                   1, 0 );
        sdesc = sopalin_data->solvmtx->starpu_desc;
    }

    /* Create the dense matrix descriptor */
    starpu_dense_matrix_init( sopalin_data->solvmtx,
                              nrhs, (char*)b, ldb,
                              sizeof(pastix_complex64_t), 1, 0 );
    ddesc = sopalin_data->solvmtx->starpu_desc_rhs;

    starpu_resume();
    starpu_zdiag_sp1dplus( pastix_data, sopalin_data );

    starpu_sparse_matrix_getoncpu( sdesc );
    starpu_dense_matrix_getoncpu( ddesc );
    starpu_task_wait_for_all();
#if defined(PASTIX_WITH_MPI)
    starpu_mpi_barrier(MPI_COMM_WORLD);
#endif
    starpu_pause();

    return;
}

/**
 *@}
 */
