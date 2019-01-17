/**
 *
 * @file starpu_zpotrf.c
 *
 * PaStiX zpotrf StarPU wrapper.
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-16
 * @precisions normal z -> s d c
 *
 * @addtogroup starpu_potrf
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
 * @brief Perform a sparse Cholesky factorization with 1D kernels.
 *
 * The function performs the Cholesky factorization of a sparse symmetric
 * positive definite (or Hermitian positive definite in the complex case) matrix
 * A.
 * The factorization has the form
 *
 *    \f[ A = L\times L^H \f]
 *
 * where L is a sparse lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[inout] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 * @param[inout] desc
 *          StarPU descriptor of the sparse matrix.
 *
 ******************************************************************************/
void
starpu_zpotrf_sp1dplus( sopalin_data_t              *sopalin_data,
                        starpu_sparse_matrix_desc_t *desc )
{
    const SolverMatrix *solvmtx = sopalin_data->solvmtx;
    SolverCblk         *cblk, *fcblk;
    SolverBlok         *blok, *lblk;
    pastix_int_t        k, m, cblknbr, cblk_n;

    cblknbr = solvmtx->cblknbr;
    cblk = solvmtx->cblktab;
    for (k=0; k<solvmtx->cblknbr; k++, cblk++){

        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            break;
        }

        starpu_task_cblk_zpotrfsp1d_panel( sopalin_data, cblk,
                                           cblknbr - k );

        blok = cblk->fblokptr + 1; /* this diagonal block */
        lblk = cblk[1].fblokptr;   /* the next diagonal block */

        /* if there are off-diagonal supernodes in the column */
        for(m=0; blok < lblk; blok++, m++ )
        {
            fcblk = (solvmtx->cblktab + blok->fcblknm);
            cblk_n = fcblk - solvmtx->cblktab;

            starpu_task_cblk_zgemmsp( PastixLCoef, PastixLCoef, PastixConjTrans,
                                      cblk, blok, fcblk, sopalin_data,
                                      cblknbr - pastix_imin( k + m, cblk_n ) );
        }
    }
    (void)desc;
}

/**
 *******************************************************************************
 *
 * @brief Perform a sparse Cholesky factorization with 1D and 2D kernels.
 *
 * The function performs the Cholesky factorization of a sparse symmetric
 * positive definite (or Hermitian positive definite in the complex case) matrix
 * A.
 * The factorization has the form
 *
 *    \f[ A = L\times L^H \f]
 *
 * where L is a sparse lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[inout] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 * @param[inout] desc
 *          StarPU descriptor of the sparse matrix.
 *
 ******************************************************************************/
void
starpu_zpotrf_sp2d( sopalin_data_t              *sopalin_data,
                    starpu_sparse_matrix_desc_t *desc )
{
    const SolverMatrix *solvmtx = sopalin_data->solvmtx;
    SolverCblk         *cblk, *fcblk;
    SolverBlok         *blok, *lblk, *blokA, *blokB;
    starpu_cblk_t      *cblkhandle;
    pastix_int_t        k, m, cblknbr, cblk_n;

    cblknbr = solvmtx->cblknbr;

    /* Let's submit all 1D tasks first */
    cblk = solvmtx->cblktab;
    for (k=0; k<=solvmtx->cblkmax1d; k++, cblk++){

        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            break;
        }

        if ( cblk->cblktype & CBLK_TASKS_2D ) {
            assert( k >= solvmtx->cblkmin2d );
            cblkhandle = desc->cblktab_handle + (k - solvmtx->cblkmin2d);
            starpu_data_partition_submit( cblk->handler[0],
                                          cblkhandle->handlenbr,
                                          cblkhandle->handletab );
            continue;
        }

        starpu_task_cblk_zpotrfsp1d_panel( sopalin_data, cblk,
                                           cblknbr - k );

        blok  = cblk->fblokptr + 1; /* this diagonal block */
        lblk = cblk[1].fblokptr;   /* the next diagonal block */

        /* if there are off-diagonal supernodes in the column */
        for(m=0; blok < lblk; blok++, m++ )
        {
            fcblk = (solvmtx->cblktab + blok->fcblknm);
            cblk_n = fcblk - solvmtx->cblktab;

            starpu_task_cblk_zgemmsp( PastixLCoef, PastixLCoef, PastixConjTrans,
                                      cblk, blok, fcblk, sopalin_data,
                                      cblknbr - pastix_imin( k + m, cblk_n ) );
        }
    }

    /* Let's submit the partitionning */
    k = solvmtx->cblkmax1d + 1;
    cblk       = solvmtx->cblktab + k;
    cblkhandle = desc->cblktab_handle + (k - solvmtx->cblkmin2d);
    for (; k<solvmtx->cblknbr; k++, cblk++, cblkhandle++) {

        if ( !(cblk->cblktype & CBLK_TASKS_2D) ) {
            continue;
        }

        starpu_data_partition_submit( cblk->handler[0],
                                      cblkhandle->handlenbr,
                                      cblkhandle->handletab );
    }

    /* Now we submit all 2D tasks */
    cblk       = solvmtx->cblktab + solvmtx->cblkmin2d;
    cblkhandle = desc->cblktab_handle;
    for (k=solvmtx->cblkmin2d; k<solvmtx->cblknbr; k++, cblk++, cblkhandle++){

        if ( !(cblk->cblktype & CBLK_TASKS_2D) ) {
            continue; /* skip 1D cblk */
        }

        if (cblk->cblktype & CBLK_IN_SCHUR)
        {
            starpu_data_unpartition_submit( cblk->handler[0],
                                            cblkhandle->handlenbr,
                                            cblkhandle->handletab, STARPU_MAIN_RAM );
            continue;
        }

        starpu_task_blok_zpotrf( sopalin_data, cblk,
                                 cblknbr - k );

        lblk = cblk[1].fblokptr;
        for(blokA=cblk->fblokptr + 1, m=0; blokA<lblk; blokA++, m++) {

            cblk_n = blokA->fcblknm;

            starpu_task_blok_ztrsmsp( PastixLCoef, PastixRight, PastixLower,
                                      PastixConjTrans, PastixNonUnit,
                                      cblk, blokA, sopalin_data,
                                      cblknbr - k );

            for(blokB=cblk->fblokptr + 1; blokB<=blokA; blokB++) {

                starpu_task_blok_zgemmsp( PastixLCoef, PastixLCoef, PastixConjTrans,
                                          cblk, solvmtx->cblktab + blokB->fcblknm,
                                          blokA, blokB, sopalin_data,
                                          cblknbr - pastix_imin( k + m, cblk_n ) );

                /* Skip B blocks facing the same cblk */
                while( (blokB < blokA) &&
                       (blokB[0].fcblknm == blokB[1].fcblknm) &&
                       (blokB[0].lcblknm == blokB[1].lcblknm) )
                {
                    blokB++;
                }
            }

            /* Skip A blocks facing the same cblk */
            while( (blokA < lblk) &&
                   (blokA[0].fcblknm == blokA[1].fcblknm) &&
                   (blokA[0].lcblknm == blokA[1].lcblknm) )
            {
                blokA++;
            }
        }
        starpu_data_unpartition_submit( cblk->handler[0],
                                        cblkhandle->handlenbr,
                                        cblkhandle->handletab, STARPU_MAIN_RAM );
    }
    (void)desc;
}

/**
 *******************************************************************************
 *
 * @brief Perform a sparse Cholesky factorization using StarPU runtime.
 *
 * The function performs the Cholesky factorization of a sparse symmetric
 * positive definite (or Hermitian positive definite in the complex case) matrix
 * A.
 * The factorization has the form
 *
 *    \f[ A = L\times L^H \f]
 *
 * where L is a sparse lower triangular matrix.
 *
 * The algorithm is automatically chosen between the 1D and 2D version based on
 * the API parameter IPARM_TASKS2D_LEVEL. If IPARM_TASKS2D_LEVEL != 0
 * the 2D scheme is applied, the 1D otherwise.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[inout] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 ******************************************************************************/
void
starpu_zpotrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    starpu_sparse_matrix_desc_t *sdesc = sopalin_data->solvmtx->starpu_desc;

    /*
     * Start StarPU if not already started
     */
    if (pastix_data->starpu == NULL) {
        int argc = 0;
        pastix_starpu_init( pastix_data, &argc, NULL, NULL );
    }

    if ( sdesc == NULL ) {
        /* Create the matrix descriptor */
        starpu_sparse_matrix_init( sopalin_data->solvmtx,
                                   sizeof( pastix_complex64_t ), PastixHermitian,
                                   1, 0 );
        sdesc = sopalin_data->solvmtx->starpu_desc;
    }

    starpu_resume();
    /*
     * Select 1D or 2D algorithm based on 2d tasks level
     */
    if ( pastix_data->iparm[IPARM_TASKS2D_LEVEL] != 0 )
    {
        starpu_zpotrf_sp2d( sopalin_data, sdesc );
    }
    else
    {
        starpu_zpotrf_sp1dplus( sopalin_data, sdesc );
    }

    starpu_sparse_matrix_getoncpu( sdesc );
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
