/**
 *
 * @file core_zgetrfsp.c
 *
 * PaStiX kernel routines for LU factorization.
 *
 * @copyright 2011-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2018-07-16
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "cblas.h"
#include "blend/solver.h"
#include "pastix_zcores.h"
#include "kernels_trace.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define MAXSIZEOFBLOCKS 64
static pastix_complex64_t zone  =  1.0;
static pastix_complex64_t mzone = -1.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup kernel_blas_lapack_null
 *
 * @brief Compute the sequential static pivoting LU factorization of the matrix
 * m-by-n A = L * U.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows and columns of the matrix A.
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[inout] A
 *          The matrix A to factorize with LU factorization. The matrix
 *          is of size lda -by- n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *
 * @param[inout] nbpivots
 *          Pointer to the number of piovting operations made during
 *          factorization. It is updated during this call
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the number of
 *          pivots is incremented.
 *
 *******************************************************************************/
static inline void
core_zgetf2sp( pastix_int_t        m,
               pastix_int_t        n,
               pastix_complex64_t *A,
               pastix_int_t        lda,
               pastix_int_t       *nbpivots,
               double              criteria )
{
    pastix_int_t k, minMN;
    pastix_complex64_t *Akk, *Aik, alpha;

    minMN = pastix_imin( m, n );

    Akk = A;
    for (k=0; k<minMN; k++) {
        Aik = Akk + 1;

        if ( cabs(*Akk) < criteria ) {
            if ( creal(*Akk) < 0. ) {
                *Akk = (pastix_complex64_t)(-criteria);
            }
            else {
                *Akk = (pastix_complex64_t)criteria;
            }
            (*nbpivots)++;
        }

        /* A_ik = A_ik / A_kk, i = k+1 .. n */
        alpha = 1.0 / (*Akk);
        cblas_zscal(m-k-1, CBLAS_SADDR( alpha ), Aik, 1 );

        if ( k+1 < minMN ) {

            /* A_ij = A_ij - A_ik * A_kj, i,j = k+1..n */
            cblas_zgeru(CblasColMajor, m-k-1, n-k-1,
                        CBLAS_SADDR(mzone),
                        Aik,        1,
                        Akk+lda, lda,
                        Aik+lda, lda);
        }

        Akk += lda+1;
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the block static pivoting LU factorization of the matrix
 * m-by-n A = L * U.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[inout] A
 *          The matrix A to factorize with LU factorization. The matrix
 *          is of size lda -by- n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *
 * @param[inout] nbpivots
 *          Pointer to the number of piovting operations made during
 *          factorization. It is updated during this call
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the number of
 *          pivots is incremented.
 *
 *******************************************************************************/
void
core_zgetrfsp( pastix_int_t        n,
               pastix_complex64_t *A,
               pastix_int_t        lda,
               pastix_int_t       *nbpivots,
               double              criteria )
{
    pastix_int_t k, blocknbr, blocksize, matrixsize, tempm;
    pastix_complex64_t *Akk, *Lik, *Ukj, *Aij;

    blocknbr = pastix_iceil( n, MAXSIZEOFBLOCKS );

    Akk = A; /* Lk,k     */

    for (k=0; k<blocknbr; k++) {

        tempm = n - k * MAXSIZEOFBLOCKS;
        blocksize = pastix_imin(MAXSIZEOFBLOCKS, tempm);
        Lik = Akk + blocksize;
        Ukj = Akk + blocksize*lda;
        Aij = Ukj + blocksize;

        /* Factorize the diagonal block Akk*/
        core_zgetf2sp( tempm, blocksize, Akk, lda, nbpivots, criteria );

        matrixsize = tempm - blocksize;
        if ( matrixsize > 0 ) {

            /* Compute the column Ukk+1 */
            cblas_ztrsm(CblasColMajor,
                        CblasLeft, CblasLower,
                        CblasNoTrans, CblasUnit,
                        blocksize, matrixsize,
                        CBLAS_SADDR(zone), Akk, lda,
                                           Ukj, lda);

            /* Update Ak+1,k+1 = Ak+1,k+1 - Lk+1,k*Uk,k+1 */
            cblas_zgemm(CblasColMajor,
                        CblasNoTrans, CblasNoTrans,
                        matrixsize, matrixsize, blocksize,
                        CBLAS_SADDR(mzone), Lik, lda,
                                            Ukj, lda,
                        CBLAS_SADDR(zone),  Aij, lda);
        }

        Akk += blocksize * (lda+1);
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the LU factorization of the diagonal block in a panel.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          Solver Matrix structure of the problem
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[inout] L
 *          The pointer to the lower matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] U
 *          The pointer to the upper matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 *******************************************************************************
 *
 * @return The number of static pivoting performed during the diagonal block
 *         factorization.
 *
 *******************************************************************************/
int
cpucblk_zgetrfsp1d_getrf( SolverMatrix       *solvmtx,
                          SolverCblk         *cblk,
                          pastix_complex64_t *L,
                          pastix_complex64_t *U )
{
    pastix_int_t ncols, stride;
    pastix_int_t nbpivots = 0;
    pastix_fixdbl_t time, flops;
    double criteria = solvmtx->diagthreshold;

    time = kernel_trace_start( PastixKernelGETRF );

    ncols  = cblk->lcolnum - cblk->fcolnum + 1;
    stride = (cblk->cblktype & CBLK_LAYOUT_2D) ? ncols : cblk->stride;

    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        assert( cblk->fblokptr->LRblock[0].rk == -1 &&
                cblk->fblokptr->LRblock[1].rk == -1 );
        L = cblk->fblokptr->LRblock[0].u;
        U = cblk->fblokptr->LRblock[1].u;
        stride = ncols;

        assert( stride == cblk->fblokptr->LRblock[0].rkmax );
        assert( stride == cblk->fblokptr->LRblock[1].rkmax );
    }

    core_zgeadd( PastixTrans, ncols, ncols,
                 1.0, U, stride,
                 1.0, L, stride );

    /* Factorize diagonal block */
    flops = FLOPS_ZGETRF( ncols, ncols );
    kernel_trace_start_lvl2( PastixKernelLvl2GETRF );
    core_zgetrfsp(ncols, L, stride, &nbpivots, criteria);
    kernel_trace_stop_lvl2( flops );

    /* Transpose Akk in ucoeftab */
    core_zgetro( ncols, ncols, L, stride, U, stride );

    kernel_trace_stop( PastixKernelGETRF, ncols, 0, 0, flops, time );

    if ( nbpivots ) {
        pastix_atomic_add_32b( &(solvmtx->nbpivots), nbpivots );
    }
    return nbpivots;
}

/**
 *******************************************************************************
 *
 * @brief Compute the LU factorization of one panel.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          Solver Matrix structure of the problem
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[inout] L
 *          The pointer to the lower matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] U
 *          The pointer to the upper matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 *******************************************************************************
 *
 * @return The number of static pivoting performed during the diagonal block
 *         factorization.
 *
 *******************************************************************************/
int
cpucblk_zgetrfsp1d_panel( SolverMatrix       *solvmtx,
                          SolverCblk         *cblk,
                          pastix_complex64_t *L,
                          pastix_complex64_t *U )
{
    pastix_int_t nbpivots;
    nbpivots = cpucblk_zgetrfsp1d_getrf( solvmtx, cblk, L, U );

    /*
     * We exploit the fact that the upper triangle is stored at the top of the L
     * column, and by transposition the L part of the diagonal block is
     * similarly stored in the U panel
     */
    cpucblk_ztrsmsp( PastixLCoef, PastixRight, PastixUpper,
                     PastixNoTrans, PastixNonUnit,
                     cblk, L, L, &(solvmtx->lowrank) );
    cpucblk_ztrsmsp( PastixUCoef, PastixRight, PastixUpper,
                     PastixNoTrans, PastixUnit,
                     cblk, U, U, &(solvmtx->lowrank) );
    return nbpivots;
}

/**
 *******************************************************************************
 *
 * @brief Perform the LU factorization of a given panel and apply all its
 * updates.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          Solver Matrix structure of the problem
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[in] work
 *          Temporary memory buffer.
 *
 * @param[in] lwork
 *          Temporary workspace dimension.
 *
 *******************************************************************************
 *
 * @return The number of static pivoting during factorization of the diagonal
 * block.
 *
 *******************************************************************************/
int
cpucblk_zgetrfsp1d( SolverMatrix       *solvmtx,
                    SolverCblk         *cblk,
                    pastix_complex64_t *work,
                    pastix_int_t        lwork )
{
    pastix_complex64_t *L = cblk->lcoeftab;
    pastix_complex64_t *U = cblk->ucoeftab;
    SolverCblk  *fcblk;
    SolverBlok  *blok, *lblk;
    pastix_int_t nbpivots;

    nbpivots = cpucblk_zgetrfsp1d_panel( solvmtx, cblk, L, U );

    blok = cblk->fblokptr + 1; /* this diagonal block */
    lblk = cblk[1].fblokptr;   /* the next diagonal block */

    /* if there are off-diagonal supernodes in the column */
    for( ; blok < lblk; blok++ )
    {
        fcblk = (solvmtx->cblktab + blok->fcblknm);

        /* Update on L */
        cpucblk_zgemmsp( PastixLCoef, PastixUCoef, PastixTrans,
                         cblk, blok, fcblk,
                         L, U, fcblk->lcoeftab,
                         work, lwork, &(solvmtx->lowrank) );

        /* Update on U */
        if ( blok+1 < lblk ) {
            cpucblk_zgemmsp( PastixUCoef, PastixLCoef, PastixTrans,
                             cblk, blok, fcblk,
                             U, L, fcblk->ucoeftab,
                             work, lwork, &(solvmtx->lowrank) );
        }
        pastix_atomic_dec_32b( &(fcblk->ctrbcnt) );
    }

    return nbpivots;
}
