/**
 * @file pastix_zcores.h
 *
 * PaStiX kernel header.
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
 */
#ifndef _pastix_zcores_h_
#define _pastix_zcores_h_

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define pastix_cblk_lock( cblk_ )    pastix_atomic_lock( &((cblk_)->lock) )
#define pastix_cblk_unlock( cblk_ )  pastix_atomic_unlock( &((cblk_)->lock) )
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @addtogroup kernel_blas_lapack
 * @{
 *    This module contains all the BLAS and LAPACK-like kernels that are working
 *    on lapack layout matrices.
 *
 *    @name PastixComplex64 BLAS kernels
 *    @{
 */
void core_zplrnt( int m, int n, pastix_complex64_t *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed );
void core_zgetro( int m, int n,
                  const pastix_complex64_t *A, int lda,
                  pastix_complex64_t *B, int ldb );
int  core_zgeadd( pastix_trans_t trans, pastix_int_t M, pastix_int_t N,
                  pastix_complex64_t alpha, const pastix_complex64_t *A, pastix_int_t LDA,
                  pastix_complex64_t beta,        pastix_complex64_t *B, pastix_int_t LDB );
int  core_zgemdm( pastix_trans_t transA, pastix_trans_t transB, int M, int N, int K,
                  pastix_complex64_t  alpha, const pastix_complex64_t *A, int LDA,
                  const pastix_complex64_t *B, int LDB,
                  pastix_complex64_t  beta, pastix_complex64_t *C, int LDC,
                  const pastix_complex64_t *D, int incD,
                  pastix_complex64_t *WORK, int LWORK );
int  core_zrrqr ( pastix_int_t m, pastix_int_t n,
                  pastix_complex64_t *A, pastix_int_t lda,
                  pastix_int_t *jpvt, pastix_complex64_t *tau,
                  pastix_complex64_t *work, pastix_int_t ldwork, double *rwork,
                  double tol, pastix_int_t nb, pastix_int_t maxrank);
int  core_ztradd( pastix_uplo_t uplo, pastix_trans_t trans, pastix_int_t M, pastix_int_t N,
                  pastix_complex64_t alpha, const pastix_complex64_t *A, pastix_int_t LDA,
                  pastix_complex64_t beta,        pastix_complex64_t *B, pastix_int_t LDB);
int  core_zscalo( pastix_trans_t trans, pastix_int_t M, pastix_int_t N,
                  const pastix_complex64_t *A, pastix_int_t lda,
                  const pastix_complex64_t *D, pastix_int_t ldd,
                  pastix_complex64_t *B, pastix_int_t ldb );

/**
 *    @}
 *    @name PastixComplex64 LAPACK kernels
 *    @{
 */
void core_zpotrfsp( pastix_int_t n, pastix_complex64_t *A, pastix_int_t lda,
                    pastix_int_t *nbpivot, double criteria );
void core_zpxtrfsp( pastix_int_t n, pastix_complex64_t *A, pastix_int_t lda,
                    pastix_int_t *nbpivot, double criteria );
void core_zgetrfsp( pastix_int_t n, pastix_complex64_t *A, pastix_int_t lda,
                    pastix_int_t *nbpivot, double criteria );
void core_zhetrfsp( pastix_int_t n, pastix_complex64_t *A, pastix_int_t lda,
                    pastix_int_t *nbpivot, double criteria );
void core_zsytrfsp( pastix_int_t n, pastix_complex64_t *A, pastix_int_t lda,
                    pastix_int_t *nbpivot, double criteria );

/**
 *     @}
 * @}
 *
 * @addtogroup kernel_fact
 * @{
 *    This module contains all the kernel working at the solver matrix structure
 *    level for the numerical factorization step.
 *
 *    @name PastixComplex64 cblk-BLAS CPU kernels
 *    @{
 */

int  cpucblk_zgeaddsp1d( const SolverCblk *cblk1, SolverCblk *cblk2,
                         const pastix_complex64_t *L1, pastix_complex64_t *L2,
                         const pastix_complex64_t *U1, pastix_complex64_t *U2 );

void cpucblk_zgemmsp( pastix_coefside_t sideA, pastix_coefside_t sideB, pastix_trans_t trans,
                      const SolverCblk *cblk, const SolverBlok *blok, SolverCblk *fcblk,
                      const pastix_complex64_t *A, const pastix_complex64_t *B, pastix_complex64_t *C,
                      pastix_complex64_t *work, pastix_int_t lwork, const pastix_lr_t *lowrank );
void cpucblk_ztrsmsp( pastix_coefside_t coef, pastix_side_t side, pastix_uplo_t uplo,
                      pastix_trans_t trans, pastix_diag_t diag, SolverCblk *cblk,
                      const pastix_complex64_t *A, pastix_complex64_t *C, const pastix_lr_t *lowrank );
void cpucblk_zscalo ( pastix_trans_t trans, SolverCblk *cblk, pastix_complex64_t *LD );

void cpublok_zgemmsp( pastix_coefside_t sideA, pastix_coefside_t sideB, pastix_trans_t trans,
                      const SolverCblk *cblk, SolverCblk *fcblk,
                      pastix_int_t blok_mk, pastix_int_t blok_nk, pastix_int_t blok_mn,
                      const pastix_complex64_t *A, const pastix_complex64_t *B, pastix_complex64_t *C,
                      const pastix_lr_t *lowrank );
void cpublok_ztrsmsp( pastix_coefside_t coef, pastix_side_t side, pastix_uplo_t uplo,
                      pastix_trans_t trans, pastix_diag_t diag,
                      SolverCblk *cblk, pastix_int_t blok_m,
                      const pastix_complex64_t *A, pastix_complex64_t *C,
                      const pastix_lr_t *lowrank );
void cpublok_zscalo ( pastix_trans_t trans,
                      SolverCblk *cblk, pastix_int_t blok_m,
                      const pastix_complex64_t *A, const pastix_complex64_t *D, pastix_complex64_t *B );

/**
 *    @}
 *    @name PastixComplex64 cblk LU kernels
 *    @{
 */
int cpucblk_zgetrfsp1d_getrf( SolverMatrix *solvmtx, SolverCblk *cblk,
                              pastix_complex64_t *L, pastix_complex64_t *U );
int cpucblk_zgetrfsp1d_panel( SolverMatrix *solvmtx, SolverCblk *cblk,
                              pastix_complex64_t *L, pastix_complex64_t *U );
int cpucblk_zgetrfsp1d      ( SolverMatrix *solvmtx, SolverCblk *cblk,
                              pastix_complex64_t *work, pastix_int_t lwork );

/**
 *    @}
 *    @name PastixComplex64 cblk Cholesky kernels
 *    @{
 */
int cpucblk_zpotrfsp1d_potrf( SolverMatrix *solvmtx, SolverCblk *cblk, pastix_complex64_t *L );
int cpucblk_zpotrfsp1d_panel( SolverMatrix *solvmtx, SolverCblk *cblk, pastix_complex64_t *L );
int cpucblk_zpotrfsp1d      ( SolverMatrix *solvmtx, SolverCblk *cblk,
                              pastix_complex64_t *work, pastix_int_t lwork );

/**
 *    @}
 *    @name PastixComplex64 cblk LDL^h kernels
 *    @{
 */
int cpucblk_zhetrfsp1d_hetrf( SolverMatrix *solvmtx, SolverCblk *cblk,
                              pastix_complex64_t *L );
int cpucblk_zhetrfsp1d_panel( SolverMatrix *solvmtx, SolverCblk *cblk,
                              pastix_complex64_t *L, pastix_complex64_t *DLh );
int cpucblk_zhetrfsp1d      ( SolverMatrix *solvmtx, SolverCblk *cblk,
                              pastix_complex64_t *work1, pastix_complex64_t *work2, pastix_int_t lwork );

/**
 *    @}
 *    @name PastixComplex64 cblk LL^t kernels
 *    @{
 */
int cpucblk_zpxtrfsp1d_pxtrf( SolverMatrix *solvmtx, SolverCblk *cblk, pastix_complex64_t *L );
int cpucblk_zpxtrfsp1d_panel( SolverMatrix *solvmtx, SolverCblk *cblk, pastix_complex64_t *L );
int cpucblk_zpxtrfsp1d      ( SolverMatrix *solvmtx, SolverCblk *cblk,
                              pastix_complex64_t *work, pastix_int_t lwork );

/**
 *    @}
 *    @name PastixComplex64 cblk LDL^t kernels
 *    @{
 */
int cpucblk_zsytrfsp1d_sytrf( SolverMatrix *solvmtx, SolverCblk *cblk,
                              pastix_complex64_t *L );
int cpucblk_zsytrfsp1d_panel( SolverMatrix *solvmtx, SolverCblk *cblk,
                              pastix_complex64_t *L, pastix_complex64_t *DLt );
int cpucblk_zsytrfsp1d      ( SolverMatrix *solvmtx, SolverCblk *cblk,
                              pastix_complex64_t *work1, pastix_complex64_t *work2, pastix_int_t lwork );

/**
 *    @}
 *    @name PastixComplex64 initialization and additionnal routines
 *    @{
 */
void cpucblk_zalloc   ( pastix_coefside_t    side,
                        SolverCblk          *cblk );
void cpucblk_zfillin  ( pastix_coefside_t    side,
                        const SolverMatrix  *solvmtx,
                        const pastix_bcsc_t *bcsc,
                        pastix_int_t         itercblk );
void cpucblk_zinit    ( pastix_coefside_t    side,
                        const SolverMatrix  *solvmtx,
                        const pastix_bcsc_t *bcsc,
                        pastix_int_t         itercblk,
                        char               **directory );
void cpucblk_zgetschur( const SolverCblk    *cblk,
                        int                  upper_part,
                        pastix_complex64_t  *S,
                        pastix_int_t         lds );
void cpucblk_zdump    ( pastix_coefside_t    side,
                        const SolverCblk    *cblk,
                        FILE                *stream );
int  cpucblk_zdiff    ( pastix_coefside_t    side,
                        const SolverCblk    *cblkA,
                        SolverCblk          *cblkB );
/**
 *    @}
 *    @name PastixComplex64 compression/uncompression routines
 *    @{
 */
pastix_int_t cpucblk_zcompress  ( pastix_coefside_t side,
                                  SolverCblk       *cblk,
                                  pastix_lr_t       lowrank );
void         cpucblk_zuncompress( pastix_coefside_t side,
                                  SolverCblk       *cblk );
void         cpucblk_zmemory    ( pastix_coefside_t   side,
                                  const SolverMatrix *solvmtx,
                                  const SolverCblk   *cblk,
                                  pastix_int_t       *gain );

/**
 *     @}
 * @}
 *
 * @addtogroup kernel_solve
 * @{
 *    This module contains all the kernel working on the solver matrix structure
 *    for the solve step.
 */

void solve_blok_ztrsm( pastix_coefside_t coefside, pastix_side_t side, pastix_uplo_t uplo,
                       pastix_trans_t trans, pastix_diag_t diag, const SolverCblk *cblk,
                       int nrhs, pastix_complex64_t *b, int ldb );
void solve_blok_zgemm( pastix_coefside_t coefside, pastix_side_t side, pastix_trans_t trans,
                       pastix_int_t nrhs, const SolverCblk *cblk, const SolverBlok *blok,
                       SolverCblk *fcbk, const pastix_complex64_t *B, pastix_int_t ldb,
                       pastix_complex64_t *C, pastix_int_t ldc );

void solve_cblk_ztrsmsp_forward( pastix_solv_mode_t mode, pastix_side_t side, pastix_uplo_t uplo,
                                 pastix_trans_t trans, pastix_diag_t diag,
                                 const SolverMatrix *datacode, const SolverCblk *cblk,
                                 int nrhs, pastix_complex64_t *b, int ldb );
void solve_cblk_ztrsmsp_backward( pastix_solv_mode_t mode, pastix_side_t side, pastix_uplo_t uplo,
                                  pastix_trans_t trans, pastix_diag_t diag,
                                  const SolverMatrix *datacode, const SolverCblk *cblk,
                                  int nrhs, pastix_complex64_t *b, int ldb );

void solve_cblk_zdiag( const SolverCblk   *cblk,
                       int                 nrhs,
                       pastix_complex64_t *b,
                       int                 ldb,
                       pastix_complex64_t *work );
/**
 * @}
 */


/**
 * @addtogroup kernel_fact_null
 * @{
 *    This module contains the three terms update functions for the LDL^t and
 *    LDL^h factorizations.
 */
void core_zhetrfsp1d_gemm( const SolverCblk *cblk, const SolverBlok *blok, SolverCblk *fcblk,
                           const pastix_complex64_t *L, pastix_complex64_t *C,
                           pastix_complex64_t *work );
void core_zsytrfsp1d_gemm( const SolverCblk *cblk, const SolverBlok *blok, SolverCblk *fcblk,
                           const pastix_complex64_t *L, pastix_complex64_t *C,
                           pastix_complex64_t *work );

/**
 * @}
 */

#endif /* _pastix_zcores_h_ */
