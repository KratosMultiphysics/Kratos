/**
 *
 * @file pastix_task_solve.c
 *
 *  PaStiX solve routines
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
#include "bcsc.h"
#include "pastix/order.h"
#include "solver.h"
#include "sopalin_data.h"

#include "bcsc_z.h"
#include "bcsc_c.h"
#include "bcsc_d.h"
#include "bcsc_s.h"

#if defined(PASTIX_DEBUG_SOLVE)
static inline void
dump_rhs( char *name, int n, double *b )
{
    int i;
    fprintf(stderr,"%s :", name );
    for (i=0; i<n; i++) {
        if (i%10 == 0)
            fprintf(stderr, "\n");
        fprintf(stderr,"%e ", b[i]);
    }
    fprintf(stderr,"\n");
}
#else
static inline void
dump_rhs( char *name, int n, double *b )
{
    (void)name;
    (void)n;
    (void)b;
}
#endif

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Apply a permutation on the right-and-side vector before the solve step.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION, IPARM_APPLYPERM_WS.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[in] flttype
 *          This arithmetic of the sparse matrix.
 *
 * @param[in] dir
 *          Forward or backword application of the permutation.
 *
 * @param[in] m
 *          Size of the right-and-side vectors.
 *
 * @param[in] n
 *          Number of right-and-side vectors.
 *
 * @param[inout] b
 *          The right-and-side vectors (can be multiple RHS).
 *
 * @param[in] ldb
 *          The leading dimension of the right-and-side vectors.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastix_subtask_applyorder( pastix_data_t *pastix_data,
                           pastix_coeftype_t flttype, pastix_dir_t dir,
                           pastix_int_t m, pastix_int_t n, void *b, pastix_int_t ldb )
{
    pastix_int_t *perm;
    int ts;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_subtask_applyorder: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (b == NULL) {
        errorPrint("pastix_subtask_applyorder: wrong b parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_CSC2BCSC) ) {
        errorPrint("pastix_subtask_applyorder: All steps from pastix_task_init() to pastix_subtask_csc2bcsc() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Make sure ordering is 0 based */
    if ( pastix_data->ordemesh->baseval != 0 ) {
        errorPrint("pastix_subtask_applyorder: ordermesh must be 0-based");
        return PASTIX_ERR_BADPARAMETER;
    }

    ts   = pastix_data->iparm[IPARM_APPLYPERM_WS];
    perm = pastix_data->ordemesh->peritab;

    /* See also xlapmr and xlapmt */
    switch( flttype ) {
    case PastixComplex64:
        bvec_zlapmr( ts, dir, m, n, b, ldb, perm );
        break;

    case PastixComplex32:
        bvec_clapmr( ts, dir, m, n, b, ldb, perm );
        break;

    case PastixFloat:
        bvec_slapmr( ts, dir, m, n, b, ldb, perm );
        break;

    case PastixDouble:
    default:
        bvec_dlapmr( ts, dir, m, n, b, ldb, perm );
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Apply a triangular solve on the right-and-side vectors.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[in] flttype
 *          This arithmetic of the sparse matrix.
 *
 * @param[in] side
 *          Left or right application.
 *
 * @param[in] uplo
 *          Upper or Lower part.
 *
 * @param[in] trans
 *          With or without transposition (or conjugate transposition).
 *
 * @param[in] diag
 *          Diagonal terms are unit or not.
 *
 * @param[in] nrhs
 *          The number of right-and-side vectors.
 *
 * @param[inout] b
 *          The right-and-side vector (can be multiple RHS).
 *          On exit, the solution is stored in place of the right-hand-side vector.
 *
 * @param[in] ldb
 *          The leading dimension of the right-and-side vectors.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastix_subtask_trsm( pastix_data_t *pastix_data,
                     pastix_coeftype_t flttype, pastix_side_t side,
                     pastix_uplo_t uplo, pastix_trans_t trans, pastix_diag_t diag,
                     pastix_int_t nrhs, void *b, pastix_int_t ldb )
{
    sopalin_data_t sopalin_data;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_subtask_trsm: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (b == NULL) {
        errorPrint("pastix_subtask_trsm: wrong b parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        errorPrint("pastix_subtask_trsm: All steps from pastix_task_init() to pastix_task_numfact() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    sopalin_data.solvmtx = pastix_data->solvmatr;

    switch (flttype) {
    case PastixComplex64:
        sopalin_ztrsm( pastix_data, side, uplo, trans, diag,
                       &sopalin_data, nrhs, (pastix_complex64_t *)b, ldb );
        break;
    case PastixComplex32:
        sopalin_ctrsm( pastix_data, side, uplo, trans, diag,
                       &sopalin_data, nrhs, (pastix_complex32_t *)b, ldb );
        break;
    case PastixDouble:
        trans = (trans == PastixConjTrans) ? PastixTrans : trans;
        sopalin_dtrsm( pastix_data, side, uplo, trans, diag,
                       &sopalin_data, nrhs, (double *)b, ldb );
        break;
    case PastixFloat:
        trans = (trans == PastixConjTrans) ? PastixTrans : trans;
        sopalin_strsm( pastix_data, side, uplo, trans, diag,
                       &sopalin_data, nrhs, (float *)b, ldb );
        break;
    default:
        fprintf(stderr, "Unknown floating point arithmetic\n" );
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Apply a diagonal operation on the right-and-side vectors.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[in] flttype
 *          This arithmetic of the sparse matrix.
 *
 * @param[in] nrhs
 *          The number of right-and-side vectors.
 *
 * @param[inout] b
 *          The right-and-side vector (can be multiple RHS).
 *          On exit, the solution is stored in place of the right-hand-side vector.
 *
 * @param[in] ldb
 *          The leading dimension of the right-and-side vectors.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastix_subtask_diag( pastix_data_t *pastix_data, pastix_coeftype_t flttype,
                     pastix_int_t nrhs, void *b, pastix_int_t ldb )
{
    sopalin_data_t sopalin_data;
    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_subtask_diag: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (b == NULL) {
        errorPrint("pastix_subtask_diag: wrong b parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        errorPrint("pastix_subtask_trsm: All steps from pastix_task_init() to pastix_task_numfact() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    sopalin_data.solvmtx = pastix_data->solvmatr;

    switch (flttype) {
    case PastixComplex64:
        sopalin_zdiag( pastix_data, &sopalin_data, nrhs, (pastix_complex64_t *)b, ldb );
        break;
    case PastixComplex32:
        sopalin_cdiag( pastix_data, &sopalin_data, nrhs, (pastix_complex32_t *)b, ldb );
        break;
    case PastixDouble:
        sopalin_ddiag( pastix_data, &sopalin_data, nrhs, (double *)b, ldb );
        break;
    case PastixFloat:
        sopalin_sdiag( pastix_data, &sopalin_data, nrhs, (float *)b, ldb );
        break;
    default:
        fprintf(stderr, "Unknown floating point arithmetic\n" );
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Solve the given problem without applying the permutation.
 *
 * @warning The input vector is considered already permuted. For a solve step
 * with permutation, see pastix_task_solve()
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[in] nrhs
 *          The number of right-and-side vectors.
 *
 * @param[inout] b
 *          The right-and-side vectors (can be multiple RHS).
 *          On exit, the solution is stored in place of the right-hand-side vector.
 *
 * @param[in] ldb
 *          The leading dimension of the right-and-side vectors.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastix_subtask_solve( pastix_data_t *pastix_data,
                      pastix_int_t nrhs, void *b, pastix_int_t ldb )
{
    pastix_int_t  *iparm;
    pastix_bcsc_t *bcsc;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_task_solve: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        errorPrint("pastix_task_solve: All steps from pastix_task_init() to pastix_task_numfact() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    iparm = pastix_data->iparm;
    bcsc  = pastix_data->bcsc;

    {
        double timer;
        pastix_trans_t trans = PastixTrans;

        clockStart(timer);
        switch ( iparm[IPARM_FACTORIZATION] ){
        case PastixFactLLH:
            trans = PastixConjTrans;

            pastix_attr_fallthrough;

        case PastixFactLLT:
            dump_rhs( "AfterPerm", bcsc->gN, b );

            /* Solve L y = P b with y = L^t P x */
            pastix_subtask_trsm( pastix_data, bcsc->flttype,
                                 PastixLeft, PastixLower,
                                 PastixNoTrans, PastixNonUnit,
                                 nrhs, b, ldb );
            dump_rhs( "AfterDown", bcsc->gN, b );

            /* Solve y = L^t (P x) */
            pastix_subtask_trsm( pastix_data, bcsc->flttype,
                                 PastixLeft, PastixLower,
                                 trans, PastixNonUnit,
                                 nrhs, b, ldb );
            dump_rhs( "AfterUp", bcsc->gN, b );
            break;

        case PastixFactLDLH:
            trans = PastixConjTrans;

            pastix_attr_fallthrough;

        case PastixFactLDLT:
            dump_rhs( "AfterPerm", bcsc->gN, b );

            /* Solve L y = P b with y = D L^t P x */
            pastix_subtask_trsm( pastix_data, bcsc->flttype,
                                 PastixLeft, PastixLower,
                                 PastixNoTrans, PastixUnit,
                                 nrhs, b, ldb );
            dump_rhs( "AfterDown", bcsc->gN, b );

            /* Solve y = D z with z = (L^t P x) */
            pastix_subtask_diag( pastix_data, bcsc->flttype, nrhs, b, ldb );
            dump_rhs( "AfterDiag", bcsc->gN, b );

            /* Solve z = L^t (P x) */
            pastix_subtask_trsm( pastix_data, bcsc->flttype,
                                 PastixLeft, PastixLower,
                                 trans, PastixUnit,
                                 nrhs, b, ldb );
            dump_rhs( "AfterUp", bcsc->gN, b );
            break;

        case PastixFactLU:
        default:
            /* Solve L y = P b with y = U P x */
            pastix_subtask_trsm( pastix_data, bcsc->flttype,
                                 PastixLeft, PastixLower,
                                 PastixNoTrans, PastixUnit,
                                 nrhs, b, ldb );

            /* Solve y = U (P x) */
            pastix_subtask_trsm( pastix_data, bcsc->flttype,
                                 PastixLeft, PastixUpper,
                                 PastixNoTrans, PastixNonUnit,
                                 nrhs, b, ldb );
            break;
        }
        clockStop(timer);

        dump_rhs( "Final", bcsc->gN, b );

        pastix_data->dparm[DPARM_SOLV_TIME] = clockVal(timer);
        if (iparm[IPARM_VERBOSE] > PastixVerboseNot) {
            pastix_print( 0, 0, OUT_TIME_SOLV,
                          pastix_data->dparm[DPARM_SOLV_TIME] );
        }
    }

    return EXIT_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_users
 *
 * @brief Solve the given problem.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[in] nrhs
 *          The number of right-and-side vectors.
 *
 * @param[inout] b
 *          The right-and-side vectors (can be multiple RHS).
 *          On exit, the solution is stored in place of the right-hand-side vector.
 *
 * @param[in] ldb
 *          The leading dimension of the right-and-side vectors.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastix_task_solve( pastix_data_t *pastix_data,
                   pastix_int_t nrhs, void *b, pastix_int_t ldb )
{
    pastix_bcsc_t *bcsc;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_task_solve: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }

    bcsc  = pastix_data->bcsc;

    /* Compute P * b */
    pastix_subtask_applyorder( pastix_data, bcsc->flttype,
                               PastixDirForward, bcsc->gN, nrhs, b, ldb );

    /* Solve A x = b */
    pastix_subtask_solve( pastix_data, nrhs, b, ldb );

    /* Compute P^t * b */
    pastix_subtask_applyorder( pastix_data, bcsc->flttype,
                               PastixDirBackward, bcsc->gN, nrhs, b, ldb );

    return EXIT_SUCCESS;
}
