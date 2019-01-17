/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @precisions normal z -> s d c z

       @author Jakub Kurzak
       @author Stan Tomov
       @author Mark Gates
       @author Azzam Haidar
       @author Ahmad Abdelfattah

*/

//#include "common.h"
#include <stdio.h>
#include <cublas_v2.h>
#include "pastix/api.h"
#include "pastix/datatypes.h"

#define PastixComplex64_PRECISION

#include "gemm_template_kernel_vbatched.cuh"

#include "gemm_config/zgemm_param_nn.h"
#include "gemm_config/zgemm_param_nt.h"
#include "gemm_config/zgemm_param_tn.h"
#include "gemm_config/zgemm_param_tt.h"

#define version(s,v) s ## _V_ ## v
/**
    Purpose
    -------
    ZGEMM performs one of the matrix-matrix operations

        C = alpha*op( A )*op( B ) + beta*C,

    where op( X ) is one of

        op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,

    alpha and beta are scalars, and A, B and C are matrices, with
    op( A ) an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.

    Parameters
    ----------
    @param[in]
    transA  CHARACTER*1.
            On entry, transA specifies the form of op( A ) to be used in
            the matrix multiplication as follows:
      -     = 'N':  op( A ) = A.
      -     = 'T':  op( A ) = A**T.
      -     = 'C':  op( A ) = A**H.

    @param[in]
    transB  CHARACTER*1.
            On entry, transB specifies the form of op( B ) to be used in
            the matrix multiplication as follows:
      -     = 'N':  op( B ) = B.
      -     = 'T':  op( B ) = B**T.
      -     = 'C':  op( B ) = B**H.

    @param[in]
    m       INTEGER.
            On entry,  M  specifies  the number  of rows  of the  matrix
            op( dA )  and of the  matrix dC.  M  must  be at least  zero.

    @param[in]
    n       INTEGER.
            On entry,  N  specifies the number  of columns of the matrix
            op( dB ) and the number of columns of the matrix dC. N must be
            at least zero.

    @param[in]
    k       INTEGER.
            On entry,  K  specifies  the number of columns of the matrix
            op( dA ) and the number of rows of the matrix op( dB ). K must
            be at least  zero.

    @param[in]
    alpha   COMPLEX_16
            On entry, ALPHA specifies the scalar alpha.

    @param[in]
    dA      COMPLEX_16 array of DIMENSION ( LDA, ka ), where ka is
            k  when  transA = PastixNoTrans,  and is  m  otherwise.
            Before entry with  transA = PastixNoTrans,  the leading  m by k
            part of the array dA must contain the matrix dA, otherwise
            the leading  k by m  part of the array dA must contain  the
            matrix dA.

    @param[in]
    ldda    INTEGER.
            On entry, LDA specifies the first dimension of A as declared
            in the calling (sub) program. When  transA = PastixNoTrans then
            LDA must be at least  max( 1, m ), otherwise  LDA must be at
            least  max( 1, k ).

    @param[in]
    dB      COMPLEX_16 array of DIMENSION ( LDB, kb ), where kb is
            n  when  transB = PastixNoTrans,  and is  k  otherwise.
            Before entry with  transB = PastixNoTrans,  the leading  k by n
            part of the array dB must contain the matrix dB, otherwise
            the leading  n by k  part of the array dB must contain  the
            matrix dB.

    @param[in]
    lddb    INTEGER.
            On entry, LDB specifies the first dimension of dB as declared
            in the calling (sub) program. When  transB = PastixNoTrans then
            LDB must be at least  max( 1, k ), otherwise  LDB must be at
            least  max( 1, n ).

    @param[in]
    beta    COMPLEX_16.
            On entry,  BETA  specifies the scalar  beta.  When  BETA  is
            supplied as zero then dC need not be set on input.

    @param[inout]
    dC      COMPLEX_16 array of DIMENSION ( LDC, n ).
            Before entry, the leading  m by n  part of the array  dC must
            contain the matrix  dC,  except when  beta  is zero, in which
            case dC need not be set on entry.
            On exit, the array  dC  is overwritten by the  m by n  matrix
            ( alpha*op( dA )*op( dB ) + beta*dC ).

    @param[in]
    lddc    INTEGER.
            On entry, LDC specifies the first dimension of dC as declared
            in  the  calling  (sub)  program.   LDC  must  be  at  least
            max( 1, m ).

    @ingroup pastix_zblas3
    ********************************************************************/
#ifdef __cplusplus
extern "C"
#endif
void
pastix_zgemm_vbatched(
    pastix_trans_t transA, pastix_trans_t transB,
    pastix_int_t* m, pastix_int_t* n, pastix_int_t* k,
    cuDoubleComplex alpha,
    cuDoubleComplex const * const * dA_array, pastix_int_t* ldda,
    cuDoubleComplex const * const * dB_array, pastix_int_t* lddb,
    cuDoubleComplex beta,
    cuDoubleComplex **dC_array, pastix_int_t* lddc,
    pastix_int_t max_m, pastix_int_t max_n, pastix_int_t max_k,
    pastix_int_t batchCount, cudaStream_t stream )
{
    pastix_int_t info = 0;
    // TODO: the checker needs modification for variable size: for example, checking lda, and checking m,n,k
    // Shall we do it on CPU or on GPU?
    /*if      ( transA != PastixNoTrans && transA != PastixTrans && transA != PastixConjTrans )
        info = -1;
    else if ( transB != PastixNoTrans && transB != PastixTrans && transB != PastixConjTrans )
        info = -2;
    else if ( m < 0 )
        info = -3;
    else if ( n < 0 )
        info = -4;
    else if ( k < 0 )
        info = -5;
    else if ( transA == PastixNoTrans ? ldda < m : ldda < k )
        info = -8;
    else if ( transB == PastixNoTrans ? lddb < k : lddb < n )
        info = -10;
    else if ( lddc < m )
        info = -13;
    */
    if (info != 0) {
        //pastix_xerbla( __func__, -(info) );
        return;  //info;
    }
    /*
    pastix_int_t arch = pastix_getdevice_arch();
    if ( arch < 200  ) {
        printf("arch < 200 not supported \n"); // TODO call cublas
        return;
    }

    if ( m <= 0 || n <= 0 || k <= 0 )
        return;
    */

    pastix_int_t shape;
    if      (transA == PastixNoTrans   && transB == PastixNoTrans)   {shape = 0;} // nn
    else if (transA == PastixNoTrans   && transB == PastixTrans)     {shape = 1;} // nt
    else if (transA == PastixNoTrans   && transB == PastixConjTrans) {shape = 2;} // nc
    else if (transA == PastixTrans     && transB == PastixNoTrans)   {shape = 3;} // tn
    else if (transA == PastixTrans     && transB == PastixTrans)     {shape = 4;} // tt
    else if (transA == PastixTrans     && transB == PastixConjTrans) {shape = 5;} // tc
    else if (transA == PastixConjTrans && transB == PastixNoTrans)   {shape = 6;} // cn
    else if (transA == PastixConjTrans && transB == PastixTrans)     {shape = 7;} // ct
    else if (transA == PastixConjTrans && transB == PastixConjTrans) {shape = 8;} // cc

    //TODO: probably the texture init code should be placed here

    size_t offsetA = 0;
    size_t offsetB = 0;
    offsetA = offsetA/sizeof(cuDoubleComplex);
    offsetB = offsetB/sizeof(cuDoubleComplex);

    switch(shape)
    {
        case 0: // nn
            {
                gemm_template_vbatched_nn<cuDoubleComplex, version(NN,18), 0, 0>
                (m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, batchCount, stream, max_m, max_n);
            }
            break;
        case 1: // nt
            {
                if(max_k <= 8)
                {
                    // version 58
                    gemm_template_vbatched_nt<cuDoubleComplex, version(NT,58), 0, 0>
                    (m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, batchCount, stream, max_m, max_n);
                }
                else
                {
                    // version 29
                    gemm_template_vbatched_nt<cuDoubleComplex, version(NT,29), 0, 0>
                    (m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, batchCount, stream, max_m, max_n);
                }
            }
            break;
        case 2: // nc
            {
                if(max_k <= 8)
                {
                    // version 58
                    gemm_template_vbatched_nt<cuDoubleComplex, version(NT,58), 0, 1>
                    (m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, batchCount, stream, max_m, max_n);
                }
                else
                {
                    // version 29
                    gemm_template_vbatched_nt<cuDoubleComplex, version(NT,29), 0, 1>
                    (m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, batchCount, stream, max_m, max_n);
                }
            }
            break;
        case 3: // tn
            {
                // version 72
                gemm_template_vbatched_tn<cuDoubleComplex, version(TN,72), 0, 0>
                (m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, batchCount, stream, max_m, max_n);
            }
            break;
        case 6: // cn
            {
                // version 72
                gemm_template_vbatched_tn<cuDoubleComplex, version(TN,72), 1, 0>
                (m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, batchCount, stream, max_m, max_n);
            }
            break;
        case 4: // tt
            {
                // version 13
                gemm_template_vbatched_tt<cuDoubleComplex, version(TT,13), 0, 0>
                (m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, batchCount, stream, max_m, max_n);
            }
            break;
        case 5: // tc
            {
                // version 13
                gemm_template_vbatched_tt<cuDoubleComplex, version(TT,13), 0, 1>
                (m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, batchCount, stream, max_m, max_n);
            }
            break;
        case 7: // ct
            {
                // version 13
                gemm_template_vbatched_tt<cuDoubleComplex, version(TT,13), 1, 0>
                (m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, batchCount, stream, max_m, max_n);
            }
            break;
        case 8: // cc
            {
                // version 13
                gemm_template_vbatched_tt<cuDoubleComplex, version(TT,13), 1, 1>
                (m, n, k, dA_array, ldda, dB_array, lddb, dC_array, lddc, alpha, beta, batchCount, stream, max_m, max_n);
            }
            break;
        default:; // propose something
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
