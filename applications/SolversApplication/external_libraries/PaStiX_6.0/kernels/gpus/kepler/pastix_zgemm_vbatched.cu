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
#include <assert.h>
#include <cublas.h>
#include "pastix.h"
#include "pastix/datatypes.h"
#include "kernels/pastix_cuda.h"
//#include "kernels/pastix_zcores.h"

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
/* pastix_zgemm_vbatched_nt( */
/*     pastix_trans_t transB, */
/*     pastix_int_t m[32], pastix_int_t n, pastix_int_t k, */
/*     cuDoubleComplex alpha, */
/*     cuDoubleComplex const * dA, pastix_int_t ldda, */
/*     cuDoubleComplex const * dB, pastix_int_t lddb, */
/*     cuDoubleComplex beta, */
/*     cuDoubleComplex       * dC_array[32], pastix_int_t lddc, */
/*     pastix_int_t max_m, pastix_int_t batchCount, const pastix_int_t Acoefind[32], cudaStream_t stream ) */
pastix_zgemm_vbatched_nt(
    pastix_trans_t transB,
    pastix_int_t n, pastix_int_t k,
    cuDoubleComplex alpha,
    const cuDoubleComplex * dB, pastix_int_t lddb,
    cuDoubleComplex beta,
    pastix_int_t max_m, pastix_int_t batchCount, cudaStream_t stream,
    gemm_params_t params )
{
    assert( transB != PastixNoTrans );

    size_t offsetA = 0;
    size_t offsetB = 0;
    offsetA = offsetA / sizeof(cuDoubleComplex);
    offsetB = offsetB / sizeof(cuDoubleComplex);

    if (transB == PastixTrans)
    {
        if(k <= 8)
        {
            // version 58
            pastix_gemm_template_vbatched_nt<cuDoubleComplex, version(NT,58), 0, 0>
                //(m, n, k, dA, ldda, dB, lddb, dC_array, lddc, alpha, beta, batchCount, Acoefind, stream, max_m);
                (n, k, dB, lddb, alpha, beta, batchCount, stream, max_m, params);
        }
        else
        {
            // version 29
            pastix_gemm_template_vbatched_nt<cuDoubleComplex, version(NT,29), 0, 0>
                (n, k, dB, lddb, alpha, beta, batchCount, stream, max_m, params);
        }
    }
    else if (transB == PastixConjTrans)
    {
        if(k <= 8)
        {
            // version 58
            pastix_gemm_template_vbatched_nt<cuDoubleComplex, version(NT,58), 0, 1>
                (n, k, dB, lddb, alpha, beta, batchCount, stream, max_m, params);
        }
        else
        {
            // version 29
            pastix_gemm_template_vbatched_nt<cuDoubleComplex, version(NT,29), 0, 1>
                (n, k, dB, lddb, alpha, beta, batchCount, stream, max_m, params);
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
