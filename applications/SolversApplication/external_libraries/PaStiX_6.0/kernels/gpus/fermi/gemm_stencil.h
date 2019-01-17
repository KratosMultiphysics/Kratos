/*
    -- MAGMA (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011
*/

#ifndef _gemm_stencil_h_
#define _gemm_stencil_h_

///////////////////////////////////////////////////////////////////////////////////////////////////
// Common parameters
///////////////////////////////////////////////////////////////////////////////////////////////////

//#define TEXTURE_1D

///////////////////////////////////////////////////////////////////////////////////////////////////

#define trans_nn 1
#define trans_nt 2
#define trans_nc 3

#define trans_tn 4
#define trans_tt 5
#define trans_tc 6

#define trans_cn 7
#define trans_ct 8
#define trans_cc 9

///////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(PastixComplex64_PRECISION)
    typedef cuDoubleComplex FloatingPoint_t;
#define GENERATE_SM_VERSION_NAME_I(func, version) fermi_z##func##_SM##version
#elif defined(PastixComplex32_PRECISION)
    typedef cuFloatComplex FloatingPoint_t;
#define GENERATE_SM_VERSION_NAME_I(func, version) fermi_c##func##_SM##version
#elif defined(PastixDouble_PRECISION)
    typedef double FloatingPoint_t;
#define GENERATE_SM_VERSION_NAME_I(func, version) fermi_d##func##_SM##version
#elif defined(PastixFloat_PRECISION)
    typedef float FloatingPoint_t;
#define GENERATE_SM_VERSION_NAME_I(func, version) fermi_s##func##_SM##version
#endif
#define GENERATE_SM_VERSION_KERNEL_NAME_I2(func, version) GENERATE_SM_VERSION_NAME_I(func, version)
#define GENERATE_SM_VERSION_KERNEL_NAME(func) GENERATE_SM_VERSION_NAME_I2(func, CUDA_SM_VERSION)

#define GENERATE_SM_VERSION_NAME_I2(func, version) GENERATE_SM_VERSION_NAME_I(func, version)
#define GENERATE_SM_VERSION_NAME(func) GENERATE_SM_VERSION_NAME_I2(func, CUDA_SM_VERSION)

///////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef TEXTURE_1D
#define fetch(A, m, n) tex_fetch(tex_ref_##A, coord_##A + n*LD##A+m)
#else
#define fetch(A, m, n) offs_d##A[n*LD##A+m]
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
#if defined(PastixComplex64_PRECISION)

#define conj(A)          cuConj(A)
#define add(A, B)        cuCadd(A, B)
#define mul(A, B)        cuCmul(A, B)
#define fma(A, B, C) C = cuCfma(A, B, C)
#define make_FloatingPoint(x, y) make_cuDoubleComplex(x, y);

#if defined(TEXTURE_1D)

static __device__
FloatingPoint_t tex_fetch(texture<int4> tex_ref, int coord)
{
    int4 v = tex1Dfetch(tex_ref, coord);
    return make_cuDoubleComplex(__hiloint2double(v.y, v.x), __hiloint2double(v.w, v.z));
}

texture<int4, 1, cudaReadModeElementType> tex_ref_A;
texture<int4, 1, cudaReadModeElementType> tex_ref_B;

#endif/* defined(TEXTURE_1D) */

///////////////////////////////////////////////////////////////////////////////////////////////////
#elif defined(PastixComplex32_PRECISION)

#define conj(A)          cuConjf(A)
#define add(A, B)        cuCaddf(A, B)
#define mul(A, B)        cuCmulf(A, B)
#define fma(A, B, C) C = cuCfmaf(A, B, C)
#define make_FloatingPoint(x, y) make_cuFloatComplex(x, y);

#if defined(TEXTURE_1D)

static __device__
FloatingPoint_t tex_fetch(texture<float2> tex_ref, int coord)
{
    return tex1Dfetch(tex_ref, coord);
}

texture<float2, 1, cudaReadModeElementType> tex_ref_A;
texture<float2, 1, cudaReadModeElementType> tex_ref_B;

#endif/* defined(TEXTURE_1D) */

///////////////////////////////////////////////////////////////////////////////////////////////////
#elif defined(PastixDouble_PRECISION)

#define conj(A)           (A)
#define add(A, B)         (A+B)
#define mul(A, B)         (A*B)
#define fma(A, B, C) C += (A*B)
#define make_FloatingPoint(x, y) (x)

#if defined(TEXTURE_1D)

static __device__
FloatingPoint_t tex_fetch(texture<int2> tex_ref, int coord)
{
    int2 v = tex1Dfetch(tex_ref, coord);
    return __hiloint2double(v.y, v.x);
}

texture<int2, 1, cudaReadModeElementType> tex_ref_A;
texture<int2, 1, cudaReadModeElementType> tex_ref_B;

#endif/* defined(TEXTURE_1D) */

///////////////////////////////////////////////////////////////////////////////////////////////////
#elif defined(PastixFloat_PRECISION)

#define conj(A)           (A)
#define add(A, B)         (A+B)
#define mul(A, B)         (A*B)
#define fma(A, B, C) C += (A*B)
#define make_FloatingPoint(x, y) (x)

#if defined(TEXTURE_1D)

static __device__
FloatingPoint_t tex_fetch(texture<float> tex_ref, int coord)
{
    return tex1Dfetch(tex_ref, coord);
}

texture<float, 1, cudaReadModeElementType> tex_ref_A;
texture<float, 1, cudaReadModeElementType> tex_ref_B;

#endif/* defined(TEXTURE_1D) */

///////////////////////////////////////////////////////////////////////////////////////////////////
#endif /* defined(PRECISION_x) */


///////////////////////////////////////////////////////////////////////////////////////////////////
//  Block sizes parameters
///////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(PastixComplex64_PRECISION)

#define DIM_X  8
#define DIM_Y  8

#define BLK_M_nn 24
#define BLK_N_nn 16
#define BLK_K_nn  8

#define DIM_XA_nn 8
#define DIM_YA_nn 8

#define DIM_XB_nn 8
#define DIM_YB_nn 8

///////////////////
#define BLK_M_nt 16
#define BLK_N_nt 24
#define BLK_K_nt  8

#define DIM_XA_nt 8
#define DIM_YA_nt 8

#define DIM_XB_nt 8
#define DIM_YB_nt 8

///////////////////
#define BLK_M_tt 16
#define BLK_N_tt 24
#define BLK_K_tt  8

#define DIM_XA_tt  4
#define DIM_YA_tt 16

#define DIM_XB_tt 8
#define DIM_YB_tt 8

///////////////////
#define BLK_M_tn 24
#define BLK_N_tn 16
#define BLK_K_tn  8

#define DIM_XA_tn 8
#define DIM_YA_tn 8

#define DIM_XB_tn 8
#define DIM_YB_tn 8

////////////////////////////////////////////////////////////////////////////
#elif defined(PastixComplex32_PRECISION)

#define DIM_X 16
#define DIM_Y 16

//////// NN ///////
#define BLK_M_nn 64
#define BLK_N_nn 64
#define BLK_K_nn 16

#define DIM_XA_nn 32
#define DIM_YA_nn  8

#define DIM_XB_nn 16
#define DIM_YB_nn 16

//////// NT ///////
#define BLK_M_nt 64
#define BLK_N_nt 64
#define BLK_K_nt 16

#define DIM_XA_nt 16
#define DIM_YA_nt 16

#define DIM_XB_nt 16
#define DIM_YB_nt 16

//////// TT ///////
#define BLK_M_tt 64
#define BLK_N_tt 64
#define BLK_K_tt 16

#define DIM_XA_tt 16
#define DIM_YA_tt 16

#define DIM_XB_tt 32
#define DIM_YB_tt  8

//////// TN ///////
#define BLK_M_tn 64
#define BLK_N_tn 64
#define BLK_K_tn 16

#define DIM_XA_tn 16
#define DIM_YA_tn 16

#define DIM_XB_tn 16
#define DIM_YB_tn 16

////////////////////////////////////////////////////////////////////////////
#elif defined(PastixDouble_PRECISION)

#define DIM_X 16
#define DIM_Y 16

//////// NN ///////
#define BLK_M_nn 64
#define BLK_N_nn 64
#define BLK_K_nn 16

#define DIM_XA_nn 16
#define DIM_YA_nn 16

#define DIM_XB_nn 16
#define DIM_YB_nn 16

//////// NT ///////
#define BLK_M_nt 64
#define BLK_N_nt 64
#define BLK_K_nt 16

#define DIM_XA_nt 16
#define DIM_YA_nt 16

#define DIM_XB_nt 16
#define DIM_YB_nt 16

//////// TT ///////
#define BLK_M_tt 64
#define BLK_N_tt 64
#define BLK_K_tt 16

#define DIM_XA_tt 16
#define DIM_YA_tt 16

#define DIM_XB_tt 16
#define DIM_YB_tt 16

//////// TN ///////
#define BLK_M_tn 64
#define BLK_N_tn 64
#define BLK_K_tn 16

#define DIM_XA_tn 16
#define DIM_YA_tn 16

#define DIM_XB_tn 16
#define DIM_YB_tn 16

////////////////////////////////////////////////////////////////////////////
#elif defined(PastixFloat_PRECISION)

#define DIM_X 16
#define DIM_Y 16

//////// NN ///////
#define BLK_M_nn 96
#define BLK_N_nn 96
#define BLK_K_nn 16

#define DIM_XA_nn 32
#define DIM_YA_nn  8

#define DIM_XB_nn  8
#define DIM_YB_nn 32

//////// NT ///////
#define BLK_M_nt 96
#define BLK_N_nt 96
#define BLK_K_nt 16

#define DIM_XA_nt 32
#define DIM_YA_nt  8

#define DIM_XB_nt 32
#define DIM_YB_nt  8

//////// TT ///////
#define BLK_M_tt 96
#define BLK_N_tt 96
#define BLK_K_tt 16

#define DIM_XA_tt 16
#define DIM_YA_tt 16

#define DIM_XB_tt 32
#define DIM_YB_tt  8

//////// TN ///////
#define BLK_M_tn 96
#define BLK_N_tn 96
#define BLK_K_tn 16

#define DIM_XA_tn 16
#define DIM_YA_tn 16

#define DIM_XB_tn 16
#define DIM_YB_tn 16

#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  NoTrans - NoTrans
//

#define BLK_M  BLK_M_nn
#define BLK_N  BLK_N_nn
#define BLK_K  BLK_K_nn
#define DIM_XA DIM_XA_nn
#define DIM_YA DIM_YA_nn
#define DIM_XB DIM_XB_nn
#define DIM_YB DIM_YB_nn

#define version trans_nn

#include "gemm_stencil_generic.cu"
#define KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#define KERNEL_RIGHT
#include "gemm_stencil_generic.cu"
#undef KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#undef KERNEL_RIGHT

#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  NoTrans - Trans
//

#define BLK_M  BLK_M_nt
#define BLK_N  BLK_N_nt
#define BLK_K  BLK_K_nt
#define DIM_XA DIM_XA_nt
#define DIM_YA DIM_YA_nt
#define DIM_XB DIM_XB_nt
#define DIM_YB DIM_YB_nt

#define version trans_nt

#include "gemm_stencil_generic.cu"
#define KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#define KERNEL_RIGHT
#include "gemm_stencil_generic.cu"
#undef KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#undef KERNEL_RIGHT

#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

#if defined(PastixComplex64_PRECISION) || defined(PastixComplex32_PRECISION)
#define BLK_M  BLK_M_nt
#define BLK_N  BLK_N_nt
#define BLK_K  BLK_K_nt
#define DIM_XA DIM_XA_nt
#define DIM_YA DIM_YA_nt
#define DIM_XB DIM_XB_nt
#define DIM_YB DIM_YB_nt

#define version trans_nc

#include "gemm_stencil_generic.cu"
#define KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#define KERNEL_RIGHT
#include "gemm_stencil_generic.cu"
#undef KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#undef KERNEL_RIGHT
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Trans - Trans
//

#define BLK_M  BLK_M_tt
#define BLK_N  BLK_N_tt
#define BLK_K  BLK_K_tt
#define DIM_XA DIM_XA_tt
#define DIM_YA DIM_YA_tt
#define DIM_XB DIM_XB_tt
#define DIM_YB DIM_YB_tt

#define version trans_tt
#include "gemm_stencil_generic.cu"
#define KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#define KERNEL_RIGHT
#include "gemm_stencil_generic.cu"
#undef KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#undef KERNEL_RIGHT
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

#if defined(PastixComplex64_PRECISION) || defined(PastixComplex32_PRECISION)
#define BLK_M  BLK_M_tt
#define BLK_N  BLK_N_tt
#define BLK_K  BLK_K_tt
#define DIM_XA DIM_XA_tt
#define DIM_YA DIM_YA_tt
#define DIM_XB DIM_XB_tt
#define DIM_YB DIM_YB_tt

#define version trans_tc

#include "gemm_stencil_generic.cu"
#define KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#define KERNEL_RIGHT
#include "gemm_stencil_generic.cu"
#undef KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#undef KERNEL_RIGHT
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

#define BLK_M  BLK_M_tt
#define BLK_N  BLK_N_tt
#define BLK_K  BLK_K_tt
#define DIM_XA DIM_XA_tt
#define DIM_YA DIM_YA_tt
#define DIM_XB DIM_XB_tt
#define DIM_YB DIM_YB_tt

#define version trans_ct

#include "gemm_stencil_generic.cu"
#define KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#define KERNEL_RIGHT
#include "gemm_stencil_generic.cu"
#undef KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#undef KERNEL_RIGHT
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

#define BLK_M  BLK_M_tt
#define BLK_N  BLK_N_tt
#define BLK_K  BLK_K_tt
#define DIM_XA DIM_XA_tt
#define DIM_YA DIM_YA_tt
#define DIM_XB DIM_XB_tt
#define DIM_YB DIM_YB_tt

#define version trans_cc

#include "gemm_stencil_generic.cu"
#define KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#define KERNEL_RIGHT
#include "gemm_stencil_generic.cu"
#undef KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#undef KERNEL_RIGHT
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Trans - NoTrans
//

#define BLK_M  BLK_M_tn
#define BLK_N  BLK_N_tn
#define BLK_K  BLK_K_tn
#define DIM_XA DIM_XA_tn
#define DIM_YA DIM_YA_tn
#define DIM_XB DIM_XB_tn
#define DIM_YB DIM_YB_tn

#define version trans_tn

#include "gemm_stencil_generic.cu"
#define KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#define KERNEL_RIGHT
#include "gemm_stencil_generic.cu"
#undef KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#undef KERNEL_RIGHT
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version

#if defined(PastixComplex64_PRECISION) || defined(PastixComplex32_PRECISION)
#define BLK_M  BLK_M_tn
#define BLK_N  BLK_N_tn
#define BLK_K  BLK_K_tn
#define DIM_XA DIM_XA_tn
#define DIM_YA DIM_YA_tn
#define DIM_XB DIM_XB_tn
#define DIM_YB DIM_YB_tn

#define version trans_cn

#include "gemm_stencil_generic.cu"
#define KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#define KERNEL_RIGHT
#include "gemm_stencil_generic.cu"
#undef KERNEL_BOTTOM
#include "gemm_stencil_generic.cu"
#undef KERNEL_RIGHT
#undef BLK_M
#undef BLK_N
#undef BLK_K
#undef DIM_XA
#undef DIM_YA
#undef DIM_XB
#undef DIM_YB
#undef version
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////

#endif /* _gemm_stencil_h_ */
