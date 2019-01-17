/*
    -- MAGMA (version 2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @author Jakub Kurzak
       @author Stan Tomov
       @author Mark Gates
       @author Ahmad Abdelfattah
       @author Azzam Haidar

*/

#ifndef GEMM_TEMPLATE_DEVICE_DEFS_H
#define GEMM_TEMPLATE_DEVICE_DEFS_H

#define mymin(a, b) (((a) < (b)) ? (a) : (b))
    
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef TEXTURE_1D

static __device__
FloatingPoint_t tex_fetch(texture<int4> tex_ref, int coord)
{
    #if (__CUDA_ARCH__ >= 200)
    int4 v = tex1Dfetch(tex_ref, coord);
    return make_cuDoubleComplex(__hiloint2double(v.y, v.x), __hiloint2double(v.w, v.z));
    #else
    return make_cuDoubleComplex( 0., 0. );  // dummy code for 1.x compile
    #endif
}

static __device__
FloatingPoint_t tex_fetch(texture<float2> tex_ref, int coord)
{
    return tex1Dfetch(tex_ref, coord);
}

static __device__
FloatingPoint_t tex_fetch(texture<int2> tex_ref, int coord)
{
    #if (__CUDA_ARCH__ >= 200)
    int2 v = tex1Dfetch(tex_ref, coord);
    return __hiloint2double(v.y, v.x);
    #else
    return 0.;  // dummy code for 1.x compile
    #endif
}

static __device__
FloatingPoint_t tex_fetch(texture<float> tex_ref, int coord)
{
    return tex1Dfetch(tex_ref, coord);
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef TEXTURE_1D
    #define fetch(A, m, n, bound) tex_fetch(tex_ref_##A, coord_##A + n*LD##A+m)
#else
    #define fetch(A, m, n, bound) offs_d##A[mymin(n*LD##A+m, bound)]
#endif


#if defined(PastixComplex64_PRECISION)
    #define conj(A)          cuConj(A)
    #define add(A, B)        cuCadd(A, B)
    #define mul(A, B)        cuCmul(A, B)
    #define fma(A, B, C) C = cuCfma(A, B, C)
#define make_FloatingPoint(x, y) make_cuDoubleComplex(x, y);
#elif defined(PastixComplex32_PRECISION)
    #define conj(A)          cuConjf(A)
    #define add(A, B)        cuCaddf(A, B)
    #define mul(A, B)        cuCmulf(A, B)
    #define fma(A, B, C) C = cuCfmaf(A, B, C)
    #define make_FloatingPoint(x, y) make_cuFloatComplex(x, y);
#else
    #define conj(A)           (A)
    #define add(A, B)         (A+B)
    #define mul(A, B)         (A*B)
    #define fma(A, B, C) C += (A*B)
    #define make_FloatingPoint(x, y) (x)
#endif


///////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef TEXTURE_1D

    #if defined(PastixComplex64_PRECISION)
	texture<int4, cudaTextureType1D, cudaReadModeElementType> tex_ref_A;
	texture<int4, cudaTextureType1D, cudaReadModeElementType> tex_ref_B;
    #elif defined(PastixComplex32_PRECISION)
	texture<float2, cudaTextureType1D, cudaReadModeElementType> tex_ref_A;
	texture<float2, cudaTextureType1D, cudaReadModeElementType> tex_ref_B;
    #elif defined(PastixDouble_PRECISION)
	texture<int2, cudaTextureType1D, cudaReadModeElementType> tex_ref_A;
	texture<int2, cudaTextureType1D, cudaReadModeElementType> tex_ref_B;
    #elif defined(PastixFloat_PRECISION)
	texture<float, cudaTextureType1D, cudaReadModeElementType> tex_ref_A;
	texture<float, cudaTextureType1D, cudaReadModeElementType> tex_ref_B;
    #endif

#endif

///////////////////////////////////////////////////////////////////////////////////////////////////

#endif //GEMM_TEMPLATE_DEVICE_DEFS_H
