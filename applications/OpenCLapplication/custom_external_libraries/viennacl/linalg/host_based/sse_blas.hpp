#ifndef VIENNACL_LINALG_HOST_BASED_SSE_BLAS_HPP_
#define VIENNACL_LINALG_HOST_BASED_SSE_BLAS_HPP_

/* =========================================================================
   Copyright (c) 2010-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennacl/linalg/host_based/sse_blas.hpp
*   @brief optimized BLAS functions using SSE2 and SSE3 intrinsic functions
*
*   Contributed by Alex Christensen.
*/

//complex BLAS functions are included but unused in this version of ViennaCL
#if defined VIENNACL_WITH_COMPLEX
#include <complex>
#endif

//defining VIENNACL_SSE3 adds a slight optimization for complex multiplication using SSE3
#if defined VIENNACL_WITH_SSE3
#include <pmmintrin.h>
#elif defined VIENNACL_WITH_SSE2
#include <emmintrin.h>
#endif

#include <cstddef>
#include <cmath>

namespace viennacl
{
  namespace linalg
  {

    namespace host_based
    {
      //saxpy, daxpy, caxpy, zaxpy
      template <class T> inline void _axpy(const T*, T*, vcl_size_t, T);

      //sdot, ddot, cdotu, zdotu
      template <class T> inline T    _dot (vcl_size_t, const T*, const T*);

      //sdot, ddot, cdotc, zdotc
      template <class T> inline T    _dotc(vcl_size_t, const T*, const T*);

      //sswap, dswap, cswap, zswap
      template <class T> inline void _swap(vcl_size_t, T*, T*);

      //scopy, dcopy, ccopy, zcopy
      template <class T> inline void _copy(vcl_size_t, T*, T*);

      //snrm2, dnrm2, euclidian norm of complex vectors
      template <class T> inline T    _nrm2(const T*, vcl_size_t);

      namespace detail
      {
        template <class T> inline T conjIfComplex(T x){return x;}
      }

      template <class T>
      inline void _axpy(const T* x, T* y, vcl_size_t n, T a)
      {
        for(vcl_size_t i=0;i<n;i++)
          y[i]+=a*x[i];
      }

      template <class T>
      inline T _dot(vcl_size_t n, const T* x, const T* y)
      {
        T sum(0);
        for(vcl_size_t i=0;i<n;i++)
          sum+=x[i]*y[i];
        return sum;
      }

      template <class T>
      inline T _dotc(vcl_size_t n, const T* x, const T* y)
      {
        T sum(0);
        for(vcl_size_t i=0;i<n;i++)
          sum+=detail::conjIfComplex(x[i])*y[i];
        return sum;
      }

      template <class T>
      inline void _swap(vcl_size_t n, T* sx, T* sy)
      {
        T t;
        for(vcl_size_t i=0;i<n;i++)
        {
          t=sx[i];
          sx[i]=sy[i];
          sy[i]=t;
        }
      }

      template <class T>
      inline void _copy(vcl_size_t n, T* cx, T* cy)
      {
        for(vcl_size_t i=0;i<n;i++)
          cx[i]=cy[i];
      }

      template <class T>
      inline T _nrm2(const T* x, vcl_size_t n)
      {
        //based on http://www.netlib.org/blas/snrm2.f, but works with std::complex

        if(n<1)
          return T(0);
        if(n==1)
          return std::abs(x[0]);
        T scale(0);
        T scaledSquareSum(1);
        for(vcl_size_t i=0;i<n;i++){
          if(x[i]!=T(0)){
            T absXi=std::abs(x[i]);
            if(std::abs(x[i])>std::abs(scale)){
              T temp=scale/absXi;
              scaledSquareSum=T(1)+scaledSquareSum*temp*temp;
              scale=absXi;
            }
            else{
              T temp=absXi/scale;
              scaledSquareSum+=temp*temp;
            }
          }
        }
        return scale*sqrt(scaledSquareSum);
      }

  #if defined VIENNACL_WITH_COMPLEX

      namespace detail
      {
        template <> inline std::complex<double> conjIfComplex(std::complex<double> x){return conj(x);}
        template <> inline std::complex<float > conjIfComplex(std::complex<float > x){return conj(x);}
      }

      template <>
      inline std::complex<double> _nrm2(const std::complex<double>* x, vcl_size_t n)
      {
        //based on http://www.netlib.org/blas/snrm2.f

        if(n<1)
          return std::complex<double>(0);
        if(n==1)
          return std::abs(x[0]);
        double scale=0.0;
        double scaledSquareSum=1.0;
        for(vcl_size_t i=0;i<n;i++){
          if(x[i].real()!=0.0){
            double absXi=std::abs(x[i].real());
            if(absXi>scale){
              double temp=scale/absXi;
              scaledSquareSum=1.0+scaledSquareSum*temp*temp;
              scale=absXi;
            }
            else{
              double temp=absXi/scale;
              scaledSquareSum+=temp*temp;
            }
          }
          if(x[i].imag()!=0.0){
            double absXi=std::abs(x[i].imag());
            if(absXi>scale){
              double temp=scale/absXi;
              scaledSquareSum=1.0+scaledSquareSum*temp*temp;
              scale=absXi;
            }
            else{
              double temp=absXi/scale;
              scaledSquareSum+=temp*temp;
            }
          }
        }
        return std::complex<double>(scale*sqrt(scaledSquareSum));
      }

      template <>
      inline std::complex<float> _nrm2(const std::complex<float>* x, vcl_size_t n)
      {
        //based on http://www.netlib.org/blas/snrm2.f

        if(n<1)
          return std::complex<float>(0);
        if(n==1)
          return std::abs(x[0]);
        float scale=0.0;
        float scaledSquareSum=1.0;
        for(vcl_size_t i=0;i<n;i++){
          if(x[i].real()!=0.0){
            float absXi=std::abs(x[i].real());
            if(absXi>scale){
              float temp=scale/absXi;
              scaledSquareSum=1.0f+scaledSquareSum*temp*temp;
              scale=absXi;
            }
            else{
              float temp=absXi/scale;
              scaledSquareSum+=temp*temp;
            }
          }
          if(x[i].imag()!=0.0){
            float absXi=std::abs(x[i].imag());
            if(absXi>scale){
              float temp=scale/absXi;
              scaledSquareSum=1.0f+scaledSquareSum*temp*temp;
              scale=absXi;
            }
            else{
              float temp=absXi/scale;
              scaledSquareSum+=temp*temp;
            }
          }
        }
        return std::complex<float>(scale*sqrt(scaledSquareSum));
      }

  #endif //defined VIENNACL_COMPLEX

  #if defined VIENNACL_WITH_SSE2

      //saxpy
      template <>
      inline void _axpy<float>(const float* x, float* y, vcl_size_t n, float a)
      {

        //if the array is short or if either array is unaligned, perform the non-SSE code
        if(n<16||((vcl_size_t)x)%16!=((vcl_size_t)y)%16||((vcl_size_t)x)%sizeof(float)!=0)
          for(vcl_size_t i=0;i<n;i++)
            y[i]+=a*x[i];
        else
        {
          //process unaligned section of arrays
          while(((vcl_size_t)x)%16)
          {
            if(n<=0)
              return;
            y[0]+=a*x[0];
            x++;
            y++;
            n--;
          }

          __m128 sum0;
          __m128 sum1;
          __m128 reg0,reg1,reg2,reg3;
          __m128 areg=_mm_set1_ps(a);
          __m128 prod;

          //add floats 8 at a time
          while(n>=8){

            //read floats into MMX registers (8 from each array)
            reg0=_mm_load_ps(x+0);
            reg1=_mm_load_ps(x+4);
            reg2=_mm_load_ps(y+0);
            reg3=_mm_load_ps(y+4);

            //add floats
            prod=_mm_mul_ps(reg0,areg);
            sum0=_mm_add_ps(prod,reg2);
            prod=_mm_mul_ps(reg1,areg);
            sum1=_mm_add_ps(prod,reg3);

            //put float sums into y
            _mm_store_ps(y+0,sum0);
            _mm_store_ps(y+4,sum1);

            x+=8;
            y+=8;
            n-=8;
          }

          //add beyond the last multiple of 8
          for(vcl_size_t i=0;i<n;i++)
            y[i]+=a*x[i];
        }
      }

      //daxpy
      template <>
      inline void _axpy<double>(const double* x, double* y, vcl_size_t n, double a)
      {
        //if the array is short or if either array is unaligned, perform the non-SSE code
        if(n<16||((vcl_size_t)x)%16!=((vcl_size_t)y)%16||((vcl_size_t)x)%sizeof(double)!=0)
          for(vcl_size_t i=0;i<n;i++)
            y[i]+=a*x[i];

        else
        {
          //process unaligned section of arrays
          while(((vcl_size_t)x)%16)
          {
            if(n<=0)
              return;
            y[0]+=a*x[0];
            x++;
            y++;
            n--;
          }

          __m128d sum0;
          __m128d sum1;
          __m128d reg0,reg1,reg2,reg3;
          __m128d areg=_mm_set1_pd(a);
          __m128d prod;

          //add doubles 4 at a time
          while(n>=8){

            //read floats into MMX registers (4 from each array)
            reg0=_mm_load_pd(x+0);
            reg1=_mm_load_pd(x+2);
            reg2=_mm_load_pd(y+0);
            reg3=_mm_load_pd(y+2);

            //add floats
            prod=_mm_mul_pd(reg0,areg);
            sum0=_mm_add_pd(prod,reg2);
            prod=_mm_mul_pd(reg1,areg);
            sum1=_mm_add_pd(prod,reg3);

            //put float sums into y
            _mm_store_pd(y+0,sum0);
            _mm_store_pd(y+2,sum1);

            x+=4;
            y+=4;
            n-=4;
          }

          //add beyond the last multiple of 4
          for(vcl_size_t i=0;i<n;i++)
            y[i]+=a*x[i];
        }
      }

      //sdot
      template <>
      inline float _dot<float>(vcl_size_t n, const float* x, const float* y)
      {

        //if the array is short or if either array is unaligned, perform the non-SSE code
        if(n<16||((vcl_size_t)x)%16!=((vcl_size_t)y)%16||((vcl_size_t)x)%sizeof(float)!=0)
        {
          float sum=0;
          for(vcl_size_t i=0;i<n;i++)
            sum+=x[i]*y[i];
          return sum;
        }
        else
        {

          //process unaligned section of array
          float sum=0;
          while(((vcl_size_t)x)%16)
          {
            if(n<=0)
              return sum;
            sum+=x[0]*y[0];
            y++;
            x++;
            n--;
          }

          __m128 sumReg=_mm_setzero_ps();
          __m128 reg0,reg1,reg2,reg3;

          //add floats 8 at a time
          while(n>=8)
          {
            //read floats into MMX registers (8 from each array)
            reg0=_mm_load_ps(x+0);
            reg1=_mm_load_ps(x+4);
            reg2=_mm_load_ps(y+0);
            reg3=_mm_load_ps(y+4);

            //multiply floats together
            reg0=_mm_mul_ps(reg0,reg2);
            reg1=_mm_mul_ps(reg1,reg3);

            //add to sums
            sumReg=_mm_add_ps(sumReg,reg0);
            sumReg=_mm_add_ps(sumReg,reg1);

            x+=8;
            y+=8;
            n-=8;
          }

          //add beyond where the inner loop stopped
          for(vcl_size_t i=0;i<n;i++)
            sum+=x[i]*y[i];

          //move the sums from the xmm registers to aligned memory on the stack
          float sums[8];
          float* pSums=(float*)((((vcl_size_t)sums)&(~15))+16);
          _mm_store_ps(pSums,sumReg);

          return sum+pSums[0]+pSums[1]+pSums[2]+pSums[3];
        }
      }

      //ddot
      template <>
      inline double _dot(vcl_size_t n, const double* x, const double* y)
      {
        //if the array is short or if either array is unaligned, perform the non-SSE code
        if(n<16||((vcl_size_t)x)%16!=((vcl_size_t)y)%16||((vcl_size_t)x)%sizeof(double)!=0)
        {
          double sum=0;
          for(vcl_size_t i=0;i<n;i++)
            sum+=x[i]*y[i];
          return sum;
        }
        else
        {
          //process unaligned section of array
          double sum=0;
          while(((vcl_size_t)x)%16)
          {
            if(n<=0)
              return sum;
            sum+=x[0]*y[0];
            y++;
            x++;
            n--;
          }

          __m128d sum0=_mm_setzero_pd();
          __m128d sum1=_mm_setzero_pd();
          __m128d reg0,reg1,reg2,reg3;

          //add doubles 4 at a time
          while(n>=4)
          {
            //read doubles into MMX registers (4 from each array)
            reg0=_mm_load_pd(x+0);
            reg1=_mm_load_pd(x+2);
            reg2=_mm_load_pd(y+0);
            reg3=_mm_load_pd(y+2);

            //multiply doubles together
            reg0=_mm_mul_pd(reg0,reg2);
            reg1=_mm_mul_pd(reg1,reg3);

            //add to sums
            sum0=_mm_add_pd(sum0,reg0);
            sum1=_mm_add_pd(sum1,reg1);

            x+=4;
            y+=4;
            n-=4;
          }

          //add beyond where the inner loop stopped
          for(vcl_size_t i=0;i<n;i++)
            sum+=x[i]*y[i];

          //move the sums from the xmm registers to aligned memory on the stack
          double sums[4];
          double* pSums=(double*)((((vcl_size_t)sums)&(~15))+16);
          sum0=_mm_add_pd(sum0,sum1);
          _mm_store_pd(pSums,sum0);

          return sum+pSums[0]+pSums[1];
        }
      }

      //conjugated dot products are the same as non-conjugated dot products for real numbers
      template <> inline float  _dotc<float >(vcl_size_t n, const float  *x, const float  *y){return _dot(n,x,y);}
      template <> inline double _dotc<double>(vcl_size_t n, const double *x, const double *y){return _dot(n,x,y);}

  #if defined VIENNACL_WITH_COMPLEX

      //caxpy
      template <>
      inline void _axpy<std::complex<float> >(const std::complex<float>* x, std::complex<float>* y, vcl_size_t n, std::complex<float> a)
      {
        //if the array is short or if either array is unaligned, perform the non-SSE code
        if(n<16||((vcl_size_t)x)%16!=((vcl_size_t)y)%16||((vcl_size_t)x)%sizeof(std::complex<float>)!=0)
          for(vcl_size_t i=0;i<n;i++)
            y[i]+=a*x[i];

        else
        {
          //process unaligned section of arrays
          while(((vcl_size_t)x)%16)
          {
            if(n<=0)
              return;
            y[0]+=a*x[0];
            x++;
            y++;
            n--;
          }

          __m128 reg0,reg1,reg2,reg3,reg4;
          __m128 areg0=_mm_set_ps(a.imag(),a.real(),a.imag(),a.real());
          __m128 areg1=_mm_set_ps(a.real(),a.imag(),a.real(),a.imag());
  #ifndef VIENNACL_WITH_SSE3
          __m128 nreg=_mm_set_ps(1.0f,-1.0f,1.0f,-1.0f);
  #endif

          //add complex floats 4 at a time
          while(n>=4)
          {
            //read floats into MMX registers (8 from each array)
            reg0=_mm_load_ps((float*)(x+0));
            reg1=_mm_load_ps((float*)(x+2));
            reg2=_mm_load_ps((float*)(y+0));
            reg3=_mm_load_ps((float*)(y+2));

            //do complex multiplication and addition
  #ifndef VIENNACL_WITH_SSE3
            reg4=_mm_shuffle_ps(reg0,reg0,0xA0);
            reg0=_mm_shuffle_ps(reg0,reg0,0xF5);
            reg4=_mm_mul_ps(reg4,areg0);
            reg0=_mm_mul_ps(reg0,areg1);
            reg0=_mm_mul_ps(reg0,nreg);
            reg0=_mm_add_ps(reg4,reg0);
            reg0=_mm_add_ps(reg0,reg2);
            reg4=_mm_shuffle_ps(reg1,reg1,0xA0);
            reg1=_mm_shuffle_ps(reg1,reg1,0xF5);
            reg4=_mm_mul_ps(reg4,areg0);
            reg1=_mm_mul_ps(reg1,areg1);
            reg1=_mm_mul_ps(reg1,nreg);
            reg1=_mm_add_ps(reg4,reg1);
            reg1=_mm_add_ps(reg1,reg3);
  #else
            reg4=_mm_moveldup_ps(reg0);
            reg0=_mm_movehdup_ps(reg0);
            reg4=_mm_mul_ps(reg4,areg0);
            reg0=_mm_mul_ps(reg0,areg1);
            reg0=_mm_addsub_ps(reg4,reg0);
            reg0=_mm_add_ps(reg0,reg2);
            reg4=_mm_moveldup_ps(reg1);
            reg1=_mm_movehdup_ps(reg1);
            reg4=_mm_mul_ps(reg4,areg0);
            reg1=_mm_mul_ps(reg1,areg1);
            reg1=_mm_addsub_ps(reg4,reg1);
            reg1=_mm_add_ps(reg1,reg3);
  #endif
            //put results into y
            _mm_store_ps((float*)(y+0),reg0);
            _mm_store_ps((float*)(y+2),reg1);

            x+=4;
            y+=4;
            n-=4;
          }

          //add beyond the last multiple of 4
          for(vcl_size_t i=0;i<n;i++)
            y[i]+=a*x[i];
        }
      }

      //zaxpy
      template <>
      inline void _axpy<std::complex<double> >(const std::complex<double>* x, std::complex<double>* y, vcl_size_t n, std::complex<double> a)
      {
        //if the array is short or if either array is unaligned, perform the non-SSE code
        if(n<16||((vcl_size_t)x)%16||((vcl_size_t)y)%16)
          for(vcl_size_t i=0;i<n;i++)
            y[i]+=a*x[i];

        else
        {
          __m128d reg0,reg1,reg2,reg3,reg4;
          __m128d areg0=_mm_set_pd(a.imag(),a.real());
          __m128d areg1=_mm_set_pd(a.real(),a.imag());
  #ifndef VIENNACL_WITH_SSE3
          __m128d nreg=_mm_set_pd(1.0,-1.0);
  #endif

          //add complex doubles 2 at a time
          while(n>=2)
          {
            //read doubles into MMX registers (4 from each array)
            reg0=_mm_load_pd((double*)(x+0));
            reg1=_mm_load_pd((double*)(x+1));
            reg2=_mm_load_pd((double*)(y+0));
            reg3=_mm_load_pd((double*)(y+1));

            //do complex multiplication and addition
  #ifndef VIENNACL_WITH_SSE3
            reg4=_mm_shuffle_pd(reg0,reg0,0x0);
            reg0=_mm_shuffle_pd(reg0,reg0,0x3);
            reg4=_mm_mul_pd(reg4,areg0);
            reg0=_mm_mul_pd(reg0,areg1);
            reg0=_mm_mul_pd(reg0,nreg);
            reg0=_mm_add_pd(reg4,reg0);
            reg0=_mm_add_pd(reg0,reg2);
            reg4=_mm_shuffle_pd(reg1,reg1,0x0);
            reg1=_mm_shuffle_pd(reg1,reg1,0x3);
            reg4=_mm_mul_pd(reg4,areg0);
            reg1=_mm_mul_pd(reg1,areg1);
            reg1=_mm_mul_pd(reg1,nreg);
            reg1=_mm_add_pd(reg4,reg1);
            reg1=_mm_add_pd(reg1,reg3);
  #else
            reg4=_mm_shuffle_pd(reg0,reg0,0x0);
            reg0=_mm_shuffle_pd(reg0,reg0,0x3);
            reg4=_mm_mul_pd(reg4,areg0);
            reg0=_mm_mul_pd(reg0,areg1);
            reg0=_mm_addsub_pd(reg4,reg0);
            reg0=_mm_add_pd(reg0,reg2);
            reg4=_mm_shuffle_pd(reg1,reg1,0x0);
            reg1=_mm_shuffle_pd(reg1,reg1,0x3);
            reg4=_mm_mul_pd(reg4,areg0);
            reg1=_mm_mul_pd(reg1,areg1);
            reg1=_mm_addsub_pd(reg4,reg1);
            reg1=_mm_add_pd(reg1,reg3);
  #endif
            //put results into y
            _mm_store_pd((double*)(y+0),reg0);
            _mm_store_pd((double*)(y+1),reg1);

            x+=2;
            y+=2;
            n-=2;
          }

          //add beyond the last multiple of 2
          if(n)
            y[0]+=a*x[0];
        }
      }

      //cdotu
      template <>
      inline std::complex<float> _dot<std::complex<float> >(vcl_size_t n, const std::complex<float>* x, const std::complex<float>* y)
      {
        //if the array is short or if either array is unaligned, perform the non-SSE code
        if(n<16||((vcl_size_t)x)%16!=((vcl_size_t)y)%16||((vcl_size_t)x)%sizeof(std::complex<float>)!=0)
        {
          std::complex<float> sum(0);
          for(vcl_size_t i=0;i<n;i++)
            sum+=x[i]*y[i];
          return sum;
        }
        else
        {
          //process unaligned section of arrays
          std::complex<float> sum(0);
          while(((vcl_size_t)x)%16)
          {
            if(n<=0)
              return sum;
            sum+=x[0]*y[0];
            y++;
            x++;
            n--;
          }

          __m128 sumReg=_mm_setzero_ps();
          __m128 reg0,reg1,reg2,reg3,reg4;
  #ifndef VIENNACL_WITH_SSE3
          __m128 nreg=_mm_set_ps(1.0f,-1.0f,1.0f,-1.0f);
  #endif

          //add complex floats 4 at a time
          while(n>=4)
          {
            //read floats into MMX registers (8 from each array)
            reg0=_mm_load_ps((float*)(x+0));
            reg1=_mm_load_ps((float*)(x+2));
            reg2=_mm_load_ps((float*)(y+0));
            reg3=_mm_load_ps((float*)(y+2));

            //multiply complex floats together
  #ifndef VIENNACL_WITH_SSE3
            reg4=_mm_shuffle_ps(reg2,reg2,0xA0);
            reg2=_mm_shuffle_ps(reg2,reg2,0xF5);
            reg4=_mm_mul_ps(reg4,reg0);
            reg2=_mm_mul_ps(reg2,reg0);
            reg2=_mm_shuffle_ps(reg2,reg2,0xB1);
            reg2=_mm_mul_ps(reg2,nreg);
            reg0=_mm_add_ps(reg4,reg2);
            reg4=_mm_shuffle_ps(reg3,reg3,0xA0);
            reg3=_mm_shuffle_ps(reg3,reg3,0xF5);
            reg4=_mm_mul_ps(reg4,reg1);
            reg3=_mm_mul_ps(reg3,reg1);
            reg3=_mm_shuffle_ps(reg3,reg3,0xB1);
            reg3=_mm_mul_ps(reg3,nreg);
            reg1=_mm_add_ps(reg4,reg3);
  #else
            reg4=_mm_moveldup_ps(reg2);
            reg2=_mm_movehdup_ps(reg2);
            reg4=_mm_mul_ps(reg4,reg0);
            reg2=_mm_mul_ps(reg2,reg0);
            reg2=_mm_shuffle_ps(reg2,reg2,0xB1);
            reg0=_mm_addsub_ps(reg4,reg2);
            reg4=_mm_moveldup_ps(reg3);
            reg3=_mm_movehdup_ps(reg3);
            reg4=_mm_mul_ps(reg4,reg1);
            reg3=_mm_mul_ps(reg3,reg1);
            reg3=_mm_shuffle_ps(reg3,reg3,0xB1);
            reg1=_mm_addsub_ps(reg4,reg3);
  #endif

            //add to sum
            sumReg=_mm_add_ps(sumReg,reg0);
            sumReg=_mm_add_ps(sumReg,reg1);

            x+=4;
            y+=4;
            n-=4;
          }

          //add beyond where the inner loop stopped
          for(vcl_size_t i=0;i<n;i++)
            sum+=x[i]*y[i];

          //move the sums from the xmm registers to aligned memory on the stack
          std::complex<float> sums[4];
          std::complex<float>* pSums=(std::complex<float>*)((((vcl_size_t)sums)&(~15))+16);
          pSums[0]=std::complex<float>(0);
          pSums[1]=std::complex<float>(0);
          _mm_store_ps((float*)pSums,sumReg);

          return sum+pSums[0]+pSums[1];
        }
      }

      //zdotu
      template <>
      inline std::complex<double> _dot<std::complex<double> >(vcl_size_t n, const std::complex<double>* x, const std::complex<double>* y)
      {
        //if the array is short or if either array is unaligned, perform the non-SSE code
        if(n<16||((vcl_size_t)x)%16||((vcl_size_t)y)%16)
        {
          std::complex<double> sum(0);
          for(vcl_size_t i=0;i<n;i++)
            sum+=x[i]*y[i];
          return sum;
        }
        else
        {
          __m128d sumReg=_mm_setzero_pd();
          __m128d reg0,reg1,reg2,reg3,reg4;
  #ifndef VIENNACL_WITH_SSE3
          __m128d nreg=_mm_set_pd(1.0,-1.0);
  #endif

          //add complex doubles 2 at a time
          while(n>=2)
          {
            //read doubles into MMX registers (4 from each array)
            reg0=_mm_load_pd((double*)(x+0));
            reg1=_mm_load_pd((double*)(x+1));
            reg2=_mm_load_pd((double*)(y+0));
            reg3=_mm_load_pd((double*)(y+1));

            //multiply complex doubles together
  #ifndef VIENNACL_WITH_SSE3
            reg4=_mm_shuffle_pd(reg2,reg2,0x0);
            reg2=_mm_shuffle_pd(reg2,reg2,0x3);
            reg4=_mm_mul_pd(reg4,reg0);
            reg2=_mm_mul_pd(reg2,reg0);
            reg2=_mm_shuffle_pd(reg2,reg2,0x1);
            reg2=_mm_mul_pd(reg2,nreg);
            reg0=_mm_add_pd(reg4,reg2);
            reg4=_mm_shuffle_pd(reg3,reg3,0x0);
            reg3=_mm_shuffle_pd(reg3,reg3,0x3);
            reg4=_mm_mul_pd(reg4,reg1);
            reg3=_mm_mul_pd(reg3,reg1);
            reg3=_mm_shuffle_pd(reg3,reg3,0x1);
            reg3=_mm_mul_pd(reg3,nreg);
            reg1=_mm_add_pd(reg4,reg3);
  #else
            reg4=_mm_shuffle_pd(reg2,reg2,0x0);
            reg2=_mm_shuffle_pd(reg2,reg2,0x3);
            reg4=_mm_mul_pd(reg4,reg0);
            reg2=_mm_mul_pd(reg2,reg0);
            reg2=_mm_shuffle_pd(reg2,reg2,0x1);
            reg0=_mm_addsub_pd(reg4,reg2);
            reg4=_mm_shuffle_pd(reg3,reg3,0x0);
            reg3=_mm_shuffle_pd(reg3,reg3,0x3);
            reg4=_mm_mul_pd(reg4,reg1);
            reg3=_mm_mul_pd(reg3,reg1);
            reg3=_mm_shuffle_pd(reg3,reg3,0x1);
            reg1=_mm_addsub_pd(reg4,reg3);
  #endif

            //add to sum
            sumReg=_mm_add_pd(sumReg,reg0);
            sumReg=_mm_add_pd(sumReg,reg1);

            x+=2;
            y+=2;
            n-=2;
          }

          //add beyond where the inner loop stopped
          std::complex<double> sum(0);
          if(n)
            sum=x[0]*y[0];

          //move the sums from the xmm registers to aligned memory on the stack
          std::complex<double> sums[2];
          std::complex<double>* pSums=(std::complex<double>*)((((vcl_size_t)sums)&(~15))+16);
          pSums[0]=std::complex<double>(0);
          _mm_store_pd((double*)pSums,sumReg);

          return sum+pSums[0];
        }
      }

      //cdotc
      template <>
      inline std::complex<float> _dotc<std::complex<float> >(vcl_size_t n, const std::complex<float>* x, const std::complex<float>* y)
      {
        //if the array is short or if either array is unaligned, perform the non-SSE code
        if(n<16||((vcl_size_t)x)%16!=((vcl_size_t)y)%16||((vcl_size_t)x)%sizeof(std::complex<float>)!=0)
        {
          std::complex<float> sum(0);
          for(vcl_size_t i=0;i<n;i++)
            sum+=conj(x[i])*y[i];
          return sum;
        }
        else
        {
          //process unaligned section of arrays
          std::complex<float> sum(0);
          while(((vcl_size_t)x)%16)
          {
            if(n<=0)
              return sum;
            sum+=conj(x[0])*y[0];
            y++;
            x++;
            n--;
          }

          __m128 sumReg=_mm_setzero_ps();
          __m128 reg0,reg1,reg2,reg3,reg4;
  #ifndef VIENNACL_WITH_SSE3
          __m128 nreg=_mm_set_ps(1.0f,-1.0f,1.0f,-1.0f);
  #endif

          //add complex floats 4 at a time
          while(n>=4)
          {
            //read floats into MMX registers (8 from each array)
            reg0=_mm_load_ps((float*)(x+0));
            reg1=_mm_load_ps((float*)(x+2));
            reg2=_mm_load_ps((float*)(y+0));
            reg3=_mm_load_ps((float*)(y+2));

            //multiply complex doubles together
  #ifndef VIENNACL_WITH_SSE3
            reg4=_mm_shuffle_ps(reg2,reg2,0xA0);
            reg2=_mm_shuffle_ps(reg2,reg2,0xF5);
            reg4=_mm_mul_ps(reg4,reg0);
            reg2=_mm_mul_ps(reg2,reg0);
            reg4=_mm_shuffle_ps(reg4,reg4,0xB1);
            reg4=_mm_mul_ps(reg4,nreg);
            reg0=_mm_add_ps(reg4,reg2);
            reg4=_mm_shuffle_ps(reg3,reg3,0xA0);
            reg3=_mm_shuffle_ps(reg3,reg3,0xF5);
            reg4=_mm_mul_ps(reg4,reg1);
            reg3=_mm_mul_ps(reg3,reg1);
            reg4=_mm_shuffle_ps(reg4,reg4,0xB1);
            reg4=_mm_mul_ps(reg4,nreg);
            reg1=_mm_add_ps(reg4,reg3);
  #else
            reg4=_mm_moveldup_ps(reg2);
            reg2=_mm_movehdup_ps(reg2);
            reg4=_mm_mul_ps(reg4,reg0);
            reg2=_mm_mul_ps(reg2,reg0);
            reg4=_mm_shuffle_ps(reg4,reg4,0xB1);
            reg0=_mm_addsub_ps(reg2,reg4);
            reg4=_mm_moveldup_ps(reg3);
            reg3=_mm_movehdup_ps(reg3);
            reg4=_mm_mul_ps(reg4,reg1);
            reg3=_mm_mul_ps(reg3,reg1);
            reg4=_mm_shuffle_ps(reg4,reg4,0xB1);
            reg1=_mm_addsub_ps(reg3,reg4);
  #endif

            //add to sum
            sumReg=_mm_add_ps(sumReg,reg0);
            sumReg=_mm_add_ps(sumReg,reg1);

            x+=4;
            y+=4;
            n-=4;
          }

          //add beyond where the inner loop stopped
          for(vcl_size_t i=0;i<n;i++)
            sum+=conj(x[i])*y[i];

          //move the sums from the xmm registers to aligned memory on the stack
          std::complex<float> sums[4];
          std::complex<float>* pSums=(std::complex<float>*)((((vcl_size_t)sums)&(~15))+16);
          sumReg=_mm_shuffle_ps(sumReg,sumReg,0xB1);//swap real and imag
          _mm_store_ps((float*)pSums,sumReg);

          return sum+pSums[0]+pSums[1];
        }
      }

      //zdotc
      template <>
      inline std::complex<double> _dotc<std::complex<double> >(vcl_size_t n, const std::complex<double>* x, const std::complex<double>* y)
      {
        //if the array is short or if either array is unaligned, perform the non-SSE code
        if(n<16||((vcl_size_t)x)%16||((vcl_size_t)y)%16)
        {
          std::complex<double> sum(0);
          for(vcl_size_t i=0;i<n;i++)
            sum+=conj(x[i])*y[i];
          return sum;
        }
        else
        {
          __m128d sumReg=_mm_setzero_pd();
          __m128d reg0,reg1,reg2,reg3,reg4;
  #ifndef VIENNACL_WITH_SSE3
          __m128d nreg=_mm_set_pd(1.0,-1.0);
  #endif

          //add complex doubles 2 at a time
          while(n>=2)
          {
            //read doubles into MMX registers (4 from each array)
            reg0=_mm_load_pd((double*)(x+0));
            reg1=_mm_load_pd((double*)(x+1));
            reg2=_mm_load_pd((double*)(y+0));
            reg3=_mm_load_pd((double*)(y+1));

            //multiply complex floats together
  #ifndef VIENNACL_WITH_SSE3
            reg4=_mm_shuffle_pd(reg2,reg2,0x0);
            reg2=_mm_shuffle_pd(reg2,reg2,0x3);
            reg4=_mm_mul_pd(reg4,reg0);
            reg2=_mm_mul_pd(reg2,reg0);
            reg4=_mm_shuffle_pd(reg4,reg4,0x1);
            reg4=_mm_mul_pd(reg4,nreg);
            reg0=_mm_add_pd(reg4,reg2);
            reg4=_mm_shuffle_pd(reg3,reg3,0x0);
            reg3=_mm_shuffle_pd(reg3,reg3,0x3);
            reg4=_mm_mul_pd(reg4,reg1);
            reg3=_mm_mul_pd(reg3,reg1);
            reg4=_mm_shuffle_pd(reg4,reg4,0x1);
            reg4=_mm_mul_pd(reg4,nreg);
            reg1=_mm_add_pd(reg4,reg3);
  #else
            reg4=_mm_shuffle_pd(reg2,reg2,0x0);
            reg2=_mm_shuffle_pd(reg2,reg2,0x3);
            reg4=_mm_mul_pd(reg4,reg0);
            reg2=_mm_mul_pd(reg2,reg0);
            reg4=_mm_shuffle_pd(reg4,reg4,0x1);
            reg0=_mm_addsub_pd(reg2,reg4);
            reg4=_mm_shuffle_pd(reg3,reg3,0x0);
            reg3=_mm_shuffle_pd(reg3,reg3,0x3);
            reg4=_mm_mul_pd(reg4,reg1);
            reg3=_mm_mul_pd(reg3,reg1);
            reg4=_mm_shuffle_pd(reg4,reg4,0x1);
            reg1=_mm_addsub_pd(reg3,reg4);

  #endif

            //add to sum
            sumReg=_mm_add_pd(sumReg,reg0);
            sumReg=_mm_add_pd(sumReg,reg1);

            x+=2;
            y+=2;
            n-=2;
          }

          //add beyond where the inner loop stopped
          std::complex<double> sum(0);
          if(n)
            sum=conj(x[0])*y[0];

          //move the sums from the xmm registers to aligned memory on the stack
          std::complex<double> sums[2];
          std::complex<double>* pSums=(std::complex<double>*)((((vcl_size_t)sums)&(~15))+16);
          sumReg=_mm_shuffle_pd(sumReg,sumReg,0x1);//swap real and imag
          _mm_store_pd((double*)pSums,sumReg);

          return sum+pSums[0];
        }
      }

  #endif //defined VIENNACL_WITH_COMPLEX

  #endif //defined VIENNACL_WITH_SSE2

    } //namespace host_based
  } //namespace linalg
} //namespace viennacl

#endif
