#ifndef VIENNACL_LINALG_QR_METHOD_COMMON_HPP
#define VIENNACL_LINALG_QR_METHOD_COMMON_HPP

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

#include <cmath>

#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/linalg/opencl/kernels/svd.hpp"
#include "viennacl/meta/result_of.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

/** @file viennacl/linalg/qr-method-common.hpp
    @brief Common routines used for the QR method and SVD. Experimental.
*/

namespace viennacl
{
  namespace linalg
  {
    const std::string SVD_BIDIAG_PACK_KERNEL = "bidiag_pack";
    const std::string SVD_HOUSEHOLDER_UPDATE_A_LEFT_KERNEL = "house_update_A_left";
    const std::string SVD_HOUSEHOLDER_UPDATE_A_RIGHT_KERNEL = "house_update_A_right";
    const std::string SVD_HOUSEHOLDER_UPDATE_QL_KERNEL = "house_update_QL";
    const std::string SVD_HOUSEHOLDER_UPDATE_QR_KERNEL = "house_update_QR";
    const std::string SVD_COPY_COL_KERNEL = "copy_col";
    const std::string SVD_COPY_ROW_KERNEL = "copy_row";
    const std::string SVD_MATRIX_TRANSPOSE_KERNEL = "transpose_inplace";
    const std::string SVD_INVERSE_SIGNS_KERNEL = "inverse_signs";
    const std::string SVD_GIVENS_PREV_KERNEL = "givens_prev";
    const std::string SVD_GIVENS_NEXT_KERNEL = "givens_next";
    const std::string SVD_FINAL_ITER_UPDATE_KERNEL = "final_iter_update";
    const std::string SVD_UPDATE_QR_COLUMN_KERNEL = "update_qr_column";

    namespace detail
    {
      //static const float EPS = 0.00001f;
      //static const vcl_size_t ITER_MAX = 50;

      static const double EPS = 1e-10;
      static const vcl_size_t ITER_MAX = 50;

      template <typename SCALARTYPE>
      SCALARTYPE pythag(SCALARTYPE a, SCALARTYPE b)
      {
        return std::sqrt(a*a + b*b);
      }

      template <typename SCALARTYPE>
      SCALARTYPE sign(SCALARTYPE val)
      {
          return (val >= 0) ? SCALARTYPE(1) : SCALARTYPE(-1);
      }

      // DEPRECATED: Replace with viennacl::linalg::norm_2
      template <typename VectorType>
      typename VectorType::value_type norm_lcl(VectorType const & x, vcl_size_t size)
      {
        typename VectorType::value_type x_norm = 0.0;
        for(vcl_size_t i = 0; i < size; i++)
          x_norm += std::pow(x[i], 2);
        return std::sqrt(x_norm);
      }

      template <typename VectorType>
      void normalize(VectorType & x, vcl_size_t size)
      {
        typename VectorType::value_type x_norm = norm_lcl(x, size);
        for(vcl_size_t i = 0; i < size; i++)
            x[i] /= x_norm;
      }



      template <typename VectorType>
      void householder_vector(VectorType & v, vcl_size_t start)
      {
        typedef typename VectorType::value_type    ScalarType;
        ScalarType x_norm = norm_lcl(v, v.size());
        ScalarType alpha = -sign(v[start]) * x_norm;
        v[start] += alpha;
        normalize(v, v.size());
      }

      template <typename MatrixType>
      void transpose(MatrixType & A)
      {
        typedef typename MatrixType::value_type                                   ScalarType;
        typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;

        viennacl::ocl::kernel & kernel = viennacl::ocl::get_kernel(viennacl::linalg::opencl::kernels::svd<CPU_ScalarType>::program_name(), SVD_MATRIX_TRANSPOSE_KERNEL);

        viennacl::ocl::enqueue(kernel(A,
                                      static_cast<cl_uint>(A.internal_size1()),
                                      static_cast<cl_uint>(A.internal_size2())
                                     )
                              );
      }



      template <typename T>
      void cdiv(T xr, T xi, T yr, T yi, T& cdivr, T& cdivi)
      {
          // Complex scalar division.
          T r;
          T d;
          if (std::fabs(yr) > std::fabs(yi))
          {
              r = yi / yr;
              d = yr + r * yi;
              cdivr = (xr + r * xi) / d;
              cdivi = (xi - r * xr) / d;
          }
          else
          {
              r = yr / yi;
              d = yi + r * yr;
              cdivr = (r * xr + xi) / d;
              cdivi = (r * xi - xr) / d;
          }
      }


      template <typename SCALARTYPE, unsigned int ALIGNMENT>
      void copy_vec(viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& A,
                    viennacl::vector<SCALARTYPE, ALIGNMENT>& V,
                    vcl_size_t row_start,
                    vcl_size_t col_start,
                    bool copy_col
      )
      {

        std::string kernel_name = copy_col ? SVD_COPY_COL_KERNEL : SVD_COPY_ROW_KERNEL;
        viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(A).context());
        viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::svd<SCALARTYPE>::program_name(), kernel_name);

        viennacl::ocl::enqueue(kernel(
                                      A,
                                      V,
                                      static_cast<cl_uint>(row_start),
                                      static_cast<cl_uint>(col_start),
                                      copy_col ? static_cast<cl_uint>(A.size1())
                                               : static_cast<cl_uint>(A.size2()),
                                      static_cast<cl_uint>(A.internal_size2())
                              ));

      }


      template<typename SCALARTYPE, unsigned int ALIGNMENT>
      void prepare_householder_vector(
                                    viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& A,
                                    viennacl::vector<SCALARTYPE, ALIGNMENT>& D,
                                    vcl_size_t size,
                                    vcl_size_t row_start,
                                    vcl_size_t col_start,
                                    vcl_size_t start,
                                    bool is_column
                                    )
      {
        boost::numeric::ublas::vector<SCALARTYPE> tmp = boost::numeric::ublas::scalar_vector<SCALARTYPE>(size, 0);

        copy_vec(A, D, row_start, col_start, is_column);
        fast_copy(D.begin(), D.begin() + vcl_ptrdiff_t(size - start), tmp.begin() + start);

        //std::cout << "1: " << tmp << "\n";

        detail::householder_vector(tmp, start);
        fast_copy(tmp, D);

        //std::cout << "2: "  << D << "\n";
      }

      template <typename SCALARTYPE, unsigned int ALIGNMENT, typename VectorType>
      void bidiag_pack(viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& A,
                       VectorType & dh,
                       VectorType & sh
                      )
      {
        viennacl::vector<SCALARTYPE, ALIGNMENT> D(dh.size());
        viennacl::vector<SCALARTYPE, ALIGNMENT> S(sh.size());

        viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(A).context());
        viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::svd<SCALARTYPE>::program_name(), SVD_BIDIAG_PACK_KERNEL);

        viennacl::ocl::enqueue(kernel(
                                      A,
                                      D,
                                      S,
                                      static_cast<cl_uint>(A.size1()),
                                      static_cast<cl_uint>(A.size2()),
                                      static_cast<cl_uint>(A.internal_size2())
                                    ));

        fast_copy(D, dh);
        fast_copy(S, sh);
      }

    }
  }
}

#endif
