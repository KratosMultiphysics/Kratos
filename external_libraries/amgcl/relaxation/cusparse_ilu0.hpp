#ifndef AMGCL_RELAXATION_CUSPARSE_ILU0_HPP
#define AMGCL_RELAXATION_CUSPARSE_ILU0_HPP

/*
The MIT License

Copyright (c) 2012-2019 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   amgcl/relaxation/cusparse_ilu0.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Implementation of ILU0 smoother for CUDA backend.
 */

#include <type_traits>

#include <thrust/device_vector.h>
#include <cusparse_v2.h>

#include <amgcl/backend/cuda.hpp>

namespace amgcl {
namespace relaxation {

template <class Backend> struct ilu0;

template <typename real>
struct ilu0< backend::cuda<real> > {
    typedef real value_type;
    typedef backend::cuda<double> Backend;

    struct params {
        /// Damping factor.
        float damping;

        params() : damping(1) {}

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, damping)
        {
            check_params(p, {"damping"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, damping);
        }
#endif
    } prm;

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    ilu0( const Matrix &A, const params &, const typename Backend::params &bprm)
        : prm(prm), handle(bprm.cusparse_handle),
          n(backend::rows(A)), nnz(backend::nonzeros(A)),
          ptr(A.ptr, A.ptr + n+1),
          col(A.col, A.col + nnz),
          val(A.val, A.val + nnz),
          y(n)
    {
        // Create matrix descriptors.
        {
            cusparseMatDescr_t descr;

            AMGCL_CALL_CUDA( cusparseCreateMatDescr(&descr) );
            AMGCL_CALL_CUDA( cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO) );
            AMGCL_CALL_CUDA( cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL) );

            descr_M.reset(descr, backend::detail::cuda_deleter());
        }
        {
            cusparseMatDescr_t descr;

            AMGCL_CALL_CUDA( cusparseCreateMatDescr(&descr) );
            AMGCL_CALL_CUDA( cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO) );
            AMGCL_CALL_CUDA( cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL) );
            AMGCL_CALL_CUDA( cusparseSetMatFillMode(descr, CUSPARSE_FILL_MODE_LOWER) );
            AMGCL_CALL_CUDA( cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_UNIT) );

            descr_L.reset(descr, backend::detail::cuda_deleter());
        }
        {
            cusparseMatDescr_t descr;

            AMGCL_CALL_CUDA( cusparseCreateMatDescr(&descr) );
            AMGCL_CALL_CUDA( cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO) );
            AMGCL_CALL_CUDA( cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL) );
            AMGCL_CALL_CUDA( cusparseSetMatFillMode(descr, CUSPARSE_FILL_MODE_UPPER) );
            AMGCL_CALL_CUDA( cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_NON_UNIT) );

            descr_U.reset(descr, backend::detail::cuda_deleter());
        }

        // Create info structures.
        {
            csrilu02Info_t info;
            AMGCL_CALL_CUDA( cusparseCreateCsrilu02Info(&info) );
            info_M.reset(info, backend::detail::cuda_deleter());
        }
        {
            csrsv2Info_t info;
            AMGCL_CALL_CUDA( cusparseCreateCsrsv2Info(&info) );
            info_L.reset(info, backend::detail::cuda_deleter());
        }
        {
            csrsv2Info_t info;
            AMGCL_CALL_CUDA( cusparseCreateCsrsv2Info(&info) );
            info_U.reset(info, backend::detail::cuda_deleter());
        }

        // Allocate scratch buffer.
        {
            const cusparseOperation_t trans_L  = CUSPARSE_OPERATION_NON_TRANSPOSE;
            const cusparseOperation_t trans_U  = CUSPARSE_OPERATION_NON_TRANSPOSE;

            int buf_size_M;
            int buf_size_L;
            int buf_size_U;

            AMGCL_CALL_CUDA(
                    cusparseXcsrilu02_bufferSize(
                        handle, n, nnz, descr_M.get(),
                        thrust::raw_pointer_cast(&val[0]),
                        thrust::raw_pointer_cast(&ptr[0]),
                        thrust::raw_pointer_cast(&col[0]),
                        info_M.get(), &buf_size_M
                        )
                    );
            AMGCL_CALL_CUDA(
                    cusparseXcsrsv2_bufferSize(
                        handle, trans_L, n, nnz, descr_L.get(),
                        thrust::raw_pointer_cast(&val[0]),
                        thrust::raw_pointer_cast(&ptr[0]),
                        thrust::raw_pointer_cast(&col[0]),
                        info_L.get(), &buf_size_L
                        )
                    );
            AMGCL_CALL_CUDA(
                    cusparseXcsrsv2_bufferSize(
                        handle, trans_U, n, nnz, descr_U.get(),
                        thrust::raw_pointer_cast(&val[0]),
                        thrust::raw_pointer_cast(&ptr[0]),
                        thrust::raw_pointer_cast(&col[0]),
                        info_U.get(), &buf_size_U
                        )
                    );

            buf.resize( std::max(buf_size_M, std::max(buf_size_L, buf_size_U)) );
        }

        // Analysis and incomplete factorization of the system matrix.
        int structural_zero;
        int numerical_zero;

        AMGCL_CALL_CUDA(
                cusparseXcsrilu02_analysis(handle, n, nnz, descr_M.get(),
                    thrust::raw_pointer_cast(&val[0]),
                    thrust::raw_pointer_cast(&ptr[0]),
                    thrust::raw_pointer_cast(&col[0]),
                    info_M.get(), policy_M,
                    thrust::raw_pointer_cast(&buf[0])
                    )
                );

        precondition(
                CUSPARSE_STATUS_ZERO_PIVOT != cusparseXcsrilu02_zeroPivot(handle, info_M.get(), &structural_zero),
                "Zero pivot in cuSPARSE ILU0"
                );

        AMGCL_CALL_CUDA(
                cusparseXcsrsv2_analysis(
                    handle, trans_L, n, nnz, descr_L.get(),
                    thrust::raw_pointer_cast(&val[0]),
                    thrust::raw_pointer_cast(&ptr[0]),
                    thrust::raw_pointer_cast(&col[0]),
                    info_L.get(), policy_L,
                    thrust::raw_pointer_cast(&buf[0])
                    )
                );

        AMGCL_CALL_CUDA(
                cusparseXcsrsv2_analysis(
                    handle, trans_U, n, nnz, descr_U.get(),
                    thrust::raw_pointer_cast(&val[0]),
                    thrust::raw_pointer_cast(&ptr[0]),
                    thrust::raw_pointer_cast(&col[0]),
                    info_U.get(), policy_U,
                    thrust::raw_pointer_cast(&buf[0])
                    )
                );

        AMGCL_CALL_CUDA(
                cusparseXcsrilu02(
                    handle, n, nnz, descr_M.get(),
                    thrust::raw_pointer_cast(&val[0]),
                    thrust::raw_pointer_cast(&ptr[0]),
                    thrust::raw_pointer_cast(&col[0]),
                    info_M.get(), policy_M,
                    thrust::raw_pointer_cast(&buf[0])
                    )
                );
        precondition(
              CUSPARSE_STATUS_ZERO_PIVOT != cusparseXcsrilu02_zeroPivot(handle, info_M.get(), &numerical_zero),
              "Zero pivot in cuSPARSE ILU0"
              );
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        solve(tmp);
        backend::axpby(prm.damping, tmp, 1, x);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        solve(tmp);
        backend::axpby(prm.damping, tmp, 1, x);
    }

    template <class Matrix, class VectorRHS, class VectorX>
    void apply(const Matrix &A, const VectorRHS &rhs, VectorX &x) const
    {
        backend::copy(rhs, x);
        solve(x);
    }

    size_t bytes() const {
        // This is incomplete, as cusparse structs are opaque.
        return
            backend::bytes(ptr) +
            backend::bytes(col) +
            backend::bytes(val) +
            backend::bytes(y) +
            backend::bytes(buf);
    }

    private:
        static const cusparseSolvePolicy_t policy_M = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
        static const cusparseSolvePolicy_t policy_L = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        static const cusparseSolvePolicy_t policy_U = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        static const cusparseOperation_t   trans_L  = CUSPARSE_OPERATION_NON_TRANSPOSE;
        static const cusparseOperation_t   trans_U  = CUSPARSE_OPERATION_NON_TRANSPOSE;

        cusparseHandle_t handle;
        int n, nnz;

        std::shared_ptr<std::remove_pointer<cusparseMatDescr_t>::type> descr_M, descr_L, descr_U;
        std::shared_ptr<std::remove_pointer<csrilu02Info_t>::type> info_M;
        std::shared_ptr<std::remove_pointer<csrsv2Info_t>::type>  info_L, info_U;

        thrust::device_vector<int> ptr, col;
        thrust::device_vector<value_type> val;
        mutable thrust::device_vector<value_type> y;
        mutable thrust::device_vector<char> buf;


        template <class VectorX>
        void solve(VectorX &x) const {
            value_type alpha = 1;

            // Solve L * y = x
            AMGCL_CALL_CUDA(
                    cusparseXcsrsv2_solve(
                        handle, trans_L, n, nnz, &alpha, descr_L.get(),
                        thrust::raw_pointer_cast(&val[0]),
                        thrust::raw_pointer_cast(&ptr[0]),
                        thrust::raw_pointer_cast(&col[0]),
                        info_L.get(),
                        thrust::raw_pointer_cast(&x[0]),
                        thrust::raw_pointer_cast(&y[0]),
                        policy_L,
                        thrust::raw_pointer_cast(&buf[0])
                        )
                    );

            // Solve U * x = y
            AMGCL_CALL_CUDA(
                    cusparseXcsrsv2_solve(
                        handle, trans_U, n, nnz, &alpha, descr_U.get(),
                        thrust::raw_pointer_cast(&val[0]),
                        thrust::raw_pointer_cast(&ptr[0]),
                        thrust::raw_pointer_cast(&col[0]),
                        info_U.get(),
                        thrust::raw_pointer_cast(&y[0]),
                        thrust::raw_pointer_cast(&x[0]),
                        policy_U,
                        thrust::raw_pointer_cast(&buf[0])
                        )
                    );
        }


        static cusparseStatus_t cusparseXcsrilu02_bufferSize(
                cusparseHandle_t handle,
                int m,
                int nnz,
                const cusparseMatDescr_t descrA,
                double *csrSortedValA,
                const int *csrSortedRowPtrA,
                const int *csrSortedColIndA,
                csrilu02Info_t info,
                int *pBufferSizeInBytes)
        {
            return cusparseDcsrilu02_bufferSize(
                handle, m, nnz, descrA, csrSortedValA, csrSortedRowPtrA,
                csrSortedColIndA, info, pBufferSizeInBytes);
        }

        static cusparseStatus_t cusparseXcsrilu02_bufferSize(
                cusparseHandle_t handle,
                int m,
                int nnz,
                const cusparseMatDescr_t descrA,
                float *csrSortedValA,
                const int *csrSortedRowPtrA,
                const int *csrSortedColIndA,
                csrilu02Info_t info,
                int *pBufferSizeInBytes)
        {
            return cusparseScsrilu02_bufferSize(
                handle, m, nnz, descrA, csrSortedValA, csrSortedRowPtrA,
                csrSortedColIndA, info, pBufferSizeInBytes);
        }

        static cusparseStatus_t cusparseXcsrsv2_bufferSize(
                cusparseHandle_t handle,
                cusparseOperation_t transA,
                int m,
                int nnz,
                const cusparseMatDescr_t descrA,
                double *csrSortedValA,
                const int *csrSortedRowPtrA,
                const int *csrSortedColIndA,
                csrsv2Info_t info,
                int *pBufferSizeInBytes
                )
        {
            return cusparseDcsrsv2_bufferSize(
                handle, transA, m, nnz, descrA, csrSortedValA,
                csrSortedRowPtrA, csrSortedColIndA, info, pBufferSizeInBytes);
        }

        static cusparseStatus_t cusparseXcsrsv2_bufferSize(
                cusparseHandle_t handle,
                cusparseOperation_t transA,
                int m,
                int nnz,
                const cusparseMatDescr_t descrA,
                float *csrSortedValA,
                const int *csrSortedRowPtrA,
                const int *csrSortedColIndA,
                csrsv2Info_t info,
                int *pBufferSizeInBytes
                )
        {
            return cusparseScsrsv2_bufferSize(
                handle, transA, m, nnz, descrA, csrSortedValA,
                csrSortedRowPtrA, csrSortedColIndA, info, pBufferSizeInBytes);
        }

        static cusparseStatus_t cusparseXcsrilu02_analysis(
                cusparseHandle_t handle,
                int m,
                int nnz,
                const cusparseMatDescr_t descrA,
                const double *csrSortedValA,
                const int *csrSortedRowPtrA,
                const int *csrSortedColIndA,
                csrilu02Info_t info,
                cusparseSolvePolicy_t policy,
                void *pBuffer
                )
        {
            return cusparseDcsrilu02_analysis(
                handle, m, nnz, descrA, csrSortedValA,
                csrSortedRowPtrA, csrSortedColIndA, info, policy, pBuffer);
        }

        static cusparseStatus_t cusparseXcsrilu02_analysis(
                cusparseHandle_t handle,
                int m,
                int nnz,
                const cusparseMatDescr_t descrA,
                const float *csrSortedValA,
                const int *csrSortedRowPtrA,
                const int *csrSortedColIndA,
                csrilu02Info_t info,
                cusparseSolvePolicy_t policy,
                void *pBuffer
                )
        {
            return cusparseScsrilu02_analysis(
                handle, m, nnz, descrA, csrSortedValA,
                csrSortedRowPtrA, csrSortedColIndA, info, policy, pBuffer);
        }

        static cusparseStatus_t cusparseXcsrsv2_analysis(
                cusparseHandle_t handle,
                cusparseOperation_t transA,
                int m,
                int nnz,
                const cusparseMatDescr_t descrA,
                const double *csrSortedValA,
                const int *csrSortedRowPtrA,
                const int *csrSortedColIndA,
                csrsv2Info_t info,
                cusparseSolvePolicy_t policy,
                void *pBuffer
                )
        {
            return cusparseDcsrsv2_analysis(
                    handle, transA, m, nnz, descrA, csrSortedValA,
                    csrSortedRowPtrA, csrSortedColIndA, info, policy, pBuffer);
        }

        static cusparseStatus_t cusparseXcsrsv2_analysis(
                cusparseHandle_t handle,
                cusparseOperation_t transA,
                int m,
                int nnz,
                const cusparseMatDescr_t descrA,
                const float *csrSortedValA,
                const int *csrSortedRowPtrA,
                const int *csrSortedColIndA,
                csrsv2Info_t info,
                cusparseSolvePolicy_t policy,
                void *pBuffer
                )
        {
            return cusparseScsrsv2_analysis(
                    handle, transA, m, nnz, descrA, csrSortedValA,
                    csrSortedRowPtrA, csrSortedColIndA, info, policy, pBuffer);
        }

        static cusparseStatus_t cusparseXcsrilu02(cusparseHandle_t handle,
                int m,
                int nnz,
                const cusparseMatDescr_t descrA,
                double *csrSortedValA_valM,
                const int *csrSortedRowPtrA,
                const int *csrSortedColIndA,
                csrilu02Info_t info,
                cusparseSolvePolicy_t policy,
                void *pBuffer
                )
        {
            return cusparseDcsrilu02(handle, m, nnz, descrA,
                    csrSortedValA_valM, csrSortedRowPtrA, csrSortedColIndA,
                    info, policy, pBuffer);
        }

        static cusparseStatus_t cusparseXcsrilu02(cusparseHandle_t handle,
                int m,
                int nnz,
                const cusparseMatDescr_t descrA,
                float *csrSortedValA_valM,
                const int *csrSortedRowPtrA,
                const int *csrSortedColIndA,
                csrilu02Info_t info,
                cusparseSolvePolicy_t policy,
                void *pBuffer
                )
        {
            return cusparseScsrilu02(handle, m, nnz, descrA,
                    csrSortedValA_valM, csrSortedRowPtrA, csrSortedColIndA,
                    info, policy, pBuffer);
        }

        static cusparseStatus_t cusparseXcsrsv2_solve(
                cusparseHandle_t handle,
                cusparseOperation_t transA,
                int m,
                int nnz,
                const double *alpha,
                const cusparseMatDescr_t descrA,
                const double *csrSortedValA,
                const int *csrSortedRowPtrA,
                const int *csrSortedColIndA,
                csrsv2Info_t info,
                const double *f,
                double *x,
                cusparseSolvePolicy_t policy,
                void *pBuffer
                )
        {
            return cusparseDcsrsv2_solve(
                    handle, transA, m,
                    nnz, alpha, descrA, csrSortedValA, csrSortedRowPtrA,
                    csrSortedColIndA, info, f, x, policy, pBuffer);
        }

        static cusparseStatus_t cusparseXcsrsv2_solve(
                cusparseHandle_t handle,
                cusparseOperation_t transA,
                int m,
                int nnz,
                const float *alpha,
                const cusparseMatDescr_t descrA,
                const float *csrSortedValA,
                const int *csrSortedRowPtrA,
                const int *csrSortedColIndA,
                csrsv2Info_t info,
                const float *f,
                float *x,
                cusparseSolvePolicy_t policy,
                void *pBuffer
                )
        {
            return cusparseScsrsv2_solve(
                    handle, transA, m,
                    nnz, alpha, descrA, csrSortedValA, csrSortedRowPtrA,
                    csrSortedColIndA, info, f, x, policy, pBuffer);
        }
};

} // namespace relaxation
} // namespace amgcl

#endif
