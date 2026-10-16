#ifndef AMGCL_RELAXATION_ROCSPARSE_ILU0_HPP
#define AMGCL_RELAXATION_ROCSPARSE_ILU0_HPP

/*
The MIT License

Copyright (c) 2012-2022 Denis Demidov <dennis.demidov@gmail.com>
Copyright (c) 2026 Advanced Micro Devices, Inc.

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
 * \file   amgcl/relaxation/rocsparse_ilu0.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \author Jeff Daily <jeff.daily@amd.com>
 * \brief  Implementation of ILU0 smoother for the ROCm/HIP backend.
 *
 * AMD GPU mirror of relaxation/cusparse_ilu0.hpp. Uses hipSPARSE (the
 * cuSPARSE-compatible ROCm interface): the csrilu02 incomplete-LU
 * factorization plus the generic SpSV triangular solves, with the same
 * persistent SpSV descriptor model as cuSPARSE.
 */

#include <type_traits>

#include <thrust/device_vector.h>

#include <hip/library_types.h>
#include <hipsparse/hipsparse.h>

#include <amgcl/backend/hip.hpp>

namespace amgcl {
namespace relaxation {

namespace detail {

// hipsparseMatDescr_t is typedef void* in hipSPARSE and collides with the
// dense-vector descriptor handled by backend::detail::hip_deleter, so the
// legacy CSR matrix descriptor used below gets its own deleter.
struct hip_mat_descr_deleter {
    void operator()(hipsparseMatDescr_t handle) {
        AMGCL_CALL_HIP( hipsparseDestroyMatDescr(handle) );
    }
};

} // namespace detail

template <class Backend> struct ilu0;

template <typename real>
struct ilu0< backend::hip<real> > {
    typedef real value_type;
    typedef backend::hip<real> Backend;

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
    ilu0( const Matrix &A, const params &prm, const typename Backend::params &bprm)
        : prm(prm), handle(bprm.hipsparse_handle),
          n(backend::rows(A)), nnz(backend::nonzeros(A)),
          ptr(A.ptr, A.ptr + n+1),
          col(A.col, A.col + nnz),
          val(A.val, A.val + nnz),
          y(n)
    {
        // LU decomposition
        std::shared_ptr<std::remove_pointer<hipsparseMatDescr_t>::type> descr_M;
        std::shared_ptr<std::remove_pointer<csrilu02Info_t>::type> info_M;

        {
            hipsparseMatDescr_t descr;
            csrilu02Info_t info;

            AMGCL_CALL_HIP( hipsparseCreateMatDescr(&descr) );
            AMGCL_CALL_HIP( hipsparseSetMatIndexBase(descr, HIPSPARSE_INDEX_BASE_ZERO) );
            AMGCL_CALL_HIP( hipsparseSetMatType(descr, HIPSPARSE_MATRIX_TYPE_GENERAL) );

            AMGCL_CALL_HIP( hipsparseCreateCsrilu02Info(&info) );

            descr_M.reset(descr, detail::hip_mat_descr_deleter());
            info_M.reset(info, backend::detail::hip_deleter());

            int buf_size;

            AMGCL_CALL_HIP(
                    hipsparseXcsrilu02_bufferSize(
                        handle, n, nnz, descr_M.get(),
                        thrust::raw_pointer_cast(&val[0]),
                        thrust::raw_pointer_cast(&ptr[0]),
                        thrust::raw_pointer_cast(&col[0]),
                        info_M.get(), &buf_size
                        )
                    );

            thrust::device_vector<char> bufLU(buf_size);

            // Analysis and incomplete factorization of the system matrix.
            int structural_zero;
            int numerical_zero;

            AMGCL_CALL_HIP(
                    hipsparseXcsrilu02_analysis(
                        handle,
                        n,
                        nnz,
                        descr_M.get(),
                        thrust::raw_pointer_cast(&val[0]),
                        thrust::raw_pointer_cast(&ptr[0]),
                        thrust::raw_pointer_cast(&col[0]),
                        info_M.get(),
                        HIPSPARSE_SOLVE_POLICY_USE_LEVEL,
                        thrust::raw_pointer_cast(&bufLU[0])
                        )
                    );

            precondition(
                    HIPSPARSE_STATUS_ZERO_PIVOT != hipsparseXcsrilu02_zeroPivot(handle, info_M.get(), &structural_zero),
                    "Zero pivot in hipSPARSE ILU0"
                    );

            AMGCL_CALL_HIP(
                    hipsparseXcsrilu02(
                        handle,
                        n,
                        nnz,
                        descr_M.get(),
                        thrust::raw_pointer_cast(&val[0]),
                        thrust::raw_pointer_cast(&ptr[0]),
                        thrust::raw_pointer_cast(&col[0]),
                        info_M.get(),
                        HIPSPARSE_SOLVE_POLICY_USE_LEVEL,
                        thrust::raw_pointer_cast(&bufLU[0])
                        )
                    );
            precondition(
                    HIPSPARSE_STATUS_ZERO_PIVOT != hipsparseXcsrilu02_zeroPivot(handle, info_M.get(), &numerical_zero),
                    "Zero pivot in hipSPARSE ILU0"
                    );

        }

        // Triangular solvers
        const real alpha = 1;
        thrust::device_vector<value_type> t(n);

        descr_y.reset(
                backend::detail::hip_vector_description(y),
                backend::detail::hip_deleter()
                );

        std::shared_ptr<std::remove_pointer<hipsparseDnVecDescr_t>::type> descr_t(
                backend::detail::hip_vector_description(t),
                backend::detail::hip_deleter()
                );

        hipsparseFillMode_t fill_lower    = HIPSPARSE_FILL_MODE_LOWER;
        hipsparseFillMode_t fill_upper    = HIPSPARSE_FILL_MODE_UPPER;
        hipsparseDiagType_t diag_unit     = HIPSPARSE_DIAG_TYPE_UNIT;
        hipsparseDiagType_t diag_non_unit = HIPSPARSE_DIAG_TYPE_NON_UNIT;

        // Triangular solver for L
        {
            descr_L.reset(
                    backend::detail::hip_matrix_description(n, n, nnz, ptr, col, val),
                    backend::detail::hip_deleter()
                    );

            AMGCL_CALL_HIP(
                    hipsparseSpMatSetAttribute(
                        descr_L.get(),
                        HIPSPARSE_SPMAT_FILL_MODE,
                        &fill_lower,
                        sizeof(fill_lower)
                        )
                    );

            AMGCL_CALL_HIP(
                    hipsparseSpMatSetAttribute(
                        descr_L.get(),
                        HIPSPARSE_SPMAT_DIAG_TYPE,
                        &diag_unit,
                        sizeof(diag_unit)
                        )
                    );


            size_t buf_size;

            hipsparseSpSVDescr_t desc;
            AMGCL_CALL_HIP( hipsparseSpSV_createDescr(&desc) );
            descr_SL.reset(desc, backend::detail::hip_deleter());

            AMGCL_CALL_HIP(
                    hipsparseSpSV_bufferSize(
                        handle,
                        HIPSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        descr_L.get(),
                        descr_t.get(),
                        descr_y.get(),
                        backend::detail::hip_datatype<real>(),
                        HIPSPARSE_SPSV_ALG_DEFAULT,
                        descr_SL.get(),
                        &buf_size
                        )
                    );

            bufL.resize(buf_size);

            AMGCL_CALL_HIP(
                    hipsparseSpSV_analysis(
                        handle,
                        HIPSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        descr_L.get(),
                        descr_t.get(),
                        descr_y.get(),
                        backend::detail::hip_datatype<real>(),
                        HIPSPARSE_SPSV_ALG_DEFAULT,
                        descr_SL.get(),
                        thrust::raw_pointer_cast(&bufL[0])
                        )
                    );
        }

        // Triangular solver for U
        {
            descr_U.reset(
                    backend::detail::hip_matrix_description(n, n, nnz, ptr, col, val),
                    backend::detail::hip_deleter()
                    );

            AMGCL_CALL_HIP(
                    hipsparseSpMatSetAttribute(
                        descr_U.get(),
                        HIPSPARSE_SPMAT_FILL_MODE,
                        &fill_upper,
                        sizeof(fill_upper)
                        )
                    );

            AMGCL_CALL_HIP(
                    hipsparseSpMatSetAttribute(
                        descr_U.get(),
                        HIPSPARSE_SPMAT_DIAG_TYPE,
                        &diag_non_unit,
                        sizeof(diag_non_unit)
                        )
                    );


            size_t buf_size;

            hipsparseSpSVDescr_t desc;
            AMGCL_CALL_HIP( hipsparseSpSV_createDescr(&desc) );
            descr_SU.reset(desc, backend::detail::hip_deleter());

            AMGCL_CALL_HIP(
                    hipsparseSpSV_bufferSize(
                        handle,
                        HIPSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        descr_U.get(),
                        descr_y.get(),
                        descr_t.get(),
                        backend::detail::hip_datatype<real>(),
                        HIPSPARSE_SPSV_ALG_DEFAULT,
                        descr_SU.get(),
                        &buf_size
                        )
                    );

            bufU.resize(buf_size);

            AMGCL_CALL_HIP(
                    hipsparseSpSV_analysis(
                        handle,
                        HIPSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        descr_U.get(),
                        descr_y.get(),
                        descr_t.get(),
                        backend::detail::hip_datatype<real>(),
                        HIPSPARSE_SPSV_ALG_DEFAULT,
                        descr_SU.get(),
                        thrust::raw_pointer_cast(&bufU[0])
                        )
                    );
        }
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
        // This is incomplete, as hipsparse structs are opaque.
        return
            backend::bytes(ptr) +
            backend::bytes(col) +
            backend::bytes(val) +
            backend::bytes(y) +
            backend::bytes(bufL) +
            backend::bytes(bufU)
            ;
    }

    private:
        hipsparseHandle_t handle;
        int n, nnz;

        thrust::device_vector<int> ptr, col;
        thrust::device_vector<value_type> val;
        mutable thrust::device_vector<value_type> y;

        std::shared_ptr<std::remove_pointer<hipsparseSpMatDescr_t>::type> descr_L, descr_U;
        std::shared_ptr<std::remove_pointer<hipsparseSpSVDescr_t>::type>  descr_SL, descr_SU;
        std::shared_ptr<std::remove_pointer<hipsparseDnVecDescr_t>::type> descr_y;
        mutable thrust::device_vector<char> bufL, bufU;

        template <class VectorX>
        void solve(VectorX &x) const {
            value_type alpha = 1;

            std::shared_ptr<std::remove_pointer<hipsparseDnVecDescr_t>::type> descr_x(
                    backend::detail::hip_vector_description(x),
                    backend::detail::hip_deleter()
                    );

            // Solve L * y = x
            AMGCL_CALL_HIP(
                    hipsparseSpSV_solve(
                        handle,
                        HIPSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        descr_L.get(),
                        descr_x.get(),
                        descr_y.get(),
                        backend::detail::hip_datatype<real>(),
                        HIPSPARSE_SPSV_ALG_DEFAULT,
                        descr_SL.get()
                        )
                    );

            // Solve U * x = y
            AMGCL_CALL_HIP(
                    hipsparseSpSV_solve(
                        handle,
                        HIPSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        descr_U.get(),
                        descr_y.get(),
                        descr_x.get(),
                        backend::detail::hip_datatype<real>(),
                        HIPSPARSE_SPSV_ALG_DEFAULT,
                        descr_SU.get()
                        )
                    );
        }

        // hipSPARSE, like cuSPARSE, only ships type-suffixed csrilu02 routines
        // (Dcsrilu02/Scsrilu02). These thin overloads dispatch on the value
        // type so the factorization code above can be written once.
        static hipsparseStatus_t hipsparseXcsrilu02_bufferSize(
                hipsparseHandle_t handle, int m, int nnz,
                const hipsparseMatDescr_t descrA, double *csrSortedValA,
                const int *csrSortedRowPtrA, const int *csrSortedColIndA,
                csrilu02Info_t info, int *pBufferSizeInBytes)
        {
            return hipsparseDcsrilu02_bufferSize(handle, m, nnz, descrA,
                csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, info,
                pBufferSizeInBytes);
        }

        static hipsparseStatus_t hipsparseXcsrilu02_bufferSize(
                hipsparseHandle_t handle, int m, int nnz,
                const hipsparseMatDescr_t descrA, float *csrSortedValA,
                const int *csrSortedRowPtrA, const int *csrSortedColIndA,
                csrilu02Info_t info, int *pBufferSizeInBytes)
        {
            return hipsparseScsrilu02_bufferSize(handle, m, nnz, descrA,
                csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, info,
                pBufferSizeInBytes);
        }

        static hipsparseStatus_t hipsparseXcsrilu02_analysis(
                hipsparseHandle_t handle, int m, int nnz,
                const hipsparseMatDescr_t descrA, const double *csrSortedValA,
                const int *csrSortedRowPtrA, const int *csrSortedColIndA,
                csrilu02Info_t info, hipsparseSolvePolicy_t policy, void *pBuffer)
        {
            return hipsparseDcsrilu02_analysis(handle, m, nnz, descrA,
                csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, info, policy,
                pBuffer);
        }

        static hipsparseStatus_t hipsparseXcsrilu02_analysis(
                hipsparseHandle_t handle, int m, int nnz,
                const hipsparseMatDescr_t descrA, const float *csrSortedValA,
                const int *csrSortedRowPtrA, const int *csrSortedColIndA,
                csrilu02Info_t info, hipsparseSolvePolicy_t policy, void *pBuffer)
        {
            return hipsparseScsrilu02_analysis(handle, m, nnz, descrA,
                csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, info, policy,
                pBuffer);
        }

        static hipsparseStatus_t hipsparseXcsrilu02(
                hipsparseHandle_t handle, int m, int nnz,
                const hipsparseMatDescr_t descrA, double *csrSortedValA_valM,
                const int *csrSortedRowPtrA, const int *csrSortedColIndA,
                csrilu02Info_t info, hipsparseSolvePolicy_t policy, void *pBuffer)
        {
            return hipsparseDcsrilu02(handle, m, nnz, descrA,
                csrSortedValA_valM, csrSortedRowPtrA, csrSortedColIndA, info,
                policy, pBuffer);
        }

        static hipsparseStatus_t hipsparseXcsrilu02(
                hipsparseHandle_t handle, int m, int nnz,
                const hipsparseMatDescr_t descrA, float *csrSortedValA_valM,
                const int *csrSortedRowPtrA, const int *csrSortedColIndA,
                csrilu02Info_t info, hipsparseSolvePolicy_t policy, void *pBuffer)
        {
            return hipsparseScsrilu02(handle, m, nnz, descrA,
                csrSortedValA_valM, csrSortedRowPtrA, csrSortedColIndA, info,
                policy, pBuffer);
        }
};

} // namespace relaxation
} // namespace amgcl

#endif
