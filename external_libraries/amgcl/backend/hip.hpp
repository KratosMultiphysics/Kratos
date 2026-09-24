#ifndef AMGCL_BACKEND_HIP_HPP
#define AMGCL_BACKEND_HIP_HPP

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
 * \file   amgcl/backend/hip.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \author Jeff Daily <jeff.daily@amd.com>
 * \brief  ROCm/HIP backend.
 *
 * AMD GPU mirror of backend/cuda.hpp. Uses hipSPARSE (the cuSPARSE-compatible
 * ROCm interface) for the CSR SpMV and rocThrust for the dense vector kernels.
 * hipSPARSE exposes the same generic SpMV API as cuSPARSE, so this header
 * follows the CUDA backend structure with only the cusparse->hipsparse symbol
 * swap; rocThrust shares the thrust:: API and header paths, so the Thrust code
 * is unchanged.
 */

#include <type_traits>
#include <memory>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/solver/skyline_lu.hpp>
#include <amgcl/util.hpp>

#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/gather.h>
#include <thrust/scatter.h>
#include <thrust/for_each.h>
#include <thrust/inner_product.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>

#include <hip/hip_runtime.h>
#include <hip/library_types.h>
#include <hipsparse/hipsparse.h>

namespace amgcl {

namespace solver {

/** Wrapper around solver::skyline_lu for use with the HIP backend.
 * Copies the rhs to the host memory, solves the problem using the host CPU,
 * then copies the solution back to the compute device(s).
 */
template <class T>
struct hip_skyline_lu : solver::skyline_lu<T> {
    typedef solver::skyline_lu<T> Base;

    mutable std::vector<T> _rhs, _x;

    template <class Matrix, class Params>
    hip_skyline_lu(const Matrix &A, const Params&)
        : Base(*A), _rhs(backend::rows(*A)), _x(backend::rows(*A))
    { }

    template <class Vec1, class Vec2>
    void operator()(const Vec1 &rhs, Vec2 &x) const {
        thrust::copy(rhs.begin(), rhs.end(), _rhs.begin());
        static_cast<const Base*>(this)->operator()(_rhs, _x);
        thrust::copy(_x.begin(), _x.end(), x.begin());
    }

    size_t bytes() const {
        return
            backend::bytes(*static_cast<const Base*>(this)) +
            backend::bytes(_rhs) +
            backend::bytes(_x);
    }
};

}

namespace backend {
namespace detail {

inline void hip_check(hipsparseStatus_t rc, const char *file, int line) {
    if (rc != HIPSPARSE_STATUS_SUCCESS) {
        std::ostringstream msg;
        msg << "HIP error " << rc << " at \"" << file << ":" << line;
        precondition(false, msg.str());
    }
}

inline void hip_check(hipError_t rc, const char *file, int line) {
    if (rc != hipSuccess) {
        std::ostringstream msg;
        msg << "HIP error " << rc << " at \"" << file << ":" << line;
        precondition(false, msg.str());
    }
}

#define AMGCL_CALL_HIP(rc)                                                     \
    amgcl::backend::detail::hip_check(rc, __FILE__, __LINE__)

struct hip_deleter {
    // hipsparseMatDescr_t and hipsparseDnVecDescr_t are both typedef'd to void*
    // in hipSPARSE (unlike cuSPARSE, where they are distinct struct pointers),
    // so they cannot be told apart by overload. The legacy hipsparseMatDescr_t
    // (used only by the ILU0 relaxation) gets its own deleter there; here the
    // void* overload destroys a dense-vector descriptor.
    void operator()(hipsparseSpMatDescr_t handle) {
        AMGCL_CALL_HIP( hipsparseDestroySpMat(handle) );
    }

    void operator()(hipsparseDnVecDescr_t handle) {
        AMGCL_CALL_HIP( hipsparseDestroyDnVec(handle) );
    }

    void operator()(hipEvent_t handle) {
        AMGCL_CALL_HIP( hipEventDestroy(handle) );
    }

    void operator()(csrilu02Info_t handle) {
        AMGCL_CALL_HIP( hipsparseDestroyCsrilu02Info(handle) );
    }

    void operator()(hipsparseSpSVDescr_t handle) {
        AMGCL_CALL_HIP( hipsparseSpSV_destroyDescr(handle) );
    }
};


template <typename real>
hipDataType hip_datatype() {
    if (sizeof(real) == sizeof(float))
        return HIP_R_32F;
    else
        return HIP_R_64F;
}

template <typename real>
hipsparseDnVecDescr_t hip_vector_description(thrust::device_vector<real> &x) {
    hipsparseDnVecDescr_t desc;
    AMGCL_CALL_HIP(
            hipsparseCreateDnVec(
                &desc,
                x.size(),
                thrust::raw_pointer_cast(&x[0]),
                hip_datatype<real>()
                )
            );
    return desc;
}

template <typename real>
hipsparseDnVecDescr_t hip_vector_description(const thrust::device_vector<real> &&x) {
    hipsparseDnVecDescr_t desc;
    AMGCL_CALL_HIP(
            hipsparseCreateDnVec(
                &desc,
                x.size(),
                thrust::raw_pointer_cast(&x[0]),
                hip_datatype<real>()
                )
            );
    return desc;
}

template <typename real>
hipsparseSpMatDescr_t hip_matrix_description(
        size_t nrows,
        size_t ncols,
        size_t nnz,
        thrust::device_vector<int> &ptr,
        thrust::device_vector<int> &col,
        thrust::device_vector<real> &val
        )
{
    hipsparseSpMatDescr_t desc;
    AMGCL_CALL_HIP(
            hipsparseCreateCsr(
                &desc,
                nrows,
                ncols,
                nnz,
                thrust::raw_pointer_cast(&ptr[0]),
                thrust::raw_pointer_cast(&col[0]),
                thrust::raw_pointer_cast(&val[0]),
                HIPSPARSE_INDEX_32I,
                HIPSPARSE_INDEX_32I,
                HIPSPARSE_INDEX_BASE_ZERO,
                detail::hip_datatype<real>()
                )
            );
    return desc;
}

} // namespace detail

/// hipSPARSE matrix in CSR format.
template <typename real>
class hip_matrix {
    public:
        typedef real value_type;

        hip_matrix(
                size_t n, size_t m,
                const ptrdiff_t *p_ptr,
                const ptrdiff_t *p_col,
                const real      *p_val,
                hipsparseHandle_t handle
                )
            : nrows(n), ncols(m), nnz(p_ptr[n]), handle(handle),
              ptr(p_ptr, p_ptr + n + 1), col(p_col, p_col + nnz), val(p_val, p_val + nnz)
        {
              desc.reset(
                      detail::hip_matrix_description(nrows, ncols, nnz, ptr, col, val),
                      backend::detail::hip_deleter()
                      );
        }

        void spmv(
                real alpha, thrust::device_vector<real> const &x,
                real beta,  thrust::device_vector<real>       &y
            ) const
        {
            std::shared_ptr<std::remove_pointer<hipsparseDnVecDescr_t>::type> xdesc(
                    detail::hip_vector_description(const_cast<thrust::device_vector<real>&>(x)),
                    backend::detail::hip_deleter()
                    );
            std::shared_ptr<std::remove_pointer<hipsparseDnVecDescr_t>::type> ydesc(
                    detail::hip_vector_description(y),
                    backend::detail::hip_deleter()
                    );

            if (!ready_for_spmv) {
                size_t buf_size;

                AMGCL_CALL_HIP(
                        hipsparseSpMV_bufferSize(
                            handle,
                            HIPSPARSE_OPERATION_NON_TRANSPOSE,
                            &alpha,
                            desc.get(),
                            xdesc.get(),
                            &beta,
                            ydesc.get(),
                            detail::hip_datatype<real>(),
                            HIPSPARSE_SPMV_CSR_ALG1,
                            &buf_size
                            )
                        );

                if (buf.size() < buf_size)
                    buf.resize(buf_size);

                AMGCL_CALL_HIP(
                        hipsparseSpMV_preprocess(
                            handle,
                            HIPSPARSE_OPERATION_NON_TRANSPOSE,
                            &alpha,
                            desc.get(),
                            xdesc.get(),
                            &beta,
                            ydesc.get(),
                            detail::hip_datatype<real>(),
                            HIPSPARSE_SPMV_CSR_ALG1,
                            thrust::raw_pointer_cast(&buf[0])
                            )
                        );

                ready_for_spmv = true;
            }

            AMGCL_CALL_HIP(
                    hipsparseSpMV(
                        handle,
                        HIPSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        desc.get(),
                        xdesc.get(),
                        &beta,
                        ydesc.get(),
                        detail::hip_datatype<real>(),
                        HIPSPARSE_SPMV_CSR_ALG1,
                        thrust::raw_pointer_cast(&buf[0])
                        )
                    );
        }

        size_t rows()     const { return nrows; }
        size_t cols()     const { return ncols; }
        size_t nonzeros() const { return nnz;   }
        size_t bytes()    const {
            return
                sizeof(int)  * (nrows + 1) +
                sizeof(int)  * nnz +
                sizeof(real) * nnz;
        }
    private:
        size_t nrows, ncols, nnz;

        hipsparseHandle_t handle;

        std::shared_ptr<std::remove_pointer<hipsparseSpMatDescr_t>::type> desc;

        thrust::device_vector<int>  ptr;
        thrust::device_vector<int>  col;
        thrust::device_vector<real> val;

        mutable thrust::device_vector<char> buf;
        mutable bool ready_for_spmv = false;

};

/// ROCm/HIP backend.
/**
 * Uses hipSPARSE for matrix operations and Thrust (rocThrust) for vector
 * operations.
 *
 * \param real Value type.
 * \ingroup backends
 */
template <typename real, class DirectSolver = solver::hip_skyline_lu<real> >
struct hip {
        static_assert(
                std::is_same<real, float>::value ||
                std::is_same<real, double>::value,
                "Unsupported value type for hip backend"
                );

    typedef real value_type;
    typedef ptrdiff_t col_type;
    typedef ptrdiff_t ptr_type;
    typedef hip_matrix<real>        matrix;
    typedef thrust::device_vector<real> vector;
    typedef thrust::device_vector<real> matrix_diagonal;
    typedef DirectSolver                direct_solver;

    struct provides_row_iterator : std::false_type {};

    /// Backend parameters.
    struct params {
        /// hipSPARSE handle.
        hipsparseHandle_t hipsparse_handle;

        params(hipsparseHandle_t handle = 0) : hipsparse_handle(handle) {}

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, hipsparse_handle)
        {
            check_params(p, {"hipsparse_handle"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, hipsparse_handle);
        }
#endif
    };

    static std::string name() { return "hip"; }

    /// Copy matrix from builtin backend.
    static std::shared_ptr<matrix>
    copy_matrix(std::shared_ptr< typename builtin<real>::matrix > A, const params &prm)
    {
        return std::make_shared<matrix>(rows(*A), cols(*A),
                A->ptr, A->col, A->val, prm.hipsparse_handle
                );
    }

    /// Copy vector from builtin backend.
    static std::shared_ptr<vector>
    copy_vector(typename builtin<real>::vector const &x, const params&)
    {
        return std::make_shared<vector>(x.data(), x.data() + x.size());
    }

    /// Copy vector from builtin backend.
    static std::shared_ptr<vector>
    copy_vector(std::shared_ptr< typename builtin<real>::vector > x, const params &prm)
    {
        return copy_vector(*x, prm);
    }

    /// Create vector of the specified size.
    static std::shared_ptr<vector>
    create_vector(size_t size, const params&)
    {
        return std::make_shared<vector>(size);
    }

    /// Create direct solver for coarse level
    static std::shared_ptr<direct_solver>
    create_solver(std::shared_ptr< typename builtin<real>::matrix > A, const params &prm)
    {
        return std::make_shared<direct_solver>(A, prm);
    }

    struct gather {
        thrust::device_vector<ptrdiff_t>  I;
        mutable thrust::device_vector<value_type> T;

        gather(size_t src_size, const std::vector<ptrdiff_t> &I, const params&)
            : I(I), T(I.size())
        { }

        void operator()(const vector &src, vector &dst) const {
            thrust::gather(I.begin(), I.end(), src.begin(), dst.begin());
        }

        void operator()(const vector &vec, std::vector<value_type> &vals) const {
            thrust::gather(I.begin(), I.end(), vec.begin(), T.begin());
            thrust::copy(T.begin(), T.end(), vals.begin());
        }
    };

    struct scatter {
        thrust::device_vector<ptrdiff_t>  I;

        scatter(size_t size, const std::vector<ptrdiff_t> &I, const params &)
            : I(I)
        { }

        void operator()(const vector &src, vector &dst) const {
            thrust::scatter(src.begin(), src.end(), I.begin(), dst.begin());
        }
    };
};

//---------------------------------------------------------------------------
// Backend interface implementation
//---------------------------------------------------------------------------
template < typename V >
struct bytes_impl< thrust::device_vector<V> > {
    static size_t get(const thrust::device_vector<V> &v) {
        return v.size() * sizeof(V);
    }
};

template < typename Alpha, typename Beta, typename V >
struct spmv_impl<
    Alpha, hip_matrix<V>, thrust::device_vector<V>,
    Beta,  thrust::device_vector<V>
    >
{
    typedef hip_matrix<V> matrix;
    typedef thrust::device_vector<V> vector;

    static void apply(Alpha alpha, const matrix &A, const vector &x,
            Beta beta, vector &y)
    {
        A.spmv(alpha, x, beta, y);
    }
};

template < typename V >
struct residual_impl<
    hip_matrix<V>,
    thrust::device_vector<V>,
    thrust::device_vector<V>,
    thrust::device_vector<V>
    >
{
    typedef hip_matrix<V> matrix;
    typedef thrust::device_vector<V> vector;

    static void apply(const vector &rhs, const matrix &A, const vector &x,
            vector &r)
    {
        thrust::copy(rhs.begin(), rhs.end(), r.begin());
        A.spmv(-1, x, 1, r);
    }
};

template < typename V >
struct clear_impl< thrust::device_vector<V> >
{
    typedef thrust::device_vector<V> vector;

    static void apply(vector &x)
    {
        thrust::fill(x.begin(), x.end(), V());
    }
};

template <class V, class T>
struct copy_impl<V, thrust::device_vector<T> >
{
    static void apply(const V &x, thrust::device_vector<T> &y)
    {
        thrust::copy(x.begin(), x.end(), y.begin());
    }
};

template <class T, class V>
struct copy_impl<thrust::device_vector<T>, V >
{
    static void apply(const thrust::device_vector<T> &x, V &y)
    {
        thrust::copy(x.begin(), x.end(), y.begin());
    }
};

template <class T1, class T2>
struct copy_impl<thrust::device_vector<T1>, thrust::device_vector<T2> >
{
    static void apply(const thrust::device_vector<T1> &x, thrust::device_vector<T2> &y)
    {
        thrust::copy(x.begin(), x.end(), y.begin());
    }
};

template < typename V >
struct inner_product_impl<
    thrust::device_vector<V>,
    thrust::device_vector<V>
    >
{
    typedef thrust::device_vector<V> vector;

    static V get(const vector &x, const vector &y)
    {
        return thrust::inner_product(x.begin(), x.end(), y.begin(), V());
    }
};

template < typename A, typename B, typename V >
struct axpby_impl<
    A, thrust::device_vector<V>,
    B, thrust::device_vector<V>
    >
{
    typedef thrust::device_vector<V> vector;

    struct functor {
        A a;
        B b;
        functor(A a, B b) : a(a), b(b) {}

        template <class Tuple>
        __host__ __device__ void operator()( Tuple t ) const {
            using thrust::get;

            if (b)
                get<1>(t) = a * get<0>(t) + b * get<1>(t);
            else
                get<1>(t) = a * get<0>(t);
        }
    };

    static void apply(A a, const vector &x, B b, vector &y)
    {
        thrust::for_each(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        x.begin(), y.begin()
                        )
                    ),
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        x.end(), y.end()
                        )
                    ),
                functor(a, b)
                );
    }
};

template < typename A, typename B, typename C, typename V >
struct axpbypcz_impl<
    A, thrust::device_vector<V>,
    B, thrust::device_vector<V>,
    C, thrust::device_vector<V>
    >
{
    typedef thrust::device_vector<V> vector;

    struct functor {
        A a;
        B b;
        C c;

        functor(A a, B b, C c) : a(a), b(b), c(c) {}

        template <class Tuple>
        __host__ __device__ void operator()( Tuple t ) const {
            using thrust::get;

            if (c)
                get<2>(t) = a * get<0>(t) + b * get<1>(t) + c * get<2>(t);
            else
                get<2>(t) = a * get<0>(t) + b * get<1>(t);
        }
    };

    static void apply(
            A a, const vector &x,
            B b, const vector &y,
            C c,       vector &z
            )
    {
        thrust::for_each(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        x.begin(), y.begin(), z.begin()
                        )
                    ),
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        x.end(), y.end(), z.end()
                        )
                    ),
                functor(a, b, c)
                );
    }
};

template < typename A, typename B, typename V >
struct vmul_impl<
    A, thrust::device_vector<V>, thrust::device_vector<V>,
    B, thrust::device_vector<V>
    >
{
    typedef thrust::device_vector<V> vector;

    struct functor {
        A a;
        B b;
        functor(A a, B b) : a(a), b(b) {}

        template <class Tuple>
        __host__ __device__ void operator()( Tuple t ) const {
            using thrust::get;

            if (b)
                get<2>(t) = a * get<0>(t) * get<1>(t) + b * get<2>(t);
            else
                get<2>(t) = a * get<0>(t) * get<1>(t);
        }
    };

    static void apply(A a, const vector &x, const vector &y, B b, vector &z)
    {
        thrust::for_each(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        x.begin(), y.begin(), z.begin()
                        )
                    ),
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        x.end(), y.end(), z.end()
                        )
                    ),
                functor(a, b)
                );
    }
};

class hip_event {
    public:
        hip_event() : e(create_event(), backend::detail::hip_deleter()) { }

        float operator-(hip_event tic) const {
            float delta;
            AMGCL_CALL_HIP( hipEventSynchronize(e.get()) );
            AMGCL_CALL_HIP( hipEventElapsedTime(&delta, tic.e.get(), e.get()) );
            return delta / 1000.0f;
        }
    private:
        std::shared_ptr<std::remove_pointer<hipEvent_t>::type> e;

        static hipEvent_t create_event() {
            hipEvent_t e;
            AMGCL_CALL_HIP( hipEventCreate(&e) );
            AMGCL_CALL_HIP( hipEventRecord(e, 0) );
            return e;
        }
};

struct hip_clock {
    typedef hip_event value_type;

    static const char* units() { return "s"; }

    hip_event current() const {
        return hip_event();
    }
};

} // namespace backend
} // namespace amgcl

#endif
