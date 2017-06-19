#ifndef AMGCL_BACKEND_CUDA_HPP
#define AMGCL_BACKEND_CUDA_HPP

/*
The MIT License

Copyright (c) 2012-2017 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/backend/cuda.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  CUDA backend.
 */

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/solver/skyline_lu.hpp>
#include <amgcl/util.hpp>

#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/for_each.h>
#include <thrust/inner_product.h>
#include <cusparse_v2.h>

namespace amgcl {

namespace solver {

/** Wrapper around solver::skyline_lu for use with the CUDA backend.
 * Copies the rhs to the host memory, solves the problem using the host CPU,
 * then copies the solution back to the compute device(s).
 */
template <class T>
struct cuda_skyline_lu : solver::skyline_lu<T> {
    typedef solver::skyline_lu<T> Base;

    mutable std::vector<T> _rhs, _x;

    template <class Matrix, class Params>
    cuda_skyline_lu(const Matrix &A, const Params&)
        : Base(*A), _rhs(backend::rows(*A)), _x(backend::rows(*A))
    { }

    template <class Vec1, class Vec2>
    void operator()(const Vec1 &rhs, Vec2 &x) const {
        thrust::copy(rhs.begin(), rhs.end(), _rhs.begin());
        static_cast<const Base*>(this)->operator()(_rhs, _x);
        thrust::copy(_x.begin(), _x.end(), x.begin());
    }
};

}

namespace backend {
namespace detail {

inline void cuda_check(cusparseStatus_t rc, const char *file, int line) {
    if (rc != CUSPARSE_STATUS_SUCCESS) {
        std::ostringstream msg;
        msg << "CUDA error " << rc << " at \"" << file << ":" << line;
        precondition(false, msg.str());
    }
}

#define AMGCL_CALL_CUDA(rc)                                                    \
    amgcl::backend::detail::cuda_check(rc, __FILE__, __LINE__)

struct cuda_deleter {
    void operator()(cusparseMatDescr_t handle) {
        AMGCL_CALL_CUDA( cusparseDestroyMatDescr(handle) );
    }

    void operator()(cusparseHybMat_t handle) {
        AMGCL_CALL_CUDA( cusparseDestroyHybMat(handle) );
    }

    void operator()(csrilu02Info_t handle) {
        AMGCL_CALL_CUDA( cusparseDestroyCsrilu02Info(handle) );
    }

    void operator()(csrsv2Info_t handle) {
        AMGCL_CALL_CUDA( cusparseDestroyCsrsv2Info(handle) );
    }
};


} // namespace detail

/// CUSPARSE matrix in Hyb format.
template <typename real>
class cuda_hyb_matrix {
    public:
        typedef real value_type;

        cuda_hyb_matrix(
                size_t n, size_t m,
                const ptrdiff_t *ptr,
                const ptrdiff_t *col,
                const real      *val,
                cusparseHandle_t handle
                )
            : nrows(n), ncols(m), nnz(ptr[n]), handle( handle ),
              desc  ( create_description(), backend::detail::cuda_deleter() ),
              mat   ( create_matrix(),      backend::detail::cuda_deleter() )
        {
            fill_matrix(n, m, ptr, col, val);
        }

        void spmv(
                real alpha, thrust::device_vector<real> const &x,
                real beta,  thrust::device_vector<real>       &y,
                boost::false_type
            ) const
        {
            AMGCL_CALL_CUDA(
                    cusparseShybmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha, desc.get(), mat.get(),
                        thrust::raw_pointer_cast(&x[0]), &beta,
                        thrust::raw_pointer_cast(&y[0])
                        )
                    );
        }

        void spmv(
                real alpha, thrust::device_vector<real> const &x,
                real beta,  thrust::device_vector<real>       &y,
                boost::true_type
            ) const
        {
            AMGCL_CALL_CUDA(
                    cusparseDhybmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha, desc.get(), mat.get(),
                        thrust::raw_pointer_cast(&x[0]), &beta,
                        thrust::raw_pointer_cast(&y[0])
                        )
                    );
        }

        size_t rows()     const { return nrows; }
        size_t cols()     const { return ncols; }
        size_t nonzeros() const { return nnz;   }
    private:
        size_t nrows, ncols, nnz;

        cusparseHandle_t handle;

        boost::shared_ptr<boost::remove_pointer<cusparseMatDescr_t>::type> desc;
        boost::shared_ptr<boost::remove_pointer<cusparseHybMat_t>::type>   mat;

        static cusparseMatDescr_t create_description() {
            cusparseMatDescr_t desc;
            AMGCL_CALL_CUDA( cusparseCreateMatDescr(&desc) );
            AMGCL_CALL_CUDA( cusparseSetMatType(desc, CUSPARSE_MATRIX_TYPE_GENERAL) );
            AMGCL_CALL_CUDA( cusparseSetMatIndexBase(desc, CUSPARSE_INDEX_BASE_ZERO) );
            return desc;
        }

        static cusparseHybMat_t create_matrix() {
            cusparseHybMat_t mat;
            AMGCL_CALL_CUDA( cusparseCreateHybMat(&mat) );
            return mat;
        }

        void fill_matrix(size_t n, size_t m,
                const ptrdiff_t *ptr, const ptrdiff_t *col, const float *val
                )
        {
            thrust::device_vector<int>   p(ptr, ptr + n + 1);
            thrust::device_vector<int>   c(col, col + ptr[n]);
            thrust::device_vector<float> v(val, val + ptr[n]);

            AMGCL_CALL_CUDA(
                    cusparseScsr2hyb(handle, n, m, desc.get(),
                        thrust::raw_pointer_cast(&v[0]),
                        thrust::raw_pointer_cast(&p[0]),
                        thrust::raw_pointer_cast(&c[0]),
                        mat.get(), 0, CUSPARSE_HYB_PARTITION_AUTO
                        )
                    );
        }

        void fill_matrix(size_t n, size_t m,
                const ptrdiff_t *ptr, const ptrdiff_t *col, const double *val
                )
        {
            thrust::device_vector<int>    p(ptr, ptr + n + 1);
            thrust::device_vector<int>    c(col, col + ptr[n]);
            thrust::device_vector<double> v(val, val + ptr[n]);

            AMGCL_CALL_CUDA(
                    cusparseDcsr2hyb(handle, n, m, desc.get(),
                        thrust::raw_pointer_cast(&v[0]),
                        thrust::raw_pointer_cast(&p[0]),
                        thrust::raw_pointer_cast(&c[0]),
                        mat.get(), 0, CUSPARSE_HYB_PARTITION_AUTO
                        )
                    );
        }
};

/// CUDA backend.
/**
 * Uses CUSPARSE for matrix operations and Thrust for vector operations.
 *
 * \param real Value type.
 * \ingroup backends
 */
template <typename real, class DirectSolver = solver::cuda_skyline_lu<real> >
struct cuda {
        BOOST_STATIC_ASSERT_MSG(
                (
                 boost::is_same<real, float>::value ||
                 boost::is_same<real, double>::value
                ),
                "Unsupported value type for cuda backend"
                );

    typedef real value_type;
    typedef cuda_hyb_matrix<real>       matrix;
    typedef thrust::device_vector<real> vector;
    typedef thrust::device_vector<real> matrix_diagonal;
    typedef DirectSolver                direct_solver;

    struct provides_row_iterator : boost::false_type {};

    /// Backend parameters.
    struct params {
        /// CUSPARSE handle.
        cusparseHandle_t cusparse_handle;

        params(cusparseHandle_t handle = 0) : cusparse_handle(handle) {}

        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, cusparse_handle)
        {
            AMGCL_PARAMS_CHECK(p, (cusparse_handle));
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, cusparse_handle);
        }
    };

    static std::string name() { return "cuda"; }

    /// Copy matrix from builtin backend.
    static boost::shared_ptr<matrix>
    copy_matrix(boost::shared_ptr< typename builtin<real>::matrix > A, const params &prm)
    {
        return boost::make_shared<matrix>(rows(*A), cols(*A),
                A->ptr, A->col, A->val, prm.cusparse_handle
                );
    }

    /// Copy vector from builtin backend.
    static boost::shared_ptr<vector>
    copy_vector(typename builtin<real>::vector const &x, const params&)
    {
        return boost::make_shared<vector>(x.data(), x.data() + x.size());
    }

    /// Copy vector from builtin backend.
    static boost::shared_ptr<vector>
    copy_vector(boost::shared_ptr< typename builtin<real>::vector > x, const params &prm)
    {
        return copy_vector(*x, prm);
    }

    /// Create vector of the specified size.
    static boost::shared_ptr<vector>
    create_vector(size_t size, const params&)
    {
        return boost::make_shared<vector>(size);
    }

    /// Create direct solver for coarse level
    static boost::shared_ptr<direct_solver>
    create_solver(boost::shared_ptr< typename builtin<real>::matrix > A, const params &prm)
    {
        return boost::make_shared<direct_solver>(A, prm);
    }
};

//---------------------------------------------------------------------------
// Backend interface implementation
//---------------------------------------------------------------------------
template < typename V >
struct rows_impl< cuda_hyb_matrix<V> > {
    static size_t get(const cuda_hyb_matrix<V> &A) {
        return A.rows();
    }
};

template < typename V >
struct cols_impl< cuda_hyb_matrix<V> > {
    static size_t get(const cuda_hyb_matrix<V> &A) {
        return A.cols();
    }
};

template < typename V >
struct nonzeros_impl< cuda_hyb_matrix<V> > {
    static size_t get(const cuda_hyb_matrix<V> &A) {
        return A.nonzeros();
    }
};

template < typename Alpha, typename Beta, typename V >
struct spmv_impl<
    Alpha, cuda_hyb_matrix<V>, thrust::device_vector<V>,
    Beta,  thrust::device_vector<V>
    >
{
    typedef cuda_hyb_matrix<V> matrix;
    typedef thrust::device_vector<V> vector;

    static void apply(Alpha alpha, const matrix &A, const vector &x,
            Beta beta, vector &y)
    {
        A.spmv(alpha, x, beta, y, typename boost::is_same<V, double>::type());
    }
};

template < typename V >
struct residual_impl<
    cuda_hyb_matrix<V>,
    thrust::device_vector<V>,
    thrust::device_vector<V>,
    thrust::device_vector<V>
    >
{
    typedef cuda_hyb_matrix<V> matrix;
    typedef thrust::device_vector<V> vector;

    static void apply(const vector &rhs, const matrix &A, const vector &x,
            vector &r)
    {
        thrust::copy(rhs.begin(), rhs.end(), r.begin());
        A.spmv(-1, x, 1, r, typename boost::is_same<V, double>::type());
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

template < typename V >
struct copy_impl<
    thrust::device_vector<V>,
    thrust::device_vector<V>
    >
{
    typedef thrust::device_vector<V> vector;

    static void apply(const vector &x, vector &y)
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

} // namespace backend
} // namespace amgcl

#endif
