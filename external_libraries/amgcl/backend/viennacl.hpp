#ifndef AMGCL_BACKEND_VIENNACL_HPP
#define AMGCL_BACKEND_VIENNACL_HPP

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
 * \file   amgcl/backend/viennacl.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  ViennaCL backend.
 */

#include <type_traits>

#include <viennacl/vector.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/ell_matrix.hpp>
#include <viennacl/hyb_matrix.hpp>
#include <viennacl/linalg/inner_prod.hpp>
#include <viennacl/linalg/prod.hpp>

#include <amgcl/util.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/solver/skyline_lu.hpp>

namespace amgcl {

namespace solver {

/** Wrapper around solver::skyline_lu for use with the ViennaCL backend.
 * Copies the rhs to the host memory, solves the problem using the host CPU,
 * then copies the solution back to the compute device(s).
 */
template <class T>
struct viennacl_skyline_lu : solver::skyline_lu<T> {
    typedef solver::skyline_lu<T> Base;

    mutable std::vector<T> _rhs, _x;

    template <class Matrix, class Params>
    viennacl_skyline_lu(const Matrix &A, const Params&)
        : Base(*A), _rhs(backend::rows(*A)), _x(backend::rows(*A))
    { }

    template <class Vec1, class Vec2>
    void operator()(const Vec1 &rhs, Vec2 &x) const {
        viennacl::fast_copy(rhs, _rhs);
        static_cast<const Base*>(this)->operator()(_rhs, _x);
        viennacl::fast_copy(_x, x);
    }
};

}

namespace backend {

/// ViennaCL backend
/**
 * This is a backend that uses types defined in the ViennaCL library
 * (http://viennacl.sourceforge.net).
 *
 * \param Matrix ViennaCL matrix to use with the backend. Possible choices are
 * viannacl::compressed_matrix<T>, viennacl::ell_matrix<T>, and
 * viennacl::hyb_matrix<T>.
 * \ingroup backends
 */
template <
    class Matrix,
    class DirectSolver = solver::viennacl_skyline_lu<typename backend::value_type<Matrix>::type>
    >
struct viennacl {
    typedef typename backend::value_type<Matrix>::type value_type;
    typedef ptrdiff_t                                  index_type;
    typedef Matrix                                     matrix;
    typedef ::viennacl::vector<value_type>             vector;
    typedef ::viennacl::vector<value_type>             matrix_diagonal;
    typedef DirectSolver                               direct_solver;

    struct provides_row_iterator : std::false_type {};

    /// Backend parameters.
    typedef amgcl::detail::empty_params params;

    static std::string name() { return "viennacl"; }

    /// Copy matrix from builtin backend.
    static std::shared_ptr<matrix>
    copy_matrix(
            std::shared_ptr< typename builtin<value_type>::matrix > A,
            const params&
            )
    {
        auto m = std::make_shared<matrix>();
        ::viennacl::copy(viennacl_matrix_adapter(*A), *m);
        return m;
    }

    /// Copy vector from builtin backend.
    static std::shared_ptr<vector>
    copy_vector(typename builtin<value_type>::vector const &x, const params&)
    {
        auto v = std::make_shared<vector>(x.size());
        ::viennacl::fast_copy(x.data(), x.data() + x.size(), v->begin());
        return v;
    }

    /// Copy vector from builtin backend.
    static std::shared_ptr<vector>
    copy_vector(
            std::shared_ptr< typename builtin<value_type>::vector > x,
            const params &prm
            )
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
    create_solver(std::shared_ptr< typename builtin<value_type>::matrix > A, const params &prm)
    {
        return std::make_shared<direct_solver>(A, prm);
    }

    private:
        struct viennacl_matrix_adapter {
            typedef ptrdiff_t index_type;
            typedef size_t    size_type;

            class const_iterator1;

            class const_iterator2 {
                public:
                    bool operator!=(const const_iterator2 &it) const {
                        return pos != it.pos;
                    }

                    const const_iterator2& operator++() {
                        ++pos;
                        return *this;
                    }

                    index_type index1() const {
                        return row;
                    }

                    index_type index2() const {
                        return col[pos];
                    }

                    value_type operator*() const {
                        return val[pos];
                    }
                private:
                    const_iterator2(index_type row, index_type pos,
                            const index_type *col, const value_type *val)
                        : row(row), pos(pos), col(col), val(val)
                    { }

                    index_type        row;
                    index_type        pos;
                    const index_type *col;
                    const value_type *val;

                    friend class const_iterator1;
            };

            class const_iterator1 {
                public:
                    bool operator!=(const const_iterator1 &it) const {
                        return pos != it.pos;
                    }

                    const const_iterator1& operator++() {
                        ++pos;
                        return *this;
                    }

                    index_type index1() const {
                        return pos;
                    }

                    const const_iterator2 begin() const {
                        return const_iterator2(pos, row[pos], col, val);
                    }

                    const const_iterator2 end() const {
                        return const_iterator2(pos, row[pos + 1], col, val);
                    }
                private:
                    const_iterator1(index_type pos,
                            const index_type *row,
                            const index_type *col,
                            const value_type *val
                            )
                        : pos(pos), row(row), col(col), val(val)
                    { }

                    index_type pos;
                    const index_type *row;
                    const index_type *col;
                    const value_type *val;

                    friend class viennacl_matrix_adapter;
            };

            viennacl_matrix_adapter(
                    const typename backend::builtin<value_type>::matrix &A)
                : rows(A.nrows), cols(A.ncols),
                  row(A.ptr), col(A.col), val(A.val)
            { }

            const_iterator1 begin1() const {
                return const_iterator1(0, row, col, val);
            }

            const_iterator1 end1() const {
                return const_iterator1(rows, row, col, val);
            }

            size_t size1() const {
                return rows;
            }

            size_t size2() const {
                return cols;
            }
            private:
                size_t rows;
                size_t cols;

                const index_type *row;
                const index_type *col;
                const value_type *val;
        };
};

template <class T>
struct is_viennacl_matrix : std::false_type {};

template <class V>
struct is_viennacl_matrix< ::viennacl::compressed_matrix<V> > : std::true_type
{};

template <class V>
struct is_viennacl_matrix< ::viennacl::hyb_matrix<V> > : std::true_type
{};

template <class V>
struct is_viennacl_matrix< ::viennacl::ell_matrix<V> > : std::true_type
{};

template <class M>
struct value_type<
    M,
    typename std::enable_if< is_viennacl_matrix<M>::value >::type
    >
{
    typedef typename M::value_type::value_type type;
};

template <class V>
struct value_type< ::viennacl::vector<V> >
{
    typedef V type;
};

template <class M>
struct rows_impl<
    M,
    typename std::enable_if< is_viennacl_matrix<M>::value >::type
    >
{
    static size_t get(const M &A) {
        return A.size1();
    }
};

template <class M>
struct cols_impl<
    M,
    typename std::enable_if< is_viennacl_matrix<M>::value >::type
    >
{
    static size_t get(const M &A) {
        return A.size2();
    }
};

template <class V>
struct nonzeros_impl< ::viennacl::compressed_matrix<V> > {
    static size_t get(const ::viennacl::compressed_matrix<V> &A) {
        return A.nnz();
    }
};

template <class V>
struct nonzeros_impl< ::viennacl::ell_matrix<V> > {
    static size_t get(const ::viennacl::ell_matrix<V> &A) {
        return A.nnz();
    }
};

template <class V>
struct nonzeros_impl< ::viennacl::hyb_matrix<V> > {
    static size_t get(const ::viennacl::hyb_matrix<V> &A) {
        return A.ell_nnz() + A.csr_nnz();
    }
};

template <class Alpha, class Mtx, class Beta, class Vec>
struct spmv_impl<
    Alpha, Mtx, Vec, Beta, Vec,
    typename std::enable_if< is_viennacl_matrix<Mtx>::value >::type
    >
{
    static void apply(Alpha alpha, const Mtx &A, const Vec &x, Beta beta, Vec &y)
    {
        if (beta)
            y = alpha * ::viennacl::linalg::prod(A, x) + beta * y;
        else
            y = alpha * ::viennacl::linalg::prod(A, x);
    }
};

template <class Mtx, class Vec>
struct residual_impl<
    Mtx, Vec, Vec, Vec,
    typename std::enable_if< is_viennacl_matrix<Mtx>::value >::type
    >
{
    typedef typename value_type<Mtx>::type V;

    static void apply(const Vec &rhs, const Mtx &A, const Vec &x, Vec &r)
    {
        r = ::viennacl::linalg::prod(A, x);
        r = rhs - r;
    }
};

template < typename V >
struct clear_impl< ::viennacl::vector<V> >
{
    static void apply(::viennacl::vector<V> &x)
    {
        x.clear();
    }
};

template < typename V >
struct inner_product_impl<
    ::viennacl::vector<V>,
    ::viennacl::vector<V>
    >
{
    static V get(const ::viennacl::vector<V> &x, const ::viennacl::vector<V> &y)
    {
        return ::viennacl::linalg::inner_prod(x, y);
    }
};

template < typename A, typename B, typename V >
struct axpby_impl<
    A, ::viennacl::vector<V>,
    B, ::viennacl::vector<V>
    >
{
    static void apply(
            A a, const ::viennacl::vector<V> &x,
            B b, ::viennacl::vector<V> &y
            )
    {
        if (b)
            y = a * x + b * y;
        else
            y = a * x;
    }
};

template < typename A, typename B, typename C, typename V >
struct axpbypcz_impl<
    A, ::viennacl::vector<V>,
    B, ::viennacl::vector<V>,
    C, ::viennacl::vector<V>
    >
{
    static void apply(
            A a, const ::viennacl::vector<V> &x,
            B b, const ::viennacl::vector<V> &y,
            C c,       ::viennacl::vector<V> &z
            )
    {
        if (c)
            z = a * x + b * y + c * z;
        else
            z = a * x + b * y;
    }
};

template < typename A, typename B, typename V >
struct vmul_impl<
    A, ::viennacl::vector<V>, ::viennacl::vector<V>,
    B, ::viennacl::vector<V>
    >
{
    static void apply(
            A a, const ::viennacl::vector<V> &x, const ::viennacl::vector<V> &y,
            B b, ::viennacl::vector<V> &z)
    {
        if (b)
            z = a * ::viennacl::linalg::element_prod(x, y) + b * z;
        else
            z = a * ::viennacl::linalg::element_prod(x, y);
    }
};

template < typename V >
struct copy_impl<
    ::viennacl::vector<V>,
    ::viennacl::vector<V>
    >
{
    static void apply(const ::viennacl::vector<V> &x, ::viennacl::vector<V> &y)
    {
        y = x;
    }
};

} // namespace backend
} // namespace amgcl

#endif
