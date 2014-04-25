#ifndef AMGCL_OPERATIONS_VIENNACL_HPP
#define AMGCL_OPERATIONS_VIENNACL_HPP

/*
The MIT License

Copyright (c) 2012-2014 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   operations_viennacl.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Adaptors for ViennaCL types.
 */

#include <amgcl/spmat.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/linalg/prod.hpp>
#include "viennacl/linalg/inner_prod.hpp"
#include <viennacl/linalg/norm_2.hpp>
#include <viennacl/traits/clear.hpp>

namespace amgcl {

template <typename T>
struct value_type< viennacl::vector<T> > {
    typedef T type;
};

template <typename T>
void clear(viennacl::vector<T> &x) {
    viennacl::traits::clear(x);
}

template <typename T>
T inner_prod(const viennacl::vector<T> &x, const viennacl::vector<T> &y) {
    return viennacl::linalg::inner_prod(x, y);
}

template <typename T>
T norm(const viennacl::vector<T> &x) {
    return viennacl::linalg::norm_2(x);
}

template <class matrix, typename real>
void residual(
        const matrix &A,
        const viennacl::vector<real> &x,
        const viennacl::vector<real> &f,
        viennacl::vector<real> &y
        )
{
    y = viennacl::linalg::prod(A, x);
    y = f - y;
}

/// Specialization of matrix-vector product for ublas types.
/** Necessary for ublas types to work with amgcl::solve() functions. */
template <class matrix, typename real>
void axpy(
        const matrix &A,
        const viennacl::vector<real> &x,
        viennacl::vector<real> &y
        )
{
    y = viennacl::linalg::prod(A, x);
}

namespace sparse {

/// Provides proxy matrix class that may be used with ViennaCL.
/** No data is copied here. The proxy class does not own the data.  */
template <class spmat>
class viennacl_matrix_adapter {
    public:
        typedef typename sparse::matrix_index<spmat>::type index_type;
        typedef typename sparse::matrix_value<spmat>::type value_type;

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

                index_type row;
                index_type pos;
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

        viennacl_matrix_adapter(const spmat &A)
            : rows(sparse::matrix_rows(A)),
              cols(sparse::matrix_cols(A)),
              row(sparse::matrix_outer_index(A)),
              col(sparse::matrix_inner_index(A)),
              val(sparse::matrix_values(A))
        { }

        const_iterator1 begin1() const {
            return const_iterator1(0, row, col, val);
        }

        const_iterator1 end1() const {
            return const_iterator1(rows, row, col, val);
        }

        index_type size1() const {
            return rows;
        }

        index_type size2() const {
            return cols;
        }
    private:
        index_type rows;
        index_type cols;

        const index_type *row;
        const index_type *col;
        const value_type *val;
};

template <class spmat>
viennacl_matrix_adapter<spmat> viennacl_map(const spmat &A) {
    return viennacl_matrix_adapter<spmat>(A);
}

} // namespace sparse

/// Wrapper around amgcl::solver that is compatible with ViennaCL solvers.
/**
 * \param vector Vector type that will be used with the preconditioner.
 * \param AMG    Type of amgcl::solver
 */
template <class vector, class AMG>
class viennacl_preconditioner {
    public:
        viennacl_preconditioner(const AMG &amg)
            : amg(amg), buf(amg.size()) {}

        void apply(vector &x) const {
            buf.swap(x);
            viennacl::traits::clear(x);
            amg.apply(buf, x);
        }
    private:
        const AMG &amg;
        mutable vector buf;
};

/// Wrapper around amgcl::solver that is compatible with ViennaCL solvers.
template <class vector, class AMG>
viennacl_preconditioner<vector, AMG> make_viennacl_precond(const AMG &amg) {
    return viennacl_preconditioner<vector, AMG>(amg);
}

} // namespace amgcl

#endif
