#ifndef AMGCL_ADAPTER_CRS_BUILDER_HPP
#define AMGCL_ADAPTER_CRS_BUILDER_HPP

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
\file    amgcl/adapter/crs_builder.hpp
\author  Denis Demidov <dennis.demidov@gmail.com>
\brief   Matrix builder that creates matrix rows as needed.
\ingroup adapters

Example:
\code
struct poisson_2d {
    typedef double val_type;
    typedef long   col_type;

    poisson_2d(size_t n) : n(n), h2i((n - 1) * (n - 1)) {}

    // Number of rows in the constructed matrix:
    size_t rows() const { return n * n; }

    // Estimated number of nonzeros in the problem:
    size_t nonzeros() const { return 5 * rows(); }

    // Fills column numbers and values of nonzero elements in the given matrix row.
    void operator()(size_t row,
            std::vector<col_type> &col,
            std::vector<val_type> &val
            ) const
    {
        size_t i = row % n;
        size_t j = row / n;

        if (j > 0) {
            col.push_back(row - n);
            val.push_back(-h2i);
        }

        if (i > 0) {
            col.push_back(row - 1);
            val.push_back(-h2i);
        }

        col.push_back(row);
        val.push_back(4 * h2i);

        if (i + 1 < n) {
            col.push_back(row + 1);
            val.push_back(-h2i);
        }

        if (j + 1 < n) {
            col.push_back(row + n);
            val.push_back(-h2i);
        }
    }

    private:
        size_t n;
        double h2i;
};

amgcl::make_solver<
    Backend, Coarsening, Relaxation, IterativeSolver
    > solve( amgcl::backend::make_matrix( poisson_2d(m) ) );
\endcode
*/

#include <amgcl/backend/interface.hpp>
#include <amgcl/backend/detail/matrix_ops.hpp>

namespace amgcl {

/// Matrix adapters.
namespace adapter {

/// Generates matrix rows as needed with help of user-provided functor.
/**
 * The generated rows are not stored anywhere.
 */
template <class RowBuilder>
struct matrix_builder {
    typedef typename RowBuilder::col_type col_type;
    typedef typename RowBuilder::val_type val_type;

    RowBuilder build_row;

    mutable std::vector<col_type> col;
    mutable std::vector<val_type> val;

    matrix_builder(const RowBuilder &row_builder) : build_row(row_builder)
    {}

    size_t rows() const { return build_row.rows(); }
    size_t cols() const { return build_row.rows(); }

    size_t nonzeros() const { return build_row.nonzeros(); }

    struct row_iterator {
        typedef typename RowBuilder::col_type col_type;
        typedef typename RowBuilder::val_type val_type;

        typedef typename std::vector<col_type>::const_iterator col_iterator;
        typedef typename std::vector<val_type>::const_iterator val_iterator;

        row_iterator(
                col_iterator col_begin,
                col_iterator col_end,
                val_iterator val_begin
                )
            : m_col(col_begin), m_end(col_end), m_val(val_begin)
        {}

        operator bool() const {
            return m_col != m_end;
        }

        row_iterator& operator++() {
            ++m_col;
            ++m_val;
            return *this;
        }

        col_type col() const {
            return *m_col;
        }

        val_type value() const {
            return *m_val;
        }

        private:
            col_iterator m_col;
            col_iterator m_end;
            val_iterator m_val;
    };

    row_iterator row_begin(size_t i) const {
        col.clear();
        val.clear();
        build_row(i, col, val);

        return row_iterator(col.begin(), col.end(), val.begin());
    }

};

/// Convenience function returning an instance of matrix_builder<RowBuilder>
template <class RowBuilder>
matrix_builder<RowBuilder> make_matrix(const RowBuilder &row_builder) {
    return matrix_builder<RowBuilder>(row_builder);
}

} // namespace adapter

namespace backend {

//---------------------------------------------------------------------------
// Specialization of matrix interface
//---------------------------------------------------------------------------
template <class RowBuilder>
struct value_type< adapter::matrix_builder<RowBuilder> >
{
    typedef typename adapter::matrix_builder<RowBuilder>::val_type type;
};

template <class RowBuilder>
struct rows_impl< adapter::matrix_builder<RowBuilder> >
{
    static size_t get(const adapter::matrix_builder<RowBuilder> &A) {
        return A.rows();
    }
};

template <class RowBuilder>
struct cols_impl< adapter::matrix_builder<RowBuilder> >
{
    static size_t get(const adapter::matrix_builder<RowBuilder> &A) {
        return A.cols();
    }
};

template <class RowBuilder>
struct nonzeros_impl< adapter::matrix_builder<RowBuilder> >
{
    static size_t get(const adapter::matrix_builder<RowBuilder> &A) {
        return A.nonzeros();
    }
};

template <class RowBuilder>
struct row_iterator< adapter::matrix_builder<RowBuilder> >
{
    typedef typename adapter::matrix_builder<RowBuilder>::row_iterator type;
};

template <class RowBuilder>
struct row_begin_impl< adapter::matrix_builder<RowBuilder> >
{
    typedef adapter::matrix_builder<RowBuilder> Matrix;
    static typename row_iterator<Matrix>::type
    get(const Matrix &matrix, size_t row) {
        return matrix.row_begin(row);
    }
};

namespace detail {

template <class RowBuilder>
struct use_builtin_matrix_ops< amgcl::adapter::matrix_builder<RowBuilder> >
    : std::true_type
{};

} // namespace detail

} // namespace backend
} // namespace amgcl

#endif
