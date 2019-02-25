#ifndef AMGCL_ADAPTER_COMPLEX_HPP
#define AMGCL_ADAPTER_COMPLEX_HPP

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
\file    amgcl/adapter/complex.hpp
\author  Denis Demidov <dennis.demidov@gmail.com>
\brief   Complex-valued matrix adapter.
\ingroup adapters
*/

#include <type_traits>
#include <boost/range/iterator_range.hpp>

#include <amgcl/backend/interface.hpp>
#include <amgcl/backend/detail/matrix_ops.hpp>

namespace amgcl {
namespace adapter {

template <class Matrix>
struct complex_adapter {
    static_assert(is_complex<typename backend::value_type<Matrix>::type>::value,
            "value type should be complex");

    typedef typename backend::value_type<Matrix>::type::value_type val_type;

    const Matrix &A;

    complex_adapter(const Matrix &A) : A(A) {}

    size_t rows() const {
        return 2 * backend::rows(A);
    }

    size_t cols() const {
        return 2 * backend::cols(A);
    }

    size_t nonzeros() const {
        return 4 * backend::nonzeros(A);
    }

    struct row_iterator {
        typedef typename backend::row_iterator<Matrix>::type Base;
        typedef typename Base::col_type col_type;

        row_iterator(const Base &base, bool row_real)
            : base(base), row_real(row_real), col_real(true) {}

        operator bool() const {
            return static_cast<bool>(base);
        }

        row_iterator& operator++() {
            col_real = !col_real;
            if (col_real) ++base;

            return *this;
        }

        col_type col() const {
            if (col_real)
                return base.col() * 2;
            else
                return base.col() * 2 + 1;
        }

        val_type value() const {
            if (row_real) {
                if (col_real)
                    return std::real(base.value());
                else
                    return -std::imag(base.value());
            } else {
                if (col_real)
                    return std::imag(base.value());
                else
                    return std::real(base.value());
            }
        }

        private:
            Base base;
            bool row_real;
            bool col_real;

    };

    row_iterator row_begin(size_t i) const {
        return row_iterator(backend::row_begin(A, i / 2), i % 2 == 0);
    }
};

template <class Matrix>
complex_adapter<Matrix> complex_matrix(const Matrix &A) {
    return complex_adapter<Matrix>(A);
}

template <class Range>
boost::iterator_range<
    typename std::add_pointer<
        typename std::conditional<
            std::is_const<Range>::value,
            typename std::add_const<
                typename boost::range_value<
                    typename std::decay<Range>::type
                    >::type::value_type
                >::type,
            typename boost::range_value<
                typename std::decay<Range>::type
                >::type::value_type
            >::type
        >::type
    >
complex_range(Range &rng) {
    typedef
        typename std::add_pointer<
            typename std::conditional<
                std::is_const<Range>::value,
                typename std::add_const<
                    typename boost::range_value<
                        typename std::decay<Range>::type
                        >::type::value_type
                    >::type,
                typename boost::range_value<
                    typename std::decay<Range>::type
                    >::type::value_type
                >::type
            >::type
        pointer_type;

    pointer_type b = reinterpret_cast<pointer_type>(&rng[0]);
    pointer_type e = b + 2 * boost::size(rng);

    return boost::iterator_range<pointer_type>(b, e);
}

} // namespace adapter

namespace backend {

//---------------------------------------------------------------------------
// Specialization of matrix interface
//---------------------------------------------------------------------------
template <class Matrix>
struct value_type< adapter::complex_adapter<Matrix> >
{
    typedef typename adapter::complex_adapter<Matrix>::val_type type;
};

template <class Matrix>
struct rows_impl< adapter::complex_adapter<Matrix> >
{
    static size_t get(const adapter::complex_adapter<Matrix> &A) {
        return A.rows();
    }
};

template <class Matrix>
struct cols_impl< adapter::complex_adapter<Matrix> >
{
    static size_t get(const adapter::complex_adapter<Matrix> &A) {
        return A.cols();
    }
};

template <class Matrix>
struct nonzeros_impl< adapter::complex_adapter<Matrix> >
{
    static size_t get(const adapter::complex_adapter<Matrix> &A) {
        return A.nonzeros();
    }
};

template <class Matrix>
struct row_iterator< adapter::complex_adapter<Matrix> >
{
    typedef typename adapter::complex_adapter<Matrix>::row_iterator type;
};

template <class Matrix>
struct row_begin_impl< adapter::complex_adapter<Matrix> >
{
    typedef adapter::complex_adapter<Matrix> CM;
    static typename row_iterator<CM>::type
    get(const CM &matrix, size_t row) {
        return matrix.row_begin(row);
    }
};

namespace detail {

template <class Matrix>
struct use_builtin_matrix_ops< amgcl::adapter::complex_adapter<Matrix> >
    : std::true_type
{};

} // namespace detail
} // namespace backend
} // namespace amgcl

#endif
