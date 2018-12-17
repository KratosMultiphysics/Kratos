#ifndef AMGCL_ADAPTER_UBLAS_HPP
#define AMGCL_ADAPTER_UBLAS_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>
Copyright (c) 2014, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)

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
\file    amgcl/adapter/ublas.hpp
\author  Denis Demidov <dennis.demidov@gmail.com>
\brief   Adapters for Boost.uBlas matrices and vectors.
\ingroup adapters
*/

#include <boost/range/iterator_range.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

namespace amgcl {
namespace backend {

// Make builtin backend recognize ublas vectors as its own:
template <typename T>
struct is_builtin_vector< boost::numeric::ublas::vector<T> >
    : std::true_type {};

/// Adapts Boost.uBlas matrix.
template <typename T>
std::tuple<
    size_t,
    boost::iterator_range<const size_t*>,
    boost::iterator_range<const size_t*>,
    boost::iterator_range<const T*>
    >
map(const boost::numeric::ublas::compressed_matrix<T, boost::numeric::ublas::row_major> &A) {
    return std::make_tuple(
            A.size1(),
            boost::make_iterator_range(
                A.index1_data().begin(), A.index1_data().end()
                ),
            boost::make_iterator_range(
                A.index2_data().begin(), A.index2_data().end()
                ),
            boost::make_iterator_range(
                A.value_data().begin(), A.value_data().end()
                )
            );
}

} // namespace backend
} // namespace amgcl

#endif
