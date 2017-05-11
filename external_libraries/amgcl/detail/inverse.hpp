#ifndef AMGCL_DETAIL_INVERSE_HPP
#define AMGCL_DETAIL_INVERSE_HPP

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
 * \file   amgcl/detail/inverse.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Inverse of a dense matrix.
 */

#include <vector>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <amgcl/util.hpp>

namespace amgcl {
namespace detail {

    template <typename value_type>
    static void inverse(size_t n, value_type *data) {
        using namespace boost::numeric::ublas;

        matrix<value_type> A(n, n);
        std::copy(data, data + n * n, &A.data()[0]);

        permutation_matrix<size_t> pm(n);

        lu_factorize(A, pm);

        matrix<value_type> I(n, n);
        I.assign( identity_matrix<value_type>(n) );

        lu_substitute(A, pm, I);

        std::copy(&I.data()[0], &I.data()[0] + n * n, data);
    }

} // namespace detail
} // namespace amgcl

#endif
