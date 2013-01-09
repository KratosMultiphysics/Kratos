#ifndef VEXCL_SPAI_HPP
#define VEXCL_SPAI_HPP

/*
The MIT License

Copyright (c) 2012 Denis Demidov <ddemidov@ksu.ru>

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
 * \file   spai.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  Compute Sparse Approximate Inverse for a given matrix.
 */

#include <amgcl/spmat.hpp>

namespace amgcl {

/// Routines computing Sparse Approximate Inverses for a given matrix.
namespace spai {

/// SPAI-0 algorithm from \ref spai_2002 "Broeker (2002)".
template <class spmat>
std::vector<typename sparse::matrix_value<spmat>::type>
level0(const spmat &A) {
    typedef typename sparse::matrix_value<spmat>::type value_t;
    typedef typename sparse::matrix_index<spmat>::type index_t;

    const index_t n = sparse::matrix_rows(A);

    BOOST_AUTO(Arow, sparse::matrix_outer_index(A));
    BOOST_AUTO(Acol, sparse::matrix_inner_index(A));
    BOOST_AUTO(Aval, sparse::matrix_values(A));

    std::vector<value_t> m(n);

#pragma omp parallel for schedule(dynamic, 1024)
    for(index_t i = 0; i < n; ++i) {
        value_t num = 0;
        value_t den = 0;

        for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j) {
            value_t v = Aval[j];
            den += v * v;
            if (Acol[j] == i) num += v;
        }

        m[i] = num / den;
    }

    return m;
}

} // namespace spai
} // namespace amgcl

#endif
