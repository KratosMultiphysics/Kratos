#ifndef AMGCL_DETAIL_SORT_ROW_HPP
#define AMGCL_DETAIL_SORT_ROW_HPP

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
 * \file   amgcl/detail/sort_row.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Sort row of CRS matrix by columns.
 */

namespace amgcl {
namespace detail {

template <typename Col, typename Val>
void sort_row(Col *col, Val *val, int n) {
    for(int j = 1; j < n; ++j) {
        Col c = col[j];
        Val v = val[j];

        int i = j - 1;

        while(i >= 0 && col[i] > c) {
            col[i + 1] = col[i];
            val[i + 1] = val[i];
            i--;
        }

        col[i + 1] = c;
        val[i + 1] = v;
    }
}

} // namespace detail
} // namespace amgcl

#endif
