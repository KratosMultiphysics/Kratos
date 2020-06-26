#ifndef AMGCL_ADAPTER_ZERO_COPY_HPP
#define AMGCL_ADAPTER_ZERO_COPY_HPP

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
\file    amgcl/adapter/zero_copy.hpp
\author  Denis Demidov <dennis.demidov@gmail.com>
\brief   Zero-copy adapter for input matrix in CRS format.
\ingroup adapters
*/

#include <type_traits>

#include <amgcl/util.hpp>
#include <amgcl/backend/builtin.hpp>

namespace amgcl {
namespace adapter {

template <typename Ptr, typename Col, typename Val>
std::shared_ptr< backend::crs<Val> >
zero_copy(size_t n, const Ptr *ptr, const Col *col, const Val *val) {
    // Check that Ptr and Col types are binary-compatible with ptrdiff_t:
    static_assert(std::is_integral<Ptr>::value, "Unsupported Ptr type");
    static_assert(std::is_integral<Col>::value, "Unsupported Col type");
    static_assert(sizeof(Ptr) == sizeof(ptrdiff_t), "Unsupported Ptr type");
    static_assert(sizeof(Col) == sizeof(ptrdiff_t), "Unsupported Col type");

    auto A = std::make_shared< backend::crs<Val> >();
    A->nrows = n;
    A->ncols = n;
    A->nnz   = ptr[n];

    A->ptr = (ptrdiff_t*)ptr;
    A->col = (ptrdiff_t*)col;
    A->val = (Val*)val;

    A->own_data = false;

    return A;
}

} // namespace adapter
} // namespace amgcl

#endif
