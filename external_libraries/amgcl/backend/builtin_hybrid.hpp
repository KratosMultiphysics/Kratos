#ifndef AMGCL_BACKEND_BUILTIN_HYBRID_HPP
#define AMGCL_BACKEND_BUILTIN_HYBRID_HPP

/*
The MIT License

Copyright (c) 2012-2022 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/backend/builtin_hybrid.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Builtin backend that uses scalar matrices to build the hierarchy, but stores the computed matrix in block format.
 */

#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/block_matrix.hpp>

namespace amgcl {
namespace backend {

// Hybrid backend uses scalar matrices to build the hierarchy,
// but stores the computed matrices in the block format.
template <typename ScalarType, typename BlockType, typename ColumnType = ptrdiff_t, typename PointerType = ColumnType>
struct builtin_hybrid : public builtin<ScalarType> {
    typedef builtin<ScalarType, ColumnType, PointerType> Base;
    typedef crs<BlockType, ColumnType, PointerType> matrix;
    struct provides_row_iterator : std::false_type {};

    static std::shared_ptr<matrix>
    copy_matrix(std::shared_ptr<typename Base::matrix> As, const typename Base::params&)
    {
        return std::make_shared<matrix>(amgcl::adapter::block_matrix<BlockType>(*As));
    }
};

template <typename T1, typename B1, typename T2, typename B2>
struct backends_compatible< builtin_hybrid<T1, B1>, builtin_hybrid<T2, B2> > : std::true_type {};

template <typename T1, typename T2, typename B2>
struct backends_compatible< builtin<T1>, builtin_hybrid<T2, B2> > : std::true_type {};

template <typename T1, typename B1, typename T2>
struct backends_compatible< builtin_hybrid<T1, B1>, builtin<T2> > : std::true_type {};

} // namespace backend
} // namespace amgcl
#endif
