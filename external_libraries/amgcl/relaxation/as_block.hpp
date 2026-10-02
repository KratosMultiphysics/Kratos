#ifndef AMGCL_RELAXATION_AS_BLOCK_HPP
#define AMGCL_RELAXATION_AS_BLOCK_HPP

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
 * \file   amgcl/relaxation/as_block.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Block matrix adapter for an amgcl smoother.
 */

#include <vector>
#include <memory>
#include <amgcl/util.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/adapter/block_matrix.hpp>

namespace amgcl {
namespace relaxation {

/// Converts input matrix to block format before constructing an amgcl smoother.
template <class BlockBackend, template <class> class Relax>
struct as_block {
    typedef typename BlockBackend::value_type BlockType;

    template <class Backend>
    class type {
        public:
            typedef Backend backend_type;

            typedef Relax<BlockBackend>       Base;

            typedef typename Backend::matrix  matrix;
            typedef typename Backend::vector  vector;
            typedef typename Base::params     params;
            typedef typename Backend::params  backend_params;

            typedef typename Backend::value_type value_type;
            typedef typename Backend::col_type   col_type;
            typedef typename Backend::ptr_type   ptr_type;
            typedef typename backend::builtin<value_type, col_type, ptr_type>::matrix build_matrix;

            template <class Matrix>
            type(
                    const Matrix &A,
                    const params &prm = params(),
                    const backend_params &bprm = backend_params()
                    )
            : base(*std::make_shared<typename backend::crs<BlockType, col_type, ptr_type>>(adapter::block_matrix<BlockType>(A)), prm, bprm),
              nrows(backend::rows(A) / math::static_rows<BlockType>::value)
            { }

            template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
            void apply_pre(
                    const Matrix &A,
                    const VectorRHS &rhs,
                    VectorX &x,
                    VectorTMP &tmp
                    ) const
            {
                auto F = backend::reinterpret_as_rhs<BlockType>(rhs);
                auto X = backend::reinterpret_as_rhs<BlockType>(x);
                auto T = backend::reinterpret_as_rhs<BlockType>(tmp);
                base.apply_pre(A, F, X, T);
            }

            template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
            void apply_post(
                    const Matrix &A,
                    const VectorRHS &rhs,
                    VectorX &x,
                    VectorTMP &tmp
                    ) const
            {
                auto F = backend::reinterpret_as_rhs<BlockType>(rhs);
                auto X = backend::reinterpret_as_rhs<BlockType>(x);
                auto T = backend::reinterpret_as_rhs<BlockType>(tmp);
                base.apply_post(A, F, X, T);
            }

            template <class Matrix, class Vec1, class Vec2>
            void apply(const Matrix &A, const Vec1 &rhs, Vec2 &&x) const {
                auto F = backend::reinterpret_as_rhs<BlockType>(rhs);
                auto X = backend::reinterpret_as_rhs<BlockType>(x);
                base.apply(A, F, X);
            }

            const matrix& system_matrix() const {
                return base.system_matrix();
            }

            std::shared_ptr<matrix> system_matrix_ptr() const {
                return base.system_matrix_ptr();
            }

            size_t bytes() const {
                return base.bytes();
            }
        private:
            Base base;
            size_t nrows;
    };
};

} // namespace relaxation
} // namespace amgcl

#endif
