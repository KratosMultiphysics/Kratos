#ifndef AMGCL_COARSENING_AS_SCALAR_HPP
#define AMGCL_COARSENING_AS_SCALAR_HPP

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
 * \file   amgcl/coarsening/as_scalar.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Scalar coarsening for block matrices.
 */

#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/interface.hpp>

namespace amgcl {
namespace coarsening {

// Takes a block matrix as input,
// converts it to the scalar format,
// applies the base coarsening,
// converts the results back to block format.
template <template <class> class Coarsening>
struct as_scalar {
    template <class Backend>
    struct type {
        typedef typename math::element_of<typename Backend::value_type>::type Scalar;
        typedef backend::builtin<Scalar, typename Backend::col_type, typename Backend::ptr_type> BaseBackend;
        typedef Coarsening<BaseBackend> Base;

        typedef typename Base::params params;
        Base base;

        type(const params &prm = params()) : base(prm) {};

        template <class Matrix>
        typename std::enable_if<
            backend::coarsening_is_supported<BaseBackend, Coarsening>::value &&
            (math::static_rows<typename backend::value_type<Matrix>::type>::value > 1),
            std::tuple<
                std::shared_ptr<Matrix>,
                std::shared_ptr<Matrix>
                >
            >::type
        transfer_operators(const Matrix &B) {
            typedef typename backend::value_type<Matrix>::type Block;
            auto T = base.transfer_operators(*adapter::unblock_matrix(B));

            auto &P = *std::get<0>(T);
            auto &R = *std::get<1>(T);

            backend::sort_rows(P);
            backend::sort_rows(R);

            return std::make_tuple(
                    std::make_shared<Matrix>(adapter::block_matrix<Block>(P)),
                    std::make_shared<Matrix>(adapter::block_matrix<Block>(R))
                    );
        }

        template <class Matrix>
        typename std::enable_if<
            backend::coarsening_is_supported<BaseBackend, Coarsening>::value &&
            (math::static_rows<typename backend::value_type<Matrix>::type>::value == 1),
            std::tuple<
                std::shared_ptr<Matrix>,
                std::shared_ptr<Matrix>
                >
            >::type
        transfer_operators(const Matrix &A) {
            return base.transfer_operators(A);
        }

        template <class Matrix>
        typename std::enable_if<
            !backend::coarsening_is_supported<BaseBackend, Coarsening>::value,
            std::tuple<
                std::shared_ptr<Matrix>,
                std::shared_ptr<Matrix>
                >
            >::type
        transfer_operators(const Matrix&) {
            throw std::logic_error("The coarsening is not supported by the backend");
        }

        template <class Matrix>
        std::shared_ptr<Matrix>
        coarse_operator(const Matrix &A, const Matrix &P, const Matrix &R) const {
            return base.coarse_operator(A, P, R);
        }
    };
};

} // namespace coarsening
} // namespace amgcl

#endif
