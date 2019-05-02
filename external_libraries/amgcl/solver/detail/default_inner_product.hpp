#ifndef AMGCL_SOLVER_DETAIL_DEFAULT_INNER_PRODUCT_HPP
#define AMGCL_SOLVER_DETAIL_DEFAULT_INNER_PRODUCT_HPP

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
 * \file   amgcl/solver/detail/default_inner_product.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Default inner product getter for iterative solvers.
 *
 * Falls through to backend::inner_product().
 */

#include <amgcl/backend/interface.hpp>

namespace amgcl {
namespace solver {
namespace detail {

struct default_inner_product {
    template <class Vec1, class Vec2>
    typename math::inner_product_impl<
        typename backend::value_type<Vec1>::type
    >::return_type
    operator()(const Vec1 &x, const Vec2 &y) const {
        return backend::inner_product(x, y);
    }
};

} // namespace detail
} // namespace solver
} // namespace amgcl


#endif
