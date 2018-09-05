#ifndef AMGCL_BACKEND_DETAIL_DEFAULT_DIRECT_SOLVER_HPP
#define AMGCL_BACKEND_DETAIL_DEFAULT_DIRECT_SOLVER_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/backend/detail/default_direct_solver.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Default direct solver for coarse level.
 */

#include <memory>
#include <amgcl/backend/builtin.hpp>

namespace amgcl {
namespace backend {
namespace detail {

template <class Backend>
struct default_direct_solver {
    typedef typename Backend::value_type   real;
    typedef typename Backend::matrix       matrix;
    typedef typename builtin<real>::matrix host_matrix;

    std::shared_ptr<matrix> Ainv;

    default_direct_solver(
            std::shared_ptr<host_matrix> A,
            typename Backend::params const &prm
            )
    {
        auto ainv = std::make_shared<host_matrix>();
        *ainv = inverse(*A);
        Ainv = Backend::copy_matrix(ainv, prm);
    }

    template <class Vec1, class Vec2>
    void operator()(const Vec1 &rhs, Vec2 &x) const {
        backend::spmv(1, *Ainv, rhs, 0, x);
    }

    static size_t coarse_enough() { return 500; }
};

} // namespace detail
} // namespace backend
} // namespace amgcl



#endif
