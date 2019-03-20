#ifndef AMGCL_SOLVER_DETAIL_GIVENS_ROTATIONS_HPP
#define AMGCL_SOLVER_DETAIL_GIVENS_ROTATIONS_HPP

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
 * \file   amgcl/solver/detail/givens_rotations.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Givens plane rotations used in GMRES variants.
 */

#include <amgcl/value_type/interface.hpp>

namespace amgcl {
namespace solver {
namespace detail {

template <class T>
inline void generate_plane_rotation(T dx, T dy, T &cs, T &sn) {
    if (math::is_zero(dy)) {
        cs = 1;
        sn = 0;
    } else if (std::abs(dy) > std::abs(dx)) {
        T tmp = dx / dy;
        sn = math::inverse(sqrt(math::identity<T>() + tmp * tmp));
        cs = tmp * sn;
    } else {
        T tmp = dy / dx;
        cs = math::inverse(sqrt(math::identity<T>() + tmp * tmp));
        sn = tmp * cs;
    }
}

template <class T>
void apply_plane_rotation(T &dx, T &dy, T cs, T sn) {
    T tmp = math::adjoint(cs) * dx + math::adjoint(sn) * dy;
    dy = -sn * dx + cs * dy;
    dx = tmp;
}

} // namespace detail
} // namespace solver
} // namespace amgcl

#endif
