#ifndef AMGCL_BACKEND_DETAIL_MIXING_HPP
#define AMGCL_BACKEND_DETAIL_MIXING_HPP

/*
The MIT License

Copyright (c) 2012-2020 Denis Demidov <dennis.demidov@gmail.com>
Copyright (c) 2016, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)

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
 * \file   amgcl/backend/detail/mixing.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Utilities for mixed-precision of mixed-block backends.
 */

#include <type_traits>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/interface.hpp>

namespace amgcl {
namespace backend {
namespace detail {

// Backend with scalar value_type of highest precision.

template <class B1, class B2, class Enable = void>
struct common_scalar_backend;

template <class B>
struct common_scalar_backend<B, B,
    typename std::enable_if<
        math::static_rows<typename B::value_type>::value == 1
        >::type >
{
    typedef B type;
};

template <class V1, class V2>
struct common_scalar_backend< backend::builtin<V1>, backend::builtin<V2>,
    typename std::enable_if<
        math::static_rows<V1>::value != 1 ||
        math::static_rows<V2>::value != 1
        >::type>
{
    typedef typename math::scalar_of<V1>::type S1;
    typedef typename math::scalar_of<V2>::type S2;

    typedef
        typename std::conditional<
            (sizeof(S1) > sizeof(S2)), backend::builtin<S1>, backend::builtin<S2>
            >::type
        type;
};

} // namespace detail
} // namespace backend
} // namespace amgcl



#endif
