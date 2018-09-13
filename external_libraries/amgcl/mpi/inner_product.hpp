#ifndef AMGCL_MPI_INNER_PRODUCT_HPP
#define AMGCL_MPI_INNER_PRODUCT_HPP

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
 * \file   amgcl/mpi/inner_product.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Inner product for distributed vectors.
 */

#include <mpi.h>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/interface.hpp>
#include <amgcl/mpi/util.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace mpi {

struct inner_product {
    communicator comm;

    inner_product(communicator comm) : comm(comm) {}

    template <class Vec1, class Vec2>
    typename math::inner_product_impl<
        typename backend::value_type<Vec1>::type
        >::return_type
    operator()(const Vec1 &x, const Vec2 &y) const {
        typedef typename backend::value_type<Vec1>::type value_type;
        typedef typename math::inner_product_impl<value_type>::return_type coef_type;

        AMGCL_TIC("inner product");
        coef_type sum = comm.reduce(MPI_SUM, backend::inner_product(x, y));
        AMGCL_TOC("inner product");

        return sum;
    }
};

} // namespace mpi
} // namespace amgcl

#endif
