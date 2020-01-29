#ifndef AMGCL_MPI_SOLVER_CG_HPP
#define AMGCL_MPI_SOLVER_CG_HPP

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
 * \file   amgcl/mpi/solver/cg.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  MPI wrapper for CG iterative method.
 */

#include <amgcl/solver/cg.hpp>
#include <amgcl/mpi/inner_product.hpp>

namespace amgcl {
namespace mpi {
namespace solver {

template <class Backend, class InnerProduct = mpi::inner_product>
class cg : public amgcl::solver::cg<Backend, InnerProduct> {
    typedef amgcl::solver::cg<Backend, InnerProduct> Base;
    public:
        using Base::Base;
};

} // namespace solver
} // namespace mpi
} // namespace amgcl


#endif
