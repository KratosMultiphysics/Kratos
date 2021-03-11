#ifndef AMGCL_MPI_RELAXATION_GAUSS_SEIDEL_HPP
#define AMGCL_MPI_RELAXATION_GAUSS_SEIDEL_HPP

/*
The MIT License

Copyright (c) 2012-2020 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/mpi/relaxation/gauss_seidel.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Distributed memory Gauss-Seidel relaxation scheme.
 */

#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>

namespace amgcl {
namespace mpi {
namespace relaxation {

template <class Backend>
struct gauss_seidel : public amgcl::relaxation::gauss_seidel<Backend> {
    typedef Backend backend_type;
    typedef amgcl::relaxation::gauss_seidel<Backend> Base;
    typedef typename Backend::params backend_params;
    typedef typename Base::params params;

    gauss_seidel(
            const distributed_matrix<Backend> &A,
            const params &prm = params(),
            const backend_params &bprm = backend_params()
         ) : Base(*A.local(), prm, bprm)
    {}

    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &t) const
    {
        Base::apply_pre(*A.local_backend(), rhs, x, t);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &t) const
    {
        Base::apply_post(*A.local_backend(), rhs, x, t);
    }

    template <class Matrix, class VectorRHS, class VectorX>
    void apply(const Matrix &A, const VectorRHS &rhs, VectorX &x) const
    {
        Base::apply(*A.local_backend(), rhs, x);
    }
};

} // namespace
} // mpi
} // amgcl

#endif
