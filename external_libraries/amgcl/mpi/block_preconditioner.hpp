#ifndef AMGCL_MPI_BLOCK_PRECONDITIONER_HPP
#define AMGCL_MPI_BLOCK_PRECONDITIONER_HPP

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
 * \file   amgcl/mpi/block_preconditioner.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Distributed block preconditioner.
 */

#include <vector>

#include <memory>

#include <mpi.h>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/inner_product.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>

namespace amgcl {
namespace mpi {

template <class Precond>
class block_preconditioner {
    public:
        typedef typename Precond::params       params;
        typedef typename Precond::backend_type backend_type;
        typedef typename backend_type::params  backend_params;

        typedef typename backend_type::value_type value_type;
        typedef typename backend_type::matrix     bmatrix;
        typedef distributed_matrix<backend_type>  matrix;

        template <class Matrix>
        block_preconditioner(
                communicator comm,
                const Matrix &Astrip,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
        {
            A = std::make_shared<matrix>(comm, Astrip, backend::rows(Astrip));
            P = std::make_shared<Precond>(A->local(), prm, bprm);
            A->set_local(P->system_matrix_ptr());
            A->move_to_backend(bprm);
        }

        block_preconditioner(
                communicator,
                std::shared_ptr<matrix> A,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                ) : A(A)
        {
            P = std::make_shared<Precond>(A->local(), prm, bprm);
            A->set_local(P->system_matrix_ptr());
            A->move_to_backend(bprm);
        }

        std::shared_ptr<matrix> system_matrix_ptr() const {
            return A;
        }

        const matrix& system_matrix() const {
            return *A;
        }

        template <class Vec1, class Vec2>
        void apply(const Vec1 &rhs, Vec2 &&x) const {
            P->apply(rhs, x);
        }
    private:
        std::shared_ptr<matrix>  A;
        std::shared_ptr<Precond> P;
};

} // namespace mpi
} // namespace amgcl

#endif
