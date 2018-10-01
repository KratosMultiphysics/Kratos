#ifndef AMGCL_MPI_RELAXATION_AS_PRECONDITIONER_HPP
#define AMGCL_MPI_RELAXATION_AS_PRECONDITIONER_HPP

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
 * \file   amgcl/mpi/relaxation/as_preconditioner.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Use a distributed amgcl smoother as a standalone preconditioner.
 */

#include <vector>
#include <memory>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/mpi/util.hpp>

namespace amgcl {
namespace mpi {
namespace relaxation {

template <class Relaxation>
struct as_preconditioner {
    typedef typename Relaxation::params                params;
    typedef typename Relaxation::backend_type          backend_type;
    typedef typename backend_type::params              backend_params;
    typedef typename backend_type::value_type          value_type;
    typedef typename math::scalar_of<value_type>::type scalar_type;
    typedef distributed_matrix<backend_type>           matrix;
    typedef typename backend_type::vector              vector;

    template <class Matrix>
    as_preconditioner(
            communicator comm,
            const Matrix &A,
            const params &prm = params(),
            const backend_params &bprm = backend_params()
       ) : A(std::make_shared>(comm, A, backend::rows(A))),
           S(A, prm, bprm)
    {
        this->A->move_to_backend(bprm);
    }

    as_preconditioner(
            communicator,
            std::shared_ptr<matrix> A,
            const params &prm = params(),
            const backend_params &bprm = backend_params()
       ) : A(A), S(*A, prm, bprm)
    {
        this->A->move_to_backend(bprm);
    }

    template <class Vec1, class Vec2>
    void apply(const Vec1 &rhs, Vec2 &&x) const {
        S.apply(*A, rhs, x);
    }

    std::shared_ptr<matrix> system_matrix_ptr() const {
        return A;
    }

    const matrix& system_matrix() const {
        return *system_matrix_ptr();
    }

    private:
        std::shared_ptr<matrix> A;
        Relaxation S;

        friend std::ostream& operator<<(std::ostream &os, const as_preconditioner &p) {
            os << "Relaxation as preconditioner" << std::endl;
            os << "  unknowns: " << p.system_matrix().glob_rows() << std::endl;
            os << "  nonzeros: " << p.system_matrix().glob_nonzeros() << std::endl;

            return os;
        }
};

} // namespace relaxation
} // namespace mpi
} // namespace amgcl

#endif
