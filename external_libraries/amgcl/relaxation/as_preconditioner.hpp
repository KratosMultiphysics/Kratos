#ifndef AMGCL_RELAXATION_AS_PRECONDITIONER_HPP
#define AMGCL_RELAXATION_AS_PRECONDITIONER_HPP

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
 * \file   amgcl/relaxation/as_preconditioner.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Use an amgcl smoother as a standalone preconditioner.
 */

#include <vector>
#include <memory>
#include <amgcl/backend/builtin.hpp>

namespace amgcl {
namespace relaxation {

/// Allows to use an amgcl smoother as standalone preconditioner.
template <class Backend, template <class> class Relax>
class as_preconditioner {
    public:
        typedef Backend backend_type;

        typedef Relax<Backend>            smoother;

        typedef typename Backend::matrix  matrix;
        typedef typename Backend::vector  vector;
        typedef typename smoother::params params;
        typedef typename Backend::params  backend_params;

        typedef typename Backend::value_type value_type;
        typedef typename backend::builtin<value_type>::matrix build_matrix;

        template <class Matrix>
        as_preconditioner(
                const Matrix &M,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
            : prm(prm)
        {
            init(std::make_shared<build_matrix>(M), bprm);
        }

        as_preconditioner(
                std::shared_ptr<build_matrix> M,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
            : prm(prm)
        {
            init(M, bprm);
        }

        template <class Vec1, class Vec2>
        void apply(const Vec1 &rhs, Vec2 &&x) const {
            S->apply(*A, rhs, x);
        }

        const matrix& system_matrix() const {
            return *A;
        }

        std::shared_ptr<matrix> system_matrix_ptr() const {
            return A;
        }

        size_t bytes() const {
            size_t b = 0;

            if (A) b += backend::bytes(*A);
            if (S) b += backend::bytes(*S);

            return b;
        }
    private:
        params prm;

        std::shared_ptr<matrix>   A;
        std::shared_ptr<smoother> S;

        void init(std::shared_ptr<build_matrix> M, const backend_params &bprm) {
            A = Backend::copy_matrix(M, bprm);
            S = std::make_shared<smoother>(*M, prm, bprm);
        }

        friend std::ostream& operator<<(std::ostream &os, const as_preconditioner &p) {
            os << "Relaxation as preconditioner" << std::endl;
            os << "  unknowns: " << backend::rows(p.system_matrix()) << std::endl;
            os << "  nonzeros: " << backend::nonzeros(p.system_matrix()) << std::endl;
            os << "  memory:   " << human_readable_memory(p.bytes()) << std::endl;

            return os;
        }
};

} // namespace relaxation
} // namespace amgcl

#endif
