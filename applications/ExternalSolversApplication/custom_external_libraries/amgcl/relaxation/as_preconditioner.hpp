#ifndef AMGCL_RELAXATION_AS_PRECONDITIONER_HPP
#define AMGCL_RELAXATION_AS_PRECONDITIONER_HPP

/*
The MIT License

Copyright (c) 2012-2015 Denis Demidov <dennis.demidov@gmail.com>

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
 * \brief  Use an amgcl smoother as standalone preconditioner.
 */

#include <vector>
#include <boost/shared_ptr.hpp>
#include <amgcl/backend/builtin.hpp>

namespace amgcl {
namespace relaxation {

/// Allows to use an amgcl smoother as standalone preconditioner.
template <class Backend, template <class> class Relax>
class as_preconditioner {
    public:
        typedef Relax<Backend>            smoother;

        typedef typename Backend::matrix  matrix;
        typedef typename Backend::vector  vector;
        typedef typename smoother::params sparams;
        typedef typename Backend::params  bparams;

        typedef typename Backend::value_type value_type;
        typedef typename backend::builtin<value_type>::matrix build_matrix;

        template <class Matrix>
        as_preconditioner(
                const Matrix &M,
                const sparams &sprm,
                const bparams &bprm
                )
            : prm(sprm)
        {
            boost::shared_ptr<build_matrix> m = boost::make_shared<build_matrix>(M);
            A = Backend::copy_matrix(m, bprm);
            S = boost::make_shared<smoother>(*m, sprm, bprm);
            tmp = Backend::create_vector(backend::rows(*m), bprm);
        }

        template <class Vec1, class Vec2>
        void apply(
                const Vec1 &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2       &x
#else
                Vec2       &&x
#endif
                ) const
        {
            backend::clear(x);
            S->apply_pre (*A, rhs, x, *tmp, prm);
            S->apply_post(*A, rhs, x, *tmp, prm);
        }

        const matrix& system_matrix() const {
            return *A;
        }
    private:
        sparams prm;

        boost::shared_ptr<matrix>   A;
        boost::shared_ptr<smoother> S;
        boost::shared_ptr<vector> tmp;
};

} // namespace relaxation
} // namespace amgcl

#endif
