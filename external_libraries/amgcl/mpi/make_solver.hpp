#ifndef AMGCL_MPI_MAKE_SOLVER_HPP
#define AMGCL_MPI_MAKE_SOLVER_HPP

/*
The MIT License

Copyright (c) 2012-2016 Denis Demidov <dennis.demidov@gmail.com>

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
 * \brief  Iterative solver wrapper for distributed linear systmes.
 */

#include <boost/property_tree/ptree.hpp>
#include <boost/shared_ptr.hpp>

#include <mpi.h>

#include <amgcl/util.hpp>
#include <amgcl/mpi/inner_product.hpp>

namespace amgcl {
namespace mpi {

template <
    class Precond,
    template <class, class> class IterativeSolver
    >
class make_solver {
    public:
        typedef typename Precond::backend_type backend_type;
        typedef typename backend_type::matrix matrix;
        typedef typename backend_type::value_type value_type;
        typedef typename backend_type::params backend_params;
        typedef typename backend::builtin<value_type>::matrix build_matrix;
        typedef typename math::scalar_of<value_type>::type scalar_type;

        typedef IterativeSolver<backend_type, mpi::inner_product> Solver;


        struct params {
            typename Precond::params precond; ///< Preconditioner parameters.
            typename Solver::params  solver;  ///< Iterative solver parameters.

            params() {}

            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_CHILD(p, precond),
                  AMGCL_PARAMS_IMPORT_CHILD(p, solver)
            {}

            void get(boost::property_tree::ptree &p, const std::string &path = "") const
            {
                AMGCL_PARAMS_EXPORT_CHILD(p, path, precond);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, solver);
            }
        } prm;

        template <class Matrix>
        make_solver(
                MPI_Comm comm, const Matrix &A,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                ) :
            prm(prm), n(backend::rows(A)),
            P(comm, A, prm.precond, bprm),
            S(backend::rows(A), prm.solver, bprm, mpi::inner_product(comm))
        {}

        make_solver(
                MPI_Comm comm, boost::shared_ptr<build_matrix> A,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                ) :
            prm(prm), n(backend::rows(*A)),
            P(comm, A, prm.precond, bprm),
            S(comm, backend::rows(*A), prm.solver, bprm, mpi::inner_product(comm))
        {}

        template <class Matrix, class Vec1, class Vec2>
        boost::tuple<size_t, scalar_type> operator()(
                Matrix  const &A,
                Vec1    const &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2          &x
#else
                Vec2          &&x
#endif
                ) const
        {
            return S(A, P, rhs, x);
        }

        template <class Vec1, class Vec2>
        boost::tuple<size_t, scalar_type> operator()(
                Vec1    const &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2          &x
#else
                Vec2          &&x
#endif
                ) const
        {
            return S(P, rhs, x);
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
            (*this)(rhs, x);
        }

        const Precond& precond() const {
            return P;
        }

        const Solver& solver() const {
            return S;
        }

        const matrix& system_matrix() const {
            return P.system_matrix();
        }

        void get_params(boost::property_tree::ptree &p) const {
            prm.get(p);
        }

        size_t size() const {
            return n;
        }
    private:
        size_t n;

        Precond P;
        Solver  S;
};

} // namespace mpi
} // namespace amgcl

#endif
