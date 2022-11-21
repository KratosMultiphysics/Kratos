#ifndef AMGCL_DEFLATED_SOLVER_HPP
#define AMGCL_DEFLATED_SOLVER_HPP

/*
The MIT License

Copyright (c) 2012-2022 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/deflated_solver.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Iterative preconditioned solver with deflation.
 */

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>
#include <amgcl/detail/inverse.hpp>

namespace amgcl {

/// Convenience class that bundles together a preconditioner and an iterative solver.
template <
    class Precond,
    class IterativeSolver
    >
class deflated_solver : public amgcl::detail::non_copyable {
    static_assert(
            backend::backends_compatible<
                typename IterativeSolver::backend_type,
                typename Precond::backend_type
            >::value,
            "Backends for preconditioner and iterative solver should be compatible"
            );
    public:
        typedef typename IterativeSolver::backend_type backend_type;
        typedef typename backend_type::matrix matrix;
        typedef typename backend_type::vector vector;

        typedef typename backend_type::value_type value_type;
        typedef typename backend_type::params backend_params;
        typedef typename backend::builtin<value_type>::matrix build_matrix;

        typedef typename math::scalar_of<value_type>::type scalar_type;

        /** Combined parameters of the bundled preconditioner and the iterative
         * solver.
         */
        struct params {
            int         nvec; ///< The number of deflation vectors
            scalar_type *vec; ///< Deflation vectors as a [nvec x n] matrix

            typename Precond::params         precond; ///< Preconditioner parameters.
            typename IterativeSolver::params solver;  ///< Iterative solver parameters.

            params() : nvec(0), vec(nullptr) {}

#ifndef AMGCL_NO_BOOST
            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, nvec),
                  AMGCL_PARAMS_IMPORT_VALUE(p, vec),
                  AMGCL_PARAMS_IMPORT_CHILD(p, precond),
                  AMGCL_PARAMS_IMPORT_CHILD(p, solver)
            {
                check_params(p, {"nvec", "vec", "precond", "solver"});
            }

            void get( boost::property_tree::ptree &p,
                    const std::string &path = ""
                    ) const
            {
                AMGCL_PARAMS_EXPORT_CHILD(p, path, nvec);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, vec);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, precond);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, solver);
            }
#endif
        } prm;

        /** Sets up the preconditioner and creates the iterative solver. */
        template <class Matrix>
        deflated_solver(
                const Matrix &A,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                ) :
            prm(prm), n(backend::rows(A)),
            P(A, prm.precond, bprm),
            S(backend::rows(A), prm.solver, bprm),
            r(backend_type::create_vector(n, bprm)),
            Z(prm.nvec),
            E(prm.nvec * prm.nvec, 0),
            d(prm.nvec)
        {
            init(A, bprm);
        }

        // Constructs the preconditioner and creates iterative solver.
        // Takes shared pointer to the matrix in internal format.
        deflated_solver(
                std::shared_ptr<build_matrix> A,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                ) :
            prm(prm), n(backend::rows(*A)),
            P(A, prm.precond, bprm),
            S(backend::rows(*A), prm.solver, bprm),
            r(backend_type::create_vector(n, bprm)),
            Z(prm.nvec),
            E(prm.nvec * prm.nvec, 0),
            d(prm.nvec)
        {
            init(*A, bprm);
        }

        template <class Matrix>
        void init(const Matrix &A, const backend_params &bprm) {
            precondition(prm.nvec > 0 && prm.vec != nullptr, "Deflation vectors are not set!");

            for(int i = 0; i < prm.nvec; ++i) {
                Z[i] = backend_type::copy_vector(
                        std::make_shared<backend::numa_vector<scalar_type>>(make_iterator_range(prm.vec + n * i, prm.vec + n * (i + 1))),
                        bprm);
            }

            std::vector<scalar_type> AZ(prm.nvec);
            std::fill(E.begin(), E.end(), math::zero<scalar_type>());
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                std::fill(AZ.begin(), AZ.end(), math::zero<scalar_type>());
                for(auto a = backend::row_begin(A, i); a; ++a) {
                    for(int j = 0; j < prm.nvec; ++j) {
                        AZ[j] += a.value() * prm.vec[j * n + a.col()];
                    }
                }

                for(int ii = 0, k = 0; ii < prm.nvec; ++ii) {
                    for(int jj = 0; jj < prm.nvec; ++jj, ++k) {
                        E[k] += prm.vec[i + ii * n] * AZ[jj];
                    }
                }
            }

            std::vector<scalar_type> t(E.size());
            std::vector<int> p(prm.nvec);
            detail::inverse(prm.nvec, E.data(), t.data(), p.data());
        }

        /** Computes the solution for the given system matrix \p A and the
         * right-hand side \p rhs.  Returns the number of iterations made and
         * the achieved residual as a ``std::tuple``. The solution vector
         * \p x provides initial approximation in input and holds the computed
         * solution on output.
         *
         * \rst
         * The system matrix may differ from the matrix used during
         * initialization. This may be used for the solution of non-stationary
         * problems with slowly changing coefficients. There is a strong chance
         * that a preconditioner built for a time step will act as a reasonably
         * good preconditioner for several subsequent time steps [DeSh12]_.
         * \endrst
         */
        template <class Matrix, class Vec1, class Vec2>
        std::tuple<size_t, scalar_type> operator()(
                const Matrix &A, const Vec1 &rhs, Vec2 &&x) const
        {
            project(rhs, x);
            return S(A, *this, rhs, x);
        }

        /** Computes the solution for the given right-hand side \p rhs.
         * Returns the number of iterations made and the achieved residual as a
         * ``std::tuple``. The solution vector \p x provides initial
         * approximation in input and holds the computed solution on output.
         */
        template <class Vec1, class Vec2>
        std::tuple<size_t, scalar_type> operator()(const Vec1 &rhs, Vec2 &&x) const {
            project(rhs, x);
            return S(*this, rhs, x);
        }

        template <class Vec1, class Vec2>
        void apply(const Vec1 &rhs, Vec2 &&x) const {
            P.apply(rhs, x);
            project(rhs, x);
        }

        template <class Vec1, class Vec2>
        void project(const Vec1 &b, Vec2 &x) const {
            // x += Z^T E^{-1} Z (b - Ax)
            backend::residual(b, P.system_matrix(), x, *r);
            std::fill(d.begin(), d.end(), math::zero<scalar_type>());
            for(int j = 0; j < prm.nvec; ++j) {
                auto fj = backend::inner_product(*Z[j], *r);
                for(int i = 0; i < prm.nvec; ++i)
                    d[i] += E[i*prm.nvec+j] * fj;
            }
            backend::lin_comb(prm.nvec, d, Z, 1, x);
        }

        /// Returns reference to the constructed preconditioner.
        const Precond& precond() const {
            return P;
        }

        /// Returns reference to the constructed preconditioner.
        Precond& precond() {
            return P;
        }

        /// Returns reference to the constructed iterative solver.
        const IterativeSolver& solver() const {
            return S;
        }

        /// Returns the system matrix in the backend format.
        std::shared_ptr<typename Precond::matrix> system_matrix_ptr() const {
            return P.system_matrix_ptr();
        }

        typename Precond::matrix const& system_matrix() const {
            return P.system_matrix();
        }

#ifndef AMGCL_NO_BOOST
        /// Stores the parameters used during construction into the property tree \p p.
        void get_params(boost::property_tree::ptree &p) const {
            prm.get(p);
        }
#endif

        /// Returns the size of the system matrix.
        size_t size() const {
            return n;
        }

        size_t bytes() const {
            return backend::bytes(S) + backend::bytes(P);
        }

        friend std::ostream& operator<<(std::ostream &os, const deflated_solver &p) {
            return os
                << "Solver\n======\n" << p.S << std::endl
                << "Preconditioner\n==============\n" << p.P;
        }
    private:
        size_t           n;
        Precond          P;
        IterativeSolver  S;
        std::shared_ptr<vector> r;
        std::vector<std::shared_ptr<vector>> Z;
        std::vector<scalar_type> E;
        mutable std::vector<scalar_type> d;
};

} // namespace amgcl


#endif
