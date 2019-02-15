#ifndef AMGCL_MAKE_SOLVER_HPP
#define AMGCL_MAKE_SOLVER_HPP

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
 * \file   amgcl/make_solver.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Tie an iterative solver and a preconditioner in a single class.
 */

#include <type_traits>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>

namespace amgcl {

/// Convenience class that bundles together a preconditioner and an iterative solver.
template <
    class Precond,
    class IterativeSolver
    >
class make_solver : public amgcl::detail::non_copyable {
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

        typedef typename backend_type::value_type value_type;
        typedef typename backend_type::params backend_params;
        typedef typename backend::builtin<value_type>::matrix build_matrix;

        typedef typename math::scalar_of<value_type>::type scalar_type;

        /** Combined parameters of the bundled preconditioner and the iterative
         * solver.
         */
        struct params {
            typename Precond::params         precond; ///< Preconditioner parameters.
            typename IterativeSolver::params solver;  ///< Iterative solver parameters.

            params() {}

#ifndef AMGCL_NO_BOOST
            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_CHILD(p, precond),
                  AMGCL_PARAMS_IMPORT_CHILD(p, solver)
            {
                check_params(p, {"precond", "solver"});
            }

            void get( boost::property_tree::ptree &p,
                    const std::string &path = ""
                    ) const
            {
                AMGCL_PARAMS_EXPORT_CHILD(p, path, precond);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, solver);
            }
#endif
        } prm;

        /** Sets up the preconditioner and creates the iterative solver. */
        template <class Matrix>
        make_solver(
                const Matrix &A,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                ) :
            prm(prm), n(backend::rows(A)),
            P(A, prm.precond, bprm),
            S(backend::rows(A), prm.solver, bprm)
        {}

        // Constructs the preconditioner and creates iterative solver.
        // Takes shared pointer to the matrix in internal format.
        make_solver(
                std::shared_ptr<build_matrix> A,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                ) :
            prm(prm), n(backend::rows(*A)),
            P(A, prm.precond, bprm),
            S(backend::rows(*A), prm.solver, bprm)
        {}

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
            return S(A, P, rhs, x);
        }

        /** Computes the solution for the given right-hand side \p rhs.
         * Returns the number of iterations made and the achieved residual as a
         * ``std::tuple``. The solution vector \p x provides initial
         * approximation in input and holds the computed solution on output.
         */
        template <class Vec1, class Vec2>
        std::tuple<size_t, scalar_type> operator()(const Vec1 &rhs, Vec2 &&x) const {
            return S(P, rhs, x);
        }

        /** Acts as a preconditioner. That is, applies the solver to the
         * right-hand side \p rhs to get the solution \p x with zero initial
         * approximation.  Iterative methods usually use estimated residual for
         * exit condition.  For some problems the value of the estimated
         * residual can get too far from the true residual due to round-off
         * errors.  Nesting iterative solvers in this way may allow to shave
         * the last bits off the error. The method should not be used directly
         * but rather allows nesting ``make_solver`` classes as in the
         * following example:
         *
         * \rst
         * .. code-block:: cpp
         *
         *   typedef amgcl::make_solver<
         *     amgcl::make_solver<
         *       amgcl::amg<
         *         Backend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::spai0
         *         >,
         *       amgcl::solver::cg<Backend>
         *       >,
         *     amgcl::solver::cg<Backend>
         *     > NestedSolver;
         * \endrst
         */
        template <class Vec1, class Vec2>
        void apply(const Vec1 &rhs, Vec2 &&x) const {
            backend::clear(x);
            (*this)(rhs, x);
        }

        /// Returns reference to the constructed preconditioner.
        const Precond& precond() const {
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

        friend std::ostream& operator<<(std::ostream &os, const make_solver &p) {
            return os
                << "Solver\n======\n" << p.S << std::endl
                << "Preconditioner\n==============\n" << p.P;
        }
    private:
        size_t           n;
        Precond          P;
        IterativeSolver  S;
};

namespace detail {

template <class Matrix>
struct scaled_matrix {
    typedef typename backend::value_type<Matrix>::type value_type;

    const Matrix &A;
    const value_type *W;

    scaled_matrix(const Matrix &A, const std::vector<value_type> &W)
        : A(A), W(&W[0])
    {}

    size_t rows()     const { return backend::rows(A);     }
    size_t cols()     const { return backend::cols(A);     }
    size_t nonzeros() const { return backend::nonzeros(A); }

    struct row_iterator : public backend::row_iterator<Matrix>::type {
        typedef typename backend::row_iterator<Matrix>::type Base;
        typedef ptrdiff_t  col_type;
        typedef value_type val_type;

        value_type Wi;
        const value_type *Wj;

        row_iterator(const Matrix &A, const value_type *W, size_t i)
            : Base(A, i), Wi(W[i]), Wj(W) {}

        val_type value() const {
            return Wi * Wj[this->col()] * static_cast<const Base*>(this)->value();
        }
    };

    row_iterator row_begin(size_t i) const {
        return row_iterator(A, W, i);
    }
};

template <class Matrix>
scaled_matrix<Matrix> make_scaled_matrix(
        const Matrix &A,
        const std::vector<typename backend::value_type<Matrix>::type> &W)
{
    return scaled_matrix<Matrix>(A, W);
}

} // namespace detail

namespace backend {

template <class Matrix>
struct value_type< amgcl::detail::scaled_matrix<Matrix> >
{
    typedef typename backend::value_type<Matrix>::type type;
};

template <class Matrix>
struct rows_impl< amgcl::detail::scaled_matrix<Matrix> >
{
    static size_t get(const amgcl::detail::scaled_matrix<Matrix> &A) {
        return A.rows();
    }
};

template <class Matrix>
struct cols_impl< amgcl::detail::scaled_matrix<Matrix> >
{
    static size_t get(const amgcl::detail::scaled_matrix<Matrix> &A) {
        return A.cols();
    }
};

template <class Matrix>
struct nonzeros_impl< amgcl::detail::scaled_matrix<Matrix> >
{
    static size_t get(const amgcl::detail::scaled_matrix<Matrix> &A) {
        return A.nonzeros();
    }
};

template <class Matrix>
struct row_iterator< amgcl::detail::scaled_matrix<Matrix> >
{
    typedef typename amgcl::detail::scaled_matrix<Matrix>::row_iterator type;
};

template <class Matrix>
struct row_begin_impl< amgcl::detail::scaled_matrix<Matrix> >
{
    typedef amgcl::detail::scaled_matrix<Matrix> M;
    static typename row_iterator<M>::type get(const M &matrix, size_t row) {
        return matrix.row_begin(row);
    }
};

} // namespace backend

/// Wrapper for make_solver that scales the matrix and the RHS so that matrix has unit diagonal
template <
    class Precond,
    class IterativeSolver
    >
class make_scaling_solver : public amgcl::detail::non_copyable {
    public:
        typedef typename Precond::backend_type             backend_type;
        typedef typename backend_type::value_type          value_type;
        typedef typename backend_type::params              backend_params;
        typedef typename backend_type::vector              vector;
        typedef typename math::scalar_of<value_type>::type scalar_type;

        typedef typename make_solver<Precond, IterativeSolver>::params params;

        template <class Matrix>
        make_scaling_solver(
                const Matrix &A,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
        {
            const ptrdiff_t n = backend::rows(A);
            std::vector<value_type> w(n);

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n; ++i) {
                scalar_type sum = math::zero<scalar_type>();

                for(auto a = backend::row_begin(A, i); a; ++a)
                    sum += a.value() * a.value();

                w[i] = 1 / sqrt(sqrt(sum));
            }

            S = std::make_shared<Solver>(amgcl::detail::scaled_matrix<Matrix>(A, w), prm, bprm);
            W = backend_type::copy_vector(w, bprm);
            t = backend_type::create_vector(n, bprm);
        }

        template <class Vec1, class Vec2>
        std::tuple<size_t, scalar_type> operator()(Vec1 const &rhs, Vec2 &&x) const {
            backend::vmul(math::identity<scalar_type>(), *W, rhs, math::zero<scalar_type>(), *t);
            std::tuple<size_t, scalar_type> c = (*S)(*t, x);
            backend::vmul(math::identity<scalar_type>(), *W, x, math::zero<scalar_type>(), x);
            return c;
        }

        template <class Vec1, class Vec2>
        void apply(const Vec1 &rhs, Vec2 &&x) const {
            backend::clear(x);
            backend::vmul(math::identity<scalar_type>(), *W, rhs, math::zero<scalar_type>(), *t);
            (*this)(*t, x);
            backend::vmul(math::identity<scalar_type>(), *W, x, math::zero<scalar_type>(), x);
        }

        /// Returns reference to the constructed preconditioner.
        const Precond& precond() const {
            return S->precond();
        }
    private:
        typedef make_solver<Precond, IterativeSolver> Solver;

        std::shared_ptr<Solver> S;
        std::shared_ptr<vector> W;
        std::shared_ptr<vector> t;
};

} // namespace amgcl

#endif
