#ifndef AMGCL_SOLVERS_CG_HPP
#define AMGCL_SOLVERS_CG_HPP

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
 * \file   amgcl/solver/cg.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Conjugate Gradient method.
 */

#include <tuple>
#include <amgcl/backend/interface.hpp>
#include <amgcl/solver/detail/default_inner_product.hpp>
#include <amgcl/util.hpp>

namespace amgcl {

/// Iterative solvers
namespace solver {

/**
 * \defgroup solvers
 * \brief Iterative solvers
 *
 * AMGCL provides several iterative solvers, but it should be easy to use it as
 * a preconditioner with a user-provided solver.  Each solver in AMGCL is a
 * class template. Its single template parameter specifies the backend to use.
 * This allows to preallocate necessary resources at class construction.
 * Obviously, the solver backend has to coincide with the AMG backend.
 */


/** Conjugate Gradients method.
 * \rst
 * An effective method for symmetric positive definite systems [Barr94]_.
 * \endrst
 */
template <
    class Backend,
    class InnerProduct = detail::default_inner_product
    >
class cg {
    public:
        typedef Backend backend_type;

        typedef typename Backend::vector     vector;
        typedef typename Backend::value_type value_type;
        typedef typename Backend::params     backend_params;

        typedef typename math::scalar_of<value_type>::type scalar_type;

        typedef typename math::inner_product_impl<
            typename math::rhs_of<value_type>::type
            >::return_type coef_type;

        /// Solver parameters.
        struct params {
            /// Maximum number of iterations.
            size_t maxiter;

            /// Target relative residual error.
            scalar_type tol;

            /// Target absolute residual error.
            scalar_type abstol;

            params()
                : maxiter(100), tol(1e-8),
                  abstol(std::numeric_limits<scalar_type>::min())
            {}

#ifndef AMGCL_NO_BOOST
            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, maxiter),
                  AMGCL_PARAMS_IMPORT_VALUE(p, tol),
                  AMGCL_PARAMS_IMPORT_VALUE(p, abstol)
            {
                check_params(p, {"maxiter", "tol", "abstol"});
            }

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, maxiter);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, tol);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, abstol);
            }
#endif
        };

        /// Preallocates necessary data structures for the system of size \p n.
        cg(
                size_t n,
                const params &prm = params(),
                const backend_params &backend_prm = backend_params(),
                const InnerProduct &inner_product = InnerProduct()
          ) : prm(prm), n(n),
              r(Backend::create_vector(n, backend_prm)),
              s(Backend::create_vector(n, backend_prm)),
              p(Backend::create_vector(n, backend_prm)),
              q(Backend::create_vector(n, backend_prm)),
              inner_product(inner_product)
        { }

        /* Computes the solution for the given system matrix \p A and the
         * right-hand side \p rhs.  Returns the number of iterations made and
         * the achieved residual as a ``std::tuple``. The solution vector
         * \p x provides initial approximation in input and holds the computed
         * solution on output.
         *
         * The system matrix may differ from the matrix used during
         * initialization. This may be used for the solution of non-stationary
         * problems with slowly changing coefficients. There is a strong chance
         * that a preconditioner built for a time step will act as a reasonably
         * good preconditioner for several subsequent time steps [DeSh12]_.
         */
        template <class Matrix, class Precond, class Vec1, class Vec2>
        std::tuple<size_t, scalar_type> operator()(
                const Matrix &A, const Precond &P, const Vec1 &rhs, Vec2 &&x) const
        {
            static const coef_type one  = math::identity<coef_type>();
            static const coef_type zero = math::zero<coef_type>();

            backend::residual(rhs, A, x, *r);
            scalar_type norm_rhs = norm(rhs);
            if (norm_rhs < amgcl::detail::eps<scalar_type>(1)) {
                backend::clear(x);
                return std::make_tuple(0, norm_rhs);
            }

            scalar_type eps  = std::max(prm.tol * norm_rhs, prm.abstol);

            coef_type rho1 = 2 * eps * one;
            coef_type rho2 = zero;
            scalar_type res_norm = norm(*r);

            size_t iter = 0;
            for(; iter < prm.maxiter && math::norm(res_norm) > eps; ++iter) {
                P.apply(*r, *s);

                rho2 = rho1;
                rho1 = inner_product(*r, *s);

                if (iter)
                    backend::axpby(one, *s, rho1 / rho2, *p);
                else
                    backend::copy(*s, *p);

                backend::spmv(one, A, *p, zero, *q);

                coef_type alpha = rho1 / inner_product(*q, *p);

                backend::axpby( alpha, *p, one,  x);
                backend::axpby(-alpha, *q, one, *r);

                res_norm = norm(*r);
            }

            return std::make_tuple(iter, res_norm / norm_rhs);
        }

        /* Computes the solution for the given right-hand side \p rhs. The
         * system matrix is the same that was used for the setup of the
         * preconditioner \p P.  Returns the number of iterations made and the
         * achieved residual as a ``std::tuple``. The solution vector \p x
         * provides initial approximation in input and holds the computed
         * solution on output.
         */
        template <class Precond, class Vec1, class Vec2>
        std::tuple<size_t, scalar_type> operator()(
                const Precond &P, const Vec1 &rhs, Vec2 &&x) const
        {
            return (*this)(P.system_matrix(), P, rhs, x);
        }

        size_t bytes() const {
            return
                backend::bytes(*r) +
                backend::bytes(*s) +
                backend::bytes(*p) +
                backend::bytes(*q);
        }

        friend std::ostream& operator<<(std::ostream &os, const cg &s) {
            return os
                << "Type:             CG"
                << "\nUnknowns:         " << s.n
                << "\nMemory footprint: " << human_readable_memory(s.bytes())
                << std::endl;
        }
    public:
        params prm;

    private:
        size_t n;

        std::shared_ptr<vector> r;
        std::shared_ptr<vector> s;
        std::shared_ptr<vector> p;
        std::shared_ptr<vector> q;

        InnerProduct inner_product;

        template <class Vec>
        scalar_type norm(const Vec &x) const {
            return sqrt(math::norm(inner_product(x, x)));
        }
};

} // namespace solver
} // namespace amgcl


#endif
