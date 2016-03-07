#ifndef AMGCL_SOLVERS_BICGSTAB_HPP
#define AMGCL_SOLVERS_BICGSTAB_HPP

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
 * \file   amgcl/solver/bicgstab.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  BiCGStab iterative method.
 */

#include <boost/tuple/tuple.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/solver/detail/default_inner_product.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace solver {

/** BiConjugate Gradient Stabilized (BiCGSTAB) method.
 * \rst
 * The BiConjugate Gradient Stabilized method (Bi-CGSTAB) was developed to
 * solve nonsymmetric linear systems while avoiding the often irregular
 * convergence patterns of the Conjugate Gradient [Barr94]_.
 * \endrst
 */
template <
    class Backend,
    class InnerProduct = detail::default_inner_product
    >
class bicgstab {
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

            /// Target residual error.
            scalar_type tol;

            params(size_t maxiter = 100, scalar_type tol = 1e-8)
                : maxiter(maxiter), tol(tol)
            {}

            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, maxiter),
                  AMGCL_PARAMS_IMPORT_VALUE(p, tol)
            {}

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, maxiter);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, tol);
            }
        };

        /// Preallocates necessary data structures for the system of size \p n.
        bicgstab(
                size_t n,
                const params &prm = params(),
                const backend_params &backend_prm = backend_params(),
                const InnerProduct &inner_product = InnerProduct()
                )
            : prm(prm), n(n),
              r ( Backend::create_vector(n, backend_prm) ),
              p ( Backend::create_vector(n, backend_prm) ),
              v ( Backend::create_vector(n, backend_prm) ),
              s ( Backend::create_vector(n, backend_prm) ),
              t ( Backend::create_vector(n, backend_prm) ),
              rh( Backend::create_vector(n, backend_prm) ),
              ph( Backend::create_vector(n, backend_prm) ),
              sh( Backend::create_vector(n, backend_prm) ),
              inner_product(inner_product)
        { }

        /* Computes the solution for the given system matrix \p A and the
         * right-hand side \p rhs.  Returns the number of iterations made and
         * the achieved residual as a ``boost::tuple``. The solution vector
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
        boost::tuple<size_t, scalar_type> operator()(
                Matrix  const &A,
                Precond const &P,
                Vec1    const &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2          &x
#else
                Vec2          &&x
#endif
                ) const
        {
            static const coef_type one  = math::identity<coef_type>();
            static const coef_type zero = math::zero<coef_type>();

            backend::residual(rhs, A, x, *r);

            scalar_type norm_rhs = norm(rhs);
            if (norm_rhs < amgcl::detail::eps<scalar_type>(n)) {
                backend::clear(x);
                return boost::make_tuple(0, norm_rhs);
            }

            scalar_type eps = norm_rhs * prm.tol;

            coef_type rho1  = zero;
            coef_type rho2  = zero;
            coef_type alpha = zero;
            coef_type omega = zero;

            size_t      iter = 0;
            scalar_type res  = norm(*r);

            backend::copy(*r, *rh);

            for(bool first = true; res > eps && iter < prm.maxiter; ++iter) {

                rho2 = rho1;
                rho1 = inner_product(*r, *rh);

                if (first) {
                    backend::copy(*r, *p);
                    first = false;
                } else {
                    precondition(!math::is_zero(rho2), "Zero rho in BiCGStab");
                    coef_type beta = (rho1 * alpha) / (rho2 * omega);
                    backend::axpbypcz(one, *r, -beta * omega, *v, beta, *p);
                }

                P.apply(*p, *ph);

                backend::spmv(one, A, *ph, zero, *v);

                alpha = rho1 / inner_product(*rh, *v);

                backend::axpbypcz(one, *r, -alpha, *v, zero, *s);

                if ((res = norm(*s)) <= eps) {
                    backend::axpby(alpha, *ph, one, x);
                } else {
                    P.apply(*s, *sh);

                    backend::spmv(one, A, *sh, zero, *t);

                    omega = inner_product(*t, *s) / inner_product(*t, *t);

                    precondition(!math::is_zero(omega), "Zero omega in BiCGStab");

                    backend::axpbypcz(alpha, *ph, omega, *sh, one, x);
                    backend::axpbypcz(one, *s, -omega, *t, zero, *r);

                    res = norm(*r);
                }
            }

            return boost::make_tuple(iter, res / norm_rhs);
        }

        /* Computes the solution for the given right-hand side \p rhs. The
         * system matrix is the same that was used for the setup of the
         * preconditioner \p P.  Returns the number of iterations made and the
         * achieved residual as a ``boost::tuple``. The solution vector \p x
         * provides initial approximation in input and holds the computed
         * solution on output.
         */
        template <class Precond, class Vec1, class Vec2>
        boost::tuple<size_t, scalar_type> operator()(
                Precond const &P,
                Vec1    const &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2          &x
#else
                Vec2          &&x
#endif
                ) const
        {
            return (*this)(P.system_matrix(), P, rhs, x);
        }


    public:
        params prm;

    private:
        size_t n;

        boost::shared_ptr<vector> r;
        boost::shared_ptr<vector> p;
        boost::shared_ptr<vector> v;
        boost::shared_ptr<vector> s;
        boost::shared_ptr<vector> t;
        boost::shared_ptr<vector> rh;
        boost::shared_ptr<vector> ph;
        boost::shared_ptr<vector> sh;

        InnerProduct inner_product;

        template <class Vec>
        scalar_type norm(const Vec &x) const {
            return sqrt(math::norm(inner_product(x, x)));
        }
};

} // namespace solver
} // namespace amgcl


#endif
