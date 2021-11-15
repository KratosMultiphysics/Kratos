#ifndef AMGCL_SOLVER_PREONLY_HPP
#define AMGCL_SOLVER_PREONLY_HPP

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
 * \file   amgcl/solver/preonly.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Only apply preconditioner once
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
 */

/** Richardson iteration */
template <
    class Backend,
    class InnerProduct = detail::default_inner_product
    >
class preonly {
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
        typedef amgcl::detail::empty_params params;

        /// Preallocates necessary data structures for the system of size \p n.
        preonly(
                size_t n,
                const params& = params(),
                const backend_params& = backend_params(),
                const InnerProduct &inner_product = InnerProduct()
          ) : n(n), inner_product(inner_product)
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
                const Matrix&, const Precond &P, const Vec1 &rhs, Vec2 &&x) const
        {
            P.apply(rhs, x);
            return std::make_tuple(0, 0);
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
            return 0;
        }

        friend std::ostream& operator<<(std::ostream &os, const preonly &s) {
            return os
                << "Type:             PreOnly"
                << "\nUnknowns:         " << s.n
                << "\nMemory footprint: " << human_readable_memory(s.bytes())
                << std::endl;
        }

    private:
        size_t n;

        InnerProduct inner_product;

        template <class Vec>
        scalar_type norm(const Vec &x) const {
            return sqrt(math::norm(inner_product(x, x)));
        }
};

} // namespace solver
} // namespace amgcl


#endif
