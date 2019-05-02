#ifndef AMGCL_RELAXATION_SPAI0_HPP
#define AMGCL_RELAXATION_SPAI0_HPP

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
 * \file   amgcl/relaxation/spai0.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Sparse approximate inverse relaxation scheme.
 */

#include <memory>
#include <amgcl/backend/interface.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace relaxation {

/// Sparse approximate interface smoother.
/**
 * The inverse matrix is approximated with diagonal matrix.
 *
 * \tparam Backend Backend for temporary structures allocation.
 * \ingroup relaxation
 * \sa \cite Broker2002
 */
template <class Backend>
struct spai0 {
    typedef typename Backend::value_type      value_type;
    typedef typename Backend::matrix_diagonal matrix_diagonal;

    typedef typename math::scalar_of<value_type>::type scalar_type;
    /// Relaxation parameters.
    typedef amgcl::detail::empty_params params;

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    spai0( const Matrix &A, const params &, const typename Backend::params &backend_prm)
    {
        const size_t n = rows(A);

        auto m = std::make_shared< backend::numa_vector<value_type> >(n, false);

#pragma omp parallel for
        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            value_type  num = math::zero<value_type>();
            scalar_type den = math::zero<scalar_type>();

            for(auto a = backend::row_begin(A, i); a; ++a) {
                value_type v = a.value();
                scalar_type norm_v = math::norm(v);
                den += norm_v * norm_v;
                if (a.col() == i) num += v;
            }

            (*m)[i] = math::inverse(den) * num;
        }

        M = Backend::copy_vector(m, backend_prm);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        static const scalar_type one = math::identity<scalar_type>();
        backend::residual(rhs, A, x, tmp);
        backend::vmul(one, *M, tmp, one, x);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        static const scalar_type one = math::identity<scalar_type>();
        backend::residual(rhs, A, x, tmp);
        backend::vmul(one, *M, tmp, one, x);
    }

    template <class Matrix, class VectorRHS, class VectorX>
    void apply( const Matrix&, const VectorRHS &rhs, VectorX &x) const
    {
        backend::vmul(math::identity<scalar_type>(), *M, rhs, math::zero<scalar_type>(), x);
    }

    size_t bytes() const {
        return backend::bytes(*M);
    }

    std::shared_ptr<matrix_diagonal> M;
};

} // namespace relaxation
} // namespace amgcl

#endif
