#ifndef AMGCL_RELAXATION_SPAI0_HPP
#define AMGCL_RELAXATION_SPAI0_HPP

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
 * \file   amgcl/relaxation/spai0.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Sparse approximate inverse relaxation scheme.
 */

#include <boost/shared_ptr.hpp>
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
    typedef typename Backend::value_type value_type;
    typedef typename Backend::vector     vector;

    /// Relaxation parameters.
    struct params {
        params() {}
        params(const boost::property_tree::ptree&) {}
        void get(boost::property_tree::ptree&, const std::string&) const {}
    };

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    spai0( const Matrix &A, const params &, const typename Backend::params &backend_prm)
    {
        typedef typename backend::row_iterator<Matrix>::type row_iterator;

        const size_t n = rows(A);

        std::vector<value_type> m(n);

#pragma omp parallel for
        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            value_type num = 0;
            value_type den = 0;

            for(row_iterator a = backend::row_begin(A, i); a; ++a) {
                value_type v = a.value();
                den += v * v;
                if (a.col() == i) num += v;
            }

            m[i] = num / den;
        }

        M = Backend::copy_vector(m, backend_prm);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
            const params&
            ) const
    {
        apply(A, rhs, x, tmp);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
            const params&
            ) const
    {
        apply(A, rhs, x, tmp);
    }

    private:
        boost::shared_ptr<vector> M;

        template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
        void apply(
                const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
                ) const
        {
            backend::residual(rhs, A, x, tmp);
            backend::vmul(1, *M, tmp, 1, x);
        }

};

} // namespace relaxation
} // namespace amgcl

#endif
