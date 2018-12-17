#ifndef AMGCL_RELAXATION_SPAI1_HPP
#define AMGCL_RELAXATION_SPAI1_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/relaxation/spai1.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Sparse approximate inverse relaxation scheme.
 */

#include <vector>

#include <memory>
#include <amgcl/backend/interface.hpp>
#include <amgcl/util.hpp>
#include <amgcl/detail/qr.hpp>

namespace amgcl {
namespace relaxation {

/// Sparse approximate interface smoother.
/**
 * Sparsity pattern of the approximate inverse matrix coincides with that of A.
 *
 * \tparam Backend Backend for temporary structures allocation.
 * \ingroup relaxation
 * \sa \cite Broker2002
 */
template <class Backend>
struct spai1 {
    typedef typename Backend::value_type value_type;
    typedef typename Backend::vector     vector;

    typedef typename math::scalar_of<value_type>::type scalar_type;

    /// Relaxation parameters.
    typedef amgcl::detail::empty_params params;

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    spai1( const Matrix &A, const params &, const typename Backend::params &backend_prm)
    {
        typedef typename backend::value_type<Matrix>::type value_type;

        const size_t n = backend::rows(A);
        const size_t m = backend::cols(A);

        auto Ainv = std::make_shared<Matrix>(A);

#pragma omp parallel
        {
            std::vector<ptrdiff_t> marker(m, -1);
            std::vector<ptrdiff_t> I, J;
            std::vector<value_type> B, ek;
            amgcl::detail::QR<value_type> qr;

#pragma omp for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                ptrdiff_t row_beg = A.ptr[i];
                ptrdiff_t row_end = A.ptr[i + 1];

                I.assign(A.col + row_beg, A.col + row_end);

                J.clear();
                for(ptrdiff_t j = row_beg; j < row_end; ++j) {
                    ptrdiff_t c = A.col[j];

                    for(ptrdiff_t jj = A.ptr[c], ee = A.ptr[c + 1]; jj < ee; ++jj) {
                        ptrdiff_t cc = A.col[jj];
                        if (marker[cc] < 0) {
                            marker[cc] = 1;
                            J.push_back(cc);
                        }
                    }
                }
                std::sort(J.begin(), J.end());
                B.assign(I.size() * J.size(), math::zero<value_type>());
                ek.assign(J.size(), math::zero<value_type>());
                for(size_t j = 0; j < J.size(); ++j) {
                    marker[J[j]] = j;
                    if (J[j] == static_cast<ptrdiff_t>(i)) ek[j] = math::identity<value_type>();
                }

                for(ptrdiff_t j = row_beg; j < row_end; ++j) {
                    ptrdiff_t c = A.col[j];

                    for(auto a = row_begin(A, c); a; ++a)
                        B[marker[a.col()] + J.size() * (j - row_beg)] = a.value();
                }

                qr.solve(J.size(), I.size(), &B[0], &ek[0], &Ainv->val[row_beg],
                        amgcl::detail::col_major);

                for(size_t j = 0; j < J.size(); ++j)
                    marker[J[j]] = -1;
            }
        }

        M = Backend::copy_matrix(Ainv, backend_prm);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        backend::spmv(math::identity<scalar_type>(), *M, tmp, math::identity<scalar_type>(), x);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        backend::spmv(math::identity<scalar_type>(), *M, tmp, math::identity<scalar_type>(), x);
    }

    template <class Matrix, class VectorRHS, class VectorX>
    void apply(const Matrix&, const VectorRHS &rhs, VectorX &x) const
    {
        backend::spmv(math::identity<scalar_type>(), *M, rhs, math::zero<scalar_type>(), x);
    }

    size_t bytes() const {
        return backend::bytes(*M);
    }

    std::shared_ptr<typename Backend::matrix> M;
};

} // namespace relaxation
} // namespace amgcl


#endif
