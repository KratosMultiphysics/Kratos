#ifndef AMGCL_RELAXATION_PARALLEL_ILU0_HPP
#define AMGCL_RELAXATION_PARALLEL_ILU0_HPP

/*
The MIT License

Copyright (c) 2012-2017 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/relaxation/parallel_ilu0.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Parallel version of ILU0 relaxation.
 */

#include <cassert>
#include <vector>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>
#include <amgcl/relaxation/detail/ilu_solve.hpp>

namespace amgcl {
namespace relaxation {

/// Parallel ILU(0) smoother.
/**
 * \note This is a parallel version of ILU0 described in \cite chow2015fine.
 *
 * \param Backend Backend for temporary structures allocation.
 * \ingroup relaxation
 */
template <class Backend>
struct parallel_ilu0 {
    typedef typename Backend::value_type      value_type;
    typedef typename Backend::vector          vector;
    typedef typename Backend::matrix          matrix;
    typedef typename Backend::params          backend_params;
    typedef typename Backend::matrix_diagonal matrix_diagonal;

    typedef typename math::scalar_of<value_type>::type scalar_type;

    /// Relaxation parameters.
    struct params {
        /// Number of factorization sweeps.
        unsigned factor_sweeps;

        /// Number of Jacobi iterations during solution
        unsigned jacobi_iters;

        /// Damping factor.
        scalar_type damping;

        params(unsigned factor_sweeps = 2, unsigned jacobi_iters = 2, scalar_type damping = 1)
            : factor_sweeps(factor_sweeps), jacobi_iters(jacobi_iters), damping(damping) {}

        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, factor_sweeps)
            , AMGCL_PARAMS_IMPORT_VALUE(p, jacobi_iters)
            , AMGCL_PARAMS_IMPORT_VALUE(p, damping)
        {
            AMGCL_PARAMS_CHECK(p, (factor_sweeps)(jacobi_iters)(damping));
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, factor_sweeps);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, jacobi_iters);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, damping);
        }
    };

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    parallel_ilu0( const Matrix &A, const params &prm, const backend_params &bprm)
    {
        typedef typename backend::row_iterator<Matrix>::type  row_iterator;
        typedef typename backend::builtin<value_type>::matrix build_matrix;

        const size_t n = backend::rows(A);

        boost::shared_ptr<build_matrix> Lh = boost::make_shared<build_matrix>();
        boost::shared_ptr<build_matrix> Uh = boost::make_shared<build_matrix>();

        Lh->set_size(n, n); std::fill(Lh->ptr, Lh->ptr + n + 1, 0);
        Uh->set_size(n, n); std::fill(Uh->ptr, Uh->ptr + n + 1, 0);

        // Create an initial approximation for L and U by copying the
        // corresponding parts of A into the matrices.
        // Transpose U in the process (convert it to CSC).
        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            for(row_iterator a = backend::row_begin(A, i); a; ++a) {
                ptrdiff_t c = a.col();
                if (c < i) {
                    ++Lh->ptr[i+1];
                } else {
                    ++Uh->ptr[c+1];
                }
            }
        }

        std::partial_sum(Lh->ptr, Lh->ptr + n + 1, Lh->ptr);
        std::partial_sum(Uh->ptr, Uh->ptr + n + 1, Uh->ptr);

        Lh->set_nonzeros(Lh->ptr[n]);
        Uh->set_nonzeros(Uh->ptr[n]);

        for(ptrdiff_t i = 0, Lhead = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            for(row_iterator a = backend::row_begin(A, i); a; ++a) {
                ptrdiff_t c = a.col();
                if (c < i) {
                    Lh->col[Lhead] = c;
                    Lh->val[Lhead] = a.value();
                    ++Lhead;
                } else {
                    ptrdiff_t head = Uh->ptr[c]++;

                    Uh->col[head] = i;
                    Uh->val[head] = a.value();
                }
            }
        }

        std::rotate(Uh->ptr, Uh->ptr + n, Uh->ptr + n + 1);
        Uh->ptr[0] = 0;

        amgcl::backend::sort_rows(*Uh);

        // Do the required number of Chow-Patel sweeps to get an approximated
        // factorization.
        value_type * Ltmp = new value_type[Lh->ptr[n]];
        value_type * Utmp = new value_type[Uh->ptr[n]];

        for (unsigned sweep = 0; sweep < prm.factor_sweeps; ++sweep) {
#pragma omp parallel for schedule(dynamic,4096)
            for(ptrdiff_t row = 0; row < static_cast<ptrdiff_t>(n); ++row) {
                {
                    // Update U
                    ptrdiff_t U_head = Uh->ptr[row];
                    ptrdiff_t U_tail = Uh->ptr[row+1];

                    for(ptrdiff_t k = U_head; k < U_tail; ++k) {
                        ptrdiff_t i = Uh->col[k];
                        value_type s = math::zero<value_type>();

                        ptrdiff_t kl = Lh->ptr[i];
                        ptrdiff_t jl = (kl < Lh->ptr[i+1]) ? Lh->col[kl] : n;

                        for(ptrdiff_t ku = U_head; ku < k; ++ku) {
                            ptrdiff_t iu = Uh->col[ku];

                            while(jl < iu) {
                                ++kl;
                                jl = Lh->col[kl];
                            }

                            if (jl == iu)
                                s -= Lh->val[kl] * Uh->val[ku];
                        }

                        row_iterator a = backend::row_begin(A, i);
                        while(a.col() < row) ++a;

                        Utmp[k] = a.value() - s;
                    }
                }

                // Update L
                {
                    row_iterator a = backend::row_begin(A, row);
                    ptrdiff_t L_head = Lh->ptr[row];
                    ptrdiff_t L_tail = Lh->ptr[row + 1];

                    for(ptrdiff_t k = L_head; k < L_tail; ++k, ++a) {
                        ptrdiff_t j = Lh->col[k];
                        value_type s = a.value();

                        ptrdiff_t ku = Uh->ptr[j];
                        ptrdiff_t iu = Uh->col[ku];

                        for(ptrdiff_t kl = L_head; kl < k; ++kl) {
                            ptrdiff_t jl = Lh->col[kl];

                            while(iu < jl) {
                                ++ku;
                                iu = Uh->col[ku];
                            }

                            if (iu == jl)
                                s -= Lh->val[kl] * Uh->val[ku];
                        }

                        Ltmp[k] = math::inverse(Uh->val[Uh->ptr[j+1]-1]) * s;
                    }
                }
            }

            std::swap(Ltmp, Lh->val);
            std::swap(Utmp, Uh->val);
        }

        delete[] Ltmp;
        delete[] Utmp;

        // Convert U to CRS, split out its diagonal.
        std::vector<value_type> Dh(n);

        ptrdiff_t  * Uptr = new ptrdiff_t [n+1];
        ptrdiff_t  * Ucol = new ptrdiff_t [Uh->ptr[n] - n];
        value_type * Uval = new value_type[Uh->ptr[n] - n];

        std::fill(Uptr, Uptr + n + 1, 0);

        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            for(ptrdiff_t k = Uh->ptr[i]; k < Uh->ptr[i+1]; ++k) {
                ptrdiff_t j = Uh->col[k];
                if (j == i) {
                    Dh[i] = math::inverse(Uh->val[k]);
                } else {
                    ++Uptr[j+1];
                }
            }
        }

        std::partial_sum(Uptr, Uptr + n + 1, Uptr);

        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            for(ptrdiff_t k = Uh->ptr[i]; k < Uh->ptr[i+1]; ++k) {
                ptrdiff_t j = Uh->col[k];
                if (j != i) {
                    ptrdiff_t head = Uptr[j]++;
                    Ucol[head] = i;
                    Uval[head] = Uh->val[k];
                }
            }
        }

        std::rotate(Uptr, Uptr + n, Uptr + n + 1);
        Uptr[0] = 0;

        std::swap(Uptr, Uh->ptr); delete[] Uptr;
        std::swap(Ucol, Uh->col); delete[] Ucol;
        std::swap(Uval, Uh->val); delete[] Uval;

        amgcl::backend::sort_rows(*Uh);

        L = Backend::copy_matrix(Lh, bprm);
        U = Backend::copy_matrix(Uh, bprm);

        D  = Backend::copy_vector(Dh, bprm);
        t1 = Backend::create_vector(n, bprm);
        t2 = Backend::create_vector(n, bprm);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
            const params &prm
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        solve(tmp, prm);
        backend::axpby(prm.damping, tmp, math::identity<scalar_type>(), x);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
            const params &prm
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        solve(tmp, prm);
        backend::axpby(prm.damping, tmp, math::identity<scalar_type>(), x);
    }

    template <class Matrix, class VectorRHS, class VectorX>
    void apply(const Matrix&, const VectorRHS &rhs, VectorX &x, const params &prm) const
    {
        backend::copy(rhs, x);
        solve(x, prm);
    }

    private:
        boost::shared_ptr<matrix> L, U;
        boost::shared_ptr<matrix_diagonal> D;
        boost::shared_ptr<vector> t1, t2;

        template <class VectorX>
        void solve(VectorX &x, const params &prm) const
        {
            relaxation::detail::parallel_ilu_solve(
                    *L, *U, *D, x, *t1, *t2, prm.jacobi_iters);
        }
};

} // namespace relaxation
} // namespace amgcl



#endif
