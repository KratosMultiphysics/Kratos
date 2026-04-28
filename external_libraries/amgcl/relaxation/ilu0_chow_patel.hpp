#ifndef AMGCL_RELAXATION_ILU0_CHOW_PATEL_HPP
#define AMGCL_RELAXATION_ILU0_CHOW_PATEL_HPP

/*
The MIT License

Copyright (c) 2026 Vicente Mataix Ferrandiz <vicente.mataix-ferrandiz@siemens.com>

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
 * \file   amgcl/relaxation/ilu0_chow_patel.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Fine-grained parallel ILU(0) factorization.
 *
 * Parallel ILU(0) factorization based on the iterative algorithm described in:
 *
 *   Edmond Chow and Aftab Patel,
 *   "Fine-Grained Parallel Incomplete LU Factorization",
 *   SIAM J. Sci. Comput., 37(2), C169-C193, 2015.
 *
 * Instead of a sequential Gaussian-elimination-based factorization, the L and
 * U factors are computed by solving the nonlinear constraint equations
 * (LU)_{ij} = a_{ij} via asynchronous fixed-point iteration sweeps.  Each
 * nonzero of L and U can be updated independently, giving fine-grained
 * parallelism that scales well regardless of matrix ordering.
 *
 * The matrix is row-scaled to have a unit diagonal before factorization,
 * as recommended by the paper for convergence.  The scaling is undone
 * when constructing the final factors for the preconditioner.
 */

#include <vector>
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>
#include <amgcl/relaxation/detail/ilu_solve.hpp>

namespace amgcl {
namespace relaxation {

/// Fine-grained parallel ILU(0) smoother (Chow-Patel algorithm).
/**
 * Uses an iterative fixed-point method to compute the ILU(0) factorization
 * in parallel.  All nonzeros of L and U can be updated simultaneously,
 * providing parallelism that is independent of the matrix ordering.
 *
 * \param Backend Backend for temporary structures allocation.
 * \ingroup relaxation
 */
template <class Backend>
struct ilu0_chow_patel {
    typedef typename Backend::value_type      value_type;
    typedef typename Backend::col_type        col_type;
    typedef typename Backend::ptr_type        ptr_type;
    typedef typename Backend::vector          vector;
    typedef typename Backend::matrix          matrix;
    typedef typename Backend::matrix_diagonal matrix_diagonal;

    typedef typename math::scalar_of<value_type>::type scalar_type;
    typedef detail::ilu_solve<Backend> ilu_solve;

    /// Relaxation parameters.
    struct params {
        /// Damping factor.
        scalar_type damping;

        /// Number of fixed-point iteration sweeps for computing L and U.
        int sweeps;

        /// Parameters for sparse triangular system solver.
        typename ilu_solve::params solve;

        params() : damping(1), sweeps(2) {}

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, damping)
            , AMGCL_PARAMS_IMPORT_VALUE(p, sweeps)
            , AMGCL_PARAMS_IMPORT_CHILD(p, solve)
        {
            check_params(p, {"damping", "sweeps", "solve"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, damping);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, sweeps);
            AMGCL_PARAMS_EXPORT_CHILD(p, path, solve);
        }
#endif
    } prm;

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    ilu0_chow_patel(
            const Matrix &A, const params &prm,
            const typename Backend::params &bprm)
      : prm(prm)
    {
        typedef typename backend::builtin<value_type, col_type, ptr_type>::matrix build_matrix;
        const ptrdiff_t n = static_cast<ptrdiff_t>(backend::rows(A));

        // -----------------------------------------------------------------
        // 1. Diagonal scaling.
        //    Row-scale A so that the diagonal is the identity:
        //      A' = diag(1/a_ii) * A.
        //    Store inv_diag[i] = 1/a_ii for unscaling later.
        // -----------------------------------------------------------------
        std::vector<value_type> inv_diag(n);
        for (ptrdiff_t i = 0; i < n; ++i) {
            value_type d = math::zero<value_type>();
            for (ptr_type j = A.ptr[i]; j < A.ptr[i + 1]; ++j) {
                if (A.col[j] == i) {
                    d = A.val[j];
                    break;
                }
            }
            precondition(!math::is_zero(d),
                         "Zero diagonal in ILU0 Chow-Patel scaling");
            inv_diag[i] = math::inverse(d);
        }

        // -----------------------------------------------------------------
        // 2. Count nonzeros per row and build prefix-sum Lptr / Uptr.
        //    Doing this in a separate pass lets the fill in step 3 run
        //    fully in parallel (each row writes to a known, disjoint range).
        //    L is strictly lower triangular (unit diagonal, not stored).
        //    U has its diagonal stored separately in Udiag.
        // -----------------------------------------------------------------
        std::vector<ptr_type> Lptr(n + 1, 0);
        std::vector<ptr_type> Uptr(n + 1, 0);

        for (ptrdiff_t i = 0; i < n; ++i) {
            for (ptr_type j = A.ptr[i]; j < A.ptr[i + 1]; ++j) {
                ptrdiff_t c = A.col[j];
                if      (c < i) ++Lptr[i + 1];
                else if (c > i) ++Uptr[i + 1];
            }
        }
        for (ptrdiff_t i = 1; i <= n; ++i) {
            Lptr[i] += Lptr[i - 1];
            Uptr[i] += Uptr[i - 1];
        }

        const size_t Lnz = static_cast<size_t>(Lptr[n]);
        const size_t Unz = static_cast<size_t>(Uptr[n]);

        // Use numa_vector (no zero-fill on allocation) so that the parallel
        // first-touch in step 3 below places each page on the NUMA node of
        // the thread that will own it during the sweep loops.
        backend::numa_vector<col_type>   Lcol(Lnz, false);
        backend::numa_vector<value_type> Lval(Lnz, false);

        backend::numa_vector<col_type>   Ucol(Unz, false);
        backend::numa_vector<value_type> Uval(Unz, false);

        backend::numa_vector<value_type> Udiag(n, false);

        // -----------------------------------------------------------------
        // 3. Standard initial guess (Section 2.3 of the paper):
        //      L^(0)_{ij} = a'_{ij}      (i > j)
        //      U^(0)_{ij} = a'_{ij}      (i < j)
        //      U^(0)_{ii} = 1            (unit diagonal after scaling)
        //    where a'_{ij} = a_{ij} / a_{ii}.
        //
        //    The loop is parallel so that each thread first-touches its own
        //    slice of Lcol/Lval, Ucol/Uval, and Udiag – matching the access
        //    pattern of the sweep loops in step 5.
        // -----------------------------------------------------------------
#ifdef _OPENMP
#  pragma omp parallel for schedule(guided, 64)
#endif
        for (ptrdiff_t i = 0; i < n; ++i) {
            ptr_type lh = Lptr[i], uh = Uptr[i];
            for (ptr_type j = A.ptr[i]; j < A.ptr[i + 1]; ++j) {
                ptrdiff_t c = A.col[j];
                value_type v = A.val[j] * inv_diag[i]; // row-scaling
                if (c < i) {
                    Lcol[lh] = static_cast<col_type>(c);
                    Lval[lh] = v;
                    ++lh;
                } else if (c > i) {
                    Ucol[uh] = static_cast<col_type>(c);
                    Uval[uh] = v;
                    ++uh;
                }
            }
            Udiag[i] = math::identity<value_type>();
        }

        // Keep a copy of the scaled matrix entries as RHS of the fixed-point
        // equations (these do not change during sweeps).  Parallel copy so
        // that a_L / a_U are first-touched on the same nodes as Lval / Uval.
        backend::numa_vector<value_type> a_L(Lnz, false);
        backend::numa_vector<value_type> a_U(Unz, false);
#ifdef _OPENMP
#  pragma omp parallel for schedule(guided, 64)
#endif
        for (ptrdiff_t i = 0; i < n; ++i) {
            for (ptr_type j = Lptr[i]; j < Lptr[i + 1]; ++j) a_L[j] = Lval[j];
            for (ptr_type j = Uptr[i]; j < Uptr[i + 1]; ++j) a_U[j] = Uval[j];
        }

        // -----------------------------------------------------------------
        // 4. Build CSC representation of U for efficient column access.
        //    The sparse-sparse inner products require rows of L (CSR) and
        //    columns of U (CSC).
        //
        //    Structure (Ucptr, Ucrow) is built sequentially once; values
        //    (Ucval) are written by a separate parallel loop so that the
        //    first-touch for Ucval happens on the threads that will read it
        //    in rebuild_ucsc_values (which uses the same row-parallel schedule).
        // -----------------------------------------------------------------
        std::vector<ptr_type>    Ucptr(n + 1, 0);
        std::vector<col_type>    Ucrow(Unz);
        backend::numa_vector<value_type> Ucval(Unz, false);

        // build_ucsc: rebuild full CSC (structure + values) from current Uval.
        // Structure is computed sequentially; values are then filled in
        // parallel (row-parallel over i) to establish NUMA locality.
        auto build_ucsc = [&]() {
            // --- structural part (sequential) ---
            std::fill(Ucptr.begin(), Ucptr.end(), ptr_type(0));
            for (size_t k = 0; k < Unz; ++k)
                ++Ucptr[Ucol[k] + 1];
            for (ptrdiff_t j = 1; j <= n; ++j)
                Ucptr[j] += Ucptr[j - 1];

            // Build Ucrow (maps CSC position → original row).
            std::vector<ptr_type> pos(Ucptr.begin(), Ucptr.end());
            for (ptrdiff_t i = 0; i < n; ++i)
                for (ptr_type k = Uptr[i]; k < Uptr[i + 1]; ++k)
                    Ucrow[pos[Ucol[k]]++] = static_cast<col_type>(i);

            // --- value part (parallel first-touch) ---
#ifdef _OPENMP
#  pragma omp parallel for schedule(guided, 64)
#endif
            for (ptrdiff_t i = 0; i < n; ++i) {
                for (ptr_type k = Uptr[i]; k < Uptr[i + 1]; ++k) {
                    ptrdiff_t c = Ucol[k];
                    for (ptr_type p = Ucptr[c]; p < Ucptr[c + 1]; ++p) {
                        if (Ucrow[p] == static_cast<col_type>(i)) {
                            Ucval[p] = Uval[k];
                            break;
                        }
                    }
                }
            }
        };

        // Parallel CSC rebuild: only the value copy needs updating each
        // sweep; the structural indexing (Ucptr, Ucrow) stays the same
        // once built.
        auto rebuild_ucsc_values = [&]() {
#ifdef _OPENMP
#  pragma omp parallel for schedule(guided, 64)
#endif
            for (ptrdiff_t i = 0; i < n; ++i) {
                for (ptr_type k = Uptr[i]; k < Uptr[i + 1]; ++k) {
                    ptrdiff_t c = Ucol[k];
                    // Find position: walk Ucptr[c]..Ucptr[c+1] for row i.
                    for (ptr_type p = Ucptr[c]; p < Ucptr[c + 1]; ++p) {
                        if (Ucrow[p] == static_cast<col_type>(i)) {
                            Ucval[p] = Uval[k];
                            break;
                        }
                    }
                }
            }
        };

        build_ucsc();

        // -----------------------------------------------------------------
        // 5. Fixed-point iteration sweeps (Algorithm 2 from the paper).
        //
        //    For each (i,j) in S, in parallel:
        //      if i > j:  l_ij = (a_ij - sum_{k<j} l_ik u_kj) / u_jj
        //      if i = j:  u_ii = a_ii - sum_{k<i} l_ik u_ki  [a_ii=1]
        //      if i < j:  u_ij = a_ij - sum_{k<i} l_ik u_kj
        //
        //    The inner products are sparse-sparse dot products between
        //    row i of L (CSR) and column j of U (CSC).
        // -----------------------------------------------------------------
        for (int sweep = 0; sweep < prm.sweeps; ++sweep) {
            // Synchronize CSC copy of U with (possibly updated) CSR values.
            // First sweep builds full CSC structure; subsequent sweeps only
            // update the values (structure is invariant).
            if (sweep == 0)
                build_ucsc();
            else
                rebuild_ucsc_values();

            // Update L entries: l_ij for i > j
#ifdef _OPENMP
#  pragma omp parallel for schedule(guided, 64)
#endif
            for (ptrdiff_t i = 0; i < n; ++i) {
                for (ptr_type jj = Lptr[i]; jj < Lptr[i + 1]; ++jj) {
                    ptrdiff_t j = Lcol[jj];

                    // Sparse dot: row i of L  dot  col j of U  (indices < j)
                    value_type s = math::zero<value_type>();
                    ptr_type lp = Lptr[i],  le = Lptr[i + 1];
                    ptr_type up = Ucptr[j], ue = Ucptr[j + 1];

                    while (lp < le && up < ue) {
                        col_type lc = Lcol[lp];
                        col_type uc = Ucrow[up];
                        if (lc >= j || uc >= j) break;
                        if      (lc < uc) { ++lp; }
                        else if (lc > uc) { ++up; }
                        else { s += Lval[lp] * Ucval[up]; ++lp; ++up; }
                    }

                    Lval[jj] = (a_L[jj] - s) * math::inverse(Udiag[j]);
                }
            }

            // Update U diagonal and upper entries
#ifdef _OPENMP
#  pragma omp parallel for schedule(guided, 64)
#endif
            for (ptrdiff_t i = 0; i < n; ++i) {
                // u_ii = 1 - sum_{k<i} l_ik u_ki
                {
                    value_type s = math::zero<value_type>();
                    ptr_type lp = Lptr[i],  le = Lptr[i + 1];
                    ptr_type up = Ucptr[i], ue = Ucptr[i + 1];

                    while (lp < le && up < ue) {
                        col_type lc = Lcol[lp];
                        col_type uc = Ucrow[up];
                        if (lc >= i || uc >= i) break;
                        if      (lc < uc) { ++lp; }
                        else if (lc > uc) { ++up; }
                        else { s += Lval[lp] * Ucval[up]; ++lp; ++up; }
                    }
                    Udiag[i] = math::identity<value_type>() - s;
                }

                // u_ij for j > i
                for (ptr_type jj = Uptr[i]; jj < Uptr[i + 1]; ++jj) {
                    ptrdiff_t j = Ucol[jj];

                    value_type s = math::zero<value_type>();
                    ptr_type lp = Lptr[i],  le = Lptr[i + 1];
                    ptr_type up = Ucptr[j], ue = Ucptr[j + 1];

                    while (lp < le && up < ue) {
                        col_type lc = Lcol[lp];
                        col_type uc = Ucrow[up];
                        if (lc >= i || uc >= i) break;
                        if      (lc < uc) { ++lp; }
                        else if (lc > uc) { ++up; }
                        else { s += Lval[lp] * Ucval[up]; ++lp; ++up; }
                    }

                    Uval[jj] = a_U[jj] - s;
                }
            }
        }

        // -----------------------------------------------------------------
        // 6. Un-scale and build output factors for ilu_solve.
        //
        //    We computed L,U for scaled system A' = diag(1/a_ii)*A such that
        //    L*U ≈ A'.  Therefore  diag(a_ii)*L*U ≈ A.
        //
        //    ilu_solve expects factors in the form
        //      M = (I + L_s)(diag(pivot) + U_s)
        //    where L_s is strictly lower triangular, U_s is strictly upper
        //    triangular, and it stores D = inv(pivot).
        //
        //    Decomposing  diag(a_ii)*(I + L_cp)*(diag(Udiag) + U_cp) yields:
        //      L_s_{ij}  = (a_ii / a_jj) * L_cp_{ij}     (i > j)
        //      pivot_i   = a_ii * Udiag_i
        //      U_s_{ij}  = a_ii * U_cp_{ij}               (i < j)
        // -----------------------------------------------------------------
        auto L_out = std::make_shared<build_matrix>();
        auto U_out = std::make_shared<build_matrix>();
        L_out->set_size(n, n); L_out->set_nonzeros(Lnz); L_out->ptr[0] = 0;
        U_out->set_size(n, n); U_out->set_nonzeros(Unz); U_out->ptr[0] = 0;

        auto D_out = std::make_shared<backend::numa_vector<value_type>>(n, false);

#ifdef _OPENMP
#  pragma omp parallel for schedule(guided, 64)
#endif
        for (ptrdiff_t i = 0; i < n; ++i) {
            value_type a_ii = math::inverse(inv_diag[i]);

            for (ptr_type jj = Lptr[i]; jj < Lptr[i + 1]; ++jj) {
                ptrdiff_t j = Lcol[jj];
                value_type a_jj = math::inverse(inv_diag[j]);
                L_out->col[jj] = Lcol[jj];
                L_out->val[jj] = (a_ii * math::inverse(a_jj)) * Lval[jj];
            }
            L_out->ptr[i + 1] = Lptr[i + 1];

            for (ptr_type jj = Uptr[i]; jj < Uptr[i + 1]; ++jj) {
                U_out->col[jj] = Ucol[jj];
                U_out->val[jj] = a_ii * Uval[jj];
            }
            U_out->ptr[i + 1] = Uptr[i + 1];

            (*D_out)[i] = math::inverse(a_ii * Udiag[i]);
        }

        ilu = std::make_shared<ilu_solve>(L_out, U_out, D_out, prm.solve, bprm);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        ilu->solve(tmp);
        backend::axpby(prm.damping, tmp, math::identity<scalar_type>(), x);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        ilu->solve(tmp);
        backend::axpby(prm.damping, tmp, math::identity<scalar_type>(), x);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply
    template <class Matrix, class VectorRHS, class VectorX>
    void apply(const Matrix&, const VectorRHS &rhs, VectorX &x) const
    {
        backend::copy(rhs, x);
        ilu->solve(x);
    }

    size_t bytes() const {
        return ilu->bytes();
    }

    private:
        std::shared_ptr<ilu_solve> ilu;
};

} // namespace relaxation
} // namespace amgcl

#endif