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
 * \author Vicente Mataix Ferrandiz <vicente.mataix-ferrandiz@siemens.com>
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
 * The matrix is scaled to have a unit diagonal before factorization,
 * as recommended by the paper for convergence.  Two scaling strategies
 * are available:
 *   - Row scaling:       A' = diag(1/a_ii) * A             (default)
 *   - Symmetric scaling: A' = diag(1/sqrt|a_ii|)*A*diag(1/sqrt|a_ii|)
 *
 * Symmetric scaling (DAD) is recommended by the paper (Section 2.3) and
 * generally gives better convergence because it balances both rows and
 * columns.  Row scaling is kept as the default for backward compatibility
 * and because it can be more robust when the diagonal entries have very
 * large magnitude variation.
 *
 * An underrelaxation factor omega can be used to stabilize the
 * fixed-point iteration for non-diagonally-dominant matrices.
 */

#include <vector>
#include <cmath>
#include <limits>
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

        /// Underrelaxation factor for the fixed-point sweeps (0 < omega <= 1).
        /// Values less than 1 dampen the update and can stabilize convergence
        /// for non-diagonally-dominant matrices at the cost of slower
        /// convergence for well-conditioned ones.
        scalar_type omega;

        /// Use symmetric diagonal scaling (DAD) instead of row scaling (DA).
        /// Symmetric scaling often improves convergence of the iterative
        /// factorization (see Section 2.3 of Chow & Patel 2015).
        bool symmetric_scaling;

        /// Parameters for sparse triangular system solver.
        typename ilu_solve::params solve;

        params() : damping(1), sweeps(5), omega(0.8), symmetric_scaling(false) {}

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, damping)
            , AMGCL_PARAMS_IMPORT_VALUE(p, sweeps)
            , AMGCL_PARAMS_IMPORT_VALUE(p, omega)
            , AMGCL_PARAMS_IMPORT_VALUE(p, symmetric_scaling)
            , AMGCL_PARAMS_IMPORT_CHILD(p, solve)
        {
            check_params(p, {"damping", "sweeps", "omega", "symmetric_scaling", "solve"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, damping);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, sweeps);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, omega);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, symmetric_scaling);
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
        //    Compute per-row scale factors so the scaled matrix has unit
        //    diagonal.  Two strategies are supported:
        //
        //    Row scaling (default):
        //      A' = diag(d_i) * A,  where d_i = 1/a_ii.
        //
        //    Symmetric scaling (Section 2.3 of the paper):
        //      A' = diag(d_i) * A * diag(d_i),  where d_i = 1/sqrt(|a_ii|).
        //      This balances both rows and columns and generally gives
        //      better convergence when the matrix is far from diagonally
        //      dominant.
        //
        //    In both cases, scale[i] stores the left (row) scale factor.
        //    For symmetric scaling col_scale[j] stores the column factor
        //    (same array – it equals scale[j]).
        // -----------------------------------------------------------------
        auto inv_diag_ptr = backend::diagonal(A, /*invert=*/true);
        auto &inv_diag = *inv_diag_ptr;  // inv_diag[i] = 1/a_ii

        // scale[i] is the left scale factor d_i for row i.
        backend::numa_vector<value_type> scale(n, false);
#ifdef _OPENMP
#  pragma omp parallel for schedule(guided, 64)
#endif
        for (ptrdiff_t i = 0; i < n; ++i) {
            if (prm.symmetric_scaling) {
                // d_i = 1 / sqrt(|a_ii|)
                scalar_type abs_aii = math::norm(math::inverse(inv_diag[i]));
                scale[i] = (abs_aii > std::numeric_limits<scalar_type>::min())
                    ? static_cast<scalar_type>(1.0 / std::sqrt(abs_aii)) * math::identity<value_type>()
                    : inv_diag[i];
            } else {
                // d_i = 1 / a_ii  (row scaling)
                scale[i] = inv_diag[i];
            }
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
        //
        //    Row scaling:       a'_{ij} = a_{ij} * scale[i]
        //    Symmetric scaling: a'_{ij} = a_{ij} * scale[i] * scale[j]
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
                // Scale: row scaling  a'_ij = a_ij * scale[i]
                //        sym scaling  a'_ij = a_ij * scale[i] * scale[c]
                value_type v = A.val[j] * scale[i];
                if (prm.symmetric_scaling) v = v * scale[c];
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
        //
        //    When omega < 1, each update is under-relaxed:
        //      x_new = omega * G(x) + (1 - omega) * x_old
        //    This can stabilize convergence for non-diagonally-dominant
        //    matrices at the cost of slower convergence.
        // -----------------------------------------------------------------
        const scalar_type omega   = prm.omega;
        const scalar_type one_m_w = static_cast<scalar_type>(1) - omega;

        // Minimum acceptable pivot norm (relative to identity).
        // Using sqrt(eps) instead of eps gives much more robust behavior
        // for ill-conditioned matrices.
        const scalar_type pivot_tol =
            std::sqrt(std::numeric_limits<scalar_type>::epsilon())
            * math::norm(math::identity<value_type>());

        for (int sweep = 0; sweep < prm.sweeps; ++sweep) {
            // Synchronize CSC copy of U with (possibly updated) CSR values.
            // First sweep builds the full CSC structure + values; subsequent
            // sweeps only update the values (structure is invariant).
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

                    // Guard against zero pivot.
                    value_type pivot = Udiag[j];
                    if (math::norm(pivot) < pivot_tol)
                        pivot = pivot_tol * math::identity<value_type>();

                    value_type l_new = (a_L[jj] - s) * math::inverse(pivot);

                    // Protect against NaN/Inf from numerical blow-up.
                    if (!std::isfinite(math::norm(l_new))) l_new = a_L[jj];

                    // Under-relaxation.
                    Lval[jj] = omega * l_new + one_m_w * Lval[jj];
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
                    value_type u_diag_new = math::identity<value_type>() - s;
                    if (!std::isfinite(math::norm(u_diag_new)))
                        u_diag_new = math::identity<value_type>();
                    Udiag[i] = omega * u_diag_new + one_m_w * Udiag[i];
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

                    value_type u_new = a_U[jj] - s;
                    if (!std::isfinite(math::norm(u_new))) u_new = a_U[jj];
                    Uval[jj] = omega * u_new + one_m_w * Uval[jj];
                }
            }
        }

        // -----------------------------------------------------------------
        // 6. Un-scale and build output factors for ilu_solve.
        //
        //    ilu_solve expects factors in the form
        //      M = (I + L_s)(diag(pivot) + U_s)
        //    where L_s is strictly lower triangular, U_s is strictly upper
        //    triangular, and it stores D = inv(pivot).
        //
        //    Row scaling:  A' = diag(1/a_ii)*A,  L'U' ≈ A'.
        //      (I + L_s)(pivot + U_s) = diag(a_ii)*(I + L')(diag(U'd) + U')
        //        L_s_{ij}  = (a_ii / a_jj) * L'_{ij}
        //        pivot_i   = a_ii * U'diag_i
        //        U_s_{ij}  = a_ii * U'_{ij}
        //
        //    Symmetric scaling: A' = D*A*D,  D = diag(1/sqrt|a_ii|).
        //      (I + L_s)(pivot + U_s) = D^{-1}(I + L')(diag(U'd) + U')D^{-1}
        //        L_s_{ij}  = sqrt(|a_ii| / |a_jj|) * L'_{ij}
        //        pivot_i   = |a_ii| * U'diag_i
        //        U_s_{ij}  = sqrt(|a_ii| * |a_jj|) * U'_{ij}
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
            // Compute un-scaling factors for row i.
            value_type a_ii = math::inverse(inv_diag[i]);
            scalar_type abs_aii = math::norm(a_ii);

            for (ptr_type jj = Lptr[i]; jj < Lptr[i + 1]; ++jj) {
                ptrdiff_t j = Lcol[jj];
                L_out->col[jj] = Lcol[jj];
                if (prm.symmetric_scaling) {
                    // L_s_{ij} = sqrt(|a_ii| / |a_jj|) * L'_{ij}
                    scalar_type abs_ajj = math::norm(math::inverse(inv_diag[j]));
                    scalar_type ratio = (abs_ajj > std::numeric_limits<scalar_type>::min())
                        ? std::sqrt(abs_aii / abs_ajj)
                        : static_cast<scalar_type>(1);
                    L_out->val[jj] = ratio * Lval[jj];
                } else {
                    // L_s_{ij} = (a_ii / a_jj) * L'_{ij}
                    value_type a_jj = math::inverse(inv_diag[j]);
                    L_out->val[jj] = (a_ii * math::inverse(a_jj)) * Lval[jj];
                }
            }
            L_out->ptr[i + 1] = Lptr[i + 1];

            for (ptr_type jj = Uptr[i]; jj < Uptr[i + 1]; ++jj) {
                U_out->col[jj] = Ucol[jj];
                if (prm.symmetric_scaling) {
                    // U_s_{ij} = sqrt(|a_ii| * |a_jj|) * U'_{ij}
                    ptrdiff_t j = Ucol[jj];
                    scalar_type abs_ajj = math::norm(math::inverse(inv_diag[j]));
                    scalar_type factor = std::sqrt(abs_aii * abs_ajj);
                    U_out->val[jj] = factor * Uval[jj];
                } else {
                    // U_s_{ij} = a_ii * U'_{ij}
                    U_out->val[jj] = a_ii * Uval[jj];
                }
            }
            U_out->ptr[i + 1] = Uptr[i + 1];

            // Compute pivot and guard against zero.
            value_type pivot;
            if (prm.symmetric_scaling) {
                // pivot_i = |a_ii| * U'diag_i
                pivot = abs_aii * Udiag[i];
            } else {
                // pivot_i = a_ii * U'diag_i
                pivot = a_ii * Udiag[i];
            }
            if (math::norm(pivot) < pivot_tol)
                pivot = pivot_tol * math::identity<value_type>();
            (*D_out)[i] = math::inverse(pivot);
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