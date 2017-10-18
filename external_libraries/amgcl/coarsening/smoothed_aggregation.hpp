#ifndef AMGCL_COARSENING_SMOOTHED_AGGREGATION_HPP
#define AMGCL_COARSENING_SMOOTHED_AGGREGATION_HPP

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
 * \file   amgcl/coarsening/smoothed_aggregation.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Smoothed aggregation coarsening scheme.
 */

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/detail/galerkin.hpp>
#include <amgcl/coarsening/pointwise_aggregates.hpp>
#include <amgcl/coarsening/tentative_prolongation.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace coarsening {

/// Smoothed aggregation coarsening.
/**
 * \ingroup coarsening
 * \sa \cite Vanek1996
 */
struct smoothed_aggregation {
    typedef pointwise_aggregates Aggregates;

    /// Coarsening parameters
    struct params {
        /// Aggregation parameters.
        Aggregates::params aggr;

        /// Near nullspace parameters.
        nullspace_params nullspace;

        /// Relaxation factor.
        /**
         * Used as a scaling for the damping factor omega.
         * When estimate_spectral_radius is set, then
         *   omega = relax * (4/3) / rho.
         * Otherwise
         *   omega = relax * (2/3).
         *
         * Piecewise constant prolongation \f$\tilde P\f$ from non-smoothed
         * aggregation is improved by a smoothing to get the final prolongation
         * matrix \f$P\f$. Simple Jacobi smoother is used here, giving the
         * prolongation matrix
         * \f[P = \left( I - \omega D^{-1} A^F \right) \tilde P.\f]
         * Here \f$A^F = (a_{ij}^F)\f$ is the filtered matrix given by
         * \f[
         * a_{ij}^F =
         * \begin{cases}
         * a_{ij} \quad \text{if} \; j \in N_i\\
         * 0 \quad \text{otherwise}
         * \end{cases}, \quad \text{if}\; i \neq j,
         * \quad a_{ii}^F = a_{ii} - \sum\limits_{j=1,j\neq i}^n
         * \left(a_{ij} - a_{ij}^F \right),
         * \f]
         * where \f$N_i\f$ is the set of variables, strongly coupled to
         * variable \f$i\f$, and \f$D\f$ denotes the diagonal of \f$A^F\f$.
         */
        float relax;

        // Use power iterations to estimate the matrix spectral radius.
        // This usually improves convergence rate and results in faster solves,
        // but costs some time during setup.
        bool estimate_spectral_radius;

        // Number of power iterations to apply for the spectral radius
        // estimation.
        int power_iters;

        params() : relax(1.0f), estimate_spectral_radius(true), power_iters(5) { }

        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_CHILD(p, aggr),
              AMGCL_PARAMS_IMPORT_CHILD(p, nullspace),
              AMGCL_PARAMS_IMPORT_VALUE(p, relax),
              AMGCL_PARAMS_IMPORT_VALUE(p, estimate_spectral_radius),
              AMGCL_PARAMS_IMPORT_VALUE(p, power_iters)
        {
            AMGCL_PARAMS_CHECK(p, (aggr)(nullspace)(relax)(estimate_spectral_radius)(power_iters));
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_CHILD(p, path, aggr);
            AMGCL_PARAMS_EXPORT_CHILD(p, path, nullspace);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, relax);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, estimate_spectral_radius);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, power_iters);
        }
    };

    /// \copydoc amgcl::coarsening::aggregation::transfer_operators
    template <class Matrix>
    static boost::tuple< boost::shared_ptr<Matrix>, boost::shared_ptr<Matrix> >
    transfer_operators(const Matrix &A, params &prm)
    {
        typedef typename backend::value_type<Matrix>::type value_type;
        typedef typename math::scalar_of<value_type>::type scalar_type;

        const size_t n = rows(A);

        AMGCL_TIC("aggregates");
        Aggregates aggr(A, prm.aggr, prm.nullspace.cols);
        prm.aggr.eps_strong *= 0.5;
        AMGCL_TOC("aggregates");

        AMGCL_TIC("interpolation");
        boost::shared_ptr<Matrix> P_tent = tentative_prolongation<Matrix>(
                n, aggr.count, aggr.id, prm.nullspace, prm.aggr.block_size
                );

        boost::shared_ptr<Matrix> P = boost::make_shared<Matrix>();
        P->set_size(rows(*P_tent), cols(*P_tent), true);

        scalar_type omega = prm.relax;
        if (prm.estimate_spectral_radius) {
            omega *= static_cast<scalar_type>(4.0/3) / spectral_radius(A, prm.power_iters);
        } else {
            omega *= static_cast<scalar_type>(2.0/3);
        }

#pragma omp parallel
        {
            std::vector<ptrdiff_t> marker(P->ncols, -1);

            // Count number of entries in P.
#pragma omp for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                for(ptrdiff_t ja = A.ptr[i], ea = A.ptr[i+1]; ja < ea; ++ja) {
                    ptrdiff_t ca = A.col[ja];

                    // Skip weak off-diagonal connections.
                    if (ca != i && !aggr.strong_connection[ja])
                        continue;

                    for(ptrdiff_t jp = P_tent->ptr[ca], ep = P_tent->ptr[ca+1]; jp < ep; ++jp) {
                        ptrdiff_t cp = P_tent->col[jp];

                        if (marker[cp] != i) {
                            marker[cp] = i;
                            ++( P->ptr[i + 1] );
                        }
                    }
                }
            }
        }

        std::partial_sum(P->ptr, P->ptr + n + 1, P->ptr);
        P->set_nonzeros();

#pragma omp parallel
        {
            std::vector<ptrdiff_t> marker(P->ncols, -1);

            // Fill the interpolation matrix.
#pragma omp for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {

                // Diagonal of the filtered matrix is the original matrix
                // diagonal minus its weak connections.
                value_type dia = math::zero<value_type>();
                for(ptrdiff_t j = A.ptr[i], e = A.ptr[i+1]; j < e; ++j) {
                    if (A.col[j] == i)
                        dia += A.val[j];
                    else if (!aggr.strong_connection[j])
                        dia -= A.val[j];
                }
                dia = math::inverse(dia);

                ptrdiff_t row_beg = P->ptr[i];
                ptrdiff_t row_end = row_beg;
                for(ptrdiff_t ja = A.ptr[i], ea = A.ptr[i + 1]; ja < ea; ++ja) {
                    ptrdiff_t ca = A.col[ja];

                    // Skip weak off-diagonal connections.
                    if (ca != i && !aggr.strong_connection[ja]) continue;

                    value_type va = (ca == i)
                        ? static_cast<value_type>(static_cast<scalar_type>(1 - omega) * math::identity<value_type>())
                        : static_cast<value_type>(static_cast<scalar_type>(-omega) * dia * A.val[ja]);

                    for(ptrdiff_t jp = P_tent->ptr[ca], ep = P_tent->ptr[ca+1]; jp < ep; ++jp) {
                        ptrdiff_t cp = P_tent->col[jp];
                        value_type vp = P_tent->val[jp];

                        if (marker[cp] < row_beg) {
                            marker[cp] = row_end;
                            P->col[row_end] = cp;
                            P->val[row_end] = va * vp;
                            ++row_end;
                        } else {
                            P->val[ marker[cp] ] += va * vp;
                        }
                    }
                }
            }
        }
        AMGCL_TOC("interpolation");

        if (prm.nullspace.cols > 0)
            prm.aggr.block_size = prm.nullspace.cols;

        return boost::make_tuple(P, transpose(*P));
    }

    /// \copydoc amgcl::coarsening::aggregation::coarse_operator
    template <class Matrix>
    static boost::shared_ptr<Matrix>
    coarse_operator(
            const Matrix &A,
            const Matrix &P,
            const Matrix &R,
            const params&
            )
    {
        return detail::galerkin(A, P, R);
    }

    // Uses power iteration to estimate spectral readius of the matrix,
    // scaled by its inverse diagonal.
    template <class Matrix>
    static
    typename math::scalar_of<typename backend::value_type<Matrix>::type>::type
    spectral_radius(const Matrix &A, int power_iters)
    {
        typedef typename backend::value_type<Matrix>::type   value_type;
        typedef typename math::rhs_of<value_type>::type      rhs_type;
        typedef typename math::scalar_of<value_type>::type   scalar_type;

        const ptrdiff_t n = backend::rows(A);

        backend::numa_vector<value_type> D(n, false);
        backend::numa_vector<rhs_type>   b0(n, false), b1(n, false);

        // Fill the initial vector with random values.
        // Also extract the inverted matrix diagonal values.
        scalar_type b0_norm = 0;
#pragma omp parallel
        {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
#else
            int tid = 0;
#endif
            boost::random::mt11213b rng(tid);
            boost::random::uniform_real_distribution<scalar_type> rnd(-1, 1);

            scalar_type loc_norm = 0;

#pragma omp for nowait
            for(ptrdiff_t i = 0; i < n; ++i) {
                rhs_type v = math::constant<rhs_type>(rnd(rng));

                b0[i] = v;
                loc_norm += math::norm(math::inner_product(v,v));

                for(ptrdiff_t j = A.ptr[i], e = A.ptr[i+1]; j < e; ++j) {
                    if (A.col[j] == i) {
                        D[i] = math::inverse(A.val[j]);
                        break;
                    }
                }
            }

#pragma omp critical
            b0_norm += loc_norm;
        }

        // Normalize b0
        b0_norm = 1 / sqrt(b0_norm);
#pragma omp parallel for
        for(ptrdiff_t i = 0; i < n; ++i) {
            b0[i] = b0_norm * b0[i];
        }

        scalar_type radius = 1;

        for(int iter = 0; iter < power_iters;) {
            // b1 = (D * A) * b0
            // b1_norm = ||b1||
            // radius = <b1,b0>
            scalar_type b1_norm = 0;
            radius = 0;
#pragma omp parallel
            {
                scalar_type loc_norm = 0;
                scalar_type loc_radi = 0;

#pragma omp for nowait
                for(ptrdiff_t i = 0; i < n; ++i) {
                    rhs_type s = math::zero<rhs_type>();

                    for(ptrdiff_t j = A.ptr[i], e = A.ptr[i+1]; j < e; ++j) {
                        s += A.val[j] * b0[A.col[j]];
                    }

                    s = D[i] * s;

                    loc_norm += math::norm(math::inner_product(s, s));
                    loc_radi += math::norm(math::inner_product(s, b0[i]));

                    b1[i] = s;
                }

#pragma omp critical
                {
                    b1_norm += loc_norm;
                    radius  += loc_radi;
                }
            }

            if (++iter < power_iters) {
                // b0 = b1 / b1_norm
                b1_norm = 1 / sqrt(b1_norm);
#pragma omp parallel for
                for(ptrdiff_t i = 0; i < n; ++i) {
                    b0[i] = b1_norm * b1[i];
                }
            }
        }

        return radius < 0 ? static_cast<scalar_type>(2) : radius;
    }
};

} // namespace coarsening
} // namespace amgcl

#endif
