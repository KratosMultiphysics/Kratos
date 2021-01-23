#ifndef AMGCL_COARSENING_SMOOTHED_AGGREGATION_HPP
#define AMGCL_COARSENING_SMOOTHED_AGGREGATION_HPP

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
 * \file   amgcl/coarsening/smoothed_aggregation.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Smoothed aggregation coarsening scheme.
 */

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <tuple>
#include <memory>

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
template <class Backend>
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

        // Estimate the matrix spectral radius.
        // This usually improves convergence rate and results in faster solves,
        // but costs some time during setup.
        bool estimate_spectral_radius;

        // Number of power iterations to apply for the spectral radius
        // estimation. Use Gershgorin disk theorem when power_iters = 0.
        int power_iters;

        params() : relax(1.0f), estimate_spectral_radius(false), power_iters(0) { }

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_CHILD(p, aggr),
              AMGCL_PARAMS_IMPORT_CHILD(p, nullspace),
              AMGCL_PARAMS_IMPORT_VALUE(p, relax),
              AMGCL_PARAMS_IMPORT_VALUE(p, estimate_spectral_radius),
              AMGCL_PARAMS_IMPORT_VALUE(p, power_iters)
        {
            check_params(p, {"aggr", "nullspace", "relax", "estimate_spectral_radius", "power_iters"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_CHILD(p, path, aggr);
            AMGCL_PARAMS_EXPORT_CHILD(p, path, nullspace);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, relax);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, estimate_spectral_radius);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, power_iters);
        }
#endif
    } prm;

    smoothed_aggregation(const params &prm = params()) : prm(prm) {}

    /// \copydoc amgcl::coarsening::aggregation::transfer_operators
    template <class Matrix>
    std::tuple< std::shared_ptr<Matrix>, std::shared_ptr<Matrix> >
    transfer_operators(const Matrix &A) {
        typedef typename backend::value_type<Matrix>::type value_type;
        typedef typename math::scalar_of<value_type>::type scalar_type;

        const size_t n = rows(A);

        AMGCL_TIC("aggregates");
        Aggregates aggr(A, prm.aggr, prm.nullspace.cols);
        prm.aggr.eps_strong *= 0.5;
        AMGCL_TOC("aggregates");

        auto P_tent = tentative_prolongation<Matrix>(
                n, aggr.count, aggr.id, prm.nullspace, prm.aggr.block_size
                );

        auto P = std::make_shared<Matrix>();
        P->set_size(rows(*P_tent), cols(*P_tent), true);

        scalar_type omega = prm.relax;
        if (prm.estimate_spectral_radius) {
            omega *= static_cast<scalar_type>(4.0/3) / backend::spectral_radius<true>(A, prm.power_iters);
        } else {
            omega *= static_cast<scalar_type>(2.0/3);
        }

        AMGCL_TIC("smoothing");
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

        P->scan_row_sizes();
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
                    if (A.col[j] == i || !aggr.strong_connection[j])
                        dia += A.val[j];
                }
                dia = -omega * math::inverse(dia);

                ptrdiff_t row_beg = P->ptr[i];
                ptrdiff_t row_end = row_beg;
                for(ptrdiff_t ja = A.ptr[i], ea = A.ptr[i + 1]; ja < ea; ++ja) {
                    ptrdiff_t ca = A.col[ja];

                    // Skip weak off-diagonal connections.
                    if (ca != i && !aggr.strong_connection[ja]) continue;

                    value_type va = (ca == i)
                        ? static_cast<value_type>(static_cast<scalar_type>(1 - omega) * math::identity<value_type>())
                        : dia * A.val[ja];

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
        AMGCL_TOC("smoothing");

        if (prm.nullspace.cols > 0)
            prm.aggr.block_size = prm.nullspace.cols;

        return std::make_tuple(P, transpose(*P));
    }

    /// \copydoc amgcl::coarsening::aggregation::coarse_operator
    template <class Matrix>
    std::shared_ptr<Matrix>
    coarse_operator(const Matrix &A, const Matrix &P, const Matrix &R) const {
        return detail::galerkin(A, P, R);
    }
};

} // namespace coarsening
} // namespace amgcl

#endif
