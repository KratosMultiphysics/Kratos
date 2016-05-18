#ifndef AMGCL_COARSENING_SMOOTHED_AGGREGATION_HPP
#define AMGCL_COARSENING_SMOOTHED_AGGREGATION_HPP

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
 * \file   amgcl/coarsening/smoothed_aggregation.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Smoothed aggregation coarsening scheme.
 */

#include <boost/typeof/typeof.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>

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

        /// Relaxation factor \f$\omega\f$.
        /**
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

        params() : relax(0.666f) { }

        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_CHILD(p, aggr),
              AMGCL_PARAMS_IMPORT_CHILD(p, nullspace),
              AMGCL_PARAMS_IMPORT_VALUE(p, relax)
        {
            AMGCL_PARAMS_CHECK(p, (aggr)(nullspace)(relax));
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_CHILD(p, path, aggr);
            AMGCL_PARAMS_EXPORT_CHILD(p, path, nullspace);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, relax);
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

        BOOST_AUTO(Aptr, A.ptr_data());
        BOOST_AUTO(Acol, A.col_data());
        BOOST_AUTO(Aval, A.val_data());

        TIC("aggregates");
        Aggregates aggr(A, prm.aggr, prm.nullspace.cols);
        prm.aggr.eps_strong *= 0.5;
        TOC("aggregates");

        TIC("interpolation");
        boost::shared_ptr<Matrix> P_tent = tentative_prolongation<Matrix>(
                n, aggr.count, aggr.id, prm.nullspace, prm.aggr.block_size
                );

        boost::shared_ptr<Matrix> P = boost::make_shared<Matrix>();
        P->nrows = rows(*P_tent);
        P->ncols = cols(*P_tent);

        P->ptr.resize(n + 1, 0);

#pragma omp parallel
        {
            std::vector<ptrdiff_t> marker(P->ncols, -1);

#ifdef _OPENMP
            int nt  = omp_get_num_threads();
            int tid = omp_get_thread_num();

            ptrdiff_t chunk_size  = (n + nt - 1) / nt;
            ptrdiff_t chunk_start = tid * chunk_size;
            ptrdiff_t chunk_end   = std::min<ptrdiff_t>(n, chunk_start + chunk_size);
#else
            ptrdiff_t chunk_start = 0;
            ptrdiff_t chunk_end   = n;
#endif

            // Count number of entries in P.
            for(ptrdiff_t i = chunk_start; i < chunk_end; ++i) {
                for(ptrdiff_t ja = Aptr[i], ea = Aptr[i+1]; ja < ea; ++ja) {
                    ptrdiff_t ca = Acol[ja];

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

            boost::fill(marker, -1);

#pragma omp barrier
#pragma omp single
            {
                boost::partial_sum(P->ptr, P->ptr.begin());
                P->col.resize(P->ptr.back());
                P->val.resize(P->ptr.back());
            }

            // Fill the interpolation matrix.
            for(ptrdiff_t i = chunk_start; i < chunk_end; ++i) {

                // Diagonal of the filtered matrix is the original matrix
                // diagonal minus its weak connections.
                value_type dia = math::zero<value_type>();
                for(ptrdiff_t j = Aptr[i], e = Aptr[i+1]; j < e; ++j) {
                    if (Acol[j] == i)
                        dia += Aval[j];
                    else if (!aggr.strong_connection[j])
                        dia -= Aval[j];
                }
                dia = math::inverse(dia);

                ptrdiff_t row_beg = P->ptr[i];
                ptrdiff_t row_end = row_beg;
                for(ptrdiff_t ja = Aptr[i], ea = Aptr[i + 1]; ja < ea; ++ja) {
                    ptrdiff_t ca = Acol[ja];

                    // Skip weak off-diagonal connections.
                    if (ca != i && !aggr.strong_connection[ja]) continue;

                    value_type va = (ca == i)
                        ? static_cast<value_type>(static_cast<scalar_type>(1 - prm.relax) * math::identity<value_type>())
                        : static_cast<value_type>(static_cast<scalar_type>(-prm.relax) * dia * Aval[ja]);

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
        TOC("interpolation");

        boost::shared_ptr<Matrix> R = boost::make_shared<Matrix>();
        *R = transpose(*P);

        if (prm.nullspace.cols > 0)
            prm.aggr.block_size = prm.nullspace.cols;

        return boost::make_tuple(P, R);
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
};

} // namespace coarsening
} // namespace amgcl

#endif
