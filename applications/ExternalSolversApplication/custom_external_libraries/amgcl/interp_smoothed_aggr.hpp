#ifndef AMGCL_INTERP_SMOOTHED_AGGR_HPP
#define AMGCL_INTERP_SMOOTHED_AGGR_HPP

/*
The MIT License

Copyright (c) 2012 Denis Demidov <ddemidov@ksu.ru>

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
 * \file   interp_smoothed_aggr.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  Interpolation scheme based on smoothed aggregation.
 */

#include <vector>
#include <algorithm>

#include <boost/typeof/typeof.hpp>

#include <amgcl/spmat.hpp>
#include <amgcl/aggr_connect.hpp>
#include <amgcl/tictoc.hpp>

namespace amgcl {

namespace interp {

/// Interpolation scheme based on smoothed aggregation.
/**
 * See \ref Vanek_1996 "Vanek (1996)"
 *
 * \param aggr_type \ref aggregation "Aggregation scheme".
 *
 * \ingroup interpolation
 */
template <class aggr_type>
struct smoothed_aggregation {

/// Parameters controlling aggregation.
struct params {
    /// Relaxation factor \f$\omega\f$.
    /**
     * See \ref Vanek_1996 "Vanek (1996)".
     * Piecewise constant prolongation \f$\tilde P\f$ from \ref
     * amgcl::interp::aggregation is improved by a smoothing to get the final
     * prolongation matrix \f$P\f$. Simple Jacobi smoother is used here, giving
     * the prolongation matrix
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
     * where \f$N_i\f$ is the set of variables, strongly coupled to variable
     * \f$i\f$, and \f$D\f$ denotes the diagonal of \f$A^F\f$.
     */
    float relax;

    /// Parameter \f$\varepsilon_{str}\f$ defining strong couplings.
    /**
     * Variable \f$i\f$ is defined to be strongly coupled to another variable,
     * \f$j\f$, if \f[|a_{ij}| \geq \varepsilon\sqrt{a_{ii} a_{jj}}\quad
     * \text{with fixed} \quad \varepsilon = \varepsilon_{str} \left(
     * \frac{1}{2} \right)^l,\f]
     * where \f$l\f$ is level number (finest level is 0).
     */
    mutable float eps_strong;

    unsigned dof_per_node;

    params() : relax(2.0f / 3.0f), eps_strong(0.08f), dof_per_node(1) {}
};

/// Constructs coarse level by aggregation.
/**
 * Returns interpolation operator, which is enough to construct system matrix
 * at coarser level.
 *
 * \param A   system matrix.
 * \param prm parameters.
 *
 * \returns interpolation operator.
 */
template < class value_t, class index_t >
static std::pair<
    sparse::matrix<value_t, index_t>,
    sparse::matrix<value_t, index_t>
    >
interp(const sparse::matrix<value_t, index_t> &A, const params &prm) {
    const index_t n = sparse::matrix_rows(A);

    std::vector<char>    S;
    std::vector<index_t> aggr;

    assert(prm.dof_per_node > 0);

    if (prm.dof_per_node == 1) {
        // Scalar system. Nothing fancy.
        TIC("connections");
        aggr::connect(A, prm.eps_strong).swap(S);
        TOC("connections");

        TIC("aggregates");
        aggr_type::aggregates(A, S).swap(aggr);
        TOC("aggregates");
    } else {
        // Non-scalar system.
        // Build reduced matrix, find connections and aggregates with it,
        // restore the vectors to full size.

        BOOST_AUTO(S_aggr, aggr::pointwise_coarsening<aggr_type>(
                    A, prm.eps_strong, prm.dof_per_node));
        S.swap(S_aggr.first);
        aggr.swap(S_aggr.second);
    }

    prm.eps_strong *= 0.5;

    index_t nc = std::max(
            static_cast<index_t>(0),
            *std::max_element(aggr.begin(), aggr.end()) + static_cast<index_t>(1)
            );

    TIC("interpolation");
    std::pair<
        sparse::matrix<value_t, index_t>,
        sparse::matrix<value_t, index_t>
        > PR;

    sparse::matrix<value_t, index_t> &P = PR.first;
    sparse::matrix<value_t, index_t> &R = PR.second;

    P.resize(n, nc);
    std::fill(P.row.begin(), P.row.end(), static_cast<index_t>(0));

    BOOST_AUTO(Arow, sparse::matrix_outer_index(A));
    BOOST_AUTO(Acol, sparse::matrix_inner_index(A));
    BOOST_AUTO(Aval, sparse::matrix_values(A));


#pragma omp parallel
    {
        std::vector<index_t> marker(nc, static_cast<index_t>(-1));

#ifdef _OPENMP
	int nt  = omp_get_num_threads();
	int tid = omp_get_thread_num();

	index_t chunk_size  = (n + nt - 1) / nt;
	index_t chunk_start = tid * chunk_size;
	index_t chunk_end   = std::min(n, chunk_start + chunk_size);
#else
	index_t chunk_start = 0;
	index_t chunk_end   = n;
#endif

        // Count number of entries in P.
        for(index_t i = chunk_start; i < chunk_end; ++i) {
            for(index_t j = Arow[i], e = Arow[i+1]; j < e; ++j) {
                index_t c = Acol[j];

                if (c != i && !S[j]) continue;

                index_t g = aggr[c];

                if (g >= 0 && marker[g] != i) {
                    marker[g] = i;
                    ++P.row[i + 1];
                }
            }
        }

        std::fill(marker.begin(), marker.end(), static_cast<index_t>(-1));

#pragma omp barrier
#pragma omp single
        {
            std::partial_sum(P.row.begin(), P.row.end(), P.row.begin());
            P.reserve(P.row.back());
        }

        // Fill the interpolation matrix.
        for(index_t i = chunk_start; i < chunk_end; ++i) {

            // Diagonal of the filtered matrix is original matrix diagonal minus
            // its weak connections.
            value_t dia = 0;
            for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j) {
                if (Acol[j] == i)
                    dia += Aval[j];
                else if (!S[j])
                    dia -= Aval[j];
            }
            dia = 1 / dia;

            index_t row_beg = P.row[i];
            index_t row_end = row_beg;
            for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j) {
                index_t c = Acol[j];

                // Skip weak couplings, ...
                if (c != i && !S[j]) continue;

                // ... and the ones not in any aggregate.
                index_t g = aggr[c];
                if (g < 0) continue;

                value_t v = (c == i) ? 1 - prm.relax : -prm.relax * dia * Aval[j];

                if (marker[g] < row_beg) {
                    marker[g] = row_end;
                    P.col[row_end] = g;
                    P.val[row_end] = v;
                    ++row_end;
                } else {
                    P.val[marker[g]] += v;
                }
            }
        }
    }
    TOC("interpolation");

    sparse::transpose(P).swap(R);
    return PR;
}

};

} // namespace interp
} // namespace amgcl



#endif
