#ifndef AMGCL_INTERP_SA_EMIN_HPP
#define AMGCL_INTERP_SA_EMIN_HPP

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
 * \file   interp_sa_emin.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  Interpolation scheme based on smoothed aggregation with energy minimization.
 */

#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>

#include <boost/typeof/typeof.hpp>

#include <amgcl/spmat.hpp>
#include <amgcl/aggr_connect.hpp>
#include <amgcl/tictoc.hpp>

namespace amgcl {

namespace interp {

/// Interpolation scheme based on smoothed aggregation with energy minimization.
/**
 * See \ref Sala_2008 "Sala (2008)"
 *
 * \param aggr_type \ref aggregation "Aggregation scheme".
 *
 * \ingroup interpolation
 */
template <class aggr_type>
struct sa_emin {

/// Parameters controlling aggregation.
struct params {
    /// Parameter \f$\varepsilon_{str}\f$ defining strong couplings.
    /**
     * Variable \f$i\f$ is defined to be strongly coupled to another variable,
     * \f$j\f$, if \f[|a_{ij}| \geq \varepsilon\sqrt{a_{ii} a_{jj}}\quad
     * \text{with fixed} \quad \varepsilon = \varepsilon_{str} \left(
     * \frac{1}{2} \right)^l,\f]
     * where \f$l\f$ is level number (finest level is 0).
     */
    mutable float eps_strong;

    /// Number of degrees of freedom (number of unknowns) per grid node.
    /**
     * Should be used with non-scalar systems of equations. Equations are
     * assumed to be ordered by nodes (not by unknowns).
     * \note This is assumed to be constant for all nodes.
     */
    unsigned dof_per_node;

    params() : eps_strong(0.08f), dof_per_node(1) {}
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

    index_t nc;
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

    nc = std::max(
            static_cast<index_t>(0),
            *std::max_element(aggr.begin(), aggr.end()) + static_cast<index_t>(1)
            );

    // Compute smoothed interpolation and restriction operators.
    static std::pair<
        sparse::matrix<value_t, index_t>,
        sparse::matrix<value_t, index_t>
    > PR;

    std::vector<value_t> omega(n);

    TIC("smoothed interpolation");
    smoothed_interpolation(A, S, aggr, nc, omega).swap(PR.first);
    TOC("smoothed interpolation");

    TIC("smoothed restriction");
    smoothed_restriction(A, S, aggr, nc, omega).swap(PR.second);
    TOC("smoothed restriction");

    return PR;
}

private:

template <typename value_t, typename index_t>
static sparse::matrix<value_t, index_t>
smoothed_interpolation(
        const sparse::matrix<value_t, index_t> &A,
        const std::vector<char> &S,
        const std::vector<index_t> &aggr,
        index_t nc, std::vector<value_t> &omega)
{
    const index_t n = sparse::matrix_rows(A);

    sparse::matrix<value_t, index_t> AP(n, nc);
    std::fill(AP.row.begin(), AP.row.end(), static_cast<index_t>(0));

    std::vector<value_t> omega_p(nc, static_cast<value_t>(0));
    std::vector<value_t> denum(nc, static_cast<value_t>(0));

    std::vector<value_t> D(n);

#pragma omp parallel
    {
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

        std::vector<index_t> marker(nc, static_cast<index_t>(-1));

        // Diagonal of filtered matrix.
        for(index_t i = chunk_start; i < chunk_end; ++i) {
            value_t dia = 0;
            for(index_t j = A.row[i], e = A.row[i + 1]; j < e; ++j) {
                index_t c = A.col[j];
                value_t v = A.val[j];

                if (c == i)
                    dia += v;
                else if (!S[j])
                    dia -= v;
            }

            D[i] = dia;
        }

#pragma omp barrier

        // Compute A * P_tent product. P_tent is stored implicitly in aggr.
        // 1. Structure of the product result:
        for(index_t i = chunk_start; i < chunk_end; ++i) {
            for(index_t j = A.row[i], e = A.row[i + 1]; j < e; ++j) {
                index_t c = A.col[j];

                if (c != i && !S[j]) continue;
                index_t g = aggr[c]; if (g < 0) continue;

                if (marker[g] != i) {
                    marker[g] = i;
                    ++AP.row[i + 1];
                }
            }
        }

        std::fill(marker.begin(), marker.end(), static_cast<index_t>(-1));

#pragma omp barrier
#pragma omp single
        {
            std::partial_sum(AP.row.begin(), AP.row.end(), AP.row.begin());
            AP.reserve(AP.row.back());
        }

        // 2. Compute the product result.
        for(index_t i = chunk_start; i < chunk_end; ++i) {
            index_t row_beg = AP.row[i];
            index_t row_end = row_beg;
            for(index_t j = A.row[i], e = A.row[i + 1]; j < e; ++j) {
                index_t c = A.col[j];

                if (c != i && !S[j]) continue;
                index_t g = aggr[c]; if (g < 0) continue;

                value_t v = (c == i ? D[i] : A.val[j]);

                if (marker[g] < row_beg) {
                    marker[g] = row_end;
                    AP.col[row_end] = g;
                    AP.val[row_end] = v;
                    ++row_end;
                } else {
                    AP.val[marker[g]] += v;
                }
            }

            // Sort columns in the new row.
            sparse::insertion_sort(&AP.col[row_beg], &AP.val[row_beg], row_end - row_beg);
        }

        std::fill(marker.begin(), marker.end(), static_cast<index_t>(-1));
        std::vector< std::pair<index_t, value_t> > adap(128);

#pragma omp barrier

        // Compute A * Dinv * AP row by row and compute columnwise scalar products
        // necessary for computation of omega_p. The actual results of
        // matrix-matrix product are not stored.
        for(index_t ia = chunk_start; ia < chunk_end; ++ia) {
            adap.clear();

            // Form current row of ADAP matrix.
            for(index_t ja = A.row[ia], ea = A.row[ia + 1]; ja < ea; ++ja) {
                index_t ca = A.col[ja];

                if (ca != ia && !S[ja]) continue;

                value_t dia = D[ca];
                value_t va = (ca == ia ? dia : A.val[ja]);

                for(index_t jb = AP.row[ca], eb = AP.row[ca + 1]; jb < eb; ++jb) {
                    index_t cb = AP.col[jb];
                    value_t vb = AP.val[jb] / dia;

                    if (marker[cb] < 0) {
                        marker[cb] = adap.size();
                        adap.push_back(std::make_pair(cb, va * vb));
                    } else {
                        adap[marker[cb]].second += va * vb;
                    }
                }
            }

            std::sort(adap.begin(), adap.end());

            // Update columnwise scalar products (AP,ADAP) and (ADAP,ADAP).
            // 1. (AP, ADAP)
            for(
                    index_t ja = AP.row[ia], ea = AP.row[ia + 1],
                    jb = 0, eb = adap.size();
                    ja < ea && jb < eb;
               )
            {
                index_t ca = AP.col[ja];
                index_t cb = adap[jb].first;

                if (ca < cb)
                    ++ja;
                else if (cb < ca)
                    ++jb;
                else /*ca == cb*/ {
#pragma omp atomic
                    omega_p[ca] += AP.val[ja] * adap[jb].second;
                    ++ja;
                    ++jb;
                }
            }

            // 2. (ADAP, ADAP) (and clear marker)
            for(index_t j = 0, e = adap.size(); j < e; ++j) {
                index_t c = adap[j].first;
                value_t v = adap[j].second;
#pragma omp atomic
                denum[c] += v * v;
                marker[c] = -1;
            }
        }
    }

    std::transform(omega_p.begin(), omega_p.end(), denum.begin(), omega_p.begin(),
            std::divides<value_t>());

    // Convert omega from (4.13) to (4.14) (Sala, Tuminaro, 2008):
#pragma omp parallel for schedule(dynamic, 1024)
    for(index_t i = 0; i < n; ++i) {
        value_t w = -1;
        for(index_t j = A.row[i], e = A.row[i + 1]; j < e; ++j) {
            index_t c = A.col[j];
            if (c != i && !S[j]) continue;
            index_t g = aggr[c]; if (g < 0) continue;
            if (omega_p[g] < w || w < 0) w = omega_p[g];
        }
        omega[i] = std::max(w, static_cast<value_t>(0));
    }

    // Update AP to obtain P.
#pragma omp parallel for schedule(dynamic, 1024)
    for(index_t i = 0; i < n; ++i) {
        value_t wd = omega[i] / D[i];
        for(index_t j = AP.row[i], e = AP.row[i + 1]; j < e; ++j)
            AP.val[j] = (AP.col[j] == aggr[i] ? 1 : 0) - wd * AP.val[j];
    }

    return AP;
}

template <typename value_t, typename index_t>
static sparse::matrix<value_t, index_t>
smoothed_restriction(
        const sparse::matrix<value_t, index_t> &A,
        const std::vector<char> &S,
        const std::vector<index_t> &aggr,
        index_t nc, const std::vector<value_t> &omega)
{
    const index_t n = sparse::matrix_rows(A);

    // Get structure of R_tent from aggr
    std::vector<index_t> R_tent_row(nc + 1, static_cast<index_t>(0));
    for(index_t i = 0; i < n; ++i) {
        index_t g = aggr[i]; if (g < 0) continue;
        ++R_tent_row[g + 1];
    }

    std::partial_sum(R_tent_row.begin(), R_tent_row.end(), R_tent_row.begin());
    std::vector<index_t> R_tent_col(R_tent_row.back());

    for(index_t i = 0; i < n; ++i) {
        index_t g = aggr[i]; if (g < 0) continue;
        R_tent_col[R_tent_row[g]++] = i;
    }

    std::rotate(R_tent_row.begin(), R_tent_row.end() - 1, R_tent_row.end());
    R_tent_row[0] = 0;

    sparse::matrix<value_t, index_t> R(nc, n);
    std::fill(R.row.begin(), R.row.end(), static_cast<index_t>(0));

    // Diagonal of filtered matrix.
    std::vector<value_t> Dinv(n);
#pragma omp parallel for schedule(dynamic, 1024)
    for(index_t i = 0; i < n; ++i) {
        value_t dia = 0;
        for(index_t j = A.row[i], e = A.row[i + 1]; j < e; ++j) {
            index_t c = A.col[j];
            value_t v = A.val[j];

            if (c == i)
                dia += v;
            else if (!S[j])
                dia -= v;
        }

        Dinv[i] = 1 / dia;
    }

    // Compute R_tent * A * Dinv.
#pragma omp parallel
    {
#ifdef _OPENMP
        int nt  = omp_get_num_threads();
        int tid = omp_get_thread_num();

        index_t chunk_size  = (nc + nt - 1) / nt;
        index_t chunk_start = tid * chunk_size;
        index_t chunk_end   = std::min(nc, chunk_start + chunk_size);
#else
        index_t chunk_start = 0;
        index_t chunk_end   = nc;
#endif

        std::vector<index_t> marker(n, static_cast<index_t>(-1));
        for(index_t ir = chunk_start; ir < chunk_end; ++ir) {
            for(index_t jr = R_tent_row[ir], er = R_tent_row[ir + 1]; jr < er; ++jr) {
                index_t cr = R_tent_col[jr];
                for(index_t ja = A.row[cr], ea = A.row[cr + 1]; ja < ea; ++ja) {
                    index_t ca = A.col[ja];
                    if (ca != cr && !S[ja]) continue;

                    if (marker[ca] != ir) {
                        marker[ca] = ir;
                        ++R.row[ir + 1];
                    }
                }
            }
        }

        std::fill(marker.begin(), marker.end(), static_cast<index_t>(-1));

#pragma omp barrier
#pragma omp single
        {
            std::partial_sum(R.row.begin(), R.row.end(), R.row.begin());
            R.reserve(R.row.back());
        }

        for(index_t ir = chunk_start; ir < chunk_end; ++ir) {
            index_t row_beg = R.row[ir];
            index_t row_end = row_beg;

            for(index_t jr = R_tent_row[ir], er = R_tent_row[ir + 1]; jr < er; ++jr) {
                index_t cr = R_tent_col[jr];

                for(index_t ja = A.row[cr], ea = A.row[cr + 1]; ja < ea; ++ja) {
                    index_t ca = A.col[ja];
                    if (ca != cr && !S[ja]) continue;
                    value_t va = (ca == cr ? 1 : A.val[ja] * Dinv[ca]);

                    if (marker[ca] < row_beg) {
                        marker[ca] = row_end;
                        R.col[row_end] = ca;
                        R.val[row_end] = va;
                        ++row_end;
                    } else {
                        R.val[marker[ca]] += va;
                    }
                }
            }
        }
    }

    // Update R.
#pragma omp parallel for schedule(dynamic, 1024)
    for(index_t i = 0; i < nc; ++i) {
        for(index_t j = R.row[i], e = R.row[i + 1]; j < e; ++j) {
            index_t c = R.col[j];
            R.val[j] = (aggr[c] == i ? 1 : 0) - omega[c] * R.val[j];
        }
    }

    return R;
}

};

} // namespace interp
} // namespace amgcl



#endif
