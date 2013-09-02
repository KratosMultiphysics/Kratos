#ifndef AMGCL_INTERP_AGGR_HPP
#define AMGCL_INTERP_AGGR_HPP

/*
The MIT License

Copyright (c) 2012-2013 Denis Demidov <ddemidov@ksu.ru>

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
 * \file   interp_aggr.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  Aggregates-based interpolation scheme.
 */

#include <vector>
#include <algorithm>

#include <amgcl/spmat.hpp>
#include <amgcl/aggr_connect.hpp>
#include <amgcl/tictoc.hpp>

namespace amgcl {

namespace interp {

/// Aggregation-based interpolation scheme.
/**
 * \param aggr_type \ref aggregation "Aggregation scheme".
 *
 * \ingroup interpolation
 */
template <class aggr_type>
struct aggregation {

/// Parameters controlling aggregation.
struct params {
    /// Over-interpolation factor \f$\alpha\f$.
    /**
     * See \ref Stuben_1999 "Stuben (1999)", Section 9.1 "Re-scaling of the
     * Galerkin operator". [In case of aggregation multigrid] Coarse-grid
     * correction of smooth error, and by this the overall convergence, can
     * often be substantially improved by using "over-interpolation", that is,
     * by multiplying the actual correction (corresponding to piecewise
     * constant interpolation) by some factor \f$\alpha>1\f$. Equivalently,
     * this means that the coarse-level Galerkin operator is re-scaled by
     * \f$1/\alpha\f$:
     * \f[I_h^HA_hI_H^h \to \frac{1}{\alpha}I_h^HA_hI_H^h.\f]
     */
    float over_interp;

    /// Parameter \f$\varepsilon_{str}\f$ defining strong couplings.
    /**
     * Variable \f$i\f$ is defined to be strongly coupled to another variable,
     * \f$j\f$, if \f[|a_{ij}| \geq \varepsilon_{str}\sqrt{a_{ii} a_{jj}}\quad
     * \text{with fixed} \quad 0 < \varepsilon_{str} < 1.\f]
     */
    float eps_strong;

    unsigned dof_per_node;

    params() : over_interp(1.5f), eps_strong(0.1f), dof_per_node(1) {}
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

    std::vector<index_t> aggr;

    assert(prm.dof_per_node > 0);

    if (prm.dof_per_node == 1) {
        // Scalar system. Nothing fancy.
        TIC("aggregates");
        aggr_type::aggregates(A, aggr::connect(A, prm.eps_strong)).swap(aggr);
        TOC("aggregates");
    } else {
        // Non-scalar system.
        // Build reduced matrix, find connections and aggregates with it,
        // restore the vectors to full size.

        std::pair<std::vector<char>, std::vector<index_t> > S_aggr = aggr::pointwise_coarsening<aggr_type>(
                    A, prm.eps_strong, prm.dof_per_node);
        aggr.swap(S_aggr.second);
    }

    index_t nc = std::max(
            static_cast<index_t>(0),
            *std::max_element(aggr.begin(), aggr.end()) + static_cast<index_t>(1)
            );

    TIC("interpolation");
    static std::pair<
        sparse::matrix<value_t, index_t>,
        sparse::matrix<value_t, index_t>
    > PR;

    sparse::matrix<value_t, index_t> &P = PR.first;
    sparse::matrix<value_t, index_t> &R = PR.second;

    P.resize(n, nc);
    P.col.reserve(n);
    P.val.reserve(n);

    P.row[0] = 0;
    for(index_t i = 0; i < n; ++i) {
        if (aggr[i] >= 0) {
            P.row[i + 1] = P.row[i] + 1;
            P.col.push_back(aggr[i]);
            P.val.push_back(static_cast<value_t>(1));
        } else {
            P.row[i + 1] = P.row[i];
        }
    }
    TOC("interpolation");

    sparse::transpose(P).swap(R);
    return PR;
}

};

/// Coarse level computing for aggregation-based AMG.
struct aggregated_operator {
    template <class spmat, class Params>
    static spmat apply(const spmat &R, const spmat &A, const spmat &P,
            const Params &prm)
    {
        typedef typename sparse::matrix_index<spmat>::type index_t;
        typedef typename sparse::matrix_value<spmat>::type value_t;

        // For now this s just a Galerking operator with possible
        // over-interpolation.
        sparse::matrix<value_t, index_t> a = sparse::prod(sparse::prod(R, A), P);

        if (prm.over_interp > 1.0f)
            for(typename std::vector<value_t>::iterator v = a.val.begin(); v != a.val.end(); ++v)
                *v /= prm.over_interp;

        return a;
    }
};

template <class T>
struct coarse_operator< aggregation<T> > {
    typedef aggregated_operator type;
};

} // namespace interp
} // namespace amgcl

#endif
