#ifndef AMGCL_AGGR_PLAIN_HPP
#define AMGCL_AGGR_PLAIN_HPP

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
 * \file   aggr_plain.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  Plain aggregation.
 */

#include <vector>
#include <algorithm>

#include <amgcl/spmat.hpp>
#include <amgcl/tictoc.hpp>

namespace amgcl {

/// Aggregation related types and functions.
namespace aggr {

/**
 * \defgroup aggregation Aggregation
 * \brief Possible aggregation strategies.
 */

/// Plain aggregation.
/**
 * Modification of a greedy aggregation scheme from \ref Vanek_1996 "Vanek (1996)".
 * Any nonzero matrix entry forms a connection. Variables without neighbours
 * (resulting, e.g., from Dirichlet conditions) are excluded from aggregation
 * process. The aggregation is completed in a single pass over variables:
 * variables adjacent to a new aggregate are temporarily marked as beloning to
 * this aggregate. Later they may be claimed by other aggregates; if nobody
 * claims them, then they just stay in their initial aggregate.
 *
 * \ingroup aggregation
 */
struct plain {

/// Constructs aggregates of variables.
/**
 * Each entry of the return vector corresponds to a variable and contains
 * number of an aggregate the variable belongs to. If an entry is negative,
 * then variable does not belong to any aggregate.
 *
 * \param A The system matrix.
 * \param S Strong couplings in A.
 *
 * \returns a vector of aggregate numbers.
 */
template <class spmat>
static std::vector< typename sparse::matrix_index<spmat>::type >
aggregates( const spmat &A, const std::vector<char> &S ) {
    typedef typename sparse::matrix_index<spmat>::type index_t;
    typedef typename sparse::matrix_value<spmat>::type value_t;

    const index_t n = sparse::matrix_rows(A);

    const index_t undefined = static_cast<index_t>(-1);
    const index_t removed   = static_cast<index_t>(-2);

    std::vector<index_t> agg(n);

    const index_t *Arow = sparse::matrix_outer_index(A);
    const index_t *Acol = sparse::matrix_inner_index(A);

    // Remove nodes without neighbours
    index_t max_neib = 0;
    for(index_t i = 0; i < n; ++i) {
        index_t j = Arow[i];
        index_t e = Arow[i + 1];

        max_neib = std::max(e - j, max_neib);

        index_t state = removed;
        for(; j < e; ++j)
            if (S[j]) {
                state = undefined;
                break;
            }

        agg[i] = state;
    }

    std::vector<index_t> neib;
    neib.reserve(max_neib);

    index_t last_g = static_cast<index_t>(-1);

    // Perform plain aggregation
    for(index_t i = 0; i < n; ++i) {
        if (agg[i] != undefined) continue;

        // The point is not adjacent to a core of any previous aggregate:
        // so its a seed of a new aggregate.
        agg[i] = ++last_g;

        neib.clear();

        // Include its neighbors as well.
        for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j) {
            index_t c = Acol[j];
            if (S[j] && agg[c] != removed) {
                agg[c] = last_g;
                neib.push_back(c);
            }
        }

        // Temporarily mark undefined points adjacent to the new aggregate as
        // beloning to the aggregate. If nobody claims them later, they will
        // stay here.
        for(typename std::vector<index_t>::const_iterator nb = neib.begin(); nb != neib.end(); ++nb)
            for(index_t j = Arow[*nb], e = Arow[*nb + 1]; j < e; ++j)
                if (S[j] && agg[Acol[j]] == undefined) agg[Acol[j]] = last_g;
    }

    assert( std::count(agg.begin(), agg.end(), undefined) == 0 );

    return agg;
}

};

} // namespace aggr
} // namespace amgcl

#endif
