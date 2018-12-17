#ifndef AMGCL_COARSENING_PLAIN_AGGREGATES_HPP
#define AMGCL_COARSENING_PLAIN_AGGREGATES_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/coarsening/plain_aggregates.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Plain aggregation.
 */

#include <vector>
#include <numeric>

#include <amgcl/util.hpp>
#include <amgcl/backend/builtin.hpp>

namespace amgcl {
namespace coarsening {

/**
 * \defgroup aggregates Aggregates
 * \brief These classes control how fine-level variables are subdivided into
 * aggregates.
 */

/// Plain aggregation.
/**
 * Modification of a greedy aggregation scheme from \cite Vanek1996.
 * Connectivity is defined in a symmetric way, that is, two variables \f$i\f$
 * and \f$j\f$ are considered to be connected to each other if
 * \f$a_{ij}^2/a_{ii}a_{jj} > \varepsilon_{strong}\f$. Variables without
 * neighbours (resulting, e.g., from Dirichlet conditions) are excluded from
 * aggregation process. The aggregation is completed in a single pass over
 * variables: variables adjacent to a new aggregate are temporarily marked as
 * beloning to this aggregate. Later they may be claimed by other aggregates;
 * if nobody claims them, then they just stay in their initial aggregate.
 *
 * \ingroup aggregates
 */
struct plain_aggregates {
    /// Aggregation parameters.
    struct params {
        /// Parameter \f$\varepsilon_{strong}\f$ defining strong couplings.
        /**
         * Connectivity is defined in a symmetric way, that is, two variables
         * \f$i\f$ and \f$j\f$ are considered to be connected to each other if
         * \f$a_{ij}^2/a_{ii}a_{jj} > \varepsilon_{strong}\f$ with fixed \f$0 <
         * \varepsilon_{strong} < 1.\f$
         */
        float eps_strong;

        params() : eps_strong(0.08f) {}

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, eps_strong)
        {
            check_params(p, {"eps_strong", "block_size"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, eps_strong);
        }
#endif
    };

    static const ptrdiff_t undefined = -1;
    static const ptrdiff_t removed   = -2;

    /// Number of aggregates.
    size_t count;

    /// Strong connectivity matrix.
    /**
     * This is just 'values' part of CRS matrix. 'col' and 'ptr' arrays are
     * borrowed from the system matrix.
     */
    std::vector<char> strong_connection;

    /// Aggerate id that each fine-level variable belongs to.
    /** When id[i] < 0, then variable i stays at the fine level (this could be
     * the case for a Dirichelt condition variable).*/
    std::vector<ptrdiff_t> id;

    /// Constructs aggregates for a given matrix.
    /**
     * \param A   The system matrix.
     * \param prm Aggregation parameters.
     */
    template <class Matrix>
    plain_aggregates(const Matrix &A, const params &prm)
        : count(0),
          strong_connection( backend::nonzeros(A) ),
          id( backend::rows(A) )
    {
        typedef typename backend::value_type<Matrix>::type value_type;
        typedef typename math::scalar_of<value_type>::type scalar_type;

        scalar_type eps_squared = prm.eps_strong * prm.eps_strong;

        const size_t n = rows(A);

        /* 1. Get strong connections */
        auto dia = diagonal(A);
#pragma omp parallel for
        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            value_type eps_dia_i = eps_squared * (*dia)[i];

            for(ptrdiff_t j = A.ptr[i], e = A.ptr[i+1]; j < e; ++j) {
                ptrdiff_t c = A.col[j];
                value_type v = A.val[j];

                strong_connection[j] = (c != i) && (eps_dia_i * (*dia)[c] < v * v);
            }
        }

        /* 2. Get aggregate ids */

        // Remove lonely nodes.
        size_t max_neib = 0;
        for(size_t i = 0; i < n; ++i) {
            ptrdiff_t j = A.ptr[i], e = A.ptr[i+1];
            max_neib    = std::max<size_t>(max_neib, e - j);

            ptrdiff_t state = removed;
            for(; j < e; ++j)
                if (strong_connection[j]) {
                    state = undefined;
                    break;
                }

            id[i] = state;
        }

        std::vector<ptrdiff_t> neib;
        neib.reserve(max_neib);

        // Perform plain aggregation
        for(size_t i = 0; i < n; ++i) {
            if (id[i] != undefined) continue;

            // The point is not adjacent to a core of any previous aggregate:
            // so its a seed of a new aggregate.
            ptrdiff_t cur_id = static_cast<ptrdiff_t>(count++);
            id[i] = cur_id;

            // (*) Include its neighbors as well.
            neib.clear();
            for(ptrdiff_t j = A.ptr[i], e = A.ptr[i+1]; j < e; ++j) {
                ptrdiff_t c = A.col[j];
                if (strong_connection[j] && id[c] != removed) {
                    id[c] = cur_id;
                    neib.push_back(c);
                }
            }

            // Temporarily mark undefined points adjacent to the new aggregate
            // as members of the aggregate.
            // If nobody claims them later, they will stay here.
            for(ptrdiff_t c : neib) {
                for(ptrdiff_t j = A.ptr[c], e = A.ptr[c+1]; j < e; ++j) {
                    ptrdiff_t cc = A.col[j];
                    if (strong_connection[j] && id[cc] == undefined)
                        id[cc] = cur_id;
                }
            }
        }

        precondition(count > 0, "Zero aggregates found.");

        // Some of the aggregates could potentially vanish during expansion
        // step (*) above. We need to exclude those and renumber the rest.
        std::vector<ptrdiff_t> cnt(count, 0);
        for(ptrdiff_t i : id)
            if (i >= 0) cnt[i] = 1;
        std::partial_sum(cnt.begin(), cnt.end(), cnt.begin());

        if (static_cast<ptrdiff_t>(count) > cnt.back()) {
            count = cnt.back();

            for(size_t i = 0; i < n; ++i)
                if (id[i] >= 0) id[i] = cnt[id[i]] - 1;
        }
    }
};

} // namespace coarsening
} // namespace amgcl

#endif
