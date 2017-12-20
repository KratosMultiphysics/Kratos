#ifndef AMGCL_COARSENING_POINTWISE_AGGREGATES_HPP
#define AMGCL_COARSENING_POINTWISE_AGGREGATES_HPP

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
 * \file   amgcl/coarsening/pointwise_aggregates.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Pointwise aggregation.
 */

#include <vector>
#include <cmath>

#include <amgcl/util.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/plain_aggregates.hpp>

namespace amgcl {
namespace coarsening {

/// Pointwise aggregation.
/**
 * The system matrix should have block structure. It is reduced to a single
 * value per block and is subjected to coarsening::plain_aggregation.
 *
 * \ingroup aggregates
 */
class pointwise_aggregates {
    public:
        /// Aggregation parameters.
        struct params : plain_aggregates::params {
            /// Block size for the system matrix.
            /**
             * When block_size=1, the scheme is equivalent to (and performs on
             * par with) plain_aggregates.
             */
            unsigned block_size;

            params() : block_size(1) {}

            params(const boost::property_tree::ptree &p)
                : plain_aggregates::params(p),
                  AMGCL_PARAMS_IMPORT_VALUE(p, block_size)
            {
                AMGCL_PARAMS_CHECK(p, (eps_strong)(block_size));
            }

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                plain_aggregates::params::get(p, path);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, block_size);
            }
        };

        static const ptrdiff_t undefined = -1;
        static const ptrdiff_t removed   = -2;

        /// \copydoc amgcl::coarsening::plain_aggregates::count
        size_t count;

        /// \copydoc amgcl::coarsening::plain_aggregates::strong_connection
        std::vector<char> strong_connection;

        /// \copydoc amgcl::coarsening::plain_aggregates::id
        std::vector<ptrdiff_t> id;

        /// \copydoc amgcl::coarsening::plain_aggregates::plain_aggregates
        template <class Matrix>
        pointwise_aggregates(const Matrix &A, const params &prm, unsigned min_aggregate)
            : count(0)
        {
            typedef typename backend::value_type<Matrix>::type value_type;
            typedef typename math::scalar_of<value_type>::type scalar_type;
            if (prm.block_size == 1) {
                plain_aggregates aggr(A, prm);

                remove_small_aggregates(A.nrows, 1, min_aggregate, aggr);

                count = aggr.count;
                strong_connection.swap(aggr.strong_connection);
                id.swap(aggr.id);
            } else {
                strong_connection.resize( nonzeros(A) );
                id.resize( rows(A) );

                backend::crs<scalar_type> Ap = pointwise_matrix(A, prm.block_size);

                plain_aggregates pw_aggr(Ap, prm);

                remove_small_aggregates(
                        Ap.nrows, prm.block_size, min_aggregate, pw_aggr);

                count = pw_aggr.count * prm.block_size;

#pragma omp parallel
                {
                    std::vector<ptrdiff_t> marker(Ap.nrows, -1);

#pragma omp for
                    for(ptrdiff_t ip = 0; ip < static_cast<ptrdiff_t>(Ap.nrows); ++ip) {
                        ptrdiff_t ia = ip * prm.block_size;
                        ptrdiff_t row_beg = Ap.ptr[ip];
                        ptrdiff_t row_end = row_beg;

                        for(unsigned k = 0; k < prm.block_size; ++k, ++ia) {
                            id[ia] = prm.block_size * pw_aggr.id[ip] + k;

                            for(ptrdiff_t ja = A.ptr[ia], ea = A.ptr[ia+1]; ja < ea; ++ja) {
                                ptrdiff_t cp = A.col[ja] / prm.block_size;

                                if (marker[cp] < row_beg) {
                                    marker[cp] = row_end;
                                    strong_connection[ja] = pw_aggr.strong_connection[row_end];
                                    ++row_end;
                                } else {
                                    strong_connection[ja] = pw_aggr.strong_connection[ marker[cp] ];
                                }
                            }
                        }
                    }
                }
            }
        }

        static void remove_small_aggregates(
                size_t n, unsigned block_size, unsigned min_aggregate,
                plain_aggregates &aggr
                )
        {
            if (min_aggregate <= 1) return; // nothing to do

            // Count entries in each of the aggregates
            std::vector<ptrdiff_t> count(aggr.count, 0);

            for(size_t i = 0; i < n; ++i) {
                ptrdiff_t id = aggr.id[i];
                if (id != removed) ++count[id];
            }

            // If any aggregate has less entries than required, remove it.
            // Renumber the rest of the aggregates to leave no gaps.
            size_t m = 0;
            for(size_t i = 0; i < aggr.count; ++i) {
                if (block_size * count[i] < min_aggregate) {
                    count[i] = removed;
                } else {
                    count[i] = m++;
                }
            }

            // Update aggregate count and aggregate ids.
            aggr.count = m;

            for(size_t i = 0; i < n; ++i) {
                ptrdiff_t id = aggr.id[i];
                if (id != removed) aggr.id[i] = count[id];
            }
        }

        template <class Matrix>
        static backend::crs<
            typename math::scalar_of<
                typename backend::value_type<Matrix>::type
                >::type
            >
        pointwise_matrix(const Matrix &A, size_t block_size) {
            typedef typename backend::value_type<Matrix>::type V;
            typedef typename math::scalar_of<V>::type S;
            typedef typename backend::row_iterator<Matrix>::type row_iterator;

            const size_t n  = backend::rows(A);
            const size_t m  = backend::cols(A);
            const size_t np = n / block_size;
            const size_t mp = m / block_size;

            precondition(n % block_size == 0 && m % block_size == 0,
                    "Matrix size should be divisible by block_size");

            backend::crs<S> Ap;
            Ap.set_size(np, mp, true);

#pragma omp parallel
            {
                std::vector<ptrdiff_t> marker(mp, -1);

                // Count number of nonzeros in block matrix.
#pragma omp for
                for(ptrdiff_t ip = 0; ip < static_cast<ptrdiff_t>(np); ++ip) {
                    ptrdiff_t ia = ip * block_size;

                    for(unsigned k = 0; k < block_size; ++k, ++ia) {
                        for(row_iterator a = backend::row_begin(A, ia); a; ++a) {
                            ptrdiff_t cp = a.col() / block_size;
                            if (marker[cp] != ip) {
                                marker[cp] = ip;
                                ++Ap.ptr[ip + 1];
                            }
                        }
                    }
                }
            }

            std::partial_sum(Ap.ptr, Ap.ptr + np + 1, Ap.ptr);
            Ap.set_nonzeros();

#pragma omp parallel
            {
                std::vector<ptrdiff_t> marker(mp, -1);

                // Fill the reduced matrix. Use max norm for blocks.
#pragma omp for
                for(ptrdiff_t ip = 0; ip < static_cast<ptrdiff_t>(np); ++ip) {
                    ptrdiff_t ia = ip * block_size;
                    ptrdiff_t row_beg = Ap.ptr[ip];
                    ptrdiff_t row_end = row_beg;

                    for(unsigned k = 0; k < block_size; ++k, ++ia) {
                        for(row_iterator a = backend::row_begin(A, ia); a; ++a) {
                            ptrdiff_t cb = a.col() / block_size;
                            S va = math::norm(a.value());

                            if (marker[cb] < row_beg) {
                                marker[cb] = row_end;
                                Ap.col[row_end] = cb;
                                Ap.val[row_end] = va;
                                ++row_end;
                            } else {
                                Ap.val[marker[cb]] = std::max(Ap.val[marker[cb]], va);
                            }
                        }
                    }
                }
            }

            return Ap;
        }
};

} // namespace coarsening
} // namespace amgcl


#endif
