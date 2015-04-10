#ifndef AMGCL_COARSENING_POINTWISE_AGGREGATES_HPP
#define AMGCL_COARSENING_POINTWISE_AGGREGATES_HPP

/*
The MIT License

Copyright (c) 2012-2015 Denis Demidov <dennis.demidov@gmail.com>

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

#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>

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
            {}

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
        pointwise_aggregates(const Matrix &A, const params &prm) : count(0)
        {
            if (prm.block_size == 1) {
                plain_aggregates aggr(A, prm);

                count = aggr.count;
                strong_connection.swap(aggr.strong_connection);
                id.swap(aggr.id);
            } else {
                strong_connection.resize( nonzeros(A) );
                id.resize( rows(A) );

                Matrix Ap = pointwise_matrix(A, prm.block_size);

                plain_aggregates pw_aggr(Ap, prm);
                count = pw_aggr.count * prm.block_size;

#pragma omp parallel
                {
                    std::vector<ptrdiff_t> marker(Ap.nrows, -1);

#ifdef _OPENMP
                    int nt  = omp_get_num_threads();
                    int tid = omp_get_thread_num();

                    size_t chunk_size  = (Ap.nrows + nt - 1) / nt;
                    size_t chunk_start = tid * chunk_size;
                    size_t chunk_end   = std::min(Ap.nrows, chunk_start + chunk_size);
#else
                    size_t chunk_start = 0;
                    size_t chunk_end   = Ap.nrows;
#endif

                    for(size_t ip = chunk_start, ia = ip * prm.block_size; ip < chunk_end; ++ip) {
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

        template <class Matrix>
        static backend::crs<typename backend::value_type<Matrix>::type>
        pointwise_matrix(const Matrix &A, size_t block_size) {
            typedef typename backend::value_type<Matrix>::type   V;
            typedef typename backend::row_iterator<Matrix>::type row_iterator;

            const size_t n  = backend::rows(A);
            const size_t m  = backend::cols(A);
            const size_t np = n / block_size;
            const size_t mp = m / block_size;

            precondition(n % block_size == 0 && m % block_size == 0,
                    "Matrix size should be divisible by block_size");

            backend::crs<V> Ap;
            Ap.nrows = np;
            Ap.ncols = mp;
            Ap.ptr.resize(np + 1, 0);

#pragma omp parallel
            {
                std::vector<ptrdiff_t> marker(mp, -1);

#ifdef _OPENMP
                int nt  = omp_get_num_threads();
                int tid = omp_get_thread_num();

                size_t chunk_size  = (np + nt - 1) / nt;
                size_t chunk_start = tid * chunk_size;
                size_t chunk_end   = std::min(np, chunk_start + chunk_size);
#else
                size_t chunk_start = 0;
                size_t chunk_end   = np;
#endif

                // Count number of nonzeros in block matrix.
                for(size_t ip = chunk_start, ia = ip * block_size; ip < chunk_end; ++ip) {
                    for(unsigned k = 0; k < block_size; ++k, ++ia) {
                        for(row_iterator a = backend::row_begin(A, ia); a; ++a) {
                            ptrdiff_t cp = a.col() / block_size;
                            if (static_cast<size_t>(marker[cp]) != ip) {
                                marker[cp] = ip;
                                ++Ap.ptr[ip + 1];
                            }
                        }
                    }
                }

                boost::fill(marker, -1);

#pragma omp barrier
#pragma omp single
                {
                    boost::partial_sum(Ap.ptr, Ap.ptr.begin());
                    Ap.col.resize(Ap.ptr.back());
                    Ap.val.resize(Ap.ptr.back());
                }

                // Fill the reduced matrix. Use max norm for blocks.
                for(size_t ip = chunk_start, ia = ip * block_size; ip < chunk_end; ++ip) {
                    ptrdiff_t row_beg = Ap.ptr[ip];
                    ptrdiff_t row_end = row_beg;

                    for(unsigned k = 0; k < block_size; ++k, ++ia) {
                        for(row_iterator a = backend::row_begin(A, ia); a; ++a) {
                            ptrdiff_t cb = a.col() / block_size;
                            V    va = fabs(a.value());

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
