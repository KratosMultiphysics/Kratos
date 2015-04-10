#ifndef AMGCL_COARSENING_TENTATIVE_PROLONGATION_HPP
#define AMGCL_COARSENING_TENTATIVE_PROLONGATION_HPP

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
 * \file   amgcl/coarsening/tentative_prolongation.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Tentative prolongation operator for aggregated AMG.
 */

#include <vector>
#include <algorithm>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/algorithm.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/detail/qr.hpp>

namespace amgcl {
namespace coarsening {
namespace detail {
    struct skip_negative {
        const std::vector<ptrdiff_t> &key;

        skip_negative(const std::vector<ptrdiff_t> &key) : key(key) { }

        bool operator()(ptrdiff_t i, ptrdiff_t j) const {
            // Cast to unsigned type to keep negative values at the end
            return
                static_cast<size_t>(key[i]) <
                static_cast<size_t>(key[j]);
        }
    };
} // namespace detail


//---------------------------------------------------------------------------
struct nullspace_params {
    /// Number of vectors in near nullspace.
    int cols;

    /// Near nullspace vectors.
    /**
     * The vectors are represented as columns of a 2D matrix stored in
     * row-major order.
     */
    std::vector<double> B;

    nullspace_params() : cols(0) {}

    nullspace_params(const boost::property_tree::ptree &p)
        : cols(p.get("cols", nullspace_params().cols))
    {
        double *b = 0;
        b = p.get("B", b);

        if (b) {
            size_t rows = 0;
            rows = p.get("rows", rows);

            precondition(cols > 0,
                    "Error in nullspace parameters: "
                    "B is set, but cols is not"
                    );

            precondition(rows > 0,
                    "Error in nullspace parameters: "
                    "B is set, but rows is not"
                    );

            B.assign(b, b + rows * cols);
        } else {
            precondition(cols == 0,
                    "Error in nullspace parameters: "
                    "cols > 0, but B is empty"
                    );
        }
    }

    void get(boost::property_tree::ptree&, const std::string&) const {}
};

/// Tentative prolongation operator
/**
 * If near nullspace vectors are not provided, returns piecewise-constant
 * prolongation operator. If user provides near nullspace vectors, those are
 * used to improve the prolongation operator.
 * \see \cite Vanek2001
 */
template <class Matrix>
boost::shared_ptr<Matrix> tentative_prolongation(
        size_t n,
        size_t naggr,
        const std::vector<ptrdiff_t> aggr,
        nullspace_params &nullspace
        )
{
    typedef typename backend::value_type<Matrix>::type value_type;

    boost::shared_ptr<Matrix> P = boost::make_shared<Matrix>();

    TIC("tentative");
    if (nullspace.cols > 0) {
        // Sort fine points by aggregate number.
        // Put points not belonging to any aggregate to the end of the list.
        std::vector<ptrdiff_t> order(
                boost::counting_iterator<ptrdiff_t>(0),
                boost::counting_iterator<ptrdiff_t>(n)
                );
        boost::sort(order, detail::skip_negative(aggr));

        // Precompute the shape of the prolongation operator.
        // Each row contains exactly nullspace.cols non-zero entries.
        // Rows that do not belong to any aggregate are empty.
        P->nrows = n;
        P->ncols = naggr * nullspace.cols;
        P->ptr.reserve(n + 1);

        P->ptr.push_back(0);
        for(size_t i = 0; i < n; ++i) {
            if (aggr[i] < 0)
                P->ptr.push_back(P->ptr.back());
            else
                P->ptr.push_back(P->ptr.back() + nullspace.cols);
        }

        P->col.resize(P->ptr.back());
        P->val.resize(P->ptr.back());

        // Compute the tentative prolongation operator and null-space vectors
        // for the coarser level.
        std::vector<double> Bnew;
        Bnew.reserve(naggr * nullspace.cols * nullspace.cols);

        size_t offset = 0;

        amgcl::detail::QR<double> qr;
        std::vector<double> Bpart;
        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(naggr); ++i) {
            int d = 0;
            Bpart.clear();
            for(size_t j = offset; j < n && aggr[order[j]] == i; ++j, ++d)
                std::copy(
                        &nullspace.B[nullspace.cols * order[j]],
                        &nullspace.B[nullspace.cols * order[j]] + nullspace.cols,
                        std::back_inserter(Bpart)
                        );

            qr.compute(d, nullspace.cols, Bpart.data());

            for(int ii = 0; ii < nullspace.cols; ++ii)
                for(int jj = 0; jj < nullspace.cols; ++jj)
                    Bnew.push_back( qr.R(ii,jj) );

            for(int ii = 0; ii < d; ++ii, ++offset) {
                ptrdiff_t  *c = &P->col[P->ptr[order[offset]]];
                value_type *v = &P->val[P->ptr[order[offset]]];

                for(int jj = 0; jj < nullspace.cols; ++jj) {
                    c[jj] = i * nullspace.cols + jj;
                    v[jj] = qr.Q(ii,jj);
                }
            }
        }

        std::swap(nullspace.B, Bnew);
    } else {
        P->nrows = n;
        P->ncols = naggr;
        P->ptr.reserve(n + 1);
        P->col.reserve(n);

        P->ptr.push_back(0);
        for(size_t i = 0; i < n; ++i) {
            if (aggr[i] >= 0) P->col.push_back(aggr[i]);
            P->ptr.push_back( static_cast<ptrdiff_t>(P->col.size()) );
        }
        P->val.resize(n, static_cast<value_type>(1));
    }
    TOC("tentative");

    return P;
}

} // namespace coarsening
} // namespace amgcl

#endif
