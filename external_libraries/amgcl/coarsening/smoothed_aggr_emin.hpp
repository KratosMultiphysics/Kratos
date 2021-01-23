#ifndef AMGCL_COARSENING_SMOOTHED_AGGR_EMIN_HPP
#define AMGCL_COARSENING_SMOOTHED_AGGR_EMIN_HPP

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
 * \file   amgcl/coarsening/smoothed_aggr_emin.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Smoothed aggregation with energy minimization coarsening.
 */

#include <limits>

#include <tuple>
#include <memory>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/detail/galerkin.hpp>
#include <amgcl/coarsening/pointwise_aggregates.hpp>
#include <amgcl/coarsening/tentative_prolongation.hpp>
#include <amgcl/util.hpp>
#include <amgcl/detail/sort_row.hpp>

namespace amgcl {
namespace coarsening {

/// Smoothed aggregation with energy minimization.
/**
 * \ingroup coarsening
 * \sa \cite Sala2008
 */
template <class Backend>
struct smoothed_aggr_emin {
    typedef pointwise_aggregates Aggregates;

    /// Coarsening parameters.
    struct params {
        /// Aggregation parameters.
        Aggregates::params aggr;

        /// Near nullspace parameters.
        nullspace_params nullspace;

        params() {}

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_CHILD(p, aggr),
              AMGCL_PARAMS_IMPORT_CHILD(p, nullspace)
        {
            check_params(p, {"aggr", "nullspace"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_CHILD(p, path, aggr);
            AMGCL_PARAMS_EXPORT_CHILD(p, path, nullspace);
        }
#endif
    } prm;

    smoothed_aggr_emin(const params &prm = params()) : prm(prm) {}

    /// \copydoc amgcl::coarsening::aggregation::transfer_operators
    template <class Matrix>
    std::tuple<
        std::shared_ptr<Matrix>,
        std::shared_ptr<Matrix>
        >
    transfer_operators(const Matrix &A) {
        typedef typename backend::value_type<Matrix>::type Val;
        typedef ptrdiff_t Idx;

        AMGCL_TIC("aggregates");
        Aggregates aggr(A, prm.aggr, prm.nullspace.cols);
        prm.aggr.eps_strong *= 0.5;
        AMGCL_TOC("aggregates");

        AMGCL_TIC("interpolation");
        auto P_tent = tentative_prolongation<Matrix>(
                rows(A), aggr.count, aggr.id, prm.nullspace, prm.aggr.block_size
                );

        // Filter the system matrix
        backend::crs<Val> Af;
        Af.set_size(rows(A), cols(A));
        Af.ptr[0] = 0;

        std::vector<Val> dia(Af.nrows);

#pragma omp parallel for
        for(Idx i = 0; i < static_cast<Idx>(Af.nrows); ++i) {
            Idx row_begin = A.ptr[i];
            Idx row_end   = A.ptr[i+1];
            Idx row_width = row_end - row_begin;

            Val D = math::zero<Val>();
            for(Idx j = row_begin; j < row_end; ++j) {
                Idx c = A.col[j];
                Val v = A.val[j];

                if (c == i)
                    D += v;
                else if (!aggr.strong_connection[j]) {
                    D += v;
                    --row_width;
                }
            }

            dia[i] = D;
            Af.ptr[i+1] = row_width;
        }

        Af.set_nonzeros(Af.scan_row_sizes());

#pragma omp parallel for
        for(Idx i = 0; i < static_cast<Idx>(Af.nrows); ++i) {
            Idx row_begin = A.ptr[i];
            Idx row_end   = A.ptr[i+1];
            Idx row_head  = Af.ptr[i];

            for(Idx j = row_begin; j < row_end; ++j) {
                Idx c = A.col[j];

                if (c == i) {
                    Af.col[row_head] = i;
                    Af.val[row_head] = dia[i];
                    ++row_head;
                } else if (aggr.strong_connection[j]) {
                    Af.col[row_head] = c;
                    Af.val[row_head] = A.val[j];
                    ++row_head;
                }
            }
        }

        std::vector<Val> omega;

        auto P = interpolation(Af, dia, *P_tent, omega);
        auto R = restriction  (Af, dia, *P_tent, omega);
        AMGCL_TOC("interpolation");

        if (prm.nullspace.cols > 0)
            prm.aggr.block_size = prm.nullspace.cols;

        return std::make_tuple(P, R);
    }

    template <class Matrix>
    std::shared_ptr<Matrix>
    coarse_operator(const Matrix &A, const Matrix &P, const Matrix &R) const {
        return detail::galerkin(A, P, R);
    }

    private:
        template <class AMatrix, typename Val, typename Col, typename Ptr>
        static std::shared_ptr< backend::crs<Val, Col, Ptr> >
        interpolation(
                const AMatrix &A, const std::vector<Val> &Adia,
                const backend::crs<Val, Col, Ptr> &P_tent,
                std::vector<Val> &omega
                )
        {
            const size_t n  = rows(P_tent);
            const size_t nc = cols(P_tent);

            auto AP = product(A, P_tent, /*sort rows: */true);

            omega.resize(nc, math::zero<Val>());
            std::vector<Val> denum(nc, math::zero<Val>());

#pragma omp parallel
            {
                std::vector<ptrdiff_t> marker(nc, -1);

                // Compute A * Dinv * AP row by row and compute columnwise
                // scalar products necessary for computation of omega. The
                // actual results of matrix-matrix product are not stored.
                std::vector<Col> adap_col(128);
                std::vector<Val> adap_val(128);

#pragma omp for
                for(ptrdiff_t ia = 0; ia < static_cast<ptrdiff_t>(n); ++ia) {
                    adap_col.clear();
                    adap_val.clear();

                    // Form current row of ADAP matrix.
                    for(auto a = A.row_begin(ia); a; ++a) {
                        Col ca  = a.col();
                        Val va  = math::inverse(Adia[ca]) * a.value();

                        for(auto p = AP->row_begin(ca); p; ++p) {
                            Col c = p.col();
                            Val v = va * p.value();

                            if (marker[c] < 0) {
                                marker[c] = adap_col.size();
                                adap_col.push_back(c);
                                adap_val.push_back(v);
                            } else {
                                adap_val[marker[c]] += v;
                            }
                        }
                    }

                    amgcl::detail::sort_row(
                            &adap_col[0], &adap_val[0], adap_col.size()
                            );

                    // Update columnwise scalar products (AP,ADAP) and (ADAP,ADAP).
                    // 1. (AP, ADAP)
                    for(
                            Ptr ja = AP->ptr[ia], ea = AP->ptr[ia + 1],
                            jb = 0, eb = adap_col.size();
                            ja < ea && jb < eb;
                       )
                    {
                        Col ca = AP->col[ja];
                        Col cb = adap_col[jb];

                        if (ca < cb)
                            ++ja;
                        else if (cb < ca)
                            ++jb;
                        else /*ca == cb*/ {
                            Val v = AP->val[ja] * adap_val[jb];
#pragma omp critical
                            omega[ca] += v;
                            ++ja;
                            ++jb;
                        }
                    }

                    // 2. (ADAP, ADAP) (and clear marker)
                    for(size_t j = 0, e = adap_col.size(); j < e; ++j) {
                        Col c = adap_col[j];
                        Val v = adap_val[j];
#pragma omp critical
                        denum[c] += v * v;
                        marker[c] = -1;
                    }
                }
            }

            for(size_t i = 0, m = omega.size(); i < m; ++i)
                omega[i] = math::inverse(denum[i]) * omega[i];

            // Update AP to obtain P: P = (P_tent - D^-1 A P Omega)
            /*
             * Here we use the fact that if P(i,j) != 0,
             * then with necessity AP(i,j) != 0:
             *
             * AP(i,j) = sum_k(A_ik P_kj), and A_ii != 0.
             */
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                Val dia = math::inverse(Adia[i]);

                for(Ptr ja = AP->ptr[i],    ea = AP->ptr[i+1],
                        jp = P_tent.ptr[i], ep = P_tent.ptr[i+1];
                        ja < ea; ++ja
                   )
                {
                    Col ca = AP->col[ja];
                    Val va = -dia * AP->val[ja] * omega[ca];

                    for(; jp < ep; ++jp) {
                        Col cp = P_tent.col[jp];
                        if (cp > ca)
                            break;

                        if (cp == ca) {
                            va += P_tent.val[jp];
                            break;
                        }
                    }

                    AP->val[ja] = va;
                }
            }

            return AP;
        }

        template <typename AMatrix, typename Val, typename Col, typename Ptr>
        static std::shared_ptr< backend::crs<Val, Col, Ptr> >
        restriction(
                const AMatrix &A, const std::vector<Val> &Adia,
                const backend::crs<Val, Col, Ptr> &P_tent,
                const std::vector<Val> &omega
                )
        {
            const size_t nc = cols(P_tent);

            auto R_tent = transpose(P_tent);
            sort_rows(*R_tent);

            auto RA = product(*R_tent, A, /*sort rows: */true);

            // Compute R = R_tent - Omega R_tent A D^-1.
            /*
             * Here we use the fact that if R(i,j) != 0,
             * then with necessity RA(i,j) != 0:
             *
             * RA(i,j) = sum_k(R_ik A_kj), and A_jj != 0.
             */
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nc); ++i) {
                Val w = omega[i];

                for(Ptr ja = RA->ptr[i],     ea = RA->ptr[i+1],
                        jr = R_tent->ptr[i], er = R_tent->ptr[i+1];
                        ja < ea; ++ja
                   )
                {
                    Col ca = RA->col[ja];
                    Val va = -w * math::inverse(Adia[ca]) * RA->val[ja];

                    for(; jr < er; ++jr) {
                        Col cr = R_tent->col[jr];
                        if (cr > ca)
                            break;

                        if (cr == ca) {
                            va += R_tent->val[jr];
                            break;
                        }
                    }

                    RA->val[ja] = va;
                }
            }

            return RA;
        }
};

} // namespace coarsening
} // namespace amgcl



#endif
