#ifndef AMGCL_COARSENING_SMOOTHED_AGGR_EMIN_HPP
#define AMGCL_COARSENING_SMOOTHED_AGGR_EMIN_HPP

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
 * \file   amgcl/coarsening/smoothed_aggr_emin.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Smoothed aggregation with energy minimization coarsening.
 */

#include <limits>

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
#include <amgcl/detail/sort_row.hpp>

namespace amgcl {
namespace coarsening {
namespace detail {

template <class Base>
struct sa_emin_filtered_matrix {
    typedef typename backend::value_type<Base>::type value_type;

    typedef value_type val_type;
    typedef ptrdiff_t  col_type;
    typedef ptrdiff_t  ptr_type;

    const Base &base;
    const std::vector<char> &strong;
    std::vector<val_type> dia;

    class row_iterator {
        public:
            row_iterator(
                    const col_type * col,
                    const col_type * end,
                    const val_type * val,
                    const char     * str,
                    col_type row,
                    val_type dia
                    )
                : m_col(col), m_end(end), m_val(val), m_str(str),
                  m_row(row), m_dia(dia)
            {}

            operator bool() const {
                return m_col < m_end;
            }

            row_iterator& operator++() {
                do {
                    ++m_col;
                    ++m_val;
                    ++m_str;
                } while(m_col < m_end && m_col[0] != m_row && !m_str[0]);

                return *this;
            }

            col_type col() const {
                return *m_col;
            }

            val_type value() const {
                if (m_col[0] == m_row)
                    return m_dia;
                else
                    return m_val[0];
            }

        private:
            const col_type * m_col;
            const col_type * m_end;
            const val_type * m_val;
            const char     * m_str;

            col_type m_row;
            val_type m_dia;
    };

    sa_emin_filtered_matrix(
            const Base &base,
            const std::vector<char> &strong
            ) : base(base), strong(strong), dia( backend::rows(base) )
    {
        const size_t n = backend::rows(base);
#pragma omp parallel for
        for(ptrdiff_t i = 0; i < ptrdiff_t(n); ++i) {
            value_type D = 0;
            for(ptrdiff_t j = base.ptr[i], e = base.ptr[i + 1]; j < e; ++j) {
                ptrdiff_t  c = base.col[j];
                value_type v = base.val[j];

                if (c == i)
                    D += v;
                else if (!strong[j])
                    D -= v;
            }

            dia[i] = D;
        }
    }

    size_t rows() const {
        return backend::rows(base);
    }

    size_t cols() const {
        return backend::cols(base);
    }

    row_iterator row_begin(size_t row) const {
        ptr_type b = base.ptr[row];
        ptr_type e = base.ptr[row + 1];

        const col_type *col = &base.col[b];
        const col_type *end = &base.col[e];
        const val_type *val = &base.val[b];
        const char     *str = &strong[b];

        return row_iterator(col, end, val, str, row, dia[row]);
    }
};

} // namespace detail
} // namespace coarsening

namespace backend {

template <class Base>
struct rows_impl< coarsening::detail::sa_emin_filtered_matrix<Base> > {
    typedef coarsening::detail::sa_emin_filtered_matrix<Base> Matrix;

    static size_t get(const Matrix &A) {
        return A.rows();
    }
};

template <class Base>
struct cols_impl< coarsening::detail::sa_emin_filtered_matrix<Base> > {
    typedef coarsening::detail::sa_emin_filtered_matrix<Base> Matrix;

    static size_t get(const Matrix &A) {
        return A.cols();
    }
};

template <class Base>
struct row_iterator< coarsening::detail::sa_emin_filtered_matrix<Base> > {
    typedef
        typename coarsening::detail::sa_emin_filtered_matrix<Base>::row_iterator
        type;
};

template <class Base>
struct row_begin_impl< coarsening::detail::sa_emin_filtered_matrix<Base> > {
    typedef coarsening::detail::sa_emin_filtered_matrix<Base> Matrix;

    static typename row_iterator<Matrix>::type
    get(const Matrix &matrix, size_t row) {
        return matrix.row_begin(row);
    }
};

} // namespace backend

namespace coarsening {

/// Smoothed aggregation with energy minimization.
/**
 * \ingroup coarsening
 * \sa \cite Sala2008
 */
struct smoothed_aggr_emin {
    typedef pointwise_aggregates Aggregates;

    /// Coarsening parameters.
    struct params {
        /// Aggregation parameters.
        Aggregates::params aggr;

        /// Near nullspace parameters.
        nullspace_params nullspace;

        params() {}

        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_CHILD(p, aggr),
              AMGCL_PARAMS_IMPORT_CHILD(p, nullspace)
        {}

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_CHILD(p, path, aggr);
            AMGCL_PARAMS_EXPORT_CHILD(p, path, nullspace);
        }
    };

    /// \copydoc amgcl::coarsening::aggregation::transfer_operators
    template <class Matrix>
    static boost::tuple<
        boost::shared_ptr<Matrix>,
        boost::shared_ptr<Matrix>
        >
    transfer_operators(const Matrix &A, params &prm)
    {
        typedef typename backend::value_type<Matrix>::type Val;

        TIC("aggregates");
        Aggregates aggr(A, prm.aggr);
        prm.aggr.eps_strong *= 0.5;
        TOC("aggregates");

        TIC("interpolation");
        boost::shared_ptr<Matrix> P_tent = tentative_prolongation<Matrix>(
                rows(A), aggr.count, aggr.id, prm.nullspace, prm.aggr.block_size
                );

        std::vector<Val> omega;
        detail::sa_emin_filtered_matrix<Matrix> Af(A, aggr.strong_connection);

        boost::shared_ptr<Matrix> P = interpolation(Af, *P_tent, omega);
        boost::shared_ptr<Matrix> R = restriction  (Af, *P_tent, omega);
        TOC("interpolation");

        if (prm.nullspace.cols > 0)
            prm.aggr.block_size = prm.nullspace.cols;

        return boost::make_tuple(P, R);
    }

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

    private:
        template <class AMatrix, typename Val, typename Col, typename Ptr>
        static boost::shared_ptr< backend::crs<Val, Col, Ptr> >
        interpolation(
                const AMatrix &A,
                const backend::crs<Val, Col, Ptr> &P_tent,
                std::vector<Val> &omega
                )
        {
            typedef backend::crs<Val, Col, Ptr> PMatrix;

            typedef typename PMatrix::row_iterator Piterator;
            typedef typename AMatrix::row_iterator Aiterator;

            const size_t n  = rows(P_tent);
            const size_t nc = cols(P_tent);

            boost::shared_ptr<PMatrix> AP = boost::make_shared<PMatrix>();

            *AP = product(A, P_tent, /*sort rows: */true);

            omega.resize(nc, Val(0));
            std::vector<Val> denum(nc, 0);

#pragma omp parallel
            {
#ifdef _OPENMP
                int nt  = omp_get_num_threads();
                int tid = omp_get_thread_num();

                size_t chunk_size  = (n + nt - 1) / nt;
                size_t chunk_start = tid * chunk_size;
                size_t chunk_end   = std::min(n, chunk_start + chunk_size);
#else
                size_t chunk_start = 0;
                size_t chunk_end   = n;
#endif

                std::vector<ptrdiff_t> marker(nc, -1);

                // Compute A * Dinv * AP row by row and compute columnwise
                // scalar products necessary for computation of omega. The
                // actual results of matrix-matrix product are not stored.
                std::vector<Col> adap_col(128);
                std::vector<Val> adap_val(128);

                for(size_t ia = chunk_start; ia < chunk_end; ++ia) {
                    adap_col.clear();
                    adap_val.clear();

                    // Form current row of ADAP matrix.
                    for(Aiterator a = A.row_begin(ia); a; ++a) {
                        Col ca  = a.col();
                        Val va  = a.value() / A.dia[ca];

                        for(Piterator p = AP->row_begin(ca); p; ++p) {
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
                            adap_col.data(), adap_val.data(), adap_col.size()
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
#pragma omp atomic
                            omega[ca] += AP->val[ja] * adap_val[jb];
                            ++ja;
                            ++jb;
                        }
                    }

                    // 2. (ADAP, ADAP) (and clear marker)
                    for(size_t j = 0, e = adap_col.size(); j < e; ++j) {
                        Col c = adap_col[j];
                        Val v = adap_val[j];
#pragma omp atomic
                        denum[c] += v * v;
                        marker[c] = -1;
                    }
                }
            }

            boost::transform(omega, denum, omega.begin(), std::divides<Val>());

            // Update AP to obtain P: P = (P_tent - D^-1 A P Omega)
            /*
             * Here we use the fact that if P(i,j) != 0,
             * then with necessity AP(i,j) != 0:
             *
             * AP(i,j) = sum_k(A_ik P_kj), and A_ii != 0.
             */
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                Val dia = 1 / A.dia[i];

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
        static boost::shared_ptr< backend::crs<Val, Col, Ptr> >
        restriction(
                const AMatrix &A,
                const backend::crs<Val, Col, Ptr> &P_tent,
                const std::vector<Val> &omega
                )
        {
            typedef backend::crs<Val, Col, Ptr> PMatrix;

            const size_t nc = cols(P_tent);

            PMatrix R_tent = transpose(P_tent);
            sort_rows(R_tent);

            boost::shared_ptr<PMatrix> RA = boost::make_shared<PMatrix>();
            *RA = product(R_tent, A, /*sort rows: */true);

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

                for(Ptr ja = RA->ptr[i],    ea = RA->ptr[i+1],
                        jr = R_tent.ptr[i], er = R_tent.ptr[i+1];
                        ja < ea; ++ja
                   )
                {
                    Col ca = RA->col[ja];
                    Val va = -w * RA->val[ja] / A.dia[ca];

                    for(; jr < er; ++jr) {
                        Col cr = R_tent.col[jr];
                        if (cr > ca)
                            break;

                        if (cr == ca) {
                            va += R_tent.val[jr];
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
