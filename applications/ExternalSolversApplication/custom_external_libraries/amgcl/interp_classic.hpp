#ifndef AMGCL_INTERP_CLASSIC_HPP
#define AMGCL_INTERP_CLASSIC_HPP

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
 * \file   interp_classic.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  Direct interpolation scheme based on \ref Stuben_1999 "Stuben (1999)".
 */

#include <vector>
#include <algorithm>

#include <boost/typeof/typeof.hpp>

#include <amgcl/spmat.hpp>
#include <amgcl/tictoc.hpp>

namespace amgcl {

namespace interp {

/**
 * \defgroup interpolation Interpolation
 * \brief Possible interpolation schemes.
 */

/// Direct interpolation scheme based on \ref Stuben_1999 "Stuben (1999)".
/** \ingroup interpolation */
class classic {
    public:
        /// Parameters controlling direct interpolation.
        /**
         * See \ref Stuben_1999 "Stuben (1999)" for detailed description of
         * these parameters.
         */
        struct params {
            /// Parameter \f$\varepsilon_{str}\f$ defining strong couplings.
            /**
             * Variable \f$i\f$ is defined to be strongly negatively coupled to
             * another variable, \f$j\f$, if \f[-a_{ij} \geq
             * \varepsilon_{str}\max\limits_{a_{ik}<0}|a_{ik}|\quad \text{with
             * fixed} \quad 0 < \varepsilon_{str} < 1.\f]
             * In practice, a value of \f$\varepsilon_{str}=0.25\f$ is usually
             * taken.
             */
            float eps_strong;

            /// Truncate prolongation operator?
            /**
             * Interpolation operators, and, hence coarse operators may
             * increase substabtially towards coarser levels. Without
             * truncation, this may become too costly. Truncation ignores all
             * interpolatory connections which are smaller (in absolute value)
             * than the largest one by a factor of \f$\varepsilon_{tr}\f$. The
             * remaining weights are rescaled so that the total sum remains
             * unchanged. In practice, a value of \f$\varepsilon_{tr}=0.2\f$ is
             * usually taken.
             */
            // TODO: is it possible to say which coefficiants of the P matrix
            // are to be dropped only looking at A matrix?
            bool  trunc_int;

            /// Truncation parameter \f$\varepsilon_{tr}\f$.
            float eps_tr;

            params() : eps_strong(0.25f), trunc_int(true), eps_tr(0.2f) {}
        };

        /// Computes prolongation operator from a system matrix.
        /**
         * \param A   The system matrix.
         * \param prm Parameters.
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

            std::vector<char> cf(n, 'U');

            TIC("conn");
            BOOST_AUTO(S, connect(A, prm, cf));
            TOC("conn");

            TIC("split");
            cfsplit(A, S, cf);
            TOC("split");

            TIC("interpolation");
            index_t nc = 0;
            std::vector<index_t> cidx(n);

            for(index_t i = 0; i < n; i++)
                if (cf[i] == 'C') cidx[i] = nc++;

            BOOST_AUTO(Arow, sparse::matrix_outer_index(A));
            BOOST_AUTO(Acol, sparse::matrix_inner_index(A));
            BOOST_AUTO(Aval, sparse::matrix_values(A));

            std::pair<
                sparse::matrix<value_t, index_t>,
                sparse::matrix<value_t, index_t>
            > PR;

            sparse::matrix<value_t, index_t> &P = PR.first;
            sparse::matrix<value_t, index_t> &R = PR.second;

            P.resize(n, nc);
            std::fill(P.row.begin(), P.row.end(), static_cast<index_t>(0));

            std::vector<value_t> Amin, Amax;

            if (prm.trunc_int) {
                Amin.resize(n);
                Amax.resize(n);
            }

#pragma omp parallel for schedule(dynamic, 1024)
            for(index_t i = 0; i < n; ++i) {
                if (cf[i] == 'C') {
                    ++P.row[i + 1];
                    continue;
                }

                if (prm.trunc_int) {
                    value_t amin = 0, amax = 0;

                    for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j) {
                        if (!S.val[j] || cf[Acol[j]] != 'C') continue;

                        amin = std::min(amin, Aval[j]);
                        amax = std::max(amax, Aval[j]);
                    }

                    Amin[i] = amin = amin * prm.eps_tr;
                    Amax[i] = amax = amax * prm.eps_tr;

                    for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j) {
                        if (!S.val[j] || cf[Acol[j]] != 'C') continue;

                        if (Aval[j] <= amin || Aval[j] >= amax)
                            ++P.row[i + 1];
                    }
                } else {
                    for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j)
                        if (S.val[j] && cf[Acol[j]] == 'C')
                            ++P.row[i + 1];
                }
            }

            std::partial_sum(P.row.begin(), P.row.end(), P.row.begin());

            P.reserve(P.row.back());

#pragma omp parallel for schedule(dynamic, 1024)
            for(index_t i = 0; i < n; ++i) {
                index_t row_head = P.row[i];

                if (cf[i] == 'C') {
                    P.col[row_head] = cidx[i];
                    P.val[row_head] = 1;
                    continue;
                }

                value_t diag  = 0;
                value_t a_num = 0, a_den = 0;
                value_t b_num = 0, b_den = 0;
                value_t d_neg = 0, d_pos = 0;

                for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j) {
                    index_t c = Acol[j];
                    value_t v = Aval[j];

                    if (c == i) {
                        diag = v;
                        continue;
                    }

                    if (v < 0) {
                        a_num += v;
                        if (S.val[j] && cf[c] == 'C') {
                            a_den += v;
                            if (prm.trunc_int && v > Amin[i]) d_neg += v;
                        }
                    } else {
                        b_num += v;
                        if (S.val[j] && cf[c] == 'C') {
                            b_den += v;
                            if (prm.trunc_int && v < Amax[i]) d_pos += v;
                        }
                    }
                }

                value_t cf_neg = 1;
                value_t cf_pos = 1;

                if (prm.trunc_int) {
                    if (fabs(a_den - d_neg) > 1e-32) cf_neg = a_den / (a_den - d_neg);
                    if (fabs(b_den - d_pos) > 1e-32) cf_pos = b_den / (b_den - d_pos);
                }

                if (b_num > 0 && fabs(b_den) < 1e-32) diag += b_num;

                value_t alpha = fabs(a_den) > 1e-32 ? -cf_neg * a_num / (diag * a_den) : 0;
                value_t beta  = fabs(b_den) > 1e-32 ? -cf_pos * b_num / (diag * b_den) : 0;

                for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j) {
                    index_t c = Acol[j];
                    value_t v = Aval[j];

                    if (!S.val[j] || cf[c] != 'C') continue;
                    if (prm.trunc_int && v > Amin[i] && v < Amax[i]) continue;

                    P.col[row_head] = cidx[c];
                    P.val[row_head] = (v < 0 ? alpha : beta) * v;
                    ++row_head;
                }
            }
            TOC("interpolation");

            sparse::transpose(P).swap(R);
            return PR;
        }

    private:
        // Extract strong connections from a system matrix.
        template < class spmat >
        static sparse::matrix<char, typename sparse::matrix_index<spmat>::type>
        connect(const spmat &A, const params &prm, std::vector<char> &cf) {
            typedef typename sparse::matrix_index<spmat>::type index_t;
            typedef typename sparse::matrix_value<spmat>::type value_t;

            const index_t n = sparse::matrix_rows(A);

            sparse::matrix<char, index_t> S(n, n);

            S.row.resize(n + 1, 0);
            S.val.resize(sparse::matrix_nonzeros(A), false);

            BOOST_AUTO(Arow, sparse::matrix_outer_index(A));
            BOOST_AUTO(Acol, sparse::matrix_inner_index(A));
            BOOST_AUTO(Aval, sparse::matrix_values(A));

#pragma omp parallel for schedule(dynamic, 1024)
            for(index_t i = 0; i < n; ++i) {
                value_t a_min = 0;

                for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j)
                    if (Acol[j] != i && Aval[j] < a_min) a_min = Aval[j];

                if (fabs(a_min) < 1e-32) {
                    cf[i] = 'F';
                    continue;
                }

                a_min *= prm.eps_strong;

                for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j)
                    if (Acol[j] != i && Aval[j] < a_min) S.val[j] = true;
            }

            for(index_t i = 0, nnz = Arow[n]; i < nnz; ++i)
                if (S.val[i]) S.row[Acol[i] + 1]++;

            std::partial_sum(S.row.begin(), S.row.end(), S.row.begin());

            S.col.resize(S.row.back());

            for(index_t i = 0; i < n; ++i)
                for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j)
                    if (S.val[j]) S.col[ S.row[ Acol[j] ]++ ] = i;

            for(index_t i = n; i > 0; --i) S.row[i] = S.row[i-1];

            return S;
        }

        // Split variables into C(oarse) and F(ine) sets.
        template < class spmat >
        static void cfsplit(
                const spmat &A,
                const sparse::matrix<char, typename sparse::matrix_index<spmat>::type> &S,
                std::vector<char> &cf
                )
        {
            typedef typename sparse::matrix_index<spmat>::type index_t;

            const index_t n = sparse::matrix_rows(A);

            BOOST_AUTO(Arow, sparse::matrix_outer_index(A));
            BOOST_AUTO(Acol, sparse::matrix_inner_index(A));

            std::vector<index_t> lambda(n);

            // Initialize lambdas:
            for(index_t i = 0; i < n; ++i) {
                index_t temp = 0;
                for(index_t j = S.row[i], e = S.row[i + 1]; j < e; ++j)
                    temp += (cf[S.col[j]] == 'U' ? 1 : 2);
                lambda[i] = temp;
            }

            // Keep track of variable groups of equal lambda values.
            // ptr - start of a group;
            // cnt - size of a group;
            // i2n - variable number;
            // n2i - vaiable position in a group.
            std::vector<index_t> ptr(n+1, static_cast<index_t>(0));
            std::vector<index_t> cnt(n, static_cast<index_t>(0));
            std::vector<index_t> i2n(n);
            std::vector<index_t> n2i(n);

            for(index_t i = 0; i < n; ++i) ptr[lambda[i] + 1]++;

            std::partial_sum(ptr.begin(), ptr.end(), ptr.begin());

            for(index_t i = 0; i < n; ++i) {
                index_t lam = lambda[i];
                index_t idx = ptr[lam] + cnt[lam]++;
                i2n[idx] = i;
                n2i[i] = idx;
            }

            // Process variables by decreasing lambda value.
            // 1. The vaiable with maximum value of lambda becomes next C-variable.
            // 2. Its neighbours from S' become F-variables.
            // 3. Keep lambda values in sync.
            for(index_t top = n - 1; top >= 0; --top) {
                index_t i = i2n[top];
                index_t lam = lambda[i];

                if (lam == 0) {
                    std::replace(cf.begin(), cf.end(), 'U', 'C');
                    break;
                }

                // Remove tne variable from its group.
                cnt[lam]--;

                if (cf[i] == 'F') continue;
                assert(cf[i] == 'U');

                // Mark the variable as 'C'.
                cf[i] = 'C';

                // Its neighbours from S' become F-variables.
                for(index_t j = S.row[i], e = S.row[i + 1]; j < e; ++j) {
                    index_t c = S.col[j];

                    if (cf[c] != 'U') continue;

                    cf[c] = 'F';

                    // Increase lambdas of the newly created F's neighbours.
                    for(index_t jj = Arow[c], ee = Arow[c + 1]; jj < ee; ++jj) {
                        if (!S.val[jj]) continue;

                        index_t cc = Acol[jj];
                        index_t lam_cc = lambda[cc];

                        if (cf[cc] != 'U' || lam_cc >= n - 1) continue;

                        index_t old_pos = n2i[cc];
                        index_t new_pos = ptr[lam_cc] + cnt[lam_cc] - 1;

                        n2i[i2n[old_pos]] = new_pos;
                        n2i[i2n[new_pos]] = old_pos;

                        std::swap(i2n[old_pos], i2n[new_pos]);

                        --cnt[lam_cc];
                        ++cnt[lam_cc + 1];
                        ptr[lam_cc + 1] = ptr[lam_cc] + cnt[lam_cc];

                        ++lambda[cc];
                    }
                }

                // Decrease lambdas of the newly create C's neighbours.
                for(index_t j = Arow[i], e = Arow[i + 1]; j < e; j++) {
                    if (!S.val[j]) continue;

                    index_t c = Acol[j];
                    index_t lam = lambda[c];

                    if (cf[c] != 'U' || lam == 0) continue;

                    index_t old_pos = n2i[c];
                    index_t new_pos = ptr[lam];

                    n2i[i2n[old_pos]] = new_pos;
                    n2i[i2n[new_pos]] = old_pos;

                    std::swap(i2n[old_pos], i2n[new_pos]);

                    --cnt[lam];
                    ++cnt[lam - 1];
                    ++ptr[lam];
                    --lambda[c];

                    assert(ptr[lam - 1] == ptr[lam] - cnt[lam - 1]);
                }
            }
        }
};

} // namespace interp

} // namespace amgcl

#endif
