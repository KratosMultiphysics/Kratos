#ifndef AMGCL_PRECONDITIONER_CPR_DRS_HPP
#define AMGCL_PRECONDITIONER_CPR_DRS_HPP

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
 * \file   amgcl/preconditioner/cpr_drs.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  CPR preconditioner with Dynamic Row Sum modification.
 */

#include <vector>
#include <algorithm>

#include <memory>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace preconditioner {

template <class PPrecond, class SPrecond>
class cpr_drs {
    static_assert(
            std::is_same<
                typename PPrecond::backend_type,
                typename SPrecond::backend_type
                >::value,
            "Backends for pressure and flow preconditioners should coinside!"
            );
    public:
        typedef typename PPrecond::backend_type backend_type;

        typedef typename backend_type::value_type value_type;
        typedef typename backend_type::matrix     matrix;
        typedef typename backend_type::vector     vector;
        typedef typename backend_type::params     backend_params;

        typedef typename backend::builtin<value_type>::matrix build_matrix;

        struct params {
            typedef typename PPrecond::params pprecond_params;
            typedef typename SPrecond::params sprecond_params;

            pprecond_params pprecond;
            sprecond_params sprecond;

            int    block_size;
            size_t active_rows;
            double eps_dd;
            double eps_ps;

            std::vector<double> weights;

            params() : block_size(2), active_rows(0), eps_dd(0.2), eps_ps(0.02) {}

#ifndef AMGCL_NO_BOOST
            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_CHILD(p, pprecond),
                  AMGCL_PARAMS_IMPORT_CHILD(p, sprecond),
                  AMGCL_PARAMS_IMPORT_VALUE(p, block_size),
                  AMGCL_PARAMS_IMPORT_VALUE(p, active_rows),
                  AMGCL_PARAMS_IMPORT_VALUE(p, eps_dd),
                  AMGCL_PARAMS_IMPORT_VALUE(p, eps_ps)
            {
                void  *ptr = 0;
                size_t n   = 0;

                ptr = p.get("weights",      ptr);
                n   = p.get("weights_size", n);

                if (ptr) {
                    precondition(n > 0,
                            "Error in cpr_wdrs parameters: "
                            "weights is set, but weights_size is not"
                            );

                    weights.assign(
                            static_cast<double*>(ptr),
                            static_cast<double*>(ptr) + n);
                }

                check_params(p, {"pprecond", "sprecond", "block_size", "active_rows", "eps_dd", "eps_ps", "weights", "weights_size"});
            }

            void get(boost::property_tree::ptree &p, const std::string &path = "") const
            {
                AMGCL_PARAMS_EXPORT_CHILD(p, path, pprecond);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, sprecond);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, block_size);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, active_rows);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, eps_dd);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, eps_ps);
            }
#endif
        } prm;

        template <class Matrix>
        cpr_drs(
                const Matrix &K,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
               ) : prm(prm), n(backend::rows(K))
        {
            init(std::make_shared<build_matrix>(K), bprm);
        }

        cpr_drs(
                std::shared_ptr<build_matrix> K,
                const params &prm = params(),
                const backend_params bprm = backend_params()
               ) : prm(prm), n(backend::rows(*K))
        {
            init(K, bprm);
        }

        template <class Vec1, class Vec2>
        void apply(const Vec1 &rhs, Vec2 &&x) const {
            S->apply(rhs, x);
            backend::residual(rhs, S->system_matrix(), x, *rs);

            backend::spmv(1, *Fpp, *rs, 0, *rp);
            P->apply(*rp, *xp);

            backend::spmv(1, *Scatter, *xp, 1, x);
        }

        std::shared_ptr<matrix> system_matrix_ptr() const {
            return S->system_matrix_ptr();
        }

        const matrix& system_matrix() const {
            return S->system_matrix();
        }

    private:
        size_t n, np;

        std::shared_ptr<PPrecond> P;
        std::shared_ptr<SPrecond> S;

        std::shared_ptr<matrix> Fpp, Scatter;
        std::shared_ptr<vector> rp, xp, rs;

        void init(std::shared_ptr<build_matrix> K, const backend_params bprm)
        {
            typedef typename backend::row_iterator<build_matrix>::type row_iterator;
            const int       B = prm.block_size;
            const ptrdiff_t N = (prm.active_rows ? prm.active_rows : n);

            precondition(
                    prm.weights.empty() || prm.weights.size() == static_cast<size_t>(N),
                    "CPR: weights size is not equal to number of active rows.");

            np = N / B;

            auto fpp = std::make_shared<build_matrix>();
            fpp->set_size(np, n);
            fpp->set_nonzeros(n);
            fpp->ptr[0] = 0;

            auto App = std::make_shared<build_matrix>();
            App->set_size(np, np, true);

#pragma omp parallel
            {
                std::vector<value_type> a_dia(B), a_off(B), a_top(B);
                std::vector<row_iterator> k; k.reserve(B);

#pragma omp for
                for(ptrdiff_t ip = 0; ip < static_cast<ptrdiff_t>(np); ++ip) {
                    ptrdiff_t ik = ip * B;
                    bool      done = true;
                    ptrdiff_t cur_col = 0;

                    std::fill(a_dia.begin(), a_dia.end(), 0);
                    std::fill(a_off.begin(), a_off.end(), 0);
                    std::fill(a_top.begin(), a_top.end(), 0);

                    k.clear();
                    for(int i = 0; i < B; ++i) {
                        k.push_back(backend::row_begin(*K, ik + i));

                        if (k.back() && k.back().col() < N) {
                            ptrdiff_t col = k.back().col() / B;
                            if (done) {
                                cur_col = col;
                                done = false;
                            } else {
                                cur_col = std::min(cur_col, col);
                            }
                        }
                    }

                    while (!done) {
                        ++App->ptr[ip+1];

                        ptrdiff_t end = (cur_col + 1) * B;

                        for(int i = 0; i < B; ++i) {
                            for(; k[i] && k[i].col() < end; ++k[i]) {
                                ptrdiff_t  c = k[i].col() % B;
                                value_type v = k[i].value();

                                if (i == 0) {
                                    a_top[c] += std::abs(v);
                                }

                                if (c == 0) {
                                    if (cur_col == ip) {
                                        a_dia[i] = v;
                                    } else {
                                        a_off[i] += std::abs(v);
                                    }
                                }
                            }
                        }

                        // Get next column number.
                        done = true;
                        for(int i = 0; i < B; ++i) {
                            if (k[i] && k[i].col() < N) {
                                ptrdiff_t col = k[i].col() / B;
                                if (done) {
                                    cur_col = col;
                                    done = false;
                                } else {
                                    cur_col = std::min(cur_col, col);
                                }
                            }
                        }
                    }

                    for(int i = 0; i < B; ++i) {
                        fpp->col[ik+i] = ik+i;
                        double delta = 1;

                        if (!prm.weights.empty())
                            delta *= prm.weights[ik+i];

                        if (i > 0) {
                            if (a_dia[i] < prm.eps_dd * a_off[i])
                                delta = 0;

                            if (a_top[i] < prm.eps_ps * std::abs(a_dia[0]))
                                delta = 0;
                        }

                        fpp->val[ik+i] = delta;
                    }

                    fpp->ptr[ip+1] = ik + B;
                }
            }

            App->set_nonzeros(App->scan_row_sizes());

            auto scatter = std::make_shared<build_matrix>();
            scatter->set_size(n, np);
            scatter->set_nonzeros(np);
            scatter->ptr[0] = 0;

#pragma omp parallel
            {
                std::vector<row_iterator> k; k.reserve(B);

#pragma omp for
                for(ptrdiff_t ip = 0; ip < static_cast<ptrdiff_t>(np); ++ip) {
                    ptrdiff_t ik = ip * B;
                    ptrdiff_t head = App->ptr[ip];
                    bool      done = true;
                    ptrdiff_t cur_col = 0;

                    value_type *d = &fpp->val[ik];

                    k.clear();
                    for(int i = 0; i < B; ++i) {
                        k.push_back(backend::row_begin(*K, ik + i));

                        if (k.back() && k.back().col() < N) {
                            ptrdiff_t col = k.back().col() / B;
                            if (done) {
                                cur_col = col;
                                done = false;
                            } else {
                                cur_col = std::min(cur_col, col);
                            }
                        }
                    }

                    while (!done) {
                        ptrdiff_t  end = (cur_col + 1) * B;
                        value_type app = 0;

                        for(int i = 0; i < B; ++i) {
                            for(; k[i] && k[i].col() < end; ++k[i]) {
                                if (k[i].col() % B == 0) {
                                    app += d[i] * k[i].value();
                                }
                            }
                        }

                        App->col[head] = cur_col;
                        App->val[head] = app;
                        ++head;

                        // Get next column number.
                        done = true;
                        for(int i = 0; i < B; ++i) {
                            if (k[i] && k[i].col() < N) {
                                ptrdiff_t col = k[i].col() / B;
                                if (done) {
                                    cur_col = col;
                                    done = false;
                                } else {
                                    cur_col = std::min(cur_col, col);
                                }
                            }
                        }
                    }

                    scatter->col[ip] = ip;
                    scatter->val[ip] = math::identity<value_type>();

                    ptrdiff_t nnz = ip;
                    for(int i = 0; i < B; ++i) {
                        if (i == 0) ++nnz;
                        scatter->ptr[ik + i + 1] = nnz;
                    }
                }
            }

            for(size_t i = N; i < n; ++i)
                scatter->ptr[i+1] = scatter->ptr[i];

            P = std::make_shared<PPrecond>(App, prm.pprecond, bprm);
            S = std::make_shared<SPrecond>(K,   prm.sprecond, bprm);

            Fpp     = backend_type::copy_matrix(fpp, bprm);
            Scatter = backend_type::copy_matrix(scatter, bprm);

            rp = backend_type::create_vector(np, bprm);
            xp = backend_type::create_vector(np, bprm);
            rs = backend_type::create_vector(n, bprm);
        }

        friend std::ostream& operator<<(std::ostream &os, const cpr_drs &p) {
            os << "CPR_DRS (two-stage preconditioner)\n"
                  "### Pressure preconditioner:\n"
               << *p.P << "\n"
                  "### Global preconditioner:\n"
               << *p.S << std::endl;
            return os;
        }
};

} // namespace preconditioner
} // namespace amgcl


#endif
