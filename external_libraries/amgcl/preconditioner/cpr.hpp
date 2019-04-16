#ifndef AMGCL_PRECONDITIONER_CPR_HPP
#define AMGCL_PRECONDITIONER_CPR_HPP

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
 * \file   amgcl/preconditioner/cpr.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Two stage preconditioner of the Constrained Pressure Residual type.
 */

#include <vector>
#include <memory>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace preconditioner {

template <class PPrecond, class SPrecond>
class cpr {
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

            params() : block_size(2), active_rows(0) {}

#ifndef AMGCL_NO_BOOST
            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_CHILD(p, pprecond),
                  AMGCL_PARAMS_IMPORT_CHILD(p, sprecond),
                  AMGCL_PARAMS_IMPORT_VALUE(p, block_size),
                  AMGCL_PARAMS_IMPORT_VALUE(p, active_rows)
            {
                check_params(p, {"pprecond", "sprecond", "block_size", "active_rows"});
            }

            void get(boost::property_tree::ptree &p, const std::string &path = "") const
            {
                AMGCL_PARAMS_EXPORT_CHILD(p, path, pprecond);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, sprecond);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, block_size);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, active_rows);
            }
#endif
        } prm;

        template <class Matrix>
        cpr(
                const Matrix &K,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
           ) : prm(prm), n(backend::rows(K))
        {
            init(std::make_shared<build_matrix>(K), bprm);
        }

        cpr(
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

            np = N / B;

            auto fpp = std::make_shared<build_matrix>();
            fpp->set_size(np, n);
            fpp->set_nonzeros(n);
            fpp->ptr[0] = 0;

            auto App = std::make_shared<build_matrix>();
            App->set_size(np, np, true);

            // Get the pressure matrix nonzero pattern,
            // extract and invert block diagonals.
#pragma omp parallel
            {
                std::vector<row_iterator> k; k.reserve(B);
                multi_array<value_type, 2> v(B, B);

#pragma omp for
                for(ptrdiff_t ip = 0; ip < static_cast<ptrdiff_t>(np); ++ip) {
                    ptrdiff_t ik = ip * B;
                    bool      done = true;
                    ptrdiff_t cur_col = 0;

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

                        fpp->col[ik + i] = ik + i;
                    }
                    fpp->ptr[ip+1] = ik + B;

                    while (!done) {
                        ++App->ptr[ip+1];

                        ptrdiff_t end = (cur_col + 1) * B;

                        if (cur_col == ip) {
                            // This is diagonal block.
                            // Capture its (transposed) value,
                            // invert it and put the relevant row into fpp.
                            for(int i = 0; i < B; ++i)
                                for(int j = 0; j < B; ++j) v(i,j) = 0;

                            for(int i = 0; i < B; ++i)
                                for(; k[i] && k[i].col() < end; ++k[i])
                                    v(k[i].col() % B, i) = k[i].value();

                            invert(v, &fpp->val[ik]);
                        } else {
                            // This is off-diagonal block.
                            // Just skip it.
                            for(int i = 0; i < B; ++i)
                                while(k[i] && k[i].col() < end) ++k[i];
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
                    ptrdiff_t cur_col;

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

        // Inverts dense matrix A;
        // Returns the first column of the inverted matrix.
        void invert(multi_array<value_type, 2> &A, value_type *y)
        {
            const int B = prm.block_size;

            // Perform LU-factorization of A in-place
            for(int k = 0; k < B; ++k) {
                value_type d = A(k,k);
                assert(!math::is_zero(d));
                for(int i = k+1; i < B; ++i) {
                    A(i,k) /= d;
                    for(int j = k+1; j < B; ++j)
                        A(i,j) -= A(i,k) * A(k,j);
                }
            }

            // Invert unit vector in-place.
            // Lower triangular solve:
            for(int i = 0; i < B; ++i) {
                value_type b = static_cast<value_type>(i == 0);
                for(int j = 0; j < i; ++j)
                    b -= A(i,j) * y[j];
                y[i] = b;
            }

            // Upper triangular solve:
            for(int i = B; i --> 0; ) {
                for(int j = i+1; j < B; ++j)
                    y[i] -= A(i,j) * y[j];
                y[i] /= A(i,i);
            }
        }

        friend std::ostream& operator<<(std::ostream &os, const cpr &p) {
            os << "CPR (two-stage preconditioner)\n"
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
