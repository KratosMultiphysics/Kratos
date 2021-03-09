#ifndef AMGCL_PRECONDITIONER_CPR_HPP
#define AMGCL_PRECONDITIONER_CPR_HPP

/*
The MIT License

Copyright (c) 2012-2020 Denis Demidov <dennis.demidov@gmail.com>

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
            math::static_rows<typename PPrecond::backend_type::value_type>::value == 1,
            "Pressure backend should have scalar value type!"
            );

    static_assert(
            backend::backends_compatible<
                typename PPrecond::backend_type,
                typename SPrecond::backend_type
                >::value,
            "Backends for pressure and flow preconditioners should coincide!"
            );

    public:
        typedef typename SPrecond::backend_type backend_type;
        typedef typename PPrecond::backend_type backend_type_p;

        typedef typename backend_type::value_type   value_type;
        typedef typename backend_type::matrix       matrix;
        typedef typename backend_type::vector       vector;
        typedef typename backend_type_p::value_type value_type_p;
        typedef typename backend_type_p::matrix     matrix_p;
        typedef typename backend_type_p::vector     vector_p;

        typedef typename backend_type::params       backend_params;

        typedef typename backend::builtin<value_type>::matrix   build_matrix;
        typedef typename backend::builtin<value_type_p>::matrix build_matrix_p;

        typedef typename math::scalar_of<value_type>::type scalar_type;

        struct params {
            typedef typename PPrecond::params pprecond_params;
            typedef typename SPrecond::params sprecond_params;

            pprecond_params pprecond;
            sprecond_params sprecond;

            int    block_size;
            size_t active_rows;

            params()
                : block_size(math::static_rows<value_type>::value == 1 ? 2 : math::static_rows<value_type>::value),
                  active_rows(0) {}

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
            init(std::make_shared<build_matrix>(K), bprm,
                    std::integral_constant<bool, math::static_rows<value_type>::value == 1>());
        }

        cpr(
                std::shared_ptr<build_matrix> K,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
           ) : prm(prm), n(backend::rows(*K))
        {
            init(K, bprm,
                    std::integral_constant<bool, math::static_rows<value_type>::value == 1>());
        }

        template <class Vec1, class Vec2>
        void apply(const Vec1 &rhs, Vec2 &&x) const {
            AMGCL_TIC("sprecond");
            S->apply(rhs, x);
            AMGCL_TOC("sprecond");
            backend::residual(rhs, S->system_matrix(), x, *rs);

            backend::spmv(1, *Fpp, *rs, 0, *rp);
            AMGCL_TIC("pprecond");
            P->apply(*rp, *xp);
            AMGCL_TOC("pprecond");

            backend::spmv(1, *Scatter, *xp, 1, x);
        }

        std::shared_ptr<matrix> system_matrix_ptr() const {
            return S->system_matrix_ptr();
        }

        const matrix& system_matrix() const {
            return S->system_matrix();
        }

        template <class Matrix>
        void partial_update(
                const Matrix &K,
                bool update_transfer_ops = true,
                const backend_params &bprm = backend_params()
              )
        {
            auto K_ptr = std::make_shared<build_matrix>(K);
            // Update global preconditioner
            S = std::make_shared<SPrecond>(K_ptr, prm.sprecond, bprm);
            if(update_transfer_ops){
              // Update transfer operator Fpp
              update_transfer(
                  K_ptr,
                  bprm,
                  std::integral_constant<bool, math::static_rows<value_type>::value == 1>()
                );
            }
        }

    private:
        size_t n, np;

        std::shared_ptr<PPrecond> P;
        std::shared_ptr<SPrecond> S;

        std::shared_ptr<matrix_p> Fpp, Scatter;
        std::shared_ptr<vector>   rs;
        std::shared_ptr<vector_p> rp, xp;

        // Returns pressure transfer operator fpp and (optionally)
        // partially constructed pressure system matrix App.
        std::tuple<std::shared_ptr<build_matrix_p>, std::shared_ptr<build_matrix_p>>
        first_scalar_pass(std::shared_ptr<build_matrix> K, bool get_app = true) {
            typedef typename backend::row_iterator<build_matrix>::type row_iterator;
            const int       B = prm.block_size;
            const ptrdiff_t N = (prm.active_rows ? prm.active_rows : n);

            np = N / B;

            auto fpp = std::make_shared<build_matrix_p>();
            fpp->set_size(np, N);
            fpp->set_nonzeros(N);
            fpp->ptr[0] = 0;

            std::shared_ptr<build_matrix_p> App;

            if (get_app) {
                App = std::make_shared<build_matrix>();
                App->set_size(np, np, true);
            }

            // Get the pressure matrix nonzero pattern,
            // extract and invert block diagonals.
#pragma omp parallel
            {
                std::vector<row_iterator> k; k.reserve(B);
                multi_array<scalar_type, 2> v(B, B);

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
                        if (get_app) ++App->ptr[ip+1];

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

                            invert(v.data(), &fpp->val[ik]);

                            if (!get_app) break;
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

            if (get_app)
                App->set_nonzeros(App->scan_row_sizes());

            return std::make_tuple(fpp, App);
        }

        // The system matrix has scalar values
        void init(std::shared_ptr<build_matrix> K, const backend_params bprm, std::true_type)
        {
            typedef typename backend::row_iterator<build_matrix>::type row_iterator;
            const int       B = prm.block_size;
            const ptrdiff_t N = (prm.active_rows ? prm.active_rows : n);

            np = N / B;

            std::shared_ptr<build_matrix_p> fpp, App;
            std::tie(fpp, App) = first_scalar_pass(K);

            auto scatter = std::make_shared<build_matrix_p>();
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

                    value_type_p *d = &fpp->val[ik];

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
                        ptrdiff_t    end = (cur_col + 1) * B;
                        value_type_p app = 0;

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
                    scatter->val[ip] = math::identity<value_type_p>();

                    ptrdiff_t nnz = ip;
                    for(int i = 0; i < B; ++i) {
                        if (i == 0) ++nnz;
                        scatter->ptr[ik + i + 1] = nnz;
                    }
                }
            }

            for(size_t i = N; i < n; ++i)
                scatter->ptr[i+1] = scatter->ptr[i];

            AMGCL_TIC("pprecond");
            P = std::make_shared<PPrecond>(App, prm.pprecond, bprm);
            AMGCL_TOC("pprecond");
            AMGCL_TIC("sprecond");
            S = std::make_shared<SPrecond>(K,   prm.sprecond, bprm);
            AMGCL_TOC("sprecond");

            Fpp     = backend_type_p::copy_matrix(fpp, bprm);
            Scatter = backend_type_p::copy_matrix(scatter, bprm);

            rp = backend_type_p::create_vector(np, bprm);
            xp = backend_type_p::create_vector(np, bprm);
            rs = backend_type::create_vector(n, bprm);
        }

        void update_transfer(std::shared_ptr<build_matrix> K, const backend_params bprm, std::true_type)
        {
            auto fpp = std::get<0>(first_scalar_pass(K, /*get_app*/false));
            Fpp = backend_type_p::copy_matrix(fpp, bprm);
        }

        // The system matrix has block values
        void init(std::shared_ptr<build_matrix> K, const backend_params bprm, std::false_type)
        {
            const int       B = math::static_rows<value_type>::value;
            const ptrdiff_t N = (prm.active_rows ? prm.active_rows : n);

            np = N;

            auto fpp = std::make_shared<build_matrix_p>();
            fpp->set_size(np, np * B);
            fpp->set_nonzeros(np * B);
            fpp->ptr[0] = 0;

            auto scatter = std::make_shared<build_matrix_p>();
            scatter->set_size(np * B, np);
            scatter->set_nonzeros(np);
            scatter->ptr[0] = 0;

            auto App = std::make_shared<build_matrix_p>();
            App->set_size(np, np, true);
            App->set_nonzeros(K->nnz);
            App->ptr[0] = 0;

#pragma omp parallel for
            for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(np); ++i) {
                ptrdiff_t ik = i * B;
                for(int k = 0; k < B; ++k, ++ik) {
                    fpp->col[ik] = ik;
                    scatter->ptr[ik + 1] = i + 1;
                }

                fpp->ptr[i + 1] = ik;
                scatter->col[i] = i;
                scatter->val[i] = math::identity<value_type_p>();

                ptrdiff_t row_beg = K->ptr[i];
                ptrdiff_t row_end = K->ptr[i + 1];
                App->ptr[i+1] = row_end;

                // Extract and invert block diagonals
                value_type_p *d = &fpp->val[i * B];
                for(ptrdiff_t j = row_beg; j < row_end; ++j) {
                    if (K->col[j] == i) {
                        value_type v = math::adjoint(K->val[j]);
                        invert(v.data(), d);
                        break;
                    }
                }

                for(ptrdiff_t j = row_beg; j < row_end; ++j) {
                    value_type_p app = 0;
                    for(int k = 0; k < B; ++k)
                        app += d[k] * K->val[j](k,0);

                    App->col[j] = K->col[j];
                    App->val[j] = app;
                }
            }

            AMGCL_TIC("pprecond");
            P = std::make_shared<PPrecond>(App, prm.pprecond, bprm);
            AMGCL_TOC("pprecond");
            AMGCL_TIC("sprecond");
            S = std::make_shared<SPrecond>(K,   prm.sprecond, bprm);
            AMGCL_TOC("sprecond");


            Fpp     = backend_type_p::copy_matrix(fpp, bprm);
            Scatter = backend_type_p::copy_matrix(scatter, bprm);

            rp = backend_type_p::create_vector(np, bprm);
            xp = backend_type_p::create_vector(np, bprm);
            rs = backend_type::create_vector(n, bprm);
        }

        void update_transfer(std::shared_ptr<build_matrix> K, const backend_params bprm, std::false_type)
        {
            const int       B = math::static_rows<value_type>::value;
            const ptrdiff_t N = (prm.active_rows ? prm.active_rows : n);

            np = N;

            auto fpp = std::make_shared<build_matrix_p>();
            fpp->set_size(np, np * B);
            fpp->set_nonzeros(np * B);
            fpp->ptr[0] = 0;

#pragma omp parallel for
            for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(np); ++i) {
                ptrdiff_t ik = i * B;
                for(int k = 0; k < B; ++k, ++ik) {
                    fpp->col[ik] = ik;
                }

                fpp->ptr[i + 1] = ik;

                ptrdiff_t row_beg = K->ptr[i];
                ptrdiff_t row_end = K->ptr[i + 1];

                // Extract and invert block diagonals
                value_type_p *d = &fpp->val[i * B];
                for(ptrdiff_t j = row_beg; j < row_end; ++j) {
                    if (K->col[j] == i) {
                        value_type v = math::adjoint(K->val[j]);
                        invert(v.data(), d);
                        break;
                    }
                }
            }
            Fpp = backend_type_p::copy_matrix(fpp, bprm);
        }

        // Inverts dense matrix A;
        // Returns the first column of the inverted matrix.
        void invert(scalar_type *A, value_type_p *y)
        {
            const int B = math::static_rows<value_type>::value == 1 ? prm.block_size : math::static_rows<value_type>::value;

            // Perform LU-factorization of A in-place
            for(int k = 0; k < B; ++k) {
                scalar_type d = A[k*B+k];
                assert(!math::is_zero(d));
                for(int i = k+1; i < B; ++i) {
                    A[i*B+k] /= d;
                    for(int j = k+1; j < B; ++j)
                        A[i*B+j] -= A[i*B+k] * A[k*B+j];
                }
            }

            // Invert unit vector in-place.
            // Lower triangular solve:
            for(int i = 0; i < B; ++i) {
                value_type_p b = static_cast<value_type_p>(i == 0);
                for(int j = 0; j < i; ++j)
                    b -= A[i*B+j] * y[j];
                y[i] = b;
            }

            // Upper triangular solve:
            for(int i = B; i --> 0; ) {
                for(int j = i+1; j < B; ++j)
                    y[i] -= A[i*B+j] * y[j];
                y[i] /= A[i*B+i];
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
