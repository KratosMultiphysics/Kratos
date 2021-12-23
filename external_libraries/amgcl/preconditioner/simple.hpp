#ifndef AMGCL_PRECONDITIONER_SIMPLE_HPP
#define AMGCL_PRECONDITIONER_SIMPLE_HPP

/*
The MIT License

Copyright (c) 2012-2020 Denis Demidov <dennis.demidov@gmail.com>
Copyright (c) 2016, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)

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
 * \file   amgcl/preconditioner/simple.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  SIMPLE[C] preconditioning scheme.
 *
 * [1] Elman, Howard, et al. "A taxonomy and comparison of parallel block
 *     multi-level preconditioners for the incompressible Navier–Stokes
 *     equations." Journal of Computational Physics 227.3 (2008): 1790-1808.
 */

#include <vector>

#include <memory>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/backend/detail/mixing.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace preconditioner {

/// Schur-complement pressure correction preconditioner
template <class USolver, class PSolver>
class simple {
    static_assert(
            backend::backends_compatible<
                typename USolver::backend_type,
                typename PSolver::backend_type
                >::value,
            "Backends for pressure and flow preconditioners should coincide!"
            );
    public:
        typedef
            typename backend::detail::common_scalar_backend<
                typename USolver::backend_type,
                typename PSolver::backend_type
                >::type
            backend_type;

        typedef typename backend_type::value_type value_type;
        typedef typename backend_type::matrix     matrix;
        typedef typename backend_type::vector     vector;
        typedef typename backend_type::params     backend_params;

        typedef typename backend::builtin<value_type>::matrix build_matrix;

        struct params {
            typedef typename USolver::params usolver_params;
            typedef typename PSolver::params psolver_params;

            usolver_params usolver;
            psolver_params psolver;

            std::vector<char> pmask;

            // Use 1/sum_j(abs(Kuu_{i,j})) instead of dia(Kuu)^-1
            // as approximation for the Kuu^-1 (as in SIMPLEC algorithm)
            bool simplec;

            params() : simplec(true) {}

#ifndef AMGCL_NO_BOOST
            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_CHILD(p, usolver),
                  AMGCL_PARAMS_IMPORT_CHILD(p, psolver),
                  AMGCL_PARAMS_IMPORT_VALUE(p, simplec)
            {
                size_t n = 0;

                n = p.get("pmask_size", n);

                precondition(n > 0, "Error in SIMPLE parameters: pmask_size is not set");

                if (p.count("pmask_pattern")) {
                    pmask.resize(n, 0);

                    std::string pattern = p.get("pmask_pattern", std::string());
                    switch (pattern[0]) {
                        case '%':
                            {
                                int start  = std::atoi(pattern.substr(1).c_str());
                                int stride = std::atoi(pattern.substr(3).c_str());
                                for(size_t i = start; i < n; i += stride) pmask[i] = 1;
                            }
                            break;
                        case '<':
                            {
                                size_t m = std::atoi(pattern.c_str()+1);
                                for(size_t i = 0; i < std::min(m, n); ++i) pmask[i] = 1;
                            }
                            break;
                        case '>':
                            {
                                size_t m = std::atoi(pattern.c_str()+1);
                                for(size_t i = m; i < n; ++i) pmask[i] = 1;
                            }
                            break;
                        default:
                            precondition(false, "Unknown pattern in pmask_pattern");
                    }
                } else if (p.count("pmask")) {
                    void *pm = 0;
                    pm = p.get("pmask", pm);
                    pmask.assign(static_cast<char*>(pm), static_cast<char*>(pm) + n);
                } else {
                    precondition(false,
                            "Error in SIMPLE parameters: "
                            "neither pmask_pattern, nor pmask is set"
                            );
                }

                check_params(p, {"usolver", "psolver", "simplec", "pmask_size"},
                        {"pmask", "pmask_pattern"});
            }

            void get(boost::property_tree::ptree &p, const std::string &path = "") const
            {
                AMGCL_PARAMS_EXPORT_CHILD(p, path, usolver);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, psolver);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, simplec);
            }
#endif
        } prm;

        template <class Matrix>
        simple(
                const Matrix &K,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
            : prm(prm), n(backend::rows(K)), np(0), nu(0)
        {
            init(std::make_shared<build_matrix>(K), bprm);
        }

        simple(
                std::shared_ptr<build_matrix> K,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
            : prm(prm), n(backend::rows(*K)), np(0), nu(0)
        {
            init(K, bprm);
        }

        template <class Vec1, class Vec2>
        void apply(const Vec1 &rhs, Vec2 &&x) const {
            backend::spmv(1, *x2u, rhs, 0, *rhs_u);
            backend::spmv(1, *x2p, rhs, 0, *rhs_p);

            // 1. Kuu u_1/2 = rhs_u
            backend::clear(*u);
            report("U", (*U)(*rhs_u, *u));

            // 2. S p = Kpu u_1/2 - rhs_p
            backend::clear(*p);
            backend::spmv(1, *Kpu, *u, -1, *rhs_p);
            report("P", (*P)(*rhs_p, *p));

            // 3. u = u_1/2 - M Kup p
            backend::spmv(1, *Kup, *p, 0, *rhs_u);
            backend::vmul(-1, *M, *rhs_u, 1, *u);

            backend::clear(x);
            backend::spmv(1, *u2x, *u, 1, x);
            backend::spmv(1, *p2x, *p, 1, x);
        }

        std::shared_ptr<matrix> system_matrix_ptr() const {
            return K;
        }

        const matrix& system_matrix() const {
            return *K;
        }

        size_t bytes() const {
            size_t b = 0;

            b += backend::bytes(*K);
            b += backend::bytes(*Kup);
            b += backend::bytes(*Kpu);
            b += backend::bytes(*x2u);
            b += backend::bytes(*x2p);
            b += backend::bytes(*u2x);
            b += backend::bytes(*p2x);
            b += backend::bytes(*rhs_u);
            b += backend::bytes(*rhs_p);
            b += backend::bytes(*u);
            b += backend::bytes(*p);
            b += backend::bytes(*U);
            b += backend::bytes(*P);
            b += backend::bytes(*M);

            return b;
        }

    private:
        size_t n, np, nu;

        std::shared_ptr<matrix> K, Kup, Kpu, x2u, x2p, u2x, p2x;
        std::shared_ptr<vector> rhs_u, rhs_p, u, p;
        std::shared_ptr<typename backend_type::matrix_diagonal> M;

        std::shared_ptr<USolver> U;
        std::shared_ptr<PSolver> P;

        void init(const std::shared_ptr<build_matrix> &K, const backend_params &bprm)
        {
            this->K = backend_type::copy_matrix(K, bprm);

            // Extract matrix subblocks.
            auto Kuu = std::make_shared<build_matrix>();
            auto Kpu = std::make_shared<build_matrix>();
            auto Kup = std::make_shared<build_matrix>();
            auto Kpp = std::make_shared<build_matrix>();

            std::vector<ptrdiff_t> idx(n);

            for(size_t i = 0; i < n; ++i)
                idx[i] = (prm.pmask[i] ? np++ : nu++);

            Kuu->set_size(nu, nu, true);
            Kup->set_size(nu, np, true);
            Kpu->set_size(np, nu, true);
            Kpp->set_size(np, np, true);

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                ptrdiff_t ci = idx[i];
                char      pi = prm.pmask[i];
                for(auto k = backend::row_begin(*K, i); k; ++k) {
                    char pj = prm.pmask[k.col()];

                    if (pi) {
                        if (pj) {
                            ++Kpp->ptr[ci+1];
                        } else {
                            ++Kpu->ptr[ci+1];
                        }
                    } else {
                        if (pj) {
                            ++Kup->ptr[ci+1];
                        } else {
                            ++Kuu->ptr[ci+1];
                        }
                    }
                }
            }

            Kuu->set_nonzeros(Kuu->scan_row_sizes());
            Kup->set_nonzeros(Kup->scan_row_sizes());
            Kpu->set_nonzeros(Kpu->scan_row_sizes());
            Kpp->set_nonzeros(Kpp->scan_row_sizes());

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                ptrdiff_t ci = idx[i];
                char      pi = prm.pmask[i];

                ptrdiff_t uu_head = 0, up_head = 0, pu_head = 0, pp_head = 0;

                if(pi) {
                    pu_head = Kpu->ptr[ci];
                    pp_head = Kpp->ptr[ci];
                } else {
                    uu_head = Kuu->ptr[ci];
                    up_head = Kup->ptr[ci];
                }

                for(auto k = backend::row_begin(*K, i); k; ++k) {
                    ptrdiff_t  j = k.col();
                    value_type v = k.value();
                    ptrdiff_t cj = idx[j];
                    char      pj = prm.pmask[j];

                    if (pi) {
                        if (pj) {
                            Kpp->col[pp_head] = cj;
                            Kpp->val[pp_head] = v;
                            ++pp_head;
                        } else {
                            Kpu->col[pu_head] = cj;
                            Kpu->val[pu_head] = v;
                            ++pu_head;
                        }
                    } else {
                        if (pj) {
                            Kup->col[up_head] = cj;
                            Kup->val[up_head] = v;
                            ++up_head;
                        } else {
                            Kuu->col[uu_head] = cj;
                            Kuu->val[uu_head] = v;
                            ++uu_head;
                        }
                    }
                }
            }

            std::shared_ptr<backend::numa_vector<value_type>> Kuu_dia;

            if (prm.simplec) {
                Kuu_dia = std::make_shared<backend::numa_vector<value_type>>(nu);
#pragma omp parallel for
                for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nu); ++i) {
                    value_type s = math::zero<value_type>();
                    for(ptrdiff_t j = Kuu->ptr[i], e = Kuu->ptr[i+1]; j < e; ++j) {
                        s += math::norm(Kuu->val[j]);
                    }
                    (*Kuu_dia)[i] = math::inverse(s);
                }
            } else {
                Kuu_dia = diagonal(*Kuu, /*invert = */true);
            }

            this->Kup = backend_type::copy_matrix(Kup, bprm);
            this->Kpu = backend_type::copy_matrix(Kpu, bprm);

            // Kpp -= Kpu * dia(Kuu)^-1 * Kup)
            {
#pragma omp parallel for
                for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(np); ++i) {
                    for(ptrdiff_t j = Kpu->ptr[i], e = Kpu->ptr[i+1]; j < e; ++j) {
                        Kpu->val[j] *= (*Kuu_dia)[Kpu->col[j]];
                    }
                }

                Kpp = backend::sum(
                        math::identity<value_type>(), *Kpp,
                       -math::identity<value_type>(), *backend::product(*Kpu, *Kup));
            }

            U = std::make_shared<USolver>(*Kuu, prm.usolver, bprm);
            P = std::make_shared<PSolver>(*Kpp, prm.psolver, bprm);
            M = backend_type::copy_vector(Kuu_dia, bprm);

            rhs_u = backend_type::create_vector(nu, bprm);
            rhs_p = backend_type::create_vector(np, bprm);

            u = backend_type::create_vector(nu, bprm);
            p = backend_type::create_vector(np, bprm);

            // Scatter/Gather matrices
            auto x2u = std::make_shared<build_matrix>();
            auto x2p = std::make_shared<build_matrix>();
            auto u2x = std::make_shared<build_matrix>();
            auto p2x = std::make_shared<build_matrix>();

            x2u->set_size(nu, n, true);
            x2p->set_size(np, n, true);
            u2x->set_size(n, nu, true);
            p2x->set_size(n, np, true);

            {
                ptrdiff_t x2u_head = 0, x2u_idx = 0;
                ptrdiff_t x2p_head = 0, x2p_idx = 0;
                ptrdiff_t u2x_head = 0, u2x_idx = 0;
                ptrdiff_t p2x_head = 0, p2x_idx = 0;

                for(size_t i = 0; i < n; ++i) {
                    if (prm.pmask[i]) {
                        x2p->ptr[++x2p_idx] = ++x2p_head;
                        ++p2x_head;
                    } else {
                        x2u->ptr[++x2u_idx] = ++x2u_head;
                        ++u2x_head;
                    }

                    p2x->ptr[++p2x_idx] = p2x_head;
                    u2x->ptr[++u2x_idx] = u2x_head;
                }
            }

            x2u->set_nonzeros();
            x2p->set_nonzeros();
            u2x->set_nonzeros();
            p2x->set_nonzeros();

            {
                ptrdiff_t x2u_head = 0;
                ptrdiff_t x2p_head = 0;
                ptrdiff_t u2x_head = 0;
                ptrdiff_t p2x_head = 0;

                for(size_t i = 0; i < n; ++i) {
                    ptrdiff_t j = idx[i];

                    if (prm.pmask[i]) {
                        x2p->col[x2p_head] = i;
                        x2p->val[x2p_head] = math::identity<value_type>();
                        ++x2p_head;

                        p2x->col[p2x_head] = j;
                        p2x->val[p2x_head] = math::identity<value_type>();
                        ++p2x_head;
                    } else {
                        x2u->col[x2u_head] = i;
                        x2u->val[x2u_head] = math::identity<value_type>();
                        ++x2u_head;

                        u2x->col[u2x_head] = j;
                        u2x->val[u2x_head] = math::identity<value_type>();
                        ++u2x_head;
                    }
                }
            }

            this->x2u = backend_type::copy_matrix(x2u, bprm);
            this->x2p = backend_type::copy_matrix(x2p, bprm);
            this->u2x = backend_type::copy_matrix(u2x, bprm);
            this->p2x = backend_type::copy_matrix(p2x, bprm);
        }

        friend std::ostream& operator<<(std::ostream &os, const simple &p) {
            os << "SIMPLE (two-stage preconditioner)" << std::endl;
            os << "  Unknowns: " << p.n << "(" << p.np << ")" << std::endl;
            os << "  Nonzeros: " << backend::nonzeros(p.system_matrix()) << std::endl;
            os << "  Memory:  " << human_readable_memory(p.bytes()) << std::endl;
            os << std::endl;
            os << "[ U ]\n" << *p.U << std::endl;
            os << "[ P ]\n" << *p.P << std::endl;

            return os;
        }

#if defined(AMGCL_DEBUG)
        template <typename I, typename E>
        static void report(const std::string &name, const std::tuple<I, E> &c) {
            std::cout << name << " (" << std::get<0>(c) << ", " << std::get<1>(c) << ")\n";
        }
#else
        template <typename I, typename E>
        static void report(const std::string&, const std::tuple<I, E>&) {
        }
#endif
};

} // namespace preconditioner
} // namespace amgcl

#endif
