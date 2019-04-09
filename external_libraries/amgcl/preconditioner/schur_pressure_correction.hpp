#ifndef AMGCL_PRECONDITIONER_SCHUR_PRESSURE_CORRECTION_HPP
#define AMGCL_PRECONDITIONER_SCHUR_PRESSURE_CORRECTION_HPP

/*
The MIT License

Copyright (c) 2012-2019 Denis Demidov <dennis.demidov@gmail.com>
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
 * \file   amgcl/preconditioner/schur_pressure_correction.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Schur-complement pressure correction preconditioning scheme.
 */

#include <vector>

#include <memory>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace preconditioner {

namespace detail {

// Same backends are always compatible
template <class B1, class B2>
struct compatible_backends
    : std::is_same<B1, B2>::type {};

// Builtin backend allows mixing backends of different value types,
// so that scalar and non-scalar backends may coexist.
template <class V1, class V2>
struct compatible_backends< backend::builtin<V1>, backend::builtin<V2> >
    : std::true_type {};

// Backend for schur complement preconditioner is selected as the one with
// lower dimensionality of its value_type.

template <class B1, class B2, class Enable = void>
struct common_backend;

template <class B>
struct common_backend<B, B> {
    typedef B type;
};

template <class V1, class V2>
struct common_backend< backend::builtin<V1>, backend::builtin<V2>,
    typename std::enable_if<!std::is_same<V1, V2>::value>::type >
{
    typedef
        typename std::conditional<
            (math::static_rows<V1>::value <= math::static_rows<V2>::value),
            backend::builtin<V1>, backend::builtin<V2>
            >::type
        type;
};

} // namespace detail

/// Schur-complement pressure correction preconditioner
template <class USolver, class PSolver>
class schur_pressure_correction {
    static_assert(
            detail::compatible_backends<
                typename USolver::backend_type,
                typename PSolver::backend_type
                >::value,
            "Backends for pressure and flow preconditioners should coincide!"
            );
    public:
        typedef
            typename detail::common_backend<
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

            // Approximate Kuu^-1 with inverted diagonal of Kuu during
            // construction of matrix-less Schur complement.
            // When false, USolver is used instead.
            bool approx_schur;

            params() : approx_schur(true) {}

#ifndef AMGCL_NO_BOOST
            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_CHILD(p, usolver),
                  AMGCL_PARAMS_IMPORT_CHILD(p, psolver),
                  AMGCL_PARAMS_IMPORT_VALUE(p, approx_schur)
            {
                size_t n = 0;

                n = p.get("pmask_size", n);

                precondition(n > 0,
                        "Error in schur_complement parameters: "
                        "pmask_size is not set");

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
                            "Error in schur_complement parameters: "
                            "neither pmask_pattern, nor pmask is set"
                            );
                }

                check_params(p, {"usolver", "psolver", "approx_schur", "pmask_size"},
                        {"pmask", "pmask_pattern"});
            }

            void get(boost::property_tree::ptree &p, const std::string &path = "") const
            {
                AMGCL_PARAMS_EXPORT_CHILD(p, path, usolver);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, psolver);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, approx_schur);
            }
#endif
        } prm;

        template <class Matrix>
        schur_pressure_correction(
                const Matrix &K,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
            : prm(prm), n(backend::rows(K)), np(0), nu(0)
        {
            init(std::make_shared<build_matrix>(K), bprm);
        }

        schur_pressure_correction(
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

            // Ai u = rhs_u
            backend::clear(*u);
            report("U1", (*U)(*rhs_u, *u));

            // rhs_p -= Kpu u
            backend::spmv(-1, *Kpu, *u, 1, *rhs_p);

            // S p = rhs_p
            backend::clear(*p);
            report("P1", (*P)(*this, *rhs_p, *p));

            // rhs_u -= Kup p
            backend::spmv(-1, *Kup, *p, 1, *rhs_u);

            // Ai u = rhs_u
            backend::clear(*u);
            report("U2", (*U)(*rhs_u, *u));

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

        template <class Alpha, class Vec1, class Beta, class Vec2>
        void spmv(Alpha alpha, const Vec1 &x, Beta beta, Vec2 &y) const {
            // y = beta y + alpha S x, where S = Kpp - Kpu Kuu^-1 Kup
            backend::spmv( alpha, P->system_matrix(), x, beta, y);

            backend::spmv(1, *Kup, x, 0, *tmp);

            if (prm.approx_schur) {
                backend::vmul(1, *M, *tmp, 0, *u);
            } else {
                backend::clear(*u);
                (*U)(*tmp, *u);
            }

            backend::spmv(-alpha, *Kpu, *u, 1, y);
        }
    private:
        size_t n, np, nu;

        std::shared_ptr<matrix> K, Kup, Kpu, x2u, x2p, u2x, p2x;
        std::shared_ptr<vector> rhs_u, rhs_p, u, p, tmp;
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

            U = std::make_shared<USolver>(*Kuu, prm.usolver, bprm);
            P = std::make_shared<PSolver>(*Kpp, prm.psolver, bprm);

            this->Kup = backend_type::copy_matrix(Kup, bprm);
            this->Kpu = backend_type::copy_matrix(Kpu, bprm);

            rhs_u = backend_type::create_vector(nu, bprm);
            rhs_p = backend_type::create_vector(np, bprm);

            u = backend_type::create_vector(nu, bprm);
            p = backend_type::create_vector(np, bprm);

            tmp = backend_type::create_vector(nu, bprm);

            if (prm.approx_schur)
                M = backend_type::copy_vector(diagonal(*Kuu, /*invert = */true), bprm);

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

        friend std::ostream& operator<<(std::ostream &os, const schur_pressure_correction &p) {
            os << "Schur complement (two-stage preconditioner)" << std::endl;
            os << "  unknowns: " << p.n << "(" << p.np << ")" << std::endl;
            os << "  nonzeros: " << backend::nonzeros(p.system_matrix()) << std::endl;

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

namespace backend {

template <class US, class PS, class Alpha, class Beta, class Vec1, class Vec2>
struct spmv_impl< Alpha, preconditioner::schur_pressure_correction<US, PS>, Vec1, Beta, Vec2>
{
    static void apply(Alpha alpha, const preconditioner::schur_pressure_correction<US, PS> &A, const Vec1 &x, Beta beta, Vec2 &y)
    {
        A.spmv(alpha, x, beta, y);
    }
};

template <class US, class PS, class Vec1, class Vec2, class Vec3>
struct residual_impl< preconditioner::schur_pressure_correction<US, PS>, Vec1, Vec2, Vec3>
{
    static void apply(const Vec1 &rhs, const preconditioner::schur_pressure_correction<US, PS> &A, const Vec2 &x, Vec3 &r)
    {
        backend::copy(rhs, r);
        A.spmv(-1, x, 1, r);
    }
};

} // namespace backend
} // namespace amgcl

#endif
