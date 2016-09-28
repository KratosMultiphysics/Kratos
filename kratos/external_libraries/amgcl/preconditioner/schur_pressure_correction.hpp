#ifndef AMGCL_PRECONDITIONER_SCHUR_PRESSURE_CORRECTION_HPP
#define AMGCL_PRECONDITIONER_SCHUR_PRESSURE_CORRECTION_HPP

/*
The MIT License

Copyright (c) 2012-2016 Denis Demidov <dennis.demidov@gmail.com>
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

#include <boost/static_assert.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace preconditioner {

namespace detail {

// Same backends are always compatible
template <class B1, class B2>
struct compatible_backends
    : boost::is_same<B1, B2>::type {};

// Builtin backend allows mixing backends of different value types,
// so that scalar and non-scalar backends may coexist.
template <class V1, class V2>
struct compatible_backends< backend::builtin<V1>, backend::builtin<V2> >
    : boost::true_type {};

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
    typename boost::disable_if<typename boost::is_same<V1, V2>::type>::type >
{
    typedef
        typename boost::conditional<
            (math::static_rows<V1>::value <= math::static_rows<V2>::value),
            backend::builtin<V1>, backend::builtin<V2>
            >::type
        type;
};

} // namespace detail

/// Schur-complement pressure correction preconditioner
template <class USolver, class PSolver>
class schur_pressure_correction {
    BOOST_STATIC_ASSERT_MSG(
            (
             detail::compatible_backends<
                 typename USolver::backend_type,
                 typename PSolver::backend_type
                 >::value
            ),
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

            params() {}

            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_CHILD(p, usolver),
                  AMGCL_PARAMS_IMPORT_CHILD(p, psolver)
            {
                void *pm = 0;
                size_t n = 0;

                pm = p.get("pmask",     pm);
                n  = p.get("pmask_size", n);

                precondition(pm,
                        "Error in schur_complement parameters: "
                        "pmask is not set");

                precondition(n > 0,
                        "Error in schur_complement parameters: "
                        "pmask is set, but pmask_size is not"
                        );

                pmask.assign(static_cast<char*>(pm), static_cast<char*>(pm) + n);

                AMGCL_PARAMS_CHECK(p, (usolver)(psolver)(pmask)(pmask_size));
            }

            void get(boost::property_tree::ptree &p, const std::string &path = "") const
            {
                AMGCL_PARAMS_EXPORT_CHILD(p, path, usolver);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, psolver);
            }
        } prm;

        template <class Matrix>
        schur_pressure_correction(
                const Matrix &K,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
            : prm(prm), n(backend::rows(K)), np(0), nu(0)
        {
            init(boost::make_shared<build_matrix>(K), bprm);
        }

        schur_pressure_correction(
                boost::shared_ptr<build_matrix> K,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
            : prm(prm), n(backend::rows(*K)), np(0), nu(0)
        {
            init(K, bprm);
        }

        template <class Vec1, class Vec2>
        void apply(
                const Vec1 &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2       &x
#else
                Vec2       &&x
#endif
                ) const
        {
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

        const matrix& system_matrix() const {
            return *K;
        }

        template <class Alpha, class Vec1, class Beta, class Vec2>
        void spmv(Alpha alpha, const Vec1 &x, Beta beta, Vec2 &y) const {
            // y = beta y + alpha S x, where S = Kpp - Kpu Kuu^-1 Kup
            backend::spmv( alpha, P->system_matrix(), x, beta, y);

            backend::spmv(1, *Kup, x, 0, *tmp);
            backend::clear(*u);
            (*U)(*tmp, *u);
            backend::spmv(-alpha, *Kpu, *u, 1, y);
        }
    private:
        size_t n, np, nu;

        boost::shared_ptr<matrix> K, Kup, Kpu, x2u, x2p, u2x, p2x;
        boost::shared_ptr<vector> rhs_u, rhs_p, u, p, tmp;

        boost::shared_ptr<USolver> U;
        boost::shared_ptr<PSolver> P;

        void init(const boost::shared_ptr<build_matrix> &K, const backend_params &bprm)
        {
            typedef typename backend::row_iterator<build_matrix>::type row_iterator;

            this->K = backend_type::copy_matrix(K, bprm);

            // Extract matrix subblocks.
            boost::shared_ptr<build_matrix> Kuu = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> Kpu = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> Kup = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> Kpp = boost::make_shared<build_matrix>();

            std::vector<ptrdiff_t> idx(n);

            for(size_t i = 0; i < n; ++i)
                idx[i] = (prm.pmask[i] ? np++ : nu++);

            boost::tie(Kuu->nrows, Kuu->ncols) = boost::make_tuple(nu, nu);
            boost::tie(Kup->nrows, Kup->ncols) = boost::make_tuple(nu, np);
            boost::tie(Kpu->nrows, Kpu->ncols) = boost::make_tuple(np, nu);
            boost::tie(Kpp->nrows, Kpp->ncols) = boost::make_tuple(np, np);

            Kuu->ptr.resize(nu + 1, 0);
            Kup->ptr.resize(nu + 1, 0);
            Kpu->ptr.resize(np + 1, 0);
            Kpp->ptr.resize(np + 1, 0);

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                ptrdiff_t ci = idx[i];
                char      pi = prm.pmask[i];
                for(row_iterator k = backend::row_begin(*K, i); k; ++k) {
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

            boost::partial_sum(Kuu->ptr, Kuu->ptr.begin());
            boost::partial_sum(Kup->ptr, Kup->ptr.begin());
            boost::partial_sum(Kpu->ptr, Kpu->ptr.begin());
            boost::partial_sum(Kpp->ptr, Kpp->ptr.begin());

            Kuu->col.resize(Kuu->ptr.back());
            Kuu->val.resize(Kuu->ptr.back());

            Kup->col.resize(Kup->ptr.back());
            Kup->val.resize(Kup->ptr.back());

            Kpu->col.resize(Kpu->ptr.back());
            Kpu->val.resize(Kpu->ptr.back());

            Kpp->col.resize(Kpp->ptr.back());
            Kpp->val.resize(Kpp->ptr.back());

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

                for(row_iterator k = backend::row_begin(*K, i); k; ++k) {
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

            U = boost::make_shared<USolver>(*Kuu, prm.usolver, bprm);
            P = boost::make_shared<PSolver>(*Kpp, prm.psolver, bprm);

            this->Kup = backend_type::copy_matrix(Kup, bprm);
            this->Kpu = backend_type::copy_matrix(Kpu, bprm);

            rhs_u = backend_type::create_vector(nu, bprm);
            rhs_p = backend_type::create_vector(np, bprm);

            u = backend_type::create_vector(nu, bprm);
            p = backend_type::create_vector(np, bprm);

            tmp = backend_type::create_vector(nu, bprm);

            // Scatter/Gather matrices
            boost::shared_ptr<build_matrix> x2u = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> x2p = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> u2x = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> p2x = boost::make_shared<build_matrix>();

            boost::tie(x2u->nrows, x2u->ncols) = boost::make_tuple(nu, n);
            boost::tie(x2p->nrows, x2p->ncols) = boost::make_tuple(np, n);
            boost::tie(u2x->nrows, u2x->ncols) = boost::make_tuple(n, nu);
            boost::tie(p2x->nrows, p2x->ncols) = boost::make_tuple(n, np);

            x2u->ptr.reserve(nu+1); x2u->ptr.push_back(0);
            x2p->ptr.reserve(np+1); x2p->ptr.push_back(0);
            u2x->ptr.reserve(n +1); u2x->ptr.push_back(0);
            p2x->ptr.reserve(n +1); p2x->ptr.push_back(0);

            x2u->col.reserve(nu);
            x2p->col.reserve(np);
            u2x->col.reserve(nu);
            p2x->col.reserve(np);

            x2u->val.resize(nu, math::identity<value_type>());
            x2p->val.resize(np, math::identity<value_type>());
            u2x->val.resize(nu, math::identity<value_type>());
            p2x->val.resize(np, math::identity<value_type>());

            for(size_t i = 0; i < n; ++i) {
                ptrdiff_t j = idx[i];

                if (prm.pmask[i]) {
                    x2p->col.push_back(i);
                    x2p->ptr.push_back(x2p->col.size());

                    p2x->col.push_back(j);
                } else {
                    x2u->col.push_back(i);
                    x2u->ptr.push_back(x2u->col.size());

                    u2x->col.push_back(j);
                }

                p2x->ptr.push_back(p2x->col.size());
                u2x->ptr.push_back(u2x->col.size());
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

        template <typename I, typename E>
        static void report(const std::string &name, const boost::tuple<I, E> &c) {
#if defined(AMGCL_DEBUG) || !defined(NDEBUG)
            std::cout << name << " (" << boost::get<0>(c) << ", " << boost::get<1>(c) << ")\n";
#endif
        }
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
