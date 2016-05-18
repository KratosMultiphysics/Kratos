#ifndef AMGCL_PRECONDITIONER_SCHUR_COMPLEMENT_HPP
#define AMGCL_PRECONDITIONER_SCHUR_COMPLEMENT_HPP

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
 * \file   amgcl/preconditioner/schur_complement.hpp
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

/// Schur-complement pressure correction preconditioner
template <class USolver, class PSolver>
class schur_complement {
    BOOST_STATIC_ASSERT_MSG(
            (
             boost::is_same<
                 typename USolver::backend_type,
                 typename PSolver::backend_type
                 >::value
            ),
            "Backends for pressure and flow preconditioners should coinside!"
            );
    public:
        typedef typename USolver::backend_type backend_type;

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
        schur_complement(
                const Matrix &K,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
            : prm(prm), n(backend::rows(K)), np(0), nu(0)
        {
            init(boost::make_shared<build_matrix>(K), bprm);
        }

        schur_complement(
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
            (*U)(*rhs_u, *u);

            // rhs_p -= Kpu u
            backend::spmv(-1, *Kpu, *u, 1, *rhs_p);

            // S p = rhs_p
            backend::clear(*p);
            (*P)(*this, *rhs_p, *p);

            // rhs_u -= Kup p
            backend::spmv(-1, *Kup, *p, 1, *rhs_u);

            // Ai u = rhs_u
            backend::clear(*u);
            (*U)(*rhs_u, *u);

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
            backend::spmv( alpha, *Kpp, x, beta, y);

            backend::spmv(1, *Kup, x, 0, *tmp);
            backend::clear(*u);
            (*U)(*tmp, *u);
            backend::spmv(-alpha, *Kpu, *u, 1, y);
        }
    private:
        size_t n, np, nu;

        boost::shared_ptr<matrix> K, Kup, Kpu, Kpp, x2u, x2p, u2x, p2x;
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
            this->Kpp = backend_type::copy_matrix(Kpp, bprm);

            rhs_u = backend_type::create_vector(nu, bprm);
            rhs_p = backend_type::create_vector(np, bprm);

            u = backend_type::create_vector(nu, bprm);
            p = backend_type::create_vector(np, bprm);

            tmp = backend_type::create_vector(nu, bprm);

            // Scatter/Gather matrices
            boost::shared_ptr<build_matrix> X2U = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> X2P = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> U2X = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> P2X = boost::make_shared<build_matrix>();

            boost::tie(X2U->nrows, X2U->ncols) = boost::make_tuple(nu, n);
            boost::tie(X2P->nrows, X2P->ncols) = boost::make_tuple(np, n);
            boost::tie(U2X->nrows, U2X->ncols) = boost::make_tuple(n, nu);
            boost::tie(P2X->nrows, P2X->ncols) = boost::make_tuple(n, np);

            X2U->ptr.reserve(nu+1); X2U->ptr.push_back(0);
            X2P->ptr.reserve(np+1); X2P->ptr.push_back(0);
            U2X->ptr.reserve(n +1); U2X->ptr.push_back(0);
            P2X->ptr.reserve(n +1); P2X->ptr.push_back(0);

            X2U->col.reserve(nu);
            X2P->col.reserve(np);
            U2X->col.reserve(nu);
            P2X->col.reserve(np);

            X2U->val.resize(nu, 1.0);
            X2P->val.resize(np, 1.0);
            U2X->val.resize(nu, 1.0);
            P2X->val.resize(np, 1.0);

            for(size_t i = 0; i < n; ++i) {
                ptrdiff_t j = idx[i];

                if (prm.pmask[i]) {
                    X2P->col.push_back(i);
                    X2P->ptr.push_back(X2P->col.size());

                    P2X->col.push_back(j);
                } else {
                    X2U->col.push_back(i);
                    X2U->ptr.push_back(X2U->col.size());

                    U2X->col.push_back(j);
                }

                P2X->ptr.push_back(P2X->col.size());
                U2X->ptr.push_back(U2X->col.size());
            }

            x2u = backend_type::copy_matrix(X2U, bprm);
            x2p = backend_type::copy_matrix(X2P, bprm);
            u2x = backend_type::copy_matrix(U2X, bprm);
            p2x = backend_type::copy_matrix(P2X, bprm);
        }

        friend std::ostream& operator<<(std::ostream &os, const schur_complement &p) {
            os << "Schur complement (two-stage preconditioner)" << std::endl;
            os << "  unknowns: " << p.n << "(" << p.np << ")" << std::endl;
            os << "  nonzeros: " << backend::nonzeros(p.system_matrix()) << std::endl;

            return os;
        }
};

} // namespace preconditioner

namespace backend {

template <class US, class PS, class Alpha, class Beta, class Vec1, class Vec2>
struct spmv_impl< Alpha, preconditioner::schur_complement<US, PS>, Vec1, Beta, Vec2>
{
    static void apply(Alpha alpha, const preconditioner::schur_complement<US, PS> &A, const Vec1 &x, Beta beta, Vec2 &y)
    {
        A.spmv(alpha, x, beta, y);
    }
};

template <class US, class PS, class Vec1, class Vec2, class Vec3>
struct residual_impl< preconditioner::schur_complement<US, PS>, Vec1, Vec2, Vec3>
{
    static void apply(const Vec1 &rhs, const preconditioner::schur_complement<US, PS> &A, const Vec2 &x, Vec3 &r)
    {
        backend::copy(rhs, r);
        A.spmv(-1, x, 1, r);
    }
};

} // namespace backend
} // namespace amgcl

#endif
