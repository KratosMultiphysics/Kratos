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
            : prm(prm), n(backend::rows(K)), np(0), nu(0),
              _K(backend_type::copy_matrix(boost::make_shared<build_matrix>(K), bprm)), idx(n)
        {
            typedef typename backend::row_iterator<Matrix>::type row_iterator;

            // Extract matrix subblocks.
            boost::shared_ptr<build_matrix> A  = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> B  = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> BT = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> C  = boost::make_shared<build_matrix>();

            for(size_t i = 0; i < n; ++i)
                idx[i] = (prm.pmask[i] ? np++ : nu++);

            boost::tie(A->nrows,  A->ncols ) = boost::make_tuple(nu, nu);
            boost::tie(BT->nrows, BT->ncols) = boost::make_tuple(nu, np);
            boost::tie(B->nrows,  B->ncols ) = boost::make_tuple(np, nu);
            boost::tie(C->nrows,  C->ncols ) = boost::make_tuple(np, np);

            A->ptr.resize(nu + 1, 0);
            BT->ptr.resize(nu + 1, 0);
            B->ptr.resize(np + 1, 0);
            C->ptr.resize(np + 1, 0);

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                ptrdiff_t ci = idx[i];
                char      pi = prm.pmask[i];
                for(row_iterator k = backend::row_begin(K, i); k; ++k) {
                    char pj = prm.pmask[k.col()];

                    if (pi) {
                        if (pj) {
                            ++C->ptr[ci+1];
                        } else {
                            ++B->ptr[ci+1];
                        }
                    } else {
                        if (pj) {
                            ++BT->ptr[ci+1];
                        } else {
                            ++A->ptr[ci+1];
                        }
                    }
                }
            }

            boost::partial_sum(A->ptr,  A->ptr.begin());
            boost::partial_sum(BT->ptr, BT->ptr.begin());
            boost::partial_sum(B->ptr,  B->ptr.begin());
            boost::partial_sum(C->ptr,  C->ptr.begin());

            A->col.resize(A->ptr.back());
            A->val.resize(A->ptr.back());

            BT->col.resize(BT->ptr.back());
            BT->val.resize(BT->ptr.back());

            B->col.resize(B->ptr.back());
            B->val.resize(B->ptr.back());

            C->col.resize(C->ptr.back());
            C->val.resize(C->ptr.back());

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                ptrdiff_t ci = idx[i];
                char      pi = prm.pmask[i];

                ptrdiff_t A_head = 0, BT_head = 0, B_head = 0, C_head = 0;
                if(pi) {
                    B_head = B->ptr[ci];
                    C_head = C->ptr[ci];
                } else {
                    A_head  = A->ptr[ci];
                    BT_head = BT->ptr[ci];
                }

                for(row_iterator k = backend::row_begin(K, i); k; ++k) {
                    ptrdiff_t  j = k.col();
                    value_type v = k.value();
                    ptrdiff_t cj = idx[j];
                    char      pj = prm.pmask[j];

                    if (pi) {
                        if (pj) {
                            C->col[C_head] = cj;
                            C->val[C_head] = -v;
                            ++C_head;
                        } else {
                            B->col[B_head] = cj;
                            B->val[B_head] = v;
                            ++B_head;
                        }
                    } else {
                        if (pj) {
                            BT->col[BT_head] = cj;
                            BT->val[BT_head] = v;
                            ++BT_head;
                        } else {
                            A->col[A_head] = cj;
                            A->val[A_head] = v;
                            ++A_head;
                        }
                    }
                }
            }

            U = boost::make_shared<USolver>(*A, prm.usolver, bprm);
            P = boost::make_shared<PSolver>(*C, prm.psolver, bprm);

            _B  = backend_type::copy_matrix(B,  bprm);
            _BT = backend_type::copy_matrix(BT, bprm);
            _C  = backend_type::copy_matrix(C,  bprm);

            rhs_u = backend_type::create_vector(nu, bprm);
            rhs_p = backend_type::create_vector(np, bprm);

            u = backend_type::create_vector(nu, bprm);
            p = backend_type::create_vector(np, bprm);

            tmp1 = backend_type::create_vector(nu, bprm);
            tmp2 = backend_type::create_vector(nu, bprm);

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
            backend::spmv( 1, *x2u, rhs, 0, *rhs_u);
            backend::spmv(-1, *x2p, rhs, 0, *rhs_p);

            size_t iters;
            double error;

            // Ai u = rhs_u
            backend::clear(*u);
            boost::tie(iters, error) = (*U)(*rhs_u, *u);
            std::cout << "rhs(p): (" << iters << ", " << error << ")" << std::endl;

            // rhs_p += B u
            backend::spmv(1, *_B, *u, 1, *rhs_p);

            // S p = rhs_p
            backend::clear(*p);
            boost::tie(iters, error) = (*P)(*this, *rhs_p, *p);
            std::cout << "S(p)=r: (" << iters << ", " << error << ")" << std::endl;

            // rhs_u -= BT p
            backend::spmv(-1, *_BT, *p, 1, *rhs_u);

            // Ai u = rhs_u
            backend::clear(*u);
            boost::tie(iters, error) = (*U)(*rhs_u, *u);
            std::cout << "A(u)=f: (" << iters << ", " << error << ")" << std::endl;

            backend::clear(x);
            backend::spmv(1, *u2x, *u, 1, x);
            backend::spmv(1, *p2x, *p, 1, x);
        }

        const matrix& system_matrix() const {
            return *_K;
        }

        template <class Alpha, class Vec1, class Beta, class Vec2>
        void spmv(Alpha alpha, const Vec1 &x, Beta beta, Vec2 &y) const {
            backend::spmv(1, *_BT, x, 0, *tmp1);
            backend::clear(*tmp2);
            (*U)(*tmp1, *tmp2);
            backend::spmv(alpha, *_B, *tmp2, beta, y);
            backend::spmv(alpha, *_C, x, 1, y);
        }
    private:
        size_t n, np, nu;

        boost::shared_ptr<matrix> _K, _B, _BT, _C, x2u, x2p, u2x, p2x;
        boost::shared_ptr<vector> rhs_u, rhs_p, u, p, tmp1, tmp2;

        std::vector<ptrdiff_t> idx;

        boost::shared_ptr<USolver> U;
        boost::shared_ptr<PSolver> P;
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
