#ifndef AMGCL_MPI_SCHUR_PRESSURE_CORRECTION_HPP
#define AMGCL_MPI_SCHUR_PRESSURE_CORRECTION_HPP

/*
The MIT License

Copyright (c) 2012-2017 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/mpi/schur_pressure_correction.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Distributed Schur complement pressure correction preconditioner.
 */

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>
#include <amgcl/mpi/inner_product.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>

namespace amgcl {
namespace mpi {

template <class USolver, class PSolver>
class schur_pressure_correction {
    BOOST_STATIC_ASSERT_MSG(
            (
             boost::is_same<
                 typename USolver::backend_type,
                 typename PSolver::backend_type
                 >::value
            ),
            "Backends for pressure and flow preconditioners should coincide!"
            );

    public:
        typedef typename USolver::backend_type backend_type;

        typedef typename backend_type::value_type value_type;
        typedef typename backend_type::matrix     bmatrix;
        typedef typename backend_type::vector     vector;
        typedef typename backend_type::params     backend_params;

        typedef distributed_matrix<backend_type> matrix;

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

                amgcl::precondition(pm,
                        "Error in schur_complement parameters: "
                        "pmask is not set");

                amgcl::precondition(n > 0,
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
                MPI_Comm mpi_comm,
                const Matrix &K,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
            : prm(prm), comm(mpi_comm)
        {
            typedef typename backend::row_iterator<Matrix>::type row_iterator;
            using boost::tie;
            using boost::make_tuple;
            using boost::shared_ptr;
            using boost::make_shared;

            // Get sizes of each domain in comm.
            ptrdiff_t n = backend::rows(K);
            std::vector<ptrdiff_t> domain = mpi::exclusive_sum(comm, n);

            ptrdiff_t loc_beg = domain[comm.rank];
            ptrdiff_t loc_end = domain[comm.rank + 1];

            // Count pressure and flow variables.
            std::vector<ptrdiff_t> idx(n);
            ptrdiff_t np = 0, nu = 0;

            for(ptrdiff_t i = 0; i < n; ++i)
                idx[i] = (prm.pmask[i] ? np++ : nu++);

            // Split the matrix into local and remote parts.
            shared_ptr<build_matrix> K_loc = make_shared<build_matrix>();
            shared_ptr<build_matrix> K_rem = make_shared<build_matrix>();

            K_loc->set_size(n, n, true);
            K_rem->set_size(n, 0, true); // number of columns is unknown at this point

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n; ++i) {
                for(row_iterator a = backend::row_begin(K, i); a; ++a) {
                    ptrdiff_t c = a.col();

                    if (loc_beg <= c && c < loc_end)
                        ++K_loc->ptr[i+1];
                    else
                        ++K_rem->ptr[i+1];
                }
            }

            std::partial_sum(K_loc->ptr, K_loc->ptr + n + 1, K_loc->ptr);
            std::partial_sum(K_rem->ptr, K_rem->ptr + n + 1, K_rem->ptr);

            K_loc->set_nonzeros(K_loc->ptr[n]);
            K_rem->set_nonzeros(K_rem->ptr[n]);

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n; ++i) {
                ptrdiff_t loc_head = K_loc->ptr[i];
                ptrdiff_t rem_head = K_rem->ptr[i];

                for(row_iterator a = backend::row_begin(K, i); a; ++a) {
                    ptrdiff_t  c = a.col();
                    value_type v = a.value();

                    if (loc_beg <= c && c < loc_end) {
                        K_loc->col[loc_head] = c - loc_beg;
                        K_loc->val[loc_head] = v;
                        ++loc_head;
                    } else {
                        K_rem->col[rem_head] = c;
                        K_rem->val[rem_head] = v;
                        ++rem_head;
                    }
                }
            }

            // Analyze communication pattern for the system matrix
            C = boost::make_shared<CommPattern>(comm, n, K_rem->nnz, K_rem->col, bprm);
            K_rem->ncols = C->renumber(K_rem->nnz, K_rem->col);

            this->K_loc = backend_type::copy_matrix(K_loc, bprm);
            this->K_rem = backend_type::copy_matrix(K_rem, bprm);
            this->K     = make_shared<matrix>(*C, *this->K_loc, *this->K_rem);

            // We know what points each of our neighbors needs from us;
            // and we know if those points are pressure or flow.
            // We can immediately provide them with our renumbering scheme.
            std::vector<ptrdiff_t> pdomain = mpi::exclusive_sum(comm, np);
            std::vector<ptrdiff_t> udomain = mpi::exclusive_sum(comm, nu);
            ptrdiff_t p_beg = pdomain[comm.rank];
            ptrdiff_t u_beg = udomain[comm.rank];

            ptrdiff_t nsend = C->send.ptr.back(), nrecv = C->recv.ptr.back();
            std::vector<char>      smask(nsend), rmask(nrecv);
            std::vector<ptrdiff_t> s_idx(nsend), r_idx(nrecv);

            for(ptrdiff_t i = 0; i < nsend; ++i) {
                ptrdiff_t c = C->send.col[i];
                smask[i] = prm.pmask[c];
                s_idx[i] = idx[c] + (smask[i] ? p_beg : u_beg);
            }

            C->exchange(&smask[0], &rmask[0]);
            C->exchange(&s_idx[0], &r_idx[0]);

            // Fill the subblocks of the system matrix.
            // Kpp and Kpp have to be constructed as whole strips, and
            // Kup and Kpu may be split to local/remote parts immediately.
            // K_rem->col may be used as direct indices into rmask and r_idx.
            shared_ptr<build_matrix> Kpp = make_shared<build_matrix>();
            shared_ptr<build_matrix> Kuu = make_shared<build_matrix>();

            shared_ptr<build_matrix> Kpu_loc = make_shared<build_matrix>();
            shared_ptr<build_matrix> Kpu_rem = make_shared<build_matrix>();
            shared_ptr<build_matrix> Kup_loc = make_shared<build_matrix>();
            shared_ptr<build_matrix> Kup_rem = make_shared<build_matrix>();

            Kpp->set_size(np, pdomain.back(), true);
            Kuu->set_size(nu, udomain.back(), true);

            Kpu_loc->set_size(np, nu, true);
            Kup_loc->set_size(nu, np, true);

            Kpu_rem->set_size(np, 0, true);
            Kup_rem->set_size(np, 0, true);

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n; ++i) {
                typedef typename backend::row_iterator<build_matrix>::type row_iterator;

                ptrdiff_t ci = idx[i];
                char      pi = prm.pmask[i];

                for(row_iterator a = row_begin(*K_loc, i); a; ++a) {
                    char pj = prm.pmask[a.col()];

                    if (pi) {
                        if (pj) {
                            ++Kpp->ptr[ci + 1];
                        } else {
                            ++Kpu_loc->ptr[ci + 1];
                        }
                    } else {
                        if (pj) {
                            ++Kup_loc->ptr[ci + 1];
                        } else {
                            ++Kuu->ptr[ci + 1];
                        }
                    }
                }

                for(row_iterator a = row_begin(*K_rem, i); a; ++a) {
                    char pj = rmask[a.col()];

                    if (pi) {
                        if (pj) {
                            ++Kpp->ptr[ci + 1];
                        } else {
                            ++Kpu_rem->ptr[ci + 1];
                        }
                    } else {
                        if (pj) {
                            ++Kup_rem->ptr[ci + 1];
                        } else {
                            ++Kuu->ptr[ci + 1];
                        }
                    }
                }
            }

            std::partial_sum(Kpp->ptr, Kpp->ptr + np + 1, Kpp->ptr);
            std::partial_sum(Kuu->ptr, Kuu->ptr + nu + 1, Kuu->ptr);

            std::partial_sum(Kpu_loc->ptr, Kpu_loc->ptr + np + 1, Kpu_loc->ptr);
            std::partial_sum(Kpu_rem->ptr, Kpu_rem->ptr + np + 1, Kpu_rem->ptr);

            std::partial_sum(Kup_loc->ptr, Kup_loc->ptr + nu + 1, Kup_loc->ptr);
            std::partial_sum(Kup_rem->ptr, Kup_rem->ptr + nu + 1, Kup_rem->ptr);

            Kpp->set_nonzeros(Kpp->ptr[np]);
            Kuu->set_nonzeros(Kuu->ptr[nu]);

            Kpu_loc->set_nonzeros(Kpu_loc->ptr[np]);
            Kpu_rem->set_nonzeros(Kpu_rem->ptr[np]);

            Kup_loc->set_nonzeros(Kup_loc->ptr[nu]);
            Kup_rem->set_nonzeros(Kup_rem->ptr[nu]);

            // Fill subblocks of the system matrix.
            // Kpp and Kuu will be fed to the solvers constructors, so the
            // columns should have global numeration.
            // Kpu and Kup are only needed as distributed matrices, so their
            // local and remote parts are constructed explicitly.
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n; ++i) {
                typedef typename backend::row_iterator<build_matrix>::type row_iterator;

                ptrdiff_t ci = idx[i];
                char      pi = prm.pmask[i];

                ptrdiff_t pp_head = 0, uu_head = 0;
                ptrdiff_t pu_loc_head = 0, pu_rem_head = 0;
                ptrdiff_t up_loc_head = 0, up_rem_head = 0;

                if(pi) {
                    pp_head = Kpp->ptr[ci];
                    pu_loc_head = Kpu_loc->ptr[ci];
                    pu_rem_head = Kpu_rem->ptr[ci];
                } else {
                    uu_head = Kuu->ptr[ci];
                    up_loc_head = Kup_loc->ptr[ci];
                    up_rem_head = Kup_rem->ptr[ci];
                }

                for(row_iterator a = row_begin(*K_loc, i); a; ++a) {
                    ptrdiff_t  j = a.col();
                    value_type v = a.value();
                    char      pj = prm.pmask[j];
                    ptrdiff_t cj = idx[j];
                    ptrdiff_t loc_beg = pj ? p_beg : u_beg;

                    if (pi) {
                        if (pj) {
                            Kpp->col[pp_head] = cj + loc_beg;
                            Kpp->val[pp_head] = v;
                            ++pp_head;
                        } else {
                            Kpu_loc->col[pu_loc_head] = cj;
                            Kpu_loc->val[pu_loc_head] = v;
                            ++pu_loc_head;
                        }
                    } else {
                        if (pj) {
                            Kup_loc->col[up_loc_head] = cj;
                            Kup_loc->val[up_loc_head] = v;
                            ++up_loc_head;
                        } else {
                            Kuu->col[uu_head] = cj + loc_beg;
                            Kuu->val[uu_head] = v;
                            ++uu_head;
                        }
                    }
                }

                for(row_iterator a = row_begin(*K_rem, i); a; ++a) {
                    ptrdiff_t  j = a.col();
                    value_type v = a.value();
                    char      pj = rmask[j];
                    ptrdiff_t cj = r_idx[j];

                    if (pi) {
                        if (pj) {
                            Kpp->col[pp_head] = cj;
                            Kpp->val[pp_head] = v;
                            ++pp_head;
                        } else {
                            Kpu_rem->col[pu_rem_head] = cj;
                            Kpu_rem->val[pu_rem_head] = v;
                            ++pu_rem_head;
                        }
                    } else {
                        if (pj) {
                            Kup_rem->col[up_rem_head] = cj;
                            Kup_rem->val[up_rem_head] = v;
                            ++up_rem_head;
                        } else {
                            Kuu->col[uu_head] = cj;
                            Kuu->val[uu_head] = v;
                            ++uu_head;
                        }
                    }
                }
            }

            Cpu = boost::make_shared<CommPattern>(comm, nu, Kpu_rem->nnz, Kpu_rem->col, bprm);
            Kpu_rem->ncols = Cpu->renumber(Kpu_rem->nnz, Kpu_rem->col);

            this->Kpu_loc = backend_type::copy_matrix(Kpu_loc, bprm);
            this->Kpu_rem = backend_type::copy_matrix(Kpu_rem, bprm);

            Kpu = make_shared<matrix>(*Cpu, *this->Kpu_loc, *this->Kpu_rem);


            Cup = boost::make_shared<CommPattern>(comm, np, Kup_rem->nnz, Kup_rem->col, bprm);
            Kup_rem->ncols = Cup->renumber(Kup_rem->nnz, Kup_rem->col);

            this->Kup_loc = backend_type::copy_matrix(Kup_loc, bprm);
            this->Kup_rem = backend_type::copy_matrix(Kup_rem, bprm);

            Kup = make_shared<matrix>(*Cup, *this->Kup_loc, *this->Kup_rem);

            U = make_shared<USolver>(mpi_comm, *Kuu, prm.usolver, bprm);
            P = make_shared<PSolver>(mpi_comm, *Kpp, prm.psolver, bprm);

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

            x2u->set_size(nu, n, true);
            x2p->set_size(np, n, true);
            u2x->set_size(n, nu, true);
            p2x->set_size(n, np, true);

            {
                ptrdiff_t x2u_head = 0, x2u_idx = 0;
                ptrdiff_t x2p_head = 0, x2p_idx = 0;
                ptrdiff_t u2x_head = 0, u2x_idx = 0;
                ptrdiff_t p2x_head = 0, p2x_idx = 0;

                for(ptrdiff_t i = 0; i < n; ++i) {
                    ptrdiff_t j = idx[i];

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

                for(ptrdiff_t i = 0; i < n; ++i) {
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
                    }
                }
            }

            this->x2u = backend_type::copy_matrix(x2u, bprm);
            this->x2p = backend_type::copy_matrix(x2p, bprm);
            this->u2x = backend_type::copy_matrix(u2x, bprm);
            this->p2x = backend_type::copy_matrix(p2x, bprm);
        }

        const matrix& system_matrix() const {
            return *K;
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
            TIC("split variables");
            backend::spmv(1, *x2u, rhs, 0, *rhs_u);
            backend::spmv(1, *x2p, rhs, 0, *rhs_p);
            TOC("split variables");

            // Ai u = rhs_u
            TIC("solve U");
            backend::clear(*u);
            report("U1", (*U)(*rhs_u, *u));
            TOC("solve U");

            // rhs_p -= Kpu u
            TIC("solve P");
            backend::spmv(-1, *Kpu, *u, 1, *rhs_p);

            // S p = rhs_p
            backend::clear(*p);
            report("P", (*P)(*this, *rhs_p, *p));
            TOC("solve P");

            // rhs_u -= Kup p
            TIC("Update U");
            backend::spmv(-1, *Kup, *p, 1, *rhs_u);

            // Ai u = rhs_u
            backend::clear(*u);
            report("U2", (*U)(*rhs_u, *u));
            TOC("Update U");

            TIC("merge variables");
            backend::clear(x);
            backend::spmv(1, *u2x, *u, 1, x);
            backend::spmv(1, *p2x, *p, 1, x);
            TOC("merge variables");
        }

        template <class Alpha, class Vec1, class Beta, class Vec2>
        void spmv(Alpha alpha, const Vec1 &x, Beta beta, Vec2 &y) const {
            // y = beta y + alpha S x, where S = Kpp - Kpu Kuu^-1 Kup
            TIC("matrix-free spmv");
            backend::spmv(alpha, P->system_matrix(), x, beta, y);

            backend::spmv(1, *Kup, x, 0, *tmp);
            backend::clear(*u);
            (*U)(*tmp, *u);
            backend::spmv(-alpha, *Kpu, *u, 1, y);
            TOC("matrix-free spmv");
        }
    private:
        typedef comm_pattern<backend_type> CommPattern;
        communicator comm;

        boost::shared_ptr<CommPattern> C, Cup, Cpu;

        boost::shared_ptr<bmatrix>  K_loc, K_rem, Kup_loc, Kup_rem, Kpu_loc, Kpu_rem;
        boost::shared_ptr<bmatrix>  x2p, x2u, p2x, u2x;
        boost::shared_ptr<matrix> K, Kpu, Kup;
        boost::shared_ptr<vector>  rhs_u, rhs_p, u, p, tmp;

        boost::shared_ptr<USolver> U;
        boost::shared_ptr<PSolver> P;

        template <typename I, typename E>
        void report(const std::string &name, const boost::tuple<I, E> &c) const {
#ifdef AMGCL_DEBUG
            if (comm.rank == 0)
                std::cout << name << " (" << boost::get<0>(c) << ", " << boost::get<1>(c) << ")\n";
#endif
        }

};

} // namespace mpi

namespace backend {

template <class US, class PS, class Alpha, class Beta, class Vec1, class Vec2>
struct spmv_impl< Alpha, mpi::schur_pressure_correction<US, PS>, Vec1, Beta, Vec2>
{
    static void apply(Alpha alpha, const mpi::schur_pressure_correction<US, PS> &A, const Vec1 &x, Beta beta, Vec2 &y)
    {
        A.spmv(alpha, x, beta, y);
    }
};

template <class US, class PS, class Vec1, class Vec2, class Vec3>
struct residual_impl< mpi::schur_pressure_correction<US, PS>, Vec1, Vec2, Vec3>
{
    static void apply(const Vec1 &rhs, const mpi::schur_pressure_correction<US, PS> &A, const Vec2 &x, Vec3 &r)
    {
        backend::copy(rhs, r);
        A.spmv(-1, x, 1, r);
    }
};

} // namespace backend
} // namespace amgcl


#endif
