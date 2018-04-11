#ifndef AMGCL_MPI_SUBDOMAIN_DEFLATION_HPP
#define AMGCL_MPI_SUBDOMAIN_DEFLATION_HPP

/*
The MIT License

Copyright (c) 2012-2017 Denis Demidov <dennis.demidov@gmail.com>
Copyright (c) 2014-2015, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)

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
 * \file   amgcl/mpi/subdomain_deflatedion.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Distributed solver based on subdomain deflation.
 */

#include <vector>
#include <algorithm>
#include <numeric>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/multi_array.hpp>
#include <boost/function.hpp>

#include <mpi.h>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/skyline_lu.hpp>
#include <amgcl/mpi/inner_product.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>

namespace amgcl {
namespace mpi {

/// Pointwise constant deflation vectors.
struct constant_deflation {
    const int block_size;
    /// Constructor
    /**
     * \param block_size Number of degrees of freedom per grid point
     */
    constant_deflation(int block_size = 1) : block_size(block_size) {}

    int dim() const {
        return block_size;
    }

    int operator()(ptrdiff_t row, int j) const {
        return row % block_size == j;
    }
};

template <class SDD, class Matrix>
struct sdd_projected_matrix {
    typedef typename SDD::value_type value_type;

    const SDD    &S;
    const Matrix &A;

    sdd_projected_matrix(const SDD &S, const Matrix &A) : S(S), A(A) {}

    template <class T, class Vec1, class Vec2>
    void mul(T alpha, const Vec1 &x, T beta, Vec2 &y) const {
        AMGCL_TIC("top/spmv");
        backend::spmv(alpha, A, x, beta, y);
        AMGCL_TOC("top/spmv");

        S.project(y);
    }

    template <class Vec1, class Vec2, class Vec3>
    void residual(const Vec1 &f, const Vec2 &x, Vec3 &r) const {
        AMGCL_TIC("top/residual");
        backend::residual(f, A, x, r);
        AMGCL_TOC("top/residual");

        S.project(r);
    }
};

template <class SDD, class Matrix>
sdd_projected_matrix<SDD, Matrix> make_sdd_projected_matrix(const SDD &S, const Matrix &A) {
    return sdd_projected_matrix<SDD, Matrix>(S, A);
}

/// Distributed solver based on subdomain deflation.
/**
 * \sa \cite Frank2001
 */
template <
    class LocalPrecond,
    template <class, class> class IterativeSolver,
    class DirectSolver = mpi::skyline_lu<typename LocalPrecond::backend_type::value_type>
    >
class subdomain_deflation {
    public:
        typedef typename LocalPrecond::backend_type backend_type;
        typedef typename backend_type::params backend_params;
        typedef IterativeSolver<backend_type, mpi::inner_product> ISolver;

        struct params {
            typename LocalPrecond::params local;
            typename ISolver::params      isolver;
            typename DirectSolver::params dsolver;

            // Number of deflation vectors.
            unsigned num_def_vec;

            // Value of deflation vector at the given row and column.
            boost::function<double(ptrdiff_t, unsigned)> def_vec;

            params() {}

            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_CHILD(p, local),
                  AMGCL_PARAMS_IMPORT_CHILD(p, isolver),
                  AMGCL_PARAMS_IMPORT_CHILD(p, dsolver),
                  AMGCL_PARAMS_IMPORT_VALUE(p, num_def_vec)
            {
                void *ptr = 0;
                ptr = p.get("def_vec", ptr);

                amgcl::precondition(ptr,
                        "Error in subdomain_deflation parameters: "
                        "def_vec is not set");

                def_vec = *static_cast<boost::function<double(ptrdiff_t, unsigned)>*>(ptr);

                AMGCL_PARAMS_CHECK(p, (local)(isolver)(dsolver)(num_def_vec)(def_vec));
            }

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_CHILD(p, path, local);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, isolver);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, dsolver);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, num_def_vec);
            }
        };

        typedef typename backend_type::value_type value_type;
        typedef typename backend_type::matrix     bmatrix;
        typedef typename backend_type::vector     vector;
        typedef distributed_matrix<backend_type>  matrix;


        template <class Matrix>
        subdomain_deflation(
                MPI_Comm mpi_comm,
                const Matrix &Astrip,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
        : comm(mpi_comm),
          nrows(backend::rows(Astrip)), ndv(prm.num_def_vec),
          dtype( datatype<value_type>() ), dv_start(comm.size + 1, 0),
          Z( ndv ), master_rank(0),
          q( backend_type::create_vector(nrows, bprm) ),
          S(nrows, prm.isolver, bprm, mpi::inner_product(mpi_comm))
        {
            AMGCL_TIC("setup deflation");
            typedef backend::crs<value_type, ptrdiff_t>                build_matrix;
            typedef typename backend::row_iterator<Matrix>::type       row_iterator1;
            typedef typename backend::row_iterator<build_matrix>::type row_iterator2;

            // Lets see how many deflation vectors are there.
            std::vector<ptrdiff_t> dv_size(comm.size);
            MPI_Allgather(&ndv, 1, datatype<ptrdiff_t>(), &dv_size[0], 1, datatype<ptrdiff_t>(), comm);
            std::partial_sum(dv_size.begin(), dv_size.end(), dv_start.begin() + 1);
            nz = dv_start.back();

            df.resize(ndv);
            dx.resize(ndv);
            dd = backend_type::create_vector(ndv, bprm);

            boost::shared_ptr<build_matrix> a_loc  = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> a_rem  = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> az_loc = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> az_rem = boost::make_shared<build_matrix>();

            // Get sizes of each domain in comm.
            std::vector<ptrdiff_t> domain = mpi::exclusive_sum(comm, nrows);
            ptrdiff_t loc_beg = domain[comm.rank];
            ptrdiff_t loc_end = domain[comm.rank + 1];

            // Fill deflation vectors.
            AMGCL_TIC("copy deflation vectors");
            {
                std::vector<value_type> z(nrows);
                for(int j = 0; j < ndv; ++j) {
#pragma omp parallel for
                    for(ptrdiff_t i = 0; i < nrows; ++i)
                        z[i] = prm.def_vec(i, j);
                    Z[j] = backend_type::copy_vector(z, bprm);
                }
            }
            AMGCL_TOC("copy deflation vectors");

            AMGCL_TIC("first pass");
            // First pass over Astrip rows:
            // 1. Count local and remote nonzeros
            a_loc->set_size(nrows, nrows, true);
            az_loc->set_size(nrows, ndv, true);
            // ncols for remote parts are not known at this point,
            // will update later:
            a_rem->set_size(nrows, 0, true);
            az_rem->set_size(nrows, 0, true);

#pragma omp parallel
            {
                std::vector<ptrdiff_t> marker(nz, -1);

#pragma omp for
                for(ptrdiff_t i = 0; i < nrows; ++i) {
                    for(row_iterator1 a = backend::row_begin(Astrip, i); a; ++a) {
                        ptrdiff_t c = a.col();

#ifdef AMGCL_SANITY_CHECK
                        precondition(comm, c >= 0 && c < domain.back(), "Column number is out of bounds");
#endif

                        ptrdiff_t d = comm.rank; // domain the column belongs to
                        bool local = (loc_beg <= c && c < loc_end);

                        if (local) {
                            ++a_loc->ptr[i+1];
                        } else {
                            ++a_rem->ptr[i+1];
                            d = std::upper_bound(domain.begin(), domain.end(), c) - domain.begin() - 1;
                        }

                        if (marker[d] != i) {
                            marker[d] = i;

                            if (local)
                                az_loc->ptr[i+1] += dv_size[d];
                            else
                                az_rem->ptr[i+1] += dv_size[d];
                        }
                    }
                }
            }
            AMGCL_TOC("first pass");

            AMGCL_TIC("second pass");
            // Second pass over Astrip rows:
            // 1. Build local and remote matrix parts.
            // 2. Build local part of AZ matrix.
            std::partial_sum(a_loc->ptr, a_loc->ptr + nrows + 1, a_loc->ptr);
            std::partial_sum(a_rem->ptr, a_rem->ptr + nrows + 1, a_rem->ptr);
            std::partial_sum(az_loc->ptr, az_loc->ptr + nrows + 1, az_loc->ptr);
            std::partial_sum(az_rem->ptr, az_rem->ptr + nrows + 1, az_rem->ptr);

            a_loc->set_nonzeros(a_loc->ptr[nrows]);
            a_rem->set_nonzeros(a_rem->ptr[nrows]);
            az_loc->set_nonzeros(az_loc->ptr[nrows]);
            az_rem->set_nonzeros(az_rem->ptr[nrows]);

#pragma omp parallel
            {
                std::vector<ptrdiff_t> marker(ndv, -1);

#pragma omp for
                for(ptrdiff_t i = 0; i < nrows; ++i) {
                    ptrdiff_t a_loc_head = a_loc->ptr[i];
                    ptrdiff_t a_rem_head = a_rem->ptr[i];

                    ptrdiff_t az_loc_head = az_loc->ptr[i];
                    ptrdiff_t az_loc_tail = az_loc_head;

                    for(row_iterator1 a = backend::row_begin(Astrip, i); a; ++a) {
                        ptrdiff_t  c = a.col();
                        value_type v = a.value();

                        if (loc_beg <= c && c < loc_end) {
                            ptrdiff_t loc_c = c - loc_beg;
                            a_loc->col[a_loc_head] = loc_c;
                            a_loc->val[a_loc_head] = v;
                            ++a_loc_head;

                            for(ptrdiff_t j = 0; j < ndv; ++j) {
                                if (marker[j] < az_loc_head) {
                                    marker[j] = az_loc_tail;
                                    assert(az_loc_tail < az_loc->nnz);
                                    az_loc->col[az_loc_tail] = j;
                                    az_loc->val[az_loc_tail] = v * prm.def_vec(loc_c, j);
                                    ++az_loc_tail;
                                } else {
                                    az_loc->val[marker[j]] += v * prm.def_vec(loc_c, j);
                                }
                            }
                        } else {
                            a_rem->col[a_rem_head] = c;
                            a_rem->val[a_rem_head] = v;
                            ++a_rem_head;
                        }
                    }
                }
            }
            AMGCL_TOC("second pass");

            // Create local preconditioner.
            P = boost::make_shared<LocalPrecond>( *a_loc, prm.local, bprm );

            // Analyze communication pattern for A, create distributed matrix.
            Acp = boost::make_shared< comm_pattern<backend_type> >(comm, nrows, a_rem->nnz, a_rem->col, bprm);
            a_rem->ncols = Acp->renumber(a_rem->nnz, a_rem->col);
            Arem = backend_type::copy_matrix(a_rem, bprm);
            A = boost::make_shared<matrix>(*Acp, P->system_matrix(), *Arem);

            AMGCL_TIC("remote(A*Z)");
            /* Construct remote part of AZ */
            // Exchange deflation vectors
            std::vector<ptrdiff_t> zrecv_ptr(Acp->recv.nbr.size() + 1, 0);
            std::vector<ptrdiff_t> zcol_ptr;
            zcol_ptr.reserve(Acp->recv.val.size() + 1);
            zcol_ptr.push_back(0);

            for(size_t i = 0; i < Acp->recv.nbr.size(); ++i) {
                ptrdiff_t ncols = Acp->recv.ptr[i + 1] - Acp->recv.ptr[i];
                ptrdiff_t nvecs = dv_size[Acp->recv.nbr[i]];
                ptrdiff_t size = nvecs * ncols;
                zrecv_ptr[i + 1] = zrecv_ptr[i] + size;

                for(ptrdiff_t j = 0; j < ncols; ++j)
                    zcol_ptr.push_back(zcol_ptr.back() + nvecs);
            }

            std::vector<value_type> zrecv(zrecv_ptr.back());
            std::vector<value_type> zsend(Acp->send.val.size() * ndv);

            for(size_t i = 0; i < Acp->recv.nbr.size(); ++i) {
                ptrdiff_t begin = zrecv_ptr[i];
                ptrdiff_t size  = zrecv_ptr[i + 1] - begin;

                MPI_Irecv(&zrecv[begin], size, dtype, Acp->recv.nbr[i],
                        tag_exc_vals, comm, &Acp->recv.req[i]);
            }

            for(size_t i = 0, k = 0; i < Acp->send.col.size(); ++i)
                for(ptrdiff_t j = 0; j < ndv; ++j, ++k)
                    zsend[k] = prm.def_vec(Acp->send.col[i], j);

            for(size_t i = 0; i < Acp->send.nbr.size(); ++i)
                MPI_Isend(
                        &zsend[ndv * Acp->send.ptr[i]], ndv * (Acp->send.ptr[i+1] - Acp->send.ptr[i]),
                        dtype, Acp->send.nbr[i], tag_exc_vals, comm, &Acp->send.req[i]);

            MPI_Waitall(Acp->recv.req.size(), &Acp->recv.req[0], MPI_STATUSES_IGNORE);
            MPI_Waitall(Acp->send.req.size(), &Acp->send.req[0], MPI_STATUSES_IGNORE);

#pragma omp parallel
            {
                std::vector<ptrdiff_t> marker(nz, -1);

                // AZ_rem = Arem * Z
#pragma omp for
                for(ptrdiff_t i = 0; i < nrows; ++i) {
                    ptrdiff_t az_rem_head = az_rem->ptr[i];
                    ptrdiff_t az_rem_tail = az_rem_head;

                    for(row_iterator2 a = backend::row_begin(*a_rem, i); a; ++a) {
                        ptrdiff_t  c = a.col();
                        value_type v = a.value();

                        // Domain the column belongs to
                        ptrdiff_t d = Acp->recv.nbr[
                            std::upper_bound(Acp->recv.ptr.begin(), Acp->recv.ptr.end(), c) -
                                Acp->recv.ptr.begin() - 1];

                        value_type *zval = &zrecv[ zcol_ptr[c] ];
                        for(ptrdiff_t j = 0, k = dv_start[d]; j < dv_size[d]; ++j, ++k) {
                            if (marker[k] < az_rem_head) {
                                marker[k] = az_rem_tail;
                                az_rem->col[az_rem_tail] = k;
                                az_rem->val[az_rem_tail] = v * zval[j];
                                ++az_rem_tail;
                            } else {
                                az_rem->val[marker[k]] += v * zval[j];
                            }
                        }
                    }
                }
            }
            AMGCL_TOC("remote(A*Z)");


            /* Build deflated matrix E. */
            AMGCL_TIC("assemble E");
            // Who is responsible for solution of coarse problem
            int nmasters = std::min(comm.size, DirectSolver::comm_size(nz, prm.dsolver));
            int nslaves  = (comm.size + nmasters - 1) / nmasters;
            int cgroup_beg = 0, cgroup_end = 0;
            std::vector<int> slaves(nmasters + 1, 0);
            for(int p = 1; p <= nmasters; ++p) {
                slaves[p] = std::min(p * nslaves, comm.size);
                if (slaves[p-1] <= comm.rank && comm.rank < slaves[p]) {
                    cgroup_beg  = slaves[p - 1];
                    cgroup_end  = slaves[p];
                    master_rank = slaves[p-1];
                    nslaves     = cgroup_end - cgroup_beg;
                }
            }

            // Communicator for masters (used to solve the coarse problem):
            MPI_Comm_split(comm,
                    comm.rank == master_rank ? 0 : MPI_UNDEFINED,
                    comm.rank, &masters_comm
                    );

            // Communicator for slaves (used to send/recv coarse data):
            MPI_Comm_split(comm, master_rank, comm.rank, &slaves_comm);

            // Count nonzeros in E.
            std::vector<int> eptr(ndv + 1, 0);
            for(int j = 0; j < comm.size; ++j) {
                if (j == comm.rank || Acp->talks_to(j)) {
                    for(int k = 0; k < ndv; ++k)
                        eptr[k + 1] += dv_size[j];
                }
            }

            sstart.resize(nslaves);
            ssize.resize(nslaves);
            for(int p = cgroup_beg, i = 0, offset = dv_start[p]; p < cgroup_end; ++p, ++i) {
                sstart[i] = dv_start[p] - offset;
                ssize[i]  = dv_start[p + 1] - dv_start[p];
            }


            std::vector<int> Eptr;
            if (comm.rank == master_rank) {
                Eptr.resize(dv_start[cgroup_end] - dv_start[cgroup_beg] + 1, 0);

                MPI_Gatherv(
                        &eptr[1], ndv, MPI_INT, &Eptr[0] + 1,
                        const_cast<int*>(&ssize[0]), const_cast<int*>(&sstart[0]),
                        MPI_INT, 0, slaves_comm
                        );
            } else {
                MPI_Gatherv(
                        &eptr[1], ndv, MPI_INT, NULL, NULL, NULL,
                        MPI_INT, 0, slaves_comm
                        );
            }

            std::partial_sum(eptr.begin(), eptr.end(), eptr.begin());
            std::partial_sum(Eptr.begin(), Eptr.end(), Eptr.begin());

            // Build local strip of E.
            boost::multi_array<value_type, 2> erow(boost::extents[ndv][nz]);
            std::fill_n(erow.data(), erow.num_elements(), 0);

            {
                ptrdiff_t loc_dv_start = dv_start[comm.rank];
                std::vector<value_type> z(ndv);
                for(ptrdiff_t i = 0; i < nrows; ++i) {
                    for(ptrdiff_t j = 0; j < ndv; ++j)
                        z[j] = prm.def_vec(i,j);

                    for(row_iterator2 a = backend::row_begin(*az_loc, i); a; ++a) {
                        ptrdiff_t  c = a.col() + loc_dv_start;
                        value_type v = a.value();

                        for(ptrdiff_t j = 0; j < ndv; ++j)
                            erow[j][c] += v * z[j];
                    }

                    for(row_iterator2 a = backend::row_begin(*az_rem, i); a; ++a) {
                        ptrdiff_t  c = a.col();
                        value_type v = a.value();

                        for(ptrdiff_t j = 0; j < ndv; ++j)
                            erow[j][c] += v * z[j];
                    }
                }
            }

            std::vector<int>        ecol(eptr.back());
            std::vector<value_type> eval(eptr.back());
            for(int i = 0; i < ndv; ++i) {
                int row_head = eptr[i];
                for(int j = 0; j < comm.size; ++j) {
                    if (j == comm.rank || Acp->talks_to(j)) {
                        for(int k = 0; k < dv_size[j]; ++k) {
                            int c = dv_start[j] + k;
                            ecol[row_head] = c;
                            eval[row_head] = erow[i][c];
                            ++row_head;
                        }
                    }
                }
            }

            // Exchange strips of E.
            std::vector<int>        Ecol;
            std::vector<value_type> Eval;
            if (comm.rank == master_rank) {
                Ecol.resize(Eptr.back());
                Eval.resize(Eptr.back());

                for(int p = cgroup_beg, i = 0, offset = dv_start[p]; p < cgroup_end; ++p, ++i) {
                    sstart[i] = Eptr[dv_start[p]     - offset];
                    ssize[i]  = Eptr[dv_start[p + 1] - offset] - sstart[i];
                }

                MPI_Gatherv(
                        &ecol[0], ecol.size(), MPI_INT, &Ecol[0],
                        const_cast<int*>(&ssize[0]), const_cast<int*>(&sstart[0]),
                        MPI_INT, 0, slaves_comm
                        );

                MPI_Gatherv(
                        &eval[0], eval.size(), dtype, &Eval[0],
                        const_cast<int*>(&ssize[0]), const_cast<int*>(&sstart[0]),
                        dtype, 0, slaves_comm
                        );
            } else {
                MPI_Gatherv(
                        &ecol[0], ecol.size(), MPI_INT, NULL, NULL, NULL,
                        MPI_INT, 0, slaves_comm
                        );

                MPI_Gatherv(
                        &eval[0], eval.size(), dtype, NULL, NULL, NULL,
                        dtype, 0, slaves_comm
                        );
            }
            AMGCL_TOC("assemble E");

            // Prepare E factorization.
            AMGCL_TIC("factorize E");
            if (comm.rank == master_rank) {
                E = boost::make_shared<DirectSolver>(
                        masters_comm, Eptr.size() - 1, Eptr, Ecol, Eval, prm.dsolver
                        );

                cf.resize(Eptr.size() - 1);
                cx.resize(Eptr.size() - 1);
            }
            AMGCL_TOC("factorize E");

            AMGCL_TIC("finish(A*Z)");
            AZcp = boost::make_shared< comm_pattern<backend_type> >(comm, ndv, az_rem->nnz, az_rem->col, bprm);
            az_rem->ncols = AZcp->renumber(az_rem->nnz, az_rem->col);
            AZloc = backend_type::copy_matrix(az_loc, bprm);
            AZrem = backend_type::copy_matrix(az_rem, bprm);
            AZ = boost::make_shared<matrix>(*AZcp, *AZloc, *AZrem);
            AMGCL_TOC("finish(A*Z)");

            // Prepare Gatherv configuration for coarse solve
            for(int p = cgroup_beg, i = 0, offset = dv_start[p]; p < cgroup_end; ++p, ++i) {
                sstart[i] = dv_start[p] - offset;
                ssize[i]  = dv_start[p + 1] - dv_start[p];
            }
            AMGCL_TOC("setup deflation");
        }

        ~subdomain_deflation() {
            E.reset();
            if (masters_comm != MPI_COMM_NULL) MPI_Comm_free(&masters_comm);
            if (slaves_comm  != MPI_COMM_NULL) MPI_Comm_free(&slaves_comm);
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
            size_t iters;
            double error;
            backend::clear(x);
            boost::tie(iters, error) = (*this)(rhs, x);
        }

        const matrix& system_matrix() const {
            return *A;
        }

        template <class Matrix, class Vec1, class Vec2>
        boost::tuple<size_t, value_type> operator()(
                Matrix  const &A,
                Vec1    const &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2          &x
#else
                Vec2          &&x
#endif
                ) const
        {
            boost::tuple<size_t, value_type> cnv = S(make_sdd_projected_matrix(*this, A), *P, rhs, x);
            postprocess(rhs, x);
            return cnv;
        }

        template <class Vec1, class Vec2>
        boost::tuple<size_t, value_type>
        operator()(
                const Vec1 &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2          &x
#else
                Vec2          &&x
#endif
                ) const
        {
            boost::tuple<size_t, value_type> cnv = S(make_sdd_projected_matrix(*this, *A), *P, rhs, x);
            postprocess(rhs, x);
            return cnv;
        }

        size_t size() const {
            return nrows;
        }

        template <class Vector>
        void project(Vector &x) const {
            AMGCL_TIC("project");

            AMGCL_TIC("local inner product");
            for(ptrdiff_t j = 0; j < ndv; ++j)
                df[j] = backend::inner_product(x, *Z[j]);
            AMGCL_TOC("local inner product");

            coarse_solve(df, dx);

            AMGCL_TIC("spmv");
            backend::copy_to_backend(dx, *dd);
            backend::spmv(-1, *AZ, *dd, 1, x);
            AMGCL_TOC("spmv");

            AMGCL_TOC("project");
        }
    private:
        static const int tag_exc_vals = 2011;
        static const int tag_exc_dmat = 3011;
        static const int tag_exc_dvec = 4011;
        static const int tag_exc_lnnz = 5011;

        communicator comm;
        ptrdiff_t nrows, ndv, nz;

        MPI_Datatype dtype;

        boost::shared_ptr< comm_pattern<backend_type> > Acp, AZcp;
        boost::shared_ptr<bmatrix> Arem, AZloc, AZrem;
        boost::shared_ptr<matrix> A, AZ;
        boost::shared_ptr<LocalPrecond> P;

        mutable std::vector<value_type> df, dx, cf, cx;
        std::vector<ptrdiff_t> dv_start;
        std::vector<int> sstart, ssize;

        std::vector< boost::shared_ptr<vector> > Z;

        MPI_Comm masters_comm, slaves_comm;
        int master_rank;
        boost::shared_ptr<DirectSolver> E;

        boost::shared_ptr<vector> q;
        boost::shared_ptr<vector> dd;

        ISolver S;

        void coarse_solve(std::vector<value_type> &f, std::vector<value_type> &x) const
        {
            AMGCL_TIC("coarse solve");
            AMGCL_TIC("exchange rhs");
            if (comm.rank == master_rank) {
                MPI_Gatherv(&f[0], ndv, dtype, &cf[0],
                        const_cast<int*>(&ssize[0]), const_cast<int*>(&sstart[0]),
                        dtype, 0, slaves_comm);
            } else {
                MPI_Gatherv(&f[0], f.size(), dtype, NULL, NULL, NULL,
                        dtype, 0, slaves_comm);
            }
            AMGCL_TOC("exchange rhs");

            if (comm.rank == master_rank) {
                AMGCL_TIC("call solver");
                (*E)(cf, cx);
                AMGCL_TOC("call solver");

                AMGCL_TIC("scatter result");
                MPI_Scatterv(&cx[0],
                        const_cast<int*>(&ssize[0]), const_cast<int*>(&sstart[0]),
                        dtype, &x[0], ndv, dtype, 0, slaves_comm);
                AMGCL_TOC("scatter result");
            } else {
                MPI_Scatterv(NULL, NULL, NULL, dtype,
                        &x[0], ndv, dtype, 0, slaves_comm);
            }

            AMGCL_TOC("coarse solve");
        }

        template <class Vec1, class Vec2>
        void postprocess(const Vec1 &rhs, Vec2 &x) const {
            AMGCL_TIC("postprocess");

            // q = rhs - Ax
            backend::copy(rhs, *q);
            backend::spmv(-1, *A, x, 1, *q);

            // df = transp(Z) * (rhs - Ax)
            AMGCL_TIC("local inner product");
            for(ptrdiff_t j = 0; j < ndv; ++j)
                df[j] = backend::inner_product(*q, *Z[j]);
            AMGCL_TOC("local inner product");

            // dx = inv(E) * df
            coarse_solve(df, dx);

            // x += Z * dx
            backend::lin_comb(ndv, dx, Z, 1, x);

            AMGCL_TOC("postprocess");
        }

};

} // namespace mpi

namespace backend {

template <
    class SDD, class Matrix,
    class Alpha, class Beta, class Vec1, class Vec2
    >
struct spmv_impl<
    Alpha, mpi::sdd_projected_matrix<SDD, Matrix>, Vec1, Beta, Vec2
    >
{
    typedef mpi::sdd_projected_matrix<SDD, Matrix> M;

    static void apply(Alpha alpha, const M &A, const Vec1 &x, Beta beta, Vec2 &y)
    {
        A.mul(alpha, x, beta, y);
    }
};

template <class SDD, class Matrix, class Vec1, class Vec2, class Vec3>
struct residual_impl<mpi::sdd_projected_matrix<SDD, Matrix>, Vec1, Vec2, Vec3>
{
    typedef mpi::sdd_projected_matrix<SDD, Matrix> M;

    static void apply(const Vec1 &rhs, const M &A, const Vec2 &x, Vec3 &r) {
        A.residual(rhs, x, r);
    }
};

} // namespace backend

} // namespace amgcl

#endif
