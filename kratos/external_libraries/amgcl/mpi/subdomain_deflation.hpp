#ifndef AMGCL_MPI_SUBDOMAIN_DEFLATION_HPP
#define AMGCL_MPI_SUBDOMAIN_DEFLATION_HPP

/*
The MIT License

Copyright (c) 2012-2016 Denis Demidov <dennis.demidov@gmail.com>
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

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/range/numeric.hpp>
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
        typedef typename backend_type::matrix     matrix;
        typedef typename backend_type::vector     vector;

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
            TIC("setup deflation");
            typedef backend::crs<value_type, ptrdiff_t>                build_matrix;
            typedef typename backend::row_iterator<Matrix>::type       row_iterator1;
            typedef typename backend::row_iterator<build_matrix>::type row_iterator2;

            // Lets see how many deflation vectors are there.
            std::vector<ptrdiff_t> dv_size(comm.size);
            MPI_Allgather(&ndv, 1, datatype<ptrdiff_t>(), &dv_size[0], 1, datatype<ptrdiff_t>(), comm);
            boost::partial_sum(dv_size, dv_start.begin() + 1);
            nz = dv_start.back();

            df.resize(ndv);
            dx.resize(nz);
            dd = backend_type::create_vector(nz, bprm);

            boost::shared_ptr<build_matrix> aloc = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> arem = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> az   = boost::make_shared<build_matrix>();

            // Get sizes of each domain in comm.
            std::vector<ptrdiff_t> domain = mpi::exclusive_sum(comm, nrows);
            ptrdiff_t loc_beg = domain[comm.rank];
            ptrdiff_t loc_end = domain[comm.rank + 1];

            // Fill deflation vectors.
            TIC("copy deflation vectors");
            {
                std::vector<value_type> z(nrows);
                for(int j = 0; j < ndv; ++j) {
                    for(ptrdiff_t i = 0; i < nrows; ++i)
                        z[i] = prm.def_vec(i, j);
                    Z[j] = backend_type::copy_vector(z, bprm);
                }
            }
            TOC("copy deflation vectors");

            // Number of nonzeros in local and remote parts of the matrix.
            ptrdiff_t loc_nnz = 0, rem_nnz = 0;

            TIC("first pass");
            // First pass over Astrip rows:
            // 1. Count local and remote nonzeros,
            // 3. Build sparsity pattern of matrix AZ.
            az->nrows = nrows;
            az->ncols = nz;
            az->ptr.resize(nrows + 1, 0);

            std::vector<ptrdiff_t> marker(nz, -1);

            for(ptrdiff_t i = 0; i < nrows; ++i) {
                for(row_iterator1 a = backend::row_begin(Astrip, i); a; ++a) {
                    ptrdiff_t c = a.col();

#ifdef AMGCL_SANITY_CHECK
                    mpi::precondition(comm,
                            c >= 0 && c < domain.back(),
                            "Column number is out of bounds"
                            );
#endif

                    ptrdiff_t d = comm.rank; // domain the column belongs to

                    if (loc_beg <= c && c < loc_end) {
                        ++loc_nnz;
                    } else {
                        ++rem_nnz;
                        d = boost::upper_bound(domain, c) - domain.begin() - 1;
                    }

                    if (marker[d] != i) {
                        marker[d] = i;
                        az->ptr[i + 1] += dv_size[d];
                    }
                }
            }
            TOC("first pass");

            TIC("second pass");
            // Second pass over Astrip rows:
            // 1. Build local and remote matrix parts.
            // 2. Build local part of AZ matrix.
            aloc->nrows = nrows;
            aloc->ncols = nrows;
            aloc->ptr.reserve(nrows + 1);
            aloc->col.reserve(loc_nnz);
            aloc->val.reserve(loc_nnz);
            aloc->ptr.push_back(0);

            arem->nrows = nrows;
            // arem->ncols = <not known at this point, does not matter>

            arem->ptr.reserve(nrows + 1);
            arem->col.reserve(rem_nnz);
            arem->val.reserve(rem_nnz);
            arem->ptr.push_back(0);

            boost::partial_sum(az->ptr, az->ptr.begin());
            az->col.resize(az->ptr.back());
            az->val.resize(az->ptr.back());
            boost::fill(marker, -1);

            for(ptrdiff_t i = 0; i < nrows; ++i) {
                ptrdiff_t az_row_beg = az->ptr[i];
                ptrdiff_t az_row_end = az_row_beg;

                for(row_iterator1 a = backend::row_begin(Astrip, i); a; ++a) {
                    ptrdiff_t  c = a.col();
                    value_type v = a.value();

                    if (loc_beg <= c && c < loc_end) {
                        ptrdiff_t loc_c = c - loc_beg;
                        aloc->col.push_back(loc_c);
                        aloc->val.push_back(v);

                        for(ptrdiff_t j = 0, k = dv_start[comm.rank]; j < ndv; ++j, ++k) {
                            if (marker[k] < az_row_beg) {
                                marker[k] = az_row_end;
                                az->col[az_row_end] = k;
                                az->val[az_row_end] = v * prm.def_vec(loc_c, j);
                                ++az_row_end;
                            } else {
                                az->val[marker[k]] += v * prm.def_vec(loc_c, j);
                            }
                        }
                    } else {
                        arem->col.push_back(c);
                        arem->val.push_back(v);
                    }
                }

                az->ptr[i] = az_row_end;

                aloc->ptr.push_back(aloc->col.size());
                arem->ptr.push_back(arem->col.size());
            }
            TOC("second pass");

            // Create local preconditioner.
            P = boost::make_shared<LocalPrecond>( *aloc, prm.local, bprm );

            // Analyze communication pattern, create distributed matrix.
            C = boost::make_shared< comm_pattern<backend_type> >(comm, nrows, arem->col, bprm);
            arem->ncols = C->renumber(arem->col);
            Arem = backend_type::copy_matrix(arem, bprm);
            A = boost::make_shared<dmatrix>(*C, P->system_matrix(), *Arem);

            TIC("A*Z");
            /* Finish construction of AZ */
            // Exchange deflation vectors
            std::vector<ptrdiff_t> zrecv_ptr(C->recv.nbr.size() + 1, 0);
            std::vector<ptrdiff_t> zcol_ptr;
            zcol_ptr.reserve(C->recv.val.size() + 1);
            zcol_ptr.push_back(0);

            for(size_t i = 0; i < C->recv.nbr.size(); ++i) {
                ptrdiff_t ncols = C->recv.ptr[i + 1] - C->recv.ptr[i];
                ptrdiff_t nvecs = dv_size[C->recv.nbr[i]];
                ptrdiff_t size = nvecs * ncols;
                zrecv_ptr[i + 1] = zrecv_ptr[i] + size;

                for(ptrdiff_t j = 0; j < ncols; ++j)
                    zcol_ptr.push_back(zcol_ptr.back() + nvecs);
            }

            std::vector<value_type> zrecv(zrecv_ptr.back());
            std::vector<value_type> zsend(C->send.val.size() * ndv);

            for(size_t i = 0; i < C->recv.nbr.size(); ++i) {
                ptrdiff_t begin = zrecv_ptr[i];
                ptrdiff_t size  = zrecv_ptr[i + 1] - begin;

                MPI_Irecv(&zrecv[begin], size, dtype, C->recv.nbr[i],
                        tag_exc_vals, comm, &C->recv.req[i]);
            }

            for(size_t i = 0, k = 0; i < C->send.col.size(); ++i)
                for(ptrdiff_t j = 0; j < ndv; ++j, ++k)
                    zsend[k] = prm.def_vec(C->send.col[i], j);

            for(size_t i = 0; i < C->send.nbr.size(); ++i)
                MPI_Isend(
                        &zsend[ndv * C->send.ptr[i]], ndv * (C->send.ptr[i+1] - C->send.ptr[i]),
                        dtype, C->send.nbr[i], tag_exc_vals, comm, &C->send.req[i]);

            MPI_Waitall(C->recv.req.size(), &C->recv.req[0], MPI_STATUSES_IGNORE);

            boost::fill(marker, -1);

            // AZ += Arem * Z
            for(ptrdiff_t i = 0; i < nrows; ++i) {
                ptrdiff_t az_row_beg = az->ptr[i];
                ptrdiff_t az_row_end = az_row_beg;

                for(row_iterator2 a = backend::row_begin(*arem, i); a; ++a) {
                    ptrdiff_t  c = a.col();
                    value_type v = a.value();

                    // Domain the column belongs to
                    ptrdiff_t d = C->recv.nbr[boost::upper_bound(C->recv.ptr, c) - C->recv.ptr.begin() - 1];

                    value_type *zval = &zrecv[ zcol_ptr[c] ];
                    for(ptrdiff_t j = 0, k = dv_start[d]; j < dv_size[d]; ++j, ++k) {
                        if (marker[k] < az_row_beg) {
                            marker[k] = az_row_end;
                            az->col[az_row_end] = k;
                            az->val[az_row_end] = v * zval[j];
                            ++az_row_end;
                        } else {
                            az->val[marker[k]] += v * zval[j];
                        }
                    }
                }

                az->ptr[i] = az_row_end;
            }

            std::rotate(az->ptr.begin(), az->ptr.end() - 1, az->ptr.end());
            az->ptr.front() = 0;
            TOC("A*Z");

            MPI_Waitall(C->send.req.size(), &C->send.req[0], MPI_STATUSES_IGNORE);

            /* Build deflated matrix E. */
            TIC("assemble E");
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
                if (j == comm.rank || C->talks_to(j)) {
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
            }

            MPI_Gatherv(
                    &eptr[1], ndv, MPI_INT, &Eptr[0] + 1,
                    const_cast<int*>(&ssize[0]), const_cast<int*>(&sstart[0]),
                    MPI_INT, 0, slaves_comm
                    );

            boost::partial_sum(eptr, eptr.begin());
            boost::partial_sum(Eptr, Eptr.begin());

            // Build local strip of E.
            boost::multi_array<value_type, 2> erow(boost::extents[ndv][nz]);
            std::fill_n(erow.data(), erow.num_elements(), 0);

            for(ptrdiff_t i = 0; i < nrows; ++i) {
                for(row_iterator2 a = backend::row_begin(*az, i); a; ++a) {
                    ptrdiff_t  c = a.col();
                    value_type v = a.value();

                    for(ptrdiff_t j = 0; j < ndv; ++j)
                        erow[j][c] += v * prm.def_vec(i, j);
                }
            }

            std::vector<int>        ecol(eptr.back());
            std::vector<value_type> eval(eptr.back());
            for(int i = 0; i < ndv; ++i) {
                int row_head = eptr[i];
                for(int j = 0; j < comm.size; ++j) {
                    if (j == comm.rank || C->talks_to(j)) {
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
            TOC("assemble E");

            // Prepare E factorization.
            TIC("factorize E");
            if (comm.rank == master_rank) {
                E = boost::make_shared<DirectSolver>(
                        masters_comm, Eptr.size() - 1, Eptr, Ecol, Eval, prm.dsolver
                        );

                cf.resize(Eptr.size() - 1);
                cx.resize(Eptr.size() - 1);
            }
            TOC("factorize E");

            TOC("setup deflation");

            // Move matrices to backend.
            AZ = backend_type::copy_matrix(az, bprm);

            // Prepare Gatherv configuration for coarse solve
            for(int p = cgroup_beg, i = 0, offset = dv_start[p]; p < cgroup_end; ++p, ++i) {
                sstart[i] = dv_start[p] - offset;
                ssize[i]  = dv_start[p + 1] - dv_start[p];
            }

            mstart.resize(nmasters);
            msize.resize(nmasters);
            for(int p = 0; p < nmasters; ++p) {
                mstart[p] = dv_start[slaves[p]];
                msize[p]  = dv_start[slaves[p+1]] - mstart[p];
            }
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
            P->apply(rhs, x);
        }

        const subdomain_deflation& system_matrix() const {
            return *this;
        }

        template <class Vec1, class Vec2>
        boost::tuple<size_t, value_type>
        operator()(const Vec1 &rhs, Vec2 &x) const {
            boost::tuple<size_t, value_type> cnv = S(*this, *P, rhs, x);
            postprocess(rhs, x);
            return cnv;
        }

        template <class Vec1, class Vec2>
        void mul(value_type alpha, const Vec1 &x, value_type beta, Vec2 &y) const {
            TIC("top/spmv");
            backend::spmv(alpha, *A, x, beta, y);
            TOC("top/spmv");

            project(y);
        }

        template <class Vec1, class Vec2, class Vec3>
        void residual(const Vec1 &f, const Vec2 &x, Vec3 &r) const {
            TIC("top/residual");
            backend::residual(f, *A, x, r);
            TOC("top/residual");

            project(r);
        }

        size_t size() const {
            return nrows;
        }
    private:
        typedef distributed_matrix<backend_type> dmatrix;

        static const int tag_exc_vals = 2011;
        static const int tag_exc_dmat = 3011;
        static const int tag_exc_dvec = 4011;
        static const int tag_exc_lnnz = 5011;

        communicator comm;
        ptrdiff_t nrows, ndv, nz;

        MPI_Datatype dtype;

        boost::shared_ptr< comm_pattern<backend_type> > C;
        boost::shared_ptr<matrix> Arem;
        boost::shared_ptr<dmatrix> A;
        boost::shared_ptr<LocalPrecond> P;

        mutable std::vector<value_type> df, dx, cf, cx;
        std::vector<ptrdiff_t> dv_start;
        std::vector<int> sstart, ssize, mstart, msize;

        std::vector< boost::shared_ptr<vector> > Z;

        MPI_Comm masters_comm, slaves_comm;
        int master_rank;
        boost::shared_ptr<DirectSolver> E;

        boost::shared_ptr<matrix> AZ;
        boost::shared_ptr<vector> q;
        boost::shared_ptr<vector> dd;

        ISolver S;

        template <class Vector>
        void project(Vector &x) const {
            TIC("project");

            TIC("local inner product");
            for(ptrdiff_t j = 0; j < ndv; ++j)
                df[j] = backend::inner_product(x, *Z[j]);
            TOC("local inner product");

            coarse_solve(df, dx);

            TIC("spmv");
            backend::copy_to_backend(dx, *dd);
            backend::spmv(-1, *AZ, *dd, 1, x);
            TOC("spmv");

            TOC("project");
        }

        void coarse_solve(std::vector<value_type> &f, std::vector<value_type> &x) const
        {
            TIC("coarse solve");
            TIC("exchange rhs");
            MPI_Gatherv(
                    &f[0], f.size(), dtype, &cf[0],
                    const_cast<int*>(&ssize[0]), const_cast<int*>(&sstart[0]),
                    dtype, 0, slaves_comm
                    );
            TOC("exchange rhs");

            if (comm.rank == master_rank) {
                TIC("call solver");
                (*E)(cf, cx);
                TOC("call solver");

                TIC("gather result");
                MPI_Gatherv(
                        &cx[0], cx.size(), dtype, &x[0],
                        const_cast<int*>(&msize[0]), const_cast<int*>(&mstart[0]),
                        dtype, 0, masters_comm
                        );
                TOC("gather result");
            }

            TIC("broadcast result");
            MPI_Bcast(&x[0], x.size(), dtype, 0, comm);
            TOC("broadcast result");
            TOC("coarse solve");
        }

        template <class Vec1, class Vec2>
        void postprocess(const Vec1 &rhs, Vec2 &x) const {
            TIC("postprocess");

            // q = Ax
            mul(1, x, 0, *q);

            // df = transp(Z) * (rhs - Ax)
            TIC("local inner product");
            for(ptrdiff_t j = 0; j < ndv; ++j)
                df[j] = backend::inner_product(rhs, *Z[j])
                      - backend::inner_product(*q,  *Z[j]);
            TOC("local inner product");

            // dx = inv(E) * df
            coarse_solve(df, dx);

            // x += Z * dx
            ptrdiff_t j = 0, k = dv_start[comm.rank];
            for(; j + 1 < ndv; j += 2, k += 2)
                backend::axpbypcz(dx[k], *Z[j], dx[k+1], *Z[j+1], 1, x);

            for(; j < ndv; ++j, ++k)
                backend::axpby(dx[k], *Z[j], 1, x);

            TOC("postprocess");
        }

};

} // namespace mpi

namespace backend {

template <
    class LocalPrecond,
    template <class, class> class ISolver,
    class DSolver,
    class Alpha, class Beta,
    class Vec1,
    class Vec2
    >
struct spmv_impl<
    Alpha,
    mpi::subdomain_deflation<LocalPrecond, ISolver, DSolver>,
    Vec1, Beta, Vec2
    >
{
    typedef mpi::subdomain_deflation<LocalPrecond, ISolver, DSolver> M;

    static void apply(Alpha alpha, const M &A, const Vec1 &x, Beta beta, Vec2 &y)
    {
        A.mul(alpha, x, beta, y);
    }
};

template <
    class LocalPrecond,
    template <class, class> class ISolver,
    class DSolver,
    class Vec1,
    class Vec2,
    class Vec3
    >
struct residual_impl<
    mpi::subdomain_deflation<LocalPrecond, ISolver, DSolver>,
    Vec1, Vec2, Vec3
    >
{
    typedef mpi::subdomain_deflation<LocalPrecond, ISolver, DSolver> M;

    typedef typename LocalPrecond::backend_type::value_type V;

    static void apply(const Vec1 &rhs, const M &A, const Vec2 &x, Vec3 &r) {
        A.residual(rhs, x, r);
    }
};

} // namespace backend

} // namespace amgcl

#endif
