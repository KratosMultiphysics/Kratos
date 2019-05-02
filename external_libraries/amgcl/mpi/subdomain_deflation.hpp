#ifndef AMGCL_MPI_SUBDOMAIN_DEFLATION_HPP
#define AMGCL_MPI_SUBDOMAIN_DEFLATION_HPP

/*
The MIT License

Copyright (c) 2012-2019 Denis Demidov <dennis.demidov@gmail.com>
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
#include <memory>
#include <functional>

#include <mpi.h>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/direct_solver/skyline_lu.hpp>
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
    class DirectSolver = mpi::direct::skyline_lu<typename LocalPrecond::backend_type::value_type>
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
            std::function<double(ptrdiff_t, unsigned)> def_vec;

            params() {}

#ifndef AMGCL_NO_BOOST
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

                def_vec = *static_cast<std::function<double(ptrdiff_t, unsigned)>*>(ptr);

                check_params(p, {"local", "isolver", "dsolver", "num_def_vec", "def_vec"});
            }

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_CHILD(p, path, local);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, isolver);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, dsolver);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, num_def_vec);
            }
#endif
        };

        typedef typename backend_type::value_type value_type;
        typedef typename backend_type::matrix     bmatrix;
        typedef typename backend_type::vector     vector;
        typedef distributed_matrix<backend_type>  matrix;


        template <class Matrix>
        subdomain_deflation(
                communicator comm,
                const Matrix &Astrip,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
        : comm(comm),
          nrows(backend::rows(Astrip)), ndv(prm.num_def_vec),
          dtype( datatype<value_type>() ), dv_start(comm.size + 1, 0),
          Z( ndv ), q( backend_type::create_vector(nrows, bprm) ),
          S(nrows, prm.isolver, bprm, mpi::inner_product(comm))
        {
            A = std::make_shared<matrix>(comm, Astrip, nrows);
            init(prm, bprm);
        }

        subdomain_deflation(
                communicator comm,
                std::shared_ptr<matrix> A,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
        : comm(comm),
          nrows(A->loc_rows()), ndv(prm.num_def_vec),
          dtype( datatype<value_type>() ), A(A), dv_start(comm.size + 1, 0),
          Z( ndv ), q( backend_type::create_vector(nrows, bprm) ),
          S(nrows, prm.isolver, bprm, mpi::inner_product(comm))
        {
            init(prm, bprm);
        }

        void init(
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
        {
            AMGCL_TIC("setup deflation");
            typedef backend::crs<value_type, ptrdiff_t>                build_matrix;

            // Lets see how many deflation vectors are there.
            std::vector<ptrdiff_t> dv_size(comm.size);
            MPI_Allgather(&ndv, 1, datatype<ptrdiff_t>(), &dv_size[0], 1, datatype<ptrdiff_t>(), comm);

            std::partial_sum(dv_size.begin(), dv_size.end(), dv_start.begin() + 1);
            nz = dv_start.back();

            df.resize(ndv);
            dx.resize(ndv);
            dd = backend_type::create_vector(ndv, bprm);

            auto az_loc = std::make_shared<build_matrix>();
            auto az_rem = std::make_shared<build_matrix>();

            auto a_loc = A->local();
            auto a_rem = A->remote();

            const comm_pattern<backend_type> &Acp = A->cpat();

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
            az_loc->set_size(nrows, ndv, true);
            az_loc->set_nonzeros(nrows * dv_size[comm.rank]);
            az_rem->set_size(nrows, 0, true);
            // 1. Build local part of AZ matrix.
            // 2. Count remote nonzeros
#pragma omp parallel
            {
                std::vector<ptrdiff_t> marker(Acp.recv.nbr.size(), -1);

#pragma omp for
                for(ptrdiff_t i = 0; i < nrows; ++i) {
                    ptrdiff_t az_loc_head = i * ndv;
                    az_loc->ptr[i+1] = az_loc_head + ndv;

                    for(ptrdiff_t j = 0; j < ndv; ++j) {
                        az_loc->col[az_loc_head + j] = j;
                        az_loc->val[az_loc_head + j] = math::zero<value_type>();
                    }

                    for(ptrdiff_t j = a_loc->ptr[i], e = a_loc->ptr[i+1]; j < e; ++j) {
                        ptrdiff_t  c = a_loc->col[j];
                        value_type v = a_loc->val[j];

                        for(ptrdiff_t j = 0; j < ndv; ++j)
                            az_loc->val[az_loc_head + j] += v * prm.def_vec(c, j);
                    }

                    for(ptrdiff_t j = a_rem->ptr[i], e = a_rem->ptr[i+1]; j < e; ++j) {
                        int d = Acp.domain(a_rem->col[j]);

                        if (marker[d] != i) {
                            marker[d] = i;
                            az_rem->ptr[i+1] += dv_size[d];
                        }
                    }
                }
            }
            az_rem->set_nonzeros(az_rem->scan_row_sizes());
            AMGCL_TOC("first pass");

            // Create local preconditioner.
            AMGCL_TIC("local preconditioner");
            P = std::make_shared<LocalPrecond>( *a_loc, prm.local, bprm );
            AMGCL_TOC("local preconditioner");

            A->set_local(P->system_matrix_ptr());
            A->move_to_backend(bprm);

            AMGCL_TIC("remote(A*Z)");
            /* Construct remote part of AZ */
            // Exchange deflation vectors
            std::vector<ptrdiff_t> zrecv_ptr(Acp.recv.nbr.size() + 1, 0);
            std::vector<ptrdiff_t> zcol_ptr;
            zcol_ptr.reserve(Acp.recv.count() + 1);
            zcol_ptr.push_back(0);

            for(size_t i = 0; i < Acp.recv.nbr.size(); ++i) {
                ptrdiff_t ncols = Acp.recv.ptr[i + 1] - Acp.recv.ptr[i];
                ptrdiff_t nvecs = dv_size[Acp.recv.nbr[i]];
                ptrdiff_t size = nvecs * ncols;
                zrecv_ptr[i + 1] = zrecv_ptr[i] + size;

                for(ptrdiff_t j = 0; j < ncols; ++j)
                    zcol_ptr.push_back(zcol_ptr.back() + nvecs);
            }

            std::vector<value_type> zrecv(zrecv_ptr.back());
            std::vector<value_type> zsend(Acp.send.count() * ndv);

            for(size_t i = 0; i < Acp.recv.nbr.size(); ++i) {
                ptrdiff_t begin = zrecv_ptr[i];
                ptrdiff_t size  = zrecv_ptr[i + 1] - begin;

                MPI_Irecv(&zrecv[begin], size, dtype, Acp.recv.nbr[i],
                        tag_exc_vals, comm, &Acp.recv.req[i]);
            }

            for(size_t i = 0, k = 0; i < Acp.send.count(); ++i)
                for(ptrdiff_t j = 0; j < ndv; ++j, ++k)
                    zsend[k] = prm.def_vec(Acp.send.col[i], j);

            for(size_t i = 0; i < Acp.send.nbr.size(); ++i)
                MPI_Isend(
                        &zsend[ndv * Acp.send.ptr[i]], ndv * (Acp.send.ptr[i+1] - Acp.send.ptr[i]),
                        dtype, Acp.send.nbr[i], tag_exc_vals, comm, &Acp.send.req[i]);

            MPI_Waitall(Acp.recv.req.size(), &Acp.recv.req[0], MPI_STATUSES_IGNORE);
            MPI_Waitall(Acp.send.req.size(), &Acp.send.req[0], MPI_STATUSES_IGNORE);

#pragma omp parallel
            {
                std::vector<ptrdiff_t> marker(nz, -1);

                // AZ_rem = Arem * Z
#pragma omp for
                for(ptrdiff_t i = 0; i < nrows; ++i) {
                    ptrdiff_t az_rem_head = az_rem->ptr[i];
                    ptrdiff_t az_rem_tail = az_rem_head;

                    for(auto a = backend::row_begin(*a_rem, i); a; ++a) {
                        ptrdiff_t  c = a.col();
                        value_type v = a.value();

                        // Domain the column belongs to
                        ptrdiff_t d = Acp.recv.nbr[
                            std::upper_bound(Acp.recv.ptr.begin(), Acp.recv.ptr.end(), c) -
                                Acp.recv.ptr.begin() - 1];

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

            /* Build solver for the deflated matrix E. */
            AMGCL_TIC("assemble E");

            // Count nonzeros in E.
            std::vector<int> nbrs; // processes we are talking to
            nbrs.reserve(1 + Acp.send.nbr.size() + Acp.recv.nbr.size());
            std::set_union(
                    Acp.send.nbr.begin(), Acp.send.nbr.end(),
                    Acp.recv.nbr.begin(), Acp.recv.nbr.end(),
                    std::back_inserter(nbrs));
            nbrs.push_back(comm.rank);

            build_matrix E;
            E.set_size(ndv, nz, false);

            {
                ptrdiff_t nnz = 0;
                for(int j : nbrs) nnz += dv_size[j];
                for(int k = 0; k <= ndv; ++k)
                    E.ptr[k] = k * nnz;
            }
            E.set_nonzeros(E.nnz = E.ptr[ndv]);

            // Build local strip of E.
#ifdef _OPENMP
            int nthreads = omp_get_max_threads();
#else
            int nthreads = 1;
#endif
            multi_array<value_type, 3> erow(nthreads, ndv, nz);
            std::fill_n(erow.data(), erow.size(), 0);

            {
                ptrdiff_t dv_offset = dv_start[comm.rank];
#pragma omp parallel
                {
#ifdef _OPENMP
                    const int tid = omp_get_thread_num();
#else
                    const int tid = 0;
#endif
                    std::vector<value_type> z(ndv);

#pragma omp for
                    for(ptrdiff_t i = 0; i < nrows; ++i) {
                        for(ptrdiff_t j = 0; j < ndv; ++j)
                            z[j] = prm.def_vec(i,j);

                        for(ptrdiff_t k = az_loc->ptr[i], e = az_loc->ptr[i+1]; k < e; ++k) {
                            ptrdiff_t  c = az_loc->col[k] + dv_offset;
                            value_type v = az_loc->val[k];

                            for(ptrdiff_t j = 0; j < ndv; ++j)
                                erow(tid, j, c) += v * z[j];
                        }

                        for(ptrdiff_t k = az_rem->ptr[i], e = az_rem->ptr[i+1]; k < e; ++k) {
                            ptrdiff_t  c = az_rem->col[k];
                            value_type v = az_rem->val[k];

                            for(ptrdiff_t j = 0; j < ndv; ++j)
                                erow(tid, j, c) += v * z[j];
                        }
                    }
                }
            }

            for(int i = 0; i < ndv; ++i) {
                int row_head = E.ptr[i];
                for(int j : nbrs) {
                    for(int k = 0; k < dv_size[j]; ++k) {
                        int c = dv_start[j] + k;
                        value_type v = math::zero<value_type>();
                        for(int t = 0; t < nthreads; ++t)
                            v += erow(t, i, c);

                        E.col[row_head] = c;
                        E.val[row_head] = v;

                        ++row_head;
                    }
                }
            }
            AMGCL_TOC("assemble E");

            AMGCL_TIC("factorize E");
            this->E = std::make_shared<DirectSolver>(comm, E, prm.dsolver);
            AMGCL_TOC("factorize E");

            AMGCL_TIC("finish(A*Z)");
            AZ = std::make_shared<matrix>(comm, az_loc, az_rem);
            AZ->move_to_backend(bprm);
            AMGCL_TOC("finish(A*Z)");
            AMGCL_TOC("setup deflation");
        }

        template <class Vec1, class Vec2>
        void apply(const Vec1 &rhs, Vec2 &&x) const {
            size_t iters;
            double error;
            backend::clear(x);
            std::tie(iters, error) = (*this)(rhs, x);
        }

        std::shared_ptr<matrix> system_matrix_ptr() const {
            return A;
        }

        const matrix& system_matrix() const {
            return *A;
        }

        template <class Matrix, class Vec1, class Vec2>
        std::tuple<size_t, value_type> operator()(
                const Matrix &A, const Vec1 &rhs, Vec2 &&x) const
        {
            std::tuple<size_t, value_type> cnv = S(make_sdd_projected_matrix(*this, A), *P, rhs, x);
            postprocess(rhs, x);
            return cnv;
        }

        template <class Vec1, class Vec2>
        std::tuple<size_t, value_type>
        operator()(const Vec1 &rhs, Vec2 &&x) const {
            std::tuple<size_t, value_type> cnv = S(make_sdd_projected_matrix(*this, *A), *P, rhs, x);
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
            backend::copy(dx, *dd);
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

        std::shared_ptr<matrix> A, AZ;
        std::shared_ptr<LocalPrecond> P;

        mutable std::vector<value_type> df, dx;
        std::vector<ptrdiff_t> dv_start;

        std::vector< std::shared_ptr<vector> > Z;

        std::shared_ptr<DirectSolver> E;

        std::shared_ptr<vector> q;
        std::shared_ptr<vector> dd;

        ISolver S;

        void coarse_solve(std::vector<value_type> &f, std::vector<value_type> &x) const
        {
            AMGCL_TIC("coarse solve");
            (*E)(f, x);
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
