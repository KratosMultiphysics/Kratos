#ifndef AMGCL_MPI_DISTRIBUTED_MATRIX_HPP
#define AMGCL_MPI_DISTRIBUTED_MATRIX_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

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

#include <vector>
#include <algorithm>

#include <memory>
#include <unordered_map>
#include <random>

#include <mpi.h>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>
#include <amgcl/mpi/util.hpp>

/**
 * \file   amgcl/mpi/distributed_matrix.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Distributed matrix implementation.
 */

namespace amgcl {
namespace mpi {

template <class Backend>
class comm_pattern {
    public:
        typedef typename Backend::value_type value_type;
        typedef typename math::rhs_of<value_type>::type rhs_type;
        typedef typename math::scalar_of<value_type>::type scalar_type;
        typedef typename Backend::matrix matrix;
        typedef typename Backend::vector vector;
        typedef typename Backend::params backend_params;

        struct {
            std::vector<ptrdiff_t> nbr;
            std::vector<ptrdiff_t> ptr;
            std::vector<ptrdiff_t> col;

            size_t count() const {
                return col.size();
            }

            mutable std::vector<rhs_type>    val;
            mutable std::vector<MPI_Request> req;
        } send;

        struct {
            std::vector<ptrdiff_t> nbr;
            std::vector<ptrdiff_t> ptr;

            size_t count() const {
                return val.size();
            }

            mutable std::vector<rhs_type>    val;
            mutable std::vector<MPI_Request> req;
        } recv;

        std::shared_ptr<vector> x_rem;

        comm_pattern(
                communicator comm,
                ptrdiff_t n_loc_cols,
                size_t n_rem_cols, const ptrdiff_t *p_rem_cols
                ) : comm(comm), loc_cols(n_loc_cols)
        {
            AMGCL_TIC("communication pattern");
            // Get domain boundaries
            std::vector<ptrdiff_t> domain = comm.exclusive_sum(n_loc_cols);
            loc_beg = domain[comm.rank];

            // Renumber remote columns,
            // find out how many remote values we need from each process.
            std::vector<ptrdiff_t> rem_cols(p_rem_cols, p_rem_cols + n_rem_cols);

            std::sort(rem_cols.begin(), rem_cols.end());
            rem_cols.erase(std::unique(rem_cols.begin(), rem_cols.end()), rem_cols.end());

            ptrdiff_t ncols = rem_cols.size();
            ptrdiff_t rnbr = 0, snbr = 0, send_size = 0;

            {
                std::vector<int> rcounts(comm.size, 0);
                std::vector<int> scounts(comm.size);

                // Build index for column renumbering;
                // count how many domains send us data and how much.
                idx.reserve(2 * ncols);
                for(int i = 0, d = 0, last = -1; i < ncols; ++i) {
                    while(rem_cols[i] >= domain[d + 1]) ++d;

                    ++rcounts[d];

                    if (last < d) {
                        last = d;
                        ++rnbr;
                    }

                    idx.insert(idx.end(), std::make_pair(
                                rem_cols[i], std::make_tuple(rnbr-1,i)));
                }

                recv.val.resize(ncols);
                recv.req.resize(rnbr);

                recv.nbr.reserve(rnbr);
                recv.ptr.reserve(rnbr + 1); recv.ptr.push_back(0);

                for(int d = 0; d < comm.size; ++d) {
                    if (rcounts[d]) {
                        recv.nbr.push_back(d);
                        recv.ptr.push_back(recv.ptr.back() + rcounts[d]);
                    }
                }

                MPI_Alltoall(&rcounts[0], 1, MPI_INT, &scounts[0], 1, MPI_INT, comm);

                for(ptrdiff_t d = 0; d < comm.size; ++d) {
                    if (scounts[d]) {
                        ++snbr;
                        send_size += scounts[d];
                    }
                }

                send.col.resize(send_size);
                send.val.resize(send_size);
                send.req.resize(snbr);

                send.nbr.reserve(snbr);
                send.ptr.reserve(snbr + 1); send.ptr.push_back(0);

                for(ptrdiff_t d = 0; d < comm.size; ++d) {
                    if (scounts[d]) {
                        send.nbr.push_back(d);
                        send.ptr.push_back(send.ptr.back() + scounts[d]);
                    }
                }
            }

            // What columns do you need from me?
            for(size_t i = 0; i < send.nbr.size(); ++i)
                MPI_Irecv(&send.col[send.ptr[i]], send.ptr[i+1] - send.ptr[i],
                        datatype<ptrdiff_t>(), send.nbr[i], tag_exc_cols, comm, &send.req[i]);

            // Here is what I need from you:
            for(size_t i = 0; i < recv.nbr.size(); ++i)
                MPI_Isend(&rem_cols[recv.ptr[i]], recv.ptr[i+1] - recv.ptr[i],
                        datatype<ptrdiff_t>(), recv.nbr[i], tag_exc_cols, comm, &recv.req[i]);

            AMGCL_TIC("MPI Wait");
            MPI_Waitall(recv.req.size(), &recv.req[0], MPI_STATUSES_IGNORE);
            MPI_Waitall(send.req.size(), &send.req[0], MPI_STATUSES_IGNORE);
            AMGCL_TOC("MPI Wait");

            // Shift columns to send to local numbering:
            for(ptrdiff_t &c : send.col) c -= loc_beg;

            AMGCL_TOC("communication pattern");
        }

        void move_to_backend(const backend_params &bprm = backend_params()) {
            x_rem  = Backend::create_vector(recv.count(), bprm);
            gather = std::make_shared<Gather>(loc_cols, send.col, bprm);
        }

        int domain(ptrdiff_t col) const {
            return std::get<0>(idx.at(col));
        }

        int local_index(ptrdiff_t col) const {
            return std::get<1>(idx.at(col));
        }

        std::tuple<int, int> remote_info(ptrdiff_t col) const {
            return idx.at(col);
        }

        size_t renumber(size_t n, ptrdiff_t *col) const {
            for(size_t i = 0; i < n; ++i)
                col[i] = std::get<1>(idx.at(col[i]));
            return recv.count();
        }

        bool needs_remote() const {
            return !recv.val.empty();
        }

        template <class Vector>
        void start_exchange(const Vector &x) const {
            // Start receiving ghost values from our neighbours.
            for(size_t i = 0; i < recv.nbr.size(); ++i)
                MPI_Irecv(&recv.val[recv.ptr[i]], recv.ptr[i+1] - recv.ptr[i],
                        datatype<rhs_type>(), recv.nbr[i], tag_exc_vals, comm, &recv.req[i]);

            // Start sending our data to neighbours.
            if (!send.val.empty()) {
                (*gather)(x, send.val);

                for(size_t i = 0; i < send.nbr.size(); ++i)
                    MPI_Isend(&send.val[send.ptr[i]], send.ptr[i+1] - send.ptr[i],
                            datatype<rhs_type>(), send.nbr[i], tag_exc_vals, comm, &send.req[i]);
            }
        }

        void finish_exchange() const {
            AMGCL_TIC("MPI Wait");
            MPI_Waitall(recv.req.size(), &recv.req[0], MPI_STATUSES_IGNORE);
            MPI_Waitall(send.req.size(), &send.req[0], MPI_STATUSES_IGNORE);
            AMGCL_TOC("MPI Wait");

            if (!recv.val.empty())
                backend::copy(recv.val, *x_rem);
        }

        template <typename T>
        void exchange(const T *send_val, T *recv_val) const {
            for(size_t i = 0; i < recv.nbr.size(); ++i)
                MPI_Irecv(&recv_val[recv.ptr[i]], recv.ptr[i+1] - recv.ptr[i],
                        datatype<T>(), recv.nbr[i], tag_exc_vals, comm, &recv.req[i]);

            for(size_t i = 0; i < send.nbr.size(); ++i)
                MPI_Isend(const_cast<T*>(&send_val[send.ptr[i]]), send.ptr[i+1] - send.ptr[i],
                        datatype<T>(), send.nbr[i], tag_exc_vals, comm, &send.req[i]);

            AMGCL_TIC("MPI Wait");
            MPI_Waitall(recv.req.size(), &recv.req[0], MPI_STATUSES_IGNORE);
            MPI_Waitall(send.req.size(), &send.req[0], MPI_STATUSES_IGNORE);
            AMGCL_TOC("MPI Wait");
        }

        communicator mpi_comm() const {
            return comm;
        }

        ptrdiff_t loc_col_shift() const {
            return loc_beg;
        }

    private:
        typedef typename Backend::gather Gather;

        static const int tag_set_comm = 1001;
        static const int tag_exc_cols = 1002;
        static const int tag_exc_vals = 1003;

        communicator comm;

        std::unordered_map<ptrdiff_t, std::tuple<int, int> > idx;
        std::shared_ptr<Gather> gather;
        ptrdiff_t loc_beg, loc_cols;
};

template <class Backend>
class distributed_matrix {
    public:
        typedef typename Backend::value_type value_type;
        typedef typename math::rhs_of<value_type>::type rhs_type;
        typedef typename math::scalar_of<value_type>::type scalar_type;
        typedef typename Backend::params backend_params;
        typedef typename Backend::matrix matrix;
        typedef comm_pattern<Backend> CommPattern;
        typedef backend::crs<value_type> build_matrix;

        distributed_matrix(
                communicator comm,
                std::shared_ptr<build_matrix> a_loc,
                std::shared_ptr<build_matrix> a_rem,
                std::shared_ptr<CommPattern>  c = std::shared_ptr<CommPattern>()
                )
            : a_loc(a_loc), a_rem(a_rem)
        {
            if (c) {
                C = c;
            } else {
                C = std::make_shared<CommPattern>(comm, a_loc->ncols, a_rem->nnz, a_rem->col);
            }

            a_rem->ncols = C->recv.count();

            n_loc_rows     = a_loc->nrows;
            n_loc_cols     = a_loc->ncols;
            n_loc_nonzeros = a_loc->nnz + a_rem->nnz;

            n_glob_rows     = comm.reduce(MPI_SUM, n_loc_rows);
            n_glob_cols     = comm.reduce(MPI_SUM, n_loc_cols);
            n_glob_nonzeros = comm.reduce(MPI_SUM, n_loc_nonzeros);
        }

        template <class Matrix>
        distributed_matrix(
                communicator comm,
                const Matrix &A,
                ptrdiff_t _n_loc_cols = -1)
            : n_loc_rows(backend::rows(A)),
              n_loc_cols(_n_loc_cols < 0 ? n_loc_rows : _n_loc_cols),
              n_loc_nonzeros(backend::nonzeros(A))
        {
            // Get sizes of each domain in comm.
            std::vector<ptrdiff_t> domain = comm.exclusive_sum(n_loc_cols);
            ptrdiff_t loc_beg = domain[comm.rank];
            ptrdiff_t loc_end = domain[comm.rank + 1];

            n_glob_cols     = domain.back();
            n_glob_rows     = comm.reduce(MPI_SUM, n_loc_rows);
            n_glob_nonzeros = comm.reduce(MPI_SUM, n_loc_nonzeros);

            // Split the matrix into local and remote parts.
            a_loc = std::make_shared<build_matrix>();
            a_rem = std::make_shared<build_matrix>();

            build_matrix &A_loc = *a_loc;
            build_matrix &A_rem = *a_rem;

            A_loc.set_size(n_loc_rows, n_loc_cols, true);
            A_rem.set_size(n_loc_rows, 0, true);

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n_loc_rows; ++i) {
                for(auto a = backend::row_begin(A, i); a; ++a) {
                    ptrdiff_t c = a.col();

                    if (loc_beg <= c && c < loc_end)
                        ++A_loc.ptr[i + 1];
                    else
                        ++A_rem.ptr[i + 1];
                }
            }

            A_loc.set_nonzeros(A_loc.scan_row_sizes());
            A_rem.set_nonzeros(A_rem.scan_row_sizes());

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n_loc_rows; ++i) {
                ptrdiff_t loc_head = A_loc.ptr[i];
                ptrdiff_t rem_head = A_rem.ptr[i];

                for(auto a = backend::row_begin(A, i); a; ++a) {
                    ptrdiff_t  c = a.col();
                    value_type v = a.value();

                    if (loc_beg <= c && c < loc_end) {
                        A_loc.col[loc_head] = c - loc_beg;
                        A_loc.val[loc_head] = v;
                        ++loc_head;
                    } else {
                        A_rem.col[rem_head] = c;
                        A_rem.val[rem_head] = v;
                        ++rem_head;
                    }
                }
            }

            C = std::make_shared<CommPattern>(comm, n_loc_cols, a_rem->nnz, a_rem->col);
            a_rem->ncols = C->recv.count();
        }

        communicator comm() const {
            return C->mpi_comm();
        }

        std::shared_ptr<build_matrix> local() const {
            return a_loc;
        }

        std::shared_ptr<build_matrix> remote() const {
            return a_rem;
        }

        std::shared_ptr<matrix> local_backend() const {
            return A_loc;
        }

        std::shared_ptr<matrix> remote_backend() const {
            return A_rem;
        }

        ptrdiff_t loc_rows() const {
            return n_loc_rows;
        }

        ptrdiff_t loc_cols() const {
            return n_loc_cols;
        }

        ptrdiff_t loc_col_shift() const {
            return C->loc_col_shift();
        }

        ptrdiff_t loc_nonzeros() const {
            return n_loc_nonzeros;
        }

        ptrdiff_t glob_rows() const {
            return n_glob_rows;
        }

        ptrdiff_t glob_cols() const {
            return n_glob_cols;
        }

        ptrdiff_t glob_nonzeros() const {
            return n_glob_nonzeros;
        }

        const comm_pattern<Backend>& cpat() const {
            return *C;
        }

        void set_local(std::shared_ptr<matrix> a) {
            A_loc = a;
        }

        void move_to_backend(const backend_params &bprm = backend_params()) {
            AMGCL_TIC("move to backend");
            if (!A_loc) {
                A_loc = Backend::copy_matrix(a_loc, bprm);
            }

            if (!A_rem) {
                C->renumber(a_rem->nnz, a_rem->col);
                A_rem = Backend::copy_matrix(a_rem, bprm);
            }

            C->move_to_backend(bprm);

            a_loc.reset();
            a_rem.reset();
            AMGCL_TOC("move to backend");
        }

        template <class A, class VecX, class B, class VecY>
        void mul(A alpha, const VecX &x, B beta, VecY &y) const {
            C->start_exchange(x);

            // Compute local part of the product.
            backend::spmv(alpha, *A_loc, x, beta, y);

            // Compute remote part of the product.
            C->finish_exchange();

            if (C->needs_remote())
                backend::spmv(alpha, *A_rem, *C->x_rem, 1, y);
        }

        template <class Vec1, class Vec2, class Vec3>
        void residual(const Vec1 &f, const Vec2 &x, Vec3 &r) const {
            C->start_exchange(x);
            backend::residual(f, *A_loc, x, r);

            C->finish_exchange();

            if (C->needs_remote())
                backend::spmv(-1, *A_rem, *C->x_rem, 1, r);
        }

    private:
        std::shared_ptr<CommPattern>  C;
        std::shared_ptr<matrix> A_loc, A_rem;
        std::shared_ptr<build_matrix> a_loc, a_rem;

        ptrdiff_t n_loc_rows, n_glob_rows;
        ptrdiff_t n_loc_cols, n_glob_cols;
        ptrdiff_t n_loc_nonzeros, n_glob_nonzeros;
};

template <class Backend>
std::shared_ptr< distributed_matrix<Backend> >
transpose(const distributed_matrix<Backend> &A) {
    AMGCL_TIC("MPI Transpose");
    typedef typename Backend::value_type value_type;
    typedef comm_pattern<Backend>        CommPattern;
    typedef backend::crs<value_type>     build_matrix;

    static const int tag_cnt = 2001;
    static const int tag_col = 2002;
    static const int tag_val = 2003;

    communicator comm = A.comm();
    const CommPattern &C = A.cpat();

    build_matrix &A_loc = *A.local();
    build_matrix &A_rem = *A.remote();

    ptrdiff_t nrows = A_loc.ncols;
    ptrdiff_t ncols = A_loc.nrows;

    std::vector<MPI_Request> recv_cnt_req(C.send.req.size());
    std::vector<MPI_Request> recv_col_req(C.send.req.size());
    std::vector<MPI_Request> recv_val_req(C.send.req.size());

    std::vector<MPI_Request> send_cnt_req(C.recv.req.size());
    std::vector<MPI_Request> send_col_req(C.recv.req.size());
    std::vector<MPI_Request> send_val_req(C.recv.req.size());

    // Our transposed remote part becomes remote part of someone else,
    // and the other way around.
    std::shared_ptr<build_matrix> t_ptr;
    {
        std::vector<ptrdiff_t> tmp_col(A_rem.col, A_rem.col + A_rem.nnz);
        C.renumber(tmp_col.size(), &tmp_col[0]);

        ptrdiff_t *a_rem_col = &tmp_col[0];
        std::swap(a_rem_col, A_rem.col);

        t_ptr = backend::transpose(A_rem);

        std::swap(a_rem_col, A_rem.col);
    }
    build_matrix &t_rem = *t_ptr;

    // Shift to global numbering:
    std::vector<ptrdiff_t> domain = comm.exclusive_sum(ncols);
    ptrdiff_t loc_beg = domain[comm.rank];
    for(size_t i = 0; i < t_rem.nnz; ++i)
        t_rem.col[i] += loc_beg;

    // Shift from row pointers to row sizes:
    std::vector<ptrdiff_t> row_size(t_rem.nrows);
    for(size_t i = 0; i < t_rem.nrows; ++i)
        row_size[i] = t_rem.ptr[i+1] - t_rem.ptr[i];

    // Sizes of transposed remote blocks:
    // 1. Exchange rem_ptr
    std::vector<ptrdiff_t> rem_ptr(C.send.count() + 1); rem_ptr[0] = 0;

    for(size_t i = 0; i < C.send.nbr.size(); ++i) {
        ptrdiff_t beg = C.send.ptr[i];
        ptrdiff_t end = C.send.ptr[i + 1];

        MPI_Irecv(&rem_ptr[beg + 1], end - beg, datatype<ptrdiff_t>(),
                C.send.nbr[i], tag_cnt, comm, &recv_cnt_req[i]);
    }

    for(size_t i = 0; i < C.recv.nbr.size(); ++i) {
        ptrdiff_t beg = C.recv.ptr[i];
        ptrdiff_t end = C.recv.ptr[i + 1];

        MPI_Isend(&row_size[beg], end - beg, datatype<ptrdiff_t>(),
                C.recv.nbr[i], tag_cnt, comm, &send_cnt_req[i]);
    }

    AMGCL_TIC("MPI Wait");
    MPI_Waitall(recv_cnt_req.size(), &recv_cnt_req[0], MPI_STATUSES_IGNORE);
    AMGCL_TOC("MPI Wait");
    std::partial_sum(rem_ptr.begin(), rem_ptr.end(), rem_ptr.begin());

    // 2. Start exchange of rem_col, rem_val
    std::vector<ptrdiff_t>  rem_col(rem_ptr.back());
    std::vector<value_type> rem_val(rem_ptr.back());

    for(size_t i = 0; i < C.send.nbr.size(); ++i) {
        ptrdiff_t rbeg = C.send.ptr[i];
        ptrdiff_t rend = C.send.ptr[i + 1];

        ptrdiff_t cbeg = rem_ptr[rbeg];
        ptrdiff_t cend = rem_ptr[rend];

        MPI_Irecv(&rem_col[cbeg], cend - cbeg, datatype<ptrdiff_t>(),
                C.send.nbr[i], tag_col, comm, &recv_col_req[i]);

        MPI_Irecv(&rem_val[cbeg], cend - cbeg, datatype<value_type>(),
                C.send.nbr[i], tag_val, comm, &recv_val_req[i]);
    }

    for(size_t i = 0; i < C.recv.nbr.size(); ++i) {
        ptrdiff_t rbeg = C.recv.ptr[i];
        ptrdiff_t rend = C.recv.ptr[i + 1];

        ptrdiff_t cbeg = t_rem.ptr[rbeg];
        ptrdiff_t cend = t_rem.ptr[rend];

        MPI_Isend(&t_rem.col[cbeg], cend - cbeg, datatype<ptrdiff_t>(),
                C.recv.nbr[i], tag_col, comm, &send_col_req[i]);

        MPI_Isend(&t_rem.val[cbeg], cend - cbeg, datatype<value_type>(),
                C.recv.nbr[i], tag_val, comm, &send_val_req[i]);
    }

    // 3. While rem_col and rem_val are in flight,
    //    start constructing our remote part:
    auto T_ptr = std::make_shared<build_matrix>();
    build_matrix &T_rem = *T_ptr;
    T_rem.set_size(nrows, 0, true);

    for(size_t i = 0; i < C.send.count(); ++i)
        T_rem.ptr[1 + C.send.col[i]] += rem_ptr[i+1] - rem_ptr[i];

    T_rem.scan_row_sizes();
    T_rem.set_nonzeros();

    // 4. Finish rem_col and rem_val exchange, and
    //    finish contruction of our remote part.
    AMGCL_TIC("MPI Wait");
    MPI_Waitall(recv_col_req.size(), &recv_col_req[0], MPI_STATUSES_IGNORE);
    MPI_Waitall(recv_val_req.size(), &recv_val_req[0], MPI_STATUSES_IGNORE);
    AMGCL_TOC("MPI Wait");

    for(size_t i = 0; i < C.send.count(); ++i) {
        ptrdiff_t row  = C.send.col[i];
        ptrdiff_t head = T_rem.ptr[row];

        for(ptrdiff_t j = rem_ptr[i]; j < rem_ptr[i+1]; ++j, ++head) {
            T_rem.col[head] = rem_col[j];
            T_rem.val[head] = rem_val[j];
        }

        T_rem.ptr[row] = head;
    }

    std::rotate(T_rem.ptr, T_rem.ptr + nrows, T_rem.ptr + nrows + 1);
    T_rem.ptr[0] = 0;

    AMGCL_TIC("MPI Wait");
    MPI_Waitall(send_cnt_req.size(), &send_cnt_req[0], MPI_STATUSES_IGNORE);
    MPI_Waitall(send_col_req.size(), &send_col_req[0], MPI_STATUSES_IGNORE);
    MPI_Waitall(send_val_req.size(), &send_val_req[0], MPI_STATUSES_IGNORE);
    AMGCL_TOC("MPI Wait");

    AMGCL_TOC("MPI Transpose");

    return std::make_shared< distributed_matrix<Backend> >(
            comm, backend::transpose(A_loc), T_ptr);
}

template <class Backend>
std::shared_ptr< backend::crs<typename Backend::value_type> >
remote_rows(
        const comm_pattern<Backend> &C,
        const distributed_matrix<Backend> &B,
        bool need_values = true
        )
{
    typedef typename Backend::value_type value_type;
    typedef backend::crs<value_type>     build_matrix;

    static const int tag_ptr = 3001;
    static const int tag_col = 3002;
    static const int tag_val = 3003;

    AMGCL_TIC("remote_rows");
    communicator comm = C.mpi_comm();

    build_matrix &B_loc = *B.local();
    build_matrix &B_rem = *B.remote();
    ptrdiff_t B_beg = B.loc_col_shift();

    size_t nrecv = C.recv.nbr.size();
    size_t nsend = C.send.nbr.size();

    // Create blocked matrix to send to each domain
    // that needs data from us:
    std::vector<MPI_Request> send_ptr_req(nsend);
    std::vector<MPI_Request> send_col_req(nsend);
    std::vector<MPI_Request> send_val_req(nsend);

    std::vector<build_matrix> send_rows(nsend);

    for(size_t k = 0; k < nsend; ++k) {
        ptrdiff_t beg = C.send.ptr[k];
        ptrdiff_t end = C.send.ptr[k + 1];

        build_matrix &m = send_rows[k];
        m.set_size(end - beg, 0, false);

        for(ptrdiff_t i = 0, ii = beg; ii < end; ++i, ++ii) {
            ptrdiff_t r = C.send.col[ii];

            ptrdiff_t w =
                (B_loc.ptr[r + 1] - B_loc.ptr[r]) +
                (B_rem.ptr[r + 1] - B_rem.ptr[r]);

            m.ptr[i] = w;
            m.nnz += w;
        }

        MPI_Isend(m.ptr, m.nrows, datatype<ptrdiff_t>(),
                C.send.nbr[k], tag_ptr, comm, &send_ptr_req[k]);

        m.set_nonzeros(m.nnz, need_values);

        for(ptrdiff_t i = 0, ii = beg, head = 0; ii < end; ++i, ++ii) {
            ptrdiff_t r = C.send.col[ii];

            // Contribution of the local part:
            for(ptrdiff_t j = B_loc.ptr[r]; j < B_loc.ptr[r+1]; ++j) {
                m.col[head] = B_loc.col[j] + B_beg;

                if (need_values) m.val[head] = B_loc.val[j];

                ++head;
            }

            // Contribution of the remote part:
            for(ptrdiff_t j = B_rem.ptr[r]; j < B_rem.ptr[r+1]; ++j) {
                m.col[head] = B_rem.col[j];

                if (need_values) m.val[head] = B_rem.val[j];

                ++head;
            }
        }

        MPI_Isend(m.col, m.nnz, datatype<ptrdiff_t>(),
                C.send.nbr[k], tag_col, comm, &send_col_req[k]);
        if (need_values)
            MPI_Isend(m.val, m.nnz, datatype<value_type>(),
                    C.send.nbr[k], tag_val, comm, &send_val_req[k]);
    }

    // Receive rows of B in block format from our neighbors:
    std::vector<MPI_Request> recv_ptr_req(nrecv);
    std::vector<MPI_Request> recv_col_req(nrecv);
    std::vector<MPI_Request> recv_val_req(nrecv);

    auto B_nbr = std::make_shared<build_matrix>();
    B_nbr->set_size(C.recv.count(), 0, false);
    B_nbr->ptr[0] = 0;

    for(size_t k = 0; k < nrecv; ++k) {
        ptrdiff_t beg = C.recv.ptr[k];
        ptrdiff_t end = C.recv.ptr[k + 1];

        MPI_Irecv(&B_nbr->ptr[beg + 1], end - beg, datatype<ptrdiff_t>(),
                C.recv.nbr[k], tag_ptr, comm, &recv_ptr_req[k]);
    }

    AMGCL_TIC("MPI Wait");
    MPI_Waitall(recv_ptr_req.size(), &recv_ptr_req[0], MPI_STATUSES_IGNORE);
    AMGCL_TOC("MPI Wait");

    B_nbr->set_nonzeros(B_nbr->scan_row_sizes(), need_values);

    for(size_t k = 0; k < nrecv; ++k) {
        ptrdiff_t rbeg = C.recv.ptr[k];
        ptrdiff_t rend = C.recv.ptr[k + 1];

        ptrdiff_t cbeg = B_nbr->ptr[rbeg];
        ptrdiff_t cend = B_nbr->ptr[rend];

        MPI_Irecv(&B_nbr->col[cbeg], cend - cbeg, datatype<ptrdiff_t>(),
                C.recv.nbr[k], tag_col, comm, &recv_col_req[k]);

        if (need_values)
            MPI_Irecv(&B_nbr->val[cbeg], cend - cbeg, datatype<value_type>(),
                    C.recv.nbr[k], tag_val, comm, &recv_val_req[k]);
    }

    AMGCL_TIC("MPI Wait");
    MPI_Waitall(send_ptr_req.size(), &send_ptr_req[0], MPI_STATUSES_IGNORE);
    MPI_Waitall(send_col_req.size(), &send_col_req[0], MPI_STATUSES_IGNORE);
    MPI_Waitall(recv_col_req.size(), &recv_col_req[0], MPI_STATUSES_IGNORE);

    if (need_values) {
        MPI_Waitall(send_val_req.size(), &send_val_req[0], MPI_STATUSES_IGNORE);
        MPI_Waitall(recv_val_req.size(), &recv_val_req[0], MPI_STATUSES_IGNORE);
    }
    AMGCL_TOC("MPI Wait");

    AMGCL_TOC("remote_rows");
    return B_nbr;
}

template <class Backend>
std::shared_ptr< distributed_matrix<Backend> >
product(const distributed_matrix<Backend> &A, const distributed_matrix<Backend> &B) {
    typedef typename Backend::value_type value_type;
    typedef backend::crs<value_type>     build_matrix;
    AMGCL_TIC("product");

    const comm_pattern<Backend> &Acp = A.cpat();

    build_matrix &A_loc = *A.local();
    build_matrix &A_rem = *A.remote();
    build_matrix &B_loc = *B.local();
    build_matrix &B_rem = *B.remote();

    ptrdiff_t A_rows = A.loc_rows();
    ptrdiff_t B_cols = B.loc_cols();

    ptrdiff_t B_beg = B.loc_col_shift();
    ptrdiff_t B_end = B_beg + B_cols;

    auto b_nbr = remote_rows(Acp, B);
    build_matrix &B_nbr = *b_nbr;

    // Build mapping from global to local column numbers in the remote part of
    // the product matrix.
    std::vector<ptrdiff_t> rem_cols(B_rem.nnz + B_nbr.nnz);

    std::copy(B_nbr.col, B_nbr.col + B_nbr.nnz,
            std::copy(B_rem.col, B_rem.col + B_rem.nnz, rem_cols.begin()));

    std::sort(rem_cols.begin(), rem_cols.end());
    rem_cols.erase(std::unique(rem_cols.begin(), rem_cols.end()), rem_cols.end());

    ptrdiff_t n_rem_cols = 0;
    std::unordered_map<ptrdiff_t, int> rem_idx(2 * rem_cols.size());
    for(ptrdiff_t c : rem_cols) {
        if (c >= B_beg && c < B_end) continue;
        rem_idx[c] = n_rem_cols++;
    }

    // Build the product.
    auto c_loc = std::make_shared<build_matrix>();
    auto c_rem = std::make_shared<build_matrix>();

    build_matrix &C_loc = *c_loc;
    build_matrix &C_rem = *c_rem;

    C_loc.set_size(A_rows, B_cols, false);
    C_rem.set_size(A_rows, 0,      false);

    C_loc.ptr[0] = 0;
    C_rem.ptr[0] = 0;

    AMGCL_TIC("analyze");
#pragma omp parallel
    {
        std::vector<ptrdiff_t> loc_marker(B_end - B_beg, -1);
        std::vector<ptrdiff_t> rem_marker(n_rem_cols,    -1);

#pragma omp for
        for(ptrdiff_t ia = 0; ia < A_rows; ++ia) {
            ptrdiff_t loc_cols = 0;
            ptrdiff_t rem_cols = 0;

            for(ptrdiff_t ja = A_loc.ptr[ia], ea = A_loc.ptr[ia + 1]; ja < ea; ++ja) {
                ptrdiff_t  ca = A_loc.col[ja];

                for(ptrdiff_t jb = B_loc.ptr[ca], eb = B_loc.ptr[ca+1]; jb < eb; ++jb) {
                    ptrdiff_t  cb = B_loc.col[jb];

                    if (loc_marker[cb] != ia) {
                        loc_marker[cb]  = ia;
                        ++loc_cols;
                    }
                }

                for(ptrdiff_t jb = B_rem.ptr[ca], eb = B_rem.ptr[ca+1]; jb < eb; ++jb) {
                    ptrdiff_t  cb = rem_idx[B_rem.col[jb]];

                    if (rem_marker[cb] != ia) {
                        rem_marker[cb]  = ia;
                        ++rem_cols;
                    }
                }
            }

            for(ptrdiff_t ja = A_rem.ptr[ia], ea = A_rem.ptr[ia + 1]; ja < ea; ++ja) {
                ptrdiff_t  ca = Acp.local_index(A_rem.col[ja]);

                for(ptrdiff_t jb = B_nbr.ptr[ca], eb = B_nbr.ptr[ca+1]; jb < eb; ++jb) {
                    ptrdiff_t  cb = B_nbr.col[jb];

                    if (cb >= B_beg && cb < B_end) {
                        cb -= B_beg;

                        if (loc_marker[cb] != ia) {
                            loc_marker[cb]  = ia;
                            ++loc_cols;
                        }
                    } else {
                        cb = rem_idx[cb];

                        if (rem_marker[cb] != ia) {
                            rem_marker[cb]  = ia;
                            ++rem_cols;
                        }
                    }
                }
            }

            C_loc.ptr[ia + 1] = loc_cols;
            C_rem.ptr[ia + 1] = rem_cols;
        }
    }
    AMGCL_TOC("analyze");

    C_loc.set_nonzeros(C_loc.scan_row_sizes());
    C_rem.set_nonzeros(C_rem.scan_row_sizes());

    AMGCL_TIC("compute");
#pragma omp parallel
    {
        std::vector<ptrdiff_t> loc_marker(B_end - B_beg, -1);
        std::vector<ptrdiff_t> rem_marker(n_rem_cols,    -1);

#pragma omp for
        for(ptrdiff_t ia = 0; ia < A_rows; ++ia) {
            ptrdiff_t loc_beg = C_loc.ptr[ia];
            ptrdiff_t rem_beg = C_rem.ptr[ia];
            ptrdiff_t loc_end = loc_beg;
            ptrdiff_t rem_end = rem_beg;

            for(ptrdiff_t ja = A_loc.ptr[ia], ea = A_loc.ptr[ia + 1]; ja < ea; ++ja) {
                ptrdiff_t  ca = A_loc.col[ja];
                value_type va = A_loc.val[ja];

                for(ptrdiff_t jb = B_loc.ptr[ca], eb = B_loc.ptr[ca+1]; jb < eb; ++jb) {
                    ptrdiff_t  cb = B_loc.col[jb];
                    value_type vb = B_loc.val[jb];

                    if (loc_marker[cb] < loc_beg) {
                        loc_marker[cb] = loc_end;

                        C_loc.col[loc_end] = cb;
                        C_loc.val[loc_end] = va * vb;

                        ++loc_end;
                    } else {
                        C_loc.val[loc_marker[cb]] += va * vb;
                    }
                }

                for(ptrdiff_t jb = B_rem.ptr[ca], eb = B_rem.ptr[ca+1]; jb < eb; ++jb) {
                    ptrdiff_t  gb = B_rem.col[jb];
                    ptrdiff_t  cb = rem_idx[gb];
                    value_type vb = B_rem.val[jb];

                    if (rem_marker[cb] < rem_beg) {
                        rem_marker[cb] = rem_end;

                        C_rem.col[rem_end] = gb;
                        C_rem.val[rem_end] = va * vb;

                        ++rem_end;
                    } else {
                        C_rem.val[rem_marker[cb]] += va * vb;
                    }
                }
            }

            for(ptrdiff_t ja = A_rem.ptr[ia], ea = A_rem.ptr[ia + 1]; ja < ea; ++ja) {
                ptrdiff_t  ca = Acp.local_index(A_rem.col[ja]);
                value_type va = A_rem.val[ja];

                for(ptrdiff_t jb = B_nbr.ptr[ca], eb = B_nbr.ptr[ca+1]; jb < eb; ++jb) {
                    ptrdiff_t  gb = B_nbr.col[jb];
                    value_type vb = B_nbr.val[jb];

                    if (gb >= B_beg && gb < B_end) {
                        ptrdiff_t cb = gb - B_beg;

                        if (loc_marker[cb] < loc_beg) {
                            loc_marker[cb] = loc_end;

                            C_loc.col[loc_end] = cb;
                            C_loc.val[loc_end] = va * vb;

                            ++loc_end;
                        } else {
                            C_loc.val[loc_marker[cb]] += va * vb;
                        }
                    } else {
                        ptrdiff_t cb = rem_idx[gb];

                        if (rem_marker[cb] < rem_beg) {
                            rem_marker[cb] = rem_end;

                            C_rem.col[rem_end] = gb;
                            C_rem.val[rem_end] = va * vb;

                            ++rem_end;
                        } else {
                            C_rem.val[rem_marker[cb]] += va * vb;
                        }
                    }
                }
            }
        }
    }
    AMGCL_TOC("compute");
    AMGCL_TOC("product");

    return std::make_shared<distributed_matrix<Backend> >(A.comm(), c_loc, c_rem);
}

template <class Backend, class T>
void scale(distributed_matrix<Backend> &A, T s) {
    typedef typename Backend::value_type value_type;
    typedef backend::crs<value_type> build_matrix;

    build_matrix &A_loc = *A.local();
    build_matrix &A_rem = *A.remote();

    ptrdiff_t n = A_loc.nrows;

#pragma omp parallel for
        for(ptrdiff_t i = 0; i < n; ++i) {
            for(ptrdiff_t j = A_loc.ptr[i], e = A_loc.ptr[i+1]; j < e; ++j)
                A_loc.val[j] *= s;
            for(ptrdiff_t j = A_rem.ptr[i], e = A_rem.ptr[i+1]; j < e; ++j)
                A_rem.val[j] *= s;
        }
}

template <class Backend>
void sort_rows(distributed_matrix<Backend> &A) {
    backend::sort_rows(*A.local());
    backend::sort_rows(*A.remote());
}

} // namespace mpi

namespace backend {

template <class Backend>
struct rows_impl< mpi::distributed_matrix<Backend> >
{
    static size_t get(const mpi::distributed_matrix<Backend> &A) {
        return A.loc_rows();
    }
};

template <class Backend, class Alpha, class Vec1, class Beta,  class Vec2>
struct spmv_impl<Alpha, mpi::distributed_matrix<Backend>, Vec1, Beta, Vec2>
{
    static void apply(
            Alpha alpha,
            const mpi::distributed_matrix<Backend> &A,
            const Vec1 &x, Beta beta, Vec2 &y)
    {
        A.mul(alpha, x, beta, y);
    }
};

template <class Backend, class Vec1, class Vec2, class Vec3>
struct residual_impl<mpi::distributed_matrix<Backend>, Vec1, Vec2, Vec3>
{
    static void apply(
            const Vec1 &rhs,
            const mpi::distributed_matrix<Backend> &A,
            const Vec2 &x, Vec3 &r)
    {
        A.residual(rhs, x, r);
    }
};

// Estimate spectral radius of the matrix.
template <bool scale, class Backend>
typename math::scalar_of<typename Backend::value_type>::type
spectral_radius(const mpi::distributed_matrix<Backend> &A, int power_iters = 0)
{
    AMGCL_TIC("spectral radius");
    typedef typename Backend::value_type               value_type;
    typedef typename math::rhs_of<value_type>::type    rhs_type;
    typedef typename math::scalar_of<value_type>::type scalar_type;
    typedef backend::crs<value_type>                   build_matrix;

    mpi::communicator comm = A.comm();

    const build_matrix &A_loc = *A.local();
    const build_matrix &A_rem = *A.remote();
    const mpi::comm_pattern<Backend> &C = A.cpat();

    const ptrdiff_t n = A_loc.nrows;
    scalar_type radius = 0;

    if (power_iters <= 0) {
#pragma omp parallel
        {
            scalar_type emax = 0;
            value_type  dia = math::identity<value_type>();

#pragma omp for nowait
            for(ptrdiff_t i = 0; i < n; ++i) {
                scalar_type s = 0;

                for(ptrdiff_t j = A_loc.ptr[i], e = A_loc.ptr[i+1]; j < e; ++j) {
                    ptrdiff_t  c = A_loc.col[j];
                    value_type v = A_loc.val[j];

                    s += math::norm(v);

                    if (scale && c == i) dia = v;
                }

                for(ptrdiff_t j = A_rem.ptr[i], e = A_rem.ptr[i+1]; j < e; ++j)
                    s += math::norm(A_rem.val[j]);

                if (scale) s *= math::norm(math::inverse(dia));

                emax = std::max(emax, s);
            }

#pragma omp critical
            radius = std::max(radius, emax);
        }
    } else {
        backend::numa_vector<rhs_type>   b0(n, false), b1(n, false);
        backend::numa_vector<ptrdiff_t>  rem_col(A_rem.nnz, false);

        // Fill the initial vector with random values.
        // Also extract the inverted matrix diagonal values.
        scalar_type b0_loc_norm = 0;

#pragma omp parallel
        {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
            int nt  = omp_get_max_threads();
#else
            int tid = 0;
            int nt  = 1;
#endif
            std::mt19937 rng(comm.size * nt + tid);
            std::uniform_real_distribution<scalar_type> rnd(-1, 1);

            scalar_type t_norm = 0;

#pragma omp for nowait
            for(ptrdiff_t i = 0; i < n; ++i) {
                rhs_type v = math::constant<rhs_type>(rnd(rng));

                b0[i] = v;
                t_norm += math::norm(math::inner_product(v,v));

                for(ptrdiff_t j = A_rem.ptr[i], e = A_rem.ptr[i+1]; j < e; ++j) {
                    rem_col[j] = C.local_index(A_rem.col[j]);
                }
            }

#pragma omp critical
            b0_loc_norm += t_norm;
        }

        scalar_type b0_norm = comm.reduce(MPI_SUM, b0_loc_norm);

        // Normalize b0
        b0_norm = 1 / sqrt(b0_norm);
#pragma omp parallel for
        for(ptrdiff_t i = 0; i < n; ++i) {
            b0[i] = b0_norm * b0[i];
        }

        std::vector<rhs_type> b0_send(C.send.count());
        std::vector<rhs_type> b0_recv(C.recv.count());

        for(size_t i = 0, m = C.send.count(); i < m; ++i)
            b0_send[i] = b0[C.send.col[i]];
        C.exchange(&b0_send[0], &b0_recv[0]);

        for(int iter = 0; iter < power_iters;) {
            // b1 = (D * A) * b0
            // b1_norm = ||b1||
            // radius = <b1,b0>
            scalar_type b1_loc_norm = 0;
            scalar_type loc_radius = 0;

#pragma omp parallel
            {
                scalar_type t_norm = 0;
                scalar_type t_radi = 0;
                value_type  dia = math::identity<value_type>();

#pragma omp for nowait
                for(ptrdiff_t i = 0; i < n; ++i) {
                    rhs_type s = math::zero<rhs_type>();

                    for(ptrdiff_t j = A_loc.ptr[i], e = A_loc.ptr[i+1]; j < e; ++j) {
                        ptrdiff_t  c = A_loc.col[j];
                        value_type v = A_loc.val[j];
                        if (scale && c == i) dia = v;
                        s += v * b0[c];
                    }

                    for(ptrdiff_t j = A_rem.ptr[i], e = A_rem.ptr[i+1]; j < e; ++j)
                        s += A_rem.val[j] * b0_recv[rem_col[j]];

                    if (scale) s = math::inverse(dia) * s;

                    t_norm += math::norm(math::inner_product(s, s));
                    t_radi += math::norm(math::inner_product(s, b0[i]));

                    b1[i] = s;
                }

#pragma omp critical
                {
                    b1_loc_norm += t_norm;
                    loc_radius  += t_radi;
                }
            }

            radius = comm.reduce(MPI_SUM, loc_radius);

            if (++iter < power_iters) {
                scalar_type b1_norm;
                b1_norm = comm.reduce(MPI_SUM, b1_loc_norm);

                // b0 = b1 / b1_norm
                b1_norm = 1 / sqrt(b1_norm);
#pragma omp parallel for
                for(ptrdiff_t i = 0; i < n; ++i) {
                    b0[i] = b1_norm * b1[i];
                }

                for(size_t i = 0, m = C.send.count(); i < m; ++i)
                    b0_send[i] = b0[C.send.col[i]];
                C.exchange(&b0_send[0], &b0_recv[0]);
            }
        }
    }
    AMGCL_TOC("spectral radius");

    return radius < 0 ? static_cast<scalar_type>(2) : radius;
}
} // namespace backend
} // namespace amgcl

#endif
