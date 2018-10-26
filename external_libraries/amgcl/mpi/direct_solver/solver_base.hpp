#ifndef AMGCL_MPI_DIRECT_SOLVER_SOLVER_BASE_HPP
#define AMGCL_MPI_DIRECT_SOLVER_SOLVER_BASE_HPP

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

/**
\file   amgcl/mpi/direct_solver/solver_base.hpp
\author Denis Demidov <dennis.demidov@gmail.com>
\brief  Basic functionality for distributed direct solvers.
*/

#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>

namespace amgcl {
namespace mpi {
namespace direct {

template <class value_type, class Solver>
class solver_base {
    public:
        typedef typename math::scalar_of<value_type>::type scalar_type;
        typedef typename math::rhs_of<value_type>::type    rhs_type;
        typedef backend::crs<value_type> build_matrix;

        solver_base() {}

        void init(communicator comm, const build_matrix &Astrip) {
            this->comm = comm;
            n = Astrip.nrows;

            std::vector<int> domain = comm.exclusive_sum(n);
            std::vector<int> active; active.reserve(comm.size);

            // Find out how many ranks are active (own non-zero matrix rows):
            int active_rank = 0;
            for(int i = 0; i < comm.size; ++i) {
                if (domain[i+1] - domain[i] > 0) {
                    if (comm.rank == i) active_rank = active.size();
                    active.push_back(i);
                }
            }

            // Consolidate the matrix on a fewer processes.
            int nmasters = std::min<int>(active.size(), solver().comm_size(domain.back()));
            int slaves_per_master = (active.size() + nmasters - 1) / nmasters;
            int group_beg = (active_rank / slaves_per_master) * slaves_per_master;

            group_master = active[group_beg];

            // Communicator for masters (used to solve the coarse problem):
            MPI_Comm_split(comm,
                    comm.rank == group_master ? 0 : MPI_UNDEFINED,
                    comm.rank, &masters_comm
                    );

            if (!n) return; // I am not active

            // Shift from row pointers to row widths:
            std::vector<ptrdiff_t> widths(n);
            for(ptrdiff_t i = 0; i < n; ++i)
                widths[i] = Astrip.ptr[i+1] - Astrip.ptr[i];

            if (comm.rank == group_master) {
                int group_end = std::min<int>(group_beg + slaves_per_master, active.size());
                group_beg += 1;
                int group_size = group_end - group_beg;

                std::vector<MPI_Request> cnt_req(group_size);
                std::vector<MPI_Request> col_req(group_size);
                std::vector<MPI_Request> val_req(group_size);

                solve_req.resize(group_size);
                slaves.reserve(group_size);
                counts.reserve(group_size);

                // Count rows in local chunk of the consolidated matrix,
                // see who is reporting to us.
                int nloc = n;
                for(int j = group_beg; j < group_end; ++j) {
                    int i = active[j];

                    int m = domain[i+1] - domain[i];
                    nloc += m;
                    counts.push_back(m);
                    slaves.push_back(i);
                }

                // Get matrix chunks from my slaves.
                build_matrix A;
                A.set_size(nloc, domain.back(), false);
                A.ptr[0] = 0;

                cons_f.resize(A.nrows);
                cons_x.resize(A.nrows);

                int shift = n+1;
                std::copy(widths.begin(), widths.end(), &A.ptr[1]);

                for(int j = 0; j < group_size; ++j) {
                    int i = slaves[j];

                    MPI_Irecv(&A.ptr[shift], counts[j], datatype<ptrdiff_t>(),
                            i, cnt_tag, comm, &cnt_req[j]);

                    shift += counts[j];
                }

                MPI_Waitall(cnt_req.size(), &cnt_req[0], MPI_STATUSES_IGNORE);

                A.set_nonzeros(A.scan_row_sizes());

                std::copy(Astrip.col, Astrip.col + Astrip.nnz, A.col);
                std::copy(Astrip.val, Astrip.val + Astrip.nnz, A.val);

                shift = Astrip.nnz;
                for(int j = 0, d0 = domain[comm.rank]; j < group_size; ++j) {
                    int i = slaves[j];

                    int nnz = A.ptr[domain[i+1] - d0] - A.ptr[domain[i] - d0];

                    MPI_Irecv(A.col + shift, nnz, datatype<ptrdiff_t>(),
                            i, col_tag, comm, &col_req[j]);

                    MPI_Irecv(A.val + shift, nnz, datatype<value_type>(),
                            i, val_tag, comm, &val_req[j]);

                    shift += nnz;
                }

                MPI_Waitall(col_req.size(), &col_req[0], MPI_STATUSES_IGNORE);
                MPI_Waitall(val_req.size(), &val_req[0], MPI_STATUSES_IGNORE);

                solver().init(masters_comm, A);
            } else {
                MPI_Send(widths.data(), n, datatype<ptrdiff_t>(),
                        group_master, cnt_tag, comm);
                MPI_Send(Astrip.col, Astrip.nnz, datatype<ptrdiff_t>(),
                        group_master, col_tag, comm);
                MPI_Send(Astrip.val, Astrip.nnz, datatype<value_type>(),
                        group_master, val_tag, comm);
            }

            host_v.resize(n);
        }

        template <class B>
        void init(communicator comm, const distributed_matrix<B> &A) {
            const build_matrix &A_loc = *A.local();
            const build_matrix &A_rem = *A.remote();

            build_matrix a;

            a.set_size(A.loc_rows(), A.glob_cols(), false);
            a.set_nonzeros(A_loc.nnz + A_rem.nnz);
            a.ptr[0] = 0;

            for(size_t i = 0, head = 0; i < A_loc.nrows; ++i) {
                ptrdiff_t shift = A.loc_col_shift();

                for(ptrdiff_t j = A_loc.ptr[i], e = A_loc.ptr[i+1]; j < e; ++j) {
                    a.col[head] = A_loc.col[j] + shift;
                    a.val[head] = A_loc.val[j];
                    ++head;
                }

                for(ptrdiff_t j = A_rem.ptr[i], e = A_rem.ptr[i+1]; j < e; ++j) {
                    a.col[head] = A_rem.col[j];
                    a.val[head] = A_rem.val[j];
                    ++head;
                }

                a.ptr[i+1] = head;
            }

            init(comm, a);
        }

        virtual ~solver_base() {
            if (masters_comm != MPI_COMM_NULL) MPI_Comm_free(&masters_comm);
        }

        Solver& solver() {
            return *static_cast<Solver*>(this);
        }

        const Solver& solver() const {
            return *static_cast<const Solver*>(this);
        }

        template <class VecF, class VecX>
        void operator()(const VecF &f, VecX &x) const {
            static const MPI_Datatype T = datatype<rhs_type>();

            if (!n) return;

            backend::copy(f, host_v);

            if (comm.rank == group_master) {
                std::copy(host_v.begin(), host_v.end(), cons_f.begin());

                int shift = n, j = 0;
                for(int i : slaves) {
                    MPI_Irecv(&cons_f[shift], counts[j], T, i, rhs_tag, comm, &solve_req[j]);
                    shift += counts[j++];
                }

                MPI_Waitall(solve_req.size(), &solve_req[0], MPI_STATUSES_IGNORE);

                solver().solve(cons_f, cons_x);

                std::copy(cons_x.begin(), cons_x.begin() + n, host_v.begin());
                shift = n;
                j = 0;

                for(int i : slaves) {
                    MPI_Isend(&cons_x[shift], counts[j], T, i, sol_tag, comm, &solve_req[j]);
                    shift += counts[j++];
                }

                MPI_Waitall(solve_req.size(), &solve_req[0], MPI_STATUSES_IGNORE);
            } else {
                MPI_Send(&host_v[0], n, T, group_master, rhs_tag, comm);
                MPI_Recv(&host_v[0], n, T, group_master, sol_tag, comm, MPI_STATUS_IGNORE);
            }

            backend::copy(host_v, x);
        }
    private:
        static const int cnt_tag = 5001;
        static const int col_tag = 5002;
        static const int val_tag = 5003;
        static const int rhs_tag = 5004;
        static const int sol_tag = 5005;

        communicator comm;
        int          n;
        int          group_master;
        MPI_Comm     masters_comm;
        std::vector<int> slaves;
        std::vector<int> counts;
        mutable std::vector<rhs_type> cons_f, cons_x, host_v;
        mutable std::vector<MPI_Request> solve_req;
};

} // namespace direct
} // namespace mpi
} // namespace amgcl

#endif
