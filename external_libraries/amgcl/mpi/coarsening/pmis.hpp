#ifndef AMGCL_MPI_COARSENING_PMIS_HPP
#define AMGCL_MPI_COARSENING_PMIS_HPP

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
 * \file   amgcl/mpi/coarsening/pmis.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Distributed PMIS aggregation.
 */

#include <tuple>
#include <memory>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>
#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>

namespace amgcl {
namespace mpi {
namespace coarsening {

template <class Backend>
struct pmis {
    typedef typename Backend::value_type value_type;
    typedef typename math::scalar_of<value_type>::type scalar_type;
    typedef distributed_matrix<Backend> matrix;
    typedef comm_pattern<Backend> CommPattern;
    typedef backend::crs<value_type> build_matrix;
    typedef backend::builtin<char> bool_backend;
    typedef backend::crs<char>     bool_matrix;


    struct params {
        // Strong connectivity threshold
        scalar_type eps_strong;

        // Block size for non-scalar problems.
        unsigned    block_size;

        params() : eps_strong(0.08), block_size(1) { }

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, eps_strong),
              AMGCL_PARAMS_IMPORT_VALUE(p, block_size)
        {
            check_params(p, {"eps_strong", "block_size"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, eps_strong);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, block_size);
        }
#endif
    };

    std::shared_ptr< distributed_matrix<bool_backend> > conn;
    std::shared_ptr< matrix > p_tent;

    pmis(const matrix &A, const params &prm = params()) {
        ptrdiff_t n = A.loc_rows();
        std::vector<ptrdiff_t> state(n);
        std::vector<int>       owner(n);

        if (prm.block_size == 1) {
            conn = conn_strength(A, prm.eps_strong);

            ptrdiff_t naggr = aggregates(*conn, state, owner);
            p_tent = tentative_prolongation(A.comm(), n, naggr, state, owner);
        } else {
            typedef typename math::scalar_of<value_type>::type scalar;
            typedef backend::builtin<scalar> sbackend;

            ptrdiff_t np = n / prm.block_size;

            assert(np * prm.block_size == n && "Matrix size should be divisible by block_size");

            distributed_matrix<sbackend> A_pw(A.comm(),
                pointwise_matrix(*A.local(),  prm.block_size),
                pointwise_matrix(*A.remote(), prm.block_size)
                );

            auto conn_pw = conn_strength(A_pw, prm.eps_strong);

            std::vector<ptrdiff_t> state_pw(np);
            std::vector<int>       owner_pw(np);

            ptrdiff_t naggr = aggregates(*conn_pw, state_pw, owner_pw);

            conn = std::make_shared< distributed_matrix<bool_backend> >(
                    A.comm(),
                    expand_conn(*A.local(),  *A_pw.local(),  *conn_pw->local(),  prm.block_size),
                    expand_conn(*A.remote(), *A_pw.remote(), *conn_pw->remote(), prm.block_size)
                    );

#pragma omp parallel for
            for(ptrdiff_t ip = 0; ip < np; ++ip) {
                ptrdiff_t i = ip * prm.block_size;
                ptrdiff_t s = state_pw[ip];
                int       o = owner_pw[ip];

                for(unsigned k = 0; k < prm.block_size; ++k) {
                    state[i + k] = (s < 0) ? s : (s * prm.block_size + k);
                    owner[i + k] = o;
                }
            }

            p_tent = tentative_prolongation(A.comm(), n, naggr * prm.block_size, state, owner);
        }
    }

    std::shared_ptr< distributed_matrix<bool_backend> >
    squared_interface(const distributed_matrix<bool_backend> &A) {
        const comm_pattern<bool_backend> &C = A.cpat();

        bool_matrix &A_loc = *A.local();
        bool_matrix &A_rem = *A.remote();

        ptrdiff_t A_rows = A.loc_rows();

        ptrdiff_t A_beg = A.loc_col_shift();
        ptrdiff_t A_end = A_beg + A_rows;

        auto a_nbr = remote_rows(C, A, false);
        bool_matrix &A_nbr = *a_nbr;

        // Build mapping from global to local column numbers in the remote part of
        // the square matrix.
        std::vector<ptrdiff_t> rem_cols(A_rem.nnz + A_nbr.nnz);

        std::copy(A_nbr.col, A_nbr.col + A_nbr.nnz,
                std::copy(A_rem.col, A_rem.col + A_rem.nnz, rem_cols.begin()));

        std::sort(rem_cols.begin(), rem_cols.end());
        rem_cols.erase(std::unique(rem_cols.begin(), rem_cols.end()), rem_cols.end());

        ptrdiff_t n_rem_cols = 0;
        std::unordered_map<ptrdiff_t, int> rem_idx(2 * rem_cols.size());
        for(ptrdiff_t c : rem_cols) {
            if (c >= A_beg && c < A_end) continue;
            rem_idx[c] = n_rem_cols++;
        }

        // Build the product.
        auto s_loc = std::make_shared<bool_matrix>();
        auto s_rem = std::make_shared<bool_matrix>();

        bool_matrix &S_loc = *s_loc;
        bool_matrix &S_rem = *s_rem;

        S_loc.set_size(A_rows, A_rows, false);
        S_rem.set_size(A_rows, 0,      false);

        S_loc.ptr[0] = 0;
        S_rem.ptr[0] = 0;

        AMGCL_TIC("analyze");
#pragma omp parallel
        {
            std::vector<ptrdiff_t> loc_marker(A_rows,     -1);
            std::vector<ptrdiff_t> rem_marker(n_rem_cols, -1);

#pragma omp for
            for(ptrdiff_t ia = 0; ia < A_rows; ++ia) {
                ptrdiff_t loc_cols = 0;
                ptrdiff_t rem_cols = 0;

                for(ptrdiff_t ja = A_rem.ptr[ia], ea = A_rem.ptr[ia + 1]; ja < ea; ++ja) {
                    ptrdiff_t  ca = C.local_index(A_rem.col[ja]);

                    for(ptrdiff_t jb = A_nbr.ptr[ca], eb = A_nbr.ptr[ca+1]; jb < eb; ++jb) {
                        ptrdiff_t  cb = A_nbr.col[jb];

                        if (cb >= A_beg && cb < A_end) {
                            cb -= A_beg;

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

                for(ptrdiff_t ja = A_loc.ptr[ia], ea = A_loc.ptr[ia + 1]; ja < ea; ++ja) {
                    ptrdiff_t  ca = A_loc.col[ja];

                    for(ptrdiff_t jb = A_rem.ptr[ca], eb = A_rem.ptr[ca+1]; jb < eb; ++jb) {
                        ptrdiff_t  cb = rem_idx[A_rem.col[jb]];

                        if (rem_marker[cb] != ia) {
                            rem_marker[cb]  = ia;
                            ++rem_cols;
                        }
                    }

                }

                if (rem_cols) {
                    for(ptrdiff_t ja = A_loc.ptr[ia], ea = A_loc.ptr[ia + 1]; ja < ea; ++ja) {
                        ptrdiff_t  ca = A_loc.col[ja];

                        for(ptrdiff_t jb = A_loc.ptr[ca], eb = A_loc.ptr[ca+1]; jb < eb; ++jb) {
                            ptrdiff_t  cb = A_loc.col[jb];

                            if (loc_marker[cb] != ia) {
                                loc_marker[cb]  = ia;
                                ++loc_cols;
                            }
                        }

                    }
                }

                S_rem.ptr[ia + 1] = rem_cols;
                S_loc.ptr[ia + 1] = rem_cols ? loc_cols : 0;
            }
        }
        AMGCL_TOC("analyze");

        S_loc.set_nonzeros(S_loc.scan_row_sizes(), false);
        S_rem.set_nonzeros(S_rem.scan_row_sizes(), false);

        AMGCL_TIC("compute");
#pragma omp parallel
        {
            std::vector<ptrdiff_t> loc_marker(A_rows,     -1);
            std::vector<ptrdiff_t> rem_marker(n_rem_cols, -1);

#pragma omp for
            for(ptrdiff_t ia = 0; ia < A_rows; ++ia) {
                ptrdiff_t loc_beg = S_loc.ptr[ia];
                ptrdiff_t rem_beg = S_rem.ptr[ia];
                ptrdiff_t loc_end = loc_beg;
                ptrdiff_t rem_end = rem_beg;

                if (rem_beg == S_rem.ptr[ia+1]) continue;

                for(ptrdiff_t ja = A_loc.ptr[ia], ea = A_loc.ptr[ia + 1]; ja < ea; ++ja) {
                    ptrdiff_t  ca = A_loc.col[ja];

                    for(ptrdiff_t jb = A_loc.ptr[ca], eb = A_loc.ptr[ca+1]; jb < eb; ++jb) {
                        ptrdiff_t  cb = A_loc.col[jb];

                        if (loc_marker[cb] < loc_beg) {
                            loc_marker[cb] = loc_end;
                            S_loc.col[loc_end] = cb;
                            ++loc_end;
                        }
                    }

                    for(ptrdiff_t jb = A_rem.ptr[ca], eb = A_rem.ptr[ca+1]; jb < eb; ++jb) {
                        ptrdiff_t  gb = A_rem.col[jb];
                        ptrdiff_t  cb = rem_idx[gb];

                        if (rem_marker[cb] < rem_beg) {
                            rem_marker[cb] = rem_end;
                            S_rem.col[rem_end] = gb;
                            ++rem_end;
                        }
                    }
                }

                for(ptrdiff_t ja = A_rem.ptr[ia], ea = A_rem.ptr[ia + 1]; ja < ea; ++ja) {
                    ptrdiff_t  ca = C.local_index(A_rem.col[ja]);

                    for(ptrdiff_t jb = A_nbr.ptr[ca], eb = A_nbr.ptr[ca+1]; jb < eb; ++jb) {
                        ptrdiff_t  gb = A_nbr.col[jb];

                        if (gb >= A_beg && gb < A_end) {
                            ptrdiff_t cb = gb - A_beg;

                            if (loc_marker[cb] < loc_beg) {
                                loc_marker[cb] = loc_end;
                                S_loc.col[loc_end] = cb;
                                ++loc_end;
                            }
                        } else {
                            ptrdiff_t cb = rem_idx[gb];

                            if (rem_marker[cb] < rem_beg) {
                                rem_marker[cb] = rem_end;
                                S_rem.col[rem_end] = gb;
                                ++rem_end;
                            }
                        }
                    }
                }
            }
        }
        AMGCL_TOC("compute");

        return std::make_shared< distributed_matrix<bool_backend> >(A.comm(), s_loc, s_rem);
    }

    template <class B>
    std::shared_ptr< distributed_matrix<bool_backend> >
    conn_strength(const distributed_matrix<B> &A, scalar_type eps_strong) {
        typedef typename B::value_type val_type;
        typedef backend::crs<val_type> B_matrix;

        AMGCL_TIC("conn_strength");
        ptrdiff_t n = A.loc_rows();

        const B_matrix &A_loc = *A.local();
        const B_matrix &A_rem = *A.remote();
        const comm_pattern<B> &C = A.cpat();

        scalar_type eps_squared = eps_strong * eps_strong;

        auto d = backend::diagonal(A_loc);
        backend::numa_vector<val_type> &D = *d;

        std::vector<val_type> D_loc(C.send.count());
        std::vector<val_type> D_rem(C.recv.count());

        for(size_t i = 0, nv = C.send.count(); i < nv; ++i)
            D_loc[i] = D[C.send.col[i]];

        C.exchange(&D_loc[0], &D_rem[0]);

        auto s_loc = std::make_shared<bool_matrix>();
        auto s_rem = std::make_shared<bool_matrix>();

        bool_matrix &S_loc = *s_loc;
        bool_matrix &S_rem = *s_rem;

        S_loc.set_size(n, n, true);
        S_rem.set_size(n, 0, true);

        S_loc.val = new char[A_loc.nnz];
        S_rem.val = new char[A_rem.nnz];

#pragma omp parallel for
        for(ptrdiff_t i = 0; i < n; ++i) {
            val_type eps_dia_i = eps_squared * D[i];

            for(ptrdiff_t j = A_loc.ptr[i], e = A_loc.ptr[i+1]; j < e; ++j) {
                ptrdiff_t c = A_loc.col[j];
                val_type  v = A_loc.val[j];

                if ((S_loc.val[j] = (c == i || (eps_dia_i * D[c] < v * v))))
                    ++S_loc.ptr[i + 1];
            }

            for(ptrdiff_t j = A_rem.ptr[i], e = A_rem.ptr[i+1]; j < e; ++j) {
                ptrdiff_t c = C.local_index(A_rem.col[j]);
                val_type  v = A_rem.val[j];

                if ((S_rem.val[j] = (eps_dia_i * D_rem[c] < v * v)))
                    ++S_rem.ptr[i + 1];
            }
        }

        S_loc.nnz = S_loc.scan_row_sizes();
        S_rem.nnz = S_rem.scan_row_sizes();

        S_loc.col = new ptrdiff_t[S_loc.nnz];
        S_rem.col = new ptrdiff_t[S_rem.nnz];

#pragma omp parallel for
        for(ptrdiff_t i = 0; i < n; ++i) {
            ptrdiff_t loc_head = S_loc.ptr[i];
            ptrdiff_t rem_head = S_rem.ptr[i];

            for(ptrdiff_t j = A_loc.ptr[i], e = A_loc.ptr[i+1]; j < e; ++j)
                if (S_loc.val[j]) S_loc.col[loc_head++] = A_loc.col[j];

            for(ptrdiff_t j = A_rem.ptr[i], e = A_rem.ptr[i+1]; j < e; ++j)
                if (S_rem.val[j]) S_rem.col[rem_head++] = A_rem.col[j];
        }
        AMGCL_TOC("conn_strength");

        return std::make_shared< distributed_matrix<bool_backend> >(
                A.comm(), s_loc, s_rem);
    }

    ptrdiff_t aggregates(
            const distributed_matrix<bool_backend> &A,
            std::vector<ptrdiff_t> &loc_state,
            std::vector<int>       &loc_owner
            )
    {
        AMGCL_TIC("PMIS");
        static const int tag_exc_cnt = 4001;
        static const int tag_exc_pts = 4002;

        const bool_matrix &A_loc = *A.local();
        const bool_matrix &A_rem = *A.remote();

        ptrdiff_t n = A_loc.nrows;

        communicator comm = A.comm();

        // 1. Get symbolic square of the connectivity matrix.
        AMGCL_TIC("symbolic square");
        auto S = squared_interface(A);
        const bool_matrix &S_loc = *S->local();
        const bool_matrix &S_rem = *S->remote();
        const comm_pattern<bool_backend> &Sp = S->cpat();
        AMGCL_TOC("symbolic square");

        // 2. Apply PMIS algorithm to the symbolic square.
        ptrdiff_t n_undone = 0;
        std::vector<ptrdiff_t> rem_state(Sp.recv.count(), undone);
        std::vector<ptrdiff_t> send_state(Sp.send.count());

        // Remove lonely nodes.
#pragma omp parallel for reduction(+:n_undone)
        for(ptrdiff_t i = 0; i < n; ++i) {
            ptrdiff_t wl = A_loc.ptr[i+1] - A_loc.ptr[i];
            ptrdiff_t wr = S_rem.ptr[i+1] - S_rem.ptr[i];

            if (wl + wr == 1) {
                loc_state[i] = deleted;
                ++n_undone;
            } else {
                loc_state[i] = undone;
            }

            loc_owner[i] = -1;
        }

        n_undone = n - n_undone;

        // Exchange state
        for(ptrdiff_t i = 0, m = Sp.send.count(); i < m; ++i)
            send_state[i] = loc_state[Sp.send.col[i]];
        Sp.exchange(&send_state[0], &rem_state[0]);

        std::vector< std::vector<ptrdiff_t> > send_pts(Sp.recv.nbr.size());
        std::vector<ptrdiff_t> recv_pts;

        std::vector<MPI_Request> send_cnt_req(Sp.recv.nbr.size());
        std::vector<MPI_Request> send_pts_req(Sp.recv.nbr.size());

        ptrdiff_t naggr = 0;

        std::vector<ptrdiff_t> nbr;

        while(true) {
            for(size_t i = 0; i < Sp.recv.nbr.size(); ++i)
                send_pts[i].clear();

            if (n_undone) {
                for(ptrdiff_t i = 0; i < n; ++i) {
                    if (loc_state[i] != undone) continue;

                    if (S_rem.ptr[i+1] > S_rem.ptr[i]) {
                        // Boundary points
                        bool selectable = true;
                        for(ptrdiff_t j = S_rem.ptr[i], e = S_rem.ptr[i+1]; j < e; ++j) {
                            int d,c;
                            std::tie(d,c) = Sp.remote_info(S_rem.col[j]);

                            if (rem_state[c] == undone && Sp.recv.nbr[d] > comm.rank) {
                                selectable = false;
                                break;
                            }
                        }

                        if (!selectable) continue;

                        ptrdiff_t id = naggr++;
                        loc_owner[i] = comm.rank;
                        loc_state[i] = id;
                        --n_undone;

                        // A gives immediate neighbors
                        for(ptrdiff_t j = A_loc.ptr[i], e = A_loc.ptr[i+1]; j < e; ++j) {
                            ptrdiff_t c = A_loc.col[j];
                            if (c != i) {
                                if (loc_state[c] == undone) --n_undone;
                                loc_owner[c] = comm.rank;
                                loc_state[c] = id;
                            }
                        }

                        for(ptrdiff_t j = A_rem.ptr[i], e = A_rem.ptr[i+1]; j < e; ++j) {
                            ptrdiff_t c = A_rem.col[j];
                            int d,k;
                            std::tie(d,k) = Sp.remote_info(c);

                            rem_state[k] = id;
                            send_pts[d].push_back(c);
                            send_pts[d].push_back(id);
                        }

                        // S gives removed neighbors
                        for(ptrdiff_t j = S_loc.ptr[i], e = S_loc.ptr[i+1]; j < e; ++j) {
                            ptrdiff_t c = S_loc.col[j];
                            if (c != i && loc_state[c] == undone) {
                                loc_owner[c] = comm.rank;
                                loc_state[c] = id;
                                --n_undone;
                            }
                        }

                        for(ptrdiff_t j = S_rem.ptr[i], e = S_rem.ptr[i+1]; j < e; ++j) {
                            ptrdiff_t c = S_rem.col[j];
                            int d,k;
                            std::tie(d,k) = Sp.remote_info(c);

                            if (rem_state[k] == undone) {
                                rem_state[k] = id;
                                send_pts[d].push_back(c);
                                send_pts[d].push_back(id);
                            }
                        }
                    } else {
                        // Inner points
                        ptrdiff_t id = naggr++;
                        loc_owner[i] = comm.rank;
                        loc_state[i] = id;
                        --n_undone;

                        nbr.clear();

                        for(ptrdiff_t j = A_loc.ptr[i], e = A_loc.ptr[i+1]; j < e; ++j) {
                            ptrdiff_t c = A_loc.col[j];

                            if (c != i && loc_state[c] != deleted) {
                                if (loc_state[c] == undone) --n_undone;
                                loc_owner[c] = comm.rank;
                                loc_state[c] = id;
                                nbr.push_back(c);
                            }
                        }

                        for(ptrdiff_t k : nbr) {
                            for(ptrdiff_t j = A_loc.ptr[k], e = A_loc.ptr[k+1]; j < e; ++j) {
                                ptrdiff_t c = A_loc.col[j];
                                if (c != k && loc_state[c] == undone) {
                                    loc_owner[c] = comm.rank;
                                    loc_state[c] = id;
                                    --n_undone;
                                }
                            }
                        }
                    }
                }
            }

            for(size_t i = 0; i < Sp.recv.nbr.size(); ++i) {
                int npts = send_pts[i].size();
                MPI_Isend(&npts, 1, MPI_INT, Sp.recv.nbr[i], tag_exc_cnt, comm, &send_cnt_req[i]);

                if (!npts) continue;
                MPI_Isend(&send_pts[i][0], npts, datatype<ptrdiff_t>(), Sp.recv.nbr[i], tag_exc_pts, comm, &send_pts_req[i]);
            }

            for(size_t i = 0; i < Sp.send.nbr.size(); ++i) {
                int npts;
                MPI_Recv(&npts, 1, MPI_INT, Sp.send.nbr[i], tag_exc_cnt, comm, MPI_STATUS_IGNORE);

                if (!npts) continue;
                recv_pts.resize(npts);
                MPI_Recv(&recv_pts[0], npts, datatype<ptrdiff_t>(), Sp.send.nbr[i], tag_exc_pts, comm, MPI_STATUS_IGNORE);

                for(int k = 0; k < npts; k += 2) {
                    ptrdiff_t c  = recv_pts[k] - Sp.loc_col_shift();
                    ptrdiff_t id = recv_pts[k+1];

                    if (loc_state[c] == undone) --n_undone;

                    loc_owner[c] = Sp.send.nbr[i];
                    loc_state[c] = id;
                }
            }

            for(size_t i = 0; i < Sp.recv.nbr.size(); ++i) {
                int npts = send_pts[i].size();
                MPI_Wait(&send_cnt_req[i], MPI_STATUS_IGNORE);
                if (!npts) continue;
                MPI_Wait(&send_pts_req[i], MPI_STATUS_IGNORE);
            }


            for(ptrdiff_t i = 0, m = Sp.send.count(); i < m; ++i)
                send_state[i] = loc_state[Sp.send.col[i]];
            Sp.exchange(&send_state[0], &rem_state[0]);

            if (0 == comm.reduce(MPI_SUM, n_undone))
                break;
        }
        AMGCL_TOC("PMIS");

        return naggr;
    }

    std::shared_ptr<matrix>
    tentative_prolongation(communicator comm, ptrdiff_t n, ptrdiff_t naggr,
            std::vector<ptrdiff_t> &state, std::vector<int> &owner)
    {
        AMGCL_TIC("tentative prolongation");
        // Form tentative prolongation operator.
        auto p_loc = std::make_shared<build_matrix>();
        auto p_rem = std::make_shared<build_matrix>();
        build_matrix &P_loc = *p_loc;
        build_matrix &P_rem = *p_rem;

        std::vector<ptrdiff_t> dom = comm.exclusive_sum(naggr);
        P_loc.set_size(n, naggr, true);
        P_rem.set_size(n, 0, true);

#pragma omp parallel for
        for(ptrdiff_t i = 0; i < n; ++i) {
            if (state[i] == deleted) continue;

            if (owner[i] == comm.rank) {
                ++P_loc.ptr[i+1];
            } else {
                ++P_rem.ptr[i+1];
            }
        }

        P_loc.set_nonzeros(P_loc.scan_row_sizes());
        P_rem.set_nonzeros(P_rem.scan_row_sizes());

#pragma omp parallel for
        for(ptrdiff_t i = 0; i < n; ++i) {
            ptrdiff_t s = state[i];
            if (s == deleted) continue;

            int d = owner[i];
            if (d == comm.rank) {
                P_loc.col[P_loc.ptr[i]] = s;
                P_loc.val[P_loc.ptr[i]] = math::identity<value_type>();
            } else {
                P_rem.col[P_rem.ptr[i]] = s + dom[d];
                P_rem.val[P_rem.ptr[i]] = math::identity<value_type>();
            }
        }
        AMGCL_TOC("tentative prolongation");

        return std::make_shared<matrix>(comm, p_loc, p_rem);
    }

    template <class pw_matrix>
    std::shared_ptr<bool_matrix>
    expand_conn(const build_matrix &A, const pw_matrix &Ap, const bool_matrix &Cp,
            unsigned block_size) const
    {
        ptrdiff_t np = Cp.nrows;
        ptrdiff_t n  = np * block_size;

        auto c = std::make_shared<bool_matrix>();
        bool_matrix &C = *c;

        C.set_size(n, n, true);
        C.val = new char[A.nnz];

#pragma omp parallel
        {
            std::vector<ptrdiff_t> j(block_size);
            std::vector<ptrdiff_t> e(block_size);

#pragma omp for
            for(ptrdiff_t ip = 0; ip < np; ++ip) {
                ptrdiff_t ia = ip * block_size;

                for(unsigned k = 0; k < block_size; ++k) {
                    j[k] = A.ptr[ia + k];
                    e[k] = A.ptr[ia + k + 1];
                }

                for(ptrdiff_t jp = Ap.ptr[ip], ep = Ap.ptr[ip+1]; jp < ep; ++jp) {
                    ptrdiff_t cp = Ap.col[jp];
                    bool      sp = Cp.val[jp];

                    ptrdiff_t col_end = (cp + 1) * block_size;

                    for(unsigned k = 0; k < block_size; ++k) {
                        ptrdiff_t beg = j[k];
                        ptrdiff_t end = e[k];

                        while(beg < end && A.col[beg] < col_end) {
                            C.val[beg++] = sp;

                            if (sp) ++C.ptr[ia + k + 1];
                        }

                        j[k] = beg;
                    }
                }
            }
        }

        C.nnz = C.scan_row_sizes();
        C.col = new ptrdiff_t[C.nnz];

#pragma omp parallel
        {
            std::vector<ptrdiff_t> j(block_size);
            std::vector<ptrdiff_t> e(block_size);
            std::vector<ptrdiff_t> h(block_size);

#pragma omp for
            for(ptrdiff_t ip = 0; ip < np; ++ip) {
                ptrdiff_t ia = ip * block_size;

                for(unsigned k = 0; k < block_size; ++k) {
                    j[k] = A.ptr[ia + k];
                    e[k] = A.ptr[ia + k + 1];
                    h[k] = C.ptr[ia + k];
                }

                for(ptrdiff_t jp = Ap.ptr[ip], ep = Ap.ptr[ip+1]; jp < ep; ++jp) {
                    ptrdiff_t cp = Ap.col[jp];
                    bool      sp = Cp.val[jp];

                    ptrdiff_t col_end = (cp + 1) * block_size;

                    for(unsigned k = 0; k < block_size; ++k) {
                        ptrdiff_t beg = j[k];
                        ptrdiff_t end = e[k];
                        ptrdiff_t hed = h[k];

                        while(beg < end && A.col[beg] < col_end) {
                            if (sp) C.col[hed++] = A.col[beg];
                            ++beg;
                        }

                        j[k] = beg;
                        h[k] = hed;
                    }
                }
            }
        }

        return c;
    }

    private:
        static const int undone = -2;
        static const int deleted = -1;
};

} // namespace coarsening
} // namespace mpi
} // namespace amgcl

#endif
