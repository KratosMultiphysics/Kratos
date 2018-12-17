#ifndef AMGCL_MPI_REPARTITION_UTIL_HPP
#define AMGCL_MPI_REPARTITION_UTIL_HPP

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
 * \file   amgcl/mpi/partition/util.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Repartitioning utils
 */

#include <vector>
#include <algorithm>
#include <numeric>

#include <tuple>

#include <amgcl/backend/interface.hpp>
#include <amgcl/value_type/interface.hpp>
#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>

namespace amgcl {
namespace mpi {
namespace partition {

template <class Backend, class Ptr, class Col>
void symm_graph(const distributed_matrix<Backend> &A,
        std::vector<Ptr> &ptr, std::vector<Col> &col)
{
    typedef typename Backend::value_type value_type;
    typedef backend::crs<value_type> build_matrix;

    AMGCL_TIC("symm graph");

    build_matrix &A_loc = *A.local();
    build_matrix &A_rem = *A.remote();

    ptrdiff_t n = A_loc.nrows;
    ptrdiff_t row_beg = A.loc_col_shift();

    auto T = transpose(A);

    build_matrix &T_loc = *T->local();
    build_matrix &T_rem = *T->remote();

    // Build symmetric graph
    ptr.resize(n + 1, 0);

#pragma omp parallel for
    for(ptrdiff_t i = 0; i < n; ++i) {
        using amgcl::detail::sort_row;

        ptrdiff_t A_loc_beg = A_loc.ptr[i];
        ptrdiff_t A_loc_end = A_loc.ptr[i+1];

        ptrdiff_t A_rem_beg = A_rem.ptr[i];
        ptrdiff_t A_rem_end = A_rem.ptr[i+1];

        ptrdiff_t T_loc_beg = T_loc.ptr[i];
        ptrdiff_t T_loc_end = T_loc.ptr[i+1];

        ptrdiff_t T_rem_beg = T_rem.ptr[i];
        ptrdiff_t T_rem_end = T_rem.ptr[i+1];

        sort_row(A_loc.col + A_loc_beg, A_loc.val + A_loc_beg, A_loc_end - A_loc_beg);
        sort_row(A_rem.col + A_rem_beg, A_rem.val + A_rem_beg, A_rem_end - A_rem_beg);

        sort_row(T_loc.col + T_loc_beg, T_loc.val + T_loc_beg, T_loc_end - T_loc_beg);
        sort_row(T_rem.col + T_rem_beg, T_rem.val + T_rem_beg, T_rem_end - T_rem_beg);

        Ptr row_width = 0;

        for(ptrdiff_t ja = A_loc_beg, jt = T_loc_beg; ja < A_loc_end || jt < T_loc_end;) {
            ptrdiff_t c;
            if (ja == A_loc_end) {
                c = T_loc.col[jt];
                ++jt;
            } else if (jt == T_loc_end) {
                c = A_loc.col[ja];
                ++ja;
            } else {
                ptrdiff_t ca = A_loc.col[ja];
                ptrdiff_t ct = T_loc.col[jt];
                if (ca < ct) {
                    c = ca;
                    ++ja;
                } else if (ca == ct) {
                    c = ca;
                    ++ja;
                    ++jt;
                } else {
                    c = ct;
                    ++jt;
                }
            }

            if (c != i) ++row_width;
        }

        for(ptrdiff_t ja = A_rem_beg, jt = T_rem_beg; ja < A_rem_end || jt < T_rem_end;) {
            if (ja == A_rem_end) {
                ++jt;
            } else if (jt == T_rem_end) {
                ++ja;
            } else {
                ptrdiff_t ca = A_rem.col[ja];
                ptrdiff_t ct = T_rem.col[jt];
                if (ca < ct) {
                    ++ja;
                } else if (ca == ct) {
                    ++ja;
                    ++jt;
                } else {
                    ++jt;
                }
            }

            ++row_width;
        }

        ptr[i+1] = row_width;
    }

    std::partial_sum(ptr.begin(), ptr.end(), ptr.begin());

    col.resize(ptr.back());
    if (col.empty()) col.reserve(1); // So that col.data() is not NULL

#pragma omp parallel for
    for(ptrdiff_t i = 0; i < n; ++i) {
        ptrdiff_t A_loc_beg = A_loc.ptr[i];
        ptrdiff_t A_loc_end = A_loc.ptr[i+1];

        ptrdiff_t A_rem_beg = A_rem.ptr[i];
        ptrdiff_t A_rem_end = A_rem.ptr[i+1];

        ptrdiff_t T_loc_beg = T_loc.ptr[i];
        ptrdiff_t T_loc_end = T_loc.ptr[i+1];

        ptrdiff_t T_rem_beg = T_rem.ptr[i];
        ptrdiff_t T_rem_end = T_rem.ptr[i+1];

        Ptr head = ptr[i];

        for(ptrdiff_t ja = A_loc_beg, jt = T_loc_beg; ja < A_loc_end || jt < T_loc_end;) {
            ptrdiff_t c;
            if (ja == A_loc_end) {
                c = T_loc.col[jt];
                ++jt;
            } else if (jt == T_loc_end) {
                c = A_loc.col[ja];
                ++ja;
            } else {
                ptrdiff_t ca = A_loc.col[ja];
                ptrdiff_t ct = T_loc.col[jt];

                if (ca < ct) {
                    c = ca;
                    ++ja;
                } else if (ca == ct) {
                    c = ca;
                    ++ja;
                    ++jt;
                } else {
                    c = ct;
                    ++jt;
                }
            }
            if (c != i) col[head++] = c + row_beg;
        }

        for(ptrdiff_t ja = A_rem_beg, jt = T_rem_beg; ja < A_rem_end || jt < T_rem_end;) {
            if (ja == A_rem_end) {
                col[head] = T_rem.col[jt];
                ++jt;
            } else if (jt == T_rem_end) {
                col[head] = A_rem.col[ja];
                ++ja;
            } else {
                ptrdiff_t ca = A_rem.col[ja];
                ptrdiff_t ct = T_rem.col[jt];

                if (ca < ct) {
                    col[head] = ca;
                    ++ja;
                } else if (ca == ct) {
                    col[head] = ca;
                    ++ja;
                    ++jt;
                } else {
                    col[head] = ct;
                    ++jt;
                }
            }
            ++head;
        }
    }

    AMGCL_TOC("symm graph");
}

template <class Idx>
std::tuple<ptrdiff_t, ptrdiff_t> graph_perm_index(
        communicator comm, int npart, const std::vector<Idx> &part,
        std::vector<ptrdiff_t> &perm)
{
    AMGCL_TIC("perm index");
    ptrdiff_t n = part.size();
    perm.resize(n);

    std::vector<ptrdiff_t> loc_part_cnt(npart, 0);
    std::vector<ptrdiff_t> loc_part_beg(npart, 0);
    std::vector<ptrdiff_t> glo_part_cnt(npart);
    std::vector<ptrdiff_t> glo_part_beg(npart + 1);

    for(Idx p : part) ++loc_part_cnt[p];

    MPI_Exscan(&loc_part_cnt[0], &loc_part_beg[0], npart, datatype<ptrdiff_t>(), MPI_SUM, comm);
    MPI_Allreduce(&loc_part_cnt[0], &glo_part_cnt[0], npart, datatype<ptrdiff_t>(), MPI_SUM, comm);

    glo_part_beg[0] = 0;
    std::partial_sum(glo_part_cnt.begin(), glo_part_cnt.end(), glo_part_beg.begin() + 1);

    std::vector<ptrdiff_t> cnt(npart, 0);
    for(ptrdiff_t i = 0; i < n; ++i) {
        Idx p = part[i];
        perm[i] = glo_part_beg[p] + loc_part_beg[p] + cnt[p]++;
    }

    AMGCL_TOC("perm index");
    return std::make_tuple(
            glo_part_beg[std::min(npart, comm.rank)],
            glo_part_beg[std::min(npart, comm.rank + 1)]
            );
}

template <class Backend, class Idx>
std::shared_ptr< distributed_matrix<Backend> > graph_perm_matrix(
        communicator comm, ptrdiff_t col_beg, ptrdiff_t col_end,
        const std::vector<Idx> &perm)
{
    typedef typename Backend::value_type value_type;
    typedef backend::crs<value_type> build_matrix;

    AMGCL_TIC("perm matrix");

    ptrdiff_t n = perm.size();
    ptrdiff_t ncols = col_end - col_beg;

    auto i_loc = std::make_shared<build_matrix>();
    auto i_rem = std::make_shared<build_matrix>();

    build_matrix &I_loc = *i_loc;
    build_matrix &I_rem = *i_rem;

    I_loc.set_size(n, ncols, false);
    I_rem.set_size(n, 0, false);

    I_loc.ptr[0] = 0;
    I_rem.ptr[0] = 0;

#pragma omp parallel for
    for(ptrdiff_t i = 0; i < n; ++i) {
        ptrdiff_t j = perm[i];

        if (col_beg <= j && j < col_end) {
            I_loc.ptr[i+1] = 1;
            I_rem.ptr[i+1] = 0;
        } else {
            I_loc.ptr[i+1] = 0;
            I_rem.ptr[i+1] = 1;
        }
    }

    I_loc.set_nonzeros(I_loc.scan_row_sizes());
    I_rem.set_nonzeros(I_rem.scan_row_sizes());

#pragma omp parallel for
    for(ptrdiff_t i = 0; i < n; ++i) {
        ptrdiff_t j = perm[i];

        if (col_beg <= j && j < col_end) {
            ptrdiff_t k = I_loc.ptr[i];
            I_loc.col[k] = j - col_beg;
            I_loc.val[k] = math::identity<value_type>();
        } else {
            ptrdiff_t k = I_rem.ptr[i];
            I_rem.col[k] = j;
            I_rem.val[k] = math::identity<value_type>();
        }
    }

    AMGCL_TOC("perm matrix");
    return std::make_shared< distributed_matrix<Backend> >(comm, i_loc, i_rem);
}

} // namespace partition
} // namespace mpi
} // namespace amgcl

#endif
