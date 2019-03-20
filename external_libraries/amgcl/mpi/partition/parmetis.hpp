#ifndef AMGCL_MPI_REPARTITION_PARMETIS_HPP
#define AMGCL_MPI_REPARTITION_PARMETIS_HPP

/*
The MIT License

Copyright (c) 2012-2019 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/mpi/partition/parmetis.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  ParMETIS partitioner.
 */

#include <memory>

#include <amgcl/backend/interface.hpp>
#include <amgcl/value_type/interface.hpp>
#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>
#include <amgcl/mpi/partition/util.hpp>

#include <parmetis.h>

namespace amgcl {
namespace mpi {
namespace partition {

template <class Backend>
struct parmetis {
    typedef typename Backend::value_type value_type;
    typedef distributed_matrix<Backend>  matrix;

    struct params {
        bool      enable;
        ptrdiff_t min_per_proc;
        int       shrink_ratio;

        params() :
            enable(false), min_per_proc(10000), shrink_ratio(8)
        {}

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, enable),
              AMGCL_PARAMS_IMPORT_VALUE(p, min_per_proc),
              AMGCL_PARAMS_IMPORT_VALUE(p, shrink_ratio)
        {
            check_params(p, {"enable", "min_per_proc", "shrink_ratio"});
        }

        void get(
                boost::property_tree::ptree &p,
                const std::string &path = ""
                ) const
        {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, enable);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, min_per_proc);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, shrink_ratio);
        }
#endif
    } prm;

    parmetis(const params &prm = params()) : prm(prm) {}

    bool is_needed(const matrix &A) const {
        if (!prm.enable) return false;

        communicator comm = A.comm();
        ptrdiff_t n = A.loc_rows();
        std::vector<ptrdiff_t> row_dom = comm.exclusive_sum(n);

        int non_empty = 0;
        ptrdiff_t min_n = std::numeric_limits<ptrdiff_t>::max();
        for(int i = 0; i < comm.size; ++i) {
            ptrdiff_t m = row_dom[i+1] - row_dom[i];
            if (m) {
                min_n = std::min(min_n, m);
                ++non_empty;
            }
        }

        return (non_empty > 1) && (min_n <= prm.min_per_proc);
    }

    std::shared_ptr<matrix> operator()(const matrix &A, unsigned block_size = 1) const {
        communicator comm = A.comm();
        idx_t n = A.loc_rows();
        ptrdiff_t row_beg = A.loc_col_shift();

        // Partition the graph.
        int active = (n > 0);
        int active_ranks = comm.reduce(MPI_SUM, active);

        idx_t npart = std::max(1, active_ranks / prm.shrink_ratio);

        if (comm.rank == 0)
            std::cout << "Partitioning[ParMETIS] " << active_ranks << " -> " << npart << std::endl;


        std::vector<ptrdiff_t> perm(n);
        ptrdiff_t col_beg, col_end;

        if (npart == 1) {
            col_beg = (comm.rank == 0) ? 0 : A.glob_rows();
            col_end = A.glob_rows();

            for(ptrdiff_t i = 0; i < n; ++i) {
                perm[i] = row_beg + i;
            }
        } else {
            if (block_size == 1) {
                std::tie(col_beg, col_end) = partition(A, npart, perm);
            } else {
                typedef typename math::scalar_of<value_type>::type scalar;
                typedef backend::builtin<scalar> sbackend;
                ptrdiff_t np = n / block_size;

                distributed_matrix<sbackend> A_pw(A.comm(),
                    pointwise_matrix(*A.local(),  block_size),
                    pointwise_matrix(*A.remote(), block_size)
                    );

                std::vector<ptrdiff_t> perm_pw(np);

                std::tie(col_beg, col_end) = partition(A_pw, npart, perm_pw);

                col_beg *= block_size;
                col_end *= block_size;

                for(ptrdiff_t ip = 0; ip < np; ++ip) {
                    ptrdiff_t i = ip * block_size;
                    ptrdiff_t j = perm_pw[ip] * block_size;

                    for(unsigned k = 0; k < block_size; ++k)
                        perm[i + k] = j + k;
                }
            }
        }

        return graph_perm_matrix<Backend>(comm, col_beg, col_end, perm);
    }

    template <class B>
    std::tuple<ptrdiff_t, ptrdiff_t>
    partition(const distributed_matrix<B> &A, idx_t npart, std::vector<ptrdiff_t> &perm) const {
        communicator comm = A.comm();
        idx_t n = A.loc_rows();
        int active = (n > 0);

        std::vector<idx_t> ptr;
        std::vector<idx_t> col;

        symm_graph(A, ptr, col);

        idx_t wgtflag = 0;
        idx_t numflag = 0;
        idx_t options = 0;
        idx_t edgecut = 0;
        idx_t ncon    = 1;

        std::vector<real_t> tpwgts(npart, 1.0 / npart);
        std::vector<real_t> ubvec(ncon, 1.05);
        std::vector<idx_t>  part(n);
        if (!n) part.reserve(1); // So that part.data() is not NULL

        MPI_Comm scomm;
        MPI_Comm_split(comm, active ? 0 : MPI_UNDEFINED, comm.rank, &scomm);

        if (active) {
            communicator sc(scomm);
            std::vector<idx_t> vtxdist = sc.exclusive_sum(n);

            sc.check(
                    METIS_OK == ParMETIS_V3_PartKway( &vtxdist[0], &ptr[0],
                        &col[0], NULL, NULL, &wgtflag, &numflag, &ncon,
                        &npart, &tpwgts[0], &ubvec[0], &options, &edgecut,
                        &part[0], &scomm),
                    "Error in ParMETIS"
                    );

            MPI_Comm_free(&scomm);
        }

        return graph_perm_index(comm, npart, part, perm);
    }
};


} // namespace partition
} // namespace mpi
} // namespace amgcl

#endif
