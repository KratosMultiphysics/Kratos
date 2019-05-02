#ifndef AMGCL_MPI_DIRECT_SOLVER_PASTIX_HPP
#define AMGCL_MPI_DIRECT_SOLVER_PASTIX_HPP

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
\file   amgcl/mpi/direct_solver/pastix.hpp
\author Denis Demidov <dennis.demidov@gmail.com>
\brief  Wrapper for PaStiX distributed sparse solver.

See http://pastix.gforge.inria.fr
*/

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <type_traits>

#include <amgcl/util.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/direct_solver/solver_base.hpp>

extern "C" {
#include <pastix.h>
}

namespace amgcl {
namespace mpi {
namespace direct {

/// Provides distributed direct solver interface for pastix solver.
/**
 * \sa http://pastix.gforge.inria.fr, \cite Henon2002
 */
template <typename value_type, bool Distrib=false>
class pastix : public solver_base< value_type, pastix<value_type, Distrib> > {
    public:
        static_assert(
                 std::is_same<value_type, float >::value ||
                 std::is_same<value_type, double>::value,
                 "Unsupported value type for pastix solver"
                );

        typedef backend::crs<value_type> build_matrix;

        struct params {
            int max_rows_per_process;

            params()
                : max_rows_per_process(50000)
            {}

#ifndef AMGCL_NO_BOOST
            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, max_rows_per_process)
            {
                check_params(p, {"max_rows_per_process"});
            }

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, max_rows_per_process);
            }
#endif
        };

        /// Constructor.
        template <class Matrix>
        pastix(communicator comm, const Matrix &A,
                const params &prm = params()) : prm(prm), pastix_data(0)
        {
            static_cast<Base*>(this)->init(comm, A);
        }

        static size_t coarse_enough() {
            return 10000;
        }

        int comm_size(int n) const {
            return Distrib ? (n + prm.max_rows_per_process - 1) / prm.max_rows_per_process : 1;
        }

        void init(communicator C, const build_matrix &A) {
            comm = C;
            nrows = A.nrows;
            ptr.assign(A.ptr, A.ptr + A.nrows + 1);
            col.assign(A.col, A.col + A.nnz);
            val.assign(A.val, A.val + A.nnz);

            row.resize(nrows);
            perm.resize(nrows);

            if (!Distrib) inv_perm.resize(nrows);

            std::vector<int> domain = comm.exclusive_sum(nrows);

            // PaStiX needs 1-based matrices:
            for(pastix_int_t &p : ptr) ++p;
            for(pastix_int_t &c : col) ++c;

            for(int i = 0, j = domain[comm.rank]; i < nrows; ++i)
                row[i] = ++j;

            // Initialize parameters with default values:
            iparm[IPARM_MODIFY_PARAMETER] = API_NO;
            call_pastix(API_TASK_INIT, API_TASK_INIT);

            // Factorize the matrix.
#ifdef NDEBUG
            iparm[IPARM_VERBOSE        ] = API_VERBOSE_NOT;
#else
            iparm[IPARM_VERBOSE        ] = API_VERBOSE_YES;
#endif
            iparm[IPARM_RHS_MAKING     ] = API_RHS_B;
            iparm[IPARM_SYM            ] = API_SYM_NO;
            iparm[IPARM_FACTORIZATION  ] = API_FACT_LU;
            iparm[IPARM_TRANSPOSE_SOLVE] = API_YES;
#ifdef _OPENMP
            iparm[IPARM_THREAD_NBR]      = omp_get_max_threads();
#endif
            call_pastix(API_TASK_ORDERING, API_TASK_NUMFACT);
        }

        /// Cleans up internal PaStiX data.
        ~pastix() {
            if(pastix_data) call_pastix(API_TASK_CLEAN, API_TASK_CLEAN);
        }

        /// Solves the problem for the given right-hand side.
        /**
         * \param rhs The right-hand side.
         * \param x   The solution.
         */
        template <class Vec1, class Vec2>
        void solve(const Vec1 &rhs, Vec2 &x) const {
            for(int i = 0; i < nrows; ++i) x[i] = rhs[i];
            call_pastix(API_TASK_SOLVE, API_TASK_SOLVE, &x[0]);
        }
    private:
        typedef solver_base< value_type, pastix<value_type, Distrib> > Base;
        params prm;
        amgcl::mpi::communicator comm;

        int nrows;

        // Pastix internal data.
        mutable pastix_data_t *pastix_data;

        // Pastix parameters
        mutable pastix_int_t   iparm[IPARM_SIZE];
        mutable double         dparm[DPARM_SIZE];

        std::vector<pastix_int_t> ptr;
        std::vector<pastix_int_t> col;
        std::vector<value_type>   val;

        // Local to global mapping
        std::vector<pastix_int_t> row;

        // Permutation array
        std::vector<pastix_int_t> perm;
        std::vector<pastix_int_t> inv_perm;

        void call_pastix(int beg, int end, value_type *x = NULL) const {
            iparm[IPARM_START_TASK] = beg;
            iparm[IPARM_END_TASK  ] = end;

            call_pastix(x);
        }

        void call_pastix(double *x) const {
            if (Distrib) {
                d_dpastix(&pastix_data, comm, nrows,
                        const_cast<pastix_int_t*>(&ptr[0]),
                        const_cast<pastix_int_t*>(&col[0]),
                        const_cast<double*      >(&val[0]),
                        const_cast<pastix_int_t*>(&row[0]),
                        const_cast<pastix_int_t*>(&perm[0]),
                        NULL, x, 1, iparm, dparm
                        );
            } else {
                d_pastix(&pastix_data, comm, nrows,
                        const_cast<pastix_int_t*>(&ptr[0]),
                        const_cast<pastix_int_t*>(&col[0]),
                        const_cast<double*      >(&val[0]),
                        const_cast<pastix_int_t*>(&perm[0]),
                        const_cast<pastix_int_t*>(&inv_perm[0]),
                        x, 1, iparm, dparm
                        );
            }
        }

        void call_pastix(float *x) const {
            if (Distrib) {
                s_dpastix(&pastix_data, comm, nrows,
                        const_cast<pastix_int_t*>(&ptr[0]),
                        const_cast<pastix_int_t*>(&col[0]),
                        const_cast<float*       >(&val[0]),
                        const_cast<pastix_int_t*>(&row[0]),
                        const_cast<pastix_int_t*>(&perm[0]),
                        NULL, x, 1, iparm, dparm
                        );
            } else {
                s_pastix(&pastix_data, comm, nrows,
                        const_cast<pastix_int_t*>(&ptr[0]),
                        const_cast<pastix_int_t*>(&col[0]),
                        const_cast<float*       >(&val[0]),
                        const_cast<pastix_int_t*>(&perm[0]),
                        const_cast<pastix_int_t*>(&inv_perm[0]),
                        x, 1, iparm, dparm
                        );
            }
        }
};

} // namespace direct
} // namespace mpi
} // namespace amgcl


#endif
