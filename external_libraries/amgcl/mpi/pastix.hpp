#ifndef AMGCL_MPI_PASTIX_HPP
#define AMGCL_MPI_PASTIX_HPP

/*
The MIT License

Copyright (c) 2012-2016 Denis Demidov <dennis.demidov@gmail.com>

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
\file   amgcl/mpi/pastix.hpp
\author Denis Demidov <dennis.demidov@gmail.com>
\brief  Wrapper for PaStiX distributed sparse solver.

See http://pastix.gforge.inria.fr
*/

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/irange.hpp>
#include <boost/foreach.hpp>

#include <amgcl/util.hpp>
#include <amgcl/mpi/util.hpp>

extern "C" {
#include <pastix.h>
}

namespace amgcl {
namespace mpi {

/// Provides distributed direct solver interface for PaStiX solver.
/**
 * \sa http://pastix.gforge.inria.fr, \cite Henon2002
 */
template <typename value_type>
class PaStiX {
    public:
        BOOST_STATIC_ASSERT_MSG( (
                 boost::is_same<value_type, float >::value ||
                 boost::is_same<value_type, double>::value
                ), "Unsupported value type for PaStiX solver"
                );

        struct params {
            int max_rows_per_process;

            params()
                : max_rows_per_process(50000)
            {}

            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, max_rows_per_process)
            {
                AMGCL_PARAMS_CHECK(p, (max_rows_per_process))
            }

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, max_rows_per_process);
            }
        };

        /// The number of processes optimal for the given problem size.
        static int comm_size(int n_global_rows, const params &prm = params()) {
            return (n_global_rows + prm.max_rows_per_process - 1) / prm.max_rows_per_process;
        }

        /// Constructor.
        /**
         * \param comm MPI communicator containing processes to participate in
         *        solution of the problem. The number of processes in
         *        communicator should be (but not necessarily) equal to the
         *        result of comm_size().
         * \param n_local_rows Number of matrix rows belonging to the calling
         *        process.
         * \param ptr Start of each row in col and val arrays.
         * \param col Column numbers of nonzero elements.
         * \param val Values of nonzero elements.
         * \param prm Solver parameters.
         */
        template <class PRng, class CRng, class VRng>
        PaStiX(
                MPI_Comm mpi_comm,
                int n_local_rows,
                const PRng &p_ptr,
                const CRng &p_col,
                const VRng &p_val,
                const params& = params()
                )
            : comm(mpi_comm), nrows(n_local_rows), pastix_data(0),
              ptr(boost::begin(p_ptr), boost::end(p_ptr)),
              col(boost::begin(p_col), boost::end(p_col)),
              val(boost::begin(p_val), boost::end(p_val)),
              row(nrows), perm(nrows)
        {
            std::vector<int> domain = mpi::exclusive_sum(comm, nrows);

            boost::copy(
                    boost::irange(domain[comm.rank], domain[comm.rank + 1]),
                    row.begin()
                    );

            // PaStiX needs 1-based matrices:
            BOOST_FOREACH(pastix_int_t &p, ptr) ++p;
            BOOST_FOREACH(pastix_int_t &c, col) ++c;
            BOOST_FOREACH(pastix_int_t &r, row) ++r;

            // Initialize parameters with default values:
            iparm[IPARM_MODIFY_PARAMETER] = API_NO;
            call_pastix(API_TASK_INIT, API_TASK_INIT);

            // Factorize the matrix.
            iparm[IPARM_VERBOSE        ] = API_VERBOSE_NOT;
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
        ~PaStiX() {
            call_pastix(API_TASK_CLEAN, API_TASK_CLEAN);
        }

        /// Solves the problem for the given right-hand side.
        /**
         * \param rhs The right-hand side.
         * \param x   The solution.
         */
        template <class Vec1, class Vec2>
        void operator()(const Vec1 &rhs, Vec2 &x) const {
            boost::copy(rhs, &x[0]);
            call_pastix(API_TASK_SOLVE, API_TASK_SOLVE, &x[0]);
        }
    private:
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

        void call_pastix(int beg, int end, value_type *x = NULL) const {
            iparm[IPARM_START_TASK] = beg;
            iparm[IPARM_END_TASK  ] = end;

            call_pastix(x);
        }

        void call_pastix(double *x) const {
            d_dpastix(&pastix_data, comm, nrows,
                    const_cast<pastix_int_t*>(&ptr[0]),
                    const_cast<pastix_int_t*>(&col[0]),
                    const_cast<double*      >(&val[0]),
                    const_cast<pastix_int_t*>(&row[0]),
                    const_cast<pastix_int_t*>(&perm[0]),
                    NULL, x, 1, iparm, dparm
                   );
        }

        void call_pastix(float *x) const {
            s_dpastix(&pastix_data, comm, nrows,
                    const_cast<pastix_int_t*>(&ptr[0]),
                    const_cast<pastix_int_t*>(&col[0]),
                    const_cast<float*       >(&val[0]),
                    const_cast<pastix_int_t*>(&row[0]),
                    const_cast<pastix_int_t*>(&perm[0]),
                    NULL, x, 1, iparm, dparm
                   );
        }
};

}
}

#endif
