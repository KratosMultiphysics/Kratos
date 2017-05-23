#ifndef AMGCL_RELAXATION_DETAIL_ILU_SOLVE_HPP
#define AMGCL_RELAXATION_DETAIL_ILU_SOLVE_HPP

#include <amgcl/backend/interface.hpp>

namespace amgcl {
namespace relaxation {
namespace detail {

template <class Backend>
class ilu_solve {
    private:
        typedef typename Backend::value_type value_type;
        typedef typename Backend::matrix matrix;
        typedef typename Backend::vector vector;
        typedef typename Backend::matrix_diagonal matrix_diagonal;
        typedef typename backend::builtin<value_type>::matrix build_matrix;
        typedef typename math::scalar_of<value_type>::type scalar_type;

        unsigned jacobi_iters;

        boost::shared_ptr<matrix> L;
        boost::shared_ptr<matrix> U;
        boost::shared_ptr<matrix_diagonal> D;
        boost::shared_ptr<vector> t1, t2;

    public:
        template <class Params>
        ilu_solve(
                boost::shared_ptr<build_matrix> L,
                boost::shared_ptr<build_matrix> U,
                boost::shared_ptr<backend::numa_vector<value_type> > D,
                const Params &prm,
                const typename Backend::params &bprm
                ) :
            jacobi_iters(prm.jacobi_iters),
            L(Backend::copy_matrix(L, bprm)),
            U(Backend::copy_matrix(U, bprm)),
            D(Backend::copy_vector(D, bprm)),
            t1(Backend::create_vector(backend::rows(*L), bprm)),
            t2(Backend::create_vector(backend::rows(*L), bprm))
        {}

        template <class Vector>
        void solve(Vector &x) {
            vector *y0 = t1.get();
            vector *y1 = t2.get();

            backend::copy(x, *y0);
            for(unsigned i = 0; i < jacobi_iters; ++i) {
                backend::residual(x, *L, *y0, *y1);
                std::swap(y0, y1);
            }

            backend::copy(*y0, x);
            for(unsigned i = 0; i < jacobi_iters; ++i) {
                backend::residual(*y0, *U, x, *y1);
                backend::vmul(math::identity<scalar_type>(), *D, *y1, math::zero<scalar_type>(), x);
            }
        }
};

template <class value_type>
class ilu_solve< backend::builtin<value_type> > {
    public:
        typedef backend::builtin<value_type> Backend;
        typedef typename Backend::matrix matrix;
        typedef typename Backend::vector vector;
        typedef typename Backend::matrix_diagonal matrix_diagonal;
        typedef typename backend::builtin<value_type>::matrix build_matrix;
        typedef typename Backend::rhs_type rhs_type;
        typedef typename math::scalar_of<value_type>::type scalar_type;

        template <class Params>
        ilu_solve(
                boost::shared_ptr<build_matrix> L,
                boost::shared_ptr<build_matrix> U,
                boost::shared_ptr<backend::numa_vector<value_type> > D,
                const Params &prm, const typename Backend::params &bprm
                ) : is_serial(prm.serial || num_threads() < 4)
        {
            if (is_serial)
                serial_init(L, U, D, bprm);
            else
                parallel_init(L, U, D, bprm);
        }

        template <class Vector>
        void solve(Vector &x) {
            if (is_serial)
                serial_solve(x);
            else
                parallel_solve(x);
        }

    private:
        bool is_serial;

        static int num_threads() {
#ifdef _OPENMP
            return omp_get_max_threads();
#else
            return 1;
#endif
        }

        static int thread_id() {
#ifdef _OPENMP
            return omp_get_thread_num();
#else
            return 0;
#endif
        }

        // copies of the input matrices for the fallback (serial)
        // implementation:
        boost::shared_ptr<matrix>          L;
        boost::shared_ptr<matrix>          U;
        boost::shared_ptr<matrix_diagonal> D;

        void serial_init(
                boost::shared_ptr<build_matrix>    L,
                boost::shared_ptr<build_matrix>    U,
                boost::shared_ptr<matrix_diagonal> D,
                const typename Backend::params&
                )
        {
            this->L = L;
            this->U = U;
            this->D = D;
        }

        template <class Vector>
        void serial_solve(Vector &x) {
            const size_t n = backend::rows(*L);

            const matrix          &L = *(this->L);
            const matrix          &U = *(this->U);
            const matrix_diagonal &D = *(this->D);

            for(size_t i = 0; i < n; i++) {
                for(ptrdiff_t j = L.ptr[i], e = L.ptr[i+1]; j < e; ++j)
                    x[i] -= L.val[j] * x[L.col[j]];
            }

            for(size_t i = n; i-- > 0;) {
                for(ptrdiff_t j = U.ptr[i], e = U.ptr[i+1]; j < e; ++j)
                    x[i] -= U.val[j] * x[U.col[j]];
                x[i] = D[i] * x[i];
            }
        }

        // OpenMP solver for sparse triangular systems.
        // The solver uses level scheduling approach.
        // Each level (a set of matrix rows that can be computed independently)
        // is split into tasks, a task per thread, and the matrix data is
        // distributed across threads to improve cache and NUMA locality.
        template <bool lower>
        struct sptr_solve {
            // a task is a set of rows that can be computed independently by a
            // single thread.
            struct task {
                ptrdiff_t beg, end; // rows to process

                task(ptrdiff_t beg, ptrdiff_t end) : beg(beg), end(end) {}
            };

            int nthreads;

            // thread-specific storage:
            std::vector< std::vector<task>       > tasks;
            std::vector< std::vector<ptrdiff_t>  > ptr;
            std::vector< std::vector<ptrdiff_t>  > col;
            std::vector< std::vector<value_type> > val;
            std::vector< std::vector<ptrdiff_t>  > ord; // rows ordered by levels
            std::vector< std::vector<value_type> > D;

            template <class Matrix>
            sptr_solve(const Matrix &A, const value_type *_D = 0)
                : nthreads(num_threads()), tasks(nthreads),
                  ptr(nthreads), col(nthreads), val(nthreads), ord(nthreads)
            {
                ptrdiff_t n    = A.nrows;
                ptrdiff_t nlev = 0;

                std::vector<ptrdiff_t> level(n, 0);
                std::vector<ptrdiff_t> order(n, 0);


                // 1. split rows into levels.
                ptrdiff_t beg = lower ? 0 : n-1;
                ptrdiff_t end = lower ? n :  -1;
                ptrdiff_t inc = lower ? 1 :  -1;

                for(ptrdiff_t i = beg; i != end; i += inc) {
                    ptrdiff_t l = level[i];

                    for(ptrdiff_t j = A.ptr[i]; j < A.ptr[i+1]; ++j)
                        l = std::max(l, level[A.col[j]]+1);

                    level[i] = l;
                    nlev = std::max(nlev, l+1);
                }


                // 2. reorder matrix rows.
                std::vector<ptrdiff_t> start(nlev+1, 0);

                for(ptrdiff_t i = 0; i < n; ++i)
                    ++start[level[i]+1];

                std::partial_sum(start.begin(), start.end(), start.begin());

                for(ptrdiff_t i = 0; i < n; ++i)
                    order[start[level[i]]++] = i;

                std::rotate(start.begin(), start.end() - 1, start.end());
                start[0] = 0;


                // 3. Organize matrix rows into tasks.
                //    Each level is split into nthreads tasks.
                std::vector<ptrdiff_t> thread_rows(nthreads, 0);
                std::vector<ptrdiff_t> thread_cols(nthreads, 0);

#pragma omp parallel
                {
                    int tid = thread_id();
                    tasks[tid].reserve(nlev);

                    for(ptrdiff_t lev = 0; lev < nlev; ++lev) {
                        // split each level into tasks.
                        ptrdiff_t lev_size = start[lev+1] - start[lev];
                        ptrdiff_t chunk_size = (lev_size + nthreads - 1) / nthreads;

                        ptrdiff_t beg = std::min(tid * chunk_size, lev_size);
                        ptrdiff_t end = std::min(beg + chunk_size, lev_size);

                        beg += start[lev];
                        end += start[lev];

                        tasks[tid].push_back(task(beg, end));

                        // count rows and nonzeros in the current task
                        thread_rows[tid] += end - beg;
                        for(ptrdiff_t i = beg; i < end; ++i) {
                            ptrdiff_t j = order[i];
                            thread_cols[tid] += A.ptr[j+1] - A.ptr[j];
                        }
                    }
                }

                // 4. reorganize matrix data for better cache and NUMA locality.
                if (!lower) D.resize(nthreads);

#pragma omp parallel
                {
                    int tid = thread_id();

                    col[tid].reserve(thread_cols[tid]);
                    val[tid].reserve(thread_cols[tid]);
                    ord[tid].reserve(thread_rows[tid]);
                    ptr[tid].reserve(thread_rows[tid] + 1);
                    ptr[tid].push_back(0);

                    if (!lower) D[tid].reserve(thread_rows[tid]);

                    BOOST_FOREACH(task &t, tasks[tid]) {
                        ptrdiff_t loc_beg = ptr[tid].size() - 1;
                        ptrdiff_t loc_end = loc_beg;

                        for(ptrdiff_t r = t.beg; r < t.end; ++r, ++loc_end) {
                            ptrdiff_t i = order[r];
                            if (!lower) D[tid].push_back(_D[i]);

                            ord[tid].push_back(i);

                            for(ptrdiff_t j = A.ptr[i]; j < A.ptr[i+1]; ++j) {
                                col[tid].push_back(A.col[j]);
                                val[tid].push_back(A.val[j]);
                            }

                            ptr[tid].push_back(col[tid].size());
                        }

                        t.beg = loc_beg;
                        t.end = loc_end;
                    }
                }
            }

            template <class Vector>
            void solve(Vector &x) const {
#pragma omp parallel
                {
                    int tid = thread_id();

                    BOOST_FOREACH(const task &t, tasks[tid]) {
                        for(ptrdiff_t r = t.beg; r < t.end; ++r) {
                            ptrdiff_t i   = ord[tid][r];
                            ptrdiff_t beg = ptr[tid][r];
                            ptrdiff_t end = ptr[tid][r+1];

                            rhs_type X = math::zero<rhs_type>();
                            for(ptrdiff_t j = beg; j < end; ++j)
                                X += val[tid][j] * x[col[tid][j]];

                            if (lower)
                                x[i] -= X;
                            else
                                x[i] = D[tid][r] * (x[i] - X);
                        }

                        // each task corresponds to a level, so we need
                        // to synchronize across threads at this point:
#pragma omp barrier
                    }
                }
            }
        };

        boost::shared_ptr< sptr_solve<true > > lower;
        boost::shared_ptr< sptr_solve<false> > upper;

        void parallel_init(
                boost::shared_ptr<build_matrix> L,
                boost::shared_ptr<build_matrix> U,
                boost::shared_ptr<backend::numa_vector<value_type> > D,
                const typename Backend::params&
                )
        {
            lower = boost::make_shared< sptr_solve<true > >(*L, D->data());
            upper = boost::make_shared< sptr_solve<false> >(*U, D->data());
        }

        template <class Vector>
        void parallel_solve(Vector &x) {
            lower->solve(x);
            upper->solve(x);
        }
};

} // namespace detail
} // namespace relaxation
} // namespace amgcl

#endif
