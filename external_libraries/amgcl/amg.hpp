#ifndef AMGCL_AMG_HPP
#define AMGCL_AMG_HPP

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
 * \file   amgcl/amg.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  An AMG preconditioner.
 */

#include <iostream>
#include <iomanip>
#include <list>
#include <memory>

#ifdef AMGCL_ASYNC_SETUP
#  include <atomic>
#  include <thread>
#  include <mutex>
#  include <condition_variable>
#endif

#include <amgcl/backend/builtin.hpp>
#include <amgcl/solver/detail/default_inner_product.hpp>
#include <amgcl/util.hpp>

/// Primary namespace.
namespace amgcl {

/// Algebraic multigrid method.
/**
 * AMG is one the most effective methods for solution of large sparse
 * unstructured systems of equations, arising, for example, from discretization
 * of PDEs on unstructured grids \cite Trottenberg2001. The method can be used
 * as a black-box solver for various computational problems, since it does not
 * require any information about the underlying geometry.
 *
 * The three template parameters allow the user to select the exact components
 * of the method:
 *  1. *Backend* to transfer the constructed hierarchy to,
 *  2. *Coarsening* strategy for hierarchy construction, and
 *  3. *Relaxation* scheme (smoother to use during the solution phase).
 *
 * Instance of the class builds the AMG hierarchy for the given system matrix
 * and is intended to be used as a preconditioner.
 */
template <
    class Backend,
    template <class> class Coarsening,
    template <class> class Relax
    >
class amg {
    public:
        typedef Backend backend_type;

        typedef typename Backend::value_type value_type;
        typedef typename Backend::matrix     matrix;
        typedef typename Backend::vector     vector;

        typedef Coarsening<Backend>          coarsening_type;
        typedef Relax<Backend>               relax_type;

        typedef typename backend::builtin<value_type>::matrix build_matrix;

        typedef typename math::scalar_of<value_type>::type scalar_type;

        /// Backend parameters.
        typedef typename Backend::params     backend_params;

        /// Parameters of the method.
        /**
         * The amgcl::amg::params struct includes parameters for each
         * component of the method as well as some universal parameters.
         */
        struct params {
            typedef typename coarsening_type::params coarsening_params;
            typedef typename relax_type::params relax_params;

            coarsening_params coarsening;   ///< Coarsening parameters.
            relax_params      relax;        ///< Relaxation parameters.

            /// Specifies when level is coarse enough to be solved directly.
            /**
             * If number of variables at a next level in the hierarchy becomes
             * lower than this threshold, then the hierarchy construction is
             * stopped and the linear system is solved directly at this level.
             */
            unsigned coarse_enough;

            /// Use direct solver at the coarsest level.
            /**
             * When set, the coarsest level is solved with a direct solver.
             * Otherwise a smoother is used as a solver.
             */
            bool direct_coarse;

            /// Maximum number of levels.
            /** If this number is reached while the size of the last level is
             * greater that `coarse_enough`, then the coarsest level will not
             * be solved exactly, but will use a smoother.
             */
            unsigned max_levels;

            /// Number of pre-relaxations.
            unsigned npre;

            /// Number of post-relaxations.
            unsigned npost;

            /// Number of cycles (1 for V-cycle, 2 for W-cycle, etc.).
            unsigned ncycle;

            /// Number of cycles to make as part of preconditioning.
            unsigned pre_cycles;

#ifdef AMGCL_ASYNC_SETUP
            /// Asynchronous setup.
            /** Starts cycling as soon as the first level is (partially)
             * constructed. May be useful for GPGPU backends as a way to split
             * the work between the host CPU and the compute device(s).
             */
            bool async_setup;
#endif

            params() :
                coarse_enough( Backend::direct_solver::coarse_enough() ),
                direct_coarse(true),
                max_levels( std::numeric_limits<unsigned>::max() ),
                npre(1), npost(1), ncycle(1), pre_cycles(1)
#ifdef AMGCL_ASYNC_SETUP
                , async_setup(false)
#endif
            {}

#ifndef AMGCL_NO_BOOST
            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_CHILD(p, coarsening),
                  AMGCL_PARAMS_IMPORT_CHILD(p, relax),
                  AMGCL_PARAMS_IMPORT_VALUE(p, coarse_enough),
                  AMGCL_PARAMS_IMPORT_VALUE(p, direct_coarse),
                  AMGCL_PARAMS_IMPORT_VALUE(p, max_levels),
                  AMGCL_PARAMS_IMPORT_VALUE(p, npre),
                  AMGCL_PARAMS_IMPORT_VALUE(p, npost),
                  AMGCL_PARAMS_IMPORT_VALUE(p, ncycle),
                  AMGCL_PARAMS_IMPORT_VALUE(p, pre_cycles)
#ifdef AMGCL_ASYNC_SETUP
                , AMGCL_PARAMS_IMPORT_VALUE(p, async_setup)
#endif
            {
                check_params(p, {"coarsening", "relax", "coarse_enough",  "direct_coarse", "max_levels", "npre", "npost", "ncycle",  "pre_cycles"
#ifdef AMGCL_ASYNC_SETUP
                        , "async_setup"
#endif
                        });

                precondition(max_levels > 0, "max_levels should be positive");
            }

            void get(
                    boost::property_tree::ptree &p,
                    const std::string &path = ""
                    ) const
            {
                AMGCL_PARAMS_EXPORT_CHILD(p, path, coarsening);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, relax);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, coarse_enough);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, direct_coarse);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, max_levels);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, npre);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, npost);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, ncycle);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, pre_cycles);
#ifdef AMGCL_ASYNC_SETUP
                AMGCL_PARAMS_EXPORT_VALUE(p, path, async_setup);
#endif
            }
#endif
        } prm;

        /// Builds the AMG hierarchy for the system matrix.
        /**
         * The input matrix is copied here and is safe to delete afterwards.
         *
         * \param A The system matrix. Should be convertible to
         *          amgcl::backend::crs<>.
         * \param p AMG parameters.
         *
         * \sa amgcl/adapter/crs_tuple.hpp
         */
        template <class Matrix>
        amg(
                const Matrix &M,
                const params &p = params(),
                const backend_params &bprm = backend_params()
           ) : prm(p)
        {
            auto A = std::make_shared<build_matrix>(M);
            sort_rows(*A);

            do_init(A, bprm);
        }

        /// Builds the AMG hierarchy for the system matrix.
        /**
         * The shared pointer to the input matrix is passed here. The matrix
         * will not be copied and should out-live the amg instance.
         * The matrix should be either in amgcl::backend::crs<T> format, or
         * inherit from the class and override its ptr(), col(), and val()
         * virtual functions.
         *
         * \param A The system matrix.
         * \param p AMG parameters.
         *
         * \sa amgcl/adapter/crs_tuple.hpp
         */
        amg(
                std::shared_ptr<build_matrix> A,
                const params &p = params(),
                const backend_params &bprm = backend_params()
           ) : prm(p)
        {
            do_init(A, bprm);
        }

#ifdef AMGCL_ASYNC_SETUP
        ~amg() {
            done = true;
            if (prm.async_setup) init_thread.join();
        }
#endif

        /// Performs single V-cycle for the given right-hand side and solution.
        /**
         * \param rhs Right-hand side vector.
         * \param x   Solution vector.
         */
        template <class Vec1, class Vec2>
        void cycle(const Vec1 &rhs, Vec2 &&x) const {
            cycle(levels.begin(), rhs, x);
        }

        /// Performs single V-cycle after clearing x.
        /**
         * This is intended for use as a preconditioning procedure.
         *
         * \param rhs Right-hand side vector.
         * \param x   Solution vector.
         */
        template <class Vec1, class Vec2>
        void apply(const Vec1 &rhs, Vec2 &&x) const {
            if (prm.pre_cycles) {
                backend::clear(x);
                for(unsigned i = 0; i < prm.pre_cycles; ++i)
                    cycle(rhs, x);
            } else {
                backend::copy(rhs, x);
            }
        }

        /// Returns the system matrix from the finest level.
        std::shared_ptr<matrix> system_matrix_ptr() const {
            return levels.front().A;
        }

        const matrix& system_matrix() const {
            return *system_matrix_ptr();
        }

        size_t bytes() const {
            size_t b = 0;
            for(const auto &lvl : levels) b += lvl.bytes();
            return b;
        }
    private:
        struct level {
            size_t m_rows, m_nonzeros;

            std::shared_ptr<vector> f;
            std::shared_ptr<vector> u;
            std::shared_ptr<vector> t;

            std::shared_ptr<matrix> A;
            std::shared_ptr<matrix> P;
            std::shared_ptr<matrix> R;

            std::shared_ptr< typename Backend::direct_solver > solve;

            std::shared_ptr<relax_type> relax;

            size_t bytes() const {
                size_t b = 0;

                if (f) b += backend::bytes(*f);
                if (u) b += backend::bytes(*u);
                if (t) b += backend::bytes(*t);

                if (A) b += backend::bytes(*A);
                if (P) b += backend::bytes(*P);
                if (R) b += backend::bytes(*R);

                if (solve) b += backend::bytes(*solve);
                if (relax) b += backend::bytes(*relax);

                return b;
            }

            level() {}

            level(std::shared_ptr<build_matrix> A,
                    params &prm, const backend_params &bprm)
                : m_rows(backend::rows(*A)), m_nonzeros(backend::nonzeros(*A))
            {
                AMGCL_TIC("move to backend");
                f = Backend::create_vector(m_rows, bprm);
                u = Backend::create_vector(m_rows, bprm);
                t = Backend::create_vector(m_rows, bprm);
                this->A = Backend::copy_matrix(A, bprm);
                AMGCL_TOC("move to backend");

                AMGCL_TIC("relaxation");
                relax = std::make_shared<relax_type>(*A, prm.relax, bprm);
                AMGCL_TOC("relaxation");
            }

            std::shared_ptr<build_matrix> step_down(
                    std::shared_ptr<build_matrix> A,
                    coarsening_type &C, const backend_params &bprm)
            {
                AMGCL_TIC("transfer operators");
                std::shared_ptr<build_matrix> P, R;

                try {
                    std::tie(P, R) = C.transfer_operators(*A);
                } catch(error::empty_level) {
                    AMGCL_TOC("transfer operators");
                    return std::shared_ptr<build_matrix>();
                }

                sort_rows(*P);
                sort_rows(*R);
                AMGCL_TOC("transfer operators");

                AMGCL_TIC("move to backend");
                this->P = Backend::copy_matrix(P, bprm);
                this->R = Backend::copy_matrix(R, bprm);
                AMGCL_TOC("move to backend");

                AMGCL_TIC("coarse operator");
                A = C.coarse_operator(*A, *P, *R);
                sort_rows(*A);
                AMGCL_TOC("coarse operator");

                return A;
            }

            void create_coarse(
                    std::shared_ptr<build_matrix> A,
                    const backend_params &bprm, bool single_level)
            {
                m_rows     = backend::rows(*A);
                m_nonzeros = backend::nonzeros(*A);

                u = Backend::create_vector(m_rows, bprm);
                f = Backend::create_vector(m_rows, bprm);

                solve = Backend::create_solver(A, bprm);
                if (single_level)
                    this->A = Backend::copy_matrix(A, bprm);
            }

            size_t rows() const {
                return m_rows;
            }

            size_t nonzeros() const {
                return m_nonzeros;
            }
        };

        typedef typename std::list<level>::const_iterator level_iterator;

        std::list<level> levels;
#ifdef AMGCL_ASYNC_SETUP
        std::thread init_thread;
        std::atomic<bool> done;
        mutable std::mutex levels_mx;
        mutable std::condition_variable ready_to_cycle;
#endif

        void init(
                std::shared_ptr<build_matrix> A,
                const backend_params &bprm = backend_params()
           )
        {
            precondition(
                    backend::rows(*A) == backend::cols(*A),
                    "Matrix should be square!"
                    );

            bool direct_coarse_solve = true;

            coarsening_type C(prm.coarsening);

            while( backend::rows(*A) > prm.coarse_enough) {
                {
#ifdef AMGCL_ASYNC_SETUP
                    if (done) break;
                    std::lock_guard<std::mutex> lock(levels_mx);
#endif
                    levels.push_back( level(A, prm, bprm) );
                }
#ifdef AMGCL_ASYNC_SETUP
                if (done) break;
                ready_to_cycle.notify_all();
#endif
                if (levels.size() >= prm.max_levels) break;

                A = levels.back().step_down(A, C, bprm);
                if (!A) {
                    // Zero-sized coarse level. Probably the system matrix on
                    // this level is diagonal, should be easily solvable with a
                    // couple of smoother iterations.
                    direct_coarse_solve = false;
                    break;
                }
            }

            if (!A || backend::rows(*A) > prm.coarse_enough) {
                // The coarse matrix is still too big to be solved directly.
                direct_coarse_solve = false;
            }

            if (direct_coarse_solve) {
                AMGCL_TIC("coarsest level");
                if (prm.direct_coarse) {
                    level l;
                    l.create_coarse(A, bprm, levels.empty());

                    {
#ifdef AMGCL_ASYNC_SETUP
                        std::lock_guard<std::mutex> lock(levels_mx);
#endif
                        levels.push_back( l );
                    }
                } else {
                    {
#ifdef AMGCL_ASYNC_SETUP
                        std::lock_guard<std::mutex> lock(levels_mx);
#endif
                        levels.push_back( level(A, prm, bprm) );
                    }
                }
#ifdef AMGCL_ASYNC_SETUP
                ready_to_cycle.notify_all();
#endif
                AMGCL_TOC("coarsest level");
            }
        }

        void do_init(
                std::shared_ptr<build_matrix> A,
                const backend_params &bprm = backend_params()
                )
        {
#ifdef AMGCL_ASYNC_SETUP
            done = false;
            if (prm.async_setup) {
                init_thread = std::thread(&amg::init, this, A, bprm);
                {
                    std::unique_lock<std::mutex> lock(levels_mx);
                    while(levels.empty()) ready_to_cycle.wait(lock);
                }
            } else
#endif
            {
                init(A, bprm);
            }
        }
        template <class Vec1, class Vec2>
        void cycle(level_iterator lvl, const Vec1 &rhs, Vec2 &x) const
        {
            level_iterator nxt = lvl, end;

            {
#ifdef AMGCL_ASYNC_SETUP
                std::lock_guard<std::mutex> lock(levels_mx);
#endif
                ++nxt;
                end = levels.end();
            }

            if (nxt == end) {
                if (lvl->solve) {
                    AMGCL_TIC("coarse");
                    (*lvl->solve)(rhs, x);
                    AMGCL_TOC("coarse");
                } else {
                    AMGCL_TIC("relax");
                    for(size_t i = 0; i < prm.npre;  ++i) lvl->relax->apply_pre(*lvl->A, rhs, x, *lvl->t);
                    for(size_t i = 0; i < prm.npost; ++i) lvl->relax->apply_post(*lvl->A, rhs, x, *lvl->t);
                    AMGCL_TOC("relax");
                }
            } else {
                for (size_t j = 0; j < prm.ncycle; ++j) {
                    AMGCL_TIC("relax");
                    for(size_t i = 0; i < prm.npre; ++i)
                        lvl->relax->apply_pre(*lvl->A, rhs, x, *lvl->t);
                    AMGCL_TOC("relax");

                    backend::residual(rhs, *lvl->A, x, *lvl->t);

                    backend::spmv(math::identity<scalar_type>(), *lvl->R, *lvl->t, math::zero<scalar_type>(), *nxt->f);

                    backend::clear(*nxt->u);
                    cycle(nxt, *nxt->f, *nxt->u);

                    backend::spmv(math::identity<scalar_type>(), *lvl->P, *nxt->u, math::identity<scalar_type>(), x);

                    AMGCL_TIC("relax");
                    for(size_t i = 0; i < prm.npost; ++i)
                        lvl->relax->apply_post(*lvl->A, rhs, x, *lvl->t);
                    AMGCL_TOC("relax");
                }
            }
        }

    template <class B, template <class> class C, template <class> class R>
    friend std::ostream& operator<<(std::ostream &os, const amg<B, C, R> &a);
};

/// Sends information about the AMG hierarchy to output stream.
template <class B, template <class> class C, template <class> class R>
std::ostream& operator<<(std::ostream &os, const amg<B, C, R> &a)
{
    typedef typename amg<B, C, R>::level level;
    std::ios_base::fmtflags ff(os.flags());
    auto fp = os.precision();

    size_t sum_dof = 0;
    size_t sum_nnz = 0;
    size_t sum_mem = 0;

    for(const level &lvl : a.levels) {
        sum_dof += lvl.rows();
        sum_nnz += lvl.nonzeros();
        sum_mem += lvl.bytes();
    }

    os << "Number of levels:    "   << a.levels.size()
        << "\nOperator complexity: " << std::fixed << std::setprecision(2)
        << 1.0 * sum_nnz / a.levels.front().nonzeros()
        << "\nGrid complexity:     " << std::fixed << std::setprecision(2)
        << 1.0 * sum_dof / a.levels.front().rows()
        << "\nMemory footprint:    " << human_readable_memory(sum_mem)
        << "\n\n"
           "level     unknowns       nonzeros      memory\n"
           "---------------------------------------------\n";

    size_t depth = 0;
    for(const level &lvl : a.levels) {
        os << std::setw(5)  << depth++
            << std::setw(13) << lvl.rows()
            << std::setw(15) << lvl.nonzeros()
            << std::setw(12) << human_readable_memory(lvl.bytes())
            << " (" << std::setw(5) << std::fixed << std::setprecision(2)
            << 100.0 * lvl.nonzeros() / sum_nnz
            << "%)" << std::endl;
    }

    os.flags(ff);
    os.precision(fp);
    return os;
}

} // namespace amgcl

#endif
