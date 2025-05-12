/* AMGCL */
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/ublas.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/preconditioner/runtime.hpp>

#ifdef AMGCL_GPGPU
#  include <amgcl/backend/vexcl.hpp>
#  include <amgcl/backend/vexcl_static_matrix.hpp>
#endif

#include "spaces/ublas_space.h"
#include "linear_solvers/amgcl_solver.h"

namespace Kratos {

#ifdef AMGCL_GPGPU
vex::Context& vexcl_context() {
    static vex::Context ctx(vex::Filter::Env);
    static bool run_once = [](){
        std::cout << "VexCL context:\n" << ctx << std::endl;
        return true;
    }();
    (void)run_once; // suppress "unused variable" warnings
    return ctx;
}

template <class TValue, int TBlockSize>
void register_vexcl_static_matrix_type() {
    static vex::scoped_program_header header(vexcl_context(),
            amgcl::backend::vexcl_static_matrix_declaration<TValue,TBlockSize>());
}
#endif

template <class TValue>
void AMGCLScalarSolve(
    typename TUblasSparseSpace<TValue>::MatrixType& rA,
    typename TUblasSparseSpace<TValue>::VectorType& rX,
    typename TUblasSparseSpace<TValue>::VectorType& rB,
    typename TUblasSparseSpace<TValue>::IndexType& rIterationNumber,
    TValue& rResidual,
    const boost::property_tree::ptree &amgclParams,
    int verbosity_level,
    bool use_gpgpu
    )
{
#ifdef AMGCL_GPGPU
    if (use_gpgpu && vexcl_context()) {
        auto &ctx = vexcl_context();

        typedef amgcl::backend::vexcl<TValue> Backend;
        typedef amgcl::make_solver<
            amgcl::runtime::preconditioner<Backend>,
            amgcl::runtime::solver::wrapper<Backend>
            > Solver;

        typename Backend::params bprm;
        bprm.q = ctx;

        Solver solve(amgcl::adapter::zero_copy(
                    TUblasSparseSpace<TValue>::Size1(rA),
                    rA.index1_data().begin(),
                    rA.index2_data().begin(),
                    rA.value_data().begin()),
                amgclParams, bprm);

        vex::vector<TValue> b(ctx, rB.size(), &rB[0]);
        vex::vector<TValue> x(ctx, rX.size(), &rX[0]);

        std::tie(rIterationNumber, rResidual) = solve(b, x);

        vex::copy(x.begin(), x.end(), rX.begin());

        if(verbosity_level > 1 )
            std::cout << "AMGCL Memory Occupation : " << amgcl::human_readable_memory(amgcl::backend::bytes(solve)) << std::endl;
    } else
#endif
    {
        typedef amgcl::backend::builtin<TValue> Backend;
        typedef amgcl::make_solver<
            amgcl::runtime::preconditioner<Backend>,
            amgcl::runtime::solver::wrapper<Backend>
            > Solver;

        Solver solve(amgcl::adapter::zero_copy(
                    TUblasSparseSpace<double>::Size1(rA),
                    rA.index1_data().begin(),
                    rA.index2_data().begin(),
                    rA.value_data().begin()),
                amgclParams);

        std::tie(rIterationNumber, rResidual) = solve(rB, rX);

        if(verbosity_level > 1 )
            std::cout << "AMGCL Memory Occupation : " << amgcl::human_readable_memory(amgcl::backend::bytes(solve)) << std::endl;
    }
}

template <class TValue, int TBlockSize>
void AMGCLBlockSolve(
    typename TUblasSparseSpace<TValue>::MatrixType & rA,
    typename TUblasSparseSpace<TValue>::VectorType& rX,
    typename TUblasSparseSpace<TValue>::VectorType& rB,
    typename TUblasSparseSpace<TValue>::IndexType& rIterationNumber,
    TValue& rResidual,
    boost::property_tree::ptree amgclParams,
    int verbosity_level,
    bool use_gpgpu
    )
{
    if(amgclParams.get<std::string>("precond.class") != "amg")
        amgclParams.erase("precond.coarsening");
    else
        amgclParams.put("precond.coarsening.aggr.block_size",1);

    typedef amgcl::static_matrix<TValue, TBlockSize, TBlockSize> value_type;
    typedef amgcl::static_matrix<TValue, TBlockSize, 1> rhs_type;

    std::size_t n = TUblasSparseSpace<TValue>::Size1(rA);
    std::size_t nb = n / TBlockSize;

#ifdef AMGCL_GPGPU
    if (use_gpgpu && vexcl_context()) {
        auto &ctx = vexcl_context();
        register_vexcl_static_matrix_type<TValue,TBlockSize>();

        typedef amgcl::backend::vexcl<value_type> Backend;

        typedef amgcl::make_solver<
            amgcl::runtime::preconditioner<Backend>,
            amgcl::runtime::solver::wrapper<Backend>
            > Solver;

        typename Backend::params bprm;
        bprm.q = ctx;

        Solver solve(
                amgcl::adapter::block_matrix<value_type>(
                    std::tie(n,rA.index1_data(),rA.index2_data(),rA.value_data() )),
                amgclParams, bprm);

        auto x_begin = reinterpret_cast<rhs_type*>(&rX[0]);
        auto b_begin = reinterpret_cast<const rhs_type*>(&rB[0]);

        vex::vector<rhs_type> x(ctx, nb, x_begin);
        vex::vector<rhs_type> b(ctx, nb, b_begin);

        std::tie(rIterationNumber, rResidual) = solve(b, x);

        vex::copy(x.begin(), x.end(), x_begin);

        if(verbosity_level > 1 )
            std::cout << "AMGCL Memory Occupation : " << amgcl::human_readable_memory(amgcl::backend::bytes(solve)) << std::endl;
    } else
#endif
    {
        typedef amgcl::backend::builtin<value_type> Backend;

        typedef amgcl::make_solver<
            amgcl::runtime::preconditioner<Backend>,
            amgcl::runtime::solver::wrapper<Backend>
            > Solver;

        Solver solve(
                amgcl::adapter::block_matrix<value_type>(
                    std::tie(n,rA.index1_data(),rA.index2_data(),rA.value_data() )),
                amgclParams);

        auto x_begin = reinterpret_cast<rhs_type*>(&rX[0]);
        boost::iterator_range<rhs_type*> x_range = boost::make_iterator_range(x_begin, x_begin + nb);

        auto b_begin = reinterpret_cast<const rhs_type*>(&rB[0]);
        boost::iterator_range<const rhs_type*> b_range = boost::make_iterator_range(b_begin, b_begin + nb);

        std::tie(rIterationNumber, rResidual) = solve(b_range, x_range);

        if(verbosity_level > 1 )
            std::cout << "AMGCL Memory Occupation : " << amgcl::human_readable_memory(amgcl::backend::bytes(solve)) << std::endl;
    }
}


namespace detail {


template <class TValue>
void AMGCLSolve(
    int block_size,
    typename TUblasSparseSpace<TValue>::MatrixType& rA,
    typename TUblasSparseSpace<TValue>::VectorType& rX,
    typename TUblasSparseSpace<TValue>::VectorType& rB,
    typename TUblasSparseSpace<TValue>::IndexType& rIterationNumber,
    TValue& rResidual,
    boost::property_tree::ptree amgclParams,
    int verbosity_level,
    bool use_gpgpu
    )
{
    if (use_gpgpu) {
        // ILU0 in a GPU backend has approximate iterative implementation.
        // Increase the default number of iterations to make ILU0 more robust.
        int ilu0_iters = 9;
        if (amgclParams.get<std::string>("precond.type", "") == "ilu0")
            amgclParams.put("precond.solve.iters", ilu0_iters);
        if (amgclParams.get<std::string>("precond.relax.type", "") == "ilu0")
            amgclParams.put("precond.relax.solve.iters", ilu0_iters);
    }

    switch (block_size) {
        case 2:
            AMGCLBlockSolve<TValue,2>(rA, rX, rB, rIterationNumber, rResidual, amgclParams, verbosity_level, use_gpgpu);
            return;
        case 3:
            AMGCLBlockSolve<TValue,3>(rA, rX, rB, rIterationNumber, rResidual, amgclParams, verbosity_level, use_gpgpu);
            return;
        case 4:
            AMGCLBlockSolve<TValue,4>(rA, rX, rB, rIterationNumber, rResidual, amgclParams, verbosity_level, use_gpgpu);
            return;
        default:
            AMGCLScalarSolve<TValue>(rA, rX, rB, rIterationNumber, rResidual, amgclParams, verbosity_level, use_gpgpu);
            return;
    }
}


} // namespace detail


void AMGCLSolve(
    int block_size,
    TUblasSparseSpace<double>::MatrixType& rA,
    TUblasSparseSpace<double>::VectorType& rX,
    TUblasSparseSpace<double>::VectorType& rB,
    TUblasSparseSpace<double>::IndexType& rIterationNumber,
    double& rResidual,
    boost::property_tree::ptree amgclParams,
    int verbosity_level,
    bool use_gpgpu
    )
{
    detail::AMGCLSolve<double>(block_size,
                               rA, rX, rB,
                               rIterationNumber, rResidual,
                               amgclParams, verbosity_level, use_gpgpu);
}


void AMGCLSolve(
    int block_size,
    TUblasSparseSpace<float>::MatrixType& rA,
    TUblasSparseSpace<float>::VectorType& rX,
    TUblasSparseSpace<float>::VectorType& rB,
    TUblasSparseSpace<float>::IndexType& rIterationNumber,
    float& rResidual,
    boost::property_tree::ptree amgclParams,
    int verbosity_level,
    bool use_gpgpu
    )
{
    detail::AMGCLSolve<float>(block_size,
                              rA, rX, rB,
                              rIterationNumber, rResidual,
                              amgclParams, verbosity_level, use_gpgpu);
}


}
