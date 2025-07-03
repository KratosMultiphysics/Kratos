/* AMGCL */
#include "includes/define.h"
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

// STL Includes
#include <variant> // std::variant


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


template <class TSparse, class TDense, class TReorderer>
struct AMGCLSolver<TSparse,TDense,TReorderer>::Impl
{
    using Scalar = typename TSparse::DataType;

    template <int StaticSize>
    using BackendMatrix = std::conditional_t<
        StaticSize == 1,
        Scalar,
        amgcl::static_matrix<Scalar,StaticSize,StaticSize>
    >;

    template <int StaticSize>
    using BackendVector = std::conditional_t<
        StaticSize == 1,
        Scalar,
        amgcl::static_matrix<Scalar,StaticSize,1>
    >;

    template <class TBackend>
    using MakeSolver = amgcl::make_solver<
        amgcl::runtime::preconditioner<TBackend>,
        amgcl::runtime::solver::wrapper<TBackend>
    >;

    std::variant<
         std::monostate
        ,MakeSolver<amgcl::backend::builtin<BackendMatrix<1>>>
        ,MakeSolver<amgcl::backend::builtin<BackendMatrix<2>>>
        ,MakeSolver<amgcl::backend::builtin<BackendMatrix<3>>>
        ,MakeSolver<amgcl::backend::builtin<BackendMatrix<4>>>
        ,MakeSolver<amgcl::backend::builtin<BackendMatrix<5>>>
        ,MakeSolver<amgcl::backend::builtin<BackendMatrix<6>>>
        #ifdef AMGCL_GPGPU
        ,MakeSolver<amgcl::backend::vexcl<BackendMatrix<1>>>
        ,MakeSolver<amgcl::backend::vexcl<BackendMatrix<2>>>
        ,MakeSolver<amgcl::backend::vexcl<BackendMatrix<3>>>
        ,MakeSolver<amgcl::backend::vexcl<BackendMatrix<4>>>
        ,MakeSolver<amgcl::backend::vexcl<BackendMatrix<5>>>
        ,MakeSolver<amgcl::backend::vexcl<BackendMatrix<6>>>
        #endif
    > mSolver;
}; // struct AMGCLSolver::Impl


template <class TSparse, class TDense, class TReorderer>
AMGCLSolver<TSparse,TDense,TReorderer>::AMGCLSolver()
    : AMGCLSolver(Parameters())
{
}


template <class TSparse, class TDense, class TReorderer>
AMGCLSolver<TSparse,TDense,TReorderer>::AMGCLSolver(Parameters Settings)
{
    KRATOS_TRY
    mpImpl.emplace(new Impl);
    this->ApplySettings(Settings);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TReorderer>
AMGCLSolver<TSparse,TDense,TReorderer>::AMGCLSolver(const std::string& rSmootherName,
                                                    const std::string& rSolverName,
                                                    double Tolerance,
                                                    int MaxIterationsNumber,
                                                    int Verbosity,
                                                    int GMRESSize)
    : mTolerance(Tolerance),
      mMaxIterationsNumber(MaxIterationsNumber),
      mVerbosity(Verbosity),
      mBlockSize(1),
      mGMRESSize(GMRESSize),
      mCoarseEnough(1000),
      mFallbackToGMRES(false),
      mProvideCoordinates(false),
      mUseBlockMatricesIfPossible(false)
{
    KRATOS_TRY
    mpImpl.emplace(new Impl);
    Parameters settings;
    settings.AddString("smoother_type", rSmootherName);
    settings.AddString("krylov_type", rSolverName);
    settings.AddString("coarsening_type", "aggregation");
    this->ApplySettings(settings);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TReorderer>
AMGCLSolver<TSparse,TDense,TReorderer>::AMGCLSolver(AMGCLSmoother Smoother,
                                                    AMGCLIterativeSolverType Solver,
                                                    AMGCLCoarseningType Coarsening,
                                                    double Tolerance,
                                                    int MaxIterationsNumber,
                                                    int Verbosity,
                                                    int GMRESSize,
                                                    bool ProvideCoordinates)
    : mTolerance(Tolerance),
      mMaxIterationsNumber(MaxIterationsNumber),
      mVerbosity(Verbosity),
      mBlockSize(1),
      mGMRESSize(GMRESSize),
      mCoarseEnough(1000),
      mFallbackToGMRES(false),
      mProvideCoordinates(ProvideCoordinates),
      mUseBlockMatricesIfPossible(false)
{
    KRATOS_TRY
    mpImpl.emplace(new Impl);
    KRATOS_INFO_IF("AMGCL Linear Solver", mVerbosity > 0) << "Setting up AMGCL for iterative solve " << std::endl;

    // Choose smoother in the list "gauss_seidel, multicolor_gauss_seidel, ilu0, parallel_ilu0, ilut, damped_jacobi, spai0, chebyshev"
    SetSmootherType(Smoother);
    // Setting iterative solver
    SetIterativeSolverType(Solver);
    // Setting coarsening type
    SetCoarseningType(Coarsening);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TReorderer>
AMGCLSolver<TSparse,TDense,TReorderer>::AMGCLSolver(AMGCLSolver&&) noexcept = default;


template <class TSparse, class TDense, class TReorderer>
AMGCLSolver<TSparse,TDense,TReorderer>::~AMGCLSolver() = default;


template <class TSparse, class TDense, class TReorderer>
void AMGCLSolver<TSparse,TDense,TReorderer>::InitializeSolutionStep(SparseMatrixType& rLhs,
                                                                    VectorType& rSolution,
                                                                    VectorType& rRhs)
{
    KRATOS_TRY
    if (mUseGPGPU) {
        // ILU0 in a GPU backend has approximate iterative implementation.
        // Increase the default number of iterations to make ILU0 more robust.
        int ilu0_iters = 9;
        if (amgclParams.get<std::string>("precond.type", "") == "ilu0")
            amgclParams.put("precond.solve.iters", ilu0_iters);
        if (amgclParams.get<std::string>("precond.relax.type", "") == "ilu0")
            amgclParams.put("precond.relax.solve.iters", ilu0_iters);
    }

    #ifdef AMGCL_GPGPU
    if (mUseGPGPU && vexcl_context()) {
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
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TReorderer>
void AMGCLSolver<TSparse,TDense,TReorderer>::ApplySettings(Parameters Settings)
{
    Parameters default_parameters = this->GetDefaultParameters();

    KRATOS_TRY
    Settings.ValidateAndAssignDefaults(default_parameters);
    KRATOS_CATCH("")

    KRATOS_TRY
    CheckIfSelectedOptionIsAvailable(Settings,
                                     "preconditioner_type",
                                     {
                                        "amg",
                                        "relaxation",
                                        "dummy"
                                     });
    KRATOS_CATCH("")

    KRATOS_TRY
    const std::string preconditioner_type = Settings["preconditioner_type"].GetString();
    mAMGCLParameters.put("precond.class", preconditioner_type);
    mUseAMGPreconditioning = preconditioner_type == "amg";
    if(preconditioner_type == "relaxation") {
        mAMGCLParameters.put("precond.type", Settings["smoother_type"].GetString());
    }

    mProvideCoordinates = Settings["provide_coordinates"].GetBool();
    mCoarseEnough = Settings["coarse_enough"].GetInt();

    mBlockSize = Settings["block_size"].GetInt(); //set the mndof to an inital number
    mTolerance = Settings["tolerance"].GetDouble();
    mMaxIterationsNumber = Settings["max_iteration"].GetInt();
    mVerbosity=Settings["verbosity"].GetInt();
    mGMRESSize = Settings["gmres_krylov_space_dimension"].GetInt();

    mFallbackToGMRES = false;
    this->SetIterativeSolverType(Settings["krylov_type"].GetString())

    //settings only needed if full AMG is used
    if(mUseAMGPreconditioning)
    {
        this->SetSmootherType(Settings["smoother_type"].GetString());
        this->SetCoarseningType(Settings["coarsening_type"].GetString());

        int max_levels = Settings["max_levels"].GetInt();
        if(max_levels >= 0)
            mAMGCLParameters.put("precond.max_levels",  max_levels);

        mAMGCLParameters.put("precond.npre",  Settings["pre_sweeps"].GetInt());
        mAMGCLParameters.put("precond.npost",  Settings["post_sweeps"].GetInt());
    }

    mUseBlockMatricesIfPossible = Settings["use_block_matrices_if_possible"].GetBool();

    mUseGPGPU = Settings["use_gpgpu"].GetBool();
    KRATOS_CATCH("")
}


} // namespace Kratos
