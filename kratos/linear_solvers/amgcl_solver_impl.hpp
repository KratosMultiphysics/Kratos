//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#pragma once

// The implementation of AMGCLSolver is split between
// - an implementation header (this file)
// - implementation sources (e.g.: "linear_solvers/amgcl_solver_impl.cpp")
//
// The reason is twofold:
// - includes from the AMGCL library are extremely heavy, so they are
//   avoided in the class declaration ("linear_solvers/amgcl_solver.h").
//   Instead, the implementation header includes them and defines logic
//   common to any matrix/vector representations. Each source file that
//   defines an instantiation of AMGCLSolver includes the implementation
//   header.
// - Shared memory and distributed memory matrix/vector representations
//   are handled in separate source files to avoid adding a Trilinos
//   dependency to core.

// External includes
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
#include <amgcl/coarsening/rigid_body_modes.hpp>

#ifdef AMGCL_GPGPU
#include <amgcl/backend/vexcl.hpp>
#include <amgcl/backend/vexcl_static_matrix.hpp>
#endif

#if defined(KRATOS_USING_MPI) && defined(KRATOS_AMGCL_MPI)
#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/make_solver.hpp>
#include <amgcl/mpi/preconditioner.hpp>
#include <amgcl/mpi/solver/runtime.hpp>
#endif

// Core includes
#include "includes/define.h"
#include "input_output/logger.h"
#include "linear_solvers/amgcl_solver.h"

// STL Includes
#include <variant> // std::variant

namespace Kratos {


template <class TSparse>
struct AMGCLAdaptor {};


#ifdef AMGCL_GPGPU
inline vex::Context& VexCLContext() {
    static vex::Context ctx(vex::Filter::Env);
    static bool run_once = [](){
        std::cout << "VexCL context:\n" << ctx << std::endl;
        return true;
    }();
    (void)run_once; // suppress "unused variable" warnings
    return ctx;
}

template <class TValue, int TBlockSize>
inline void register_vexcl_static_matrix_type() {
    static vex::scoped_program_header header(VexCLContext(),
            amgcl::backend::vexcl_static_matrix_declaration<TValue,TBlockSize>());
}
#endif


template <class TSparse, class TDense>
struct AMGCLSolver<TSparse,TDense>::Impl
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
    using MakeSharedMemorySolver = std::unique_ptr<amgcl::make_solver<
        amgcl::runtime::preconditioner<TBackend>,
        amgcl::runtime::solver::wrapper<TBackend>
    >>;

    #if defined(KRATOS_USING_MPI) && defined(KRATOS_AMGCL_MPI)
    template <class TBackend>
    using MakeDistributedSolver = std::unique_ptr<amgcl::mpi::make_solver<
        amgcl::runtime::mpi::preconditioner<TBackend>,
        amgcl::runtime::mpi::solver::wrapper<TBackend>
    >>;
    #endif

    std::conditional_t<
        TSparse::IsDistributed(),
        std::variant<
            #if defined(KRATOS_USING_MPI) && defined(KRATOS_AMGCL_MPI)
                 MakeDistributedSolver<amgcl::backend::builtin<BackendMatrix<1>>>
                ,MakeDistributedSolver<amgcl::backend::builtin<BackendMatrix<2>>>
                ,MakeDistributedSolver<amgcl::backend::builtin<BackendMatrix<3>>>
                ,MakeDistributedSolver<amgcl::backend::builtin<BackendMatrix<4>>>
                ,MakeDistributedSolver<amgcl::backend::builtin<BackendMatrix<5>>>
                ,MakeDistributedSolver<amgcl::backend::builtin<BackendMatrix<6>>>
                #ifdef AMGCL_GPGPU
                    ,MakeDistributedSolver<amgcl::backend::vexcl<BackendMatrix<1>>>
                    ,MakeDistributedSolver<amgcl::backend::vexcl<BackendMatrix<2>>>
                    ,MakeDistributedSolver<amgcl::backend::vexcl<BackendMatrix<3>>>
                    ,MakeDistributedSolver<amgcl::backend::vexcl<BackendMatrix<4>>>
                    ,MakeDistributedSolver<amgcl::backend::vexcl<BackendMatrix<5>>>
                    ,MakeDistributedSolver<amgcl::backend::vexcl<BackendMatrix<6>>>
                #endif
            #endif
        >,
        std::variant<
             MakeSharedMemorySolver<amgcl::backend::builtin<BackendMatrix<1>>>
            ,MakeSharedMemorySolver<amgcl::backend::builtin<BackendMatrix<2>>>
            ,MakeSharedMemorySolver<amgcl::backend::builtin<BackendMatrix<3>>>
            ,MakeSharedMemorySolver<amgcl::backend::builtin<BackendMatrix<4>>>
            ,MakeSharedMemorySolver<amgcl::backend::builtin<BackendMatrix<5>>>
            ,MakeSharedMemorySolver<amgcl::backend::builtin<BackendMatrix<6>>>
            #ifdef AMGCL_GPGPU
                ,MakeSharedMemorySolver<amgcl::backend::vexcl<BackendMatrix<1>>>
                ,MakeSharedMemorySolver<amgcl::backend::vexcl<BackendMatrix<2>>>
                ,MakeSharedMemorySolver<amgcl::backend::vexcl<BackendMatrix<3>>>
                ,MakeSharedMemorySolver<amgcl::backend::vexcl<BackendMatrix<4>>>
                ,MakeSharedMemorySolver<amgcl::backend::vexcl<BackendMatrix<5>>>
                ,MakeSharedMemorySolver<amgcl::backend::vexcl<BackendMatrix<6>>>
            #endif
        >
    > mpSolver;
}; // struct AMGCLSolver::Impl


template <class TSparse, class TDense>
AMGCLSolver<TSparse,TDense>::AMGCLSolver()
    : AMGCLSolver(Parameters())
{
}


template <class TSparse, class TDense>
AMGCLSolver<TSparse,TDense>::AMGCLSolver(Parameters Settings)
    : mTolerance(0),
      mMaxIterationsNumber(0),
      mVerbosity(0),
      mBlockSize(),
      mGMRESSize(0),
      mCoarseEnough(0),
      mFallbackToGMRES(false),
      mProvideCoordinates(false),
      mUseBlockMatricesIfPossible(false)
{
    KRATOS_TRY
    mpImpl.reset(new Impl);
    this->ApplySettings(Settings);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
AMGCLSolver<TSparse,TDense>::AMGCLSolver(const std::string& rSmootherName,
                                         const std::string& rSolverName,
                                         double Tolerance,
                                         int MaxIterationsNumber,
                                         int Verbosity,
                                         int GMRESSize)
    : mTolerance(Tolerance),
      mMaxIterationsNumber(MaxIterationsNumber),
      mVerbosity(Verbosity),
      mBlockSize(),
      mGMRESSize(GMRESSize),
      mCoarseEnough(1000),
      mFallbackToGMRES(false),
      mProvideCoordinates(false),
      mUseBlockMatricesIfPossible(false)
{
    KRATOS_TRY
    mpImpl.reset(new Impl);
    Parameters settings;
    settings.AddString("smoother_type", rSmootherName);
    settings.AddString("krylov_type", rSolverName);
    settings.AddString("coarsening_type", "aggregation");
    this->ApplySettings(settings);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
AMGCLSolver<TSparse,TDense>::AMGCLSolver(const std::string& rSmootherName,
                                         const std::string& rSolverName,
                                         const std::string& rCoarseningName,
                                         double Tolerance,
                                         int MaxIterationsNumber,
                                         int Verbosity,
                                         int GMRESSize,
                                         bool ProvideCoordinates)
    : mTolerance(Tolerance),
      mMaxIterationsNumber(MaxIterationsNumber),
      mVerbosity(Verbosity),
      mBlockSize(),
      mGMRESSize(GMRESSize),
      mCoarseEnough(1000),
      mFallbackToGMRES(false),
      mProvideCoordinates(ProvideCoordinates),
      mUseBlockMatricesIfPossible(false)
{
    KRATOS_TRY
    mpImpl.reset(new Impl);
    Parameters settings;
    settings.AddString("smoother_type", rSmootherName);
    settings.AddString("krylov_type", rSolverName);
    settings.AddString("coarsening_type", rCoarseningName);
    this->ApplySettings(settings);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
AMGCLSolver<TSparse,TDense>::AMGCLSolver(AMGCLSolver&&) noexcept = default;


template <class TSparse, class TDense>
AMGCLSolver<TSparse,TDense>::~AMGCLSolver() = default;


template <class TSparse, class TDense>
void AMGCLSolver<TSparse,TDense>::ApplySettings(Parameters Settings)
{
    Parameters default_parameters = this->GetDefaultParameters();

    // Optionally set the a user-defined block size.
    KRATOS_TRY
    if (Settings.Has("block_size")) {
        Parameters block_size_settings = Settings["block_size"].Clone();
        default_parameters.RemoveValue("block_size");

        if (block_size_settings.IsInt()) {
            mBlockSize = Settings["block_size"].GetInt();
            default_parameters.AddInt("block_size", 1);
        } else if (block_size_settings.IsString()) {
            KRATOS_ERROR_IF_NOT(block_size_settings.GetString() == "auto")
                << "Invalid value for \"block_size\": \"" << block_size_settings.GetString() << "\". "
                << "Expecting a positive integer, or \"auto\".";
            default_parameters.AddString("block_size", "auto");
        } else {
            KRATOS_ERROR << "Invalid type for \"block_size\". Expecting a positive integer, or \"auto\".";
        }
    }
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
    if (preconditioner_type == "relaxation")
        this->SetSmootherType(Settings["smoother_type"].GetString());

    mVerbosity=Settings["verbosity"].GetInt();
    mGMRESSize = Settings["gmres_krylov_space_dimension"].GetInt();
    mUseBlockMatricesIfPossible = Settings["use_block_matrices_if_possible"].GetBool();

    // Termination criteria.
    mTolerance = Settings["tolerance"].GetDouble();
    mMaxIterationsNumber = Settings["max_iteration"].GetInt();
    mAMGCLParameters.put("solver.tol", mTolerance);
    mAMGCLParameters.put("solver.maxiter", mMaxIterationsNumber);

    mFallbackToGMRES = false;
    this->SetIterativeSolverType(Settings["krylov_type"].GetString());

    // Multigrid preconditioner settings.
    mProvideCoordinates = Settings["provide_coordinates"].GetBool();
    mCoarseEnough = Settings["coarse_enough"].GetInt();
    mUseAMGPreconditioning = preconditioner_type == "amg";

    if(mUseAMGPreconditioning) {
        this->SetSmootherType(Settings["smoother_type"].GetString());
        this->SetCoarseningType(Settings["coarsening_type"].GetString());

        const int max_levels = Settings["max_levels"].GetInt();
        if(max_levels >= 0) mAMGCLParameters.put("precond.max_levels",  max_levels);

        mAMGCLParameters.put("precond.npre",  Settings["pre_sweeps"].GetInt());
        mAMGCLParameters.put("precond.npost",  Settings["post_sweeps"].GetInt());
    } // is mUseAMGPreconditioning

    // GPU settings.
    mUseGPGPU = Settings["use_gpgpu"].GetBool();
    if (mUseGPGPU) {
        // ILU0 in a GPU backend has approximate iterative implementation.
        // Increase the default number of iterations to make ILU0 more robust.
        int ilu0_iters = 9;
        if (mAMGCLParameters.get<std::string>("precond.type", "") == "ilu0")
            mAMGCLParameters.put("precond.solve.iters", ilu0_iters);
        if (mAMGCLParameters.get<std::string>("precond.relax.type", "") == "ilu0")
            mAMGCLParameters.put("precond.relax.solve.iters", ilu0_iters);
    }
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void AMGCLSolver<TSparse,TDense>::ProvideAdditionalData(SparseMatrixType& rA,
                                                        VectorType& rX,
                                                        VectorType& rB,
                                                        DofsArrayType& rDofSet,
                                                        ModelPart& rModelPart)
{
    int detected_block_size = 1;
    int old_ndof = -1;
    int ndof=0;

    if (!rModelPart.IsDistributed())
    {
        unsigned int old_node_id = rDofSet.size() ? rDofSet.begin()->Id() : 0;
        for (auto it = rDofSet.begin(); it!=rDofSet.end(); it++) {
            if(it->EquationId() < TSparse::Size1(rA) ) {
                IndexType id = it->Id();
                if(id != old_node_id) {
                    old_node_id = id;
                    if(old_ndof == -1) old_ndof = ndof;
                    else if(old_ndof != ndof) { //if it is different than the block size is 1
                        old_ndof = -1;
                        break;
                    }

                    ndof=1;
                } else {
                    ndof++;
                }
            }
        }

        if(old_ndof == -1 || old_ndof != ndof)
            detected_block_size = 1;
        else
            detected_block_size = ndof;

    }
    else //distribute
    {
        const std::size_t system_size = TSparse::Size1(rA);
        int current_rank = rModelPart.GetCommunicator().GetDataCommunicator().Rank();
        unsigned int old_node_id = rDofSet.size() ? rDofSet.begin()->Id() : 0;
        for (auto it = rDofSet.begin(); it!=rDofSet.end(); it++) {
            if(it->EquationId() < system_size  && it->GetSolutionStepValue(PARTITION_INDEX) == current_rank) {
                IndexType id = it->Id();
                if(id != old_node_id) {
                    old_node_id = id;
                    if(old_ndof == -1) old_ndof = ndof;
                    else if(old_ndof != ndof) { //if it is different than the block size is 1
                        old_ndof = -1;
                        break;
                    }

                    ndof=1;
                } else {
                    ndof++;
                }
            }
        }

        if(old_ndof != -1)
            detected_block_size = ndof;

        int max_block_size = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(detected_block_size);

        if( old_ndof == -1) {
            detected_block_size = max_block_size;
        }

        KRATOS_ERROR_IF(detected_block_size != max_block_size)
            << "Block size is not consistent. Local: " << detected_block_size
            << " Max: " << max_block_size << std::endl;
    }

    if (mBlockSize.has_value()) {
        KRATOS_WARNING_IF("AMGCLSolver", mBlockSize.value() != detected_block_size)
            << "The user-specified block size (" << mBlockSize.value() << ") "
            << "does not match the detected block size (" << detected_block_size << "). "
            << "Continuing with the user-specified block size (" << mBlockSize.value() << ").";
    } else {
        mBlockSize = detected_block_size;
    }

    KRATOS_INFO_IF("AMGCL Linear Solver", mVerbosity > 1)
        << "mndof: " << mBlockSize.value() << std::endl;

    KRATOS_ERROR_IF(TSparse::Size1(rA) % mBlockSize.value())
        << "LHS matrix (" << TSparse::Size1(rA) << "," << TSparse::Size2(rA) << ") "
        << "is incompatible with a \"block_size\" of " << mBlockSize.value();

    if(mProvideCoordinates) {
        mCoordinates.resize(TSparse::Size1(rA) / mBlockSize.value());
        unsigned int i=0;
        for (auto it_dof = rDofSet.begin(); it_dof!=rDofSet.end(); it_dof += mBlockSize.value()) {
            if(it_dof->EquationId() < TSparse::Size1(rA) ) {
                auto it_node = rModelPart.Nodes().find(it_dof->Id());
                mCoordinates[ i ] = it_node->Coordinates();
                i++;
            }
        }
    }
}


template <class TSparse, class TDense>
void AMGCLSolver<TSparse,TDense>::Clear()
{
    mpImpl->mpSolver = decltype(mpImpl->mpSolver)();
    mCoordinates = decltype(mCoordinates)();
}


template <class TSparse, class TDense>
void AMGCLSolver<TSparse,TDense>::InitializeSolutionStep(SparseMatrixType& rLhs,
                                                         VectorType& rSolution,
                                                         VectorType& rRhs)
{
    // Sanity checks.
    KRATOS_ERROR_IF(   TSparse::Size1(rLhs) != TSparse::Size2(rLhs)
                    || TSparse::Size(rSolution)  != TSparse::Size2(rLhs)
                    || TSparse::Size(rRhs) != TSparse::Size1(rLhs))
        << "inconsistent linear system shapes: "
        << "(" << TSparse::Size1(rLhs) << "x" << TSparse::Size2(rLhs) << ") "
        << "@"
        << "(" << TSparse::Size(rSolution) << "x1) = "
        << "(" << TSparse::Size(rRhs) << "x1)";

    // Default block size to 1 if it has not been previously set.
    // The only way this can happen is that the "block_size" setting
    // was not set during the construction of the solver, nor was
    // ProvideAdditionalData called.
    // Honestly, this should result in an error but some parts of the
    // code do it anyway and I don't have the patience to babysit PRs
    // fixing each and every one of them.
    // @matekelemen
    if (!mBlockSize.has_value()) {
        KRATOS_WARNING_IF("AMGCLSolver", 0 < mVerbosity)
            << "System solution requested without choosing a block size or calling AMGCLSolver::ProvideAdditionalData. "
            << "Defaulting to a block size of 1."
            << std::endl;
        mBlockSize = 1;
    }

    KRATOS_TRY
    // Apply settings that depend on the LHS matrix.
    // Use rigid body modes or set block size
    int static_block_size = mUseBlockMatricesIfPossible ? mBlockSize.value() : 1;
    std::vector<double> B;
    if(mUseAMGPreconditioning && mProvideCoordinates && (mBlockSize.value() == 2 || mBlockSize.value() == 3)) {
        int nmodes = amgcl::coarsening::rigid_body_modes(mBlockSize.value(),
                boost::make_iterator_range(
                    &mCoordinates[0][0],
                    &mCoordinates[0][0] + TSparse::Size1(rLhs)),
                B);

        if (static_block_size != 1 && static_block_size != 3) {
            KRATOS_WARNING("AMGCL Linear Solver") << "Can only combine block matrices with coordinates in 3D. Falling back to scalar matrices" << std::endl;
            static_block_size = 1;
        }

        mAMGCLParameters.put("precond.coarsening.aggr.eps_strong", 0.0);
        mAMGCLParameters.put("precond.coarsening.aggr.block_size", 1);
        mAMGCLParameters.put("precond.coarsening.nullspace.cols",  nmodes);
        mAMGCLParameters.put("precond.coarsening.nullspace.rows",  TSparse::Size1(rLhs));
        mAMGCLParameters.put("precond.coarsening.nullspace.B",     &B[0]);
    } else if(mUseAMGPreconditioning && mAMGCLParameters.get<std::string>("precond.coarsening.type") != std::string("ruge_stuben")) {
        mAMGCLParameters.put("precond.coarsening.aggr.eps_strong", 0.0);
        mAMGCLParameters.put("precond.coarsening.aggr.block_size", 1);
    }

    if(mUseAMGPreconditioning)
        mAMGCLParameters.put("precond.coarse_enough",mCoarseEnough / mBlockSize.value());

    if (mVerbosity > 2) {
        write_json(std::cout, mAMGCLParameters);
    }

    if(mVerbosity == 4) {
        //output to matrix market
        TSparse::WriteMatrixMarketMatrix("A.mm", rLhs, false);
        TSparse::WriteMatrixMarketVector("b.mm", rRhs);

        if(mProvideCoordinates) {
            //output of coordinates
            std::ofstream coordsfile;
            coordsfile.open ("coordinates.csv");
            for(unsigned int i=0; i<mCoordinates.size(); i++) {
                coordsfile << mCoordinates[i][0] << "," << mCoordinates[i][1] << "," << mCoordinates[i][2] << "\n";
            }
        }

        KRATOS_ERROR << " Verbosity = 4 prints the matrix and exits" << std::endl;
    }

    if(mUseBlockMatricesIfPossible) {
        KRATOS_ERROR_IF(TSparse::Size1(rLhs) % mBlockSize.value())
            << "The requested block size (" << mBlockSize.value() << ") "
            << "is not an exact multiple of the matrix size (" << TSparse::Size1(rLhs) << ")";
    }

    // Compute the multigrid hierarchy.
    if (mUseGPGPU) {
        #ifdef AMGCL_GPGPU

        // Initialize the static GPU environment.
        auto &vexcl_context = VexCLContext();
        KRATOS_ERROR_IF_NOT(vexcl_context) << "failed to initialize VexCL context";

        // Construct a GPU-bound shared-memory solver.
        #define KRATOS_MAKE_SHARED_MEMORY_AMGCL_SOLVER(BLOCK_SIZE)                                          \
            using BackendType = amgcl::backend::vexcl<typename Impl::template BackendMatrix<BLOCK_SIZE>>;   \
            using SolverType = typename Impl::template MakeSharedMemorySolver<BackendType>::element_type;   \
                                                                                                            \
            typename BackendType::params backend_parameters;                                                \
            backend_parameters.q = vexcl_context;                                                           \
                                                                                                            \
            KRATOS_TRY                                                                                      \
            mpImpl->mpSolver = std::unique_ptr<SolverType>(new SolverType(                                  \
                AMGCLAdaptor<TSparse>().template MakeMatrixAdaptor<BLOCK_SIZE>(rLhs),                       \
                mAMGCLParameters,                                                                           \
                backend_parameters));                                                                       \
            KRATOS_CATCH("")

        // Construct a GPU-bound distributed-memory solver.
        #ifdef KRATOS_USING_MPI
            #define KRATOS_MAKE_DISTRIBUTED_AMGCL_SOLVER(BLOCK_SIZE)                                            \
                using BackendType = amgcl::backend::vexcl<typename Impl::template BackendMatrix<BLOCK_SIZE>>;   \
                using SolverType = typename Impl::template MakeDistributedSolver<BackendType>::element_type;    \
                                                                                                                \
                typename BackendType::params backend_parameters;                                                \
                backend_parameters.q = vexcl_context;                                                           \
                                                                                                                \
                KRATOS_TRY                                                                                      \
                mpImpl->mpSolver = std::unique_ptr<SolverType>(new SolverType(                                  \
                    AMGCLAdaptor<TSparse>().GetCommunicator(rLhs),                                              \
                    AMGCLAdaptor<TSparse>().template MakeMatrixAdaptor<BLOCK_SIZE>(rLhs),                       \
                    mAMGCLParameters,                                                                           \
                    backend_parameters));                                                                       \
                KRATOS_CATCH("")
        #else
            #define KRATOS_MAKE_DISTRIBUTED_AMGCL_SOLVER(BLOCK_SIZE)                                            \
                KRATOS_ERROR << "Requesting a distributed memory AMGCL solver but Kratos was compiled "         \
                             << "without MPI support.";
        #endif

        #define KRATOS_MAKE_AMGCL_SOLVER(BLOCK_SIZE)                                                        \
            if constexpr (TSparse::IsDistributed()) {                                                       \
                KRATOS_MAKE_DISTRIBUTED_AMGCL_SOLVER(BLOCK_SIZE)                                            \
            } else {                                                                                        \
                KRATOS_MAKE_SHARED_MEMORY_AMGCL_SOLVER(BLOCK_SIZE)                                          \
            }

        switch (mBlockSize.value()) {
            case 1: {KRATOS_MAKE_AMGCL_SOLVER(1) break;}
            case 2: {KRATOS_MAKE_AMGCL_SOLVER(2) break;}
            case 3: {KRATOS_MAKE_AMGCL_SOLVER(3) break;}
            case 4: {KRATOS_MAKE_AMGCL_SOLVER(4) break;}
            case 5: {KRATOS_MAKE_AMGCL_SOLVER(5) break;}
            case 6: {KRATOS_MAKE_AMGCL_SOLVER(6) break;}
            default: KRATOS_ERROR << "Unsupported block size for AMGCLSolver: " << mBlockSize.value() << ". "
                                  << "Options are 1, 2, 3, 4, 5, 6.";
        } // switch mBlockSize

        #undef KRATOS_MAKE_AMGCL_SOLVER
        #undef KRATOS_MAKE_DISTRIBUTED_AMGCL_SOLVER
        #undef KRATOS_MAKE_SHARED_MEMORY_AMGCL_SOLVER

        #else
        KRATOS_ERROR << "Requested a GPU-bound AMGCL solver, but Kratos was compiled without VexCL support!";
        #endif
    } /*if mUseGPGPU*/ else {
        // Construct a CPU-bound shared-memory solver.
        #define KRATOS_MAKE_SHARED_MEMORY_AMGCL_SOLVER(BLOCK_SIZE)                                          \
            using BackendType = amgcl::backend::builtin<typename Impl::template BackendMatrix<BLOCK_SIZE>>; \
            using SolverType = typename Impl::template MakeSharedMemorySolver<BackendType>::element_type;   \
                                                                                                            \
            KRATOS_TRY                                                                                      \
            mpImpl->mpSolver = std::unique_ptr<SolverType>(new SolverType(                                  \
                AMGCLAdaptor<TSparse>().template MakeMatrixAdaptor<BLOCK_SIZE>(rLhs),                       \
                mAMGCLParameters));                                                                         \
            KRATOS_CATCH("")

        // Construct a CPU-bound distributed-memory solver.
        #ifdef KRATOS_USING_MPI
            #define KRATOS_MAKE_DISTRIBUTED_AMGCL_SOLVER(BLOCK_SIZE)                                            \
                using BackendType = amgcl::backend::builtin<typename Impl::template BackendMatrix<BLOCK_SIZE>>; \
                using SolverType = typename Impl::template MakeDistributedSolver<BackendType>::element_type;    \
                                                                                                                \
                KRATOS_TRY                                                                                      \
                mpImpl->mpSolver = std::unique_ptr<SolverType>(new SolverType(                                  \
                    AMGCLAdaptor<TSparse>().GetCommunicator(rLhs),                                              \
                    AMGCLAdaptor<TSparse>().template MakeMatrixAdaptor<BLOCK_SIZE>(rLhs),                       \
                    mAMGCLParameters));                                                                         \
                KRATOS_CATCH("")
        #else
            #define KRATOS_MAKE_DISTRIBUTED_AMGCL_SOLVER(BLOCK_SIZE)                                            \
                KRATOS_ERROR << "Requesting a distributed memory AMGCL solver but Kratos was compiled "         \
                             << "without MPI support.";
        #endif

        #define KRATOS_MAKE_AMGCL_SOLVER(BLOCK_SIZE)                                                        \
            if constexpr (TSparse::IsDistributed()) {                                                       \
                KRATOS_MAKE_DISTRIBUTED_AMGCL_SOLVER(BLOCK_SIZE)                                            \
            } else {                                                                                        \
                KRATOS_MAKE_SHARED_MEMORY_AMGCL_SOLVER(BLOCK_SIZE)                                          \
            }

        switch (mBlockSize.value()) {
            case 1: {KRATOS_MAKE_AMGCL_SOLVER(1) break;}
            case 2: {KRATOS_MAKE_AMGCL_SOLVER(2) break;}
            case 3: {KRATOS_MAKE_AMGCL_SOLVER(3) break;}
            case 4: {KRATOS_MAKE_AMGCL_SOLVER(4) break;}
            case 5: {KRATOS_MAKE_AMGCL_SOLVER(5) break;}
            case 6: {KRATOS_MAKE_AMGCL_SOLVER(6) break;}
            default: KRATOS_ERROR << "Unsupported block size for AMGCLSolver: " << mBlockSize.value() << ". "
                                  << "Options are 1, 2, 3, 4, 5, 6.";
        } // switch mBlockSize

        #undef KRATOS_MAKE_AMGCL_SOLVER
        #undef KRATOS_MAKE_DISTRIBUTED_AMGCL_SOLVER
        #undef KRATOS_MAKE_SHARED_MEMORY_AMGCL_SOLVER
    } /*if mUseGPGPU else*/

    // Issue message about memory footprint.
    if (1 < mVerbosity) {
        std::visit(
            [] (const auto& rp_solver) -> void {
                KRATOS_INFO("AMGCLSolver")
                    << amgcl::human_readable_memory(amgcl::backend::bytes(*rp_solver))
                    << " memory footprint"
                    << std::endl;
            },
            mpImpl->mpSolver
        );
    }
    KRATOS_CATCH("")
}


template <class T>
struct AMGCLStaticVectorTraits {};


template <>
struct AMGCLStaticVectorTraits<double>
{
    static constexpr int value = 1;
    using type = double;
};


template <>
struct AMGCLStaticVectorTraits<float>
{
    static constexpr int value = 1;
    using type = float;
};


template <class TValue, int RowCount, int ColumnCount>
struct AMGCLStaticVectorTraits<amgcl::static_matrix<TValue,RowCount,ColumnCount>>
{
    static constexpr int value = RowCount;
    using type = amgcl::static_matrix<TValue,RowCount,1>;
};


template <class TSparse, class TDense>
bool AMGCLSolver<TSparse,TDense>::PerformSolutionStep(SparseMatrixType& rLhs,
                                                      VectorType& rSolution,
                                                      VectorType& rRhs)
{
    KRATOS_TRY

    const auto [iteration_count, residual_norm] = std::visit(
        [&rSolution, &rRhs, &rLhs] (auto& rp_solver) -> std::pair<std::size_t,typename TSparse::DataType> {
            using ElementType = std::remove_reference_t<decltype(rp_solver)>;
            KRATOS_ERROR_IF_NOT(rp_solver) << "AMGCL solver is uninitialized";

            auto& r_solver = *rp_solver;
            using SolverType = typename std::pointer_traits<ElementType>::element_type;
            using BackendType = typename SolverType::backend_type;
            using StaticMatrixType = typename BackendType::value_type;
            using StaticVectorType = typename AMGCLStaticVectorTraits<StaticMatrixType>::type;

            if constexpr (std::is_same_v<BackendType,amgcl::backend::builtin<StaticMatrixType>>) {
                const std::size_t block_system_size = AMGCLAdaptor<TSparse>().template BlockSystemSize<StaticMatrixType>(rLhs);

                auto it_solution_begin = reinterpret_cast<StaticVectorType*>(AMGCLAdaptor<TSparse>().MakeVectorIterator(rSolution));
                const auto it_rhs_begin = reinterpret_cast<StaticVectorType*>(AMGCLAdaptor<TSparse>().MakeVectorIterator(rRhs));
                const auto [iteration_count, residual_norm] = r_solver(
                    boost::make_iterator_range(it_rhs_begin, it_rhs_begin + block_system_size),
                    boost::make_iterator_range(it_solution_begin, it_solution_begin + block_system_size));
                return std::pair<std::size_t,typename TSparse::DataType>(iteration_count, residual_norm);

            #ifdef AMGCL_GPGPU
            } else if constexpr (std::is_same_v<BackendType,amgcl::backend::vexcl<StaticMatrixType>>) {
                using SolverType = typename std::pointer_traits<ElementType>::element_type;

                auto& vexcl_context = VexCLContext();
                KRATOS_ERROR_IF_NOT(vexcl_context) << "invalid VexCL context";

                using StaticMatrixType = typename SolverType::backend_type::value_type;
                using StaticVectorType = typename AMGCLStaticVectorTraits<StaticMatrixType>::type;
                const std::size_t block_system_size = AMGCLAdaptor<TSparse>().template BlockSystemSize<StaticMatrixType>(rLhs);

                auto it_solution_begin = reinterpret_cast<StaticVectorType*>(AMGCLAdaptor<TSparse>().MakeVectorIterator(rSolution));
                const auto it_rhs_begin = reinterpret_cast<const StaticVectorType*>(AMGCLAdaptor<TSparse>().MakeVectorIterator(rRhs));
                vex::vector<StaticVectorType> solution(vexcl_context, block_system_size, it_solution_begin),
                                              rhs     (vexcl_context, block_system_size, it_rhs_begin);

                const auto [iteration_count, residual_norm] = r_solver(rhs, solution);
                vex::copy(solution.begin(), solution.end(), it_solution_begin);
                return std::pair<std::size_t,typename TSparse::DataType>(iteration_count, residual_norm);
            #endif

            } else {
                static_assert(std::is_same_v<ElementType,std::monostate>, "unhandled solver type");
                return {};
            }
        },
        mpImpl->mpSolver
    );

    // Failed to solve the system; retry with GMRES.
    // This results in a lasting changes to the original
    // configuration of this solver instance.
    // (this functionality should be handled outside this class)
    if (mFallbackToGMRES && mTolerance < residual_norm) {
        const std::string iterative_solver_name = mAMGCLParameters.get<std::string>("solver.type");
        if (iterative_solver_name != "gmres") {
            KRATOS_INFO("AMGCLSolver")
                << "Failed to solve the system using " << iterative_solver_name << " "
                << "(" << residual_norm << " residual in " << iteration_count << " iterations) "
                << ". Falling back to \"gmres\".";

            // Override the user-provided configuration and retry with
            // GMRES and a block size of 1.
            // => requires a refactorization.
            mAMGCLParameters.put("solver.type", "gmres");
            mAMGCLParameters.put("solver.M",  mGMRESSize);
            this->mBlockSize = 1;
            this->InitializeSolutionStep(rLhs, rSolution, rRhs);

            return this->PerformSolutionStep(rLhs, rSolution, rRhs);
        }
    }

    return residual_norm < mTolerance;

    KRATOS_CATCH("")
}


} // namespace Kratos
