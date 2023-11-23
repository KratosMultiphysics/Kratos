//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#ifndef AMGCL_PARAM_UNKNOWN
    #include "input_output/logger.h"
    #define AMGCL_PARAM_UNKNOWN(NAME)                       \
        KRATOS_ERROR                                        \
            << KRATOS_CODE_LOCATION                         \
            << "Unknown parameter " << NAME << std::endl
#endif

// Project includes
#include "amgcl_raw_solver.h"
#include "spaces/ublas_space.h"
#include "utilities/profiler.h"

// External includes
#include "amgcl/adapter/ublas.hpp"
#include "amgcl/adapter/zero_copy.hpp"
#include "amgcl/backend/builtin.hpp"
#include "amgcl/value_type/static_matrix.hpp"
#include "amgcl/make_solver.hpp"
#include "amgcl/make_block_solver.hpp"
#include "amgcl/solver/runtime.hpp"
#include "amgcl/preconditioner/runtime.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

#ifdef AMGCL_GPGPU
#include <amgcl/backend/vexcl.hpp>
#include <amgcl/backend/vexcl_static_matrix.hpp>
#endif

// STL includes
#include <sstream> // stringstream
#include <optional> // optional
#include <variant> // variant


namespace Kratos {


#ifdef AMGCL_GPGPU
vex::Context& GetVexCLContext() {
    static vex::Context ctx(vex::Filter::Env);
    [[maybe_unused]] static bool run_once = [](){
        std::cout << "VexCL context:\n" << ctx << std::endl;
        return true;
    }();
    return ctx;
}

template <int TBlockSize>
void RegisterVexCLStaticMatrixType() {
    static vex::scoped_program_header header(GetVexCLContext(),
            amgcl::backend::vexcl_static_matrix_declaration<double,TBlockSize>());
}
#endif


// -------------------------------------
// --- Template Nightmare Land Begin ---
// -------------------------------------


namespace Detail {


template <class TValue, unsigned BlockSize>
using AMGCLBlock = amgcl::static_matrix<TValue,BlockSize,BlockSize>;


// Trait class defining an AMGCL backend for the provided template arguments,
// and its associated preconditioner, wrapper, and solver.
template <unsigned BlockSize>
struct AMGCLTraits
{
    template <template <class, class ...> class TBackend, class TValue, class ...TArgs>
    struct Impl
    {
        // Indicates whether the system matrix and vectors must be transferred
        // to/from another device (GPU) before/after invoking the solver.
        constexpr static inline bool IsGPUBound =
        #ifdef AMGCL_GPGPU
            std::is_same_v<TBackend<TValue,TArgs...>,amgcl::backend::vexcl<TValue,TArgs...>>;
        #else
            false;
        #endif

        constexpr static inline bool IsScalar = BlockSize == 1;

        using Scalar = TValue;

        using Value = std::conditional_t<
            IsScalar,
            Scalar,
            amgcl::static_matrix<Scalar,BlockSize,BlockSize>
        >;

        using RHSValue = std::conditional_t<
            IsScalar,
            Scalar,
            amgcl::static_matrix<TValue,BlockSize,1>
        >;

        using MatrixAdapter = std::conditional_t<
            IsGPUBound,
            amgcl::adapter::block_matrix_adapter<typename TUblasSparseSpace<Scalar>::MatrixType,Value>,
            amgcl::backend::crs<Scalar>
        >;

        using Backend = std::conditional_t<
            IsScalar,
            TBackend<TValue>,
            TBackend<AMGCLBlock<TValue,BlockSize>>
        >;

        using Preconditioner = amgcl::runtime::preconditioner<Backend>;

        using SolverWrapper = amgcl::runtime::solver::wrapper<Backend>;

        using Solver = std::conditional_t<
            IsScalar,
            amgcl::make_solver<Preconditioner,SolverWrapper>,
            amgcl::make_block_solver<Preconditioner,SolverWrapper>
        >;
    }; // struct Impl

    template <template <class, class ...> class TBackend, class TValue, class ...TArgs>
    using Solver = typename Impl<TBackend,TValue,TArgs...>::Solver;
}; // struct AMGCLTraits


// Wrapper class bundling a solver with an associated matrix adapter.
template <class TAMGCLTraits>
struct AMGCLBundle
{
    using Traits = TAMGCLTraits;

    std::unique_ptr<typename Traits::Solver> mpSolver;

    std::shared_ptr<typename Traits::MatrixAdapter> mpMatrixAdapter;
}; // struct AMGCLBundle


} // namespace Detail


// -------------------------------------
// --- Template Nightmare Land End -----
// -------------------------------------


enum class AMGCLBackendType
{
    CPU
    #ifdef AMGCL_GPGPU
    ,SinglePrecisionGPU
    ,DoublePrecisionGPU
    #endif
}; // enum class VexCLBackendType


template <class TSparseSpace,
          class TDenseSpace,
          class TReorderer>
struct AMGCLRawSolver<TSparseSpace,TDenseSpace,TReorderer>::Impl
{
    int mVerbosity;

    double mTolerance;

    // AMGCL backend selected by the user.
    AMGCLBackendType mBackendType = AMGCLBackendType::CPU;

    // Block size computed in AMGCLRawSolver::ProvideAdditionalData and
    // set in the settings passed to AMGCL.
    std::size_t mDoFCount;

    // Settings to be translated from Kratos::Parameters
    // to something that AMGCL uses.
    boost::property_tree::ptree mAMGCLSettings;

    // Pointer to the system matrix to be set in AMGCLRawSolver::ProvideAdditionalData.
    // This matrix is then required to reside at the same address when invoking
    // AMGCLRawSolver::Solve.
    const typename TSparseSpace::MatrixType* mpA = nullptr;

    // A variant for grouping members related to the wrapped AMGCL solver. It's meant
    // to bundle solvers using different backends as well as its associated matrix wrapper.
    // The wrapped types are composed of the permutations of the following attributes:
    // - block size [1, 2, 3, 4]
    // - backend type [bultin, vexcl]
    // - value type [double, float (vexcl only)]
    std::variant<
        // Double precision CPU backend with a block size of 1.
        Detail::AMGCLBundle<typename Detail::AMGCLTraits<1>::Impl<amgcl::backend::builtin,double>>,

        // Double precision CPU backend with a block size of 2.
        Detail::AMGCLBundle<typename Detail::AMGCLTraits<2>::Impl<amgcl::backend::builtin,double>>,

        // Double precision CPU backend with a block size of 3.
        Detail::AMGCLBundle<typename Detail::AMGCLTraits<3>::Impl<amgcl::backend::builtin,double>>,

        // Double precision CPU backend with a block size of 4.
        Detail::AMGCLBundle<typename Detail::AMGCLTraits<4>::Impl<amgcl::backend::builtin,double>>,

        // Double precision GPU backend with a block size of 1.
        Detail::AMGCLBundle<typename Detail::AMGCLTraits<1>::Impl<amgcl::backend::vexcl,double>>,

        // Double precision GPU backend with a block size of 2.
        Detail::AMGCLBundle<typename Detail::AMGCLTraits<2>::Impl<amgcl::backend::vexcl,double>>,

        // Double precision GPU backend with a block size of 3.
        Detail::AMGCLBundle<typename Detail::AMGCLTraits<3>::Impl<amgcl::backend::vexcl,double>>,

        // Double precision GPU backend with a block size of 4.
        Detail::AMGCLBundle<typename Detail::AMGCLTraits<4>::Impl<amgcl::backend::vexcl,double>>,

        // Double precision GPU backend with a block size of 1.
        Detail::AMGCLBundle<typename Detail::AMGCLTraits<1>::Impl<amgcl::backend::vexcl,float>>,

        // Double precision GPU backend with a block size of 2.
        Detail::AMGCLBundle<typename Detail::AMGCLTraits<2>::Impl<amgcl::backend::vexcl,float>>,

        // Double precision GPU backend with a block size of 3.
        Detail::AMGCLBundle<typename Detail::AMGCLTraits<3>::Impl<amgcl::backend::vexcl,float>>,

        // Double precision GPU backend with a block size of 4.
        Detail::AMGCLBundle<typename Detail::AMGCLTraits<4>::Impl<amgcl::backend::vexcl,float>>
    > mSolverVariant;
}; // struct AMGCLRawSolver::Impl



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
AMGCLRawSolver<TSparseSpace,TDenseSpace,TReorderer>::AMGCLRawSolver(Parameters parameters)
    : mpImpl(new Impl)
{
    KRATOS_TRY
    Parameters default_parameters = this->GetDefaultParameters();
    parameters.ValidateAndAssignDefaults(default_parameters);

    KRATOS_ERROR_IF_NOT(parameters["solver_type"].GetString() == "amgcl_raw")
        << "Requested a(n) '" << parameters["solver_type"].GetString() << "' solver,"
        << " but constructing an AMGCLRawSolver";

    mpImpl->mVerbosity = parameters["verbosity"].Get<int>();
    mpImpl->mTolerance = parameters["tolerance"].Get<double>();

    // Get requested backend type. Supported arguments:
    // - ""         : empty string means no GPU backend => default to the built-in CPU backend
    // - "float"    : single precision floating point GPU backend
    // - "double"   : double precision floating point GPU backend
    const std::string requested_backend_type = parameters["gpgpu_backend"].GetString();
    if (requested_backend_type.empty()) {
        // CPU backend by default
        mpImpl->mBackendType = AMGCLBackendType::CPU;
    } else if (requested_backend_type == "float") {
        #ifdef AMGCL_GPGPU
            mpImpl->mBackendType = AMGCLBackendType::SinglePrecisionGPU;
        #else
            KRATOS_ERROR << "the requested gpgpu backend 'float' is not available because Kratos was compiled without GPU support\n";
        #endif
    } else if (requested_backend_type == "double") {
        #ifdef AMGCL_GPGPU
            mpImpl->mBackendType = AMGCLBackendType::DoublePrecisionGPU;
        #else
            KRATOS_ERROR << "the requested gpgpu backend 'double' is not available because Kratos was compiled without GPU support\n";
        #endif
    } else {
        KRATOS_ERROR << "unsupported argument for 'gpgpu_backend': "
                     << requested_backend_type
                     << ". Available options are '', 'float', or 'double'.\n";
    }

    // Convert parameters to AMGCL settings
    std::stringstream json_stream;
    json_stream << parameters["amgcl_settings"].PrettyPrintJsonString();
    boost::property_tree::read_json(
        json_stream,
        mpImpl->mAMGCLSettings
    );
    KRATOS_CATCH("")
}



// Necessary for PIMPL
template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
AMGCLRawSolver<TSparseSpace,TDenseSpace,TReorderer>::~AMGCLRawSolver()
{
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool AMGCLRawSolver<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                Vector& rX,
                                                                Vector& rB)
{
    KRATOS_TRY

    // Override AMGCL settings
    //mpImpl->mAMGCLSettings.put("solver.verbose", 1 < mpImpl->mVerbosity);

    KRATOS_ERROR_IF_NOT(&rA == mpImpl->mSolverStruct.mpMatrix)
        << "solver got a different matrix than it was initialized with "
        << &rA << " != " << mpImpl->mSolverStruct.mpMatrix << "\n";

    const auto [iteration_count, residual] = std::visit(
        [&rB, &rX] (auto& rpSolver) -> std::tuple<std::size_t,double> {
            #ifdef AMGCL_GPGPU
            if constexpr (std::is_same_v<
                            typename std::remove_reference_t<decltype(rpSolver)>::element_type,
                            VexCLDoublePrecisionSolver
                          >) {
                auto& vexcl_context = GetVexCLContext();
                KRATOS_ERROR_IF_NOT(vexcl_context) << "invalid VexCL context state\n";

                vex::vector<double> x(vexcl_context, rX.size(), &*rX.begin());
                vex::vector<double> b(vexcl_context, rB.size(), &*rB.begin());
                const auto results = rpSolver->operator()(b, x);

                vex::copy(x.begin(), x.end(), rX.begin());
                return results;
            } else if constexpr(std::is_same_v<
                            typename std::remove_reference_t<decltype(rpSolver)>::element_type,
                            VexCLDoublePrecisionSolver
                          >) {
                auto& vexcl_context = GetVexCLContext();
                KRATOS_ERROR_IF_NOT(vexcl_context) << "invalid VexCL context state\n";

                std::vector<float> x_float(rX.begin(), rX.end());
                std::vector<float> b_float(rB.begin(), rB.end());
                vex::vector<float> x(vexcl_context, x_float.size(), x_float.data());
                vex::vector<float> b(vexcl_context, b_float.size(), b_float.data());
                const auto results = rpSolver->operator()(b, x);

                vex::copy(GetVexCLContext(), x.begin(), x.end(), x_float.begin());
                std::copy(x_float.begin(), x_float.end(), rX.begin());
                return results;
            } else
            #endif
            {
                const auto results = rpSolver->operator()(rB, rX);
                return results;
            }
        },
        mpImpl->mSolverStruct.mSolver
    );

    KRATOS_WARNING_IF("AMGCLRawSolver", 1 <= mpImpl->mVerbosity && mpImpl->mTolerance <= residual)
        << "Failed to converge. Residual: " << residual << "\n";

    if(1 < mpImpl->mVerbosity) {
        std::cout << "Iterations: " << iteration_count << "\n"
                  << "Error: " << residual << "\n"
                  << "\n";
    }

    // Construct solver
    KRATOS_INFO_IF("AMGCLRawSolver", 2 <= mpImpl->mVerbosity)
        << "Solver memory usage: "
        << amgcl::human_readable_memory(std::visit([](auto& rpSolver){return amgcl::backend::bytes(*rpSolver);},
                                                   mpImpl->mSolverStruct.mSolver))
        << "\n";

    return residual < mpImpl->mTolerance ? true : false;
    KRATOS_CATCH("")
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void AMGCLRawSolver<TSparseSpace,TDenseSpace,TReorderer>::ProvideAdditionalData(SparseMatrix& rA,
                                                                                Vector& rX,
                                                                                Vector& rB,
                                                                                ModelPart::DofsArrayType& rDofs,
                                                                                ModelPart& rModelPart)
{
    KRATOS_TRY
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);

    std::optional<std::size_t> old_dof_count;
    std::size_t dof_count = 0;

    // Compute block size
    if (!rModelPart.IsDistributed()) {
        std::size_t old_dof_id = rDofs.empty() ? 0ul : rDofs.begin()->Id();
        for (const auto& rDof : rDofs) {
            if (rDof.EquationId() < TSparseSpace::Size1(rA)) {
                const auto dof_id = rDof.Id();
                if(dof_id != old_dof_id) {
                    old_dof_id = dof_id;
                    if (!old_dof_count.has_value()) {
                        old_dof_count = dof_count;
                    } else if (old_dof_count != dof_count) { //if it is different than the block size is 1
                        old_dof_count = -1;
                        break;
                    }
                    dof_count = 1;
                } else { // else (dof_id != old_dof_id)
                    ++dof_count;
                }
            }
        }
        mpImpl->mDoFCount = old_dof_count.has_value() ? 1 : dof_count;
    } else { //distributed
        const std::size_t system_size = TSparseSpace::Size1(rA);
        int current_rank = rModelPart.GetCommunicator().GetDataCommunicator().Rank();
        unsigned int old_node_id = rDofs.size() ? rDofs.begin()->Id() : 0;
        for (auto it = rDofs.begin(); it!=rDofs.end(); it++) {
            if(it->EquationId() < system_size  && it->GetSolutionStepValue(PARTITION_INDEX) == current_rank) {
                IndexType id = it->Id();
                if(id != old_node_id) {
                    old_node_id = id;
                    if(!old_dof_count.has_value()) {
                        old_dof_count = dof_count;
                    } else if (old_dof_count.value() != dof_count) { //if it is different than the block size is 1
                        old_dof_count = -1;
                        break;
                    }

                    dof_count=1;
                } else {
                    dof_count++;
                }
            }
        }

        if(old_dof_count.has_value()) {
            mpImpl->mDoFCount = dof_count;
        }

        const std::size_t max_block_size = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(mpImpl->mDoFCount);

        if(!old_dof_count.has_value()) {
            mpImpl->mDoFCount = max_block_size;
        }

        KRATOS_ERROR_IF(mpImpl->mDoFCount != max_block_size) << "Block size is not consistent. Local: " << mpImpl->mDoFCount  << " Max: " << max_block_size << std::endl;
    }

    KRATOS_INFO_IF("AMGCLRawSolver", 1 <= mpImpl->mVerbosity)
        << "block size: " << mpImpl->mDoFCount << "\n";

    // Construct solver and matrix adapter
    mpImpl->mpA = &rA;
    std::visit([&rA, this](auto& rAMGCLBundle){
            using Traits = typename std::remove_reference_t<decltype(rAMGCLBundle)>::Traits;
            using ScalarType = typename Traits::Scalar;
            constexpr bool is_gpu_bound = Traits::IsGPUBound;

            if constexpr (std::is_same_v<ScalarType,double>) {
                if constexpr (is_gpu_bound) {
                    rAMGCLBundle.mpMatrixAdapter = std::shared_ptr<typename Traits::MatrixAdapter>(
                        new typename Traits::MatrixAdapter(rA));
                } else {
                    rAMGCLBundle.mpMatrixAdapter = amgcl::adapter::zero_copy(TSparseSpace::Size1(rA),
                                                                             rA.index1_data().begin(),
                                                                             rA.index2_data().begin(),
                                                                             rA.value_data().begin());
                }
            } else if constexpr (std::is_same_v<ScalarType,float>) {
                if constexpr (is_gpu_bound) {
                    rAMGCLBundle.mpMatrixAdapter = std::shared_ptr<typename Traits::MatrixAdapter>(
                        new typename Traits::MatrixAdapter(rA));
                } else {
                    rAMGCLBundle.mpMatrixAdapter = amgcl::adapter::zero_copy(TSparseSpace::Size1(rA),
                                                                             rA.index1_data().begin(),
                                                                             rA.index2_data().begin(),
                                                                             rA.value_data().begin());
                }
            } else {
                static_assert(std::is_same_v<ScalarType,void>, "unhandled variant");
            }
        },
        mpImpl->mSolverVariant);

    // Construct the solver
    switch (mpImpl->mBackendType) {
        #ifdef AMGCL_GPGPU
            case AMGCLBackendType::SinglePrecisionGPU: {
                VexCLDoublePrecisionBackend::params vexcl_parameters;
                vexcl_parameters.q = GetVexCLContext();

                RegisterVexCLStaticMatrixType<2>();
                RegisterVexCLStaticMatrixType<3>();
                RegisterVexCLStaticMatrixType<4>();

                mpImpl->mSolverStruct.mSolver = std::make_unique<VexCLSinglePrecisionSolver>(
                    *mpImpl->mSolverStruct.mpAdapter,
                    mpImpl->mAMGCLSettings
                );
                break;
            } // case AMGCLBackendType::SinglePrecisionGPU
            case AMGCLBackendType::DoublePrecisionGPU: {
                VexCLDoublePrecisionBackend::params vexcl_parameters;
                vexcl_parameters.q = GetVexCLContext();

                RegisterVexCLStaticMatrixType<2>();
                RegisterVexCLStaticMatrixType<3>();
                RegisterVexCLStaticMatrixType<4>();

                mpImpl->mSolverStruct.mSolver = std::make_unique<VexCLDoublePrecisionSolver>(
                    *mpImpl->mSolverStruct.mpAdapter,
                    mpImpl->mAMGCLSettings
                );
                break;
            } // case AMGCLBackendType::DoublePredicisonGPU
        #endif
        case AMGCLBackendType::CPU: {
            switch (mpImpl->mDoFCount) {
                case 1:
                    mpImpl->mSolverStruct.mSolver = std::make_unique<ScalarSolver>(
                        *mpImpl->mSolverStruct.mpAdapter,
                        mpImpl->mAMGCLSettings
                    );
                    break;
                case 2:
                    mpImpl->mSolverStruct.mSolver = std::make_unique<BlockSolver<2>>(
                        *mpImpl->mSolverStruct.mpAdapter,
                        mpImpl->mAMGCLSettings
                    );
                    break;
                case 3:
                    mpImpl->mSolverStruct.mSolver = std::make_unique<BlockSolver<3>>(
                        *mpImpl->mSolverStruct.mpAdapter,
                        mpImpl->mAMGCLSettings
                    );
                    break;
                case 4:
                    mpImpl->mSolverStruct.mSolver = std::make_unique<BlockSolver<4>>(
                        *mpImpl->mSolverStruct.mpAdapter,
                        mpImpl->mAMGCLSettings
                    );
                    break;
                default:
                    KRATOS_ERROR << "unsupported block size: " << mpImpl->mDoFCount << "\n";
            } // switch mpImpl->mDoFCount
            break;
        } // case AMGCLBackendType::CPU
        default: KRATOS_ERROR << "unhandled AMGCL backend enum: " << (int)mpImpl->mBackendType << "\n";
    } // switch mpImpl->mBackendType

    KRATOS_INFO_IF("AMGCLRawSolver", 1 < mpImpl->mVerbosity)
        << "Block DoFs: " << mpImpl->mDoFCount << "\n";
    KRATOS_CATCH("")
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
Parameters
AMGCLRawSolver<TSparseSpace,TDenseSpace,TReorderer>::GetDefaultParameters()
{
    return Parameters(R"(
{
    "solver_type" : "amgcl_raw",
    "verbosity" : 0,
    "tolerance" : 1e-6,
    "gpgpu_backend" : "",
    "amgcl_settings" : {
        "precond" : {
            "class" : "amg",
            "relax" : {
                "type" : "ilu0"
            },
            "coarsening" : {
                "type" : "aggregation",
                "aggr" : {
                    "eps_strong" : 0.08,
                    "block_size" : 1
                }
            },
            "coarse_enough" : 333,
            "npre" : 1,
            "npost" : 1
        },
        "solver" : {
            "type" : "cg",
            "maxiter" : 555,
            "tol" : 1e-6
        }
    }
}
    )");
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool AMGCLRawSolver<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                DenseMatrix& rX,
                                                                DenseMatrix& rB)
{
    return false;
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void AMGCLRawSolver<TSparseSpace,TDenseSpace,TReorderer>::PrintInfo(std::ostream& rStream) const
{
    rStream << "AMGCLRawSolver";
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void AMGCLRawSolver<TSparseSpace,TDenseSpace,TReorderer>::PrintData(std::ostream& rStream) const
{
    rStream
        << "tolerance     : " << mpImpl->mTolerance << "\n"
        << "DoF size      : " << mpImpl->mDoFCount << "\n"
        << "verbosity     : " << mpImpl->mVerbosity << "\n"
        << "AMGCL settings: "
        ;
    boost::property_tree::json_parser::write_json(rStream, mpImpl->mAMGCLSettings);
}



template
class AMGCLRawSolver<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>,
    Reorderer<
        TUblasSparseSpace<double>,
        TUblasDenseSpace<double>
    >
>;


} // namespace Kratos
