//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

// System includes

// External includes

// Project includes
#include "factories/standard_linear_solver_factory.h"
#include "linear_solvers/cg_solver.h"
#include "linear_solvers/deflated_cg_solver.h"
#include "linear_solvers/bicgstab_solver.h"
#include "linear_solvers/tfqmr_solver.h"
#include "linear_solvers/amgcl_solver.h"
#include "linear_solvers/amgcl_ns_solver.h"
#include "linear_solvers/scaling_solver.h"
#include "linear_solvers/fallback_linear_solver.h"
#include "linear_solvers/monotonicity_preserving_solver.h"
#include "linear_solvers/skyline_lu_custom_scalar_solver.h"
#include "spaces/ublas_space.h"

namespace Kratos {


namespace detail {


template <class TSparseDataType,
          class TDenseDataType>
void RegisterLinearSolvers()
{
    using SpaceType = TUblasSparseSpace<TSparseDataType>;
    using LocalSpaceType = TUblasDenseSpace<TDenseDataType>;
    using ComplexSpaceType = TUblasSparseSpace<std::complex<TSparseDataType>>;
    using ComplexLocalSpaceType = TUblasDenseSpace<std::complex<TDenseDataType>>;

    using CGSolverType = CGSolver<SpaceType,  LocalSpaceType>;
    static auto CGSolverFactory = StandardLinearSolverFactory<SpaceType,LocalSpaceType,CGSolverType>();
    KratosComponents<LinearSolverFactory<SpaceType,LocalSpaceType>>::Add("cg", CGSolverFactory);

    using BICGSTABSolverType = BICGSTABSolver<SpaceType,  LocalSpaceType>;
    static auto BICGSTABSolverFactory = StandardLinearSolverFactory<SpaceType,LocalSpaceType,BICGSTABSolverType>();
    KratosComponents<LinearSolverFactory<SpaceType,LocalSpaceType>>::Add("bicgstab", BICGSTABSolverFactory);

    using SkylineLUFactorizationSolverType = SkylineLUFactorizationSolver<SpaceType,  LocalSpaceType>;
    static auto SkylineLUFactorizationSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,SkylineLUFactorizationSolverType>();
    KratosComponents<LinearSolverFactory<SpaceType,LocalSpaceType>>::Add("skyline_lu_factorization",SkylineLUFactorizationSolverFactory );

    using TFQMRSolverType = TFQMRSolver<SpaceType,  LocalSpaceType>;
    static auto TFQMRSolverFactory = StandardLinearSolverFactory<SpaceType,LocalSpaceType,TFQMRSolverType>();
    KratosComponents<LinearSolverFactory<SpaceType,LocalSpaceType>>::Add("tfqmr", TFQMRSolverFactory);

    using AMGCLSolverType = AMGCLSolver<SpaceType,  LocalSpaceType>;
    static auto AMGCLSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,AMGCLSolverType>();
    KratosComponents<LinearSolverFactory<SpaceType,LocalSpaceType>>::Add("amgcl", AMGCLSolverFactory);

    using AMGCL_NS_SolverType = AMGCL_NS_Solver<SpaceType,  LocalSpaceType>;
    static auto AMGCL_NS_SolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,AMGCL_NS_SolverType>();
    KratosComponents<LinearSolverFactory<SpaceType,LocalSpaceType>>::Add("amgcl_ns",AMGCL_NS_SolverFactory );

    using ScalingSolverType = ScalingSolver<SpaceType, LocalSpaceType>;
    static auto ScalingSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,ScalingSolverType>();
    KratosComponents<LinearSolverFactory<SpaceType,LocalSpaceType>>::Add("scaling",ScalingSolverFactory );

    using FallbackLinearSolverType = FallbackLinearSolver<SpaceType, LocalSpaceType>;
    static auto FallbackLinearSolverFactory = StandardLinearSolverFactory<SpaceType, LocalSpaceType, FallbackLinearSolverType>();
    KratosComponents<LinearSolverFactory<SpaceType,LocalSpaceType>>::Add("fallback_linear_solver", FallbackLinearSolverFactory);

    using MonotonicityPreservingSolverType = MonotonicityPreservingSolver<SpaceType,  LocalSpaceType>;
    static auto MonotonicityPreservingSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,MonotonicityPreservingSolverType>();
    KratosComponents<LinearSolverFactory<SpaceType,LocalSpaceType>>::Add("monotonicity_preserving",MonotonicityPreservingSolverFactory );

    using SkylineLUComplexSolverType = SkylineLUCustomScalarSolver<ComplexSpaceType, ComplexLocalSpaceType>;
    static auto SkylineLUComplexSolverFactory = StandardLinearSolverFactory<ComplexSpaceType, ComplexLocalSpaceType, SkylineLUComplexSolverType>();
    KratosComponents<LinearSolverFactory<ComplexSpaceType,ComplexLocalSpaceType>>::Add("skyline_lu_complex", SkylineLUComplexSolverFactory);

    // using LinearSolverType = LinearSolver<SpaceType,  LocalSpaceType>;
    // using IterativeSolverType = IterativeSolver<SpaceType,  LocalSpaceType>;
//         KRATOS_REGISTER_LINEAR_SOLVER("LinearSolver", StandardLinearSolverFactory<SpaceType,LocalSpaceType,LinearSolverType>());

    // DeflatedCGSolver's utilities are not templated on the space type (and thus the floating point type),
    // so they can only be defined for double precision floats.
    if constexpr (std::is_same_v<TSparseDataType,double>) {
        using DeflatedCGSolverType = DeflatedCGSolver<SpaceType,  LocalSpaceType>;
        static auto DeflatedCGSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,DeflatedCGSolverType>();
        KratosComponents<LinearSolverFactory<SpaceType,LocalSpaceType>>::Add("deflated_cg", DeflatedCGSolverFactory);
    }
}


} // namespace detail


void RegisterLinearSolvers()
{
    detail::RegisterLinearSolvers</*TSparseDataType=*/double,/*TDenseDataType*/double>();
    detail::RegisterLinearSolvers</*TSparseDataType=*/float,/*TDenseDataType*/double>();
};
} // Namespace Kratos
