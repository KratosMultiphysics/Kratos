// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

// External includes
#include "spaces/ublas_space.h"

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "includes/kratos_parameters.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_strategy.hpp"
#include "custom_strategies/strategies/geo_mechanics_ramm_arc_length_strategy.hpp"
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_erosion_process_strategy.hpp"

//builders and solvers
#include "custom_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_mass_and_damping.h"

//schemes
#include "custom_strategies/schemes/backward_euler_T_scheme.hpp"
#include "custom_strategies/schemes/backward_euler_quasistatic_Pw_scheme.hpp"
#include "custom_strategies/schemes/backward_euler_quasistatic_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/generalized_newmark_T_scheme.hpp"
#include "custom_strategies/schemes/newmark_dynamic_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_quasistatic_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_quasistatic_damped_U_Pw_scheme.hpp"

//linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos::Python {

void AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType >;
    using BaseSolvingStrategyType = ImplicitSolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >;
    using BaseSchemeType = Scheme< SparseSpaceType, LocalSpaceType >;
    using BuilderAndSolverType = BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >;
    using ConvergenceCriteriaType = ConvergenceCriteria< SparseSpaceType, LocalSpaceType >;

    using NewmarkQuasistaticUPwSchemeType = NewmarkQuasistaticUPwScheme< SparseSpaceType, LocalSpaceType >;
    using NewmarkQuasistaticDampedUPwSchemeType = NewmarkQuasistaticDampedUPwScheme< SparseSpaceType, LocalSpaceType >;
    using NewmarkDynamicUPwSchemeType = NewmarkDynamicUPwScheme< SparseSpaceType, LocalSpaceType >;
    using NewmarkQuasistaticPwSchemeType = NewmarkQuasistaticPwScheme< SparseSpaceType, LocalSpaceType >;
    using NewmarkQuasistaticTSchemeType =
        GeneralizedNewmarkTScheme< SparseSpaceType, LocalSpaceType >;

    using BackwardEulerQuasistaticUPwSchemeType = BackwardEulerQuasistaticUPwScheme< SparseSpaceType, LocalSpaceType >;
    using BackwardEulerQuasistaticPwSchemeType = BackwardEulerQuasistaticPwScheme< SparseSpaceType, LocalSpaceType >;
    using BackwardEulerQuasistaticTSchemeType =
        BackwardEulerTScheme< SparseSpaceType, LocalSpaceType >;

    using GeoMechanicsNewtonRaphsonStrategyType = GeoMechanicsNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >;
    using GeoMechanicsRammArcLengthStrategyType = GeoMechanicsRammArcLengthStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >;
    using GeoMechanicsNewtonRaphsonErosionProcessStrategyType = GeoMechanicsNewtonRaphsonErosionProcessStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    py::class_< NewmarkQuasistaticUPwSchemeType, typename NewmarkQuasistaticUPwSchemeType::Pointer, BaseSchemeType >
    (m, "NewmarkQuasistaticUPwScheme", py::module_local())
    .def(py::init<  double, double, double >());

    py::class_< NewmarkQuasistaticDampedUPwSchemeType, typename NewmarkQuasistaticDampedUPwSchemeType::Pointer, BaseSchemeType >
    (m, "NewmarkQuasistaticDampedUPwScheme", py::module_local())
    .def(py::init<  double, double, double >());

    py::class_< NewmarkDynamicUPwSchemeType,typename NewmarkDynamicUPwSchemeType::Pointer, BaseSchemeType >
    (m, "NewmarkDynamicUPwScheme", py::module_local())
    .def(py::init<  double, double, double >());

    py::class_< NewmarkQuasistaticPwSchemeType, typename NewmarkQuasistaticPwSchemeType::Pointer, BaseSchemeType >
    (m, "NewmarkQuasistaticPwScheme")
    .def(py::init<  double >());

    py::class_<NewmarkQuasistaticTSchemeType, typename NewmarkQuasistaticTSchemeType::Pointer, BaseSchemeType>
        (m, "GeneralizedNewmarkTScheme")
            .def(py::init<double>());

    py::class_< BackwardEulerQuasistaticUPwSchemeType, typename BackwardEulerQuasistaticUPwSchemeType::Pointer, BaseSchemeType >
    (m, "BackwardEulerQuasistaticUPwScheme")
    .def(py::init< >());

    py::class_< BackwardEulerQuasistaticPwSchemeType, typename BackwardEulerQuasistaticPwSchemeType::Pointer, BaseSchemeType >
    (m, "BackwardEulerQuasistaticPwScheme")
    .def(py::init< >());

    py::class_< BackwardEulerQuasistaticTSchemeType, typename BackwardEulerQuasistaticTSchemeType::Pointer, BaseSchemeType >
        (m, "BackwardEulerTScheme")
            .def(py::init< >());

    py::class_< GeoMechanicsNewtonRaphsonStrategyType, typename GeoMechanicsNewtonRaphsonStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "GeoMechanicsNewtonRaphsonStrategy")
    .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
        BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >());

    py::class_< GeoMechanicsNewtonRaphsonErosionProcessStrategyType, typename GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "GeoMechanicsNewtonRaphsonErosionProcessStrategy")
    .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
        BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >());

    py::class_< GeoMechanicsRammArcLengthStrategyType, typename GeoMechanicsRammArcLengthStrategyType::Pointer, BaseSolvingStrategyType >
    (m, "GeoMechanicsRammArcLengthStrategy")
    .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
        BuilderAndSolverType::Pointer, Parameters&, int, bool, bool, bool >())
    .def("UpdateLoads",&GeoMechanicsRammArcLengthStrategyType::UpdateLoads);



    using ResidualBasedBlockBuilderAndSolverWithMassAndDampingType = ResidualBasedBlockBuilderAndSolverWithMassAndDamping< SparseSpaceType, LocalSpaceType, LinearSolverType >;
    py::class_< ResidualBasedBlockBuilderAndSolverWithMassAndDampingType, ResidualBasedBlockBuilderAndSolverWithMassAndDampingType::Pointer, BuilderAndSolverType>(m, "ResidualBasedBlockBuilderAndSolverWithMassAndDamping")
        .def(py::init< LinearSolverType::Pointer >())
        .def(py::init< LinearSolverType::Pointer, Parameters >())
        ;
}

} // Namespace Kratos::Python
