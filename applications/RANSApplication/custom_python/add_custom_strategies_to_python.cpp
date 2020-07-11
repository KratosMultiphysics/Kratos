//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "spaces/ublas_space.h"

// strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/rans_fractional_step_strategy.h"
#include "custom_utilities/solver_settings.h"

// schemes
#include "custom_strategies/generic_residual_based_bossak_velocity_scalar_scheme.h"
#include "custom_strategies/generic_residualbased_simple_steady_scalar_scheme.h"
#include "custom_strategies/algebraic_flux_corrected_scalar_steady_scheme.h"

// convergence criterians
#include "custom_strategies/generic_convergence_criteria.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

// Include base h
#include "custom_python/add_custom_strategies_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
    using BaseSchemeType = Scheme<SparseSpaceType, LocalSpaceType>;
    using BaseConvergenceCriteriaType = ConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    using BaseSolvingStrategyType = SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;

    // Convergence criteria
    using GenericConvergenceCriteriaType = GenericConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    py::class_<GenericConvergenceCriteriaType, typename GenericConvergenceCriteriaType::Pointer, BaseConvergenceCriteriaType>(m, "GenericScalarConvergenceCriteria")
        .def(py::init<double, double>());

    // add schemes
    using GenericResidualBasedBossakVelocityScalarSchemeType = GenericResidualBasedBossakVelocityScalarScheme<SparseSpaceType, LocalSpaceType>;
    py::class_<GenericResidualBasedBossakVelocityScalarSchemeType, typename GenericResidualBasedBossakVelocityScalarSchemeType::Pointer, BaseSchemeType>(m, "GenericResidualBasedBossakVelocityDynamicScalarScheme")
        .def(py::init<const double, const double, const Variable<double>&, const Variable<double>&, const Variable<double>&>());

    using GenericResidualBasedSimpleSteadyScalarSchemeType = GenericResidualBasedSimpleSteadyScalarScheme<SparseSpaceType, LocalSpaceType>;
    py::class_<GenericResidualBasedSimpleSteadyScalarSchemeType, typename GenericResidualBasedSimpleSteadyScalarSchemeType::Pointer, BaseSchemeType>(m, "GenericResidualBasedSimpleSteadyScalarScheme")
        .def(py::init<const double>());

    using AlgebraicFluxCorrectedScalarSteadySchemeType = AlgebraicFluxCorrectedScalarSteadyScheme<SparseSpaceType, LocalSpaceType>;
    py::class_<AlgebraicFluxCorrectedScalarSteadySchemeType, typename AlgebraicFluxCorrectedScalarSteadySchemeType::Pointer, BaseSchemeType>(
        m, "AlgebraicFluxCorrectedScalarSteadyScheme")
        .def(py::init<const double, const Flags&>())
        .def(py::init<const double, const Flags&, const Variable<int>&>());

    // strategies
    using RansFractionalStepStrategyType = RansFractionalStepStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    py::class_<RansFractionalStepStrategyType, typename RansFractionalStepStrategyType::Pointer, BaseSolvingStrategyType>(m, "RansFractionalStepStrategy")
        .def(py::init<ModelPart&, SolverSettings<SparseSpaceType, LocalSpaceType, LinearSolverType>&, bool, bool>())
        .def(py::init<ModelPart&, SolverSettings<SparseSpaceType, LocalSpaceType, LinearSolverType>&, bool, bool, const Kratos::Variable<int>&>());

}

} // namespace Python.
} // Namespace Kratos
