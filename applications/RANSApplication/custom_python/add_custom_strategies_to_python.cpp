//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
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
#include "custom_strategies/bossak_relaxation_scalar_scheme.h"
#include "custom_strategies/steady_scalar_scheme.h"
#include "custom_strategies/algebraic_flux_corrected_steady_scalar_scheme.h"

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

    // strategies
    using RansFractionalStepStrategyType = RansFractionalStepStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    py::class_<RansFractionalStepStrategyType, typename RansFractionalStepStrategyType::Pointer, BaseSolvingStrategyType>(m, "RansFractionalStepStrategy")
        .def(py::init<ModelPart&, SolverSettings<SparseSpaceType, LocalSpaceType, LinearSolverType>&, bool, bool>())
        .def(py::init<ModelPart&, SolverSettings<SparseSpaceType, LocalSpaceType, LinearSolverType>&, bool, bool, const Kratos::Variable<int>&>());

    // add schemes
    using SteadyScalarSchemeType = SteadyScalarScheme<SparseSpaceType, LocalSpaceType>;
    py::class_<SteadyScalarSchemeType, typename SteadyScalarSchemeType::Pointer, BaseSchemeType>(m, "SteadyScalarScheme")
        .def(py::init<const double>());

    using AlgebraicFluxCorrectedSteadyScalarSchemeType = AlgebraicFluxCorrectedSteadyScalarScheme<SparseSpaceType, LocalSpaceType>;
    py::class_<AlgebraicFluxCorrectedSteadyScalarSchemeType, typename AlgebraicFluxCorrectedSteadyScalarSchemeType::Pointer, BaseSchemeType>(m, "AlgebraicFluxCorrectedSteadyScalarScheme")
        .def(py::init<const double, const Flags&>())
        .def(py::init<const double, const Flags&, const Variable<int>&>());

    using BossakRelaxationScalarSchemeType = BossakRelaxationScalarScheme<SparseSpaceType, LocalSpaceType>;
    py::class_<BossakRelaxationScalarSchemeType, typename BossakRelaxationScalarSchemeType::Pointer, BaseSchemeType>(m, "BossakRelaxationScalarScheme")
        .def(py::init<const double, const double, const Variable<double>&, const Variable<double>&, const Variable<double>&>());

    // Convergence criteria
    using GenericConvergenceCriteriaType = GenericConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    py::class_<GenericConvergenceCriteriaType, typename GenericConvergenceCriteriaType::Pointer, BaseConvergenceCriteriaType>(m, "GenericScalarConvergenceCriteria")
        .def(py::init<double, double>());


}

} // namespace Python.
} // Namespace Kratos
