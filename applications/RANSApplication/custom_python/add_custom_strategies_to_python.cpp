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
#include "custom_python/add_custom_strategies_to_python.h"
#include "includes/define_python.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

// strategies
#include "custom_strategies/coupled_strategy.h"
#include "custom_strategies/coupled_strategy_item.h"

// schemes
#include "custom_strategies/generic_residual_based_bossak_velocity_scalar_scheme.h"
#include "custom_strategies/generic_residualbased_simple_steady_scalar_scheme.h"

// convergence criterians
#include "custom_strategies/generic_convergence_criteria.h"

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
    using BaseSolvingStrategyType = SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    using BaseConvergenceCriteriaType = ConvergenceCriteria<SparseSpaceType, LocalSpaceType>;

    // Add coupled strategy item
    using CoupledStrategyItemType = CoupledStrategyItem<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    py::class_<CoupledStrategyItemType, typename CoupledStrategyItemType::Pointer>(m, "CoupledStrategyItem")
        .def(py::init<BaseSolvingStrategyType::Pointer, std::string, int>())
        .def(py::init<BaseSolvingStrategyType::Pointer, std::string, int, std::vector<int>>())
        .def("AddAuxiliaryProcess", &CoupledStrategyItemType::AddAuxiliaryProcess)
        .def("GetName", &CoupledStrategyItemType::GetName)
        .def("GetStrategy", &CoupledStrategyItemType::GetStrategy)
        .def("GetAuxiliaryProcessList", &CoupledStrategyItemType::GetStrategy)
        .def("GetStrategyInfo", &CoupledStrategyItemType::GetStrategyInfo)
        .def("GetStrategySolvabilityPattern", &CoupledStrategyItemType::GetStrategySolvabilityPattern)
        .def("SetStrategySolvabilityPattern", &CoupledStrategyItemType::SetStrategySolvabilityPattern);

    // Add strtegies
    using CoupledStrategyType = CoupledStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    py::class_<CoupledStrategyType, typename CoupledStrategyType::Pointer, BaseSolvingStrategyType>(m, "CoupledStrategy")
        .def(py::init<ModelPart&, bool, bool, bool, int>())
        .def("AddStrategyItem", &CoupledStrategyType::AddStrategyItem)
        .def("AddConvergenceCheckVariable", &CoupledStrategyType::AddConvergenceCheckVariable);

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

}

} // namespace Python.
} // Namespace Kratos
