//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// Trilinos includes
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_MpiComm.h"

// KratosCore dependencies
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "spaces/ublas_space.h"

// TrilinosApplication dependencies
#include "trilinos_space.h"

// strategies
#include "custom_strategies/coupled_strategy.h"
#include "custom_strategies/coupled_strategy_item.h"

// schemes
#include "custom_strategies/generic_convergence_criteria.h"
#include "custom_strategies/generic_residual_based_bossak_velocity_scalar_scheme.h"
#include "custom_strategies/generic_residualbased_simple_steady_scalar_scheme.h"

// Include base h
#include "add_trilinos_strategies_to_python.h"

namespace Kratos
{
namespace Python
{
void AddTrilinosStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using MPISparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using MPILinearSolverType = LinearSolver<MPISparseSpaceType, LocalSpaceType>;
    using MPIBaseSchemeType = Scheme<MPISparseSpaceType, LocalSpaceType>;
    using MPIBaseConvergenceCriteriaType = ConvergenceCriteria<MPISparseSpaceType, LocalSpaceType>;
    using MPIBaseSolvingStrategyType = SolvingStrategy<MPISparseSpaceType, LocalSpaceType, MPILinearSolverType>;

    // Add coupled strategy item
    using MPICoupledStrategyItemType = CoupledStrategyItem<MPISparseSpaceType, LocalSpaceType, MPILinearSolverType>;
    py::class_<MPICoupledStrategyItemType, typename MPICoupledStrategyItemType::Pointer>(m, "MPICoupledStrategyItem")
        .def(py::init<MPIBaseSolvingStrategyType::Pointer, std::string, int>())
        .def(py::init<MPIBaseSolvingStrategyType::Pointer, std::string, int, std::vector<int>>())
        .def("AddAuxiliaryProcess", &MPICoupledStrategyItemType::AddAuxiliaryProcess)
        .def("GetName", &MPICoupledStrategyItemType::GetName)
        .def("GetStrategy", &MPICoupledStrategyItemType::GetStrategy)
        .def("GetAuxiliaryProcessList", &MPICoupledStrategyItemType::GetStrategy)
        .def("GetStrategyInfo", &MPICoupledStrategyItemType::GetStrategyInfo)
        .def("GetStrategySolvabilityPattern", &MPICoupledStrategyItemType::GetStrategySolvabilityPattern)
        .def("SetStrategySolvabilityPattern", &MPICoupledStrategyItemType::SetStrategySolvabilityPattern);;

    // Add strtegies
    using MPICoupledStrategyType = CoupledStrategy<MPISparseSpaceType, LocalSpaceType, MPILinearSolverType>;
    py::class_<MPICoupledStrategyType, typename MPICoupledStrategyType::Pointer, MPIBaseSolvingStrategyType>(m, "MPICoupledStrategy")
        .def(py::init<ModelPart&, bool, bool, bool, int>())
        .def("AddStrategyItem", &MPICoupledStrategyType::AddStrategyItem)
        .def("AddConvergenceCheckVariable", &MPICoupledStrategyType::AddConvergenceCheckVariable);

    // Convergence criteria
    using MPIGenericConvergenceCriteriaType = GenericConvergenceCriteria<MPISparseSpaceType, LocalSpaceType>;
    py::class_<MPIGenericConvergenceCriteriaType, typename MPIGenericConvergenceCriteriaType::Pointer, MPIBaseConvergenceCriteriaType>(m, "MPIGenericScalarConvergenceCriteria")
        .def(py::init<double, double>());

    // add schemes
    using MPIGenericResidualBasedBossakVelocityScalarSchemeType = GenericResidualBasedBossakVelocityScalarScheme<MPISparseSpaceType, LocalSpaceType>;
    py::class_<MPIGenericResidualBasedBossakVelocityScalarSchemeType, typename MPIGenericResidualBasedBossakVelocityScalarSchemeType::Pointer, MPIBaseSchemeType>(m, "MPIGenericResidualBasedBossakVelocityDynamicScalarScheme")
        .def(py::init<const double, const double, const Variable<double>&, const Variable<double>&, const Variable<double>&>());

    using MPIGenericResidualBasedSimpleSteadyScalarSchemeType = GenericResidualBasedSimpleSteadyScalarScheme<MPISparseSpaceType, LocalSpaceType>;
    py::class_<MPIGenericResidualBasedSimpleSteadyScalarSchemeType, typename MPIGenericResidualBasedSimpleSteadyScalarSchemeType::Pointer, MPIBaseSchemeType>(m, "MPIGenericResidualBasedSimpleSteadyScalarScheme")
        .def(py::init<const double>());
}

} // namespace Python
} // namespace Kratos
