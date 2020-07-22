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

// RANS trilinos extensions
// strategies
#include "custom_strategies/algebraic_flux_corrected_steady_scalar_scheme.h"
#include "custom_strategies/generic_convergence_criteria.h"
#include "custom_strategies/bossak_scalar_scheme.h"
#include "custom_strategies/steady_scalar_scheme.h"
#include "custom_strategies/rans_fractional_step_strategy.h"
#include "custom_utilities/solver_settings.h"
#include "solving_strategies/strategies/solving_strategy.h"

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

    // strategies
    using MPIRansFractionalStepStrategyType = RansFractionalStepStrategy<MPISparseSpaceType, LocalSpaceType, MPILinearSolverType>;
    py::class_<MPIRansFractionalStepStrategyType, typename MPIRansFractionalStepStrategyType::Pointer, MPIBaseSolvingStrategyType>(m, "MPIRansFractionalStepStrategy")
        .def(py::init<ModelPart&, SolverSettings<MPISparseSpaceType, LocalSpaceType, MPILinearSolverType>&, bool, bool>())
        .def(py::init<ModelPart&, SolverSettings<MPISparseSpaceType, LocalSpaceType, MPILinearSolverType>&, bool, bool, const Kratos::Variable<int>&>());

    // add schemes
    using MPIAlgebraicFluxCorrectedSteadyScalarSchemeType = AlgebraicFluxCorrectedSteadyScalarScheme<MPISparseSpaceType, LocalSpaceType>;
    py::class_<MPIAlgebraicFluxCorrectedSteadyScalarSchemeType, typename MPIAlgebraicFluxCorrectedSteadyScalarSchemeType::Pointer, MPIBaseSchemeType>(m, "MPIAlgebraicFluxCorrectedSteadyScalarScheme")
        .def(py::init<const double, const Flags&>())
        .def(py::init<const double, const Flags&, const Variable<int>&>());

    using MPISteadyScalarSchemeType = SteadyScalarScheme<MPISparseSpaceType, LocalSpaceType>;
    py::class_<MPISteadyScalarSchemeType, typename MPISteadyScalarSchemeType::Pointer, MPIBaseSchemeType>(m, "MPISteadyScalarScheme")
        .def(py::init<const double>());

    using MPIBossakScalarSchemeType = BossakScalarScheme<MPISparseSpaceType, LocalSpaceType>;
    py::class_<MPIBossakScalarSchemeType, typename MPIBossakScalarSchemeType::Pointer, MPIBaseSchemeType>(m, "MPIBossakScalarScheme")
        .def(py::init<const double, const double, const Variable<double>&, const Variable<double>&, const Variable<double>&>());

    // Convergence criteria
    using MPIGenericConvergenceCriteriaType = GenericConvergenceCriteria<MPISparseSpaceType, LocalSpaceType>;
    py::class_<MPIGenericConvergenceCriteriaType, typename MPIGenericConvergenceCriteriaType::Pointer, MPIBaseConvergenceCriteriaType>(m, "MPIGenericScalarConvergenceCriteria")
        .def(py::init<double, double>());
}

} // namespace Python
} // namespace Kratos
