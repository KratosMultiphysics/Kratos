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

#include "add_trilinos_utilities_to_python.h"

// Trilinos includes
#include "Epetra_FEVector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_MpiComm.h"

// KratosCore dependencies
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "spaces/ublas_space.h"

// TrilinosApplication dependencies
#include "trilinos_space.h"

// FluidDynamics trilinos extensions
#include "custom_strategies/strategies/fractional_step_strategy.h"
#include "custom_utilities/solver_settings.h"

// adjoint schemes
#include "custom_strategies/schemes/simple_steady_adjoint_scheme.h"
#include "custom_strategies/schemes/velocity_bossak_adjoint_scheme.h"

// sensitivity builder schemes
#include "custom_strategies/schemes/simple_steady_sensitivity_builder_scheme.h"
#include "custom_strategies/schemes/velocity_bossak_sensitivity_builder_scheme.h"

namespace Kratos {
namespace Python {

void AddTrilinosStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using TrilinosSparseSpace = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using UblasLocalSpace = UblasSpace<double, Matrix, Vector>;
    using TrilinosLinearSolver = LinearSolver<TrilinosSparseSpace, UblasLocalSpace>;

    using TrilinosBaseSolvingStrategy = ImplicitSolvingStrategy< TrilinosSparseSpace, UblasLocalSpace, TrilinosLinearSolver >;
    using BaseSolverSettings = SolverSettings<TrilinosSparseSpace, UblasLocalSpace, TrilinosLinearSolver>;
    using BaseSchemeType = Scheme<TrilinosSparseSpace, UblasLocalSpace>;

    using TrilinosFractionalStepStrategy = FractionalStepStrategy< TrilinosSparseSpace, UblasLocalSpace, TrilinosLinearSolver>;
    py::class_< TrilinosFractionalStepStrategy, typename TrilinosFractionalStepStrategy::Pointer, TrilinosBaseSolvingStrategy >(m,"TrilinosFractionalStepStrategy")
    .def(py::init< ModelPart&, BaseSolverSettings&, bool >())
    .def(py::init< ModelPart&, BaseSolverSettings&, bool, bool >())
    .def(py::init< ModelPart&, BaseSolverSettings&, bool, const Kratos::Variable<int>& >())
    .def(py::init< ModelPart&, BaseSolverSettings&, bool, bool, const Kratos::Variable<int>& >())
    .def("CalculateReactions", [](TrilinosFractionalStepStrategy& self) {
        KRATOS_WARNING("TrilinosFractionalStepStrategy") << "\'CalculateReactions()\' exposure is deprecated. Use the constructor with the \'CalculateReactionsFlag\' instead." << std::endl;
        self.CalculateReactions();})
    .def("AddIterationStep",&TrilinosFractionalStepStrategy::AddIterationStep)
    .def("ClearExtraIterationSteps",&TrilinosFractionalStepStrategy::ClearExtraIterationSteps)
    ;

    using TrilinosSimpleSteadyAdjointSchemeType = SimpleSteadyAdjointScheme<TrilinosSparseSpace, UblasLocalSpace>;
    py::class_<TrilinosSimpleSteadyAdjointSchemeType, typename TrilinosSimpleSteadyAdjointSchemeType::Pointer, BaseSchemeType>
        (m, "TrilinosSimpleSteadyAdjointScheme")
        .def(py::init<AdjointResponseFunction::Pointer, const std::size_t, const std::size_t>())
        .def(py::init<AdjointResponseFunction::Pointer, const std::size_t, const std::size_t, ElementRefinementProcess::Pointer>())
        ;

    using TrilinosVelocityBossakAdjointSchemeType = VelocityBossakAdjointScheme<TrilinosSparseSpace, UblasLocalSpace>;
    py::class_<TrilinosVelocityBossakAdjointSchemeType, typename TrilinosVelocityBossakAdjointSchemeType::Pointer, BaseSchemeType>
        (m, "TrilinosVelocityBossakAdjointScheme")
        .def(py::init<Parameters, AdjointResponseFunction::Pointer, const std::size_t, const std::size_t>())
        .def(py::init<Parameters, AdjointResponseFunction::Pointer, const std::size_t, const std::size_t, ElementRefinementProcess::Pointer>())
        .def(py::init<Parameters, AdjointResponseFunction::Pointer, const std::size_t, const std::size_t, ElementRefinementProcess::Pointer, const std::size_t>())
        ;

}

}
}
