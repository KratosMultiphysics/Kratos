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

// System includes

// External includes
#ifdef KRATOS_USE_AMATRIX
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#endif                                    // KRATOS_USE_AMATRIX

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

#include "processes/process.h"
#include "custom_utilities/solver_settings.h"

#include "spaces/ublas_space.h"

// builder_and_solvers
#include "solving_strategies/builder_and_solvers/explicit_builder.h"
#include "custom_strategies/builder_and_solvers/residualbased_block_builder_and_solver_periodic.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/fractional_step_strategy.h"

//schemes
#include "custom_strategies/schemes/bdf2_turbulent_scheme.h"
#include "custom_strategies/schemes/residualbased_simple_steady_scheme.h"
#include "custom_strategies/schemes/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent.h"
#include "custom_strategies/strategies/compressible_navier_stokes_explicit_solving_strategy_runge_kutta_4.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
namespace Python
{

void AddCustomStrategiesToPython(pybind11::module &m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> BaseSolvingStrategyType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    //********************************************************************
    //********************************************************************

    py::class_<
        ResidualBasedBlockBuilderAndSolverPeriodic<SparseSpaceType, LocalSpaceType, LinearSolverType>,
        typename ResidualBasedBlockBuilderAndSolverPeriodic<SparseSpaceType, LocalSpaceType, LinearSolverType>::Pointer,
        ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>>(m, "ResidualBasedBlockBuilderAndSolverPeriodic")
    .def(py::init<LinearSolverType::Pointer, const Variable<int> &>());

    py::class_<
        CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4<SparseSpaceType, LocalSpaceType>,
        typename CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4<SparseSpaceType, LocalSpaceType>::Pointer,
        ExplicitSolvingStrategyRungeKutta4<SparseSpaceType, LocalSpaceType>>(m, "CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4")
    .def(py::init<ModelPart&, bool, int>())
    .def(py::init<ModelPart&, Parameters>())
    .def(py::init<ModelPart&, ExplicitBuilder<SparseSpaceType, LocalSpaceType>::Pointer, bool, int>());

    py::class_<
        FractionalStepStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>,
        typename FractionalStepStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>::Pointer,
        BaseSolvingStrategyType>(m, "FractionalStepStrategy")
    .def(py::init<ModelPart &, SolverSettings<SparseSpaceType, LocalSpaceType, LinearSolverType> &, bool>())
    .def(py::init<ModelPart &, SolverSettings<SparseSpaceType, LocalSpaceType, LinearSolverType> &, bool, bool>())
    .def(py::init<ModelPart &, SolverSettings<SparseSpaceType, LocalSpaceType, LinearSolverType> &, bool, const Kratos::Variable<int> &>())
    .def(py::init<ModelPart &, SolverSettings<SparseSpaceType, LocalSpaceType, LinearSolverType> &, bool, bool, const Kratos::Variable<int> &>())
    .def("CalculateReactions", [](FractionalStepStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>& self) {
        KRATOS_WARNING("FractionalStepStrategy") << "\'CalculateReactions()\' exposure is deprecated. Use the constructor with the \'CalculateReactionsFlag\' instead." << std::endl;
        self.CalculateReactions();})
    .def("AddIterationStep", &FractionalStepStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>::AddIterationStep)
    .def("ClearExtraIterationSteps", &FractionalStepStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>::ClearExtraIterationSteps);

    py::class_<
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<SparseSpaceType, LocalSpaceType>,
        typename ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<SparseSpaceType, LocalSpaceType>::Pointer,
        BaseSchemeType>(m, "ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent")
    .def(py::init<double, double, unsigned int, Process::Pointer>())
    .def(py::init<double, double, unsigned int, double, Process::Pointer>())
    .def(py::init<double, double, unsigned int>())                        // constructor without a turbulence model
    .def(py::init<double, unsigned int, const Kratos::Variable<int> &>()) // constructor without a turbulence model for periodic boundary conditions
    ;

    typedef ResidualBasedSimpleSteadyScheme<SparseSpaceType, LocalSpaceType> ResidualBasedSimpleSteadySchemeType;
    py::class_<
        ResidualBasedSimpleSteadySchemeType,
        typename ResidualBasedSimpleSteadySchemeType::Pointer,
        BaseSchemeType>(m, "ResidualBasedSimpleSteadyScheme")
    .def(py::init<double, double, unsigned int, Process::Pointer>())
    .def(py::init<double, double, unsigned int>()) // constructor without a turbulence model
    .def("GetVelocityRelaxationFactor", &ResidualBasedSimpleSteadySchemeType::GetVelocityRelaxationFactor)
    .def("SetVelocityRelaxationFactor", &ResidualBasedSimpleSteadySchemeType::SetVelocityRelaxationFactor)
    .def("GetPressureRelaxationFactor", &ResidualBasedSimpleSteadySchemeType::GetPressureRelaxationFactor)
    .def("SetPressureRelaxationFactor", &ResidualBasedSimpleSteadySchemeType::SetPressureRelaxationFactor)
    ;

    py::class_<
        BDF2TurbulentScheme<SparseSpaceType, LocalSpaceType>,
        typename BDF2TurbulentScheme<SparseSpaceType, LocalSpaceType>::Pointer,
        BaseSchemeType>(m, "BDF2TurbulentScheme")
    .def(py::init<>())                 // default constructor
    .def(py::init<Process::Pointer>()) // constructor passing a turbulence model
    ;

}

} // namespace Python.

} // Namespace Kratos
