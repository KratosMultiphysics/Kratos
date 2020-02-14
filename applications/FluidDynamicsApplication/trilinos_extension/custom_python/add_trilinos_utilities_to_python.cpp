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
#include "containers/variable.h"
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "processes/process.h"
#include "spaces/ublas_space.h"

// TrilinosApplication dependencies
#include "trilinos_space.h"

#include "custom_utilities/trilinos_fractional_step_settings.h"
#include "custom_utilities/trilinos_fractional_step_settings_periodic.h"

namespace Kratos {
namespace Python {

void AddTrilinosUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using TrilinosSparseSpace = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using UblasLocalSpace = UblasSpace<double, Matrix, Vector>;
    using TrilinosLinearSolver = LinearSolver<TrilinosSparseSpace, UblasLocalSpace>;

    using BaseSolverSettings = SolverSettings<TrilinosSparseSpace, UblasLocalSpace, TrilinosLinearSolver>;
    typedef void (BaseSolverSettings::*BuildTurbulenceModel)(BaseSolverSettings::TurbulenceModelLabel const&, typename TrilinosLinearSolver::Pointer, const double, const unsigned int);
    typedef void (BaseSolverSettings::*PassTurbulenceModel)(Process::Pointer);
    BuildTurbulenceModel set_turbulence_model_by_build_overload = &BaseSolverSettings::SetTurbulenceModel;
    PassTurbulenceModel set_turbulence_model_by_pass_overload = &BaseSolverSettings::SetTurbulenceModel;

    // Note: this class is just here to provide a basis for derived classes. It has no constructor and should not be creable from python.
    py::class_< BaseSolverSettings, typename BaseSolverSettings::Pointer >(m,"BaseSolverSettings" )
    .def("SetTurbulenceModel",set_turbulence_model_by_build_overload)
    .def("SetTurbulenceModel",set_turbulence_model_by_pass_overload)
    ;

    py::enum_<BaseSolverSettings::StrategyLabel>(m,"TrilinosStrategyLabel")
    .value("Velocity",BaseSolverSettings::Velocity)
    .value("Pressure",BaseSolverSettings::Pressure)
    ;

    py::enum_<BaseSolverSettings::TurbulenceModelLabel>(m,"TrilinosTurbulenceModelLabel")
    .value("SpalartAllmaras",BaseSolverSettings::SpalartAllmaras)
    ;

    using FractionalStepSettings = TrilinosFractionalStepSettings<TrilinosSparseSpace,UblasLocalSpace,TrilinosLinearSolver>;
    typedef void (FractionalStepSettings::*SetStrategyByParamsType)(FractionalStepSettings::StrategyLabel const&,TrilinosLinearSolver::Pointer,const double,const unsigned int);
    SetStrategyByParamsType ThisSetStrategyOverload = &FractionalStepSettings::SetStrategy;

    py::class_< FractionalStepSettings, typename FractionalStepSettings::Pointer, BaseSolverSettings>(m,"TrilinosFractionalStepSettings")
    .def(py::init<Epetra_MpiComm&,ModelPart&,unsigned int,unsigned int,bool,bool,bool>())
    .def("SetStrategy",ThisSetStrategyOverload)
    .def("GetStrategy",&FractionalStepSettings::pGetStrategy)
    .def("SetEchoLevel",&FractionalStepSettings::SetEchoLevel)
    ;

    using FractionalStepSettingsPeriodic = TrilinosFractionalStepSettingsPeriodic<TrilinosSparseSpace,UblasLocalSpace,TrilinosLinearSolver>;
    typedef void (FractionalStepSettingsPeriodic::*SetStrategyByParamsPeriodicType)(BaseSolverSettings::StrategyLabel const&,TrilinosLinearSolver::Pointer,const double,const unsigned int);
    SetStrategyByParamsPeriodicType ThatSetStrategyOverload = &FractionalStepSettingsPeriodic::SetStrategy;

    py::class_< FractionalStepSettingsPeriodic, typename FractionalStepSettingsPeriodic::Pointer, BaseSolverSettings>(m,"TrilinosFractionalStepSettingsPeriodic")
    .def(py::init<Epetra_MpiComm&,ModelPart&,unsigned int,unsigned int,bool,bool,bool,const Kratos::Variable<int>&>())
    .def("SetStrategy",ThatSetStrategyOverload)
    .def("GetStrategy",&FractionalStepSettingsPeriodic::pGetStrategy)
    .def("SetEchoLevel",&FractionalStepSettingsPeriodic::SetEchoLevel)
    ;
}

}
}
