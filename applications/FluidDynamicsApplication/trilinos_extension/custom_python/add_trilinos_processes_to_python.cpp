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
#include "processes/process.h"
#include "spaces/ublas_space.h"

// TrilinosApplication dependencies
#include "trilinos_space.h"

// FluidDynamics trilinos extensions
#include "custom_processes/trilinos_spalart_allmaras_turbulence_model.h"
#include "custom_processes/trilinos_stokes_initialization_process.h"

namespace Kratos {
namespace Python {

void AddTrilinosProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using TrilinosSparseSpace = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using UblasLocalSpace = UblasSpace<double, Matrix, Vector>;
    using TrilinosLinearSolver = LinearSolver<TrilinosSparseSpace, UblasLocalSpace>;

    // Turbulence models
    using SpalartAllmarasProcess = TrilinosSpalartAllmarasTurbulenceModel<TrilinosSparseSpace, UblasLocalSpace, TrilinosLinearSolver>;
    py::class_< SpalartAllmarasProcess, typename SpalartAllmarasProcess::Pointer, Process>(m, "TrilinosSpalartAllmarasTurbulenceModel")
    .def(py::init < Epetra_MpiComm&, ModelPart&, typename TrilinosLinearSolver::Pointer, unsigned int, double, unsigned int, bool, unsigned int>())
    .def("ActivateDES", &SpalartAllmarasProcess::ActivateDES)
    .def("AdaptForFractionalStep", &SpalartAllmarasProcess::AdaptForFractionalStep)
    ;

    // Stokes initialization processes
    using StokesInitializationProcess = TrilinosStokesInitializationProcess<TrilinosSparseSpace, UblasLocalSpace, TrilinosLinearSolver>;
    py::class_< StokesInitializationProcess, typename StokesInitializationProcess::Pointer, Process>(m,"TrilinosStokesInitializationProcess")
    .def(py::init<Epetra_MpiComm&, ModelPart&, typename TrilinosLinearSolver::Pointer, unsigned int, const Kratos::Variable<int>& >())
    .def("SetConditions",&StokesInitializationProcess::SetConditions)
    ;
}

}
}
