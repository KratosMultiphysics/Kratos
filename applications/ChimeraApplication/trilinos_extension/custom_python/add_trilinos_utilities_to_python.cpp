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

#include "custom_utilities/trilinos_chimera_fractional_step_settings.h"


namespace Kratos {
namespace Python {

void AddTrilinosUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using TrilinosSparseSpace = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using UblasLocalSpace = UblasSpace<double, Matrix, Vector>;
    using TrilinosLinearSolver = LinearSolver<TrilinosSparseSpace, UblasLocalSpace>;

    using BaseSolverSettings = SolverSettings<TrilinosSparseSpace, UblasLocalSpace, TrilinosLinearSolver>;
    typedef ChimeraTrilinosFractionalStepSettings<TrilinosSparseSpace,UblasLocalSpace,TrilinosLinearSolver> ChimeraTrilinosFractionalStepSettingsType;

    // // TrilinosChimeraFractionalStepSettings
    // py::class_< ChimeraTrilinosFractionalStepSettingsType, ChimeraTrilinosFractionalStepSettingsType::Pointer, BaseSolverSettings>(m,"TrilinosFractionalStepSettings")
    // .def(py::init<Epetra_MpiComm&, ModelPart&,unsigned int,unsigned int,bool,bool,bool>())
    // .def("SetStrategy",&ChimeraTrilinosFractionalStepSettingsType::SetStrategy)
    // .def("GetStrategy",&ChimeraTrilinosFractionalStepSettingsType::pGetStrategy)
    // .def("SetEchoLevel",&ChimeraTrilinosFractionalStepSettingsType::SetEchoLevel)
    // ;
}

}
}
