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
#include "pybind11/pybind11.h"

// Project includes
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

#include "Epetra_FEVector.h"
#include "trilinos_space.h"

// Application includes
#include "custom_processes/rans_formulation_process.h"
#include "trilinos_extension/custom_processes/rans_trilinos_wall_distance_calculation_process.h"

// Include base h
#include "add_trilinos_processes_to_python.h"

namespace Kratos
{

namespace Python
{
void AddTrilinosProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<TrilinosRansWallDistanceCalculationProcess, TrilinosRansWallDistanceCalculationProcess::Pointer, RansFormulationProcess>(m, "TrilinosRansWallDistanceCalculationProcess")
        .def(py::init<Model&, Parameters&>());
}

} // namespace Python.
} // Namespace Kratos
