//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Carlos A. Roig
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "mpi/utilities/debug_utilities.h"

namespace Kratos {
namespace Python {

// template<class TVariableTpe>
// void CheckScalarVariables(pybind11::module& m) {

// }


// CheckNonHistoricalNodeVariableConsistency()

void AddMPIDebugUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<MpiDebugUtilities>(m,"MPIDebugUtilities");

    // CheckScalarVariables<int>(m);
    // CheckScalarVariables<double>(m);
    // CheckScalarVariables<bool>(m);

    // .def("SetMPICommunicator",&MpiDebugUtilities::SetMPICommunicator)
    // ;


}

} // namespace Python
} // namespace Kratos

