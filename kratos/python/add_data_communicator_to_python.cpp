//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/data_communicator.h"

namespace Kratos {

namespace Python {

void AddDataCommunicatorToPython(pybind11::module &m)
{
    namespace py = pybind11;

    py::class_<DataCommunicator, DataCommunicator::Pointer>(m,"DataCommunicator");
}

} // namespace Python.

} // Namespace Kratos
