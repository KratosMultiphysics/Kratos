//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//                   Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "operations/operation.h"
#include "add_operations_to_python.h"

namespace Kratos::Python
{

void AddOperationsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<Operation, Operation::Pointer>(m,"Operation")
    .def(py::init<>())
    .def("Create", &Operation::Create)
    .def("Execute", &Operation::Execute)
    .def("GetDefaultParameters", &Operation::GetDefaultParameters)
    .def("__str__", PrintObject<Operation>)
    ;
}

} // Namespace Kratos::Python
