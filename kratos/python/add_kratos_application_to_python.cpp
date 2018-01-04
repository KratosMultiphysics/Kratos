//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "python/add_kratos_application_to_python.h"

namespace Kratos {
namespace Python {
using namespace pybind11;

void AddKratosApplicationToPython(pybind11::module& m) {
    class_<KratosApplication, KratosApplication::Pointer>(m,"KratosApplication")
        .def(init<std::string>())
        .def("Register", &KratosApplication::Register)
        .def("__repr__", &KratosApplication::Info);
}

}  // namespace Python.

}  // Namespace Kratos
