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
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "python/add_kratos_application_to_python.h"

namespace Kratos {
namespace Python {
using namespace boost::python;

void AddKratosApplicationToPython() {
    class_<KratosApplication, KratosApplication::Pointer, boost::noncopyable>(
        "KratosApplication", init<std::string>())
        .def("Register", &KratosApplication::Register)
        //.def("",&Kernel::Initialize)
        .def(self_ns::str(self));
}

}  // namespace Python.

}  // Namespace Kratos
