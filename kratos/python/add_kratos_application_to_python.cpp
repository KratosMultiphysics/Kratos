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
#include "includes/define_python.h"
#include "includes/kratos_application.h"
#include "python/add_kratos_application_to_python.h"

namespace Kratos {
namespace Python {
namespace py = pybind11;

void RegisterToPythonApplicationVariables(std::string ApplicationName)
{
    auto comp = KratosComponents<VariableData>::GetComponents();
    auto m = pybind11::module::import(ApplicationName.c_str()); //Note that this is added to KratosMultiphysics not to

    for(auto item = comp.begin(); item!=comp.end(); item++)
    {
        auto& var = (item->second);
        std::string name = item->first;

        m.attr(name.c_str()) = var;
    }
}


void AddKratosApplicationToPython(pybind11::module& m) {
    py::class_<KratosApplication, KratosApplication::Pointer>(m,"KratosApplication")
        .def(py::init<std::string>())
        .def("Register", [](KratosApplication& self){
            std::cout << "*************************************" << std::endl;
            std::cout << "application name = " << self.Name() << std::endl;
            self.Register();
            RegisterToPythonApplicationVariables(self.Name());
        }
        )
        .def("__str__", PrintObject<KratosApplication>)
        ;
}

}  // namespace Python.

}  // Namespace Kratos
