//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//                   Carlos Roig
//                   Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/registry.h"
#include "add_registry_to_python.h"
#include "includes/registry_value_item.h"
#include "operations/operation.h"
#include "processes/process.h"

namespace Kratos
{

namespace Python
{

void AddRegistryToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<Registry, Registry::Pointer>(m, "Registry")
        .def_static("HasItem", &Registry::HasItem)
        .def_static("GetItem", &Registry::GetItem, py::return_value_policy::reference)
        .def_static("GetOperation", &Registry::GetValue<Operation>, py::return_value_policy::reference)
        .def_static("GetProcess", &Registry::GetValue<Process>, py::return_value_policy::reference)
        .def_static("RemoveItem", &Registry::RemoveItem)
        ;
}

}  // namespace Python.

} // Namespace Kratos
