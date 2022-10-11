//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/registry.h"
#include "add_registry_to_python.h"
#include "processes/process.h"
#include "includes/registry_value_item.h"

namespace Kratos
{

namespace Python
{

void AddRegistryToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<RegistryItem, RegistryItem::Pointer>(m, "RegistryItem")
        .def(py::init<std::string>())
        .def("GetProcess", &RegistryItem::GetValue<Process>, py::return_value_policy::reference)
        ;

    py::class_<RegistryValueItem<Process>, RegistryValueItem<Process>::Pointer, RegistryItem>(m, "RegistryProcessItem")
        .def(py::init<std::string, Process>())
        ;

    py::class_<Registry, Registry::Pointer>(m, "Registry")
        .def_static("HasItem", Registry::HasItem)
        .def_static("GetItem", Registry::GetItem, py::return_value_policy::reference)
        .def_static("RemoveItem", Registry::RemoveItem)
        ;
}

}  // namespace Python.

} // Namespace Kratos
