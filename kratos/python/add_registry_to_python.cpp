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

namespace Kratos::Python
{

namespace
{
    pybind11::list registry_item_keys(const std::string& rItemName)
    {
        pybind11::list t;
        auto& r_registry_item = Kratos::Registry::GetItem(rItemName);
        for (auto it = r_registry_item.begin(); it != r_registry_item.end(); ++it) {
            t.append(it->first);
        }
        return t;
    }

    pybind11::list registry_keys()
    {
        pybind11::list t;
        for (auto it = Kratos::Registry::begin(); it != Kratos::Registry::end(); ++it) {
            t.append(it->first);
        }
        return t;
    }
}

void AddRegistryToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<RegistryItem, RegistryItem::Pointer>(m, "RegistryItem")
        .def("Name", &RegistryItem::Name)
        // .def("__iter__", [](RegistryItem& self){return py::make_iterator(self.begin(), self.end());}, py::keep_alive<0,1>()) // Keep RegistryItem alive while the iterator is used
        ;

    py::class_<Registry, Registry::Pointer>(m, "CppRegistry")
        .def_static("HasItem", &Registry::HasItem)
        .def_static("GetItem", &Registry::GetItem, py::return_value_policy::reference)
        .def_static("GetOperation", &Registry::GetValue<Operation>, py::return_value_policy::reference)
        .def_static("GetProcess", &Registry::GetValue<Process>, py::return_value_policy::reference)
        .def_static("RemoveItem", &Registry::RemoveItem)
        .def_static("keys", &registry_keys)
        ;
}

}  // namespace Kratos::Python
