//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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
#include "operations/operation.h"
#include "processes/process.h"

namespace Kratos::Python
{

namespace
{
    pybind11::list GetRegistryItemKeys(const RegistryItem& rSelf)
    {
        pybind11::list t;
        KRATOS_ERROR_IF(!rSelf.HasItems()) << "Asking for the keys of " << rSelf.Name() << "which has no subitems." << std::endl;
        for (auto it = rSelf.cbegin(); it != rSelf.cend(); ++it) {
            t.append(it->first);
        }
        return t;
    }

    pybind11::list GetRegistryKeys()
    {
        pybind11::list t;
        for (auto it = Kratos::Registry::cbegin(); it != Kratos::Registry::cend(); ++it) {
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
        .def("HasItems", &RegistryItem::HasItems)
        .def("HasValue", &RegistryItem::HasValue)
        .def("keys", &GetRegistryItemKeys)
        .def("size", &RegistryItem::size)
        .def("__iter__", [](RegistryItem& rSelf){return py::make_iterator(rSelf.KeyConstBegin(), rSelf.KeyConstEnd());}, py::keep_alive<0,1>())
        .def("__str__", PrintObject<RegistryItem>)
        ;

    py::class_<Registry, Registry::Pointer>(m, "CppRegistry")
        .def_static("HasItem", &Registry::HasItem)
        .def_static("HasItems", &Registry::HasItems)
        .def_static("HasValue", &Registry::HasValue)
        .def_static("GetItem", &Registry::GetItem, py::return_value_policy::reference)
        .def_static("GetOperation", &Registry::GetValue<Operation>, py::return_value_policy::reference)
        .def_static("GetProcess", &Registry::GetValue<Process>, py::return_value_policy::reference)
        .def_static("RemoveItem", &Registry::RemoveItem)
        .def_static("keys", &GetRegistryKeys)
        .def_static("size", &Registry::size)
        ;
}

}  // namespace Kratos::Python
