//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Pooyan Dadvand
//

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/registry.h"
#include "python/add_registry_to_python.h"
#include "operations/operation.h"
#include "processes/process.h"

// System includes


namespace Kratos
{

namespace Python
{

void AddRegistryToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<Registry, Registry::Pointer>(m,"Registry")
        .def_static("GetItem", &Registry::GetItem)
        .def_static("GetOperation", &Registry::GetValue<Operation>, py::return_value_policy::reference)
        .def_static("GetProcess", &Registry::GetValue<Process>, py::return_value_policy::reference)
        .def_static("RemoveItem", &Registry::RemoveItem)
        .def_static("HasItem", &Registry::HasItem)
        // .def("__str__", PrintObject<Kernel>)
    ;
}

}  // namespace Python.
}  // Namespace Kratos
