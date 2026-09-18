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
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/schema.h"
#include "add_schema_to_python.h"

namespace Kratos::Python
{

void AddSchemaToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<Schema, std::shared_ptr<Schema>>(m, "Schema")
        .def_static("Get", &Schema::Get)
        .def_static("FindInstances", &Schema::FindInstances)
        .def_static("Keys", []() {const auto& map = Schema::GetByName(); return py::make_key_iterator(
            map.begin(),
            map.end()
        );})
        .def_static("Types", []() {const auto& map = Schema::GetByType(); return py::make_key_iterator(
            map.begin(),
            map.end()
        );})
        ;
}

}  // namespace Kratos::Python
