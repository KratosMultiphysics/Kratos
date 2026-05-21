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
//                   Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes
#include "add_accessors_to_python.h"
#include "includes/define_python.h"
#include "includes/accessor.h"
#include "includes/properties.h"

using AccessorBindType = std::unique_ptr<Kratos::Accessor>;

PYBIND11_MAKE_OPAQUE(AccessorBindType);

namespace Kratos::Python
{

namespace py = pybind11;

void AddAccessorsToPython(py::module& m)
{
    py::class_<AccessorBindType>(m, "Accessor")
        .def_static("Create", []() { 
            return std::make_unique<Accessor>();
        })
        ;
}

}  // namespace Kratos::Python.

