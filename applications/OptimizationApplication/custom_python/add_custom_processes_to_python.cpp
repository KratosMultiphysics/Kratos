//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes

// // External includes
#include <pybind11/stl.h>

// Project includes

// Application includes
#include "custom_processes/entity_specific_properties_process.h"

// Include base h
#include "add_custom_processes_to_python.h"

namespace Kratos {
namespace Python {

void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<EntitySpecificPropertiesProcess, EntitySpecificPropertiesProcess::Pointer, Process>(m, "EntitySpecificPropertiesProcess")
        .def(py::init<Model&, Parameters&>());
}

} // namespace Python.
} // Namespace Kratos
