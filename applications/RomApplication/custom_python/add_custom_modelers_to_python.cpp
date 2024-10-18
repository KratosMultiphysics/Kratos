//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Ruben Zorrilla
//
//


// System includes

// External includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Project includes
#include "includes/define.h"
#include "modeler/modeler.h"

// Application includes
#include "custom_python/add_custom_modelers_to_python.h"
#include "custom_modelers/hrom_visualization_mesh_modeler.h"

namespace Kratos {

namespace Python {

namespace py = pybind11;

void AddCustomModelersToPython(py::module& m)
{

    py::class_<HRomVisualizationMeshModeler, typename HRomVisualizationMeshModeler::Pointer, Modeler>(m, "HRomVisualizationMeshModeler")
        .def(py::init<Model&, Parameters>())
    ;

}

} // namespace Python.
} // Namespace Kratos
