//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_modelers_to_python.h"
#include "custom_modelers/voxel_mesh_generator_modeler.h"
#include "modeler/modeler.h"

namespace Kratos {
namespace Python {


void AddCustomModelersToPython(pybind11::module& m)
{
    namespace py = pybind11;
    
    py::class_<VoxelMeshGeneratorModeler, VoxelMeshGeneratorModeler::Pointer, Modeler>(m, "VoxelMeshGeneratorModeler")
        .def(py::init<Model &, Parameters>())
        ;
}


} // namespace Python.
} // Namespace Kratos
