//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

// Internal includes
#include "add_custom_processes_to_python.h"

// Project inlcudes
#include "custom_processes/impose_mesh_motion_process.h"

namespace Kratos
{
namespace Python
{

void AddCustomProcessesToPython(pybind11::module& rModule)
{
    pybind11::class_<ImposeMeshMotionProcess,ImposeMeshMotionProcess::Pointer,Process>(rModule, "ImposeMeshMotionProcess")
        .def(pybind11::init<ModelPart&, Parameters>())
        ;
}

} // namespace Python
} // namespace Kratos