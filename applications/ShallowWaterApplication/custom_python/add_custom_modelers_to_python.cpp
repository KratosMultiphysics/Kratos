//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes


// Project includes
#include "includes/define_python.h"
#include "add_custom_modelers_to_python.h"
#include "custom_modelers/mesh_moving_modeler.h"

namespace Kratos
{

namespace Python
{

void  AddCustomModelersToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<MeshMovingModeler, MeshMovingModeler::Pointer, Modeler>(m,"MeshMovingModeler")
        .def(py::init<>())
        .def(py::init<Model&, Parameters>())
    ;
}

}  // namespace Python.

} // Namespace Kratos
