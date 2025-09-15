//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/material_point_erase_process.h"
#include "includes/node.h"

namespace Kratos{
namespace Python{

    void  AddCustomProcessesToPython(pybind11::module& m)
    {
        namespace py = pybind11;

        py::class_<MaterialPointEraseProcess, MaterialPointEraseProcess::Pointer, Process>(m,"MaterialPointEraseProcess")
        .def(py::init<ModelPart&>());

    }

}  // namespace Python.

} // Namespace Kratos


