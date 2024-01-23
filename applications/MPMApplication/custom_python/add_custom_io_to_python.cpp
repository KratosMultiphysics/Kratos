//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Crescenzio
//
//


// System includes

// External includes

// Project includes
#include "custom_python/add_custom_io_to_python.h"
#include "includes/define.h"
#include "processes/process.h"
#include "custom_io/mpm_vtk_output.h"

namespace Kratos::Python{

    void  AddCustomIOToPython(pybind11::module& m)
    {
        namespace py = pybind11;

        py::class_<MPMVtkOutput, MPMVtkOutput::Pointer, IO>(m, "MPMVtkOutput")
            .def(py::init< ModelPart&, Parameters >())
            .def("PrintOutput", &MPMVtkOutput::PrintOutput, py::arg("output_filename")="")
            .def_static("GetDefaultParameters", &MPMVtkOutput::GetDefaultParameters);
    }

}  // namespace Kratos::Python.
