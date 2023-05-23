//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Fabian Meister
//

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// // ------------------------------------------------------------------------------
// // External includes
// // ------------------------------------------------------------------------------
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <vector>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_python/add_custom_input_output_to_python.h"
#include "custom_input_output/generic_vtk_output.h"

// ==============================================================================

namespace Kratos
{
    namespace Python
    {

        // ==============================================================================
        void AddCustomInputOutputToPython(pybind11::module &m)
        {
            namespace py = pybind11;

            py::class_<GenericVtkOutput>(m, "GenericVtkOutput")
                .def(py::init<>())
                .def("outputStructuredPoints", &GenericVtkOutput::outputStructuredPoints,
                      py::arg("path") ,
                      py::arg("numPointsX") ,
                      py::arg("numPointsY") ,
                      py::arg("numPointsZ") ,
                      py::arg("spaceBetweenPointsX") ,
                      py::arg("spaceBetweenPointsY") ,
                      py::arg("spaceBetweenPointsZ") ,
                      py::arg("scalarDataName"),
                      py::arg("vectorDataName"),
                      py::arg("tensorDataName"),
                      py::arg("scalarData"),
                      py::arg("vectorData"),
                      py::arg("tensorData")
                );

            // ================================================================
            //
            // ================================================================
            // py::class_<ContainerVariableDataVtkOutput>(m, "ContainerVariableDataVtkOutput")
            // .def(py::init<ModelPart&,Parameters>())
            // ;
            // py::class_<ContainerVariableDataVtkOutput>(m, "ContainerVariableDataVtkOutput")
            //     .def(py::init<ModelPart &, Parameters>())
            //     .def("TestFunction", &ContainerVariableDataVtkOutput::TestFunction)
            //     .def("printNumberType", &ContainerVariableDataVtkOutput::printNumberType<int>)
            //     .def("printNumberType", &ContainerVariableDataVtkOutput::printNumberType<double>);
            // .def("WriteContainerDataToFile", &ContainerVariableDataVtkOutput::WriteContainerDataToFile<ModelPart::NodesContainerType>)
            // .def("WriteContainerDataToFile", &ContainerVariableDataVtkOutput::WriteContainerDataToFile<ModelPart::ConditionsContainerType>)
            // .def("WriteContainerDataToFile", &ContainerVariableDataVtkOutput::WriteContainerDataToFile<ModelPart::ElementsContainerType>);
            // py::class_<ContainerVariableDataVtkOutput>(m, "GenericVtkOutput")
            //     .def(py::init<ModelPart &, Parameters>())
            //     .def("outputStructuredGrid", &GenericVtkOutput::outputStructuredGrid);
        }

    } // namespace Python.
} // Namespace Kratos
