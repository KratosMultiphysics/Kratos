//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// // ------------------------------------------------------------------------------
// // External includes
// // ------------------------------------------------------------------------------
#include <pybind11/stl.h>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_python/add_custom_input_output_to_python.h"
#include "custom_input_output/container_variable_data_vtk_output.h"

// ==============================================================================

namespace Kratos
{
    namespace Python
    {

        // ==============================================================================
        void AddCustomInputOutputToPython(pybind11::module &m)
        {
            namespace py = pybind11;

            // ================================================================
            //
            // ================================================================
            py::class_<ContainerVariableDataVtkOutput>(m, "ContainerVariableDataVtkOutput")
                .def(py::init<ModelPart &, Parameters>())
                .def("WriteContainerDataToFile", &ContainerVariableDataVtkOutput::WriteContainerDataToFile<ModelPart::NodesContainerType>)
                .def("WriteContainerDataToFile", &ContainerVariableDataVtkOutput::WriteContainerDataToFile<ModelPart::ConditionsContainerType>)
                .def("WriteContainerDataToFile", &ContainerVariableDataVtkOutput::WriteContainerDataToFile<ModelPart::ElementsContainerType>);
        }

    } // namespace Python.
} // Namespace Kratos
