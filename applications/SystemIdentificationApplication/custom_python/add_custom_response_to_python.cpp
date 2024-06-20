// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Talhah Ansari
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes
#include "custom_python/add_custom_response_to_python.h"

// Processes

// Response Functions

#include "custom_response/eigenfreq_response_function_utility.h"
#include "custom_response/mac_response_function_utility.h"


namespace Kratos::Python {

void  AddCustomResponseToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Response Functions

    py::class_<EigenfreqResponseFunctionUtility, EigenfreqResponseFunctionUtility::Pointer >
        (m, "EigenfreqResponseFunctionUtility")
        .def(py::init<ModelPart&, Parameters>(),  py::arg("model_part"), py::arg("sensitivity_settings"))
        .def("Initialize", &EigenfreqResponseFunctionUtility::Initialize)
        .def("CalculateValue", &EigenfreqResponseFunctionUtility::CalculateValue)
        .def("CalculateGradient", &EigenfreqResponseFunctionUtility::CalculateGradient);

    py::class_<MacResponseFunctionUtility, MacResponseFunctionUtility::Pointer >
        (m, "MacResponseFunctionUtility")
        .def(py::init<ModelPart&, Parameters>(),  py::arg("model_part"), py::arg("sensitivity_settings"))
        .def("Initialize", &MacResponseFunctionUtility::Initialize)
        .def("CalculateValue", &MacResponseFunctionUtility::CalculateValue)
        .def("CalculateGradient", &MacResponseFunctionUtility::CalculateGradient);

    // Processes
    
    // Response Functions
   
}
}  // namespace Kratos::Python.

