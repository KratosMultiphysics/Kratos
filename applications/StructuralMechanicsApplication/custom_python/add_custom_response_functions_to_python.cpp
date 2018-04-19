// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:        BSD License
//	                license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_response_functions_to_python.h"

//Utilities

//Response Functions
#include "custom_response_functions/response_utilities/response_function.h"
#include "custom_response_functions/response_utilities/strain_energy_response_function.h"
#include "custom_response_functions/response_utilities/mass_response_function.h"
#include "custom_response_functions/response_utilities/eigenfrequency_response_function.h"
#include "custom_response_functions/response_utilities/eigenfrequency_response_function_lin_scal.h"


namespace Kratos
{
namespace Python
{

using namespace pybind11;

void  AddCustomResponseFunctionUtilitiesToPython(pybind11::module& m)
{

    // Response Functions
    class_<ResponseFunctionUtility, ResponseFunctionUtility::Pointer>(m, "ResponseFunctionUtility")
      .def("Initialize", &ResponseFunctionUtility::Initialize)
      .def("CalculateValue", &ResponseFunctionUtility::CalculateValue)
      .def("CalculateGradient", &ResponseFunctionUtility::CalculateGradient);

    class_<StrainEnergyResponseFunctionUtility, StrainEnergyResponseFunctionUtility::Pointer, ResponseFunctionUtility >
      (m, "StrainEnergyResponseFunctionUtility")
      .def(init<ModelPart&, Parameters>());

    class_<MassResponseFunctionUtility, MassResponseFunctionUtility::Pointer, ResponseFunctionUtility >
      (m, "MassResponseFunctionUtility")
      .def(init<ModelPart&, Parameters>());

    class_<EigenfrequencyResponseFunctionUtility, EigenfrequencyResponseFunctionUtility::Pointer, ResponseFunctionUtility >
      (m, "EigenfrequencyResponseFunctionUtility")
      .def(init<ModelPart&, Parameters>());

    class_<EigenfrequencyResponseFunctionLinScalUtility, EigenfrequencyResponseFunctionLinScalUtility::Pointer, ResponseFunctionUtility >
      (m, "EigenfrequencyResponseFunctionLinScalUtility")
      .def(init<ModelPart&, Parameters>());
}

}  // namespace Python.

} // Namespace Kratos

