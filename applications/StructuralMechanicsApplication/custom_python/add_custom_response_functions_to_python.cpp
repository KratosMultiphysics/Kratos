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
#include "custom_python/add_custom_response_functions_to_python.h"

//Response Functions
#include "custom_response_functions/response_utilities/strain_energy_response_function_utility.h"
#include "custom_response_functions/response_utilities/mass_response_function_utility.h"
#include "custom_response_functions/response_utilities/eigenfrequency_response_function_utility.h"


namespace Kratos
{
namespace Python
{

using namespace pybind11;

void  AddCustomResponseFunctionUtilitiesToPython(pybind11::module& m)
{

    // Response Functions
    class_<StrainEnergyResponseFunctionUtility, StrainEnergyResponseFunctionUtility::Pointer >
      (m, "StrainEnergyResponseFunctionUtility")
      .def(init<ModelPart&, Parameters>())
      .def("Initialize", &StrainEnergyResponseFunctionUtility::Initialize)
      .def("CalculateValue", &StrainEnergyResponseFunctionUtility::CalculateValue)
      .def("CalculateGradient", &StrainEnergyResponseFunctionUtility::CalculateGradient);

    class_<MassResponseFunctionUtility, MassResponseFunctionUtility::Pointer >
      (m, "MassResponseFunctionUtility")
      .def(init<ModelPart&, Parameters>())
      .def("Initialize", &MassResponseFunctionUtility::Initialize)
      .def("CalculateValue", &MassResponseFunctionUtility::CalculateValue)
      .def("CalculateGradient", &MassResponseFunctionUtility::CalculateGradient);

    class_<EigenfrequencyResponseFunctionUtility, EigenfrequencyResponseFunctionUtility::Pointer >
      (m, "EigenfrequencyResponseFunctionUtility")
      .def(init<ModelPart&, Parameters>())
      .def("Initialize", &EigenfrequencyResponseFunctionUtility::Initialize)
      .def("CalculateValue", &EigenfrequencyResponseFunctionUtility::CalculateValue)
      .def("CalculateGradient", &EigenfrequencyResponseFunctionUtility::CalculateGradient);
}

}  // namespace Python.

} // Namespace Kratos

