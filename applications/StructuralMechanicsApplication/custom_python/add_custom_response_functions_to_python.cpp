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

void  AddCustomResponseFunctionsToPython(pybind11::module& m)
{

    // Response Functions
    class_<ResponseFunction, ResponseFunction::Pointer>(m, "ResponseFunction")
      .def("Initialize", &ResponseFunction::Initialize)
      .def("CalculateValue", &ResponseFunction::CalculateValue)
      .def("CalculateGradient", &ResponseFunction::CalculateGradient);

    class_<StrainEnergyResponseFunction, StrainEnergyResponseFunction::Pointer, ResponseFunction >
      (m, "StrainEnergyResponseFunction")
      .def(init<ModelPart&, Parameters>());

    class_<MassResponseFunction, MassResponseFunction::Pointer, ResponseFunction >
      (m, "MassResponseFunction")
      .def(init<ModelPart&, Parameters>());

    class_<EigenfrequencyResponseFunction, EigenfrequencyResponseFunction::Pointer, ResponseFunction >
      (m, "EigenfrequencyResponseFunction")
      .def(init<ModelPart&, Parameters>());

    class_<EigenfrequencyResponseFunctionLinScal, EigenfrequencyResponseFunctionLinScal::Pointer, ResponseFunction >
      (m, "EigenfrequencyResponseFunctionLinScal")
      .def(init<ModelPart&, Parameters>());
}

}  // namespace Python.

} // Namespace Kratos

