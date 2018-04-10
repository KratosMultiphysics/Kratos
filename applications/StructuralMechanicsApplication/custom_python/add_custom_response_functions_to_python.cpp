// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
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

void  AddCustomResponseFunctionsToPython()
{
    using namespace boost::python;

    // Response Functions
    class_<ResponseFunction, boost::noncopyable >
      ("ResponseFunction", no_init)
      .def("Initialize", &ResponseFunction::Initialize)
      .def("CalculateValue", &ResponseFunction::CalculateValue)
      .def("CalculateGradient", &ResponseFunction::CalculateGradient);

    class_<StrainEnergyResponseFunction, bases<ResponseFunction>, boost::noncopyable >
      ("StrainEnergyResponseFunction", init<ModelPart&, Parameters>());

    class_<MassResponseFunction, bases<ResponseFunction>, boost::noncopyable >
      ("MassResponseFunction", init<ModelPart&, Parameters>());

    class_<EigenfrequencyResponseFunction, bases<ResponseFunction>, boost::noncopyable >
      ("EigenfrequencyResponseFunction", init<ModelPart&, Parameters&>());

    class_<EigenfrequencyResponseFunctionLinScal, bases<ResponseFunction>, boost::noncopyable >
      ("EigenfrequencyResponseFunctionLinScal", init<ModelPart&, Parameters&>());
}

}  // namespace Python.

} // Namespace Kratos

