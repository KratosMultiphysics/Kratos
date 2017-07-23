// External includes
#include <boost/python.hpp>

// Project includes

// Application includes
#include "custom_utilities/response_function.h"
#include "custom_utilities/drag_response_function.h"
#include "custom_python/add_custom_response_functions_to_python.h"

namespace Kratos
{
namespace Python
{
using namespace boost::python;

void AddCustomResponseFunctionsToPython()
{
  class_<ResponseFunction, boost::noncopyable>("ResponseFunction", init<ModelPart&, Parameters&>())
        .def("Initialize", &ResponseFunction::Initialize)
        .def("InitializeSolutionStep", &ResponseFunction::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &ResponseFunction::FinalizeSolutionStep)
        .def("Check", &ResponseFunction::Check)
        .def("CalculateFirstDerivativesGradient",
             &ResponseFunction::CalculateFirstDerivativesGradient)
        .def("CalculateSecondDerivativesGradient",
             &ResponseFunction::CalculateSecondDerivativesGradient)
        .def("CalculateSensitivityContribution", &ResponseFunction::CalculateSensitivityContribution)
        .def("CalculateValue", &ResponseFunction::CalculateValue);

    class_<DragResponseFunction<2>, bases<ResponseFunction>, boost::noncopyable>
      ("DragResponseFunction2D", init<ModelPart&, Parameters&>());

    class_<DragResponseFunction<3>, bases<ResponseFunction>, boost::noncopyable>
      ("DragResponseFunction3D", init<ModelPart&, Parameters&>());

}

} // namespace Python

} // namespace Kratos
