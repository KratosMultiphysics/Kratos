// External includes
#include "pybind11/pybind11.h"

// Project includes

// Application includes
#include "custom_utilities/response_function.h"
#include "custom_utilities/drag_response_function.h"
#include "custom_python/add_custom_response_functions_to_python.h"

namespace Kratos
{
namespace Python
{
using namespace pybind11;

void AddCustomResponseFunctionsToPython(pybind11::module& m)
{
  class_<ResponseFunction, ResponseFunction::Pointer>(m,"ResponseFunction")
        .def(init<ModelPart&, Parameters&>())
        .def("Initialize", &ResponseFunction::Initialize)
        .def("InitializeSolutionStep", &ResponseFunction::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &ResponseFunction::FinalizeSolutionStep)
        .def("Check", &ResponseFunction::Check)
        .def("Clear", &ResponseFunction::Clear)
        .def("CalculateGradient",
             &ResponseFunction::CalculateGradient)
        .def("CalculateFirstDerivativesGradient",
             &ResponseFunction::CalculateFirstDerivativesGradient)
        .def("CalculateSecondDerivativesGradient",
             &ResponseFunction::CalculateSecondDerivativesGradient)
        .def("CalculateValue", &ResponseFunction::CalculateValue)
        .def("UpdateSensitivities", &ResponseFunction::UpdateSensitivities);

    class_<
        DragResponseFunction<2>,
        typename DragResponseFunction<2>::Pointer,
        ResponseFunction>(m,"DragResponseFunction2D")
        .def(init<ModelPart&, Parameters&>());

    class_<
        DragResponseFunction<3>,
        typename DragResponseFunction<3>::Pointer,
        ResponseFunction>(m,"DragResponseFunction3D")
        .def(init<ModelPart&, Parameters&>());

}

} // namespace Python

} // namespace Kratos
