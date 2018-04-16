// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "python/add_response_functions_to_python.h"
#include "solving_strategies/response_functions/response_function.h"

namespace Kratos
{
namespace Python
{
using namespace pybind11;

void AddResponseFunctionsToPython(pybind11::module& m)
{
      class_<ResponseFunction, ResponseFunction::Pointer>(m,"ResponseFunction")
        .def("Initialize", &ResponseFunction::Initialize)
        .def("InitializeSolutionStep", &ResponseFunction::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &ResponseFunction::FinalizeSolutionStep)
        .def("Check", &ResponseFunction::Check)
        .def("Clear", &ResponseFunction::Clear)
        .def("CalculateGradient", (void (ResponseFunction::*)(Element const&,Matrix const&,Vector&,ProcessInfo const& ) const )  &ResponseFunction::CalculateGradient)
        .def("CalculateGradient", (void (ResponseFunction::*)(Condition const&,Matrix const&,Vector&,ProcessInfo const& ) const )  &ResponseFunction::CalculateGradient)
        .def("CalculateFirstDerivativesGradient", (void (ResponseFunction::*)(Element const&,Matrix const&,Vector&,ProcessInfo const& ) const )  &ResponseFunction::CalculateFirstDerivativesGradient)
        .def("CalculateFirstDerivativesGradient", (void (ResponseFunction::*)(Condition const&,Matrix const&,Vector&,ProcessInfo const& ) const )  &ResponseFunction::CalculateFirstDerivativesGradient)
        .def("CalculateSecondDerivativesGradient", (void (ResponseFunction::*)(Element const&,Matrix const&,Vector&,ProcessInfo const& ) const )  &ResponseFunction::CalculateSecondDerivativesGradient)
        .def("CalculateSecondDerivativesGradient", (void (ResponseFunction::*)(Condition const&,Matrix const&,Vector&,ProcessInfo const& ) const )  &ResponseFunction::CalculateSecondDerivativesGradient)
        .def("CalculateValue", &ResponseFunction::CalculateValue)
        .def("UpdateSensitivities", &ResponseFunction::UpdateSensitivities);
}

} // namespace Python

} // namespace Kratos
