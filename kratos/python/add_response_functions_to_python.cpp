//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "add_response_functions_to_python.h"
#include "response_functions/adjoint_response_function.h"
// Deprecated!!! To be removed after porting to new response function.
#include "solving_strategies/response_functions/response_function.h"

namespace Kratos
{

namespace Python
{

class PyAdjointResponseFunction : public AdjointResponseFunction
{
    public:
        using AdjointResponseFunction::AdjointResponseFunction;
        void Initialize() override {
            PYBIND11_OVERLOAD(void, AdjointResponseFunction, Initialize, );
        }
        void InitializeSolutionStep() override {
            PYBIND11_OVERLOAD(void, AdjointResponseFunction, InitializeSolutionStep, );
        }
        void FinalizeSolutionStep() override {
            PYBIND11_OVERLOAD(void, AdjointResponseFunction, FinalizeSolutionStep, );
        }
        double CalculateValue() override {
            PYBIND11_OVERLOAD_PURE(double, AdjointResponseFunction, CalculateValue, );
        }
};

void AddResponseFunctionsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<AdjointResponseFunction, PyAdjointResponseFunction> adjoint_response_function(
        m, "AdjointResponseFunction");
    adjoint_response_function
        .def(py::init<>())
        .def("Initialize", &AdjointResponseFunction::Initialize)
        .def("InitializeSolutionStep", &AdjointResponseFunction::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &AdjointResponseFunction::FinalizeSolutionStep)
        .def("CalculateValue", &AdjointResponseFunction::CalculateValue);

      // Deprecated!!! To be removed after porting to new response function.
      py::class_<ResponseFunction, ResponseFunction::Pointer>(m,"ResponseFunction")
        .def("Initialize", &ResponseFunction::Initialize)
        .def("InitializeSolutionStep", &ResponseFunction::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &ResponseFunction::FinalizeSolutionStep)
        .def("CalculateGradient", (void (ResponseFunction::*)(Element const&,Matrix const&,Vector&,ProcessInfo const& ) const )  &ResponseFunction::CalculateGradient)
        .def("CalculateGradient", (void (ResponseFunction::*)(Condition const&,Matrix const&,Vector&,ProcessInfo const& ) const )  &ResponseFunction::CalculateGradient)
        .def("CalculateFirstDerivativesGradient", (void (ResponseFunction::*)(Element const&,Matrix const&,Vector&,ProcessInfo const& ) const )  &ResponseFunction::CalculateFirstDerivativesGradient)
        .def("CalculateFirstDerivativesGradient", (void (ResponseFunction::*)(Condition const&,Matrix const&,Vector&,ProcessInfo const& ) const )  &ResponseFunction::CalculateFirstDerivativesGradient)
        .def("CalculateSecondDerivativesGradient", (void (ResponseFunction::*)(Element const&,Matrix const&,Vector&,ProcessInfo const& ) const )  &ResponseFunction::CalculateSecondDerivativesGradient)
        .def("CalculateSecondDerivativesGradient", (void (ResponseFunction::*)(Condition const&,Matrix const&,Vector&,ProcessInfo const& ) const )  &ResponseFunction::CalculateSecondDerivativesGradient)
        .def("CalculateValue", &ResponseFunction::CalculateValue)
        .def("UpdateSensitivities", &ResponseFunction::UpdateSensitivities);
}

}  // namespace Python.

} // Namespace Kratos
