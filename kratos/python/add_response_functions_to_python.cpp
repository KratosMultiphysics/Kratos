//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "add_response_functions_to_python.h"
#include "response_functions/adjoint_response_function.h"
#include "includes/model_part.h"

namespace Kratos::Python
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
        double CalculateValue(ModelPart& rModelPart) override {
            PYBIND11_OVERLOAD_PURE(double, AdjointResponseFunction, CalculateValue, rModelPart);
        }
};

void AddResponseFunctionsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<AdjointResponseFunction, PyAdjointResponseFunction, AdjointResponseFunction::Pointer> adjoint_response_function(
        m, "AdjointResponseFunction");
    adjoint_response_function
        .def(py::init<>())
        .def("Initialize", &AdjointResponseFunction::Initialize)
        .def("InitializeSolutionStep", &AdjointResponseFunction::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &AdjointResponseFunction::FinalizeSolutionStep)
        .def("CalculateValue", &AdjointResponseFunction::CalculateValue)
        .def("CalculateGradient", static_cast<void(AdjointResponseFunction::*)(const Condition&, const Matrix&, Vector&, const ProcessInfo&)>(&AdjointResponseFunction::CalculateGradient))
        .def("CalculateGradient", static_cast<void(AdjointResponseFunction::*)(const Element&, const Matrix&, Vector&, const ProcessInfo&)>(&AdjointResponseFunction::CalculateGradient))
        .def("CalculateFirstDerivativesGradient", static_cast<void(AdjointResponseFunction::*)(const Condition&, const Matrix&, Vector&, const ProcessInfo&)>(&AdjointResponseFunction::CalculateFirstDerivativesGradient))
        .def("CalculateFirstDerivativesGradient", static_cast<void(AdjointResponseFunction::*)(const Element&, const Matrix&, Vector&, const ProcessInfo&)>(&AdjointResponseFunction::CalculateFirstDerivativesGradient))
        .def("CalculateSecondDerivativesGradient", static_cast<void(AdjointResponseFunction::*)(const Condition&, const Matrix&, Vector&, const ProcessInfo&)>(&AdjointResponseFunction::CalculateSecondDerivativesGradient))
        .def("CalculateSecondDerivativesGradient", static_cast<void(AdjointResponseFunction::*)(const Element&, const Matrix&, Vector&, const ProcessInfo&)>(&AdjointResponseFunction::CalculateSecondDerivativesGradient))
        .def("CalculatePartialSensitivity", static_cast<void(AdjointResponseFunction::*)(Condition&, const Variable<double>&, const Matrix&, Vector&, const ProcessInfo&)>(&AdjointResponseFunction::CalculatePartialSensitivity))
        .def("CalculatePartialSensitivity", static_cast<void(AdjointResponseFunction::*)(Element&, const Variable<double>&, const Matrix&, Vector&, const ProcessInfo&)>(&AdjointResponseFunction::CalculatePartialSensitivity))
        .def("CalculatePartialSensitivity", static_cast<void(AdjointResponseFunction::*)(Condition&, const Variable<array_1d<double, 3>>&, const Matrix&, Vector&, const ProcessInfo&)>(&AdjointResponseFunction::CalculatePartialSensitivity))
        .def("CalculatePartialSensitivity", static_cast<void(AdjointResponseFunction::*)(Element&, const Variable<array_1d<double, 3>>&, const Matrix&, Vector&, const ProcessInfo&)>(&AdjointResponseFunction::CalculatePartialSensitivity))
        ;
}

}  // namespace Kratos::Python.
