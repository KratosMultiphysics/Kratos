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

    py::class_<AdjointResponseFunction, PyAdjointResponseFunction, AdjointResponseFunction::Pointer> adjoint_response_function(
        m, "AdjointResponseFunction");
    adjoint_response_function
        .def(py::init<>())
        .def("Initialize", &AdjointResponseFunction::Initialize)
        .def("InitializeSolutionStep", &AdjointResponseFunction::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &AdjointResponseFunction::FinalizeSolutionStep)
        .def("CalculateValue", &AdjointResponseFunction::CalculateValue);
}

}  // namespace Python.

} // Namespace Kratos
