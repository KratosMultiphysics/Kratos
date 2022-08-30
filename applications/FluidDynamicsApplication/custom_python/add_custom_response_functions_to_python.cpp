// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "response_functions/adjoint_response_function.h"

// Application includes
#include "custom_response_functions/drag_response_function.h"
#include "custom_response_functions/velocity_pressure_norm_square_response_function.h"
#include "custom_python/add_custom_response_functions_to_python.h"

namespace Kratos
{

template<unsigned TDim>
class PyDragResponseFunction : public DragResponseFunction<TDim>
{
    public:
        using DragResponseFunction<TDim>::DragResponseFunction;
        void Initialize() override {
            PYBIND11_OVERLOAD(void, DragResponseFunction<TDim>, Initialize, );
        }
        void InitializeSolutionStep() override {
            PYBIND11_OVERLOAD(void, DragResponseFunction<TDim>, InitializeSolutionStep, );
        }
        void FinalizeSolutionStep() override {
            PYBIND11_OVERLOAD(void, DragResponseFunction<TDim>, FinalizeSolutionStep, );
        }
        double CalculateValue(ModelPart& rModelPart) override {
            PYBIND11_OVERLOAD_PURE(double, DragResponseFunction<TDim>, CalculateValue, rModelPart);
        }
};


namespace Python
{

void AddCustomResponseFunctionsToPython(pybind11::module& m)
{
    namespace py = pybind11;
    py::class_<
        DragResponseFunction<2>,
        PyDragResponseFunction<2>,
        DragResponseFunction<2>::Pointer,
        AdjointResponseFunction>(m,"DragResponseFunction2D")
        .def(py::init<Parameters, ModelPart&>());

    py::class_<
        DragResponseFunction<3>,
        PyDragResponseFunction<3>,
        DragResponseFunction<3>::Pointer,
        AdjointResponseFunction>(m,"DragResponseFunction3D")
        .def(py::init<Parameters, ModelPart&>());

    py::class_<
        VelocityPressureNormSquareResponseFunction,
        VelocityPressureNormSquareResponseFunction::Pointer,
        AdjointResponseFunction>(m,"VelocityPressureNormSquareResponseFunction")
        .def(py::init<Parameters, Model&>());

}

} // namespace Python

} // namespace Kratos
