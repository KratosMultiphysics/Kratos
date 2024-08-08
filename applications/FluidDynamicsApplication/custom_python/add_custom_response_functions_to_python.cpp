// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "response_functions/adjoint_response_function.h"

// Application includes
#include "custom_response_functions/drag_response_function.h"
#include "custom_response_functions/drag_frequency_response_function.h"
#include "custom_response_functions/velocity_pressure_norm_square_response_function.h"
#include "custom_python/add_custom_response_functions_to_python.h"
#include "custom_response_functions/residual_response_function.h"
#include "custom_response_functions/moment_response_function.h"
#include "custom_response_functions/domain_integrated_response_function.h"
#include "custom_response_functions/domain_integrated_3d_vector_magnitude_square_p_mean_response_function.h"

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
        DragFrequencyResponseFunction<2>,
        DragFrequencyResponseFunction<2>::Pointer,
        AdjointResponseFunction>(m,"DragFrequencyResponseFunction2D")
        .def(py::init<Parameters, ModelPart&>())
        ;

    py::class_<
        DragFrequencyResponseFunction<3>,
        DragFrequencyResponseFunction<3>::Pointer,
        AdjointResponseFunction>(m,"DragFrequencyResponseFunction3D")
        .def(py::init<Parameters, ModelPart&>())
        ;

    py::class_<
        MomentResponseFunction<2>,
        MomentResponseFunction<2>::Pointer,
        AdjointResponseFunction>(m,"MomentResponseFunction2D")
        .def(py::init<Parameters, ModelPart&>())
        ;

    py::class_<
        MomentResponseFunction<3>,
        MomentResponseFunction<3>::Pointer,
        AdjointResponseFunction>(m,"MomentResponseFunction3D")
        .def(py::init<Parameters, ModelPart&>())
        ;

    py::class_<
        DomainIntegratedResponseFunction,
        DomainIntegratedResponseFunction::Pointer,
        AdjointResponseFunction>(m,"DomainIntegratedResponseFunction")
        .def(py::init<Parameters, ModelPart&>());

    py::class_<
        DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<2>,
        DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<2>::Pointer,
        AdjointResponseFunction>(m,"DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction2D")
        .def(py::init<Parameters, ModelPart&>());

    py::class_<
        DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<3>,
        DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<3>::Pointer,
        AdjointResponseFunction>(m,"DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction3D")
        .def(py::init<Parameters, ModelPart&>());

    py::class_<
        VelocityPressureNormSquareResponseFunction,
        VelocityPressureNormSquareResponseFunction::Pointer,
        AdjointResponseFunction>(m,"VelocityPressureNormSquareResponseFunction")
        .def(py::init<Parameters, Model&>());

    py::class_<
        ResidualResponseFunction<2>,
        ResidualResponseFunction<2>::Pointer,
        AdjointResponseFunction>(m,"ResidualResponseFunction2D")
        .def(py::init<Parameters, ModelPart&>());

    py::class_<
        ResidualResponseFunction<3>,
        ResidualResponseFunction<3>::Pointer,
        AdjointResponseFunction>(m,"ResidualResponseFunction3D")
        .def(py::init<Parameters, ModelPart&>());

}

} // namespace Python

} // namespace Kratos
