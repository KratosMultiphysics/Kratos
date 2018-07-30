// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "response_functions/adjoint_response_function.h"

// Application includes
#include "custom_response_functions/drag_response_function.h"
#include "custom_python/add_custom_response_functions_to_python.h"

namespace Kratos
{
namespace Python
{
using namespace pybind11;

void AddCustomResponseFunctionsToPython(pybind11::module& m)
{
    class_<SensitivityBuilder>(m, "SensitivityBuilder")
        .def(init<Parameters&, ModelPart&, AdjointResponseFunction::Pointer>())
        .def("Initialize", &SensitivityBuilder::Initialize)
        .def("UpdateSensitivities", &SensitivityBuilder::UpdateSensitivities);

    class_<
        DragResponseFunction<2>,
        DragResponseFunction<2>::Pointer,
        AdjointResponseFunction>(m,"DragResponseFunction2D")
        .def(init<Parameters&, ModelPart&>());

    class_<
        DragResponseFunction<3>,
        DragResponseFunction<3>::Pointer,
        AdjointResponseFunction>(m,"DragResponseFunction3D")
        .def(init<Parameters&, ModelPart&>());

}

} // namespace Python

} // namespace Kratos
