// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "solving_strategies/response_functions/response_function.h"

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
    class_<
        DragResponseFunction<2>,
        typename DragResponseFunction<2>::Pointer,
        ResponseFunction>(m,"DragResponseFunction2D")
        .def(init<Parameters&>());

    class_<
        DragResponseFunction<3>,
        typename DragResponseFunction<3>::Pointer,
        ResponseFunction>(m,"DragResponseFunction3D")
        .def(init<Parameters&>());

}

} // namespace Python

} // namespace Kratos
