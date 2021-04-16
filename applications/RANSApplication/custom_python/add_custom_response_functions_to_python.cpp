//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "response_functions/adjoint_response_function.h"

// Application includes
#include "custom_response_functions/windowed_frequency_bin_component_response_function.h"

// Include base h
#include "custom_python/add_custom_response_functions_to_python.h"


namespace Kratos
{
namespace Python
{

void AddCustomResponseFunctionsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<
        WindowedFrequencyBinComponentResponseFunction<2>,
        WindowedFrequencyBinComponentResponseFunction<2>::Pointer,
        AdjointResponseFunction>(m,"WindowedFrequencyBinComponentResponseFunction2D")
        .def(py::init<Parameters, ModelPart&>());

    py::class_<
        WindowedFrequencyBinComponentResponseFunction<3>,
        WindowedFrequencyBinComponentResponseFunction<3>::Pointer,
        AdjointResponseFunction>(m,"WindowedFrequencyBinComponentResponseFunction3D")
        .def(py::init<Parameters, ModelPart&>());

}

} // namespace Python

} // namespace Kratos
