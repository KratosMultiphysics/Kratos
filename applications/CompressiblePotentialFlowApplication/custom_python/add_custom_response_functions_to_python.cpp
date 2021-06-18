//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//
//  Main authors:    Marc Nuñez, based on Armin Geiser and Martin Fusseder work
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_response_functions_to_python.h"

// Response Functions
#include "custom_response_functions/adjoint_lift_response_function_coordinates_jump.h"
#include "custom_response_functions/adjoint_far_field_lift_response_function.h"


namespace Kratos {
namespace Python {

void  AddCustomResponseFunctionUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Response Functions
    py::class_<AdjointLiftJumpCoordinatesResponseFunction, AdjointLiftJumpCoordinatesResponseFunction::Pointer, AdjointResponseFunction>
        (m, "AdjointLiftJumpCoordinatesResponseFunction")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<AdjointLiftFarFieldResponseFunction, AdjointLiftFarFieldResponseFunction::Pointer, AdjointResponseFunction>
        (m, "AdjointLiftFarFieldResponseFunction")
        .def(py::init<ModelPart&, Parameters>());

}

}  // namespace Python.
} // Namespace Kratos

