// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:        BSD License
//	                license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes
#include "custom_python/add_custom_response_functions_to_python.h"

#include "boost/numeric/ublas/vector.hpp"
#include "spaces/ublas_space.h"

// Response Functions
// #include "custom_response_functions/adjoint_structural_response_function.h"
#include "custom_response_functions/adjoint_lift_response_function_coordinates_jump.h"


namespace Kratos {
namespace Python {

void  AddCustomResponseFunctionUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Response Functions
    py::class_<AdjointLiftJumpCoordinatesResponseFunction, AdjointLiftJumpCoordinatesResponseFunction::Pointer, AdjointResponseFunction>
        (m, "AdjointLiftJumpCoordinatesResponseFunction")
        .def(py::init<ModelPart&, Parameters>());

}

}  // namespace Python.
} // Namespace Kratos

