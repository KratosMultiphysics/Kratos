// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// // ------------------------------------------------------------------------------
// // External includes
// // ------------------------------------------------------------------------------
#include <pybind11/stl.h>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "response_functions/adjoint_response_function.h"

// ------------------------------------------------------------------------------
// Application includes
// ------------------------------------------------------------------------------
#include "custom_response_functions/gauss_point_kreisselmeier_aggregation_response_function.h"

// ------------------------------------------------------------------------------
// Base header include
// ------------------------------------------------------------------------------
#include "add_custom_response_functions_to_python.h"

// ==============================================================================

namespace Kratos {
namespace Python {

// ==============================================================================
void  AddCustomResponseFunctionsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<
        GaussPointKreisselmeierAggregationResponseFunction,
        GaussPointKreisselmeierAggregationResponseFunction::Pointer,
        AdjointResponseFunction>(m,"GaussPointKreisselmeierAggregationResponseFunction")
        .def(py::init<Parameters, ModelPart&>());

}

}  // namespace Python.
} // Namespace Kratos

