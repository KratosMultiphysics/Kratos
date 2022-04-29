// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// =================================================================================

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
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_python/add_custom_controls_to_python.h"
#include "custom_responses/shape_responses/plane_symmetry.h"


// ==============================================================================

namespace Kratos {
namespace Python {



// ==============================================================================
void  AddCustomResponsesToPython(pybind11::module& m)
{
    namespace py = pybind11;
    // ================================================================
    // 
    // ================================================================
    py::class_<PlaneSymmetry >(m, "PlaneSymmetry")
        .def(py::init<std::string, Model&, Parameters>())
        .def("Initialize", &PlaneSymmetry::Initialize)
        .def("CalculateValue", &PlaneSymmetry::CalculateValue)
        .def("CalculateGradient", &PlaneSymmetry::CalculateGradient)        
        ;               
 
}

}  // namespace Python.
} // Namespace Kratos

