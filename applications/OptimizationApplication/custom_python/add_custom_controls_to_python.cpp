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
#include "custom_controls/shape_controls/vertex_morphing/explicit_vertex_morphing.h"


// ==============================================================================

namespace Kratos {
namespace Python {



// ==============================================================================
void  AddCustomControlsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // ================================================================
    // 
    // ================================================================
    py::class_<ExplicitVertexMorphing >(m, "ExplicitVertexMorphing")
        .def(py::init<std::string, Model&, Parameters>())
        .def("Initialize", &ExplicitVertexMorphing::Initialize)
        .def("Update", &ExplicitVertexMorphing::Update)
        .def("MapControlUpdate", &ExplicitVertexMorphing::MapControlUpdate)
        .def("MapFirstDerivative", &ExplicitVertexMorphing::MapFirstDerivative)
        ;
 
}

}  // namespace Python.
} // Namespace Kratos

