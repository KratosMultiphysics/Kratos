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
#include "custom_controls/shape_controls/vertex_morphing/implicit_vertex_morphing.h"


// ==============================================================================

namespace Kratos {
namespace Python {



// ==============================================================================
void  AddCustomControlsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;    
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    // ================================================================
    // 
    // ================================================================
    py::class_<ImplicitVertexMorphing >(m, "ImplicitVertexMorphing")
        .def(py::init<std::string, Model&, LinearSolverType::Pointer, Parameters>())
        .def("Initialize", &ImplicitVertexMorphing::Initialize)
        .def("Update", &ImplicitVertexMorphing::Update)
        .def("MapControlUpdate", &ImplicitVertexMorphing::MapControlUpdate)
        .def("MapFirstDerivative", &ImplicitVertexMorphing::MapFirstDerivative)
        ;
 
}

}  // namespace Python.
} // Namespace Kratos

