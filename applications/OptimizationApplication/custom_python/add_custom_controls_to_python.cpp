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
#include "custom_controls/thickness_controls/helmholtz_thickness.h"
#include "custom_controls/topology_controls/helmholtz_topology.h"


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
        .def(py::init<std::string, Model&, std::vector<LinearSolverType::Pointer>&, Parameters>())
        .def("Initialize", &ImplicitVertexMorphing::Initialize)
        .def("Update", &ImplicitVertexMorphing::Update)
        .def("MapControlUpdate", &ImplicitVertexMorphing::MapControlUpdate)
        .def("MapFirstDerivative", &ImplicitVertexMorphing::MapFirstDerivative)
        .def("Finalize", &ImplicitVertexMorphing::Finalize)
        ;

    py::class_<HelmholtzThickness >(m, "HelmholtzThickness")
        .def(py::init<std::string, Model&, std::vector<LinearSolverType::Pointer>&, Parameters>())
        .def("Initialize", &HelmholtzThickness::Initialize)
        .def("Update", &HelmholtzThickness::Update)
        .def("MapControlUpdate", &HelmholtzThickness::MapControlUpdate)
        .def("MapFirstDerivative", &HelmholtzThickness::MapFirstDerivative)
        .def("Finalize", &HelmholtzThickness::Finalize)
        ;

    py::class_<HelmholtzTopology >(m, "HelmholtzTopology")
        .def(py::init<std::string, Model&, std::vector<LinearSolverType::Pointer>&, Parameters>())
        .def("Initialize", &HelmholtzTopology::Initialize)
        .def("Update", &HelmholtzTopology::Update)
        .def("MapControlUpdate", &HelmholtzTopology::MapControlUpdate)
        .def("MapFirstDerivative", &HelmholtzTopology::MapFirstDerivative)
        .def("Finalize", &HelmholtzTopology::Finalize)
        ;                  
 
}

}  // namespace Python.
} // Namespace Kratos

