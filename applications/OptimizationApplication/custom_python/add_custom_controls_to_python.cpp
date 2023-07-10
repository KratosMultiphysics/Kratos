//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//

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
#include "custom_controls/material_controls/helmholtz_material.h"
#include "custom_controls/material_controls/helmholtz_partition.h"

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

    py::class_<HelmholtzMaterial >(m, "HelmholtzMaterial")
        .def(py::init<std::string, Model&, std::vector<LinearSolverType::Pointer>&, Parameters>())
        .def("Initialize", &HelmholtzMaterial::Initialize)
        .def("Update", &HelmholtzMaterial::Update)
        .def("MapControlUpdate", &HelmholtzMaterial::MapControlUpdate)
        .def("MapFirstDerivative", &HelmholtzMaterial::MapFirstDerivative)
        .def("Finalize", &HelmholtzMaterial::Finalize)
        ;                  
 
    py::class_<HelmholtzPartition >(m, "HelmholtzPartition")
        .def(py::init<std::string, Model&, std::vector<LinearSolverType::Pointer>&, Parameters>())
        .def("Initialize", &HelmholtzPartition::Initialize)
        .def("Update", &HelmholtzPartition::Update)
        .def("MapControlUpdate", &HelmholtzPartition::MapControlUpdate)
        .def("MapFirstDerivative", &HelmholtzPartition::MapFirstDerivative)
        .def("Finalize", &HelmholtzPartition::Finalize)
        ; 

}

}  // namespace Python.
} // Namespace Kratos

