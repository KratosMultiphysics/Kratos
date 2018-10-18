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

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/optimization_utilities.h"
#include "custom_utilities/geometry_utilities.h"
#include "custom_utilities/mapping/mapper_vertex_morphing.h"
#include "custom_utilities/mapping/mapper_vertex_morphing_matrix_free.h"
#include "custom_utilities/mapping/mapper_vertex_morphing_improved_integration.h"
#include "custom_utilities/damping/damping_utilities.h"
#include "custom_utilities/mesh_controller_utilities.h"
#include "custom_utilities/input_output/universal_file_io.h"
#include "custom_utilities/input_output/vtk_file_io.h"
// ==============================================================================

namespace Kratos
{

namespace Python
{

// Overloaded functions
template<typename TMapper>
inline void MapScalar(TMapper& mapper,
                const Variable< array_1d<double, 3> >& origin_variable,
                const Variable< array_1d<double, 3> >& destination_variable)
{
    mapper.Map(origin_variable, destination_variable);
}

template<typename TMapper>
inline void MapVector(TMapper& mapper,
                const Variable< double >& origin_variable,
                const Variable< double >& destination_variable)
{
    mapper.Map(origin_variable, destination_variable);
}

template<typename TMapper>
inline void InverseMapScalar(TMapper& mapper,
                       const Variable< array_1d<double, 3> >& origin_variable,
                       const Variable< array_1d<double, 3> >& destination_variable)
{
    mapper.InverseMap(origin_variable, destination_variable);
}
template<typename TMapper>
inline void InverseMapVector(TMapper& mapper,
                       const Variable< double >& origin_variable,
                       const Variable< double >& destination_variable)
{
    mapper.InverseMap(origin_variable, destination_variable);
}

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;

    // ================================================================
    // For perfoming the mapping according to Vertex Morphing
    // ================================================================
    class_<MapperVertexMorphing >(m, "MapperVertexMorphing")
        .def(init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphing::Initialize)
        .def("Map", MapScalar<MapperVertexMorphing>)
        .def("Map", MapVector<MapperVertexMorphing>)
        .def("InverseMap", InverseMapScalar<MapperVertexMorphing>)
        .def("InverseMap", InverseMapVector<MapperVertexMorphing>)
        ;
    class_<MapperVertexMorphingMatrixFree >(m, "MapperVertexMorphingMatrixFree")
        .def(init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphingMatrixFree::Initialize)
        .def("Map", MapScalar<MapperVertexMorphingMatrixFree>)
        .def("Map", MapVector<MapperVertexMorphingMatrixFree>)
        .def("InverseMap", InverseMapScalar<MapperVertexMorphingMatrixFree>)
        .def("InverseMap", InverseMapVector<MapperVertexMorphingMatrixFree>)
        ;
    class_<MapperVertexMorphingImprovedIntegration >(m, "MapperVertexMorphingImprovedIntegration")
        .def(init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphingImprovedIntegration::Initialize)
        .def("Map", MapScalar<MapperVertexMorphingImprovedIntegration>)
        .def("Map", MapVector<MapperVertexMorphingImprovedIntegration>)
        .def("InverseMap", InverseMapScalar<MapperVertexMorphingImprovedIntegration>)
        .def("InverseMap", InverseMapVector<MapperVertexMorphingImprovedIntegration>)
        ;

    // ================================================================
    // For a possible damping of nodal variables
    // ================================================================
    class_<DampingUtilities >(m, "DampingUtilities")
        .def(init<ModelPart&, pybind11::dict, Parameters>())
        .def("DampNodalVariable", &DampingUtilities::DampNodalVariable)
        ;

    // ========================================================================
    // For performing individual steps of an optimization algorithm
    // ========================================================================
    class_<OptimizationUtilities >(m, "OptimizationUtilities")
        .def(init<ModelPart&, Parameters>())
        // ----------------------------------------------------------------
        // For running unconstrained descent methods
        // ----------------------------------------------------------------
        .def("ComputeSearchDirectionSteepestDescent", &OptimizationUtilities::ComputeSearchDirectionSteepestDescent)
        // ----------------------------------------------------------------
        // For running penalized projection method
        // ----------------------------------------------------------------
        .def("ComputeProjectedSearchDirection", &OptimizationUtilities::ComputeProjectedSearchDirection)
        .def("CorrectProjectedSearchDirection", &OptimizationUtilities::CorrectProjectedSearchDirection)
        // ----------------------------------------------------------------
        // General optimization operations
        // ----------------------------------------------------------------
        .def("ComputeControlPointUpdate", &OptimizationUtilities::ComputeControlPointUpdate)
        .def("AddFirstVariableToSecondVariable", &OptimizationUtilities::AddFirstVariableToSecondVariable)
        ;

    // ========================================================================
    // For pre- and post-processing of geometry data
    // ========================================================================
    class_<GeometryUtilities >(m, "GeometryUtilities")
        .def(init<ModelPart&>())
        .def("ComputeUnitSurfaceNormals", &GeometryUtilities::ComputeUnitSurfaceNormals)
        .def("ProjectNodalVariableOnUnitSurfaceNormals", &GeometryUtilities::ProjectNodalVariableOnUnitSurfaceNormals)
        .def("ExtractBoundaryNodes", &GeometryUtilities::ExtractBoundaryNodes)
        ;

    // ========================================================================
    // For mesh handling
    // ========================================================================
    class_<MeshControllerUtilities >(m, "MeshControllerUtilities")
        .def(init<ModelPart&>())
        .def("UpdateMeshAccordingInputVariable", &MeshControllerUtilities::UpdateMeshAccordingInputVariable)
        .def("LogMeshChangeAccordingInputVariable", &MeshControllerUtilities::LogMeshChangeAccordingInputVariable)
        .def("SetMeshToReferenceMesh", &MeshControllerUtilities::SetMeshToReferenceMesh)
        .def("SetReferenceMeshToMesh", &MeshControllerUtilities::SetReferenceMeshToMesh)
        .def("SetDeformationVariablesToZero", &MeshControllerUtilities::SetDeformationVariablesToZero)
        ;

    // ========================================================================
    // For input / output
    // ========================================================================
    class_<UniversalFileIO >(m, "UniversalFileIO")
        .def(init<ModelPart&, std::string, std::string, Parameters>())
        .def("InitializeLogging", &UniversalFileIO::InitializeLogging)
        .def("LogNodalResults", &UniversalFileIO::LogNodalResults)
        ;
    class_<VTKFileIO >(m, "VTKFileIO")
        .def(init<ModelPart&, std::string, std::string, Parameters>())
        .def("InitializeLogging", &VTKFileIO::InitializeLogging)
        .def("LogNodalResults", &VTKFileIO::LogNodalResults)
        ;
}


}  // namespace Python.

} // Namespace Kratos

