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
#include "custom_utilities/response_functions/strain_energy_response_function.h"
#include "custom_utilities/response_functions/mass_response_function.h"
#include "custom_utilities/input_output/universal_file_io.h"
#include "custom_utilities/input_output/vtk_file_io.h"
#include "custom_utilities/response_functions/eigenfrequency_response_function.h"
#include "custom_utilities/response_functions/eigenfrequency_response_function_lin_scal.h"
// ==============================================================================

namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;

    // ================================================================
    // For perfoming the mapping according to Vertex Morphing
    // ================================================================
    class_<MapperVertexMorphing >(m, "MapperVertexMorphing")
        .def(init<ModelPart&, Parameters>())
        .def("MapToDesignSpace", &MapperVertexMorphing::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphing::MapToGeometrySpace)
        ;
    class_<MapperVertexMorphingMatrixFree >(m, "MapperVertexMorphingMatrixFree")
        .def(init<ModelPart&, Parameters>())
        .def("MapToDesignSpace", &MapperVertexMorphingMatrixFree::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphingMatrixFree::MapToGeometrySpace)
        ;
    class_<MapperVertexMorphingImprovedIntegration >(m, "MapperVertexMorphingImprovedIntegration")
        .def(init<ModelPart&, Parameters>())
        .def("MapToDesignSpace", &MapperVertexMorphingImprovedIntegration::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphingImprovedIntegration::MapToGeometrySpace)
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
        .def("ExtractSurfaceNodes", &GeometryUtilities::ExtractSurfaceNodes)
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
    // For calculations related to response functions
    // ========================================================================
    class_<StrainEnergyResponseFunction >(m, "StrainEnergyResponseFunction")
        .def(init<ModelPart&, Parameters>())
        .def("Initialize", &StrainEnergyResponseFunction::Initialize)
        .def("CalculateValue", &StrainEnergyResponseFunction::CalculateValue)
        .def("CalculateGradient", &StrainEnergyResponseFunction::CalculateGradient)
        .def("GetValue", &StrainEnergyResponseFunction::GetValue)
        .def("GetGradient", &StrainEnergyResponseFunction::GetGradient)
        ;
    class_<MassResponseFunction >(m, "MassResponseFunction")
        .def(init<ModelPart&, Parameters>())
        .def("Initialize", &MassResponseFunction::Initialize)
        .def("CalculateValue", &MassResponseFunction::CalculateValue)
        .def("CalculateGradient", &MassResponseFunction::CalculateGradient)
        .def("GetValue", &MassResponseFunction::GetValue)
        .def("GetGradient", &MassResponseFunction::GetGradient)
        ;

   class_<EigenfrequencyResponseFunction >(m, "EigenfrequencyResponseFunction")
        .def(init<ModelPart&, Parameters&>())
        .def("Initialize", &EigenfrequencyResponseFunction::Initialize)
        .def("CalculateValue", &EigenfrequencyResponseFunction::CalculateValue)
        .def("CalculateGradient", &EigenfrequencyResponseFunction::CalculateGradient)
        .def("GetValue", &EigenfrequencyResponseFunction::GetValue)
        .def("GetGradient", &EigenfrequencyResponseFunction::GetGradient)
        ;

    class_<EigenfrequencyResponseFunctionLinScal >(m, "EigenfrequencyResponseFunctionLinScal")
        .def(init<ModelPart&, Parameters&>())
        .def("Initialize", &EigenfrequencyResponseFunctionLinScal::Initialize)
        .def("CalculateValue", &EigenfrequencyResponseFunctionLinScal::CalculateValue)
        .def("CalculateGradient", &EigenfrequencyResponseFunctionLinScal::CalculateGradient)
        .def("GetValue", &EigenfrequencyResponseFunctionLinScal::GetValue)
        .def("GetGradient", &EigenfrequencyResponseFunctionLinScal::GetGradient)
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

