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

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>

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

#if defined(EMPIRE_NURBS_VERTEX_MORPHING)
    #include "custom_utilities/mapping/mapper_empire_nurbs.h"
#endif

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


void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;

    // ================================================================
    // For perfoming the mapping according to Vertex Morphing
    // ================================================================
    class_<MapperVertexMorphing, bases<Process> >("MapperVertexMorphing", init<ModelPart&, Parameters>())
        .def("MapToDesignSpace", &MapperVertexMorphing::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphing::MapToGeometrySpace)
        ;
 
    #if defined(EMPIRE_NURBS_VERTEX_MORPHING)
        class_<MapperEmpireNURBS, bases<Process> >("MapperEmpireNURBS", init<ModelPart&, Parameters&>())
            .def("MapToDesignSpace", &MapperEmpireNURBS::MapToDesignSpace)
            .def("MapToGeometrySpace", &MapperEmpireNURBS::MapToGeometrySpace)
            ;
    #endif

    class_<MapperVertexMorphingMatrixFree, bases<Process> >("MapperVertexMorphingMatrixFree", init<ModelPart&, Parameters>())
        .def("MapToDesignSpace", &MapperVertexMorphingMatrixFree::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphingMatrixFree::MapToGeometrySpace)
        ;
    class_<MapperVertexMorphingImprovedIntegration, bases<Process> >("MapperVertexMorphingImprovedIntegration", init<ModelPart&, Parameters>())
        .def("MapToDesignSpace", &MapperVertexMorphingImprovedIntegration::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphingImprovedIntegration::MapToGeometrySpace)
        ;

    class_<MapperVertexMorphing, bases<Process> >("MapperVertexMorphing", init<ModelPart&, Parameters&>())
        .def("MapToDesignSpace", &MapperVertexMorphing::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphing::MapToGeometrySpace)
        ;        

    // ================================================================
    // For a possible damping of nodal variables
    // ================================================================
    class_<DampingUtilities, bases<Process> >("DampingUtilities", init<ModelPart&, boost::python::dict, Parameters>())
        .def("DampNodalVariable", &DampingUtilities::DampNodalVariable)
        ;

    // ========================================================================
    // For performing individual steps of an optimization algorithm
    // ========================================================================
    class_<OptimizationUtilities, bases<Process> >("OptimizationUtilities", init<ModelPart&, Parameters>())
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
        .def("UpdateControlPointChangeByInputVariable", &OptimizationUtilities::UpdateControlPointChangeByInputVariable)        
        // ----------------------------------------------------------------
        // Adjoint extract design surface shape sensitivities
        // ----------------------------------------------------------------
        .def("GetAdjointDesignSurfaceSensitivities", &OptimizationUtilities::GetAdjointDesignSurfaceShapeSensitivities)
        ;

    // ========================================================================
    // For pre- and post-processing of geometry data
    // ========================================================================
    class_<GeometryUtilities, bases<Process> >("GeometryUtilities", init<ModelPart&>())
        .def("ComputeUnitSurfaceNormals", &GeometryUtilities::ComputeUnitSurfaceNormals)
        .def("ProjectNodalVariableOnUnitSurfaceNormals", &GeometryUtilities::ProjectNodalVariableOnUnitSurfaceNormals)
        .def("ExtractSurfaceNodes", &GeometryUtilities::ExtractSurfaceNodes)
        ;

    // ========================================================================
    // For mesh handling
    // ========================================================================
    class_<MeshControllerUtilities, bases<Process> >("MeshControllerUtilities", init<ModelPart&>())
        .def("UpdateMeshAccordingInputVariable", &MeshControllerUtilities::UpdateMeshAccordingInputVariable)
        .def("LogMeshChangeAccordingInputVariable", &MeshControllerUtilities::LogMeshChangeAccordingInputVariable)
        .def("SetMeshToReferenceMesh", &MeshControllerUtilities::SetMeshToReferenceMesh)
        .def("SetReferenceMeshToMesh", &MeshControllerUtilities::SetReferenceMeshToMesh)
        .def("SetDeformationVariablesToZero", &MeshControllerUtilities::SetDeformationVariablesToZero)
        ;

    // ========================================================================
    // For calculations related to response functions
    // ========================================================================
    class_<StrainEnergyResponseFunction, bases<Process> >("StrainEnergyResponseFunction", init<ModelPart&, Parameters>())
        .def("Initialize", &StrainEnergyResponseFunction::Initialize)
        .def("CalculateValue", &StrainEnergyResponseFunction::CalculateValue)
        .def("CalculateGradient", &StrainEnergyResponseFunction::CalculateGradient)
        .def("GetValue", &StrainEnergyResponseFunction::GetValue)
        .def("GetGradient", &StrainEnergyResponseFunction::GetGradient)
        ;
    class_<MassResponseFunction, bases<Process> >("MassResponseFunction", init<ModelPart&, Parameters>())
        .def("Initialize", &MassResponseFunction::Initialize)
        .def("CalculateValue", &MassResponseFunction::CalculateValue)
        .def("CalculateGradient", &MassResponseFunction::CalculateGradient)
        .def("GetValue", &MassResponseFunction::GetValue)
        .def("GetGradient", &MassResponseFunction::GetGradient)
        ;

   class_<EigenfrequencyResponseFunction, bases<Process> >("EigenfrequencyResponseFunction", init<ModelPart&, Parameters&>())
        .def("Initialize", &EigenfrequencyResponseFunction::Initialize)
        .def("CalculateValue", &EigenfrequencyResponseFunction::CalculateValue)
        .def("CalculateGradient", &EigenfrequencyResponseFunction::CalculateGradient)
        .def("GetValue", &EigenfrequencyResponseFunction::GetValue)
        .def("GetGradient", &EigenfrequencyResponseFunction::GetGradient)
        ;

    class_<EigenfrequencyResponseFunctionLinScal, bases<Process> >("EigenfrequencyResponseFunctionLinScal", init<ModelPart&, Parameters&>())
        .def("Initialize", &EigenfrequencyResponseFunctionLinScal::Initialize)
        .def("CalculateValue", &EigenfrequencyResponseFunctionLinScal::CalculateValue)
        .def("CalculateGradient", &EigenfrequencyResponseFunctionLinScal::CalculateGradient)
        .def("GetValue", &EigenfrequencyResponseFunctionLinScal::GetValue)
        .def("GetGradient", &EigenfrequencyResponseFunctionLinScal::GetGradient)
        ;

    // ========================================================================
    // For input / output
    // ========================================================================
    class_<UniversalFileIO, bases<Process> >("UniversalFileIO", init<ModelPart&, std::string, std::string, Parameters>())
        .def("InitializeLogging", &UniversalFileIO::InitializeLogging)
        .def("LogNodalResults", &UniversalFileIO::LogNodalResults)
        ;
    class_<VTKFileIO, bases<Process> >("VTKFileIO", init<ModelPart&, std::string, std::string, Parameters>())
        .def("InitializeLogging", &VTKFileIO::InitializeLogging)
        .def("LogNodalResults", &VTKFileIO::LogNodalResults)
        ;
}


}  // namespace Python.

} // Namespace Kratos

