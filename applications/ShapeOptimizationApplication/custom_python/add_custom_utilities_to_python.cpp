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
#include "custom_utilities/damping/damping_utilities.h"
#include "custom_utilities/response_functions/strain_energy_response_function.h"
#include "custom_utilities/response_functions/mass_response_function.h"
#include "custom_utilities/input_output/universal_file_io.h"
#include "custom_utilities/input_output/vtk_file_io.h"
//#include "custom_utilities/response_functions/eigenfrequency_response_function.h"
//#include "custom_utilities/response_functions/eigenfrequency_response_function_lin_scal.h"
//#include "custom_utilities/response_functions/eigenfrequency_response_function_KS.h"
#include "custom_utilities/response_functions/local_stress_response_function.h"
#include "custom_utilities/response_functions/nodal_displacement_response_function.h"
#include "custom_utilities/response_functions/rework_strain_energy_response_function.h" //TODO: fusseder rename it after finishing
#include "custom_utilities/response_functions/rework_eigenfrequency_response_function.h" //TODO: fusseder rename it after finishing
#include "custom_utilities/finite_differences_utilities.h" //MFusseder
// ==============================================================================

namespace Kratos
{

namespace Python
{

inline
void CalculateGradient1(
        StructuralResponseFunction& rThisUtil,
        const Condition& rAdjointCondition,
        const Matrix& rAdjointMatrix,
        Vector& rResponseGradient,
        ProcessInfo& rProcessInfo)
{
    rThisUtil.CalculateGradient(rAdjointCondition,rAdjointMatrix,rResponseGradient,rProcessInfo);
}

inline
void CalculateGradient2(
        StructuralResponseFunction& rThisUtil,
        const Element& rAdjointElem,
        const Matrix& rAdjointMatrix,
        Vector& rResponseGradient,
        ProcessInfo& rProcessInfo)
{
    rThisUtil.CalculateGradient(rAdjointElem,rAdjointMatrix,rResponseGradient,rProcessInfo);
}
   

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

    class_<MapperVertexMorphingMatrixFree, bases<Process> >("MapperVertexMorphingMatrixFree", init<ModelPart&, Parameters>())
        .def("MapToDesignSpace", &MapperVertexMorphingMatrixFree::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphingMatrixFree::MapToGeometrySpace)
        ;

    class_<MapperVertexMorphingImprovedIntegration, bases<Process> >("MapperVertexMorphingImprovedIntegration", init<ModelPart&, Parameters>())
        .def("MapToDesignSpace", &MapperVertexMorphingImprovedIntegration::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphingImprovedIntegration::MapToGeometrySpace)
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
        ;

    // ========================================================================
    // For pre- and post-processing of geometry data
    // ========================================================================
    class_<GeometryUtilities, bases<Process> >("GeometryUtilities", init<ModelPart&>())
        .def("ComputeUnitSurfaceNormals", &GeometryUtilities::ComputeUnitSurfaceNormals)
        .def("ProjectNodalVariableOnUnitSurfaceNormals", &GeometryUtilities::ProjectNodalVariableOnUnitSurfaceNormals)
        .def("UpdateShapeChangeByInputVariable", &GeometryUtilities::UpdateShapeChangeByInputVariable)
        .def("ExtractSurfaceNodes", &GeometryUtilities::ExtractSurfaceNodes)
        ;

    // ========================================================================
    // For calculations related to response functions
    // ========================================================================
    class_<StrainEnergyResponseFunction, bases<Process> >("StrainEnergyResponseFunction", init<ModelPart&, Parameters>())
        .def("Initialize", &StrainEnergyResponseFunction::Initialize)
        .def("CalculateValue", &StrainEnergyResponseFunction::CalculateValue)
        .def("CalculateGradient", &StrainEnergyResponseFunction::CalculateGradient)
        .def("GetValue", &StrainEnergyResponseFunction::GetValue)
        .def("GetInitialValue", &StrainEnergyResponseFunction::GetInitialValue)
        .def("GetGradient", &StrainEnergyResponseFunction::GetGradient)
        ;

    class_<MassResponseFunction, bases<Process> >("MassResponseFunction", init<ModelPart&, Parameters>())
        .def("Initialize", &MassResponseFunction::Initialize)
        .def("CalculateValue", &MassResponseFunction::CalculateValue)
        .def("CalculateGradient", &MassResponseFunction::CalculateGradient)
        .def("GetValue", &MassResponseFunction::GetValue)
        .def("GetInitialValue", &MassResponseFunction::GetInitialValue)
        .def("GetGradient", &MassResponseFunction::GetGradient)
        ;
  
   /* class_<EigenfrequencyResponseFunction, bases<Process> >("EigenfrequencyResponseFunction", init<ModelPart&, Parameters&>())
        .def("initialize", &EigenfrequencyResponseFunction::initialize)
        .def("calculate_value", &EigenfrequencyResponseFunction::calculate_value)
        .def("calculate_gradient", &EigenfrequencyResponseFunction::calculate_gradient) 
        .def("get_value", &EigenfrequencyResponseFunction::get_value)
        .def("get_initial_value", &EigenfrequencyResponseFunction::get_initial_value)  
        .def("get_gradient", &EigenfrequencyResponseFunction::get_gradient)   
        ;   

    class_<EigenfrequencyResponseFunctionLinScal, bases<Process> >("EigenfrequencyResponseFunctionLinScal", init<ModelPart&, Parameters&>())
        .def("initialize", &EigenfrequencyResponseFunctionLinScal::initialize)
        .def("calculate_value", &EigenfrequencyResponseFunctionLinScal::calculate_value)
        .def("calculate_gradient", &EigenfrequencyResponseFunctionLinScal::calculate_gradient) 
        .def("get_value", &EigenfrequencyResponseFunctionLinScal::get_value)
        .def("get_initial_value", &EigenfrequencyResponseFunctionLinScal::get_initial_value)  
        .def("get_gradient", &EigenfrequencyResponseFunctionLinScal::get_gradient)   
        ;  

    class_<EigenfrequencyResponseFunctionKS, bases<Process> >("EigenfrequencyResponseFunctionKS", init<ModelPart&, Parameters&>())
        .def("initialize", &EigenfrequencyResponseFunctionKS::initialize)
        .def("calculate_value", &EigenfrequencyResponseFunctionKS::calculate_value)
        .def("calculate_gradient", &EigenfrequencyResponseFunctionKS::calculate_gradient) 
        .def("get_value", &EigenfrequencyResponseFunctionKS::get_value)
        .def("get_initial_value", &EigenfrequencyResponseFunctionKS::get_initial_value)  
        .def("get_gradient", &EigenfrequencyResponseFunctionKS::get_gradient)   
        ;*/

    class_<StructuralResponseFunction, boost::noncopyable>("StructuralResponseFunction", init<ModelPart&, Parameters&>())
        .def("Initialize", &StructuralResponseFunction::Initialize)
        .def("InitializeSolutionStep", &StructuralResponseFunction::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &StructuralResponseFunction::FinalizeSolutionStep)
        .def("Check", &StructuralResponseFunction::Check)
        .def("Clear", &StructuralResponseFunction::Clear)
        .def("CalculateGradient", CalculateGradient1)
        .def("CalculateGradient", CalculateGradient2)
        .def("CalculateFirstDerivativesGradient",
             &StructuralResponseFunction::CalculateFirstDerivativesGradient)
        .def("CalculateSecondDerivativesGradient",
             &StructuralResponseFunction::CalculateSecondDerivativesGradient)
        .def("CalculateValue", &StructuralResponseFunction::CalculateValue)
        .def("UpdateSensitivities", &StructuralResponseFunction::UpdateSensitivities);  

    class_<LocalStressResponseFunction, bases<StructuralResponseFunction>, boost::noncopyable>
      ("LocalStressResponseFunction", init<ModelPart&, Parameters&>()); 

    class_<NodalDisplacementResponseFunction, bases<StructuralResponseFunction>, boost::noncopyable>
      ("NodalDisplacementResponseFunction", init<ModelPart&, Parameters&>());    

    class_<ReworkStrainEnergyResponseFunction, bases<StructuralResponseFunction>, boost::noncopyable>
      ("ReworkStrainEnergyResponseFunction", init<ModelPart&, Parameters&>());  

    class_<ReworkEigenfrequencyResponseFunction, bases<StructuralResponseFunction>, boost::noncopyable>
      ("ReworkEigenfrequencyResponseFunction", init<ModelPart&, Parameters&>()); 

    // ================================================================
    // For Finite Differences TODO: is this needed?
    // ================================================================
    class_<FiniteDifferencesUtilities, boost::noncopyable>("FiniteDifferencesUtilities", init< >())
        .def("SetDesignVariable", &FiniteDifferencesUtilities::SetDesignVariable)
        .def("GetDesignVariable", &FiniteDifferencesUtilities::GetDesignVariable)
        .def("SetDerivedObject", &FiniteDifferencesUtilities::SetDerivedObject)
        .def("GetDerivedObject", &FiniteDifferencesUtilities::GetDerivedObject)
        .def("DisturbElementDesignVariable", &FiniteDifferencesUtilities::DisturbElementDesignVariable)
        .def("UndisturbElementDesignVariable", &FiniteDifferencesUtilities::UndisturbElementDesignVariable)
        .def("GetStressResultantBeam", &FiniteDifferencesUtilities::GetStressResultantBeam)
        .def("GetStressResultantShell", &FiniteDifferencesUtilities::GetStressResultantShell)
        .def("GetNodalDisplacement", &FiniteDifferencesUtilities::GetNodalDisplacement)
        .def("GetStrainEnergy", &FiniteDifferencesUtilities::GetStrainEnergy)
        
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

