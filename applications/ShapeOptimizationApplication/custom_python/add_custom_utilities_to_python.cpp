// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
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
#include "custom_utilities/damping/damping_utilities.h"
#include "custom_utilities/response_functions/strain_energy_response_function.h"
#include "custom_utilities/response_functions/mass_response_function.h"
#include "custom_utilities/input_output/universal_file_io.h"
#include "custom_utilities/input_output/vtk_file_io.h"
#include "custom_utilities/response_functions/eigenfrequency_response_function.h"
#include "custom_utilities/response_functions/eigenfrequency_response_function_lin_scal.h"
#include "custom_utilities/response_functions/eigenfrequency_response_function_KS.h"
#include "custom_utilities/response_functions/local_stress_response_function.h"
#include "custom_utilities/response_functions/rework_strain_energy_response_function.h" //fusseder rename it after finishing

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
    class_<MapperVertexMorphing, bases<Process> >("MapperVertexMorphing", init<ModelPart&, Parameters&>())
        .def("MapToDesignSpace", &MapperVertexMorphing::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphing::MapToGeometrySpace)
        ;

    class_<MapperVertexMorphingMatrixFree, bases<Process> >("MapperVertexMorphingMatrixFree", init<ModelPart&, Parameters&>())
        .def("MapToDesignSpace", &MapperVertexMorphingMatrixFree::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphingMatrixFree::MapToGeometrySpace)
        ;
    
    // ================================================================
    // For a possible damping of nodal variables
    // ================================================================
    class_<DampingUtilities, bases<Process> >("DampingUtilities", init<ModelPart&, boost::python::dict, Parameters&>())
        .def("DampNodalVariable", &DampingUtilities::DampNodalVariable)
        ;
 
    // ========================================================================
    // For performing individual steps of an optimization algorithm
    // ========================================================================
    class_<OptimizationUtilities, bases<Process> >("OptimizationUtilities", init<ModelPart&, Parameters::Pointer>())
        // ----------------------------------------------------------------
        // For running unconstrained descent methods
        // ----------------------------------------------------------------
        .def("compute_search_direction_steepest_descent", &OptimizationUtilities::compute_search_direction_steepest_descent)
        // ----------------------------------------------------------------
        // For running penalized projection method
        // ----------------------------------------------------------------
        .def("compute_projected_search_direction", &OptimizationUtilities::compute_projected_search_direction)
        .def("correct_projected_search_direction", &OptimizationUtilities::correct_projected_search_direction)
        // ----------------------------------------------------------------
        // General optimization operations
        // ----------------------------------------------------------------
        .def("compute_design_update", &OptimizationUtilities::compute_design_update)
        ;

    // ========================================================================
    // For pre- and post-processing of geometry data
    // ========================================================================
    class_<GeometryUtilities, bases<Process> >("GeometryUtilities", init<ModelPart&>())
        .def("compute_unit_surface_normals", &GeometryUtilities::compute_unit_surface_normals)
        .def("project_nodal_variable_on_unit_surface_normals", &GeometryUtilities::project_nodal_variable_on_unit_surface_normals)
        .def("update_coordinates_according_to_input_variable", &GeometryUtilities::update_coordinates_according_to_input_variable)
        .def("extract_surface_nodes", &GeometryUtilities::extract_surface_nodes)
        ;

    // ========================================================================
    // For calculations related to response functions
    // ========================================================================
    class_<StrainEnergyResponseFunction, bases<Process> >("StrainEnergyResponseFunction", init<ModelPart&, Parameters&>())
        .def("initialize", &StrainEnergyResponseFunction::initialize)
        .def("calculate_value", &StrainEnergyResponseFunction::calculate_value)
        .def("calculate_gradient", &StrainEnergyResponseFunction::calculate_gradient) 
        .def("get_value", &StrainEnergyResponseFunction::get_value)
        .def("get_initial_value", &StrainEnergyResponseFunction::get_initial_value)  
        .def("get_gradient", &StrainEnergyResponseFunction::get_gradient)                              
        ; 
    class_<MassResponseFunction, bases<Process> >("MassResponseFunction", init<ModelPart&, Parameters&>())
        .def("initialize", &MassResponseFunction::initialize)
        .def("calculate_value", &MassResponseFunction::calculate_value)
        .def("calculate_gradient", &MassResponseFunction::calculate_gradient)  
        .def("get_value", &MassResponseFunction::get_value)
        .def("get_initial_value", &MassResponseFunction::get_initial_value) 
        .def("get_gradient", &MassResponseFunction::get_gradient)                              
        ;   

    class_<EigenfrequencyResponseFunction, bases<Process> >("EigenfrequencyResponseFunction", init<ModelPart&, Parameters&>())
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
        ;

    class_<LocalStressResponseFunction, bases<Process> >("LocalStressResponseFunction", init<ModelPart&, Parameters&>())
        .def("Initialize", &LocalStressResponseFunction::Initialize)
        .def("CalculateValue", &LocalStressResponseFunction::CalculateValue)
        .def("UpdateSensitivities", &LocalStressResponseFunction::UpdateSensitivities) 
        .def("GetValue", &LocalStressResponseFunction::GetValue)
        .def("GetInitialValue", &LocalStressResponseFunction::GetInitialValue)                             
        ;  

    class_<ReworkStrainEnergyResponseFunction, bases<Process> >("ReworkStrainEnergyResponseFunction", init<ModelPart&, Parameters&>())
        .def("Initialize", &ReworkStrainEnergyResponseFunction::Initialize)
        .def("CalculateValue", &ReworkStrainEnergyResponseFunction::CalculateValue)
        .def("UpdateSensitivities", &ReworkStrainEnergyResponseFunction::UpdateSensitivities) 
        .def("GetValue", &ReworkStrainEnergyResponseFunction::GetValue)
        .def("GetInitialValue", &ReworkStrainEnergyResponseFunction::GetInitialValue)                             
        ;           

    // ========================================================================
    // For input / output
    // ======================================================================== 
    class_<UniversalFileIO, bases<Process> >("UniversalFileIO", init<ModelPart&, Parameters&>())
        .def("initializeLogging", &UniversalFileIO::initializeLogging)
        .def("logNodalResults", &UniversalFileIO::logNodalResults)
        ;           
     
    class_<VTKFileIO, bases<Process> >("VTKFileIO", init<ModelPart&, Parameters&>())
        .def("initializeLogging", &VTKFileIO::initializeLogging)
        .def("logNodalResults", &VTKFileIO::logNodalResults)
        ;           
}


}  // namespace Python.

} // Namespace Kratos

