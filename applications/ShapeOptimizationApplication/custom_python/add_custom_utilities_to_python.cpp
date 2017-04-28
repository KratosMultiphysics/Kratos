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
#include "custom_utilities/mapping/mapper_vertex_morphing_iterative.h"
#include "custom_utilities/response_functions/strain_energy_response_function.h"
#include "custom_utilities/response_functions/mass_response_function.h"
#include "custom_utilities/input_output/universal_file_io.h"

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
    class_<MapperVertexMorphing, bases<Process> >("MapperVertexMorphing", init<ModelPart&, boost::python::dict, Parameters&>())
        .def("map_to_design_space", &MapperVertexMorphing::map_to_design_space)
        .def("map_to_geometry_space", &MapperVertexMorphing::map_to_geometry_space)
        ;

    // ================================================================
    // For perfoming the mapping according to Vertex Morphing iteratively
    // ================================================================
    class_<MapperVertexMorphingIterative, bases<Process> >("MapperVertexMorphingIterative", init<ModelPart&, boost::python::dict, Parameters&>())
        .def("map_to_design_space", &MapperVertexMorphingIterative::map_to_design_space)
        .def("map_to_geometry_space", &MapperVertexMorphingIterative::map_to_geometry_space)
        ;
 
    // ========================================================================
    // For performing individual steps of an optimization algorithm
    // ========================================================================
    class_<OptimizationUtilities, bases<Process> >("OptimizationUtilities", init<ModelPart&, Parameters&>())
        // ----------------------------------------------------------------
        // For running unconstrained descent methods
        // ----------------------------------------------------------------
        .def("compute_search_direction_steepest_descent", &OptimizationUtilities::compute_search_direction_steepest_descent)
        // ----------------------------------------------------------------
        // For running penalized projection method
        // ----------------------------------------------------------------
        .def("compute_projected_search_direction", &OptimizationUtilities::compute_projected_search_direction)
        .def("correct_projected_search_direction", &OptimizationUtilities::correct_projected_search_direction)
        .def("get_correction_scaling", &OptimizationUtilities::get_correction_scaling)
        .def("set_correction_scaling", &OptimizationUtilities::set_correction_scaling)
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
        .def("project_grad_on_unit_surface_normal", &GeometryUtilities::project_grad_on_unit_surface_normal)
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

    // ========================================================================
    // For input / output
    // ======================================================================== 
    class_<UniversalFileIO, bases<Process> >("UniversalFileIO", init<ModelPart&, Parameters&>())
        .def("initializeLogging", &UniversalFileIO::initializeLogging)
        .def("logNodalResults", &UniversalFileIO::logNodalResults)
        ;            
}


}  // namespace Python.

} // Namespace Kratos

