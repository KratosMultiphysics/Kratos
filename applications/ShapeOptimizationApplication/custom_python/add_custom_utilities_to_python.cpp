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
#include "custom_utilities/response_functions/eigenfrequency_response_function.h"
// #include "custom_utilities/response_functions/eigenfrequency_response_function_lin_scal.h"
// #include "custom_utilities/response_functions/eigenfrequency_response_function_KS.h"
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

   class_<EigenfrequencyResponseFunction, bases<Process> >("EigenfrequencyResponseFunction", init<ModelPart&, Parameters&>())
        .def("Initialize", &EigenfrequencyResponseFunction::Initialize)
        .def("CalculateValue", &EigenfrequencyResponseFunction::CalculateValue)
        .def("CalculateGradient", &EigenfrequencyResponseFunction::CalculateGradient)
        .def("GetValue", &EigenfrequencyResponseFunction::GetValue)
        .def("GetInitialValue", &EigenfrequencyResponseFunction::GetInitialValue)
        .def("GetGradient", &EigenfrequencyResponseFunction::GetGradient)
        ;

    // class_<EigenfrequencyResponseFunctionLinScal, bases<Process> >("EigenfrequencyResponseFunctionLinScal", init<ModelPart&, Parameters&>())
    //     .def("initialize", &EigenfrequencyResponseFunctionLinScal::initialize)
    //     .def("calculate_value", &EigenfrequencyResponseFunctionLinScal::calculate_value)
    //     .def("calculate_gradient", &EigenfrequencyResponseFunctionLinScal::calculate_gradient)
    //     .def("get_value", &EigenfrequencyResponseFunctionLinScal::get_value)
    //     .def("get_initial_value", &EigenfrequencyResponseFunctionLinScal::get_initial_value)
    //     .def("get_gradient", &EigenfrequencyResponseFunctionLinScal::get_gradient)
    //     ;

    // class_<EigenfrequencyResponseFunctionKS, bases<Process> >("EigenfrequencyResponseFunctionKS", init<ModelPart&, Parameters&>())
    //     .def("initialize", &EigenfrequencyResponseFunctionKS::initialize)
    //     .def("calculate_value", &EigenfrequencyResponseFunctionKS::calculate_value)
    //     .def("calculate_gradient", &EigenfrequencyResponseFunctionKS::calculate_gradient)
    //     .def("get_value", &EigenfrequencyResponseFunctionKS::get_value)
    //     .def("get_initial_value", &EigenfrequencyResponseFunctionKS::get_initial_value)
    //     .def("get_gradient", &EigenfrequencyResponseFunctionKS::get_gradient)
    //     ;

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

