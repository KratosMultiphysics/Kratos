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
#include "custom_utilities/mapping/mapper_vertex_morphing_symmetric.h"
#include "custom_utilities/damping/damping_utilities.h"
#include "custom_utilities/mesh_controller_utilities.h"
#include "custom_utilities/input_output/universal_file_io.h"
#include "custom_utilities/search_based_functions.h"
#include "custom_utilities/response_functions/face_angle_response_function_utility.h"

// ==============================================================================

namespace Kratos {
namespace Python {

// Overloaded functions
template<typename TMapper>
inline void MapVector(TMapper& mapper,
                const Variable< array_1d<double, 3> >& origin_variable,
                const Variable< array_1d<double, 3> >& destination_variable)
{
    mapper.Map(origin_variable, destination_variable);
}

template<typename TMapper>
inline void MapScalar(TMapper& mapper,
                const Variable< double >& origin_variable,
                const Variable< double >& destination_variable)
{
    mapper.Map(origin_variable, destination_variable);
}

template<typename TMapper>
inline void InverseMapVector(TMapper& mapper,
                       const Variable< array_1d<double, 3> >& origin_variable,
                       const Variable< array_1d<double, 3> >& destination_variable)
{
    mapper.InverseMap(origin_variable, destination_variable);
}
template<typename TMapper>
inline void InverseMapScalar(TMapper& mapper,
                       const Variable< double >& origin_variable,
                       const Variable< double >& destination_variable)
{
    mapper.InverseMap(origin_variable, destination_variable);
}

inline void AssembleMatrixForVariableList(
    ModelPart& rModelPart,
    Matrix& rMatrix,
    pybind11::list& rVariables)
{
    std::size_t list_length = pybind11::len(rVariables);
    std::vector<Variable<OptimizationUtilities::array_3d>*> variables_vector(list_length);
    for (std::size_t i = 0; i < list_length; i++)
    {
        variables_vector[i] = (rVariables[i]).cast<Variable<OptimizationUtilities::array_3d>*>();
    }
    return OptimizationUtilities::AssembleMatrix(rModelPart, rMatrix, variables_vector);
}

// ==============================================================================
void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // ================================================================
    // For perfoming the mapping according to Vertex Morphing
    // ================================================================
    py::class_<MapperVertexMorphing >(m, "MapperVertexMorphing")
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphing::Initialize)
        .def("Update", &MapperVertexMorphing::Update)
        .def("Map", MapScalar<MapperVertexMorphing>)
        .def("Map", MapVector<MapperVertexMorphing>)
        .def("InverseMap", InverseMapScalar<MapperVertexMorphing>)
        .def("InverseMap", InverseMapVector<MapperVertexMorphing>)
        ;
    py::class_<MapperVertexMorphingMatrixFree >(m, "MapperVertexMorphingMatrixFree")
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphingMatrixFree::Initialize)
        .def("Update", &MapperVertexMorphingMatrixFree::Update)
        .def("Map", MapScalar<MapperVertexMorphingMatrixFree>)
        .def("Map", MapVector<MapperVertexMorphingMatrixFree>)
        .def("InverseMap", InverseMapScalar<MapperVertexMorphingMatrixFree>)
        .def("InverseMap", InverseMapVector<MapperVertexMorphingMatrixFree>)
        ;
    py::class_<MapperVertexMorphingImprovedIntegration >(m, "MapperVertexMorphingImprovedIntegration")
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphingImprovedIntegration::Initialize)
        .def("Update", &MapperVertexMorphingImprovedIntegration::Update)
        .def("Map", MapScalar<MapperVertexMorphingImprovedIntegration>)
        .def("Map", MapVector<MapperVertexMorphingImprovedIntegration>)
        .def("InverseMap", InverseMapScalar<MapperVertexMorphingImprovedIntegration>)
        .def("InverseMap", InverseMapVector<MapperVertexMorphingImprovedIntegration>)
        ;
    py::class_<MapperVertexMorphingSymmetric >(m, "MapperVertexMorphingSymmetric")
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphingSymmetric::Initialize)
        .def("Update", &MapperVertexMorphingSymmetric::Update)
        .def("Map", MapScalar<MapperVertexMorphingSymmetric>) // TODO
        .def("Map", MapVector<MapperVertexMorphingSymmetric>)
        .def("InverseMap", InverseMapScalar<MapperVertexMorphingSymmetric>) // TODO
        .def("InverseMap", InverseMapVector<MapperVertexMorphingSymmetric>)
        ;

    // ================================================================
    // For a possible damping of nodal variables
    // ================================================================
    py::class_<DampingUtilities >(m, "DampingUtilities")
        .def(py::init<ModelPart&, pybind11::dict, Parameters>())
        .def("DampNodalVariable", &DampingUtilities::DampNodalVariable)
        ;

    // ========================================================================
    // For performing individual steps of an optimization algorithm
    // ========================================================================
    py::class_<OptimizationUtilities >(m, "OptimizationUtilities")
        // ----------------------------------------------------------------
        // For running unconstrained descent methods
        // ----------------------------------------------------------------
        .def_static("ComputeSearchDirectionSteepestDescent", &OptimizationUtilities::ComputeSearchDirectionSteepestDescent)
        // ----------------------------------------------------------------
        // For running penalized projection method
        // ----------------------------------------------------------------
        .def_static("ComputeProjectedSearchDirection", &OptimizationUtilities::ComputeProjectedSearchDirection)
        .def_static("CorrectProjectedSearchDirection", &OptimizationUtilities::CorrectProjectedSearchDirection)
        // ----------------------------------------------------------------
        // General optimization operations
        // ----------------------------------------------------------------
        .def_static("ComputeControlPointUpdate", &OptimizationUtilities::ComputeControlPointUpdate)
        .def_static("AddFirstVariableToSecondVariable", &OptimizationUtilities::AddFirstVariableToSecondVariable)
        .def_static("ComputeL2NormOfNodalVariable", [](ModelPart& rModelPart, const Variable< double >& rVariable){
                                                        return OptimizationUtilities::ComputeL2NormOfNodalVariable(rModelPart, rVariable);
                                                    })
        .def_static("ComputeL2NormOfNodalVariable", [](ModelPart& rModelPart, const Variable< array_1d<double, 3> >& rVariable){
                                                        return OptimizationUtilities::ComputeL2NormOfNodalVariable(rModelPart, rVariable);
                                                    })
        .def_static("ComputeMaxNormOfNodalVariable", [](ModelPart& rModelPart, const Variable< double >& rVariable){
                                                        return OptimizationUtilities::ComputeMaxNormOfNodalVariable(rModelPart, rVariable);
                                                        })
        .def_static("ComputeMaxNormOfNodalVariable", [](ModelPart& rModelPart, const Variable< array_1d<double, 3> >& rVariable){
                                                        return OptimizationUtilities::ComputeMaxNormOfNodalVariable(rModelPart, rVariable);
                                                        })
        .def_static("AssembleVector", &OptimizationUtilities::AssembleVector)
        .def_static("AssignVectorToVariable", &OptimizationUtilities::AssignVectorToVariable)
        .def_static("AssembleMatrix", [](ModelPart& rModelPart, Matrix& rMatrix, pybind11::list& rVariables){
                                            std::size_t list_length = pybind11::len(rVariables);
                                            std::vector<Variable<OptimizationUtilities::array_3d>*> variables_vector(list_length);
                                            for (std::size_t i = 0; i < list_length; i++)
                                            {
                                                variables_vector[i] = (rVariables[i]).cast<Variable<OptimizationUtilities::array_3d>*>();
                                            }
                                            return OptimizationUtilities::AssembleMatrix(rModelPart, rMatrix, variables_vector);
                                        })
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &OptimizationUtilities::CalculateProjectedSearchDirectionAndCorrection)
        ;

    // ========================================================================
    // For pre- and post-processing of geometry data
    // ========================================================================
    py::class_<GeometryUtilities >(m, "GeometryUtilities")
        .def(py::init<ModelPart&>())
        .def("ComputeUnitSurfaceNormals", &GeometryUtilities::ComputeUnitSurfaceNormals)
        .def("ProjectNodalVariableOnUnitSurfaceNormals", &GeometryUtilities::ProjectNodalVariableOnUnitSurfaceNormals)
        .def("ProjectNodalVariableOnDirection", &GeometryUtilities::ProjectNodalVariableOnDirection)
        .def("ProjectNodalVariableOnTangentPlane", &GeometryUtilities::ProjectNodalVariableOnTangentPlane)
        .def("ExtractBoundaryNodes", &GeometryUtilities::ExtractBoundaryNodes)
        .def("ComputeDistancesToBoundingModelPart", &GeometryUtilities::ComputeDistancesToBoundingModelPart)
        .def("CalculateLength",&GeometryUtilities::CalculateLength<ModelPart::ElementsContainerType>)
        .def("CalculateLength",&GeometryUtilities::CalculateLength<ModelPart::ConditionsContainerType>)
        ;

    // ========================================================================
    // For mesh handling
    // ========================================================================
    py::class_<MeshControllerUtilities >(m, "MeshControllerUtilities")
        .def(py::init<ModelPart&>())
        .def("UpdateMeshAccordingInputVariable", &MeshControllerUtilities::UpdateMeshAccordingInputVariable)
        .def("RevertMeshUpdateAccordingInputVariable", &MeshControllerUtilities::RevertMeshUpdateAccordingInputVariable)
        .def("LogMeshChangeAccordingInputVariable", &MeshControllerUtilities::LogMeshChangeAccordingInputVariable)
        .def("SetMeshToReferenceMesh", &MeshControllerUtilities::SetMeshToReferenceMesh)
        .def("SetReferenceMeshToMesh", &MeshControllerUtilities::SetReferenceMeshToMesh)
        .def("SetDeformationVariablesToZero", &MeshControllerUtilities::SetDeformationVariablesToZero)
        .def("WriteCoordinatesToVariable", &MeshControllerUtilities::WriteCoordinatesToVariable)
        .def("SubtractCoordinatesFromVariable", &MeshControllerUtilities::SubtractCoordinatesFromVariable)
        .def("AddFirstVariableToSecondVariable", &MeshControllerUtilities::AddFirstVariableToSecondVariable)
        ;

    // ========================================================================
    // For input / output
    // ========================================================================
    py::class_<UniversalFileIO >(m, "UniversalFileIO")
        .def(py::init<ModelPart&, std::string, std::string, Parameters>())
        .def("InitializeLogging", &UniversalFileIO::InitializeLogging)
        .def("LogNodalResults", &UniversalFileIO::LogNodalResults)
        ;

    // ========================================================================
    // For geometric response functions
    // ========================================================================
    py::class_<FaceAngleResponseFunctionUtility >(m, "FaceAngleResponseFunctionUtility")
        .def(py::init<ModelPart&, Parameters>())
        .def("Initialize", &FaceAngleResponseFunctionUtility::Initialize)
        .def("CalculateValue", &FaceAngleResponseFunctionUtility::CalculateValue)
        .def("CalculateGradient", &FaceAngleResponseFunctionUtility::CalculateGradient)
        ;

    // ========================================================================
    // Additional operations
    // ========================================================================
    py::class_<SearchBasedFunctions >(m, "SearchBasedFunctions")
        .def(py::init<ModelPart&>())
        .def("FlagNodesInRadius", &SearchBasedFunctions::FlagNodesInRadius)
        ;

}

}  // namespace Python.
} // Namespace Kratos

