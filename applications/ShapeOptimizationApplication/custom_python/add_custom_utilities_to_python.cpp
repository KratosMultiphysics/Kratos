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
#include <pybind11/stl.h>

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
#include "custom_utilities/mapping/mapper_vertex_morphing_adaptive_radius.h"
#include "custom_utilities/damping/damping_utilities.h"
#include "custom_utilities/damping/direction_damping_utilities.h"
#include "custom_utilities/damping/thickness_damping_utilities.h"
#include "custom_utilities/mesh_controller_utilities.h"
#include "custom_utilities/input_output/universal_file_io.h"
#include "custom_utilities/search_based_functions.h"
#include "custom_utilities/response_functions/face_angle_response_function_utility.h"
#include "custom_utilities/response_functions/water_drain_response_function_utility.h"
#include "custom_utilities/response_functions/directional_derivative_response_function_utility.h"

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
    py::class_<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphing> >(m, "MapperVertexMorphingAdaptiveRadius")
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphingAdaptiveRadius<MapperVertexMorphing>::Initialize)
        .def("Update", &MapperVertexMorphingAdaptiveRadius<MapperVertexMorphing>::Update)
        .def("Map", MapScalar<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphing>>)
        .def("Map", MapVector<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphing>>)
        .def("InverseMap", InverseMapScalar<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphing>>)
        .def("InverseMap", InverseMapVector<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphing>>)
        ;
    py::class_<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingMatrixFree> >(m, "MapperVertexMorphingMatrixFreeAdaptiveRadius")
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingMatrixFree>::Initialize)
        .def("Update", &MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingMatrixFree>::Update)
        .def("Map", MapScalar<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingMatrixFree>>)
        .def("Map", MapVector<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingMatrixFree>>)
        .def("InverseMap", InverseMapScalar<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingMatrixFree>>)
        .def("InverseMap", InverseMapVector<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingMatrixFree>>)
        ;
    py::class_<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingImprovedIntegration> >(m, "MapperVertexMorphingImprovedIntegrationAdaptiveRadius")
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingImprovedIntegration>::Initialize)
        .def("Update", &MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingImprovedIntegration>::Update)
        .def("Map", MapScalar<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingImprovedIntegration>>)
        .def("Map", MapVector<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingImprovedIntegration>>)
        .def("InverseMap", InverseMapScalar<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingImprovedIntegration>>)
        .def("InverseMap", InverseMapVector<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingImprovedIntegration>>)
        ;
    py::class_<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingSymmetric> >(m, "MapperVertexMorphingSymmetricAdaptiveRadius")
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingSymmetric>::Initialize)
        .def("Update", &MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingSymmetric>::Update)
        .def("Map", MapScalar<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingSymmetric>>)
        .def("Map", MapVector<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingSymmetric>>)
        .def("InverseMap", InverseMapScalar<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingSymmetric>>)
        .def("InverseMap", InverseMapVector<MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingSymmetric>>)
        ;

    // ================================================================
    // For a possible damping of nodal variables
    // ================================================================
    py::class_<DampingUtilities >(m, "DampingUtilities")
        .def(py::init<ModelPart&, Parameters>())
        .def("DampNodalVariable", &DampingUtilities::DampNodalVariable)
        ;

    py::class_<DirectionDampingUtilities >(m, "DirectionDampingUtilities")
        .def(py::init<ModelPart&, Parameters>())
        .def("DampNodalVariable", &DirectionDampingUtilities::DampNodalVariable)
        ;

    py::class_<ThicknessDampingUtilities >(m, "ThicknessDampingUtilities")
        .def(py::init<ModelPart&, Parameters>())
        .def("DampNodalVariable", &ThicknessDampingUtilities::DampNodalVariable)
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
        .def_static("AssembleVector", [](ModelPart& rModelPart, Vector& rVector, const Variable< double >& rVariable){
                                         OptimizationUtilities::AssembleVector(rModelPart, rVector, rVariable);
                                         })
        .def_static("AssembleVector", [](ModelPart& rModelPart, Vector& rVector, const Variable< array_1d<double, 3> >& rVariable){
                                         OptimizationUtilities::AssembleVector(rModelPart, rVector, rVariable);
                                         })
        .def_static("AssignVectorToVariable", [](ModelPart& rModelPart, Vector& rVector, const Variable< double >& rVariable) {
                                                 OptimizationUtilities::AssignVectorToVariable(rModelPart, rVector, rVariable);
                                                 })
        .def_static("AssignVectorToVariable", [](ModelPart& rModelPart, Vector& rVector, const Variable< array_1d<double, 3> >& rVariable){
                                                 OptimizationUtilities::AssignVectorToVariable(rModelPart, rVector, rVariable);
                                                 })
        .def_static("AssembleMatrix", [](ModelPart& rModelPart, Matrix& rMatrix, pybind11::list& rVariables){
                                            std::size_t list_length = pybind11::len(rVariables);
                                            std::vector<Variable<OptimizationUtilities::array_3d>*> variables_vector(list_length);
                                            for (std::size_t i = 0; i < list_length; i++)
                                            {
                                                variables_vector[i] = (rVariables[i]).cast<Variable<OptimizationUtilities::array_3d>*>();
                                            }
                                            return OptimizationUtilities::AssembleMatrix(rModelPart, rMatrix, variables_vector);
                                        })
        .def_static("AssembleMatrixScalarVariables", [](ModelPart& rModelPart, Matrix& rMatrix, pybind11::list& rVariables){
                                    std::size_t list_length = pybind11::len(rVariables);
                                    std::vector<Variable<double>*> variables_vector(list_length);
                                    for (std::size_t i = 0; i < list_length; i++)
                                    {
                                        variables_vector[i] = (rVariables[i]).cast<Variable<double>*>();
                                    }
                                    return OptimizationUtilities::AssembleMatrix(rModelPart, rMatrix, variables_vector);
                                })
        .def_static("AssembleMatrixFromGradientVectors", [](ModelPart& rModelPart, Matrix& rMatrix, pybind11::list& rGradientVectors){
                                            std::size_t list_length = pybind11::len(rGradientVectors);
                                            std::vector<Vector*> gradient_vectors(list_length);
                                            for (std::size_t i = 0; i < list_length; i++)
                                            {
                                                gradient_vectors[i] = (rGradientVectors[i]).cast<Vector*>();
                                            }
                                            return OptimizationUtilities::AssembleMatrix(rModelPart, rMatrix, gradient_vectors);
                                        })
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &OptimizationUtilities::CalculateProjectedSearchDirectionAndCorrection)
        .def_static("AssembleBufferMatrix", &OptimizationUtilities::AssembleBufferMatrix)
        .def_static("CalculateRelaxedProjectedSearchDirectionAndCorrection", &OptimizationUtilities::CalculateRelaxedProjectedSearchDirectionAndCorrection)
        ;

    // ========================================================================
    // For pre- and post-processing of geometry data
    // ========================================================================
    py::class_<GeometryUtilities >(m, "GeometryUtilities")
        .def(py::init<ModelPart&>())
        .def("CalculateNodalAreasFromConditions", &GeometryUtilities::CalculateNodalAreasFromConditions)
        .def("ComputeUnitSurfaceNormals", &GeometryUtilities::ComputeUnitSurfaceNormals)
        .def("ProjectNodalVariableOnUnitSurfaceNormals", &GeometryUtilities::ProjectNodalVariableOnUnitSurfaceNormals)
        .def("ProjectNodalVariableOnDirection", &GeometryUtilities::ProjectNodalVariableOnDirection)
        .def("ProjectNodalVariableOnTangentPlane", &GeometryUtilities::ProjectNodalVariableOnTangentPlane)
        .def("ExtractBoundaryNodes", &GeometryUtilities::ExtractBoundaryNodes)
        .def("ExtractEdgeNodes", &GeometryUtilities::ExtractEdgeNodes)
        .def("ComputeDistancesToBoundingModelPart", &GeometryUtilities::ComputeDistancesToBoundingModelPart)
        .def("CalculateLength",&GeometryUtilities::CalculateLength<ModelPart::ElementsContainerType>)
        .def("CalculateLength",&GeometryUtilities::CalculateLength<ModelPart::ConditionsContainerType>)
        .def("ComputeVolume", &GeometryUtilities::ComputeVolume)
        .def("ComputeVolumeShapeDerivatives", &GeometryUtilities::ComputeVolumeShapeDerivatives)
        .def("CalculateAverageElementSize", &GeometryUtilities::CalculateAverageElementSize)
        ;

    // ========================================================================
    // For mesh handling
    // ========================================================================
    py::class_<MeshControllerUtilities >(m, "MeshControllerUtilities")
        .def(py::init<ModelPart&>())
        .def("UpdateMeshAccordingInputVariable", &MeshControllerUtilities::UpdateMeshAccordingInputVariable)
        .def("UpdateThicknessAccordingInputVariable", &MeshControllerUtilities::UpdateThicknessAccordingInputVariable)
        .def("UpdateThicknessAccordingInitialAndInputVariable", &MeshControllerUtilities::UpdateThicknessAccordingInitialAndInputVariable)
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

    py::class_<WaterDrainResponseFunctionUtility >(m, "WaterDrainResponseFunctionUtility")
        .def(py::init<ModelPart&, Parameters>())
        .def("Initialize", &WaterDrainResponseFunctionUtility::Initialize)
        .def("InitializeSolutionStep", &WaterDrainResponseFunctionUtility::InitializeSolutionStep)
        .def("CalculateValue", &WaterDrainResponseFunctionUtility::CalculateValue)
        .def("CalculateGradient", &WaterDrainResponseFunctionUtility::CalculateGradient)

    py::class_<DirectionalDerivativeResponseFunctionUtility >(m, "DirectionalDerivativeResponseFunctionUtility")
        .def(py::init<ModelPart&, Parameters>())
        .def("Initialize", &DirectionalDerivativeResponseFunctionUtility::Initialize)
        .def("CalculateValue", &DirectionalDerivativeResponseFunctionUtility::CalculateValue)
        .def("CalculateGradient", &DirectionalDerivativeResponseFunctionUtility::CalculateGradient)
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

