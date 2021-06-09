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
#include "custom_utilities/mapping/mapper_vertex_morphing_mesh_independent.h"
#include "custom_utilities/mapping/mapper_vertex_morphing_normal.h"
#include "custom_utilities/mapping/mapper_vertex_morphing_mesh_independent_normal.h"
#include "custom_utilities/mapping/mapper_vertex_morphing_symmetric.h"
#include "custom_utilities/damping/damping_utilities.h"
#include "custom_utilities/mesh_controller_utilities.h"
#include "custom_utilities/input_output/universal_file_io.h"
#include "custom_utilities/search_based_functions.h"
#include "custom_utilities/lumped_integration_utility.h"
#include "custom_utilities/response_functions/surface_area_response_function_utility.h"
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
inline void MapMixed(TMapper& mapper,
                const Variable< double >& origin_variable,
                const Variable< array_1d<double, 3> >& destination_variable)
{
    mapper.Map(origin_variable, destination_variable);
}


template<typename TMapper>
inline void InverseMapVector(TMapper& mapper,
                       const Variable< array_1d<double, 3> >& destination_variable,
                       const Variable< array_1d<double, 3> >& origin_variable)
{
    mapper.InverseMap(destination_variable, origin_variable);
}
template<typename TMapper>
inline void InverseMapScalar(TMapper& mapper,
                       const Variable< double >& destination_variable,
                       const Variable< double >& origin_variable)
{
    mapper.InverseMap(destination_variable, origin_variable);
}
template<typename TMapper>
inline void InverseMapMixed(TMapper& mapper,
                       const Variable< array_1d<double, 3> >& destination_variable,
                       const Variable< double >& origin_variable)
{
    mapper.InverseMap(destination_variable, origin_variable);
}

inline void ComputeSearchDirectionSteepestDescentScalar(OptimizationUtilities& utils, const Variable<double>& rSearchDirection, const Variable<double>& rObjectiveGradient)
{
    return utils.ComputeSearchDirectionSteepestDescent(rSearchDirection, rObjectiveGradient);
}

inline void ComputeSearchDirectionSteepestDescentVector(OptimizationUtilities& utils, const Variable<array_1d<double, 3>>& rSearchDirection, const Variable<array_1d<double, 3>>& rObjectiveGradient)
{
    return utils.ComputeSearchDirectionSteepestDescent(rSearchDirection, rObjectiveGradient);
}

inline void ComputeProjectedSearchDirectionScalar(OptimizationUtilities& utils, const Variable<double>& rSearchDirection, const Variable<double>& rObjectiveGradient, const Variable<double>& rConstraintGradient)
{
    return utils.ComputeProjectedSearchDirection(rSearchDirection, rObjectiveGradient, rConstraintGradient);
}

inline void ComputeProjectedSearchDirectionVector(OptimizationUtilities& utils, const Variable<array_1d<double, 3>>& rSearchDirection, const Variable<array_1d<double, 3>>& rObjectiveGradient, const Variable<array_1d<double, 3>>& rConstraintGradient)
{
    return utils.ComputeProjectedSearchDirection(rSearchDirection, rObjectiveGradient, rConstraintGradient);
}

inline void CorrectProjectedSearchDirectionScalar(OptimizationUtilities& utils, double ConstraintValue, const Variable<double>& rSearchDirection, const Variable<double>& rConstraintGradient)
{
    return utils.CorrectProjectedSearchDirection(ConstraintValue, rSearchDirection, rConstraintGradient);
}

inline void CorrectProjectedSearchDirectionVector(OptimizationUtilities& utils, double ConstraintValue, const Variable<array_1d<double, 3>>& rSearchDirection, const Variable<array_1d<double, 3>>& rConstraintGradient)
{
    return utils.CorrectProjectedSearchDirection(ConstraintValue, rSearchDirection, rConstraintGradient);
}

inline void ComputeControlPointUpdateScalar(OptimizationUtilities& utils, const double StepSize, const Variable<double>& rSearchDirection, const Variable<double>& rControlUpdate)
{
    return utils.ComputeControlPointUpdate(StepSize, rSearchDirection, rControlUpdate);
}

inline void ComputeControlPointUpdateVector(OptimizationUtilities& utils, const double StepSize, const Variable<array_1d<double, 3>>& rSearchDirection, const Variable<array_1d<double, 3>>& rControlUpdate)
{
    return utils.ComputeControlPointUpdate(StepSize, rSearchDirection, rControlUpdate);
}

inline void AddFirstVariableToSecondVariableScalar(OptimizationUtilities& utils, const Variable<double>& rVariable1, const Variable<double>& rVariable2)
{
    return utils.AddFirstVariableToSecondVariable(rVariable1, rVariable2);
}

inline void AddFirstVariableToSecondVariableVector(OptimizationUtilities& utils, const Variable<array_1d<double, 3>>& rVariable1, const Variable<array_1d<double, 3>>& rVariable2)
{
    return utils.AddFirstVariableToSecondVariable(rVariable1, rVariable2);
}

inline double ComputeL2NormScalar(OptimizationUtilities& utils, const Variable< double >& variable)
{
    return utils.ComputeL2NormOfNodalVariable(variable);
}

inline double ComputeL2NormVector(OptimizationUtilities& utils, const Variable< array_1d<double, 3> >& variable)
{
    return utils.ComputeL2NormOfNodalVariable(variable);
}

inline double ComputeMaxNormScalar(OptimizationUtilities& utils, const Variable< double >& variable)
{
    return utils.ComputeMaxNormOfNodalVariable(variable);
}

inline double ComputeMaxNormVector(OptimizationUtilities& utils, const Variable< array_1d<double, 3> >& variable)
{
    return utils.ComputeMaxNormOfNodalVariable(variable);
}

inline void AssembleMatrixForScalarVariableList(
    OptimizationUtilities& utils,
    Matrix& rMatrix,
    pybind11::list& rVariables)
{
    std::size_t list_length = pybind11::len(rVariables);
    std::vector<Variable<double>*> variables_vector(list_length);
    for (std::size_t i = 0; i < list_length; i++)
    {
        variables_vector[i] = (rVariables[i]).cast<Variable<double>*>();
    }
    return utils.AssembleMatrix(rMatrix, variables_vector);
}

inline void AssembleMatrixForVectorVariableList(
    OptimizationUtilities& utils,
    Matrix& rMatrix,
    pybind11::list& rVariables)
{
    std::size_t list_length = pybind11::len(rVariables);
    std::vector<Variable<OptimizationUtilities::array_3d>*> variables_vector(list_length);
    for (std::size_t i = 0; i < list_length; i++)
    {
        variables_vector[i] = (rVariables[i]).cast<Variable<OptimizationUtilities::array_3d>*>();
    }
    return utils.AssembleMatrix(rMatrix, variables_vector);
}

inline void IntegrateScalarVariable(LumpedIntegrationUtility& util, const Variable< double >& variable)
{
    util.Integrate(variable);
}

inline void IntegrateVectorVariable(LumpedIntegrationUtility& util, const Variable< array_1d<double, 3> >& variable)
{
    util.Integrate(variable);
}

inline void AssembleScalarToVector(OptimizationUtilities& util, Vector& rVector, const Variable<double> &rVariable)
{
    util.AssembleVector(rVector, rVariable);
}
inline void AssembleVectorToVector(OptimizationUtilities& util, Vector& rVector, const Variable<array_1d<double, 3>> &rVariable)
{
    util.AssembleVector(rVector, rVariable);
}

inline void AssignVectorToScalarVariable(OptimizationUtilities& util, const Vector& rVector, const Variable<double> &rVariable)
{
    util.AssignVectorToVariable(rVector, rVariable);
}
inline void AssignVectorToVectorVariable(OptimizationUtilities& util, const Vector& rVector, const Variable<array_1d<double, 3>> &rVariable)
{
    util.AssignVectorToVariable(rVector, rVariable);
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
    py::class_<MapperVertexMorphingMeshIndependent >(m, "MapperVertexMorphingMeshIndependent")
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphingMeshIndependent::Initialize)
        .def("Update", &MapperVertexMorphingMeshIndependent::Update)
        .def("Map", MapScalar<MapperVertexMorphingMeshIndependent>)
        .def("Map", MapVector<MapperVertexMorphingMeshIndependent>)
        .def("InverseMap", InverseMapScalar<MapperVertexMorphingMeshIndependent>)
        .def("InverseMap", InverseMapVector<MapperVertexMorphingMeshIndependent>)
        ;
    py::class_<MapperVertexMorphingNormal >(m, "MapperVertexMorphingNormal")
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphingNormal::Initialize)
        .def("Update", &MapperVertexMorphingNormal::Update)
        .def("Map", MapMixed<MapperVertexMorphingNormal>)
        .def("InverseMap", InverseMapMixed<MapperVertexMorphingNormal>)
        ;
    py::class_<MapperVertexMorphingMeshIndependentNormal >(m, "MapperVertexMorphingMeshIndependentNormal")
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("Initialize", &MapperVertexMorphingMeshIndependentNormal::Initialize)
        .def("Update", &MapperVertexMorphingMeshIndependentNormal::Update)
        .def("Map", MapMixed<MapperVertexMorphingMeshIndependentNormal>)
        .def("InverseMap", InverseMapMixed<MapperVertexMorphingMeshIndependentNormal>)
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
        .def(py::init<ModelPart&, Parameters>())
        // ----------------------------------------------------------------
        // For running unconstrained descent methods
        // ----------------------------------------------------------------
        .def("ComputeSearchDirectionSteepestDescent", &ComputeSearchDirectionSteepestDescentScalar)
        .def("ComputeSearchDirectionSteepestDescent", &ComputeSearchDirectionSteepestDescentVector)
        // ----------------------------------------------------------------
        // For running penalized projection method
        // ----------------------------------------------------------------
        .def("ComputeProjectedSearchDirection", &ComputeProjectedSearchDirectionScalar)
        .def("ComputeProjectedSearchDirection", &ComputeProjectedSearchDirectionVector)
        .def("CorrectProjectedSearchDirection", &CorrectProjectedSearchDirectionScalar)
        .def("CorrectProjectedSearchDirection", &CorrectProjectedSearchDirectionVector)
        .def("GetCorrectionScaling", &OptimizationUtilities::GetCorrectionScaling)
        // ----------------------------------------------------------------
        // General optimization operations
        // ----------------------------------------------------------------
        .def("ComputeControlPointUpdate", &ComputeControlPointUpdateScalar)
        .def("ComputeControlPointUpdate", &ComputeControlPointUpdateVector)
        .def("AddFirstVariableToSecondVariable", &AddFirstVariableToSecondVariableScalar)
        .def("AddFirstVariableToSecondVariable", &AddFirstVariableToSecondVariableVector)
        .def("ComputeL2NormOfNodalVariable", ComputeL2NormScalar)
        .def("ComputeL2NormOfNodalVariable", ComputeL2NormVector)
        .def("ComputeMaxNormOfNodalVariable", ComputeMaxNormScalar)
        .def("ComputeMaxNormOfNodalVariable", ComputeMaxNormVector)
        .def("AssembleVector", &AssembleVectorToVector)
        .def("AssignVectorToVariable", &AssignVectorToVectorVariable)
        .def("AssembleMatrixForVector", &AssembleMatrixForVectorVariableList)
        .def("AssembleVector", &AssembleScalarToVector)
        .def("AssignVectorToVariable", &AssignVectorToScalarVariable)
        .def("AssembleMatrixForScalar", &AssembleMatrixForScalarVariableList)
        .def("CalculateProjectedSearchDirectionAndCorrection", &OptimizationUtilities::CalculateProjectedSearchDirectionAndCorrection)
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

    py::class_<LumpedIntegrationUtility>(m, "LumpedIntegrationUtility")
        .def(py::init<ModelPart&>())
        .def("CalculateLumpedAreas", &LumpedIntegrationUtility::CalculateLumpedAreas)
        .def("Integrate", IntegrateVectorVariable)
        .def("Integrate", IntegrateScalarVariable)
        ;

    py::class_<SurfaceAreaResponseFunctionUtility>(m, "SurfaceAreaResponseFunctionUtility")
        .def(py::init<ModelPart&, Parameters>())
        .def("CalculateValue", &SurfaceAreaResponseFunctionUtility::CalculateValue)
        .def("CalculateGradient", &SurfaceAreaResponseFunctionUtility::CalculateGradient)
        ;
}

}  // namespace Python.
} // Namespace Kratos

