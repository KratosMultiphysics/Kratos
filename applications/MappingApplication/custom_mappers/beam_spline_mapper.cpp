//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
// Main authors:
// Contributor:


// System includes
#include <algorithm>
#include <array>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

// External includes
#include "utilities/math_utils.h"
#include "geometries/line_3d_2.h"

// Project includes
#include "beam_spline_mapper.h"
#include "mappers/mapper_define.h"
#include "mapping_application_variables.h"

namespace Kratos
{

void BeamSplineMapperInterfaceInfo::ProcessSearchResult(const InterfaceObject& rInterfaceObject)
{
    SaveSearchResult(rInterfaceObject, true);
}

void BeamSplineMapperInterfaceInfo::ProcessSearchResultForApproximation(const InterfaceObject& rInterfaceObject)
{
    const auto p_geom = rInterfaceObject.pGetBaseGeometry();

    double proj_dist;

    const Point point_to_proj(this->Coordinates());
    Point projected_point;

    mCoordinates = point_to_proj;

    Vector linear_shape_function_values;
    std::vector<int> eq_ids;

    for (auto& r_node : p_geom->Points()) {
        r_node.X() = r_node.X0();
        r_node.Y() = r_node.Y0();
        r_node.Z() = r_node.Z0();
    }

    // First project the surface point onto the infinite line of this beam candidate.
    std::ignore = GeometricalProjectionUtilities::FastProjectOnLine((*p_geom), point_to_proj, projected_point);

    const double distance_to_node_1 = norm_2(projected_point - (*p_geom)[0]);
    const double distance_to_node_2 = norm_2(projected_point - (*p_geom)[1]);
    const double proj_dist_nodes = std::min(distance_to_node_1, distance_to_node_2);

    if (proj_dist_nodes < mClosestProjectionDistance) {
        SetIsApproximation();

        mPairingIndex = ProjectionUtilities::PairingIndex::Closest_Point;
        p_geom->IsInside(projected_point, local_coords, 1e-14);
        p_geom->PointLocalCoordinates(local_coords, projected_point);

        // Approximate projections can fall outside the beam segment; clamp them to the valid line element range.
        local_coords[0] = std::clamp(local_coords[0], -1.0, 1.0);

        const bool compute_approximation = false;
        std::ignore = ProjectionUtilities::ProjectOnLine(
            *p_geom,
            point_to_proj,
            mLocalCoordTol,
            linear_shape_function_values,
            eq_ids,
            proj_dist,
            compute_approximation);

        if (eq_ids.size() != p_geom->PointsNumber()) {
            eq_ids.resize(p_geom->PointsNumber());
            for (IndexType i = 0; i < p_geom->PointsNumber(); ++i) {
                KRATOS_DEBUG_ERROR_IF_NOT((*p_geom)[i].Has(INTERFACE_EQUATION_ID))
                    << (*p_geom)[i] << " does not have an \"INTERFACE_EQUATION_ID\"" << std::endl;
                eq_ids[i] = (*p_geom)[i].GetValue(INTERFACE_EQUATION_ID);
            }
        }

        // Recompute the interpolation weights from the clamped local coordinate used by the approximation.
        p_geom->ShapeFunctionsValues(linear_shape_function_values, local_coords);

        mClosestProjectionDistance = proj_dist_nodes;
        mNodeIds = eq_ids;

        const IndexType num_values_linear = linear_shape_function_values.size();
        if (mLinearShapeFunctionValues.size() != num_values_linear) {
            mLinearShapeFunctionValues.resize(num_values_linear);
        }

        for (IndexType i = 0; i < num_values_linear; ++i) {
            mLinearShapeFunctionValues[i] = linear_shape_function_values[i];
        }

        mProjectionOfPoint = projected_point;
        mpInterfaceObject = make_shared<InterfaceGeometryObject>(rInterfaceObject.pGetBaseGeometry());
    }

}

void BeamSplineMapperInterfaceInfo::SaveSearchResult(
    const InterfaceObject& rInterfaceObject,
    const bool ComputeApproximation)
{
    const auto p_geom = rInterfaceObject.pGetBaseGeometry();

    double proj_dist;

    const Point point_to_proj(this->Coordinates());
    Point projected_point;

    mCoordinates = point_to_proj;

    Vector linear_shape_function_values;
    std::vector<int> eq_ids;

    for (auto& r_node : p_geom->Points()) {
        r_node.X() = r_node.X0();
        r_node.Y() = r_node.Y0();
        r_node.Z() = r_node.Z0();
    }

    const auto geom_family = p_geom->GetGeometryFamily();
    KRATOS_ERROR_IF(geom_family != GeometryData::KratosGeometryFamily::Kratos_Linear) << "Invalid geometry of the Origin! The geometry should be a beam!";

    // Compute the line projection weights and classify how good this candidate beam is.
    const auto pairing_index_linear = ProjectionUtilities::ProjectOnLine(
        *p_geom,
        point_to_proj,
        mLocalCoordTol,
        linear_shape_function_values,
        eq_ids,
        proj_dist,
        ComputeApproximation);

    std::ignore = GeometricalProjectionUtilities::FastProjectOnLine(
        *p_geom,
        point_to_proj,
        projected_point);

    const bool is_full_projection = (pairing_index_linear == ProjectionUtilities::PairingIndex::Line_Inside);
    if (is_full_projection) {
        SetLocalSearchWasSuccessful();
    } else {
        if (!ComputeApproximation) {
            return;
        }
        SetIsApproximation();
    }

    // Several beam candidates can be visited; keep the best pairing and, for ties, the nearest projection.
    if (pairing_index_linear > mPairingIndex ||
        (pairing_index_linear == mPairingIndex && proj_dist < mClosestProjectionDistance)) {
        mPairingIndex = pairing_index_linear;
        mClosestProjectionDistance = proj_dist;
        mNodeIds = eq_ids;

        const IndexType num_values_linear = linear_shape_function_values.size();
        if (mLinearShapeFunctionValues.size() != num_values_linear) {
            mLinearShapeFunctionValues.resize(num_values_linear);
        }

        for (IndexType i = 0; i < num_values_linear; ++i) {
            mLinearShapeFunctionValues[i] = linear_shape_function_values[i];
        }

        mProjectionOfPoint = projected_point;
        mpInterfaceObject = make_shared<InterfaceGeometryObject>(rInterfaceObject.pGetBaseGeometry());
    }

}

void BeamSplineMapperInterfaceInfo::ComputeRotationMatrix()
{
    std::vector<double> axis_X;
    std::vector<double> axis_Y;
    std::vector<double> axis_Z;

    axis_X.resize(3);
    axis_Y.resize(3);
    axis_Z.resize(3);
    
    const auto p_geom = mpInterfaceObject->pGetBaseGeometry();

    auto temp_v = (*p_geom)[1].Coordinates() - (*p_geom)[0].Coordinates();
    double length_X = sqrt(temp_v[0]*temp_v[0] + temp_v[1]*temp_v[1] + temp_v[2]*temp_v[2]);
    
    KRATOS_ERROR_IF(length_X < 0.000001) << "Lenght of the beam is 0.0" << std::endl;
    
    axis_X[0] = temp_v[0] / length_X;
    axis_X[1] = temp_v[1] / length_X;
    axis_X[2] = temp_v[2] / length_X;   

    if (axis_X[0] == 1.0 && axis_X[1] == 0.0 && axis_X[2] == 0.0 ){
        axis_Y[0] = 0.0;
        axis_Y[1] = 1.0;
        axis_Y[2] = 0.0;
        axis_Z[0] = 0.0;
        axis_Z[1] = 0.0;
        axis_Z[2] = 1.0;
    }
    else if (axis_X[0] == 0.0 && axis_X[1] == 1.0 && axis_X[2] == 0.0 ){
        axis_Y[0] = 0.0;
        axis_Y[1] = 0.0;
        axis_Y[2] = 1.0;
        axis_Z[0] = 1.0;
        axis_Z[1] = 0.0;
        axis_Z[2] = 0.0;
    }
    else if (axis_X[0] == 0.0 && axis_X[1] == 0.0 && axis_X[2] == 1.0 ){
        axis_Y[0] = 0.0;
        axis_Y[1] = 1.0;
        axis_Y[2] = 0.0;
        axis_Z[0] = 1.0;
        axis_Z[1] = 0.0;
        axis_Z[2] = 0.0;
    }
    else if (axis_X[0] != 0.0 && axis_X[1] != 0.0 && axis_X[2] == 0.0 ){
        axis_Y[0] = -axis_X[1];
        axis_Y[1] =  axis_X[0];
        axis_Y[2] =  0.0;
        axis_Z[0] = axis_X[1]*axis_Y[2] - axis_X[2]*axis_Y[1]; 
        axis_Z[1] = axis_X[2]*axis_Y[0] - axis_X[0]*axis_Y[2];
        axis_Z[2] = axis_X[0]*axis_Y[1] - axis_X[1]*axis_Y[0];
    }
    else if (axis_X[0] != 0.0 && axis_X[1] == 0.0 && axis_X[2] != 0.0 ){
        axis_Y[0] = -axis_X[2];
        axis_Y[1] =  0;
        axis_Y[2] =  axis_X[0];
        axis_Z[0] = axis_X[1]*axis_Y[2] - axis_X[2]*axis_Y[1]; 
        axis_Z[1] = axis_X[2]*axis_Y[0] - axis_X[0]*axis_Y[2];
        axis_Z[2] = axis_X[0]*axis_Y[1] - axis_X[1]*axis_Y[0];
    }
    else if (axis_X[0] == 0.0 && axis_X[1] != 0.0 && axis_X[2] != 0.0){
        axis_Y[0] =  0;
        axis_Y[1] = -axis_X[2];
        axis_Y[2] =  axis_X[1];
        axis_Z[0] = axis_X[1]*axis_Y[2] - axis_X[2]*axis_Y[1]; 
        axis_Z[1] = axis_X[2]*axis_Y[0] - axis_X[0]*axis_Y[2];
        axis_Z[2] = axis_X[0]*axis_Y[1] - axis_X[1]*axis_Y[0];
    }
    else{
        axis_Y[0] = 1;
        axis_Y[1] = 1;
        axis_Y[2] = (-axis_X[0] - axis_X[1]) / axis_X[2];
        double length_Y = sqrt(axis_Y[0]*axis_Y[0] + axis_Y[1]*axis_Y[1] + axis_Y[2]*axis_Y[2]);
        axis_Y[0] = axis_Y[0]/length_Y;
        axis_Y[1] = axis_Y[1]/length_Y;
        axis_Y[2] = axis_Y[2]/length_Y;
        
        axis_Z[0] = axis_X[1]*axis_Y[2] - axis_X[2]*axis_Y[1]; 
        axis_Z[1] = axis_X[2]*axis_Y[0] - axis_X[0]*axis_Y[2];
        axis_Z[2] = axis_X[0]*axis_Y[1] - axis_X[1]*axis_Y[0];
    }

    MatrixType rotation_matrix(3, 3, 0.0);

    for(IndexType j = 0; j < 3; j++)
    {
        rotation_matrix(j, 0) = axis_X[j];
        rotation_matrix(j, 1) = axis_Y[j];
        rotation_matrix(j, 2) = axis_Z[j];
    }

    mRotationMatrix_G_L = rotation_matrix;

}

void BeamSplineMapperInterfaceInfo::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MapperInterfaceInfo);
}

void BeamSplineMapperInterfaceInfo::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, MapperInterfaceInfo);
}

void BeamSplineMapperLocalSystem::PairingInfo(std::ostream& rOStream, const int EchoLevel) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not initialized!" << std::endl;
    rOStream << "BeamSplineMapperLocalSystem based on " << mpNode->Info();
}

BeamSplineMapperLocalSystem::ProjectionData BeamSplineMapperLocalSystem::CalculateProjectionData()
{
    ProjectionData projection_data;
    projection_data.pNode = mpNode;
    return projection_data;
}

template<class TSparseSpace, class TDenseSpace>
BeamSplineMapper<TSparseSpace, TDenseSpace>::BeamSplineMapper(
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    Parameters JsonParameters)
    : mrModelPartOrigin(rModelPartOrigin),
      mrModelPartDestination(rModelPartDestination),
      mMapperSettings(JsonParameters)
{
    KRATOS_TRY

    mpInterfaceVectorContainerOrigin = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartOrigin);
    mpInterfaceVectorContainerDestination = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartDestination);

    ValidateInput();
    mLocalCoordTol = JsonParameters["local_coord_tolerance"].GetDouble();

    Initialize();

    KRATOS_CATCH("")
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::UpdateInterface(
    Kratos::Flags MappingOptions,
    double SearchRadius)
{
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::Map(
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
    const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement,
    Kratos::Flags MappingOptions)
{
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::Map(
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement,
    Kratos::Flags MappingOptions)
{
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::InverseMap(
    const Variable<array_1d<double, 3>>& rOriginVariablesForces,
    const Variable<array_1d<double, 3>>& rOriginVariablesMoments,
    const Variable<array_1d<double, 3>>& rDestinationVariableForces,
    Kratos::Flags MappingOptions)
{
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::MapperUniquePointerType
BeamSplineMapper<TSparseSpace, TDenseSpace>::Clone(
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    Parameters JsonParameters) const
{
    return Kratos::make_unique<BeamSplineMapper<TSparseSpace, TDenseSpace>>(
        rModelPartOrigin,
        rModelPartDestination,
        JsonParameters);
}

template<class TSparseSpace, class TDenseSpace>
std::string BeamSplineMapper<TSparseSpace, TDenseSpace>::Info() const
{
    return "BeamSplineMapper";
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "BeamSplineMapper";
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::PrintData(std::ostream& rOStream) const
{
    BaseType::PrintData(rOStream);
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::ValidateInput()
{
    Parameters mapper_default_settings(GetMapperDefaultSettings());
    mMapperSettings.ValidateAndAssignDefaults(mapper_default_settings);
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::Initialize()
{
    mBeamChainCache.clear();
    mNodeIdToBeamChainKey.clear();
    InitializeInterfaceCommunicator();
    InitializeInterface();
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::InitializeInterfaceCommunicator()
{
    mpInterfaceCommunicator = Kratos::make_unique<InterfaceCommunicator>(
        mrModelPartOrigin,
        mMapperLocalSystems,
        mMapperSettings["search_settings"]);
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    mBeamChainCache.clear();
    mNodeIdToBeamChainKey.clear();
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildProblem(Kratos::Flags MappingOptions)
{
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::CreateMapperLocalSystems(
    const Communicator& rModelPartCommunicator,
    std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems)
{
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::MapDisplacements(
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
    const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement,
    Kratos::Flags MappingOptions)
{
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::InitializeInformationBeams(
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
    const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement)
{
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::InitializeInformationBeamsInverse(
    const Variable<array_1d<double, 3>>& rOriginVariablesForces,
    const Variable<array_1d<double, 3>>& rOriginVariablesMoments,
    const Variable<array_1d<double, 3>>& rDestinationVariableForces,
    const Kratos::Flags& rMappingOptions)
{
}

template<class TSparseSpace, class TDenseSpace>
double BeamSplineMapper<TSparseSpace, TDenseSpace>::EvaluateKernelValue(const double Distance) const
{
    return 0.0;
}

template<class TSparseSpace, class TDenseSpace>
double BeamSplineMapper<TSparseSpace, TDenseSpace>::EvaluateKernelFirstDerivative(const double Distance) const
{
    return 0.0;
}

template<class TSparseSpace, class TDenseSpace>
double BeamSplineMapper<TSparseSpace, TDenseSpace>::EvaluateKernelSecondDerivative(const double Distance) const
{
    return 0.0;
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildSplineSystemMatrices(
    const std::vector<double>& rSourceCoordinates,
    MatrixType& rMss,
    MatrixType& rMssFirstDerivative,
    MatrixType& rMssSecondDerivative,
    MatrixType& rPs,
    MatrixType& rPsFirstDerivative) const
{
    rMss.resize(0, 0, false);
    rMssFirstDerivative.resize(0, 0, false);
    rMssSecondDerivative.resize(0, 0, false);
    rPs.resize(0, 0, false);
    rPsFirstDerivative.resize(0, 0, false);
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildEvaluationRow(
    const std::vector<double>& rSourceCoordinates,
    const double ProjectionCoordinate) const
{
    return VectorType();
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildRightHandSide(
    const std::vector<double>& rDisplacements,
    const std::vector<double>& rRotations) const
{
    return VectorType();
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildLocalSourceData(
    const BeamChainCacheData& rBeamChainCacheData,
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
    typename BeamSplineMapper<TSparseSpace, TDenseSpace>::ComponentArrayType& rLocalDisplacements,
    typename BeamSplineMapper<TSparseSpace, TDenseSpace>::ComponentArrayType& rLocalRotations) const
{
    for (IndexType i = 0; i < 3; ++i) {
        rLocalDisplacements[i].clear();
        rLocalRotations[i].clear();
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::PrepareBeamChainCacheData()
{
}

template<class TSparseSpace, class TDenseSpace>
std::unordered_map<std::string, typename BeamSplineMapper<TSparseSpace, TDenseSpace>::BeamChainSourceStateData>
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildAllBeamChainSourceStates(
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rOriginVariablesRotations) const
{
    return {};
}

template<class TSparseSpace, class TDenseSpace>
const typename BeamSplineMapper<TSparseSpace, TDenseSpace>::BeamChainCacheData&
BeamSplineMapper<TSparseSpace, TDenseSpace>::GetOrCreateBeamChainCacheData(
    const GeometryType& rBeamGeometry)
{
    return mBeamChainCache.emplace(std::string(), BeamChainCacheData()).first->second;
}

template<class TSparseSpace, class TDenseSpace>
std::string BeamSplineMapper<TSparseSpace, TDenseSpace>::CreateBeamChainKey(
    const std::vector<IndexType>& rSupportNodeIds) const
{
    return std::string();
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::MatrixType
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildSplineSystemMatrix(
    const std::vector<double>& rSourceCoordinates) const
{
    return MatrixType();
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::ComputeSupportReferenceData(
    const std::vector<IndexType>& rSupportNodeIds,
    std::vector<double>& rLocalXCoordinates,
    std::vector<MatrixType>& rSupportFramesLocalToGlobal,
    std::vector<array_1d<double, 3>>& rSupportReferenceCoordinates) const
{
    rLocalXCoordinates.clear();
    rSupportFramesLocalToGlobal.clear();
    rSupportReferenceCoordinates.clear();
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::MatrixType
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildEvaluationFrameLocalToGlobal(
    const std::vector<IndexType>& rSupportNodeIds,
    const std::unordered_map<IndexType, IndexType>& rSupportNodeIdToLocalIndex,
    const std::vector<MatrixType>& rSupportFramesLocalToGlobal,
    const std::vector<array_1d<double, 3>>& rSupportReferenceCoordinates,
    const GeometryType& rBeamGeometry,
    const VectorType& rLinearShapeValues) const
{
    return MatrixType();
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::ComputeBeamChainSupport(
    const GeometryType& rBeamGeometry,
    std::vector<IndexType>& rSupportNodeIds,
    std::unordered_map<IndexType, IndexType>& rSupportNodeIdToLocalIndex) const
{
    rSupportNodeIds.clear();
    rSupportNodeIdToLocalIndex.clear();
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::InitializeOriginForcesAndMoments(
    const Variable<array_1d<double, 3>>& rOriginVariablesForces,
    const Variable<array_1d<double, 3>>& rOriginVariablesMoments)
{
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::SolveSplineCoefficients(
    const MatrixType& rSplineSystemMatrix,
    const VectorType& rRightHandSide) const
{
    return VectorType();
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::TransformGlobalToLocal(
    const MatrixType& rRotationMatrixGlobalToLocal,
    const array_1d<double, 3>& rReferencePoint,
    const array_1d<double, 3>& rCoordinates) const
{
    return VectorType();
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::TransformVectorToLocal(
    const MatrixType& rRotationMatrixGlobalToLocal,
    const array_1d<double, 3>& rVector) const
{
    return VectorType();
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::TransformVectorToGlobal(
    const MatrixType& rRotationMatrixLocalToGlobal,
    const VectorType& rLocalVector) const
{
    return VectorType();
}

template<class TSparseSpace, class TDenseSpace>
array_1d<double, 3> BeamSplineMapper<TSparseSpace, TDenseSpace>::GetReferenceCoordinates(
    const Node& rNode) const
{
    return ZeroVector(3);
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::MatrixType
BeamSplineMapper<TSparseSpace, TDenseSpace>::ComputeInitialFrameLocalToGlobal(
    const array_1d<double, 3>& rTangent) const
{
    return MatrixType();
}

template<class TSparseSpace, class TDenseSpace>
array_1d<double, 3> BeamSplineMapper<TSparseSpace, TDenseSpace>::RotateVectorAroundAxis(
    const array_1d<double, 3>& rVector,
    const array_1d<double, 3>& rAxis,
    const double Angle) const
{
    return ZeroVector(3);
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::EvaluatePointDisplacementLocal(
    const VectorType& rEvaluationRow,
    const VectorType& rSplineCoefficientsY,
    const VectorType& rSplineCoefficientsZ,
    const double AxialDisplacement,
    const double TorsionalRotation,
    const VectorType& rOffsetVectorLocal) const
{
    return VectorType();
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::MapperInterfaceInfoUniquePointerType
BeamSplineMapper<TSparseSpace, TDenseSpace>::GetMapperInterfaceInfo() const
{
    return Kratos::make_unique<BeamSplineMapperInterfaceInfo>(mLocalCoordTol);
}

template<class TSparseSpace, class TDenseSpace>
Parameters BeamSplineMapper<TSparseSpace, TDenseSpace>::GetMapperDefaultSettings() const
{
    return Parameters(R"({
        "search_settings"              : {},
        "search_radius"                : -1.0,
        "search_iterations"            : 3,
        "local_coord_tolerance"        : 0.25,
        "echo_level"                   : 0
    })");
}

template class BeamSplineMapper<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>;

}  // namespace Kratos
