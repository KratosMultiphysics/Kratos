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

namespace
{
constexpr double PolynomialTolerance = 1.0e-12;
constexpr double FrameTolerance = 1.0e-10;
constexpr std::array<const char*, 3> ComponentSuffixes{{"_X", "_Y", "_Z"}};

template<class TContainerType>
void AddBeamConnectivity(
    const TContainerType& rContainer,
    std::unordered_map<IndexType, std::unordered_set<IndexType>>& rAdjacency)
{
    for (const auto& r_entity : rContainer) {
        const auto& r_geometry = r_entity.GetGeometry();
        if (r_geometry.size() != 2) {
            continue;
        }

        const IndexType node_id_i = r_geometry[0].Id();
        const IndexType node_id_j = r_geometry[1].Id();
        rAdjacency[node_id_i].insert(node_id_j);
        rAdjacency[node_id_j].insert(node_id_i);
    }
}

template<class TMatrixType>
TMatrixType TransposeRotationMatrix(const TMatrixType& rRotationMatrixLocalToGlobal)
{
    TMatrixType rotation_matrix_global_to_local(3, 3);
    for (IndexType i = 0; i < 3; ++i) {
        for (IndexType j = 0; j < 3; ++j) {
            rotation_matrix_global_to_local(i, j) = rRotationMatrixLocalToGlobal(j, i);
        }
    }
    return rotation_matrix_global_to_local;
}

array_1d<double, 3> CrossProduct(
    const array_1d<double, 3>& rA,
    const array_1d<double, 3>& rB)
{
    array_1d<double, 3> cross_product = ZeroVector(3);
    cross_product[0] = rA[1] * rB[2] - rA[2] * rB[1];
    cross_product[1] = rA[2] * rB[0] - rA[0] * rB[2];
    cross_product[2] = rA[0] * rB[1] - rA[1] * rB[0];
    return cross_product;
}

double DotProduct(
    const array_1d<double, 3>& rA,
    const array_1d<double, 3>& rB)
{
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}

array_1d<double, 3> NormalizeVector(const array_1d<double, 3>& rVector)
{
    const double norm = norm_2(rVector);
    KRATOS_ERROR_IF(norm < FrameTolerance) << "Cannot normalize vector with near-zero norm." << std::endl;

    array_1d<double, 3> normalized_vector = rVector;
    normalized_vector /= norm;
    return normalized_vector;
}

struct ReferenceLineProjection
{
    Point ProjectionPoint;
    Vector ShapeValues{2};
    double LocalCoordinate = 0.0;
    double Distance = 0.0;
    bool IsInside = false;
};

ReferenceLineProjection ProjectOnReferenceLine(
    const Geometry<Node>& rGeometry,
    const Point& rPoint,
    const bool ClampToSegment)
{
    KRATOS_ERROR_IF(rGeometry.size() != 2)
        << "BeamSplineMapper expects a two-node line geometry." << std::endl;

    array_1d<double, 3> x0{{rGeometry[0].X0(), rGeometry[0].Y0(), rGeometry[0].Z0()}};
    array_1d<double, 3> x1{{rGeometry[1].X0(), rGeometry[1].Y0(), rGeometry[1].Z0()}};
    const array_1d<double, 3> direction = x1 - x0;
    const double length_squared = inner_prod(direction, direction);
    KRATOS_ERROR_IF(length_squared <= FrameTolerance * FrameTolerance)
        << "Cannot project on a zero-length reference beam." << std::endl;

    array_1d<double, 3> relative = rPoint.Coordinates() - x0;
    double segment_coordinate = inner_prod(relative, direction) / length_squared;
    const double unclamped_coordinate = segment_coordinate;
    if (ClampToSegment) {
        segment_coordinate = std::clamp(segment_coordinate, 0.0, 1.0);
    }

    ReferenceLineProjection result;
    result.ShapeValues(0) = 1.0 - segment_coordinate;
    result.ShapeValues(1) = segment_coordinate;
    result.LocalCoordinate = 2.0 * segment_coordinate - 1.0;
    result.IsInside = unclamped_coordinate >= 0.0 && unclamped_coordinate <= 1.0;
    noalias(result.ProjectionPoint.Coordinates()) =
        result.ShapeValues(0) * x0 + result.ShapeValues(1) * x1;
    result.Distance = norm_2(rPoint.Coordinates() - result.ProjectionPoint.Coordinates());
    return result;
}

Matrix BuildReferenceFrame(const Geometry<Node>& rGeometry)
{
    array_1d<double, 3> tangent{{
        rGeometry[1].X0() - rGeometry[0].X0(),
        rGeometry[1].Y0() - rGeometry[0].Y0(),
        rGeometry[1].Z0() - rGeometry[0].Z0()}};
    tangent = NormalizeVector(tangent);

    // Select the Cartesian axis least aligned with the tangent.  Projection
    // followed by cross products is continuous away from the unavoidable
    // equal-alignment ties and always produces a right-handed frame.
    IndexType reference_index = 0;
    if (std::abs(tangent[1]) < std::abs(tangent[reference_index])) {
        reference_index = 1;
    }
    if (std::abs(tangent[2]) < std::abs(tangent[reference_index])) {
        reference_index = 2;
    }
    array_1d<double, 3> normal = ZeroVector(3);
    normal[reference_index] = 1.0;
    normal -= DotProduct(normal, tangent) * tangent;
    normal = NormalizeVector(normal);
    array_1d<double, 3> binormal = NormalizeVector(CrossProduct(tangent, normal));
    normal = NormalizeVector(CrossProduct(binormal, tangent));

    Matrix frame(3, 3, 0.0);
    for (IndexType i = 0; i < 3; ++i) {
        frame(i, 0) = tangent[i];
        frame(i, 1) = normal[i];
        frame(i, 2) = binormal[i];
    }
    return frame;
}

array_1d<double, 3> MakeArrayFromVector(const Vector& rVector)
{
    array_1d<double, 3> coordinates = ZeroVector(3);
    coordinates[0] = rVector(0);
    coordinates[1] = rVector(1);
    coordinates[2] = rVector(2);
    return coordinates;
}

const Variable<double>& GetComponentVariable(
    const Variable<array_1d<double, 3>>& rVariable,
    const IndexType ComponentIndex)
{
    return KratosComponents<Variable<double>>::Get(rVariable.Name() + ComponentSuffixes[ComponentIndex]);
}
}

void BeamSplineMapperInterfaceInfo::ProcessSearchResult(const InterfaceObject& rInterfaceObject)
{
    SaveSearchResult(rInterfaceObject, false);
}

void BeamSplineMapperInterfaceInfo::ProcessSearchResultForApproximation(const InterfaceObject& rInterfaceObject)
{
    const auto p_geom = rInterfaceObject.pGetBaseGeometry();
    const Point point_to_project(this->Coordinates());
    mCoordinates = point_to_project;
    const auto projection = ProjectOnReferenceLine(*p_geom, point_to_project, true);

    if (ProjectionUtilities::PairingIndex::Closest_Point > mPairingIndex ||
        (ProjectionUtilities::PairingIndex::Closest_Point == mPairingIndex &&
         projection.Distance < mClosestProjectionDistance)) {
        SetIsApproximation();
        mPairingIndex = ProjectionUtilities::PairingIndex::Closest_Point;
        mClosestProjectionDistance = projection.Distance;
        mProjectionOfPoint = projection.ProjectionPoint;
        mNodeIds.resize(2);
        mLinearShapeFunctionValues.resize(2);
        for (IndexType i = 0; i < 2; ++i) {
            KRATOS_DEBUG_ERROR_IF_NOT((*p_geom)[i].Has(INTERFACE_EQUATION_ID))
                << (*p_geom)[i] << " does not have an INTERFACE_EQUATION_ID." << std::endl;
            mNodeIds[i] = (*p_geom)[i].GetValue(INTERFACE_EQUATION_ID);
            mLinearShapeFunctionValues[i] = projection.ShapeValues(i);
        }
        mpInterfaceObject = make_shared<InterfaceGeometryObject>(rInterfaceObject.pGetBaseGeometry());
    }
}

void BeamSplineMapperInterfaceInfo::SaveSearchResult(
    const InterfaceObject& rInterfaceObject,
    const bool ComputeApproximation)
{
    const auto p_geom = rInterfaceObject.pGetBaseGeometry();

    const Point point_to_project(this->Coordinates());
    mCoordinates = point_to_project;
    const auto projection = ProjectOnReferenceLine(*p_geom, point_to_project, ComputeApproximation);
    const bool is_full_projection = projection.IsInside;
    if (!is_full_projection && !ComputeApproximation) {
        return;
    }

    const auto pairing_index = is_full_projection
        ? ProjectionUtilities::PairingIndex::Line_Inside
        : ProjectionUtilities::PairingIndex::Closest_Point;

    if (is_full_projection) {
        SetLocalSearchWasSuccessful();
    } else {
        SetIsApproximation();
    }

    if (pairing_index > mPairingIndex ||
        (pairing_index == mPairingIndex && projection.Distance < mClosestProjectionDistance)) {
        mPairingIndex = pairing_index;
        mClosestProjectionDistance = projection.Distance;
        mProjectionOfPoint = projection.ProjectionPoint;
        mNodeIds.resize(2);
        mLinearShapeFunctionValues.resize(2);
        for (IndexType i = 0; i < 2; ++i) {
            KRATOS_DEBUG_ERROR_IF_NOT((*p_geom)[i].Has(INTERFACE_EQUATION_ID))
                << (*p_geom)[i] << " does not have an INTERFACE_EQUATION_ID." << std::endl;
            mNodeIds[i] = (*p_geom)[i].GetValue(INTERFACE_EQUATION_ID);
            mLinearShapeFunctionValues[i] = projection.ShapeValues(i);
        }
        mpInterfaceObject = make_shared<InterfaceGeometryObject>(rInterfaceObject.pGetBaseGeometry());
    }
}

void BeamSplineMapperInterfaceInfo::ComputeRotationMatrix()
{
    KRATOS_ERROR_IF_NOT(mpInterfaceObject)
        << "Cannot build a beam frame without an interface object." << std::endl;
    const auto p_geom = mpInterfaceObject->pGetBaseGeometry();
    mRotationMatrixLocalToGlobal = BuildReferenceFrame(*p_geom);
}

void BeamSplineMapperInterfaceInfo::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MapperInterfaceInfo);
    rSerializer.save("NodeIds", mNodeIds);
    rSerializer.save("LinearShapeFunctionValues", mLinearShapeFunctionValues);
    rSerializer.save("ProjectionOfPoint", mProjectionOfPoint);
    rSerializer.save("ClosestProjectionDistance", mClosestProjectionDistance);
    rSerializer.save("PairingIndex", static_cast<int>(mPairingIndex));
}

void BeamSplineMapperInterfaceInfo::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, MapperInterfaceInfo);
    rSerializer.load("NodeIds", mNodeIds);
    rSerializer.load("LinearShapeFunctionValues", mLinearShapeFunctionValues);
    rSerializer.load("ProjectionOfPoint", mProjectionOfPoint);
    rSerializer.load("ClosestProjectionDistance", mClosestProjectionDistance);
    int pairing_index = 0;
    rSerializer.load("PairingIndex", pairing_index);
    mPairingIndex = static_cast<ProjectionUtilities::PairingIndex>(pairing_index);
}

void BeamSplineMapperLocalSystem::PairingInfo(std::ostream& rOStream, const int EchoLevel) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not initialized!" << std::endl;

    rOStream << "BeamSplineMapperLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 1) {
        rOStream << " at Coordinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    }
}

BeamSplineMapperLocalSystem::ProjectionData BeamSplineMapperLocalSystem::CalculateProjectionData()
{
    ProjectionData projection_data;
    projection_data.pNode = mpNode;

    for (auto& r_interface_info : mInterfaceInfos) {
        auto beam_interface_info = std::dynamic_pointer_cast<BeamSplineMapperInterfaceInfo>(r_interface_info);
        KRATOS_ERROR_IF_NOT(beam_interface_info)
            << "Expected BeamSplineMapperInterfaceInfo in mInterfaceInfos but got nullptr or wrong type." << std::endl;

        beam_interface_info->ComputeRotationMatrixInterfaceObject();
        beam_interface_info->GetValue(
            projection_data.RotationMatrixLocalToGlobal,
            projection_data.ProjectionPoint,
            projection_data.LinearShapeValues,
            projection_data.BeamGeometry);
    }

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
    KRATOS_ERROR_IF(mLocalCoordTol < 0.0) << "The local_coord_tolerance cannot be negative." << std::endl;

    Initialize();

    KRATOS_CATCH("")
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::UpdateInterface(
    Kratos::Flags MappingOptions,
    double SearchRadius)
{
    if (SearchRadius > 0.0) {
        mMapperSettings["search_settings"]["search_radius"].SetDouble(SearchRadius);
    }

    if (MappingOptions.Is(MapperFlags::REMESHED)) {
        InitializeInterface(MappingOptions);
    } else {
        BuildProblem(MappingOptions);
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::Map(
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
    const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement,
    Kratos::Flags MappingOptions)
{
    MapDisplacements(
        rOriginVariablesDisplacements,
        rOriginVariablesRotations,
        rDestinationVariableDisplacement,
        MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::Map(
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement,
    Kratos::Flags MappingOptions)
{
    MapDisplacements(
        rOriginVariablesDisplacements,
        ROTATION,
        rDestinationVariableDisplacement,
        MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::InverseMap(
    const Variable<array_1d<double, 3>>& rOriginVariablesForces,
    const Variable<array_1d<double, 3>>& rOriginVariablesMoments,
    const Variable<array_1d<double, 3>>& rDestinationVariableForces,
    Kratos::Flags MappingOptions)
{
    InitializeOriginForcesAndMoments(rOriginVariablesForces, rOriginVariablesMoments);
    InitializeInformationBeamsInverse(
        rOriginVariablesForces,
        rOriginVariablesMoments,
        rDestinationVariableForces,
        MappingOptions);
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

    if (mMapperSettings["search_radius"].GetDouble() < 0.0) {
        const double search_radius = MapperUtilities::ComputeSearchRadius(
            mrModelPartOrigin,
            mrModelPartDestination,
            mMapperSettings["echo_level"].GetInt());
        mMapperSettings["search_radius"].SetDouble(search_radius);
    }
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
    CreateMapperLocalSystems(mrModelPartDestination.GetCommunicator(), mMapperLocalSystems);
    BuildProblem(MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildProblem(Kratos::Flags MappingOptions)
{
    MapperUtilities::AssignInterfaceEquationIdsToNodes(mrModelPartOrigin.GetCommunicator());
    MapperUtilities::AssignInterfaceEquationIdsToNodes(mrModelPartDestination.GetCommunicator());

    KRATOS_ERROR_IF_NOT(mpInterfaceCommunicator) << "mpInterfaceCommunicator is a nullptr." << std::endl;
    const MapperInterfaceInfoUniquePointerType p_ref_interface_info = GetMapperInterfaceInfo();
    mpInterfaceCommunicator->ExchangeInterfaceData(
        mrModelPartDestination.GetCommunicator(),
        p_ref_interface_info);

    PrepareBeamChainCacheData();
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::CreateMapperLocalSystems(
    const Communicator& rModelPartCommunicator,
    std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems)
{
    MapperUtilities::CreateMapperLocalSystemsFromNodes(
        BeamSplineMapperLocalSystem(nullptr),
        rModelPartCommunicator,
        rLocalSystems);
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::MapDisplacements(
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
    const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement,
    Kratos::Flags MappingOptions)
{
    InitializeInformationBeams(
        rOriginVariablesDisplacements,
        rOriginVariablesRotations,
        rDestinationVariableDisplacement);
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::InitializeInformationBeams(
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
    const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement)
{
    const auto beam_chain_source_states = BuildAllBeamChainSourceStates(
        rOriginVariablesDisplacements,
        rOriginVariablesRotations);

    for (auto& r_local_sys : mMapperLocalSystems) {
        if (!r_local_sys->HasInterfaceInfo()) {
            continue;
        }

        auto beam_sys = dynamic_cast<BeamSplineMapperLocalSystem*>(r_local_sys.get());
        KRATOS_ERROR_IF_NOT(beam_sys) << "Expected BeamSplineMapperLocalSystem!" << std::endl;

        auto projection_data = beam_sys->CalculateProjectionData();

        KRATOS_ERROR_IF_NOT(projection_data.pNode) << "Destination node is a nullptr." << std::endl;

        const auto& r_beam_chain_cache_data = GetOrCreateBeamChainCacheData(projection_data.BeamGeometry);
        const auto source_state_iterator = beam_chain_source_states.find(r_beam_chain_cache_data.Key);
        KRATOS_ERROR_IF(source_state_iterator == beam_chain_source_states.end())
            << "Missing source-state data for beam chain '" << r_beam_chain_cache_data.Key << "'." << std::endl;
        const auto& r_source_state_data = source_state_iterator->second;

        const MatrixType evaluation_rotation_matrix_local_to_global =
            BuildEvaluationFrameLocalToGlobal(
                r_beam_chain_cache_data.SupportNodeIds,
                r_beam_chain_cache_data.SupportNodeIdToLocalIndex,
                r_beam_chain_cache_data.SupportFramesLocalToGlobal,
                r_beam_chain_cache_data.SupportReferenceCoordinates,
                projection_data.BeamGeometry,
                projection_data.LinearShapeValues);
        const MatrixType evaluation_rotation_matrix_global_to_local =
            TransposeRotationMatrix(evaluation_rotation_matrix_local_to_global);

        const IndexType first_segment_node_index =
            r_beam_chain_cache_data.SupportNodeIdToLocalIndex.at(projection_data.BeamGeometry[0].Id());
        const IndexType second_segment_node_index =
            r_beam_chain_cache_data.SupportNodeIdToLocalIndex.at(projection_data.BeamGeometry[1].Id());
        const double projection_coordinate =
            projection_data.LinearShapeValues(0) * r_beam_chain_cache_data.LocalXCoordinates[first_segment_node_index] +
            projection_data.LinearShapeValues(1) * r_beam_chain_cache_data.LocalXCoordinates[second_segment_node_index];
        const VectorType evaluation_row = BuildEvaluationRow(
            r_beam_chain_cache_data.LocalXCoordinates,
            projection_coordinate);
        const VectorType evaluation_derivative_row = BuildEvaluationDerivativeRow(
            r_beam_chain_cache_data.LocalXCoordinates,
            projection_coordinate,
            r_beam_chain_cache_data.CoordinateHalfLength);

        const double axial_displacement =
            projection_data.LinearShapeValues(0) * r_source_state_data.LocalDisplacements[0][first_segment_node_index] +
            projection_data.LinearShapeValues(1) * r_source_state_data.LocalDisplacements[0][second_segment_node_index];

        const double torsional_rotation =
            projection_data.LinearShapeValues(0) * r_source_state_data.LocalRotations[0][first_segment_node_index] +
            projection_data.LinearShapeValues(1) * r_source_state_data.LocalRotations[0][second_segment_node_index];
        VectorType evaluation_rotation_vector(3);
        evaluation_rotation_vector(0) = torsional_rotation;
        evaluation_rotation_vector(1) = -inner_prod(
            evaluation_derivative_row,
            r_source_state_data.SplineCoefficientsZ);
        evaluation_rotation_vector(2) = inner_prod(
            evaluation_derivative_row,
            r_source_state_data.SplineCoefficientsY);
        beam_sys->SaveEvaluationRotationVector(evaluation_rotation_vector);

        const array_1d<double, 3> projection_point_reference = MakeArrayFromVector(projection_data.ProjectionPoint);
        const array_1d<double, 3> destination_point_reference = GetReferenceCoordinates(*projection_data.pNode);
        const VectorType offset_vector_local = TransformGlobalToLocal(
            evaluation_rotation_matrix_global_to_local,
            projection_point_reference,
            destination_point_reference);

        const VectorType local_displacement = EvaluatePointDisplacementLocal(
            evaluation_row,
            evaluation_derivative_row,
            r_source_state_data.SplineCoefficientsY,
            r_source_state_data.SplineCoefficientsZ,
            axial_displacement,
            torsional_rotation,
            offset_vector_local);

        const VectorType global_displacement = TransformVectorToGlobal(
            evaluation_rotation_matrix_local_to_global,
            local_displacement);

        for (IndexType i = 0; i < 3; ++i) {
            projection_data.pNode->FastGetSolutionStepValue(GetComponentVariable(rDestinationVariableDisplacement, i)) =
                global_displacement(i);
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::InitializeInformationBeamsInverse(
    const Variable<array_1d<double, 3>>& rOriginVariablesForces,
    const Variable<array_1d<double, 3>>& rOriginVariablesMoments,
    const Variable<array_1d<double, 3>>& rDestinationVariableForces,
    const Kratos::Flags& rMappingOptions)
{
    const double factor = rMappingOptions.Is(MapperFlags::SWAP_SIGN) ? -1.0 : 1.0;

    for (auto& r_local_sys : mMapperLocalSystems) {
        if (!r_local_sys->HasInterfaceInfo()) {
            continue;
        }

        auto beam_sys = dynamic_cast<BeamSplineMapperLocalSystem*>(r_local_sys.get());
        KRATOS_ERROR_IF_NOT(beam_sys) << "Expected BeamSplineMapperLocalSystem!" << std::endl;

        auto projection_data = beam_sys->CalculateProjectionData();

        KRATOS_ERROR_IF_NOT(projection_data.pNode) << "Destination node is a nullptr." << std::endl;
        KRATOS_ERROR_IF(projection_data.BeamGeometry.size() != 2)
            << "BeamSplineMapper inverse load transfer expects a 2-node beam segment." << std::endl;
        KRATOS_ERROR_IF(projection_data.LinearShapeValues.size() != 2)
            << "BeamSplineMapper inverse load transfer expects two linear shape function values." << std::endl;

        const auto& r_beam_chain_cache_data = GetOrCreateBeamChainCacheData(projection_data.BeamGeometry);

        const MatrixType evaluation_rotation_matrix_local_to_global =
            BuildEvaluationFrameLocalToGlobal(
                r_beam_chain_cache_data.SupportNodeIds,
                r_beam_chain_cache_data.SupportNodeIdToLocalIndex,
                r_beam_chain_cache_data.SupportFramesLocalToGlobal,
                r_beam_chain_cache_data.SupportReferenceCoordinates,
                projection_data.BeamGeometry,
                projection_data.LinearShapeValues);
        const MatrixType evaluation_rotation_matrix_global_to_local =
            TransposeRotationMatrix(evaluation_rotation_matrix_local_to_global);

        const IndexType first_segment_node_index =
            r_beam_chain_cache_data.SupportNodeIdToLocalIndex.at(projection_data.BeamGeometry[0].Id());
        const IndexType second_segment_node_index =
            r_beam_chain_cache_data.SupportNodeIdToLocalIndex.at(projection_data.BeamGeometry[1].Id());
        const double projection_coordinate =
            projection_data.LinearShapeValues(0) * r_beam_chain_cache_data.LocalXCoordinates[first_segment_node_index] +
            projection_data.LinearShapeValues(1) * r_beam_chain_cache_data.LocalXCoordinates[second_segment_node_index];
        const VectorType evaluation_row = BuildEvaluationRow(
            r_beam_chain_cache_data.LocalXCoordinates,
            projection_coordinate);
        const VectorType evaluation_derivative_row = BuildEvaluationDerivativeRow(
            r_beam_chain_cache_data.LocalXCoordinates,
            projection_coordinate,
            r_beam_chain_cache_data.CoordinateHalfLength);

        array_1d<double, 3> surface_force_global = ZeroVector(3);
        for (IndexType i = 0; i < 3; ++i) {
            surface_force_global[i] =
                projection_data.pNode->FastGetSolutionStepValue(GetComponentVariable(rDestinationVariableForces, i));
        }
        const VectorType surface_force_local = TransformVectorToLocal(
            evaluation_rotation_matrix_global_to_local,
            surface_force_global);

        const array_1d<double, 3> projection_point_reference = MakeArrayFromVector(projection_data.ProjectionPoint);
        const array_1d<double, 3> destination_point_reference = GetReferenceCoordinates(*projection_data.pNode);
        const VectorType offset_vector_local = TransformGlobalToLocal(
            evaluation_rotation_matrix_global_to_local,
            projection_point_reference,
            destination_point_reference);

        VectorType evaluation_rotation_vector(3);
        beam_sys->GetEvaluationRotationVector(evaluation_rotation_vector);
        const double theta_x = evaluation_rotation_vector(0);
        const double theta_y = evaluation_rotation_vector(1);
        const double theta_z = evaluation_rotation_vector(2);

        const double cos_x = std::cos(theta_x);
        const double sin_x = std::sin(theta_x);
        const double cos_y = std::cos(theta_y);
        const double sin_y = std::sin(theta_y);
        const double cos_z = std::cos(theta_z);
        const double sin_z = std::sin(theta_z);

        MatrixType rotation_x(3, 3, 0.0);
        rotation_x(0, 0) = 1.0;
        rotation_x(1, 1) = cos_x;
        rotation_x(1, 2) = -sin_x;
        rotation_x(2, 1) = sin_x;
        rotation_x(2, 2) = cos_x;

        MatrixType rotation_y(3, 3, 0.0);
        rotation_y(0, 0) = cos_y;
        rotation_y(0, 2) = sin_y;
        rotation_y(1, 1) = 1.0;
        rotation_y(2, 0) = -sin_y;
        rotation_y(2, 2) = cos_y;

        MatrixType rotation_z(3, 3, 0.0);
        rotation_z(0, 0) = cos_z;
        rotation_z(0, 1) = -sin_z;
        rotation_z(1, 0) = sin_z;
        rotation_z(1, 1) = cos_z;
        rotation_z(2, 2) = 1.0;

        MatrixType derivative_rotation_x(3, 3, 0.0);
        derivative_rotation_x(1, 1) = -sin_x;
        derivative_rotation_x(1, 2) = -cos_x;
        derivative_rotation_x(2, 1) = cos_x;
        derivative_rotation_x(2, 2) = -sin_x;

        MatrixType derivative_rotation_y(3, 3, 0.0);
        derivative_rotation_y(0, 0) = -sin_y;
        derivative_rotation_y(0, 2) = cos_y;
        derivative_rotation_y(2, 0) = -cos_y;
        derivative_rotation_y(2, 2) = -sin_y;

        MatrixType derivative_rotation_z(3, 3, 0.0);
        derivative_rotation_z(0, 0) = -sin_z;
        derivative_rotation_z(0, 1) = -cos_z;
        derivative_rotation_z(1, 0) = cos_z;
        derivative_rotation_z(1, 1) = -sin_z;

        const MatrixType rotation_z_y = prod(rotation_z, rotation_y);
        const MatrixType rotation_x_z = prod(rotation_x, rotation_z);
        const MatrixType derivative_rotation_z_y = prod(derivative_rotation_z, rotation_y);
        const MatrixType derivative_matrix_x = prod(derivative_rotation_x, rotation_z_y);
        const MatrixType derivative_matrix_y = prod(rotation_x_z, derivative_rotation_y);
        const MatrixType derivative_matrix_z = prod(rotation_x, derivative_rotation_z_y);

        VectorType derivative_offset_x(3, 0.0);
        VectorType derivative_offset_y(3, 0.0);
        VectorType derivative_offset_z(3, 0.0);
        TDenseSpace::Mult(derivative_matrix_x, offset_vector_local, derivative_offset_x);
        TDenseSpace::Mult(derivative_matrix_y, offset_vector_local, derivative_offset_y);
        TDenseSpace::Mult(derivative_matrix_z, offset_vector_local, derivative_offset_z);

        VectorType rotation_work_conjugate(3);
        rotation_work_conjugate(0) = inner_prod(surface_force_local, derivative_offset_x);
        rotation_work_conjugate(1) = inner_prod(surface_force_local, derivative_offset_y);
        rotation_work_conjugate(2) = inner_prod(surface_force_local, derivative_offset_z);

        const IndexType number_of_support_nodes = r_beam_chain_cache_data.SupportNodeIds.size();
        ComponentArrayType local_forces;
        ComponentArrayType local_moments;
        for (IndexType component_index = 0; component_index < 3; ++component_index) {
            local_forces[component_index].assign(number_of_support_nodes, 0.0);
            local_moments[component_index].assign(number_of_support_nodes, 0.0);
        }

        local_forces[0][first_segment_node_index] +=
            projection_data.LinearShapeValues(0) * surface_force_local(0);
        local_forces[0][second_segment_node_index] +=
            projection_data.LinearShapeValues(1) * surface_force_local(0);
        local_moments[0][first_segment_node_index] +=
            projection_data.LinearShapeValues(0) * rotation_work_conjugate(0);
        local_moments[0][second_segment_node_index] +=
            projection_data.LinearShapeValues(1) * rotation_work_conjugate(0);

        const IndexType spline_system_size = r_beam_chain_cache_data.SplineSystemMatrix.size1();
        VectorType right_hand_side_y(spline_system_size, 0.0);
        VectorType right_hand_side_z(spline_system_size, 0.0);
        for (IndexType i = 0; i < spline_system_size; ++i) {
            right_hand_side_y(i) =
                evaluation_row(i) * surface_force_local(1) +
                evaluation_derivative_row(i) * rotation_work_conjugate(2);
            right_hand_side_z(i) =
                evaluation_row(i) * surface_force_local(2) -
                evaluation_derivative_row(i) * rotation_work_conjugate(1);
        }

        VectorType adjoint_coefficients_y(spline_system_size, 0.0);
        VectorType adjoint_coefficients_z(spline_system_size, 0.0);
        MathUtils<double>::Solve(
            r_beam_chain_cache_data.TransposedSplineSystemMatrix,
            adjoint_coefficients_y,
            right_hand_side_y);
        MathUtils<double>::Solve(
            r_beam_chain_cache_data.TransposedSplineSystemMatrix,
            adjoint_coefficients_z,
            right_hand_side_z);

        for (IndexType i = 0; i < number_of_support_nodes; ++i) {
            local_forces[1][i] += adjoint_coefficients_y(i);
            local_moments[2][i] += r_beam_chain_cache_data.CoordinateHalfLength *
                adjoint_coefficients_y(number_of_support_nodes + i);
            local_forces[2][i] += adjoint_coefficients_z(i);
            local_moments[1][i] -= r_beam_chain_cache_data.CoordinateHalfLength *
                adjoint_coefficients_z(number_of_support_nodes + i);
        }

        for (IndexType i = 0; i < number_of_support_nodes; ++i) {
            VectorType nodal_force_local(3);
            VectorType nodal_moment_local(3);
            for (IndexType component_index = 0; component_index < 3; ++component_index) {
                nodal_force_local(component_index) = factor * local_forces[component_index][i];
                nodal_moment_local(component_index) = factor * local_moments[component_index][i];
            }

            const VectorType nodal_force_global = TransformVectorToGlobal(
                r_beam_chain_cache_data.SupportFramesLocalToGlobal[i],
                nodal_force_local);
            const VectorType nodal_moment_global = TransformVectorToGlobal(
                r_beam_chain_cache_data.SupportFramesLocalToGlobal[i],
                nodal_moment_local);

            auto& r_node = mrModelPartOrigin.GetNode(r_beam_chain_cache_data.SupportNodeIds[i]);
            for (IndexType component_index = 0; component_index < 3; ++component_index) {
                const auto& r_origin_force_variable =
                    GetComponentVariable(rOriginVariablesForces, component_index);
                const auto& r_origin_moment_variable =
                    GetComponentVariable(rOriginVariablesMoments, component_index);

                r_node.FastGetSolutionStepValue(r_origin_force_variable) +=
                    nodal_force_global(component_index);
                r_node.FastGetSolutionStepValue(r_origin_moment_variable) +=
                    nodal_moment_global(component_index);
            }
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
double BeamSplineMapper<TSparseSpace, TDenseSpace>::EvaluateKernelValue(const double Distance) const
{
    return std::abs(Distance) * Distance * Distance;
}

template<class TSparseSpace, class TDenseSpace>
double BeamSplineMapper<TSparseSpace, TDenseSpace>::EvaluateKernelFirstDerivative(const double Distance) const
{
    return 3.0 * Distance * std::abs(Distance);
}

template<class TSparseSpace, class TDenseSpace>
double BeamSplineMapper<TSparseSpace, TDenseSpace>::EvaluateKernelSecondDerivative(const double Distance) const
{
    if (std::abs(Distance) < PolynomialTolerance) {
        return 0.0;
    }
    return 6.0 * std::abs(Distance);
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
    const IndexType number_of_source_nodes = rSourceCoordinates.size();

    rMss.resize(number_of_source_nodes, number_of_source_nodes, false);
    rMssFirstDerivative.resize(number_of_source_nodes, number_of_source_nodes, false);
    rMssSecondDerivative.resize(number_of_source_nodes, number_of_source_nodes, false);
    rPs.resize(number_of_source_nodes, 2, false);
    rPsFirstDerivative.resize(number_of_source_nodes, 2, false);

    for (IndexType i = 0; i < number_of_source_nodes; ++i) {
        rPs(i, 0) = 1.0;
        rPs(i, 1) = rSourceCoordinates[i];
        rPsFirstDerivative(i, 0) = 0.0;
        rPsFirstDerivative(i, 1) = 1.0;

        for (IndexType j = 0; j < number_of_source_nodes; ++j) {
            const double distance = rSourceCoordinates[i] - rSourceCoordinates[j];
            rMss(i, j) = EvaluateKernelValue(distance);
            rMssFirstDerivative(i, j) = EvaluateKernelFirstDerivative(distance);
            rMssSecondDerivative(i, j) = EvaluateKernelSecondDerivative(distance);
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildEvaluationRow(
    const std::vector<double>& rSourceCoordinates,
    const double ProjectionCoordinate) const
{
    const IndexType number_of_source_nodes = rSourceCoordinates.size();
    VectorType evaluation_row(2 * number_of_source_nodes + 2, 0.0);

    for (IndexType j = 0; j < number_of_source_nodes; ++j) {
        const double distance = ProjectionCoordinate - rSourceCoordinates[j];
        evaluation_row(j) = EvaluateKernelValue(distance);
        evaluation_row(number_of_source_nodes + j) = EvaluateKernelFirstDerivative(distance);
    }
    evaluation_row(2 * number_of_source_nodes + 0) = 1.0;
    evaluation_row(2 * number_of_source_nodes + 1) = ProjectionCoordinate;

    return evaluation_row;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildEvaluationDerivativeRow(
    const std::vector<double>& rSourceCoordinates,
    const double ProjectionCoordinate,
    const double CoordinateHalfLength) const
{
    const IndexType number_of_source_nodes = rSourceCoordinates.size();
    VectorType evaluation_derivative_row(2 * number_of_source_nodes + 2, 0.0);

    for (IndexType j = 0; j < number_of_source_nodes; ++j) {
        const double distance = ProjectionCoordinate - rSourceCoordinates[j];
        evaluation_derivative_row(j) = EvaluateKernelFirstDerivative(distance) / CoordinateHalfLength;
        evaluation_derivative_row(number_of_source_nodes + j) =
            EvaluateKernelSecondDerivative(distance) / CoordinateHalfLength;
    }
    evaluation_derivative_row(2 * number_of_source_nodes + 1) = 1.0 / CoordinateHalfLength;

    return evaluation_derivative_row;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildRightHandSide(
    const std::vector<double>& rDisplacements,
    const std::vector<double>& rRotations,
    const double CoordinateHalfLength) const
{
    const IndexType number_of_source_nodes = rDisplacements.size();
    VectorType right_hand_side(2 * number_of_source_nodes + 2, 0.0);

    for (IndexType i = 0; i < number_of_source_nodes; ++i) {
        right_hand_side(i) = rDisplacements[i];
        right_hand_side(number_of_source_nodes + i) = CoordinateHalfLength * rRotations[i];
    }

    return right_hand_side;
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildLocalSourceData(
    const BeamChainCacheData& rBeamChainCacheData,
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
    typename BeamSplineMapper<TSparseSpace, TDenseSpace>::ComponentArrayType& rLocalDisplacements,
    typename BeamSplineMapper<TSparseSpace, TDenseSpace>::ComponentArrayType& rLocalRotations) const
{
    const IndexType number_of_source_nodes = rBeamChainCacheData.SupportNodeIds.size();
    for (IndexType component_index = 0; component_index < 3; ++component_index) {
        rLocalDisplacements[component_index].assign(number_of_source_nodes, 0.0);
        rLocalRotations[component_index].assign(number_of_source_nodes, 0.0);
    }

    for (IndexType i = 0; i < number_of_source_nodes; ++i) {
        const auto& r_node = mrModelPartOrigin.GetNode(rBeamChainCacheData.SupportNodeIds[i]);
        array_1d<double, 3> displacement_global = ZeroVector(3);
        array_1d<double, 3> rotation_global = ZeroVector(3);
        for (IndexType j = 0; j < 3; ++j) {
            displacement_global[j] =
                r_node.FastGetSolutionStepValue(GetComponentVariable(rOriginVariablesDisplacements, j));
            rotation_global[j] =
                r_node.FastGetSolutionStepValue(GetComponentVariable(rOriginVariablesRotations, j));
        }

        const MatrixType rotation_matrix_global_to_local =
            TransposeRotationMatrix(rBeamChainCacheData.SupportFramesLocalToGlobal[i]);
        const VectorType displacement_local = TransformVectorToLocal(rotation_matrix_global_to_local, displacement_global);
        const VectorType rotation_local = TransformVectorToLocal(rotation_matrix_global_to_local, rotation_global);

        for (IndexType j = 0; j < 3; ++j) {
            rLocalDisplacements[j][i] = displacement_local(j);
            rLocalRotations[j][i] = rotation_local(j);
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::PrepareBeamChainCacheData()
{
    for (auto& r_local_sys : mMapperLocalSystems) {
        if (!r_local_sys->HasInterfaceInfo()) {
            continue;
        }

        auto beam_sys = dynamic_cast<BeamSplineMapperLocalSystem*>(r_local_sys.get());
        KRATOS_ERROR_IF_NOT(beam_sys) << "Expected BeamSplineMapperLocalSystem!" << std::endl;

        GetOrCreateBeamChainCacheData(beam_sys->CalculateProjectionData().BeamGeometry);
    }
}

template<class TSparseSpace, class TDenseSpace>
std::unordered_map<std::string, typename BeamSplineMapper<TSparseSpace, TDenseSpace>::BeamChainSourceStateData>
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildAllBeamChainSourceStates(
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rOriginVariablesRotations) const
{
    std::unordered_map<std::string, BeamChainSourceStateData> beam_chain_source_states;
    beam_chain_source_states.reserve(mBeamChainCache.size());

    for (const auto& r_key_cache_pair : mBeamChainCache) {
        const auto& r_beam_chain_cache_data = r_key_cache_pair.second;

        BeamChainSourceStateData source_state_data;
        BuildLocalSourceData(
            r_beam_chain_cache_data,
            rOriginVariablesDisplacements,
            rOriginVariablesRotations,
            source_state_data.LocalDisplacements,
            source_state_data.LocalRotations);

        const VectorType right_hand_side_y = BuildRightHandSide(
            source_state_data.LocalDisplacements[1],
            source_state_data.LocalRotations[2],
            r_beam_chain_cache_data.CoordinateHalfLength);

        std::vector<double> negative_local_rotation_y(
            source_state_data.LocalRotations[1].size(),
            0.0);
        for (IndexType i = 0; i < source_state_data.LocalRotations[1].size(); ++i) {
            negative_local_rotation_y[i] = -source_state_data.LocalRotations[1][i];
        }

        const VectorType right_hand_side_z = BuildRightHandSide(
            source_state_data.LocalDisplacements[2],
            negative_local_rotation_y,
            r_beam_chain_cache_data.CoordinateHalfLength);

        source_state_data.SplineCoefficientsY = SolveSplineCoefficients(
            r_beam_chain_cache_data.SplineSystemMatrix,
            right_hand_side_y);
        source_state_data.SplineCoefficientsZ = SolveSplineCoefficients(
            r_beam_chain_cache_data.SplineSystemMatrix,
            right_hand_side_z);

        beam_chain_source_states.emplace(r_beam_chain_cache_data.Key, std::move(source_state_data));
    }

    return beam_chain_source_states;
}

template<class TSparseSpace, class TDenseSpace>
const typename BeamSplineMapper<TSparseSpace, TDenseSpace>::BeamChainCacheData&
BeamSplineMapper<TSparseSpace, TDenseSpace>::GetOrCreateBeamChainCacheData(
    const GeometryType& rBeamGeometry)
{
    const IndexType first_node_id = rBeamGeometry[0].Id();
    const auto key_iterator = mNodeIdToBeamChainKey.find(first_node_id);
    if (key_iterator != mNodeIdToBeamChainKey.end()) {
        return mBeamChainCache.at(key_iterator->second);
    }

    BeamChainCacheData beam_chain_cache_data;
    ComputeBeamChainSupport(
        rBeamGeometry,
        beam_chain_cache_data.SupportNodeIds,
        beam_chain_cache_data.SupportNodeIdToLocalIndex);
    ComputeSupportReferenceData(
        beam_chain_cache_data.SupportNodeIds,
        beam_chain_cache_data.LocalXCoordinates,
        beam_chain_cache_data.CoordinateHalfLength,
        beam_chain_cache_data.SupportFramesLocalToGlobal,
        beam_chain_cache_data.SupportReferenceCoordinates);
    beam_chain_cache_data.Key = CreateBeamChainKey(beam_chain_cache_data.SupportNodeIds);
    beam_chain_cache_data.SplineSystemMatrix = BuildSplineSystemMatrix(
        beam_chain_cache_data.LocalXCoordinates);
    beam_chain_cache_data.TransposedSplineSystemMatrix =
        BuildTransposedMatrix(beam_chain_cache_data.SplineSystemMatrix);

    for (const IndexType node_id : beam_chain_cache_data.SupportNodeIds) {
        mNodeIdToBeamChainKey[node_id] = beam_chain_cache_data.Key;
    }

    return mBeamChainCache.emplace(
        beam_chain_cache_data.Key,
        std::move(beam_chain_cache_data)).first->second;
}

template<class TSparseSpace, class TDenseSpace>
std::string BeamSplineMapper<TSparseSpace, TDenseSpace>::CreateBeamChainKey(
    const std::vector<IndexType>& rSupportNodeIds) const
{
    std::ostringstream key_stream;
    for (const IndexType node_id : rSupportNodeIds) {
        key_stream << node_id << '-';
    }
    return key_stream.str();
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::MatrixType
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildSplineSystemMatrix(
    const std::vector<double>& rSourceCoordinates) const
{
    MatrixType mss;
    MatrixType mss_first_derivative;
    MatrixType mss_second_derivative;
    MatrixType ps;
    MatrixType ps_first_derivative;

    BuildSplineSystemMatrices(
        rSourceCoordinates,
        mss,
        mss_first_derivative,
        mss_second_derivative,
        ps,
        ps_first_derivative);

    const IndexType number_of_source_nodes = rSourceCoordinates.size();
    const IndexType system_size = 2 * number_of_source_nodes + 2;
    MatrixType spline_system(system_size, system_size, 0.0);

    for (IndexType i = 0; i < number_of_source_nodes; ++i) {
        for (IndexType j = 0; j < number_of_source_nodes; ++j) {
            spline_system(i, j) = mss(i, j);
            spline_system(i, number_of_source_nodes + j) = mss_first_derivative(i, j);
            spline_system(number_of_source_nodes + i, j) = mss_first_derivative(i, j);
            spline_system(number_of_source_nodes + i, number_of_source_nodes + j) =
                mss_second_derivative(i, j);
        }

        spline_system(i, 2 * number_of_source_nodes + 0) = ps(i, 0);
        spline_system(i, 2 * number_of_source_nodes + 1) = ps(i, 1);
        spline_system(number_of_source_nodes + i, 2 * number_of_source_nodes + 0) =
            ps_first_derivative(i, 0);
        spline_system(number_of_source_nodes + i, 2 * number_of_source_nodes + 1) =
            ps_first_derivative(i, 1);
        spline_system(2 * number_of_source_nodes + 0, i) = ps(i, 0);
        spline_system(2 * number_of_source_nodes + 1, i) = ps(i, 1);
        spline_system(2 * number_of_source_nodes + 0, number_of_source_nodes + i) =
            ps_first_derivative(i, 0);
        spline_system(2 * number_of_source_nodes + 1, number_of_source_nodes + i) =
            ps_first_derivative(i, 1);
    }

    return spline_system;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::MatrixType
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildTransposedMatrix(
    const MatrixType& rMatrix) const
{
    const IndexType matrix_size = rMatrix.size1();
    KRATOS_ERROR_IF(matrix_size == 0)
        << "Cannot transpose an empty beam-spline matrix." << std::endl;
    KRATOS_ERROR_IF(rMatrix.size2() != matrix_size)
        << "Cannot transpose a non-square beam-spline matrix." << std::endl;

    MatrixType transposed_matrix(matrix_size, matrix_size, 0.0);
    for (IndexType i = 0; i < matrix_size; ++i) {
        for (IndexType j = 0; j < matrix_size; ++j) {
            transposed_matrix(i, j) = rMatrix(j, i);
        }
    }

    return transposed_matrix;
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::ComputeSupportReferenceData(
    const std::vector<IndexType>& rSupportNodeIds,
    std::vector<double>& rLocalXCoordinates,
    double& rCoordinateHalfLength,
    std::vector<MatrixType>& rSupportFramesLocalToGlobal,
    std::vector<array_1d<double, 3>>& rSupportReferenceCoordinates) const
{
    const IndexType number_of_source_nodes = rSupportNodeIds.size();
    KRATOS_ERROR_IF(number_of_source_nodes < 2)
        << "BeamSplineMapper requires at least two support nodes." << std::endl;

    rLocalXCoordinates.assign(number_of_source_nodes, 0.0);
    rSupportFramesLocalToGlobal.assign(number_of_source_nodes, MatrixType(3, 3, 0.0));
    rSupportReferenceCoordinates.assign(number_of_source_nodes, ZeroVector(3));

    for (IndexType i = 0; i < number_of_source_nodes; ++i) {
        const auto& r_node = mrModelPartOrigin.GetNode(rSupportNodeIds[i]);
        rSupportReferenceCoordinates[i] = GetReferenceCoordinates(r_node);
        if (i > 0) {
            rLocalXCoordinates[i] = rLocalXCoordinates[i - 1] +
                norm_2(rSupportReferenceCoordinates[i] - rSupportReferenceCoordinates[i - 1]);
        }
    }

    const double total_length = rLocalXCoordinates.back();
    KRATOS_ERROR_IF(total_length <= 2.0 * FrameTolerance)
        << "BeamSplineMapper cannot normalize a zero-length beam chain." << std::endl;
    rCoordinateHalfLength = 0.5 * total_length;
    for (double& r_coordinate : rLocalXCoordinates) {
        r_coordinate = r_coordinate / rCoordinateHalfLength - 1.0;
    }

    std::vector<array_1d<double, 3>> tangents(number_of_source_nodes, ZeroVector(3));
    tangents.front() = NormalizeVector(
        rSupportReferenceCoordinates[1] - rSupportReferenceCoordinates[0]);
    tangents.back() = NormalizeVector(
        rSupportReferenceCoordinates[number_of_source_nodes - 1] -
        rSupportReferenceCoordinates[number_of_source_nodes - 2]);

    for (IndexType i = 1; i + 1 < number_of_source_nodes; ++i) {
        tangents[i] = NormalizeVector(
            rSupportReferenceCoordinates[i + 1] - rSupportReferenceCoordinates[i - 1]);
    }

    rSupportFramesLocalToGlobal[0] = ComputeInitialFrameLocalToGlobal(tangents[0]);

    for (IndexType i = 1; i < number_of_source_nodes; ++i) {
        array_1d<double, 3> previous_normal = ZeroVector(3);
        for (IndexType k = 0; k < 3; ++k) {
            previous_normal[k] = rSupportFramesLocalToGlobal[i - 1](k, 1);
        }

        array_1d<double, 3> current_normal = previous_normal;
        array_1d<double, 3> transport_axis = CrossProduct(tangents[i - 1], tangents[i]);
        const double axis_norm = norm_2(transport_axis);
        const double tangent_dot = std::max(-1.0, std::min(1.0, DotProduct(tangents[i - 1], tangents[i])));

        if (axis_norm > FrameTolerance) {
            transport_axis /= axis_norm;
            const double angle = std::atan2(axis_norm, tangent_dot);
            current_normal = RotateVectorAroundAxis(previous_normal, transport_axis, angle);
        }

        current_normal -= DotProduct(current_normal, tangents[i]) * tangents[i];
        if (norm_2(current_normal) < FrameTolerance) {
            const MatrixType fallback_frame = ComputeInitialFrameLocalToGlobal(tangents[i]);
            for (IndexType k = 0; k < 3; ++k) {
                current_normal[k] = fallback_frame(k, 1);
            }
        }
        current_normal = NormalizeVector(current_normal);
        array_1d<double, 3> current_binormal = NormalizeVector(CrossProduct(tangents[i], current_normal));
        current_normal = NormalizeVector(CrossProduct(current_binormal, tangents[i]));

        for (IndexType k = 0; k < 3; ++k) {
            rSupportFramesLocalToGlobal[i](k, 0) = tangents[i][k];
            rSupportFramesLocalToGlobal[i](k, 1) = current_normal[k];
            rSupportFramesLocalToGlobal[i](k, 2) = current_binormal[k];
        }
    }
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
    const IndexType first_segment_node_index =
        rSupportNodeIdToLocalIndex.at(rBeamGeometry[0].Id());
    const IndexType second_segment_node_index =
        rSupportNodeIdToLocalIndex.at(rBeamGeometry[1].Id());

    array_1d<double, 3> tangent = NormalizeVector(
        rSupportReferenceCoordinates[second_segment_node_index] -
        rSupportReferenceCoordinates[first_segment_node_index]);

    array_1d<double, 3> normal = ZeroVector(3);
    for (IndexType k = 0; k < 3; ++k) {
        normal[k] =
            rLinearShapeValues(0) * rSupportFramesLocalToGlobal[first_segment_node_index](k, 1) +
            rLinearShapeValues(1) * rSupportFramesLocalToGlobal[second_segment_node_index](k, 1);
    }

    normal -= DotProduct(normal, tangent) * tangent;
    if (norm_2(normal) < FrameTolerance) {
        for (IndexType k = 0; k < 3; ++k) {
            normal[k] = rSupportFramesLocalToGlobal[first_segment_node_index](k, 1);
        }
        normal -= DotProduct(normal, tangent) * tangent;
    }
    normal = NormalizeVector(normal);

    array_1d<double, 3> binormal = NormalizeVector(CrossProduct(tangent, normal));
    normal = NormalizeVector(CrossProduct(binormal, tangent));

    MatrixType evaluation_frame_local_to_global(3, 3, 0.0);
    for (IndexType k = 0; k < 3; ++k) {
        evaluation_frame_local_to_global(k, 0) = tangent[k];
        evaluation_frame_local_to_global(k, 1) = normal[k];
        evaluation_frame_local_to_global(k, 2) = binormal[k];
    }

    return evaluation_frame_local_to_global;
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::ComputeBeamChainSupport(
    const GeometryType& rBeamGeometry,
    std::vector<IndexType>& rSupportNodeIds,
    std::unordered_map<IndexType, IndexType>& rSupportNodeIdToLocalIndex) const
{
    std::unordered_map<IndexType, std::unordered_set<IndexType>> adjacency;
    AddBeamConnectivity(mrModelPartOrigin.Elements(), adjacency);
    AddBeamConnectivity(mrModelPartOrigin.Conditions(), adjacency);

    const IndexType first_geometry_node_id = rBeamGeometry[0].Id();
    const IndexType second_geometry_node_id = rBeamGeometry[1].Id();

    KRATOS_ERROR_IF(adjacency.find(first_geometry_node_id) == adjacency.end())
        << "Beam node " << first_geometry_node_id << " is not connected to any beam support entity." << std::endl;
    KRATOS_ERROR_IF(adjacency.find(second_geometry_node_id) == adjacency.end())
        << "Beam node " << second_geometry_node_id << " is not connected to any beam support entity." << std::endl;

    std::unordered_set<IndexType> support_node_id_set;
    std::queue<IndexType> pending_nodes;
    pending_nodes.push(first_geometry_node_id);
    support_node_id_set.insert(first_geometry_node_id);

    while (!pending_nodes.empty()) {
        const IndexType current_node_id = pending_nodes.front();
        pending_nodes.pop();

        for (const IndexType neighbour_node_id : adjacency[current_node_id]) {
            if (support_node_id_set.insert(neighbour_node_id).second) {
                pending_nodes.push(neighbour_node_id);
            }
        }
    }

    KRATOS_ERROR_IF(support_node_id_set.find(second_geometry_node_id) == support_node_id_set.end())
        << "Beam geometry nodes " << first_geometry_node_id << " and " << second_geometry_node_id
        << " do not belong to the same beam chain." << std::endl;

    IndexType start_node_id = 0;
    for (const IndexType node_id : support_node_id_set) {
        IndexType local_degree = 0;
        for (const IndexType neighbour_node_id : adjacency[node_id]) {
            if (support_node_id_set.find(neighbour_node_id) != support_node_id_set.end()) {
                ++local_degree;
            }
        }

        KRATOS_ERROR_IF(local_degree > 2)
            << "BeamSplineMapper only supports non-branching beam chains." << std::endl;

        if (local_degree == 1) {
            start_node_id = node_id;
        }
    }

    KRATOS_ERROR_IF(start_node_id == 0)
        << "BeamSplineMapper requires an open beam chain with two end nodes." << std::endl;

    rSupportNodeIds.clear();
    rSupportNodeIds.reserve(support_node_id_set.size());

    IndexType previous_node_id = 0;
    IndexType current_node_id = start_node_id;
    while (true) {
        rSupportNodeIds.push_back(current_node_id);

        IndexType next_node_id = 0;
        for (const IndexType neighbour_node_id : adjacency[current_node_id]) {
            if (support_node_id_set.find(neighbour_node_id) == support_node_id_set.end()) {
                continue;
            }

            if (neighbour_node_id != previous_node_id) {
                next_node_id = neighbour_node_id;
                break;
            }
        }

        if (next_node_id == 0) {
            break;
        }

        previous_node_id = current_node_id;
        current_node_id = next_node_id;
    }

    KRATOS_ERROR_IF(rSupportNodeIds.size() != support_node_id_set.size())
        << "Failed to order all nodes in the current beam chain." << std::endl;

    auto first_geometry_position = std::find(
        rSupportNodeIds.begin(),
        rSupportNodeIds.end(),
        first_geometry_node_id);
    auto second_geometry_position = std::find(
        rSupportNodeIds.begin(),
        rSupportNodeIds.end(),
        second_geometry_node_id);

    KRATOS_ERROR_IF(
        first_geometry_position == rSupportNodeIds.end() ||
        second_geometry_position == rSupportNodeIds.end())
        << "Failed to recover the current beam geometry inside the ordered support." << std::endl;

    if (first_geometry_position > second_geometry_position) {
        std::reverse(rSupportNodeIds.begin(), rSupportNodeIds.end());
    }

    rSupportNodeIdToLocalIndex.clear();
    for (IndexType i = 0; i < rSupportNodeIds.size(); ++i) {
        rSupportNodeIdToLocalIndex[rSupportNodeIds[i]] = i;
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapper<TSparseSpace, TDenseSpace>::InitializeOriginForcesAndMoments(
    const Variable<array_1d<double, 3>>& rOriginVariablesForces,
    const Variable<array_1d<double, 3>>& rOriginVariablesMoments)
{
    for (auto& r_node : mrModelPartOrigin.Nodes()) {
        for (IndexType i = 0; i < 3; ++i) {
            r_node.FastGetSolutionStepValue(GetComponentVariable(rOriginVariablesForces, i)) = 0.0;
            r_node.FastGetSolutionStepValue(GetComponentVariable(rOriginVariablesMoments, i)) = 0.0;
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::SolveSplineCoefficients(
    const MatrixType& rSplineSystemMatrix,
    const VectorType& rRightHandSide) const
{
    const IndexType system_size = rSplineSystemMatrix.size1();
    KRATOS_ERROR_IF(system_size == 0)
        << "Cannot solve an empty spline coefficient system." << std::endl;
    KRATOS_ERROR_IF(rSplineSystemMatrix.size2() != system_size)
        << "Spline coefficient system matrix must be square." << std::endl;
    KRATOS_ERROR_IF(rRightHandSide.size() != system_size)
        << "Spline coefficient right-hand side has the wrong size." << std::endl;

    VectorType coefficients(system_size, 0.0);
    MathUtils<double>::Solve(rSplineSystemMatrix, coefficients, rRightHandSide);

    return coefficients;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::TransformGlobalToLocal(
    const MatrixType& rRotationMatrixGlobalToLocal,
    const array_1d<double, 3>& rReferencePoint,
    const array_1d<double, 3>& rCoordinates) const
{
    VectorType global_relative_vector(3);
    global_relative_vector(0) = rCoordinates[0] - rReferencePoint[0];
    global_relative_vector(1) = rCoordinates[1] - rReferencePoint[1];
    global_relative_vector(2) = rCoordinates[2] - rReferencePoint[2];

    VectorType local_vector(3, 0.0);
    TDenseSpace::Mult(rRotationMatrixGlobalToLocal, global_relative_vector, local_vector);
    return local_vector;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::TransformVectorToLocal(
    const MatrixType& rRotationMatrixGlobalToLocal,
    const array_1d<double, 3>& rVector) const
{
    VectorType local_vector(3, 0.0);
    TDenseSpace::Mult(rRotationMatrixGlobalToLocal, rVector, local_vector);
    return local_vector;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::TransformVectorToGlobal(
    const MatrixType& rRotationMatrixLocalToGlobal,
    const VectorType& rLocalVector) const
{
    VectorType global_vector(3, 0.0);
    TDenseSpace::Mult(rRotationMatrixLocalToGlobal, rLocalVector, global_vector);
    return global_vector;
}

template<class TSparseSpace, class TDenseSpace>
array_1d<double, 3> BeamSplineMapper<TSparseSpace, TDenseSpace>::GetReferenceCoordinates(
    const Node& rNode) const
{
    array_1d<double, 3> reference_coordinates = ZeroVector(3);
    reference_coordinates[0] = rNode.X0();
    reference_coordinates[1] = rNode.Y0();
    reference_coordinates[2] = rNode.Z0();
    return reference_coordinates;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::MatrixType
BeamSplineMapper<TSparseSpace, TDenseSpace>::ComputeInitialFrameLocalToGlobal(
    const array_1d<double, 3>& rTangent) const
{
    array_1d<double, 3> reference_axis = ZeroVector(3);
    reference_axis[2] = 1.0;
    if (std::abs(DotProduct(rTangent, reference_axis)) > 0.9) {
        reference_axis = ZeroVector(3);
        reference_axis[1] = 1.0;
    }

    array_1d<double, 3> normal = CrossProduct(reference_axis, rTangent);
    if (norm_2(normal) < FrameTolerance) {
        reference_axis = ZeroVector(3);
        reference_axis[0] = 1.0;
        normal = CrossProduct(reference_axis, rTangent);
    }
    normal = NormalizeVector(normal);
    array_1d<double, 3> binormal = NormalizeVector(CrossProduct(rTangent, normal));
    normal = NormalizeVector(CrossProduct(binormal, rTangent));

    MatrixType frame_local_to_global(3, 3, 0.0);
    for (IndexType i = 0; i < 3; ++i) {
        frame_local_to_global(i, 0) = rTangent[i];
        frame_local_to_global(i, 1) = normal[i];
        frame_local_to_global(i, 2) = binormal[i];
    }
    return frame_local_to_global;
}

template<class TSparseSpace, class TDenseSpace>
array_1d<double, 3> BeamSplineMapper<TSparseSpace, TDenseSpace>::RotateVectorAroundAxis(
    const array_1d<double, 3>& rVector,
    const array_1d<double, 3>& rAxis,
    const double Angle) const
{
    const double cosine = std::cos(Angle);
    const double sine = std::sin(Angle);

    array_1d<double, 3> rotated_vector = cosine * rVector;
    rotated_vector += sine * CrossProduct(rAxis, rVector);
    rotated_vector += (1.0 - cosine) * DotProduct(rAxis, rVector) * rAxis;
    return rotated_vector;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::EvaluatePointDisplacementLocal(
    const VectorType& rEvaluationRow,
    const VectorType& rEvaluationDerivativeRow,
    const VectorType& rSplineCoefficientsY,
    const VectorType& rSplineCoefficientsZ,
    const double AxialDisplacement,
    const double TorsionalRotation,
    const VectorType& rOffsetVectorLocal) const
{
    const double transverse_displacement_y = inner_prod(rEvaluationRow, rSplineCoefficientsY);
    const double transverse_displacement_z = inner_prod(rEvaluationRow, rSplineCoefficientsZ);
    const double transverse_slope_y = inner_prod(rEvaluationDerivativeRow, rSplineCoefficientsY);
    const double transverse_slope_z = inner_prod(rEvaluationDerivativeRow, rSplineCoefficientsZ);

    VectorType centerline_displacement(3);
    centerline_displacement(0) = AxialDisplacement;
    centerline_displacement(1) = transverse_displacement_y;
    centerline_displacement(2) = transverse_displacement_z;

    const double theta_x = TorsionalRotation;
    const double theta_y = -transverse_slope_z;
    const double theta_z = transverse_slope_y;

    const double cos_x = std::cos(theta_x);
    const double sin_x = std::sin(theta_x);
    const double cos_y = std::cos(theta_y);
    const double sin_y = std::sin(theta_y);
    const double cos_z = std::cos(theta_z);
    const double sin_z = std::sin(theta_z);

    MatrixType rotation_x(3, 3, 0.0);
    rotation_x(0, 0) = 1.0;
    rotation_x(1, 1) = cos_x;
    rotation_x(1, 2) = -sin_x;
    rotation_x(2, 1) = sin_x;
    rotation_x(2, 2) = cos_x;

    MatrixType rotation_y(3, 3, 0.0);
    rotation_y(0, 0) = cos_y;
    rotation_y(0, 2) = sin_y;
    rotation_y(1, 1) = 1.0;
    rotation_y(2, 0) = -sin_y;
    rotation_y(2, 2) = cos_y;

    MatrixType rotation_z(3, 3, 0.0);
    rotation_z(0, 0) = cos_z;
    rotation_z(0, 1) = -sin_z;
    rotation_z(1, 0) = sin_z;
    rotation_z(1, 1) = cos_z;
    rotation_z(2, 2) = 1.0;

    const MatrixType rotation_x_z = prod(rotation_x, rotation_z);
    const MatrixType rotation_matrix = prod(rotation_x_z, rotation_y);
    VectorType rotated_offset(3, 0.0);
    TDenseSpace::Mult(rotation_matrix, rOffsetVectorLocal, rotated_offset);

    VectorType local_displacement(3);
    for (IndexType i = 0; i < 3; ++i) {
        local_displacement(i) =
            centerline_displacement(i) + rotated_offset(i) - rOffsetVectorLocal(i);
    }

    return local_displacement;
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
