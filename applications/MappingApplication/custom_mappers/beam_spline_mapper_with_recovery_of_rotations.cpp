//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
// Main authors:    Qinfei Ran - Initial development
// Contributor:

/*
Compilation Commands:

make -C build/applications/MappingApplication custom_mappers/beam_spline_mapper.o
make -f applications/MappingApplication/CMakeFiles/KratosMappingCore.dir/build.make applications/MappingApplication/libKratosMappingCore.so
make -f applications/MappingApplication/CMakeFiles/KratosMappingApplication.dir/build.make applications/MappingApplication/KratosMappingApplication.cpython-310-x86_64-linux-gnu.so


*/



// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

// External includes
#include "utilities/math_utils.h"
#include "geometries/line_3d_2.h"

// Project includes
#include "beam_spline_mapper_with_recovery_of_rotations.h"
#include "factories/linear_solver_factory.h"
#include "mappers/mapper_define.h"
#include "mapping_application_variables.h"

namespace Kratos
{

namespace
{
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

const Variable<double>& GetComponentVariable(
    const Variable<array_1d<double, 3>>& rVariable,
    const IndexType ComponentIndex)
{
    return KratosComponents<Variable<double>>::Get(rVariable.Name() + ComponentSuffixes[ComponentIndex]);
}

}

void BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo::ProcessSearchResult(const InterfaceObject& rInterfaceObject)
{
    SaveSearchResult(rInterfaceObject, true);
}

void BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo::ProcessSearchResultForApproximation(const InterfaceObject& rInterfaceObject)
{
    const auto p_geom = rInterfaceObject.pGetBaseGeometry();

    double proj_dist;

    const Point point_to_proj(this->Coordinates());
    Point projected_point;
    Point local_coords;

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

    p_geom->IsInside(projected_point, local_coords, 1e-14);
    p_geom->PointLocalCoordinates(local_coords, projected_point);

    // Approximate projections can fall outside the beam segment; clamp them to the valid line element range.
    local_coords[0] = std::clamp(local_coords[0], -1.0, 1.0);
    p_geom->ShapeFunctionsValues(linear_shape_function_values, local_coords);

    Point clamped_projected_point;
    noalias(clamped_projected_point.Coordinates()) =
        linear_shape_function_values[0] * (*p_geom)[0].Coordinates() +
        linear_shape_function_values[1] * (*p_geom)[1].Coordinates();

    const double projection_distance = norm_2(point_to_proj - clamped_projected_point);

    if (projection_distance < mClosestProjectionDistance) {
        SetIsApproximation();

        mPairingIndex = ProjectionUtilities::PairingIndex::Closest_Point;

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

        mClosestProjectionDistance = projection_distance;
        mNodeIds = eq_ids;

        const IndexType num_values_linear = linear_shape_function_values.size();
        if (mLinearShapeFunctionValues.size() != num_values_linear) {
            mLinearShapeFunctionValues.resize(num_values_linear);
        }

        for (IndexType i = 0; i < num_values_linear; ++i) {
            mLinearShapeFunctionValues[i] = linear_shape_function_values[i];
        }

        mProjectionOfPoint = clamped_projected_point;
        mpInterfaceObject = make_shared<InterfaceGeometryObject>(rInterfaceObject.pGetBaseGeometry());
    }

}

void BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo::SaveSearchResult(
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

void BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo::ComputeRotationMatrix()
{
    // Known limitations of this frame construction:
    // - It compares normalized floating-point axis components with exact 0.0/1.0 values, so nearly axis-aligned
    //   beams can enter an unintended branch due to round-off.
    // - Negative axis-aligned beams, e.g. (-1, 0, 0), can fall through to the generic branch and divide by
    //   axis_X[2] when it is zero.
    // - The special case axis_X = (0, 0, 1) sets axis_Z = (1, 0, 0), while axis_X x axis_Y = (-1, 0, 0);
    //   this creates a left-handed local frame.
    // - A more robust construction would choose a non-parallel reference vector and build the frame with
    //   normalized cross products.
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

    mRotationMatrix_L_G = rotation_matrix;

}

void BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MapperInterfaceInfo);
    rSerializer.save("NodeIds", mNodeIds);
    rSerializer.save("LinearShapeFunctionValues", mLinearShapeFunctionValues);
    rSerializer.save("ProjectionOfPoint", mProjectionOfPoint);
    rSerializer.save("ClosestProjectionDistance", mClosestProjectionDistance);
    rSerializer.save("PairingIndex", static_cast<int>(mPairingIndex));
}

void BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo::load(Serializer& rSerializer)
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

void BeamSplineMapperWithRecoveryOfRotationsLocalSystem::PairingInfo(std::ostream& rOStream, const int EchoLevel) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not initialized!" << std::endl;
    rOStream << "BeamSplineMapperWithRecoveryOfRotationsLocalSystem based on " << mpNode->Info();
}

BeamSplineMapperWithRecoveryOfRotationsLocalSystem::ProjectionData BeamSplineMapperWithRecoveryOfRotationsLocalSystem::CalculateProjectionData()
{
    ProjectionData projection_data;
    projection_data.pNode = mpNode;

    KRATOS_ERROR_IF_NOT(mpNode) << "Destination node is a nullptr." << std::endl;
    KRATOS_ERROR_IF(mInterfaceInfos.empty())
        << "Cannot calculate BeamSpline projection data without interface information." << std::endl;

    for (auto& r_interface_info : mInterfaceInfos) {
        auto beam_interface_info = std::dynamic_pointer_cast<BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo>(r_interface_info);
        KRATOS_ERROR_IF_NOT(beam_interface_info)
            << "Expected BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo in mInterfaceInfos but got nullptr or wrong type." << std::endl;

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
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::BeamSplineMapperWithRecoveryOfRotations(
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
    mKernelType = JsonParameters["kernel_type"].GetString();
    mKernelRadius = JsonParameters["kernel_radius"].GetDouble();
    mRegularization = JsonParameters["regularization"].GetDouble();
    mPolynomialBasis = JsonParameters["polynomial_basis"].GetString();

    std::transform(mKernelType.begin(), mKernelType.end(), mKernelType.begin(), ::tolower);
    std::transform(mPolynomialBasis.begin(), mPolynomialBasis.end(), mPolynomialBasis.begin(), ::tolower);
    if (mPolynomialBasis == "auto") {
        mPolynomialBasis = mKernelType == "cubic" ? "line_adapted" : "full_3d";
    }

    KRATOS_ERROR_IF(
        mKernelType != "cubic" &&
        mKernelType != "inverse_multiquadric" &&
        mKernelType != "multiquadric" &&
        mKernelType != "gaussian" &&
        mKernelType != "thin_plate_spline" &&
        mKernelType != "wendland_c2")
        << "BeamSplineMapperWithRecoveryOfRotations supports the following kernels: "
        << "'cubic', 'inverse_multiquadric', 'multiquadric', 'gaussian', "
        << "'thin_plate_spline', 'wendland_c2'." << std::endl;
    KRATOS_ERROR_IF(mKernelType != "cubic" && mKernelType != "thin_plate_spline" && mKernelRadius <= 0.0)
        << "The kernel_radius must be positive for kernel_type '" << mKernelType << "'." << std::endl;
    KRATOS_ERROR_IF(mRegularization < 0.0)
        << "The regularization must be non-negative." << std::endl;
    KRATOS_ERROR_IF(
        mPolynomialBasis != "full_3d" &&
        mPolynomialBasis != "line_adapted" &&
        mPolynomialBasis != "enriched_line_adapted" &&
        mPolynomialBasis != "high_order_line_adapted")
        << "BeamSplineMapperWithRecoveryOfRotations supports the following polynomial bases: "
        << "'full_3d', 'line_adapted', 'enriched_line_adapted', "
        << "'high_order_line_adapted' and 'auto'." << std::endl;

    Initialize();

    KRATOS_CATCH("")
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::UpdateInterface(
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
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::Map(
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
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::Map(
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
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::InverseMap(
    const Variable<array_1d<double, 3>>& rOriginVariablesForces,
    const Variable<array_1d<double, 3>>& rOriginVariablesMoments,
    const Variable<array_1d<double, 3>>& rDestinationVariableForces,
    Kratos::Flags MappingOptions)
{
    InitializeOriginForcesAndMoments(
        rOriginVariablesForces,
        rOriginVariablesMoments);

    InitializeInformationBeamsInverse(
        rOriginVariablesForces,
        rOriginVariablesMoments,
        rDestinationVariableForces,
        MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::MapperUniquePointerType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::Clone(
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    Parameters JsonParameters) const
{
    return Kratos::make_unique<BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>>(
        rModelPartOrigin,
        rModelPartDestination,
        JsonParameters);
}

template<class TSparseSpace, class TDenseSpace>
std::string BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::Info() const
{
    return "BeamSplineMapperWithRecoveryOfRotations";
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "BeamSplineMapperWithRecoveryOfRotations";
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::PrintData(std::ostream& rOStream) const
{
    BaseType::PrintData(rOStream);
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::ValidateInput()
{
    Parameters mapper_default_settings(GetMapperDefaultSettings());
    mMapperSettings.ValidateAndAssignDefaults(mapper_default_settings);
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::Initialize()
{
    mBeamChainCache.clear();
    mNodeIdToBeamChainKey.clear();
    CreateLinearSolver();
    InitializeInterfaceCommunicator();
    InitializeInterface();
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::InitializeInterfaceCommunicator()
{
    mpInterfaceCommunicator = Kratos::make_unique<InterfaceCommunicator>(
        mrModelPartOrigin,
        mMapperLocalSystems,
        mMapperSettings["search_settings"]);
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::CreateLinearSolver()
{
    LinearSolverFactory<TSparseSpace, TDenseSpace> solver_factory;

    if (mMapperSettings["linear_solver_settings"].Has("solver_type")) {
        const std::string solver_type =
            mMapperSettings["linear_solver_settings"]["solver_type"].GetString();

        KRATOS_INFO("BeamSplineMapperWithRecoveryOfRotations")
            << "Using specified linear solver: " << solver_type << std::endl;

        mpLinearSolver = solver_factory.Create(mMapperSettings["linear_solver_settings"]);
    } else if (solver_factory.Has("pardiso_lu")) {
        KRATOS_INFO("BeamSplineMapperWithRecoveryOfRotations")
            << "No linear solver specified. Using PardisoLU." << std::endl;

        Parameters default_settings(R"({
            "solver_type" : "pardiso_lu"
        })");
        mpLinearSolver = solver_factory.Create(default_settings);
    } else if (solver_factory.Has("pardiso_ldlt")) {
        KRATOS_INFO("BeamSplineMapperWithRecoveryOfRotations")
            << "No linear solver specified. Using PardisoLDLT." << std::endl;

        Parameters default_settings(R"({
            "solver_type" : "pardiso_ldlt"
        })");
        mpLinearSolver = solver_factory.Create(default_settings);
    } else if (solver_factory.Has("sparse_lu")) {
        KRATOS_INFO("BeamSplineMapperWithRecoveryOfRotations")
            << "No linear solver specified. Using SparseLU." << std::endl;

        Parameters default_settings(R"({
            "solver_type" : "sparse_lu"
        })");
        mpLinearSolver = solver_factory.Create(default_settings);
    } else if (solver_factory.Has("skyline_lu_factorization")) {
        KRATOS_INFO("BeamSplineMapperWithRecoveryOfRotations")
            << "No linear solver specified. Using SkylineLUFactorization." << std::endl;

        Parameters default_settings(R"({
            "solver_type" : "skyline_lu_factorization"
        })");
        mpLinearSolver = solver_factory.Create(default_settings);
    } else {
        KRATOS_WARNING("BeamSplineMapperWithRecoveryOfRotations")
            << "No pivoting direct linear solver is registered. Using AMGCL GMRES fallback." << std::endl;

        Parameters default_settings(R"({
            "solver_type"                  : "amgcl",
            "preconditioner_type"          : "dummy",
            "krylov_type"                  : "gmres",
            "tolerance"                    : 1.0e-12,
            "max_iteration"                : 5000,
            "gmres_krylov_space_dimension" : 500,
            "verbosity"                    : 0,
            "scaling"                      : true,
            "block_size"                   : 1
        })");
        mpLinearSolver = solver_factory.Create(default_settings);
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    mBeamChainCache.clear();
    mNodeIdToBeamChainKey.clear();
    CreateMapperLocalSystems(mrModelPartDestination.GetCommunicator(), mMapperLocalSystems);
    BuildProblem(MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::BuildProblem(Kratos::Flags MappingOptions)
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
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::CreateMapperLocalSystems(
    const Communicator& rModelPartCommunicator,
    std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems)
{
    MapperUtilities::CreateMapperLocalSystemsFromNodes(
        BeamSplineMapperWithRecoveryOfRotationsLocalSystem(nullptr),
        rModelPartCommunicator,
        rLocalSystems);
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::MapDisplacements(
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
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::InitializeInformationBeams(
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

        auto beam_sys = dynamic_cast<BeamSplineMapperWithRecoveryOfRotationsLocalSystem*>(r_local_sys.get());
        KRATOS_ERROR_IF_NOT(beam_sys) << "Expected BeamSplineMapperWithRecoveryOfRotationsLocalSystem!" << std::endl;

        auto projection_data = beam_sys->CalculateProjectionData();

        KRATOS_ERROR_IF_NOT(projection_data.pNode) << "Destination node is a nullptr." << std::endl;

        const auto& r_beam_chain_cache_data = GetOrCreateBeamChainCacheData(projection_data.BeamGeometry);
        const auto source_state_iterator = beam_chain_source_states.find(r_beam_chain_cache_data.Key);
        KRATOS_ERROR_IF(source_state_iterator == beam_chain_source_states.end())
            << "Missing source-state data for beam chain '" << r_beam_chain_cache_data.Key << "'." << std::endl;
        const auto& r_source_state_data = source_state_iterator->second;

        const MatrixType evaluation_matrix = BuildEvaluationMatrix(
            r_beam_chain_cache_data,
            GetReferenceCoordinates(*projection_data.pNode));
        const VectorType global_displacement = EvaluateDisplacement(
            evaluation_matrix,
            r_source_state_data.Coefficients);

        for (IndexType i = 0; i < 3; ++i) {
            projection_data.pNode->FastGetSolutionStepValue(GetComponentVariable(rDestinationVariableDisplacement, i)) =
                global_displacement(i);
        }
    }

}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::InitializeInformationBeamsInverse(
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

        auto beam_sys = dynamic_cast<BeamSplineMapperWithRecoveryOfRotationsLocalSystem*>(r_local_sys.get());
        KRATOS_ERROR_IF_NOT(beam_sys) << "Expected BeamSplineMapperWithRecoveryOfRotationsLocalSystem!" << std::endl;

        auto projection_data = beam_sys->CalculateProjectionData();

        KRATOS_ERROR_IF_NOT(projection_data.pNode) << "Destination node is a nullptr." << std::endl;
        KRATOS_ERROR_IF(projection_data.BeamGeometry.size() != 2)
            << "BeamSplineMapperWithRecoveryOfRotations inverse load transfer expects a 2-node beam segment." << std::endl;
        KRATOS_ERROR_IF(projection_data.LinearShapeValues.size() != 2)
            << "BeamSplineMapperWithRecoveryOfRotations inverse load transfer expects two linear shape function values." << std::endl;

        const auto& r_beam_chain_cache_data = GetOrCreateBeamChainCacheData(projection_data.BeamGeometry);

        array_1d<double, 3> surface_force_global = ZeroVector(3);
        for (IndexType i = 0; i < 3; ++i) {
            surface_force_global[i] =
                projection_data.pNode->FastGetSolutionStepValue(GetComponentVariable(rDestinationVariableForces, i));
        }

        const IndexType number_of_support_nodes = r_beam_chain_cache_data.SupportNodeIds.size();
        const IndexType system_size = r_beam_chain_cache_data.SystemMatrix.size1();
        const MatrixType evaluation_matrix = BuildEvaluationMatrix(
            r_beam_chain_cache_data,
            GetReferenceCoordinates(*projection_data.pNode));

        VectorType adjoint_rhs(system_size, 0.0);
        for (IndexType row = 0; row < 3; ++row) {
            for (IndexType col = 0; col < system_size; ++col) {
                adjoint_rhs(col) += evaluation_matrix(row, col) * surface_force_global[row];
            }
        }

        VectorType generalized_loads(system_size, 0.0);
        mpLinearSolver->Solve(
            r_beam_chain_cache_data.TransposedSparseSystemMatrix,
            generalized_loads,
            adjoint_rhs);

        for (IndexType i = 0; i < number_of_support_nodes; ++i) {
            auto& r_node = mrModelPartOrigin.GetNode(r_beam_chain_cache_data.SupportNodeIds[i]);
            for (IndexType component_index = 0; component_index < 3; ++component_index) {
                const auto& r_origin_force_variable =
                    GetComponentVariable(rOriginVariablesForces, component_index);
                const auto& r_origin_moment_variable =
                    GetComponentVariable(rOriginVariablesMoments, component_index);

                r_node.FastGetSolutionStepValue(r_origin_force_variable) +=
                    factor * generalized_loads(3 * i + component_index);
                r_node.FastGetSolutionStepValue(r_origin_moment_variable) +=
                    factor * 2.0 * generalized_loads(3 * number_of_support_nodes + 3 * i + component_index);
            }
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
double BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::KernelValue(
    const array_1d<double, 3>& rDifference) const
{
    const double rho = norm_2(rDifference);
    if (mKernelType == "cubic") {
        return rho * rho * rho;
    }
    if (mKernelType == "inverse_multiquadric") {
        const double q = mKernelRadius * rho;
        return 1.0 / std::sqrt(1.0 + q * q);
    }
    if (mKernelType == "multiquadric") {
        const double q = rho / mKernelRadius;
        return std::sqrt(1.0 + q * q);
    }
    if (mKernelType == "gaussian") {
        const double q = rho / mKernelRadius;
        return std::exp(-0.5 * q * q);
    }
    if (mKernelType == "thin_plate_spline") {
        if (rho < 1.0e-12) {
            return 0.0;
        }
        const double rho_squared = rho * rho;
        return rho_squared * std::log(rho_squared);
    }
    if (mKernelType == "wendland_c2") {
        const double q = rho / mKernelRadius;
        if (q >= 1.0) {
            return 0.0;
        }
        return std::pow(1.0 - q, 4) * (4.0 * q + 1.0);
    }

    KRATOS_ERROR << "Unsupported kernel_type '" << mKernelType << "'." << std::endl;
}

template<class TSparseSpace, class TDenseSpace>
array_1d<double, 3> BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::KernelGradient(
    const array_1d<double, 3>& rDifference) const
{
    array_1d<double, 3> gradient = ZeroVector(3);
    const double rho = norm_2(rDifference);
    if (rho < std::numeric_limits<double>::epsilon()) {
        return gradient;
    }

    double radial_derivative = 0.0;
    if (mKernelType == "cubic") {
        radial_derivative = 3.0 * rho * rho;
    } else if (mKernelType == "inverse_multiquadric") {
        const double h_squared = mKernelRadius * mKernelRadius;
        const double a = 1.0 + h_squared * rho * rho;
        radial_derivative = -h_squared * rho * std::pow(a, -1.5);
    } else if (mKernelType == "multiquadric") {
        const double h_squared = mKernelRadius * mKernelRadius;
        const double a = 1.0 + rho * rho / h_squared;
        radial_derivative = rho / (h_squared * std::sqrt(a));
    } else if (mKernelType == "gaussian") {
        const double h_squared = mKernelRadius * mKernelRadius;
        radial_derivative = -rho / h_squared * KernelValue(rDifference);
    } else if (mKernelType == "thin_plate_spline") {
        if (rho < 1.0e-12) {
            return gradient;
        }
        radial_derivative = 2.0 * rho * (std::log(rho * rho) + 1.0);
    } else if (mKernelType == "wendland_c2") {
        const double q = rho / mKernelRadius;
        if (q >= 1.0) {
            return gradient;
        }
        radial_derivative = -20.0 * q * std::pow(1.0 - q, 3) / mKernelRadius;
    } else {
        KRATOS_ERROR << "Unsupported kernel_type '" << mKernelType << "'." << std::endl;
    }

    const double factor = radial_derivative / rho;
    for (IndexType i = 0; i < 3; ++i) {
        gradient[i] = factor * rDifference[i];
    }
    return gradient;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::MatrixType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::KernelHessian(
    const array_1d<double, 3>& rDifference) const
{
    const double rho = norm_2(rDifference);
    MatrixType hessian(3, 3, 0.0);
    if (rho < std::numeric_limits<double>::epsilon()) {
        return hessian;
    }

    double radial_derivative = 0.0;
    double radial_second_derivative = 0.0;
    if (mKernelType == "cubic") {
        radial_derivative = 3.0 * rho * rho;
        radial_second_derivative = 6.0 * rho;
    } else if (mKernelType == "inverse_multiquadric") {
        const double h_squared = mKernelRadius * mKernelRadius;
        const double h_fourth = h_squared * h_squared;
        const double a = 1.0 + h_squared * rho * rho;
        radial_derivative = -h_squared * rho * std::pow(a, -1.5);
        radial_second_derivative =
            -h_squared * std::pow(a, -1.5) +
            3.0 * h_fourth * rho * rho * std::pow(a, -2.5);
    } else if (mKernelType == "multiquadric") {
        const double h_squared = mKernelRadius * mKernelRadius;
        const double a = 1.0 + rho * rho / h_squared;
        radial_derivative = rho / (h_squared * std::sqrt(a));
        radial_second_derivative = 1.0 / (h_squared * std::pow(a, 1.5));
    } else if (mKernelType == "gaussian") {
        const double h_squared = mKernelRadius * mKernelRadius;
        const double h_fourth = h_squared * h_squared;
        const double phi = KernelValue(rDifference);
        radial_derivative = -rho / h_squared * phi;
        radial_second_derivative = (rho * rho / h_fourth - 1.0 / h_squared) * phi;
    } else if (mKernelType == "thin_plate_spline") {
        if (rho < 1.0e-12) {
            return hessian;
        }
        radial_derivative = 2.0 * rho * (std::log(rho * rho) + 1.0);
        radial_second_derivative = 2.0 * std::log(rho * rho) + 6.0;
    } else if (mKernelType == "wendland_c2") {
        const double q = rho / mKernelRadius;
        if (q >= 1.0) {
            return hessian;
        }
        radial_derivative = -20.0 * q * std::pow(1.0 - q, 3) / mKernelRadius;
        radial_second_derivative =
            20.0 * std::pow(1.0 - q, 2) * (4.0 * q - 1.0) /
            (mKernelRadius * mKernelRadius);
    } else {
        KRATOS_ERROR << "Unsupported kernel_type '" << mKernelType << "'." << std::endl;
    }

    const double first_over_r = radial_derivative / rho;
    const double tensor_factor =
        radial_second_derivative / (rho * rho) -
        radial_derivative / (rho * rho * rho);
    for (IndexType i = 0; i < 3; ++i) {
        for (IndexType j = 0; j < 3; ++j) {
            hessian(i, j) = tensor_factor * rDifference[i] * rDifference[j];
            if (i == j) {
                hessian(i, j) += first_over_r;
            }
        }
    }
    return hessian;
}

template<class TSparseSpace, class TDenseSpace>
double BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::ApplyFunctionalsToKernel(
    const FunctionalType RowFunctional,
    const FunctionalType ColumnFunctional,
    const array_1d<double, 3>& rRowCoordinates,
    const array_1d<double, 3>& rColumnCoordinates) const
{
    struct ScalarDifferentialOperator
    {
        IndexType ComponentIndex;
        int DerivativeDirection;
        double Coefficient;
    };

    const auto get_operators = [](const FunctionalType Functional) {
        switch (Functional) {
            case FunctionalType::ValueX:
                return std::array<ScalarDifferentialOperator, 2>{{{0, -1, 1.0}, {0, -1, 0.0}}};
            case FunctionalType::ValueY:
                return std::array<ScalarDifferentialOperator, 2>{{{1, -1, 1.0}, {0, -1, 0.0}}};
            case FunctionalType::ValueZ:
                return std::array<ScalarDifferentialOperator, 2>{{{2, -1, 1.0}, {0, -1, 0.0}}};
            case FunctionalType::CurlX:
                return std::array<ScalarDifferentialOperator, 2>{{{1, 2, -1.0}, {2, 1, 1.0}}};
            case FunctionalType::CurlY:
                return std::array<ScalarDifferentialOperator, 2>{{{0, 2, 1.0}, {2, 0, -1.0}}};
            case FunctionalType::CurlZ:
                return std::array<ScalarDifferentialOperator, 2>{{{0, 1, -1.0}, {1, 0, 1.0}}};
        }

        KRATOS_ERROR << "Unsupported functional type." << std::endl;
    };

    array_1d<double, 3> difference = rRowCoordinates - rColumnCoordinates;
    const double phi = KernelValue(difference);
    const array_1d<double, 3> gradient = KernelGradient(difference);
    const MatrixType hessian = KernelHessian(difference);
    const auto row_operators = get_operators(RowFunctional);
    const auto column_operators = get_operators(ColumnFunctional);

    double value = 0.0;
    for (const auto& r_row_operator : row_operators) {
        if (std::abs(r_row_operator.Coefficient) == 0.0) {
            continue;
        }
        for (const auto& r_column_operator : column_operators) {
            if (std::abs(r_column_operator.Coefficient) == 0.0 ||
                r_row_operator.ComponentIndex != r_column_operator.ComponentIndex) {
                continue;
            }

            double derivative_value = 0.0;
            if (r_row_operator.DerivativeDirection < 0 && r_column_operator.DerivativeDirection < 0) {
                derivative_value = phi;
            } else if (r_row_operator.DerivativeDirection >= 0 && r_column_operator.DerivativeDirection < 0) {
                derivative_value = gradient[r_row_operator.DerivativeDirection];
            } else if (r_row_operator.DerivativeDirection < 0 && r_column_operator.DerivativeDirection >= 0) {
                derivative_value = -gradient[r_column_operator.DerivativeDirection];
            } else {
                derivative_value =
                    -hessian(r_row_operator.DerivativeDirection, r_column_operator.DerivativeDirection);
            }

            value += r_row_operator.Coefficient * r_column_operator.Coefficient * derivative_value;
        }
    }

    return value;
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::BuildAijBlock(
    const array_1d<double, 3>& rRowCoordinates,
    const array_1d<double, 3>& rColumnCoordinates,
    MatrixType& rBlock) const
{
    rBlock.resize(6, 6, false);
    const std::array<FunctionalType, 6> functionals{{
        FunctionalType::ValueX,
        FunctionalType::ValueY,
        FunctionalType::ValueZ,
        FunctionalType::CurlX,
        FunctionalType::CurlY,
        FunctionalType::CurlZ}};

    for (IndexType i = 0; i < 6; ++i) {
        for (IndexType j = 0; j < 6; ++j) {
            rBlock(i, j) = ApplyFunctionalsToKernel(
                functionals[i],
                functionals[j],
                rRowCoordinates,
                rColumnCoordinates);
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::BuildPolynomialMatrix(
    const RecoveryOfRotationsCacheData& rCacheData,
    MatrixType& rPolynomialMatrix) const
{
    const auto& r_support_reference_coordinates = rCacheData.SupportReferenceCoordinates;
    const IndexType number_of_support_nodes = r_support_reference_coordinates.size();
    const IndexType polynomial_size = GetVectorPolynomialSize(rCacheData.PolynomialBasis);
    rPolynomialMatrix.resize(6 * number_of_support_nodes, polynomial_size, false);
    rPolynomialMatrix.clear();

    for (IndexType i = 0; i < number_of_support_nodes; ++i) {
        const auto& r_x = r_support_reference_coordinates[i];

        if (UsesLineAdaptedPolynomialBasis(rCacheData.PolynomialBasis)) {
            std::vector<array_1d<double, 3>> values;
            std::vector<array_1d<double, 3>> curls;
            EvaluateLineAdaptedPolynomialModes(rCacheData, r_x, values, curls);

            for (IndexType mode_index = 0; mode_index < polynomial_size; ++mode_index) {
                for (IndexType component_index = 0; component_index < 3; ++component_index) {
                    rPolynomialMatrix(3 * i + component_index, mode_index) = values[mode_index][component_index];
                    rPolynomialMatrix(3 * number_of_support_nodes + 3 * i + component_index, mode_index) =
                        curls[mode_index][component_index];
                }
            }
        } else {
            const std::array<double, Full3DPolynomialBasisSize> p{{1.0, r_x[0], r_x[1], r_x[2]}};
            const std::array<array_1d<double, 3>, Full3DPolynomialBasisSize> grad{{
                ZeroVector(3),
                array_1d<double, 3>{1.0, 0.0, 0.0},
                array_1d<double, 3>{0.0, 1.0, 0.0},
                array_1d<double, 3>{0.0, 0.0, 1.0}}};

            for (IndexType basis_index = 0; basis_index < Full3DPolynomialBasisSize; ++basis_index) {
                for (IndexType component_index = 0; component_index < 3; ++component_index) {
                    rPolynomialMatrix(3 * i + component_index, 3 * basis_index + component_index) =
                        p[basis_index];
                }

                const auto& r_grad = grad[basis_index];
                rPolynomialMatrix(3 * number_of_support_nodes + 3 * i, 3 * basis_index + 1) = -r_grad[2];
                rPolynomialMatrix(3 * number_of_support_nodes + 3 * i, 3 * basis_index + 2) =  r_grad[1];
                rPolynomialMatrix(3 * number_of_support_nodes + 3 * i + 1, 3 * basis_index) =  r_grad[2];
                rPolynomialMatrix(3 * number_of_support_nodes + 3 * i + 1, 3 * basis_index + 2) = -r_grad[0];
                rPolynomialMatrix(3 * number_of_support_nodes + 3 * i + 2, 3 * basis_index) = -r_grad[1];
                rPolynomialMatrix(3 * number_of_support_nodes + 3 * i + 2, 3 * basis_index + 1) =  r_grad[0];
            }
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::EvaluateLineAdaptedPolynomialModes(
    const RecoveryOfRotationsCacheData& rCacheData,
    const array_1d<double, 3>& rCoordinates,
    std::vector<array_1d<double, 3>>& rValues,
    std::vector<array_1d<double, 3>>& rCurls) const
{
    const IndexType number_of_modes = GetVectorPolynomialSize(rCacheData.PolynomialBasis);
    rValues.assign(number_of_modes, ZeroVector(3));
    rCurls.assign(number_of_modes, ZeroVector(3));

    const auto& r_reference = rCacheData.PolynomialReferencePoint;
    const auto& r_tangent = rCacheData.PolynomialTangent;
    const array_1d<double, 3> relative_coordinates = rCoordinates - r_reference;
    const double s =
        relative_coordinates[0] * r_tangent[0] +
        relative_coordinates[1] * r_tangent[1] +
        relative_coordinates[2] * r_tangent[2];

    rValues[0][0] = 1.0;
    rValues[1][1] = 1.0;
    rValues[2][2] = 1.0;

    rValues[3][1] = -relative_coordinates[2];
    rValues[3][2] =  relative_coordinates[1];
    rCurls[3][0] = 2.0;

    rValues[4][0] =  relative_coordinates[2];
    rValues[4][2] = -relative_coordinates[0];
    rCurls[4][1] = 2.0;

    rValues[5][0] = -relative_coordinates[1];
    rValues[5][1] =  relative_coordinates[0];
    rCurls[5][2] = 2.0;

    rValues[6][0] = s;
    rCurls[6][1] =  r_tangent[2];
    rCurls[6][2] = -r_tangent[1];

    rValues[7][1] = s;
    rCurls[7][0] = -r_tangent[2];
    rCurls[7][2] =  r_tangent[0];

    rValues[8][2] = s;
    rCurls[8][0] =  r_tangent[1];
    rCurls[8][1] = -r_tangent[0];

    if (rCacheData.PolynomialBasis == "enriched_line_adapted" ||
        rCacheData.PolynomialBasis == "high_order_line_adapted") {
        const std::array<array_1d<double, 3>, 3> axes{{
            array_1d<double, 3>{1.0, 0.0, 0.0},
            array_1d<double, 3>{0.0, 1.0, 0.0},
            array_1d<double, 3>{0.0, 0.0, 1.0}}};

        for (IndexType axis_index = 0; axis_index < 3; ++axis_index) {
            const array_1d<double, 3> rigid_rotation =
                CrossProduct(axes[axis_index], relative_coordinates);
            const IndexType mode_index = 9 + axis_index;

            rValues[mode_index] = s * rigid_rotation;
            rCurls[mode_index] = CrossProduct(r_tangent, rigid_rotation);
            rCurls[mode_index] += 2.0 * s * axes[axis_index];
        }

        if (rCacheData.PolynomialBasis == "high_order_line_adapted") {
            IndexType mode_index = EnrichedLineAdaptedVectorPolynomialSize;

            // Higher-order line modes improve reproduction of smoothly varying
            // beam rotations and centerline displacements from collinear support.
            for (IndexType degree = 2; degree <= 3; ++degree) {
                const double s_power = std::pow(s, static_cast<int>(degree));
                const double curl_factor =
                    static_cast<double>(degree) * std::pow(s, static_cast<int>(degree - 1));
                for (const auto& r_axis : axes) {
                    rValues[mode_index] = s_power * r_axis;
                    rCurls[mode_index] = curl_factor * CrossProduct(r_tangent, r_axis);
                    ++mode_index;
                }
            }

            for (IndexType degree = 2; degree <= 5; ++degree) {
                const double s_power = std::pow(s, static_cast<int>(degree));
                const double derivative_factor =
                    static_cast<double>(degree) * std::pow(s, static_cast<int>(degree - 1));
                for (const auto& r_axis : axes) {
                    const array_1d<double, 3> rigid_rotation =
                        CrossProduct(r_axis, relative_coordinates);

                    rValues[mode_index] = s_power * rigid_rotation;
                    rCurls[mode_index] = derivative_factor * CrossProduct(r_tangent, rigid_rotation);
                    rCurls[mode_index] += 2.0 * s_power * r_axis;
                    ++mode_index;
                }
            }
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
std::string BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::GetEffectivePolynomialBasis(
    const IndexType NumberOfSupportNodes) const
{
    if (mPolynomialBasis == "high_order_line_adapted") {
        if (NumberOfSupportNodes < 4) {
            return "line_adapted";
        }
        if (NumberOfSupportNodes < 5) {
            return "enriched_line_adapted";
        }
    }

    if (mPolynomialBasis == "enriched_line_adapted" && NumberOfSupportNodes < 3) {
        return "line_adapted";
    }

    return mPolynomialBasis;
}

template<class TSparseSpace, class TDenseSpace>
IndexType BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::GetVectorPolynomialSize(
    const std::string& rPolynomialBasis) const
{
    if (rPolynomialBasis == "line_adapted") {
        return LineAdaptedVectorPolynomialSize;
    }
    if (rPolynomialBasis == "enriched_line_adapted") {
        return EnrichedLineAdaptedVectorPolynomialSize;
    }
    if (rPolynomialBasis == "high_order_line_adapted") {
        return HighOrderLineAdaptedVectorPolynomialSize;
    }
    return Full3DVectorPolynomialSize;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::BuildRightHandSide(
    const RecoveryOfRotationsCacheData& rCacheData,
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rOriginVariablesRotations) const
{
    const IndexType number_of_support_nodes = rCacheData.SupportNodeIds.size();
    VectorType right_hand_side(
        6 * number_of_support_nodes + GetVectorPolynomialSize(rCacheData.PolynomialBasis),
        0.0);

    for (IndexType i = 0; i < number_of_support_nodes; ++i) {
        const auto& r_node = mrModelPartOrigin.GetNode(rCacheData.SupportNodeIds[i]);
        for (IndexType component_index = 0; component_index < 3; ++component_index) {
            right_hand_side(3 * i + component_index) =
                r_node.FastGetSolutionStepValue(GetComponentVariable(rOriginVariablesDisplacements, component_index));
            right_hand_side(3 * number_of_support_nodes + 3 * i + component_index) =
                2.0 * r_node.FastGetSolutionStepValue(GetComponentVariable(rOriginVariablesRotations, component_index));
        }
    }

    return right_hand_side;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::MatrixType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::BuildEvaluationMatrix(
    const RecoveryOfRotationsCacheData& rCacheData,
    const array_1d<double, 3>& rEvaluationCoordinates) const
{
    const IndexType number_of_support_nodes = rCacheData.SupportNodeIds.size();
    const IndexType polynomial_size = GetVectorPolynomialSize(rCacheData.PolynomialBasis);
    const IndexType system_size = 6 * number_of_support_nodes + polynomial_size;
    MatrixType evaluation_matrix(3, system_size, 0.0);
    const std::array<FunctionalType, 3> row_functionals{{
        FunctionalType::ValueX,
        FunctionalType::ValueY,
        FunctionalType::ValueZ}};
    const std::array<FunctionalType, 6> column_functionals{{
        FunctionalType::ValueX,
        FunctionalType::ValueY,
        FunctionalType::ValueZ,
        FunctionalType::CurlX,
        FunctionalType::CurlY,
        FunctionalType::CurlZ}};

    const auto coefficient_index = [number_of_support_nodes](
        const IndexType NodeIndex,
        const IndexType FunctionalIndex) {
        if (FunctionalIndex < 3) {
            return 3 * NodeIndex + FunctionalIndex;
        }
        return 3 * number_of_support_nodes + 3 * NodeIndex + FunctionalIndex - 3;
    };

    for (IndexType i = 0; i < 3; ++i) {
        for (IndexType j = 0; j < number_of_support_nodes; ++j) {
            for (IndexType functional_index = 0; functional_index < 6; ++functional_index) {
                evaluation_matrix(i, coefficient_index(j, functional_index)) = ApplyFunctionalsToKernel(
                    row_functionals[i],
                    column_functionals[functional_index],
                    rEvaluationCoordinates,
                    rCacheData.SupportReferenceCoordinates[j]);
            }
        }

        if (UsesLineAdaptedPolynomialBasis(rCacheData.PolynomialBasis)) {
            std::vector<array_1d<double, 3>> values;
            std::vector<array_1d<double, 3>> curls;
            EvaluateLineAdaptedPolynomialModes(rCacheData, rEvaluationCoordinates, values, curls);
            for (IndexType mode_index = 0; mode_index < polynomial_size; ++mode_index) {
                evaluation_matrix(i, 6 * number_of_support_nodes + mode_index) =
                    values[mode_index][i];
            }
        } else {
            const std::array<double, Full3DPolynomialBasisSize> p{{
                1.0,
                rEvaluationCoordinates[0],
                rEvaluationCoordinates[1],
                rEvaluationCoordinates[2]}};
            for (IndexType basis_index = 0; basis_index < Full3DPolynomialBasisSize; ++basis_index) {
                evaluation_matrix(i, 6 * number_of_support_nodes + 3 * basis_index + i) =
                    p[basis_index];
            }
        }
    }

    return evaluation_matrix;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::EvaluateDisplacement(
    const MatrixType& rEvaluationMatrix,
    const VectorType& rCoefficients) const
{
    KRATOS_ERROR_IF(rEvaluationMatrix.size1() != 3)
        << "Expected a 3-row evaluation matrix." << std::endl;
    KRATOS_ERROR_IF(rEvaluationMatrix.size2() != rCoefficients.size())
        << "Evaluation matrix column count (" << rEvaluationMatrix.size2()
        << ") does not match coefficient size (" << rCoefficients.size() << ")." << std::endl;

    VectorType displacement(3, 0.0);
    TDenseSpace::Mult(rEvaluationMatrix, rCoefficients, displacement);
    return displacement;
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::PrepareBeamChainCacheData()
{
    for (auto& r_local_sys : mMapperLocalSystems) {
        if (!r_local_sys->HasInterfaceInfo()) {
            continue;
        }

        auto beam_sys = dynamic_cast<BeamSplineMapperWithRecoveryOfRotationsLocalSystem*>(r_local_sys.get());
        KRATOS_ERROR_IF_NOT(beam_sys) << "Expected BeamSplineMapperWithRecoveryOfRotationsLocalSystem!" << std::endl;

        GetOrCreateBeamChainCacheData(beam_sys->CalculateProjectionData().BeamGeometry);
    }
}

template<class TSparseSpace, class TDenseSpace>
std::unordered_map<std::string, typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::RecoveryOfRotationsSourceStateData>
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::BuildAllBeamChainSourceStates(
    const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
    const Variable<array_1d<double, 3>>& rOriginVariablesRotations) const
{
    std::unordered_map<std::string, RecoveryOfRotationsSourceStateData> beam_chain_source_states;
    beam_chain_source_states.reserve(mBeamChainCache.size());

    for (const auto& r_key_cache_pair : mBeamChainCache) {
        const auto& r_beam_chain_cache_data = r_key_cache_pair.second;

        RecoveryOfRotationsSourceStateData source_state_data;
        const VectorType right_hand_side = BuildRightHandSide(
            r_beam_chain_cache_data,
            rOriginVariablesDisplacements,
            rOriginVariablesRotations);
        source_state_data.Coefficients = SolveSplineCoefficients(
            r_beam_chain_cache_data.SparseSystemMatrix,
            right_hand_side);

        beam_chain_source_states.emplace(r_beam_chain_cache_data.Key, std::move(source_state_data));
    }

    return beam_chain_source_states;
}

template<class TSparseSpace, class TDenseSpace>
const typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::RecoveryOfRotationsCacheData&
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::GetOrCreateBeamChainCacheData(
    const GeometryType& rBeamGeometry)
{
    const IndexType first_node_id = rBeamGeometry[0].Id();
    const auto key_iterator = mNodeIdToBeamChainKey.find(first_node_id);
    if (key_iterator != mNodeIdToBeamChainKey.end()) {
        return mBeamChainCache.at(key_iterator->second);
    }

    RecoveryOfRotationsCacheData beam_chain_cache_data;
    ComputeBeamChainSupport(
        rBeamGeometry,
        beam_chain_cache_data.SupportNodeIds,
        beam_chain_cache_data.SupportNodeIdToLocalIndex);
    ComputeSupportReferenceData(
        beam_chain_cache_data.SupportNodeIds,
        beam_chain_cache_data.SupportReferenceCoordinates);
    beam_chain_cache_data.PolynomialBasis = GetEffectivePolynomialBasis(
        beam_chain_cache_data.SupportNodeIds.size());
    if (beam_chain_cache_data.PolynomialBasis != mPolynomialBasis) {
        KRATOS_INFO_IF("BeamSplineMapperWithRecoveryOfRotations", mMapperSettings["echo_level"].GetInt() > 0)
            << "Using '" << beam_chain_cache_data.PolynomialBasis
            << "' polynomial basis for a beam chain with "
            << beam_chain_cache_data.SupportNodeIds.size()
            << " support nodes instead of requested '" << mPolynomialBasis
            << "'." << std::endl;
    }
    ComputeLineAdaptedPolynomialData(beam_chain_cache_data);
    beam_chain_cache_data.Key = CreateBeamChainKey(beam_chain_cache_data.SupportNodeIds);
    beam_chain_cache_data.SystemMatrix = BuildSplineSystemMatrix(beam_chain_cache_data);
    beam_chain_cache_data.SparseSystemMatrix =
        BuildSparseMatrix(beam_chain_cache_data.SystemMatrix);
    beam_chain_cache_data.TransposedSparseSystemMatrix =
        BuildTransposedSparseMatrix(beam_chain_cache_data.SystemMatrix);

    for (const IndexType node_id : beam_chain_cache_data.SupportNodeIds) {
        mNodeIdToBeamChainKey[node_id] = beam_chain_cache_data.Key;
    }

    return mBeamChainCache.emplace(
        beam_chain_cache_data.Key,
        std::move(beam_chain_cache_data)).first->second;
}

template<class TSparseSpace, class TDenseSpace>
std::string BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::CreateBeamChainKey(
    const std::vector<IndexType>& rSupportNodeIds) const
{
    std::ostringstream key_stream;
    for (const IndexType node_id : rSupportNodeIds) {
        key_stream << node_id << '-';
    }

    return key_stream.str();
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::MatrixType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::BuildSplineSystemMatrix(
    const RecoveryOfRotationsCacheData& rCacheData) const
{
    const auto& r_support_reference_coordinates = rCacheData.SupportReferenceCoordinates;
    const IndexType number_of_source_nodes = r_support_reference_coordinates.size();
    KRATOS_ERROR_IF(number_of_source_nodes < 2)
        << "BeamSplineMapperWithRecoveryOfRotations requires at least two source coordinates to build a system." << std::endl;

    const IndexType interpolation_size = 6 * number_of_source_nodes;
    const IndexType polynomial_size = GetVectorPolynomialSize(rCacheData.PolynomialBasis);
    const IndexType system_size = interpolation_size + polynomial_size;
    MatrixType spline_system(system_size, system_size, 0.0);

    MatrixType block(6, 6);
    const auto interpolation_index = [number_of_source_nodes](
        const IndexType NodeIndex,
        const IndexType FunctionalIndex) {
        if (FunctionalIndex < 3) {
            return 3 * NodeIndex + FunctionalIndex;
        }
        return 3 * number_of_source_nodes + 3 * NodeIndex + FunctionalIndex - 3;
    };

    for (IndexType i = 0; i < number_of_source_nodes; ++i) {
        for (IndexType j = 0; j < number_of_source_nodes; ++j) {
            BuildAijBlock(
                r_support_reference_coordinates[i],
                r_support_reference_coordinates[j],
                block);
            for (IndexType row = 0; row < 6; ++row) {
                for (IndexType col = 0; col < 6; ++col) {
                    spline_system(
                        interpolation_index(i, row),
                        interpolation_index(j, col)) = block(row, col);
                }
            }
        }
    }

    MatrixType polynomial_matrix;
    BuildPolynomialMatrix(rCacheData, polynomial_matrix);
    for (IndexType i = 0; i < interpolation_size; ++i) {
        for (IndexType j = 0; j < polynomial_size; ++j) {
            spline_system(i, interpolation_size + j) = polynomial_matrix(i, j);
            spline_system(interpolation_size + j, i) = polynomial_matrix(i, j);
        }
    }

    const IndexType regularized_size =
        UsesLineAdaptedPolynomialBasis(rCacheData.PolynomialBasis) ? interpolation_size : system_size;
    for (IndexType i = 0; i < regularized_size; ++i) {
        spline_system(i, i) += mRegularization;
    }

    return spline_system;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::SparseMatrixType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::BuildSparseMatrix(
    const MatrixType& rDenseMatrix) const
{
    const IndexType number_of_rows = rDenseMatrix.size1();
    const IndexType number_of_columns = rDenseMatrix.size2();
    SparseMatrixType sparse_matrix(number_of_rows, number_of_columns);

    for (IndexType i = 0; i < number_of_rows; ++i) {
        for (IndexType j = 0; j < number_of_columns; ++j) {
            const double value = rDenseMatrix(i, j);
            if (std::abs(value) > 0.0) {
                sparse_matrix(i, j) = value;
            }
        }
    }

    return sparse_matrix;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::SparseMatrixType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::BuildTransposedSparseMatrix(
    const MatrixType& rDenseMatrix) const
{
    const IndexType number_of_rows = rDenseMatrix.size1();
    const IndexType number_of_columns = rDenseMatrix.size2();
    SparseMatrixType transposed_sparse_matrix(number_of_columns, number_of_rows);

    for (IndexType i = 0; i < number_of_rows; ++i) {
        for (IndexType j = 0; j < number_of_columns; ++j) {
            const double value = rDenseMatrix(i, j);
            if (std::abs(value) > 0.0) {
                transposed_sparse_matrix(j, i) = value;
            }
        }
    }

    return transposed_sparse_matrix;
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::ComputeSupportReferenceData(
    const std::vector<IndexType>& rSupportNodeIds,
    std::vector<array_1d<double, 3>>& rSupportReferenceCoordinates) const
{
    const IndexType number_of_source_nodes = rSupportNodeIds.size();
    KRATOS_ERROR_IF(number_of_source_nodes < 2)
        << "BeamSplineMapperWithRecoveryOfRotations requires at least two support nodes." << std::endl;

    rSupportReferenceCoordinates.assign(number_of_source_nodes, ZeroVector(3));

    for (IndexType i = 0; i < number_of_source_nodes; ++i) {
        const auto& r_node = mrModelPartOrigin.GetNode(rSupportNodeIds[i]);
        rSupportReferenceCoordinates[i] = GetReferenceCoordinates(r_node);
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::ComputeLineAdaptedPolynomialData(
    RecoveryOfRotationsCacheData& rCacheData) const
{
    const auto& r_support_reference_coordinates = rCacheData.SupportReferenceCoordinates;
    KRATOS_ERROR_IF(r_support_reference_coordinates.size() < 2)
        << "BeamSplineMapperWithRecoveryOfRotations requires at least two support coordinates." << std::endl;

    rCacheData.PolynomialReferencePoint = r_support_reference_coordinates.front();
    array_1d<double, 3> tangent = r_support_reference_coordinates.back() - r_support_reference_coordinates.front();
    const double tangent_norm = std::sqrt(
        tangent[0] * tangent[0] +
        tangent[1] * tangent[1] +
        tangent[2] * tangent[2]);
    KRATOS_ERROR_IF(tangent_norm <= FrameTolerance)
        << "Cannot build a line-adapted polynomial basis from coincident support endpoints." << std::endl;

    tangent /= tangent_norm;
    rCacheData.PolynomialTangent = tangent;
}

template<class TSparseSpace, class TDenseSpace>
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::ComputeBeamChainSupport(
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
            << "BeamSplineMapperWithRecoveryOfRotations only supports non-branching beam chains." << std::endl;

        if (local_degree == 1) {
            start_node_id = node_id;
        }
    }

    KRATOS_ERROR_IF(start_node_id == 0)
        << "BeamSplineMapperWithRecoveryOfRotations requires an open beam chain with two end nodes." << std::endl;

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
void BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::InitializeOriginForcesAndMoments(
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
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::SolveSplineCoefficients(
    SparseMatrixType& rSystemMatrix,
    const VectorType& rRightHandSide) const
{
    const IndexType system_size = rSystemMatrix.size1();
    KRATOS_ERROR_IF(system_size == 0)
        << "Cannot solve an empty spline coefficient system." << std::endl;
    KRATOS_ERROR_IF(rSystemMatrix.size2() != system_size)
        << "Spline coefficient system matrix must be square. Got "
        << rSystemMatrix.size1() << " x " << rSystemMatrix.size2() << "." << std::endl;
    KRATOS_ERROR_IF(rRightHandSide.size() != system_size)
        << "Spline coefficient right-hand side size (" << rRightHandSide.size()
        << ") does not match system size (" << system_size << ")." << std::endl;

    VectorType coefficients(system_size, 0.0);

    KRATOS_ERROR_IF_NOT(mpLinearSolver)
        << "BeamSplineMapperWithRecoveryOfRotations linear solver was not created before solving the spline system." << std::endl;

    VectorType right_hand_side = rRightHandSide;
    mpLinearSolver->Solve(
        rSystemMatrix,
        coefficients,
        right_hand_side);

    return coefficients;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::TransformGlobalToLocal(
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
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::TransformVectorToLocal(
    const MatrixType& rRotationMatrixGlobalToLocal,
    const array_1d<double, 3>& rVector) const
{
    VectorType local_vector(3, 0.0);
    TDenseSpace::Mult(rRotationMatrixGlobalToLocal, rVector, local_vector);
    return local_vector;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::TransformVectorToGlobal(
    const MatrixType& rRotationMatrixLocalToGlobal,
    const VectorType& rLocalVector) const
{
    VectorType global_vector(3, 0.0);
    TDenseSpace::Mult(rRotationMatrixLocalToGlobal, rLocalVector, global_vector);
    return global_vector;
}

template<class TSparseSpace, class TDenseSpace>
array_1d<double, 3> BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::GetReferenceCoordinates(
    const Node& rNode) const
{
    array_1d<double, 3> reference_coordinates = ZeroVector(3);
    reference_coordinates[0] = rNode.X0();
    reference_coordinates[1] = rNode.Y0();
    reference_coordinates[2] = rNode.Z0();
    return reference_coordinates;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::MatrixType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::ComputeInitialFrameLocalToGlobal(
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
array_1d<double, 3> BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::RotateVectorAroundAxis(
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
typename BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::MapperInterfaceInfoUniquePointerType
BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::GetMapperInterfaceInfo() const
{
    return Kratos::make_unique<BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo>(mLocalCoordTol);
}

template<class TSparseSpace, class TDenseSpace>
Parameters BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::GetMapperDefaultSettings() const
{
    return Parameters(R"({
        "search_settings"              : {},
        "search_radius"                : -1.0,
        "search_iterations"            : 3,
        "local_coord_tolerance"        : 0.25,
        "kernel_type"                  : "cubic",
        "kernel_radius"                : 1.0,
        "regularization"               : 1.0e-8,
        "polynomial_basis"             : "auto",
        "linear_solver_settings"       : {},
        "echo_level"                   : 0
    })");
}

template<class TSparseSpace, class TDenseSpace>
bool BeamSplineMapperWithRecoveryOfRotations<TSparseSpace, TDenseSpace>::UsesLineAdaptedPolynomialBasis(
    const std::string& rPolynomialBasis) const
{
    return rPolynomialBasis == "line_adapted" ||
           rPolynomialBasis == "enriched_line_adapted" ||
           rPolynomialBasis == "high_order_line_adapted";
}

template class BeamSplineMapperWithRecoveryOfRotations<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>;

}  // namespace Kratos
