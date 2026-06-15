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
#include <iostream>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

// External includes
#include "utilities/math_utils.h"
#include "geometries/line_3d_2.h"

// Project includes
#include "beam_spline_mapper.h"
#include "factories/linear_solver_factory.h"
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

// TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
template<class TVectorType>
void PrintDebugVector(
    const std::string& rName,
    const TVectorType& rVector)
{
    (void)rName;
    (void)rVector;
    // Debug output intentionally disabled. Keep this helper as a no-op so
    // temporary diagnostic call sites do not affect the mapper calculation.
}

// TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
template<class TMatrixType>
void PrintDebugMatrix(
    const std::string& rName,
    const TMatrixType& rMatrix)
{
    (void)rName;
    (void)rMatrix;
    // Debug output intentionally disabled. Keep this helper as a no-op so
    // temporary diagnostic call sites do not affect the mapper calculation.
}
}

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
}

BeamSplineMapperLocalSystem::ProjectionData BeamSplineMapperLocalSystem::CalculateProjectionData()
{
    ProjectionData projection_data;
    projection_data.pNode = mpNode;

    KRATOS_ERROR_IF_NOT(mpNode) << "Destination node is a nullptr." << std::endl;
    KRATOS_ERROR_IF(mInterfaceInfos.empty())
        << "Cannot calculate BeamSpline projection data without interface information." << std::endl;

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
    CreateLinearSolver();
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
void BeamSplineMapper<TSparseSpace, TDenseSpace>::CreateLinearSolver()
{
    LinearSolverFactory<TSparseSpace, TDenseSpace> solver_factory;

    if (mMapperSettings["linear_solver_settings"].Has("solver_type")) {
        const std::string solver_type =
            mMapperSettings["linear_solver_settings"]["solver_type"].GetString();

        KRATOS_INFO("BeamSplineMapper")
            << "Using specified linear solver: " << solver_type << std::endl;

        mpLinearSolver = solver_factory.Create(mMapperSettings["linear_solver_settings"]);
    } else if (solver_factory.Has("pardiso_lu")) {
        KRATOS_INFO("BeamSplineMapper")
            << "No linear solver specified. Using PardisoLU." << std::endl;

        Parameters default_settings(R"({
            "solver_type" : "pardiso_lu"
        })");
        mpLinearSolver = solver_factory.Create(default_settings);
    } else if (solver_factory.Has("pardiso_ldlt")) {
        KRATOS_INFO("BeamSplineMapper")
            << "No linear solver specified. Using PardisoLDLT." << std::endl;

        Parameters default_settings(R"({
            "solver_type" : "pardiso_ldlt"
        })");
        mpLinearSolver = solver_factory.Create(default_settings);
    } else if (solver_factory.Has("sparse_lu")) {
        KRATOS_INFO("BeamSplineMapper")
            << "No linear solver specified. Using SparseLU." << std::endl;

        Parameters default_settings(R"({
            "solver_type" : "sparse_lu"
        })");
        mpLinearSolver = solver_factory.Create(default_settings);
    } else {
        KRATOS_WARNING("BeamSplineMapper")
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
            projection_coordinate);

        const double axial_displacement =
            projection_data.LinearShapeValues(0) * r_source_state_data.LocalDisplacements[0][first_segment_node_index] +
            projection_data.LinearShapeValues(1) * r_source_state_data.LocalDisplacements[0][second_segment_node_index];

        const double torsional_rotation =
            projection_data.LinearShapeValues(0) * r_source_state_data.LocalRotations[0][first_segment_node_index] +
            projection_data.LinearShapeValues(1) * r_source_state_data.LocalRotations[0][second_segment_node_index];
        VectorType evaluation_rotation_vector(3);
        evaluation_rotation_vector(0) = torsional_rotation;
        evaluation_rotation_vector(1) = -inner_prod(evaluation_derivative_row, r_source_state_data.SplineCoefficientsZ);
        evaluation_rotation_vector(2) = inner_prod(evaluation_derivative_row, r_source_state_data.SplineCoefficientsY);
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

        // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
        VectorType initial_destination_displacement(3, 0.0);
        for (IndexType i = 0; i < 3; ++i) {
            initial_destination_displacement(i) =
                projection_data.pNode->FastGetSolutionStepValue(GetComponentVariable(rDestinationVariableDisplacement, i));
        }
        // TEMP DEBUG OUTPUT DISABLED - REMOVE AFTER VERIFICATION
        // std::cout << "[BeamSplineMapper DEBUG] MapDisplacements::destination_node_id = "
        //           << projection_data.pNode->Id() << std::endl;
        // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
        PrintDebugVector("MapDisplacements::initial_destination_displacement_global", initial_destination_displacement);
        // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
        PrintDebugVector("MapDisplacements::final_destination_displacement_local", local_displacement);
        // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
        PrintDebugVector("MapDisplacements::final_destination_displacement_global", global_displacement);

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
            projection_coordinate);

        array_1d<double, 3> surface_force_global = ZeroVector(3);
        for (IndexType i = 0; i < 3; ++i) {
            surface_force_global[i] =
                projection_data.pNode->FastGetSolutionStepValue(GetComponentVariable(rDestinationVariableForces, i));
        }

        const VectorType surface_force_local = TransformVectorToLocal(
            evaluation_rotation_matrix_global_to_local,
            surface_force_global);

        const array_1d<double, 3> projection_point_reference =
            MakeArrayFromVector(projection_data.ProjectionPoint);
        const array_1d<double, 3> destination_point_reference =
            GetReferenceCoordinates(*projection_data.pNode);
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
        TDenseSpace::Mult(
            r_beam_chain_cache_data.InverseTransposedSplineSystemMatrix,
            right_hand_side_y,
            adjoint_coefficients_y);
        TDenseSpace::Mult(
            r_beam_chain_cache_data.InverseTransposedSplineSystemMatrix,
            right_hand_side_z,
            adjoint_coefficients_z);

        for (IndexType i = 0; i < number_of_support_nodes; ++i) {
            local_forces[1][i] += adjoint_coefficients_y(i);
            local_moments[2][i] += adjoint_coefficients_y(number_of_support_nodes + i);
            local_forces[2][i] += adjoint_coefficients_z(i);
            local_moments[1][i] -= adjoint_coefficients_z(number_of_support_nodes + i);
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
    double abs_distance = std::abs(Distance);
    double KernelValue = abs_distance * abs_distance * abs_distance;

    return KernelValue;
}

template<class TSparseSpace, class TDenseSpace>
double BeamSplineMapper<TSparseSpace, TDenseSpace>::EvaluateKernelFirstDerivative(const double Distance) const
{
    double abs_distance = std::abs(Distance);
    double KernelFirstDerivative = 3.0 * Distance * abs_distance;

    return KernelFirstDerivative;
}

template<class TSparseSpace, class TDenseSpace>
double BeamSplineMapper<TSparseSpace, TDenseSpace>::EvaluateKernelSecondDerivative(const double Distance) const
{
    const double abs_distance = std::abs(Distance);
    if (abs_distance < PolynomialTolerance) {
        return 0.0;
    }
    return 6.0 * abs_distance;
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
    KRATOS_ERROR_IF(number_of_source_nodes < 2)
        << "BeamSplineMapper requires at least two source coordinates to build a spline system." << std::endl;

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
    KRATOS_ERROR_IF(number_of_source_nodes < 2)
        << "BeamSplineMapper requires at least two source coordinates to build an evaluation row." << std::endl;

    VectorType evaluation_row(2 * number_of_source_nodes + 2, 0.0);

    for (IndexType j = 0; j < number_of_source_nodes; ++j) {
        const double distance = ProjectionCoordinate - rSourceCoordinates[j];
        evaluation_row(j) = EvaluateKernelValue(distance);
        evaluation_row(number_of_source_nodes + j) = EvaluateKernelFirstDerivative(distance);
    }

    evaluation_row(2 * number_of_source_nodes) = 1.0;
    evaluation_row(2 * number_of_source_nodes + 1) = ProjectionCoordinate;

    /*
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("BuildEvaluationRow::source_coordinates", rSourceCoordinates);
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    std::cout << "[BeamSplineMapper DEBUG] BuildEvaluationRow::projection_coordinate = "
              << ProjectionCoordinate << std::endl;
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("BuildEvaluationRow::C_f_[Mfs_Mprimefs_Pf]", evaluation_row);
    */

    return evaluation_row;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildEvaluationDerivativeRow(
    const std::vector<double>& rSourceCoordinates,
    const double ProjectionCoordinate) const
{
    const IndexType number_of_source_nodes = rSourceCoordinates.size();
    KRATOS_ERROR_IF(number_of_source_nodes < 2)
        << "BeamSplineMapper requires at least two source coordinates to build an evaluation derivative row." << std::endl;

    VectorType evaluation_derivative_row(2 * number_of_source_nodes + 2, 0.0);

    for (IndexType j = 0; j < number_of_source_nodes; ++j) {
        const double distance = ProjectionCoordinate - rSourceCoordinates[j];
        evaluation_derivative_row(j) = EvaluateKernelFirstDerivative(distance);
        evaluation_derivative_row(number_of_source_nodes + j) = EvaluateKernelSecondDerivative(distance);
    }

    evaluation_derivative_row(2 * number_of_source_nodes) = 0.0;
    evaluation_derivative_row(2 * number_of_source_nodes + 1) = 1.0;

    return evaluation_derivative_row;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::VectorType
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildRightHandSide(
    const std::vector<double>& rDisplacements,
    const std::vector<double>& rRotations) const
{
    const IndexType number_of_source_nodes = rDisplacements.size();
    KRATOS_ERROR_IF(number_of_source_nodes < 2)
        << "BeamSplineMapper requires at least two displacement values to build the spline right-hand side." << std::endl;
    KRATOS_ERROR_IF(rRotations.size() != number_of_source_nodes)
        << "The number of spline rotation values (" << rRotations.size()
        << ") does not match the number of displacement values (" << number_of_source_nodes << ")." << std::endl;

    VectorType right_hand_side(2 * number_of_source_nodes + 2, 0.0);

    for (IndexType i = 0; i < number_of_source_nodes; ++i) {
        right_hand_side(i) = rDisplacements[i];
        right_hand_side(number_of_source_nodes + i) = rRotations[i];
    }

    /*
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("BuildRightHandSide::displacements", rDisplacements);
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("BuildRightHandSide::rotations", rRotations);
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("BuildRightHandSide::right_hand_side", right_hand_side);
    */

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

    /*
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("BuildLocalSourceData::local_displacements_x", rLocalDisplacements[0]);
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("BuildLocalSourceData::local_displacements_y", rLocalDisplacements[1]);
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("BuildLocalSourceData::local_displacements_z", rLocalDisplacements[2]);
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("BuildLocalSourceData::local_rotations_x", rLocalRotations[0]);
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("BuildLocalSourceData::local_rotations_y", rLocalRotations[1]);
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("BuildLocalSourceData::local_rotations_z", rLocalRotations[2]);
    */
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
            source_state_data.LocalRotations[2]);

        std::vector<double> negative_local_rotation_y(
            source_state_data.LocalRotations[1].size(),
            0.0);
        for (IndexType i = 0; i < source_state_data.LocalRotations[1].size(); ++i) {
            negative_local_rotation_y[i] = -source_state_data.LocalRotations[1][i];
        }

        const VectorType right_hand_side_z = BuildRightHandSide(
            source_state_data.LocalDisplacements[2],
            negative_local_rotation_y);

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
        beam_chain_cache_data.SupportFramesLocalToGlobal,
        beam_chain_cache_data.SupportReferenceCoordinates);
    beam_chain_cache_data.Key = CreateBeamChainKey(beam_chain_cache_data.SupportNodeIds);
    beam_chain_cache_data.SplineSystemMatrix = BuildSplineSystemMatrix(
        beam_chain_cache_data.LocalXCoordinates);
    beam_chain_cache_data.InverseTransposedSplineSystemMatrix =
        BuildInverseTransposedMatrix(beam_chain_cache_data.SplineSystemMatrix);

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

    /*
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("CreateBeamChainKey::support_node_ids", rSupportNodeIds);
    */

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

        spline_system(i, 2 * number_of_source_nodes) = ps(i, 0);
        spline_system(i, 2 * number_of_source_nodes + 1) = ps(i, 1);
        spline_system(number_of_source_nodes + i, 2 * number_of_source_nodes) =
            ps_first_derivative(i, 0);
        spline_system(number_of_source_nodes + i, 2 * number_of_source_nodes + 1) =
            ps_first_derivative(i, 1);

        spline_system(2 * number_of_source_nodes, i) = ps(i, 0);
        spline_system(2 * number_of_source_nodes + 1, i) = ps(i, 1);
        spline_system(2 * number_of_source_nodes, number_of_source_nodes + i) =
            ps_first_derivative(i, 0);
        spline_system(2 * number_of_source_nodes + 1, number_of_source_nodes + i) =
            ps_first_derivative(i, 1);
    }

    /*
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("BuildSplineSystemMatrix::source_coordinates", rSourceCoordinates);
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugMatrix("BuildSplineSystemMatrix::A", spline_system);
    */

    return spline_system;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::MatrixType
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildInverseTransposedMatrix(
    const MatrixType& rMatrix) const
{
    const IndexType matrix_size = rMatrix.size1();
    KRATOS_ERROR_IF(matrix_size == 0)
        << "Cannot invert an empty matrix." << std::endl;
    KRATOS_ERROR_IF(rMatrix.size2() != matrix_size)
        << "Cannot invert the transpose of a non-square matrix. Got "
        << rMatrix.size1() << " x " << rMatrix.size2() << "." << std::endl;

    MatrixType transposed_matrix(matrix_size, matrix_size, 0.0);
    for (IndexType i = 0; i < matrix_size; ++i) {
        for (IndexType j = 0; j < matrix_size; ++j) {
            transposed_matrix(i, j) = rMatrix(j, i);
        }
    }

    MatrixType inverse_transposed_matrix(matrix_size, matrix_size, 0.0);
    SparseMatrixType sparse_transposed_matrix = BuildSparseMatrix(transposed_matrix);
    VectorType right_hand_side(matrix_size, 0.0);
    VectorType solution(matrix_size, 0.0);
    for (IndexType column_index = 0; column_index < matrix_size; ++column_index) {
        std::fill(right_hand_side.begin(), right_hand_side.end(), 0.0);
        std::fill(solution.begin(), solution.end(), 0.0);
        right_hand_side(column_index) = 1.0;

        mpLinearSolver->Solve(
            sparse_transposed_matrix,
            solution,
            right_hand_side);

        for (IndexType row_index = 0; row_index < matrix_size; ++row_index) {
            inverse_transposed_matrix(row_index, column_index) = solution(row_index);
        }
    }

    return inverse_transposed_matrix;
}

template<class TSparseSpace, class TDenseSpace>
typename BeamSplineMapper<TSparseSpace, TDenseSpace>::SparseMatrixType
BeamSplineMapper<TSparseSpace, TDenseSpace>::BuildSparseMatrix(
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
void BeamSplineMapper<TSparseSpace, TDenseSpace>::ComputeSupportReferenceData(
    const std::vector<IndexType>& rSupportNodeIds,
    std::vector<double>& rLocalXCoordinates,
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

    /*
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("ComputeSupportReferenceData::support_node_ids", rSupportNodeIds);
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("ComputeSupportReferenceData::local_x_coordinates", rLocalXCoordinates);
    for (IndexType i = 0; i < rSupportFramesLocalToGlobal.size(); ++i) {
        // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
        PrintDebugMatrix("ComputeSupportReferenceData::support_frame_local_to_global[" + std::to_string(i) + "]",
            rSupportFramesLocalToGlobal[i]);
    }
    */
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
        << "Spline coefficient system matrix must be square. Got "
        << rSplineSystemMatrix.size1() << " x " << rSplineSystemMatrix.size2() << "." << std::endl;
    KRATOS_ERROR_IF(rRightHandSide.size() != system_size)
        << "Spline coefficient right-hand side size (" << rRightHandSide.size()
        << ") does not match system size (" << system_size << ")." << std::endl;

    VectorType coefficients(system_size, 0.0);

    KRATOS_ERROR_IF_NOT(mpLinearSolver)
        << "BeamSplineMapper linear solver was not created before solving the spline system." << std::endl;

    SparseMatrixType sparse_spline_system = BuildSparseMatrix(rSplineSystemMatrix);
    VectorType right_hand_side = rRightHandSide;
    mpLinearSolver->Solve(
        sparse_spline_system,
        coefficients,
        right_hand_side);

    /*
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugMatrix("SolveSplineCoefficients::A", rSplineSystemMatrix);
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("SolveSplineCoefficients::right_hand_side", rRightHandSide);
    // TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    PrintDebugVector("SolveSplineCoefficients::coefficients", coefficients);
    */

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

    MatrixType rotation_matrix(3, 3);
    MatrixType tmp_rotation(3, 3);
    tmp_rotation = prod(rotation_x, rotation_z);
    rotation_matrix = prod(tmp_rotation, rotation_y);

    // get R * rOffset 
    VectorType rotated_offset(3, 0.0);
    TDenseSpace::Mult(rotation_matrix, rOffsetVectorLocal, rotated_offset);


    VectorType local_displacement(3);
    // Δr = R * rOffset - rOffset
    //       = (R - I) * rOffset
    // u_local = u_centerline + Δr
    //         = u_centerline + (R - I) * rOffset
    local_displacement(0) = centerline_displacement(0) + rotated_offset(0) - rOffsetVectorLocal(0);
    local_displacement(1) = centerline_displacement(1) + rotated_offset(1) - rOffsetVectorLocal(1);
    local_displacement(2) = centerline_displacement(2) + rotated_offset(2) - rOffsetVectorLocal(2);

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
        "linear_solver_settings"       : {},
        "echo_level"                   : 0
    })");
}

template class BeamSplineMapper<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>;

}  // namespace Kratos
