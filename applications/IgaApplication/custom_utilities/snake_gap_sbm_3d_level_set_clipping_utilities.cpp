//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#include "custom_utilities/snake_gap_sbm_3d_level_set_clipping_utilities.h"

// System includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

// Project includes
#include "geometries/coupling_geometry.h"
#include "spatial_containers/geometrical_objects_bins.h"
#include "geometries/quadrature_point_geometry.h"
#include "geometries/geometry_shape_function_container.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_3d_3.h"
#include "integration/integration_point.h"
#include "utilities/quadrature_points_utility.h"

namespace Kratos
{
namespace
{

using NodeType = Node;
using NodePointerType = NodeType::Pointer;
using GeometryType = Geometry<NodeType>;
using GeometryPointerType = GeometryType::Pointer;

struct SkinTriangle
{
    std::array<NodePointerType, 3> Nodes;
    array_1d<double, 3> Normal = ZeroVector(3);
    const Condition* pSourceCondition = nullptr;
};

struct SurrogateFaceOrientation
{
    array_1d<double, 3> Center = ZeroVector(3);
    array_1d<double, 3> Normal = ZeroVector(3);
};

struct ProjectionResult
{
    array_1d<double, 3> Point = ZeroVector(3);
    array_1d<double, 3> Normal = ZeroVector(3);
    NodePointerType pNearestSkinNode;
    double DistanceSquared = std::numeric_limits<double>::max();
};

struct Vertex
{
    NodePointerType pNode;
    double Phi = 0.0;
};

struct EdgeKey
{
    const NodeType* pFirst = nullptr;
    const NodeType* pSecond = nullptr;

    EdgeKey() = default;

    EdgeKey(const NodePointerType& pNode0, const NodePointerType& pNode1)
    {
        pFirst = pNode0.get();
        pSecond = pNode1.get();
        if (std::less<const NodeType*>{}(pSecond, pFirst)) {
            std::swap(pFirst, pSecond);
        }
    }

    bool operator<(const EdgeKey& rOther) const
    {
        return std::tie(pFirst, pSecond) < std::tie(rOther.pFirst, rOther.pSecond);
    }
};

struct FaceKey
{
    std::array<const NodeType*, 3> Nodes = {{nullptr, nullptr, nullptr}};

    FaceKey() = default;

    FaceKey(
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2)
    {
        Nodes = {{pNode0.get(), pNode1.get(), pNode2.get()}};
        std::sort(Nodes.begin(), Nodes.end(), std::less<const NodeType*>{});
    }

    bool operator<(const FaceKey& rOther) const
    {
        return Nodes < rOther.Nodes;
    }
};

struct FaceOccurrence
{
    std::array<NodePointerType, 3> Nodes;
    GeometryPointerType pNeighbourGeometry;
    array_1d<double, 3> OwnerSubcellCenter = ZeroVector(3);
    bool IsCutFacet = false;
};

using CutEdgeRegistry = std::map<EdgeKey, NodePointerType>;
using SubcellFaceRegistry = std::map<FaceKey, std::vector<FaceOccurrence>>;
using GridVertexRegistry = std::map<std::array<int, 3>, NodePointerType>;

constexpr double RelativeTolerance = 1.0e-12;

const array_1d<double, 3>& Coordinates(const NodePointerType& pNode)
{
    return pNode->Coordinates();
}

double SquaredNorm(const array_1d<double, 3>& rVector)
{
    return inner_prod(rVector, rVector);
}

double TetraVolume(
    const NodePointerType& pNode0,
    const NodePointerType& pNode1,
    const NodePointerType& pNode2,
    const NodePointerType& pNode3)
{
    const auto edge01 = Coordinates(pNode1) - Coordinates(pNode0);
    const auto edge02 = Coordinates(pNode2) - Coordinates(pNode0);
    const auto edge03 = Coordinates(pNode3) - Coordinates(pNode0);
    array_1d<double, 3> cross_product = ZeroVector(3);
    MathUtils<double>::CrossProduct(cross_product, edge02, edge03);
    return std::abs(inner_prod(edge01, cross_product)) / 6.0;
}

array_1d<double, 3> ClosestPointOnTriangle(
    const array_1d<double, 3>& rPoint,
    const array_1d<double, 3>& rPoint0,
    const array_1d<double, 3>& rPoint1,
    const array_1d<double, 3>& rPoint2)
{
    const auto edge01 = rPoint1 - rPoint0;
    const auto edge02 = rPoint2 - rPoint0;
    const auto point0 = rPoint - rPoint0;
    const double dot00 = inner_prod(edge01, point0);
    const double dot01 = inner_prod(edge02, point0);
    if (dot00 <= 0.0 && dot01 <= 0.0) {
        return rPoint0;
    }

    const auto point1 = rPoint - rPoint1;
    const double dot10 = inner_prod(edge01, point1);
    const double dot11 = inner_prod(edge02, point1);
    if (dot10 >= 0.0 && dot11 <= dot10) {
        return rPoint1;
    }

    const double edge01_region = dot00 * dot11 - dot10 * dot01;
    if (edge01_region <= 0.0 && dot00 >= 0.0 && dot10 <= 0.0) {
        return rPoint0 + edge01 * (dot00 / (dot00 - dot10));
    }

    const auto point2 = rPoint - rPoint2;
    const double dot20 = inner_prod(edge01, point2);
    const double dot21 = inner_prod(edge02, point2);
    if (dot21 >= 0.0 && dot20 <= dot21) {
        return rPoint2;
    }

    const double edge02_region = dot20 * dot01 - dot00 * dot21;
    if (edge02_region <= 0.0 && dot01 >= 0.0 && dot21 <= 0.0) {
        return rPoint0 + edge02 * (dot01 / (dot01 - dot21));
    }

    const double edge12_region = dot10 * dot21 - dot20 * dot11;
    if (edge12_region <= 0.0 && dot11 - dot10 >= 0.0 && dot20 - dot21 >= 0.0) {
        return rPoint1 + (rPoint2 - rPoint1) *
            ((dot11 - dot10) / ((dot11 - dot10) + (dot20 - dot21)));
    }

    const double denominator = 1.0 / (edge01_region + edge02_region + edge12_region);
    const double barycentric1 = edge02_region * denominator;
    const double barycentric2 = edge01_region * denominator;
    return rPoint0 + edge01 * barycentric1 + edge02 * barycentric2;
}

ProjectionResult ProjectToSkin(
    const array_1d<double, 3>& rPoint,
    const std::vector<SkinTriangle>& rSkinTriangles,
    GeometricalObjectsBins* pSkinBins = nullptr,
    const std::unordered_map<const Condition*, std::size_t>* pTriangleByCondition = nullptr,
    const double InitialSearchRadius = 0.0)
{
    ProjectionResult result;
    const auto update_projection = [&rPoint, &result, &rSkinTriangles](const std::size_t TriangleIndex) {
        const auto& r_triangle = rSkinTriangles[TriangleIndex];
        const auto projection = ClosestPointOnTriangle(
            rPoint,
            Coordinates(r_triangle.Nodes[0]),
            Coordinates(r_triangle.Nodes[1]),
            Coordinates(r_triangle.Nodes[2]));
        const double distance_squared = SquaredNorm(rPoint - projection);
        if (distance_squared >= result.DistanceSquared) {
            return;
        }

        result.Point = projection;
        result.Normal = r_triangle.Normal;
        result.DistanceSquared = distance_squared;
        result.pNearestSkinNode = r_triangle.Nodes[0];
        for (const auto& p_node : r_triangle.Nodes) {
            if (SquaredNorm(rPoint - Coordinates(p_node)) <
                SquaredNorm(rPoint - Coordinates(result.pNearestSkinNode))) {
                result.pNearestSkinNode = p_node;
            }
        }
    };

    // For the normal (triangular) skin conditions used by this workflow the
    // geometrical-object bins returns the exact closest condition, avoiding a
    // full scan of all skin triangles for every generated condition QP.
    if (pSkinBins != nullptr && pTriangleByCondition != nullptr) {
        Point search_point(rPoint[0], rPoint[1], rPoint[2]);
        GeometricalObjectsBins::ResultType search_result;
        double radius = std::max(InitialSearchRadius, 1.0e-12);
        for (int attempt = 0; attempt < 8 && !search_result.GetIsObjectFound(); ++attempt) {
            search_result = pSkinBins->SearchNearestInRadius(search_point, radius);
            radius *= 2.0;
        }
        if (!search_result.GetIsObjectFound()) {
            search_result = pSkinBins->SearchNearest(search_point);
        }
        if (search_result.GetIsObjectFound()) {
            const auto* p_condition = dynamic_cast<const Condition*>(search_result.Get().get());
            const auto it = p_condition != nullptr ? pTriangleByCondition->find(p_condition) :
                pTriangleByCondition->end();
            if (it != pTriangleByCondition->end()) {
                update_projection(it->second);
                return result;
            }
        }
    }
    for (std::size_t i = 0; i < rSkinTriangles.size(); ++i) {
        update_projection(i);
    }
    return result;
}

bool ContainsNode(
    const std::vector<NodePointerType>& rNodes,
    const NodePointerType& pNode)
{
    return std::any_of(rNodes.begin(), rNodes.end(), [&pNode](const auto& pCandidate) {
        return pCandidate.get() == pNode.get();
    });
}

void AddUniqueNode(std::vector<NodePointerType>& rNodes, const NodePointerType& pNode)
{
    if (pNode && !ContainsNode(rNodes, pNode)) {
        rNodes.push_back(pNode);
    }
}

void RotateToCanonicalFirstNode(std::vector<NodePointerType>& rNodes)
{
    const auto it = std::min_element(
        rNodes.begin(), rNodes.end(),
        [](const auto& pNode0, const auto& pNode1) {
            return std::less<const NodeType*>{}(pNode0.get(), pNode1.get());
        });
    std::rotate(rNodes.begin(), it, rNodes.end());
}

GeometryData::IntegrationMethod SelectTetraIntegrationMethod(const std::size_t Order);

std::vector<GeometryPointerType> CreateFacetQuadraturePointGeometries(
    const std::array<NodePointerType, 3>& rNodes,
    const std::size_t IntegrationOrder)
{
    const auto edge01 = Coordinates(rNodes[1]) - Coordinates(rNodes[0]);
    const auto edge02 = Coordinates(rNodes[2]) - Coordinates(rNodes[0]);
    array_1d<double, 3> normal = ZeroVector(3);
    MathUtils<double>::CrossProduct(normal, edge01, edge02);
    const double area = 0.5 * norm_2(normal);
    KRATOS_ERROR_IF(area <= std::numeric_limits<double>::epsilon())
        << "Cannot create a quadrature point from a degenerate cut facet.\n";

    PointerVector<NodeType> points;
    points.push_back(rNodes[0]);
    points.push_back(rNodes[1]);
    points.push_back(rNodes[2]);

    Triangle3D3<NodeType> triangle(rNodes[0], rNodes[1], rNodes[2]);
    const auto& integration_points = triangle.IntegrationPoints(
        SelectTetraIntegrationMethod(IntegrationOrder));

    std::vector<GeometryPointerType> result;
    result.reserve(integration_points.size());
    for (const auto& r_integration_point : integration_points) {
        array_1d<double, 3> local_coordinates = ZeroVector(3);
        local_coordinates[0] = r_integration_point[0];
        local_coordinates[1] = r_integration_point[1];

        // Triangle3D3 reference weights sum to 1/2.  The physical
        // determinant is 2*area, hence these weights sum exactly to area.
        IntegrationPoint<3> integration_point(
            local_coordinates,
            r_integration_point.Weight() * 2.0 * area);

        Matrix shape_functions(1, 3);
        shape_functions(0, 0) = 1.0 - local_coordinates[0] - local_coordinates[1];
        shape_functions(0, 1) = local_coordinates[0];
        shape_functions(0, 2) = local_coordinates[1];
        Matrix shape_function_gradients(3, 2);
        shape_function_gradients(0, 0) = -1.0;
        shape_function_gradients(0, 1) = -1.0;
        shape_function_gradients(1, 0) = 1.0;
        shape_function_gradients(1, 1) = 0.0;
        shape_function_gradients(2, 0) = 0.0;
        shape_function_gradients(2, 1) = 1.0;

        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
            SelectTetraIntegrationMethod(IntegrationOrder),
            integration_point,
            shape_functions,
            shape_function_gradients);

        result.push_back(CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePoint(
            3,
            2,
            data_container,
            points));
    }
    return result;
}

GeometryPointerType CreateVolumeQuadraturePointGeometry(
    const array_1d<double, 3>& rPoint,
    const double Weight)
{
    PointerVector<NodeType> points;
    points.push_back(NodePointerType(new NodeType(0, rPoint[0], rPoint[1], rPoint[2])));
    array_1d<double, 3> local_coordinates = ZeroVector(3);
    IntegrationPoint<3> integration_point(local_coordinates, Weight);
    Matrix shape_functions(1, 1);
    shape_functions(0, 0) = 1.0;
    Matrix shape_function_gradients(1, 3);
    noalias(shape_function_gradients) = ZeroMatrix(1, 3);
    GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
        GeometryData::IntegrationMethod::GI_GAUSS_1,
        integration_point,
        shape_functions,
        shape_function_gradients);
    return CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePoint(
        3,
        3,
        data_container,
        points);
}

GeometryData::IntegrationMethod SelectTetraIntegrationMethod(const std::size_t Order)
{
    switch (std::min<std::size_t>(std::max<std::size_t>(Order, 1), 5)) {
        case 1: return GeometryData::IntegrationMethod::GI_GAUSS_1;
        case 2: return GeometryData::IntegrationMethod::GI_GAUSS_2;
        case 3: return GeometryData::IntegrationMethod::GI_GAUSS_3;
        case 4: return GeometryData::IntegrationMethod::GI_GAUSS_4;
        default: return GeometryData::IntegrationMethod::GI_GAUSS_5;
    }
}

std::vector<GeometryPointerType> CreateNativeTetraQuadraturePointGeometries(
    const std::array<NodePointerType, 4>& rNodes,
    const std::size_t IntegrationOrder)
{
    auto p_tetrahedron = Kratos::make_shared<Tetrahedra3D4<NodeType>>(
        rNodes[0], rNodes[1], rNodes[2], rNodes[3]);
    const auto& integration_points =
        p_tetrahedron->IntegrationPoints(SelectTetraIntegrationMethod(IntegrationOrder));
    std::vector<GeometryPointerType> result;
    result.reserve(integration_points.size());
    Matrix jacobian;
    for (const auto& integration_point : integration_points) {
        array_1d<double, 3> local_coordinates = ZeroVector(3);
        local_coordinates[0] = integration_point[0];
        local_coordinates[1] = integration_point[1];
        local_coordinates[2] = integration_point[2];
        array_1d<double, 3> global_coordinates = ZeroVector(3);
        p_tetrahedron->GlobalCoordinates(global_coordinates, local_coordinates);
        p_tetrahedron->Jacobian(jacobian, local_coordinates);
        const double weight = integration_point.Weight() * std::abs(MathUtils<double>::Det(jacobian));
        if (weight > std::numeric_limits<double>::min()) {
            result.push_back(CreateVolumeQuadraturePointGeometry(global_coordinates, weight));
        }
    }
    return result;
}

void SetCommonConditionData(
    Condition& rCondition,
    const Vector& rKnotSpanSizes,
    const std::vector<GeometryPointerType>& rNeighbourGeometries,
    const double CharacteristicLength)
{
    rCondition.SetValue(KNOT_SPAN_SIZES, rKnotSpanSizes);
    rCondition.SetValue(NEIGHBOUR_GEOMETRIES, rNeighbourGeometries);
    array_1d<double, 3> characteristic_length = ZeroVector(3);
    characteristic_length[0] = CharacteristicLength;
    rCondition.SetValue(CHARACTERISTIC_GEOMETRY_LENGTH, characteristic_length);
}

} // unnamed namespace

SnakeGapSbm3DLevelSetClippingUtilities::Counters
SnakeGapSbm3DLevelSetClippingUtilities::Create(
    SnakeGapSbm3DUtilities& rGapUtilities,
    ModelPart& rRootModelPart,
    ModelPart& rSkinSubModelPart,
    const ModelPart& rSurrogateSubModelPart,
    const SnakeGapSbm3DUtilities::KnotSpanGridInfo& rGridInfo,
    const std::set<SnakeGapSbm3DUtilities::SpanKey3D>& rActiveBackgroundSpans,
    const Vector& rKnotSpanSizes,
    ModelPart& rGapElementsModelPart,
    ModelPart& rGapConditionsModelPart,
    ModelPart& rGapInterfacesModelPart,
    Properties::Pointer pProperties,
    const Settings& rSettings) const
{
    KRATOS_ERROR_IF_NOT(pProperties)
        << "Level-set clipping requires valid element/condition properties.\n";
    KRATOS_ERROR_IF(rSkinSubModelPart.NumberOfConditions() == 0)
        << "Level-set clipping requires a non-empty skin model part.\n";

    Counters counters;
    std::vector<SkinTriangle> skin_triangles;
    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        const auto& r_geometry = r_condition.GetGeometry();
        for (IndexType i = 1; i + 1 < r_geometry.PointsNumber(); ++i) {
            SkinTriangle triangle{{r_geometry.pGetPoint(0), r_geometry.pGetPoint(i), r_geometry.pGetPoint(i + 1)}};
            triangle.pSourceCondition = &r_condition;
            const auto edge01 = Coordinates(triangle.Nodes[1]) - Coordinates(triangle.Nodes[0]);
            const auto edge02 = Coordinates(triangle.Nodes[2]) - Coordinates(triangle.Nodes[0]);
            MathUtils<double>::CrossProduct(triangle.Normal, edge01, edge02);
            const double norm = norm_2(triangle.Normal);
            if (norm <= RelativeTolerance * std::max(SquaredNorm(edge01), SquaredNorm(edge02))) {
                continue;
            }
            triangle.Normal /= norm;
            skin_triangles.push_back(std::move(triangle));
        }
    }
    KRATOS_ERROR_IF(skin_triangles.empty())
        << "Level-set clipping found no non-degenerate skin triangles.\n";

    GeometricalObjectsBins skin_condition_bins(rSkinSubModelPart.Conditions());
    std::unordered_map<const Condition*, std::size_t> triangle_by_condition;
    for (std::size_t i = 0; i < skin_triangles.size(); ++i) {
        // A linear surface condition contributes exactly one triangle. Keep a
        // conservative fallback to the exhaustive search for other geometry
        // types, where a condition can contribute more than one fan triangle.
        const auto* p_condition = skin_triangles[i].pSourceCondition;
        if (p_condition != nullptr && p_condition->GetGeometry().PointsNumber() == 3) {
            triangle_by_condition.emplace(p_condition, i);
        }
    }

    // The winding of imported STL triangles is not a reliable definition of
    // the Gap-SBM exterior.  Match the orientation convention of the legacy
    // workflow: outer skin normals agree with the nearest surrogate-face
    // normal; inner skin normals oppose it.  Without this, KeepNegativePhi
    // can retain the volume outside the true skin.
    std::vector<SurrogateFaceOrientation> surrogate_face_orientations;
    surrogate_face_orientations.reserve(rSurrogateSubModelPart.NumberOfConditions());
    for (const auto& r_condition : rSurrogateSubModelPart.Conditions()) {
        KRATOS_ERROR_IF_NOT(r_condition.Has(NORMAL))
            << "Level-set clipping requires NORMAL on surrogate condition #"
            << r_condition.Id() << ".\n";
        const Vector& r_normal_vector = r_condition.GetValue(NORMAL);
        KRATOS_ERROR_IF(r_normal_vector.size() < 3)
            << "Invalid NORMAL on surrogate condition #" << r_condition.Id() << ".\n";
        SurrogateFaceOrientation orientation;
        orientation.Center = r_condition.GetGeometry().Center();
        orientation.Normal[0] = r_normal_vector[0];
        orientation.Normal[1] = r_normal_vector[1];
        orientation.Normal[2] = r_normal_vector[2];
        const double normal_norm = norm_2(orientation.Normal);
        KRATOS_ERROR_IF(normal_norm <= RelativeTolerance)
            << "Zero NORMAL on surrogate condition #" << r_condition.Id() << ".\n";
        orientation.Normal /= normal_norm;
        surrogate_face_orientations.push_back(orientation);
    }
    KRATOS_ERROR_IF(surrogate_face_orientations.empty())
        << "Level-set clipping requires surrogate face orientations.\n";

    for (auto& r_triangle : skin_triangles) {
        array_1d<double, 3> triangle_center = ZeroVector(3);
        for (const auto& p_node : r_triangle.Nodes) {
            triangle_center += Coordinates(p_node);
        }
        triangle_center /= 3.0;

        const SurrogateFaceOrientation* p_nearest_orientation = nullptr;
        double closest_distance_squared = std::numeric_limits<double>::max();
        for (const auto& r_orientation : surrogate_face_orientations) {
            const double distance_squared = SquaredNorm(
                triangle_center - r_orientation.Center);
            if (distance_squared < closest_distance_squared) {
                closest_distance_squared = distance_squared;
                p_nearest_orientation = &r_orientation;
            }
        }
        KRATOS_ERROR_IF_NOT(p_nearest_orientation)
            << "Level-set clipping could not orient a skin triangle.\n";

        const bool must_reverse = rSettings.KeepNegativePhi
            ? inner_prod(r_triangle.Normal, p_nearest_orientation->Normal) < 0.0
            : inner_prod(r_triangle.Normal, p_nearest_orientation->Normal) > 0.0;
        if (must_reverse) {
            r_triangle.Normal *= -1.0;
        }
    }

    IndexType next_projection_node_id = 1;
    const ModelPart& r_skin_root_model_part = rSkinSubModelPart.GetRootModelPart();
    for (const auto& r_node : r_skin_root_model_part.Nodes()) {
        next_projection_node_id = std::max(
            next_projection_node_id,
            r_node.Id() + 1);
    }

    const auto create_projection_node = [&rSkinSubModelPart,
                                         &next_projection_node_id](
        const ProjectionResult& rProjection) {
        auto p_projection_node = rSkinSubModelPart.CreateNewNode(
            next_projection_node_id++,
            rProjection.Point[0],
            rProjection.Point[1],
            rProjection.Point[2]);
        p_projection_node->SetValue(NORMAL, rProjection.Normal);
        p_projection_node->SetValue(IDENTIFIER, "LEVEL_SET_CLIPPING_PROJECTION");
        return p_projection_node;
    };

    // The neighbour is a
    // BREP/NURBS quadrature-point geometry at the centre of a surrogate face,
    // not the (linear) surrogate condition geometry.  In particular this is
    // what gives elements and conditions the correct displacement DOFs.
    const auto surrogate_neighbour_geometries =
        rGapUtilities.CreateSurrogateFaceCentreNeighbourGeometries(
            rSurrogateSubModelPart);

    GridVertexRegistry grid_vertex_registry;
    CutEdgeRegistry cut_edge_registry;
    SubcellFaceRegistry subcell_face_registry;
    std::set<FaceKey> cut_facet_keys;
    std::unordered_map<const GeometryType*, std::vector<GeometryPointerType>> volume_quadrature_by_neighbour;
    std::unordered_map<const GeometryType*, GeometryPointerType> neighbour_geometry_by_key;
    std::unordered_map<const GeometryType*, std::size_t> clipped_tetrahedra_per_neighbour;

    const double characteristic_length = std::max({
        rGridInfo.SpanSizeX,
        rGridInfo.SpanSizeY,
        rGridInfo.SpanSizeZ});
    const double distance_tolerance = RelativeTolerance * characteristic_length;
    const double volume_tolerance = RelativeTolerance * characteristic_length * characteristic_length * characteristic_length;
    const int subdivisions = static_cast<int>(std::max<std::size_t>(rSettings.SubdivisionsPerSpan, 1));

    auto get_grid_vertex = [&](const int I, const int J, const int K) {
        const std::array<int, 3> key = {{I, J, K}};
        const auto it = grid_vertex_registry.find(key);
        if (it != grid_vertex_registry.end()) {
            return it->second;
        }
        const double x = rGridInfo.MinU + static_cast<double>(I) * rGridInfo.SpanSizeX / subdivisions;
        const double y = rGridInfo.MinV + static_cast<double>(J) * rGridInfo.SpanSizeY / subdivisions;
        const double z = rGridInfo.MinW + static_cast<double>(K) * rGridInfo.SpanSizeZ / subdivisions;
        auto p_node = NodePointerType(new NodeType(0, x, y, z));
        grid_vertex_registry.emplace(key, p_node);
        return p_node;
    };

    auto signed_distance = [&](const NodePointerType& pNode) {
        const auto projection = ProjectToSkin(
            Coordinates(pNode), skin_triangles, &skin_condition_bins, &triangle_by_condition,
            characteristic_length);
        const double phi = inner_prod(Coordinates(pNode) - projection.Point, projection.Normal);
        return std::abs(phi) <= distance_tolerance ? 0.0 : phi;
    };

    const auto is_kept = [&rSettings, distance_tolerance](const double Phi) {
        return rSettings.KeepNegativePhi ? Phi <= distance_tolerance : Phi >= -distance_tolerance;
    };

    auto get_intersection = [&](const Vertex& rVertex0, const Vertex& rVertex1) {
        if (std::abs(rVertex0.Phi) <= distance_tolerance) {
            return rVertex0.pNode;
        }
        if (std::abs(rVertex1.Phi) <= distance_tolerance) {
            return rVertex1.pNode;
        }
        const EdgeKey edge_key(rVertex0.pNode, rVertex1.pNode);
        const auto existing = cut_edge_registry.find(edge_key);
        if (existing != cut_edge_registry.end()) {
            return existing->second;
        }
        const double interpolation = rVertex0.Phi / (rVertex0.Phi - rVertex1.Phi);
        array_1d<double, 3> point = Coordinates(rVertex0.pNode);
        point += interpolation *
            (Coordinates(rVertex1.pNode) - Coordinates(rVertex0.pNode));
        auto p_node = NodePointerType(new NodeType(0, point[0], point[1], point[2]));
        cut_edge_registry.emplace(edge_key, p_node);
        return p_node;
    };

    auto clip_face = [&](const std::array<Vertex, 3>& rFace) {
        std::vector<Vertex> result;
        Vertex start = rFace.back();
        bool start_is_kept = is_kept(start.Phi);
        for (const auto& end : rFace) {
            const bool end_is_kept = is_kept(end.Phi);
            if (start_is_kept != end_is_kept) {
                result.push_back({get_intersection(start, end), 0.0});
            }
            if (end_is_kept) {
                result.push_back(end);
            }
            start = end;
            start_is_kept = end_is_kept;
        }
        std::vector<NodePointerType> nodes;
        for (const auto& r_vertex : result) {
            AddUniqueNode(nodes, r_vertex.pNode);
        }
        return nodes;
    };

    auto register_tetra_faces = [&](const std::array<NodePointerType, 4>& rTetra,
                                     const GeometryPointerType& pNeighbour) {
        array_1d<double, 3> subcell_center = ZeroVector(3);
        for (const auto& p_node : rTetra) {
            subcell_center += Coordinates(p_node);
        }
        subcell_center /= 4.0;
        constexpr std::array<std::array<int, 3>, 4> local_faces = {{
            {{0, 2, 1}}, {{0, 1, 3}}, {{1, 2, 3}}, {{2, 0, 3}}}};
        for (const auto& r_local_face : local_faces) {
            const std::array<NodePointerType, 3> face_nodes = {{
                rTetra[r_local_face[0]], rTetra[r_local_face[1]], rTetra[r_local_face[2]]}};
            const FaceKey face_key(face_nodes[0], face_nodes[1], face_nodes[2]);
            subcell_face_registry[face_key].push_back({
                face_nodes,
                pNeighbour,
                subcell_center,
                cut_facet_keys.find(face_key) != cut_facet_keys.end()});
        }
    };

    constexpr std::array<std::array<int, 4>, 6> hex_tetrahedra = {{
        {{0, 1, 2, 6}}, {{0, 2, 3, 6}}, {{0, 3, 7, 6}},
        {{0, 7, 4, 6}}, {{0, 4, 5, 6}}, {{0, 5, 1, 6}}}};
    constexpr std::array<std::array<int, 3>, 4> tetra_faces = {{
        {{0, 2, 1}}, {{0, 1, 3}}, {{1, 2, 3}}, {{2, 0, 3}}}};
    constexpr std::array<std::array<int, 2>, 6> tetra_edges = {{
        {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}, {{1, 3}}, {{2, 3}}}};

    // Only inactive cells adjacent to the active background region can
    // contribute to the gap.  This is a direct active/inactive-topology test,
    // not a type1/type2/type3 classification, and avoids traversing the full
    // exterior background volume at every refinement level.
    const auto is_valid_span = [&rGridInfo](const int I, const int J, const int K) {
        return I >= 0 && J >= 0 && K >= 0 &&
            I < static_cast<int>(rGridInfo.NumberOfSpansX) &&
            J < static_cast<int>(rGridInfo.NumberOfSpansY) &&
            K < static_cast<int>(rGridInfo.NumberOfSpansZ);
    };
    const auto is_gap_candidate = [&rActiveBackgroundSpans, &is_valid_span](
        const int I, const int J, const int K) {
        if (rActiveBackgroundSpans.find({I, J, K}) != rActiveBackgroundSpans.end()) {
            return false;
        }
        for (int di = -1; di <= 1; ++di) {
            for (int dj = -1; dj <= 1; ++dj) {
                for (int dk = -1; dk <= 1; ++dk) {
                    if (di == 0 && dj == 0 && dk == 0) continue;
                    if (is_valid_span(I + di, J + dj, K + dk) &&
                        rActiveBackgroundSpans.find({I + di, J + dj, K + dk}) !=
                            rActiveBackgroundSpans.end()) return true;
                }
            }
        }
        return false;
    };

    std::vector<GeometryPointerType> all_neighbour_geometries;
    all_neighbour_geometries.reserve(surrogate_neighbour_geometries.size());
    for (const auto& r_entry : surrogate_neighbour_geometries) {
        if (r_entry.second) {
            all_neighbour_geometries.push_back(r_entry.second);
        }
    }
    KRATOS_ERROR_IF(all_neighbour_geometries.empty())
        << "Level-set clipping could not create surrogate NEIGHBOUR_GEOMETRIES.\n";

    const auto select_neighbour_geometry =
        [&all_neighbour_geometries](const array_1d<double, 3>& rPoint) {
            GeometryPointerType p_result;
            double closest_distance_squared = std::numeric_limits<double>::max();
            for (const auto& p_candidate : all_neighbour_geometries) {
                const double distance_squared = SquaredNorm(rPoint - p_candidate->Center());
                if (distance_squared < closest_distance_squared) {
                    closest_distance_squared = distance_squared;
                    p_result = p_candidate;
                }
            }
            return p_result;
        };

    for (int I = 0; I < static_cast<int>(rGridInfo.NumberOfSpansX); ++I) {
        for (int J = 0; J < static_cast<int>(rGridInfo.NumberOfSpansY); ++J) {
            for (int K = 0; K < static_cast<int>(rGridInfo.NumberOfSpansZ); ++K) {
        if (!is_gap_candidate(I, J, K)) {
            continue;
        }
        for (int a = 0; a < subdivisions; ++a) {
            for (int b = 0; b < subdivisions; ++b) {
                for (int c = 0; c < subdivisions; ++c) {
        const int i = I * subdivisions + a;
        const int j = J * subdivisions + b;
        const int k = K * subdivisions + c;
        const std::array<NodePointerType, 8> hex_nodes = {{
            get_grid_vertex(i,     j,     k),     get_grid_vertex(i + 1, j,     k),
            get_grid_vertex(i + 1, j + 1, k),     get_grid_vertex(i,     j + 1, k),
            get_grid_vertex(i,     j,     k + 1), get_grid_vertex(i + 1, j,     k + 1),
            get_grid_vertex(i + 1, j + 1, k + 1), get_grid_vertex(i,     j + 1, k + 1)}};

        std::array<double, 8> hex_phi;
        for (std::size_t i = 0; i < hex_nodes.size(); ++i) {
            hex_phi[i] = signed_distance(hex_nodes[i]);
        }

        ++counters.NumberOfProcessedHexCells;
        for (const auto& r_local_tetra : hex_tetrahedra) {
            ++counters.NumberOfGeneratedTetrahedraBeforeClipping;
            std::array<Vertex, 4> tetra;
            for (std::size_t i = 0; i < tetra.size(); ++i) {
                tetra[i] = {hex_nodes[r_local_tetra[i]], hex_phi[r_local_tetra[i]]};
            }

            std::vector<NodePointerType> clipped_nodes;
            for (const auto& r_vertex : tetra) {
                if (is_kept(r_vertex.Phi)) {
                    AddUniqueNode(clipped_nodes, r_vertex.pNode);
                }
            }
            // Include vertices that lie on the zero level as well as strict
            // edge crossings.  The latter alone misses a valid cut facet
            // whenever phi=0 passes through a background-grid vertex/edge.
            std::vector<NodePointerType> cut_nodes;
            for (const auto& r_vertex : tetra) {
                if (std::abs(r_vertex.Phi) <= distance_tolerance) {
                    AddUniqueNode(cut_nodes, r_vertex.pNode);
                }
            }
            for (const auto& r_edge : tetra_edges) {
                const auto& r_vertex0 = tetra[r_edge[0]];
                const auto& r_vertex1 = tetra[r_edge[1]];
                if ((r_vertex0.Phi < -distance_tolerance && r_vertex1.Phi > distance_tolerance) ||
                    (r_vertex0.Phi > distance_tolerance && r_vertex1.Phi < -distance_tolerance)) {
                    const auto p_cut_node = get_intersection(r_vertex0, r_vertex1);
                    AddUniqueNode(clipped_nodes, p_cut_node);
                    AddUniqueNode(cut_nodes, p_cut_node);
                }
            }
            if (clipped_nodes.size() < 4) {
                continue;
            }

            std::vector<std::pair<std::vector<NodePointerType>, bool>> polyhedron_faces;
            for (const auto& r_local_face : tetra_faces) {
                const std::array<Vertex, 3> face = {{
                    tetra[r_local_face[0]], tetra[r_local_face[1]], tetra[r_local_face[2]]}};
                auto face_nodes = clip_face(face);
                if (face_nodes.size() >= 3) {
                    RotateToCanonicalFirstNode(face_nodes);
                    polyhedron_faces.emplace_back(std::move(face_nodes), false);
                }
            }
            if (cut_nodes.size() >= 3) {
                const auto cut_center = [&cut_nodes]() {
                    array_1d<double, 3> center = ZeroVector(3);
                    for (const auto& p_node : cut_nodes) {
                        center += Coordinates(p_node);
                    }
                    center /= static_cast<double>(cut_nodes.size());
                    return center;
                }();
                const auto cut_projection = ProjectToSkin(
                    cut_center, skin_triangles, &skin_condition_bins, &triangle_by_condition,
                    characteristic_length);
                array_1d<double, 3> desired_normal = cut_projection.Normal;
                if (!rSettings.KeepNegativePhi) {
                    desired_normal *= -1.0;
                }
                array_1d<double, 3> reference = std::abs(desired_normal[0]) < 0.8
                    ? array_1d<double, 3>{{1.0, 0.0, 0.0}}
                    : array_1d<double, 3>{{0.0, 1.0, 0.0}};
                array_1d<double, 3> tangent0 = ZeroVector(3);
                MathUtils<double>::CrossProduct(tangent0, desired_normal, reference);
                tangent0 /= norm_2(tangent0);
                array_1d<double, 3> tangent1 = ZeroVector(3);
                MathUtils<double>::CrossProduct(tangent1, desired_normal, tangent0);
                std::sort(cut_nodes.begin(), cut_nodes.end(), [&cut_center, &tangent0, &tangent1](const auto& pNode0, const auto& pNode1) {
                    const auto vector0 = Coordinates(pNode0) - cut_center;
                    const auto vector1 = Coordinates(pNode1) - cut_center;
                    return std::atan2(inner_prod(vector0, tangent1), inner_prod(vector0, tangent0)) <
                           std::atan2(inner_prod(vector1, tangent1), inner_prod(vector1, tangent0));
                });
                const auto face_normal = [&cut_nodes]() {
                    array_1d<double, 3> normal = ZeroVector(3);
                    MathUtils<double>::CrossProduct(
                        normal,
                        Coordinates(cut_nodes[1]) - Coordinates(cut_nodes[0]),
                        Coordinates(cut_nodes[2]) - Coordinates(cut_nodes[0]));
                    return normal;
                }();
                if (inner_prod(face_normal, desired_normal) < 0.0) {
                    std::reverse(cut_nodes.begin() + 1, cut_nodes.end());
                }
                RotateToCanonicalFirstNode(cut_nodes);
                polyhedron_faces.emplace_back(cut_nodes, true);
            }

            array_1d<double, 3> centroid = ZeroVector(3);
            for (const auto& p_node : clipped_nodes) {
                centroid += Coordinates(p_node);
            }
            centroid /= static_cast<double>(clipped_nodes.size());
            const auto p_centroid = NodePointerType(new NodeType(0, centroid[0], centroid[1], centroid[2]));
            // One parent extension point per clipped original tetrahedron.
            // All fan tetrahedra are only an integration tessellation of the
            // same clipped domain and must not create artificial changes of
            // NEIGHBOUR_GEOMETRY (and hence artificial interfaces).
            const auto p_clipped_tetra_neighbour = select_neighbour_geometry(centroid);
            KRATOS_ERROR_IF_NOT(p_clipped_tetra_neighbour)
                << "Level-set clipping failed to select a parent "
                << "NEIGHBOUR_GEOMETRY for a clipped tetrahedron.\n";
            neighbour_geometry_by_key.emplace(
                p_clipped_tetra_neighbour.get(),
                p_clipped_tetra_neighbour);
            ++clipped_tetrahedra_per_neighbour[p_clipped_tetra_neighbour.get()];

            double clipped_polyhedron_volume = 0.0;
            double generated_tetra_volume = 0.0;

            for (const auto& r_face_data : polyhedron_faces) {
                auto face_nodes = r_face_data.first;
                const bool is_cut_facet = r_face_data.second;
                RotateToCanonicalFirstNode(face_nodes);
                for (std::size_t i = 1; i + 1 < face_nodes.size(); ++i) {
                    std::array<NodePointerType, 3> facet = {{
                        face_nodes[0], face_nodes[i], face_nodes[i + 1]}};
                    const FaceKey facet_key(facet[0], facet[1], facet[2]);
                    if (is_cut_facet) {
                        cut_facet_keys.insert(facet_key);
                    }
                    const std::array<NodePointerType, 4> subcell = {{
                        p_centroid, facet[0], facet[1], facet[2]}};
                    const double subcell_volume =
                        TetraVolume(subcell[0], subcell[1], subcell[2], subcell[3]);
                    if (subcell_volume <= volume_tolerance) {
                        continue;
                    }

                    // Independent polyhedron-volume audit.  Orient the
                    // boundary triangle away from the interior centroid and
                    // use the signed pyramid volume.  This detects missing
                    // or duplicated face triangles in the fan construction.
                    array_1d<double, 3> edge01 = Coordinates(facet[1]) - Coordinates(facet[0]);
                    array_1d<double, 3> edge02 = Coordinates(facet[2]) - Coordinates(facet[0]);
                    array_1d<double, 3> facet_normal = ZeroVector(3);
                    MathUtils<double>::CrossProduct(facet_normal, edge01, edge02);
                    const array_1d<double, 3> facet_center =
                        (Coordinates(facet[0]) + Coordinates(facet[1]) + Coordinates(facet[2])) / 3.0;
                    if (inner_prod(facet_normal, facet_center - centroid) < 0.0) {
                        std::swap(facet[1], facet[2]);
                        edge01 = Coordinates(facet[1]) - Coordinates(facet[0]);
                        edge02 = Coordinates(facet[2]) - Coordinates(facet[0]);
                        MathUtils<double>::CrossProduct(facet_normal, edge01, edge02);
                    }
                    clipped_polyhedron_volume += inner_prod(
                        Coordinates(facet[0]) - centroid,
                        facet_normal) / 6.0;
                    generated_tetra_volume += subcell_volume;

                    auto quadrature_geometries =
                        CreateNativeTetraQuadraturePointGeometries(
                            subcell,
                            rSettings.VolumeIntegrationOrder);
                    auto& r_batch = volume_quadrature_by_neighbour[p_clipped_tetra_neighbour.get()];
                    for (const auto& p_quadrature_geometry : quadrature_geometries) {
                        r_batch.push_back(p_quadrature_geometry);
                    }
                    register_tetra_faces(subcell, p_clipped_tetra_neighbour);
                    ++counters.NumberOfClippedVolumeSubcells;
                }
            }

            counters.TotalClippedPolyhedronVolume += clipped_polyhedron_volume;
            counters.TotalGeneratedTetraVolume += generated_tetra_volume;
            const double audit_tolerance = 1.0e-10 * std::max(
                {1.0, clipped_polyhedron_volume, generated_tetra_volume});
            if (std::abs(clipped_polyhedron_volume - generated_tetra_volume) >
                audit_tolerance) {
                ++counters.NumberOfVolumeAuditFailures;
            }
        }
                }
            }
        }
            }
        }
    }

    IndexType next_element_id = 1;
    for (const auto& r_element : rRootModelPart.Elements()) {
        next_element_id = std::max(next_element_id, r_element.Id() + 1);
    }
    const Element& r_reference_element = KratosComponents<Element>::Get(rSettings.VolumeElementName);
    for (auto& r_batch_entry : volume_quadrature_by_neighbour) {
        auto& r_quadrature_geometries = r_batch_entry.second;
        if (r_quadrature_geometries.empty()) {
            continue;
        }
        std::vector<GeometryPointerType> geometries;
        geometries.swap(r_quadrature_geometries);
        auto p_coupling_geometry = Kratos::make_shared<CouplingGeometry<NodeType>>(std::move(geometries));
        auto p_element = r_reference_element.Create(next_element_id++, p_coupling_geometry, pProperties);
        const auto p_neighbour = neighbour_geometry_by_key.at(r_batch_entry.first);
        p_element->SetValue(NEIGHBOUR_GEOMETRIES, NeighbourGeometriesVectorType{p_neighbour});

        // CouplingGeometry compacts its geometry parts after initialization.
        // Preserve the physical QP positions explicitly for diagnostics and
        // post-processing; these are the exact points used by the volume rule.
        std::vector<Vector> integration_points;
        integration_points.reserve(p_coupling_geometry->NumberOfGeometryParts());
        for (std::size_t point_index = 0;
             point_index < p_coupling_geometry->NumberOfGeometryParts();
             ++point_index) {
            const auto point = p_coupling_geometry->GetGeometryPart(point_index).Center();
            Vector coordinates(3);
            coordinates[0] = point[0];
            coordinates[1] = point[1];
            coordinates[2] = point[2];
            integration_points.push_back(std::move(coordinates));
        }
        p_element->SetValue(INTEGRATION_POINTS, integration_points);
        rGapElementsModelPart.AddElement(p_element);
        ++counters.NumberOfQuadraturePointElementsCreated;
    }

    IndexType next_condition_id = 1;
    for (const auto& r_condition : rRootModelPart.Conditions()) {
        next_condition_id = std::max(next_condition_id, r_condition.Id() + 1);
    }
    const Condition& r_boundary_condition = KratosComponents<Condition>::Get(rSettings.BoundaryConditionName);
    const Condition& r_interface_condition = KratosComponents<Condition>::Get(rSettings.InterfaceConditionName);

    for (const auto& r_face_entry : subcell_face_registry) {
        const auto& r_occurrences = r_face_entry.second;
        if (r_occurrences.size() > 2) {
            ++counters.NumberOfNonManifoldFaces;
            continue;
        }
        if (r_occurrences.size() == 1) {
            const auto& r_occurrence = r_occurrences.front();
            if (!r_occurrence.IsCutFacet) {
                continue;
            }
            const auto quadrature_geometries =
                CreateFacetQuadraturePointGeometries(
                    r_occurrence.Nodes, rSettings.BoundaryIntegrationOrder);
            for (const auto& p_quadrature_geometry : quadrature_geometries) {
                auto p_condition = r_boundary_condition.Create(
                    next_condition_id++,
                    p_quadrature_geometry,
                    pProperties);
                SetCommonConditionData(
                    *p_condition,
                    rKnotSpanSizes,
                    NeighbourGeometriesVectorType{r_occurrence.pNeighbourGeometry},
                    characteristic_length);
                const auto projection = ProjectToSkin(
                    p_quadrature_geometry->Center().Coordinates(),
                    skin_triangles,
                    &skin_condition_bins,
                    &triangle_by_condition,
                    characteristic_length);
                const auto p_projection_node = create_projection_node(projection);
                p_condition->SetValue(PROJECTION_NODE, p_projection_node);
                p_condition->SetValue(
                    PROJECTION_NODE_ID,
                    static_cast<int>(p_projection_node->Id()));
                p_condition->SetValue(NORMAL, projection.Normal);
                rGapConditionsModelPart.AddCondition(p_condition);
            }
            ++counters.NumberOfCutBoundaryFacets;
            continue;
        }

        const auto& r_occurrence0 = r_occurrences[0];
        const auto& r_occurrence1 = r_occurrences[1];
        if (!rSettings.CreateInterfaceConditions ||
            r_occurrence0.pNeighbourGeometry.get() == r_occurrence1.pNeighbourGeometry.get()) {
            continue;
        }
        // The first neighbour is the plus/outgoing side.  Orient the facet
        // normal from its owning subcell (occurrence 0) towards the other
        // subcell, independently of surrogate-face centres or insertion order.
        auto interface_nodes = r_occurrence0.Nodes;
        const auto edge01 = Coordinates(interface_nodes[1]) - Coordinates(interface_nodes[0]);
        const auto edge02 = Coordinates(interface_nodes[2]) - Coordinates(interface_nodes[0]);
        array_1d<double, 3> interface_normal = ZeroVector(3);
        MathUtils<double>::CrossProduct(interface_normal, edge01, edge02);
        const array_1d<double, 3> outgoing_direction =
            r_occurrence1.OwnerSubcellCenter - r_occurrence0.OwnerSubcellCenter;
        if (inner_prod(interface_normal, outgoing_direction) < 0.0) {
            std::swap(interface_nodes[1], interface_nodes[2]);
        }
        const auto quadrature_geometries =
            CreateFacetQuadraturePointGeometries(
                interface_nodes, rSettings.BoundaryIntegrationOrder);
        const auto interface_edge01 = Coordinates(interface_nodes[1]) - Coordinates(interface_nodes[0]);
        const auto interface_edge02 = Coordinates(interface_nodes[2]) - Coordinates(interface_nodes[0]);
        array_1d<double, 3> interface_area_normal = ZeroVector(3);
        MathUtils<double>::CrossProduct(
            interface_area_normal, interface_edge01, interface_edge02);
        counters.TotalInterfaceArea += 0.5 * norm_2(interface_area_normal);
        for (const auto& p_quadrature_geometry : quadrature_geometries) {
            auto p_condition = r_interface_condition.Create(
                next_condition_id++,
                p_quadrature_geometry,
                pProperties);
            SetCommonConditionData(
                *p_condition,
                rKnotSpanSizes,
                NeighbourGeometriesVectorType{
                    r_occurrence0.pNeighbourGeometry,
                    r_occurrence1.pNeighbourGeometry},
                characteristic_length);
            rGapInterfacesModelPart.AddCondition(p_condition);
        }
        ++counters.NumberOfInterfaceFaces;
    }

    KRATOS_ERROR_IF(counters.NumberOfNonManifoldFaces > 0)
        << "Level-set clipping generated " << counters.NumberOfNonManifoldFaces
        << " non-manifold subcell faces.\n";

    KRATOS_ERROR_IF(counters.NumberOfVolumeAuditFailures > 0)
        << "Level-set clipping volume audit failed for "
        << counters.NumberOfVolumeAuditFailures
        << " clipped tetrahedra. Expected clipped-polyhedron volume="
        << counters.TotalClippedPolyhedronVolume
        << ", generated tetra volume="
        << counters.TotalGeneratedTetraVolume << ".\n";

    counters.NumberOfDistinctNeighbourGeometries = clipped_tetrahedra_per_neighbour.size();
    if (!clipped_tetrahedra_per_neighbour.empty()) {
        counters.MinimumClippedTetrahedraPerNeighbour = std::numeric_limits<std::size_t>::max();
        for (const auto& r_entry : clipped_tetrahedra_per_neighbour) {
            counters.MinimumClippedTetrahedraPerNeighbour = std::min(
                counters.MinimumClippedTetrahedraPerNeighbour, r_entry.second);
            counters.MaximumClippedTetrahedraPerNeighbour = std::max(
                counters.MaximumClippedTetrahedraPerNeighbour, r_entry.second);
        }
    }

    KRATOS_INFO_IF("SnakeGapSbm3DLevelSetClipping", rSettings.EchoLevel > 0)
        << "3D Gap-SBM level-set clipping summary:\n"
        << "  processed hex cells: " << counters.NumberOfProcessedHexCells << "\n"
        << "  tetrahedra before clipping: " << counters.NumberOfGeneratedTetrahedraBeforeClipping << "\n"
        << "  clipped volume subcells: " << counters.NumberOfClippedVolumeSubcells << "\n"
        << "  cut boundary facets: " << counters.NumberOfCutBoundaryFacets << "\n"
        << "  interface faces: " << counters.NumberOfInterfaceFaces << "\n"
        << "  non-manifold faces: " << counters.NumberOfNonManifoldFaces << "\n"
        << "  quadrature-point elements: " << counters.NumberOfQuadraturePointElementsCreated << "\n"
        << "  clipped-polyhedron volume: " << counters.TotalClippedPolyhedronVolume << "\n"
        << "  generated tetra volume: " << counters.TotalGeneratedTetraVolume << "\n"
        << "  volume-audit failures: " << counters.NumberOfVolumeAuditFailures << "\n"
        << "  distinct neighbour geometries: " << counters.NumberOfDistinctNeighbourGeometries << "\n"
        << "  clipped tetrahedra per neighbour [min,max]: ["
        << counters.MinimumClippedTetrahedraPerNeighbour << ", "
        << counters.MaximumClippedTetrahedraPerNeighbour << "]\n"
        << "  total interface area: " << counters.TotalInterfaceArea << "\n";

    return counters;
}

} // namespace Kratos
