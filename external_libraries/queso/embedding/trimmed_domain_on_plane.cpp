//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

//// STL includes

//// Project includes
#include "queso/embedding/trimmed_domain_on_plane.h"
#include "queso/embedding/brep_operator.h"
#include "queso/embedding/trimmed_domain.h"
#include "queso/embedding/polygon.h"

namespace queso
{

typedef TrimmedDomainOnPlane::Point2DType Point2DType;
typedef TrimmedDomainOnPlane::TriangleMeshPtrType TriangleMeshPtrType;
typedef TrimmedDomainOnPlane::Edge2D Edge2D;
typedef TrimmedDomainOnPlane::Point2DType Point2DType;
typedef TrimmedDomainOnPlane::Point2DSetType Point2DSetType;

void TrimmedDomainOnPlane::CollectEdgesOnPlane(const TriangleMeshInterface &rTriangleMesh)
{
    const auto& edges_on_planes = rTriangleMesh.GetEdgesOnPlanes();
    const IndexType plane_index = DIRINDEX3*2UL + static_cast<IndexType>(mUpperBoundary);

    const auto& r_vertices = rTriangleMesh.GetVertices();
    const auto& edges_on_plane = edges_on_planes[plane_index];
    // Should only require 2UL*edges_on_plane.size(), however small buffer is used.
    // Note, initial capacity must be large enough. New allocation is not allowed and will crash.
    Reserve( std::max<IndexType>(3UL*edges_on_plane.size(), 10UL) );
    for( const auto& edge : edges_on_plane ){
        const IndexType vertex_index_1 = std::get<0>(edge);
        const IndexType vertex_index_2 = std::get<1>(edge);
        const IndexType normal_index = std::get<2>(edge);

        const auto& P1 = r_vertices[vertex_index_1];
        const auto& P2 = r_vertices[vertex_index_2];
        const auto& normal = rTriangleMesh.Normal(normal_index);
        InsertEdge(P1, P2, normal);
    }
}

void TrimmedDomainOnPlane::CloseContourEdges(const BRepOperator* pOperator) {
    // 1. Insert edges such that positive oriented edges span the entire AABB.
    // Find all intersections of Edges with upper bound of AABB.
    //     UP2 x1-------x2---x3
    //         | out   /     |      Point  x2 will be found.
    //         |------/      |      Points x1 and x3 are also added.
    //     LB2 |  inside     |
    //        LB1           UP1
    // LB1 - lower bound of AABB in DIRINDEX1
    // UP1 - upper bound of AABB in DIRINDEX1
    // We store only the value in DIRINDEX1 in intersected_vertices, since the position in
    // DIRINDEX2 is always mUpperBound[DIRINDEX2].
    // We also store if vertex is left (true) or right (false) boundary,
    // x1-left, x2-left, x3-right.
    // Only, x2 and x3 are connected with one edge, since x1 and x2 are both left boundaries.
    // Get all edged that intersect upper boundary
    // These are edges that have at least one vertex located on the upper boundary: mUpperBoundary[DIRINDEX2]).
    std::vector<Edge2D> intersect_edges_positive{};
    std::vector<Edge2D> intersect_edges_negative{};
    std::vector<Edge2D> intersect_edges_vertical{};
    FindIntersectingEdgesWithUpperBound(intersect_edges_positive, Orientation::Positive);
    FindIntersectingEdgesWithUpperBound(intersect_edges_negative, Orientation::Negative);
    FindIntersectingEdgesWithUpperBound(intersect_edges_vertical, Orientation::Vertical);

    // Remove dublicates
    RemoveDublicateVerticalEdges(intersect_edges_vertical, mVerticesVertical);
    RemoveDublicateEdges(intersect_edges_positive, mVerticesPositive);
    RemoveDublicateEdges(intersect_edges_negative, mVerticesNegative);

    // double: position in DIRINDEX1, bool = true: left boundary, bool = false: right boundary.
    std::vector<std::pair<double, bool>> intersected_vertices{};

    // Check if corner points need to be inserted.
    Point2DType corner_left = { mLowerBound[DIRINDEX1], mUpperBound[DIRINDEX2] };
    Point2DType corner_right = { mUpperBound[DIRINDEX1], mUpperBound[DIRINDEX2] };
    bool add_corner_left = !PointExists(corner_left, *mVerticesSetPositive)
        && !PointExists(corner_left, *mVerticesSetNegative) && !PointExists(corner_left, *mVerticesSetVertical);
    bool add_corner_right = !PointExists(corner_right, *mVerticesSetPositive)
        && !PointExists(corner_right, *mVerticesSetNegative) && !PointExists(corner_right, *mVerticesSetVertical);

    // Add corner point left if neccesary.
    if( add_corner_left ){
        intersected_vertices.push_back( std::make_pair(mLowerBound[DIRINDEX1], true) );
    }

    // Loop until all intersected edges are visited.
    const IndexType size_positive = intersect_edges_positive.size();
    const IndexType size_negative = intersect_edges_negative.size();
    const IndexType size_vertical = intersect_edges_vertical.size();
    IndexType pos_positive = 0;
    IndexType pos_negative = 0;
    IndexType pos_vertical = 0;
    while( pos_vertical < size_vertical || pos_negative < size_negative || pos_positive < size_positive){
        // Get ptr to edges.
        Edge2D* edge_positive = nullptr;
        Edge2D* edge_negative = nullptr;
        Edge2D* edge_vertical = nullptr;
        if( pos_positive < size_positive )
            edge_positive = &intersect_edges_positive[pos_positive];
        if( pos_negative < size_negative )
            edge_negative = &intersect_edges_negative[pos_negative];
        if( pos_vertical < size_vertical)
            edge_vertical = &intersect_edges_vertical[pos_vertical];

        // IndexType size = intersected_vertices.size();
        // DIRINDEX1-value of current position.
        const double left_bound = mLowerBound[DIRINDEX1];

        // Get distance to current positive edge.
        double distance_pos = MAXD;
        if( edge_positive ){
            const auto status = edge_positive->IsVertexOnUpperBoundary();
            if( status.first ){
                distance_pos = std::abs(mVerticesPositive[edge_positive->V1()][0] - left_bound);
            }
            else {
                distance_pos = std::abs(mVerticesPositive[edge_positive->V2()][0] - left_bound);
            }
        }

        // Get distance to current negative edge.
        double distance_neg = MAXD;
        if( edge_negative ){
            const auto status = edge_negative->IsVertexOnUpperBoundary();
            if( status.first ) {
                distance_neg = std::abs(mVerticesNegative[edge_negative->V1()][0] - left_bound);
            } else {
                distance_neg = std::abs(mVerticesNegative[edge_negative->V2()][0] - left_bound);
            }
        }

        // Get distance to current vertical edge.
        double distance_ver = MAXD;
        if( edge_vertical ){
            const auto status = edge_vertical->IsVertexOnUpperBoundary();
            if( status.first ){
                distance_ver = std::abs(mVerticesVertical[edge_vertical->V1()][0] - left_bound);
            }
            else {
                distance_ver = std::abs(mVerticesVertical[edge_vertical->V2()][0] - left_bound);
            }
        }

        const double min_distance = std::min( {distance_pos, distance_neg, distance_ver} );
        const bool pos_found =  distance_pos < (min_distance + mSnapTolerance);
        const bool neg_found =  distance_neg < (min_distance + mSnapTolerance);
        const bool ver_found = distance_ver < (min_distance + mSnapTolerance);

        const IndexType found_count = static_cast<int>(pos_found) + static_cast<int>(neg_found) + static_cast<int>(ver_found);

        // Ignore if double vertex.
        if( found_count > 1 ) {
            pos_positive += static_cast<int>(pos_found);
            pos_negative += static_cast<int>(neg_found);
            pos_vertical += static_cast<int>(ver_found);
        } // ADD POSITIVE
        else if (pos_found) {
            // Add and increment positive
            AddIntersectedVertexPositive(edge_positive, intersected_vertices);
            ++pos_positive;
        } // ADD Negative
        else if( neg_found ){
            // Add and increment negative
            AddIntersectedVertexNegative(edge_negative, intersected_vertices);
            ++pos_negative;
        } // ADD VERTICAL
        else if( ver_found ){
            // Add and increment vertical
            AddIntersectedVertexVertical(edge_vertical, intersected_vertices);
            ++pos_vertical;
        }
    }

    // Add corner right
    if( add_corner_right ){
        IndexType size = intersected_vertices.size();
        // size > 0, makes sure that intersected_vertices[size-1] does not result in segmentation-fault.
        if( size > 0 && std::abs(intersected_vertices[size-1].first - mUpperBound[DIRINDEX1]) > mSnapTolerance ){
            intersected_vertices.push_back( std::make_pair(mUpperBound[DIRINDEX1], false) );
        } else {
            add_corner_right = false;
        }
    }

    const double y_max_negative = GetMaxDIRINDEX2(mVerticesNegative);
    const double y_max_positive = GetMaxDIRINDEX2(mVerticesPositive);
    const double y_max = std::max( std::max( y_max_negative, y_max_positive), mLowerBound[DIRINDEX2]);

    // Insert edges
    const double plane_position = GetPlanePosition();
    if( add_corner_left && add_corner_right && intersected_vertices.size() == 2 ){
        // Corner case, if only the corner points are part of intersected_vertices.
        // In this case, we have to check if this plane is inside or outside.
        const double v_left = intersected_vertices[0].first;
        const double v_right = intersected_vertices[1].first;
        double center = v_left + 0.5 * (v_right - v_left);
        Point3DType new_point{}; // We need 3D point, to use IsInsideTrimmedDomain().
        new_point[DIRINDEX1] = center;
        new_point[DIRINDEX2] = 0.5*(mUpperBound[DIRINDEX2]+y_max);
        new_point[DIRINDEX3] = plane_position;
        Point2DType normal = {0, 1};
        bool Success = true;
        if ( mpTrimmedDomain->IsInsideTrimmedDomain(new_point, Success) ) {
            InsertEdge({v_left, mUpperBound[DIRINDEX2]}, {v_right, mUpperBound[DIRINDEX2]}, normal, Orientation::Positive);
        }
        if( !Success ){ // Safety Test.
            if( pOperator->IsInside(new_point) ){
                    InsertEdge({v_left, mUpperBound[DIRINDEX2]}, {v_right, mUpperBound[DIRINDEX2]}, normal, Orientation::Positive);
            }
        }
    }
    else if( intersected_vertices.size() > 1 ){
        // Insert edge, if two consecutive intersected_vertices are left (true) and right (false) boundary.
        for (IndexType i = 0; i < intersected_vertices.size() - 1; ++i) {
            const double v_left = intersected_vertices[i].first;
            const double v_right = intersected_vertices[i + 1].first;
            Point2DType normal = {0, 1};
            if( intersected_vertices[i].second && !intersected_vertices[i + 1].second  ){
                InsertEdge({v_left, mUpperBound[DIRINDEX2]}, {v_right, mUpperBound[DIRINDEX2]}, normal, Orientation::Positive);
            }
        }
    }

    // 2. We split all positive edges at all negative oriented vertices and vice versa.
    //                 x--^-----x---^--x   positive oriented edges
    //                 |  |     |   |                                     DIRINDEX2
    //                 |  |     |   |                                         ^
    //            x----v--x-----v---x      negative oriented edges            |---> DIRINDEX1

    // Set split points at negative oriented vertices.
    for (IndexType i = 0; i < mVerticesPositive.size(); ++i) {
        SetSplitPoint(mVerticesPositive[i], Orientation::Negative);
    }

    // Set split points at positive oriented vertices.
    for (IndexType i = 0; i < mVerticesNegative.size(); ++i) {
        SetSplitPoint(mVerticesNegative[i], Orientation::Positive);
    }

    // Split edges
    SplitEdgesAtSplitPoint(Orientation::Negative);
    SplitEdgesAtSplitPoint(Orientation::Positive);
}

TriangleMeshPtrType TrimmedDomainOnPlane::TriangulateDomain() const
{
    Point3DType normal = {0.0, 0.0, 0.0};
    double plane_position = GetPlanePosition();
    if (mUpperBoundary) {
        normal[DIRINDEX3] = 1.0;
    }
    else {
        normal[DIRINDEX3] = -1.0;
    }

    auto orientation_origin = Orientation::Positive;

    auto &r_edges_origin = GetEdges(orientation_origin);

    //Instantiate new mesh ptr.
    auto p_new_mesh = MakeUnique<TriangleMesh>();
    p_new_mesh->Reserve(5 * r_edges_origin.size());

    // Loop over positive oriented edges.
    for (IndexType i = 0; i < r_edges_origin.size(); ++i) {
        const auto &v1_up = V1byEdgeId(i, orientation_origin);
        const auto &v2_up = V2byEdgeId(i, orientation_origin);
        const auto& r_edge = r_edges_origin[i];

        int edge_id_dest = FindNegativePartnerEdge(v1_up, v2_up, r_edge.Normal() );

        bool skip = false;
        if (edge_id_dest > -1) { // Lower edge is found.
            const auto &v1_low = V1byEdgeId(edge_id_dest, Orientation::Negative);
            const auto &v2_low = V2byEdgeId(edge_id_dest, Orientation::Negative);

            // If v1 of lower and upper edge coincide.
            // Add:    /|
            //        / |
            //         \|
            if (std::abs(v1_low[1] - v1_up[1]) < mSnapTolerance) {
                Point3DType tmp_point = {0.0, 0.0, 0.0};
                tmp_point[DIRINDEX3] = plane_position;
                tmp_point[DIRINDEX1] = v1_low[0];
                tmp_point[DIRINDEX2] = v1_low[1];
                IndexType v1 = p_new_mesh->AddVertex(tmp_point);
                tmp_point[DIRINDEX1] = v2_low[0];
                tmp_point[DIRINDEX2] = v2_low[1];
                IndexType v2 = p_new_mesh->AddVertex(tmp_point);
                tmp_point[DIRINDEX1] = v2_up[0];
                tmp_point[DIRINDEX2] = v2_up[1];
                IndexType v3 = p_new_mesh->AddVertex(tmp_point);
                if (mUpperBoundary ^ mSwitchOrientation)
                    p_new_mesh->AddTriangle(Vector3i{v1, v2, v3});
                else
                    p_new_mesh->AddTriangle(Vector3i{v2, v1, v3});
                p_new_mesh->AddNormal(normal);

                skip = true; // Skip polygon construction.
            }

            else if (std::abs(v2_low[1] - v2_up[1]) < mSnapTolerance) {
                /* If v2 of lower and upper edge coincide.
                   Add:  |\
                         | \
                         |/                              */
                Point3DType tmp_point = {0.0, 0.0, 0.0};
                tmp_point[DIRINDEX3] = plane_position;
                tmp_point[DIRINDEX1] = v1_low[0];
                tmp_point[DIRINDEX2] = v1_low[1];
                IndexType v1 = p_new_mesh->AddVertex(tmp_point);
                tmp_point[DIRINDEX1] = v2_low[0];
                tmp_point[DIRINDEX2] = v2_low[1];
                IndexType v2 = p_new_mesh->AddVertex(tmp_point);
                tmp_point[DIRINDEX1] = v1_up[0];
                tmp_point[DIRINDEX2] = v1_up[1];
                IndexType v3 = p_new_mesh->AddVertex(tmp_point);
                if (mUpperBoundary ^ mSwitchOrientation)
                    p_new_mesh->AddTriangle(Vector3i{v1, v2, v3});
                else
                    p_new_mesh->AddTriangle(Vector3i{v2, v1, v3});
                p_new_mesh->AddNormal(normal);

                skip = true; // Skip polygon construction.
            }
        }

        if (!skip) {
            // Counter clock-wise orientation. Get first and last point of polygon.
            std::array<Point3DType, 4> corner_points{};
            corner_points[0][DIRINDEX1] = v1_up[0];
            corner_points[0][DIRINDEX2] = v1_up[1];
            corner_points[0][DIRINDEX3] = plane_position;

            corner_points[3][DIRINDEX1] = v2_up[0];
            corner_points[3][DIRINDEX2] = v2_up[1];
            corner_points[3][DIRINDEX3] = plane_position;

            // Second and third point depends if lower edge is found or not.
            if (edge_id_dest > -1) {
                const auto &v1_low = V1byEdgeId(edge_id_dest, Orientation::Negative);
                const auto &v2_low = V2byEdgeId(edge_id_dest, Orientation::Negative);

                corner_points[1][DIRINDEX1] = v1_low[0];
                corner_points[1][DIRINDEX2] = v1_low[1];
                corner_points[1][DIRINDEX3] = plane_position;

                corner_points[2][DIRINDEX1] = v2_low[0];
                corner_points[2][DIRINDEX2] = v2_low[1];
                corner_points[2][DIRINDEX3] = plane_position;
            }
            else {
                corner_points[1][DIRINDEX1] = v1_up[0];
                corner_points[1][DIRINDEX2] = mLowerBound[DIRINDEX2];
                corner_points[1][DIRINDEX3] = plane_position;

                corner_points[2][DIRINDEX1] = v2_up[0];
                corner_points[2][DIRINDEX2] = mLowerBound[DIRINDEX2];
                corner_points[2][DIRINDEX3] = plane_position;
            }

            // Instantiate new polygon.
            Polygon<4> polygon(normal);
            // Orientation of triangle depends whether we are on upper or lower bound of AABB.
            if (mUpperBoundary ^ mSwitchOrientation) {
                polygon.AddVertex(corner_points[0]);
                polygon.AddVertex(corner_points[1]);
                polygon.AddVertex(corner_points[2]);
                polygon.AddVertex(corner_points[3]);
            }
            else {
                polygon.AddVertex(corner_points[0]);
                polygon.AddVertex(corner_points[3]);
                polygon.AddVertex(corner_points[2]);
                polygon.AddVertex(corner_points[1]);
            }
            polygon.AddToTriangleMesh(*p_new_mesh.get());
        }
    }
    //MeshUtilities::Refine(*p_new_mesh, 10);
    return p_new_mesh;
}

void TrimmedDomainOnPlane::FindIntersectingEdgesWithUpperBound(std::vector<Edge2D>& rEdges, OrientationType Orientation ) {
    auto& r_edges = GetEdges(Orientation);
    for (IndexType edge_id = 0; edge_id < GetNumberEdges(Orientation); ++edge_id) {
        const auto &v1 = V1byEdgeId(edge_id, Orientation);
        const auto &v2 = V2byEdgeId(edge_id, Orientation);
        bool v1_on_edge = std::abs(v1[1] - mUpperBound[DIRINDEX2]) < 10.0*mSnapTolerance;
        bool v2_on_edge = std::abs(v2[1] - mUpperBound[DIRINDEX2]) < 10.0*mSnapTolerance;
        if( v1_on_edge ^ v2_on_edge ){
            auto new_edge = r_edges[edge_id];
            new_edge.SetVerticesOnUpperBoundary(v1_on_edge, v2_on_edge);
            rEdges.push_back( new_edge );
        }
    }

    const auto& r_vertices = GetVertices(Orientation);
    // Sort Edges from left to right (DIRINDEX1)
    std::sort( rEdges.begin(), rEdges.end(), [&r_vertices](const auto& rLHs, const auto& rRHs){
        const auto status_1 = rLHs.IsVertexOnUpperBoundary();
        const auto status_2 = rRHs.IsVertexOnUpperBoundary();

        // Sort by vertex on upper boundary from left to right.
        const IndexType index_v_left = status_1.first ? rLHs.V1() : rLHs.V2();
        const IndexType index_v_right = status_2.first ? rRHs.V1() : rRHs.V2();

        double value_left = r_vertices[index_v_left][0];
        double value_right = r_vertices[index_v_right][0];

        // If two edges share a vertex, sort by center of edges.
        const bool same_vertices = (index_v_left == index_v_right);
        if( same_vertices) {
            value_left = 0.5*(r_vertices[rLHs.V1()][0] + r_vertices[rLHs.V2()][0]);
            value_right = 0.5*(r_vertices[rRHs.V1()][0] + r_vertices[rRHs.V2()][0]);
        }

        return value_left < value_right;
    }  );
}

int TrimmedDomainOnPlane::FindNegativePartnerEdge(const Point2DType &rV1, const Point2DType &rV2, const Point2DType &rNormal) const {
    // Get center
    const Point2DType c_positive = {0.5 * (rV1[0] + rV2[0]), 0.5 * (rV1[1] + rV2[1])};
    double min_distance = MAXD;
    IndexType found_id = -1;

    for (IndexType i = 0; i < GetNumberEdges(Orientation::Negative); ++i) {
        const auto &v1 = V1byEdgeId(i, Orientation::Negative);
        const auto &v2 = V2byEdgeId(i, Orientation::Negative);
        // Check if center is containerd.
        if (c_positive[0] > v1[0] && c_positive[0] < v2[0] ) {
            // Check if left and right vertices (DIRINDEX1) are the same.
            if( std::abs(rV1[0] - v1[0]) < 10*mSnapTolerance && std::abs(rV2[0] - v2[0]) < 10*mSnapTolerance) {
                Point2DType c_negative = {0.5 * (v1[0] + v2[0]), 0.5 * (v1[1] + v2[1])};
                // Vector from center to center.
                Point2DType c_pos_c_neg = {c_negative[0] - c_positive[0], c_negative[1] - c_positive[1]  };
                // Dot produced along normal.
                double value = c_pos_c_neg[0]*rNormal[0] + c_pos_c_neg[1]*rNormal[1];

                if( std::abs(value) <= ZEROTOL ){
                    // This means, partner edge is aligned with initial edge.
                    // Two options are possible.
                    //    | Either this is inside or outside.
                    //    V
                    //   \ |  <- We are currently lokking from here.
                    //    \| x-Point   We need to check if x (test_point) is inside or outside.
                    //    /|           If x is not inside, return.
                    //   / |
                    Point3DType test_point{};
                    test_point[DIRINDEX1] = c_positive[0];
                    test_point[DIRINDEX2] = c_positive[1];
                    double plane_position = GetPlanePosition();
                    test_point[DIRINDEX3] = mUpperBoundary ? plane_position + 100.0*mSnapTolerance : plane_position - 100.0*mSnapTolerance;
                    bool is_inside = mpTrimmedDomain->IsInsideTrimmedDomain(test_point);
                    if ( !is_inside ){
                        const double distance = c_positive[1] - 0.5*(v1[1] + v2[1]);
                        // Distance must be larger than 0.0 and large than alreade found min_distance.
                        if (distance > -ZEROTOL && distance < min_distance) {
                            found_id = i;
                            min_distance = distance;
                        }
                    }
                }
                else if( value < ZEROTOL ){
                    const double distance = c_positive[1] - 0.5*(v1[1] + v2[1]);
                    // Distance must be larger than 0.0 and large than alreade found min_distance.
                    if (distance >= -ZEROTOL && distance < min_distance) {
                        found_id = i;
                        min_distance = distance;
                    }
                }
            }

        }
    }
    return found_id;
}

void TrimmedDomainOnPlane::SetSplitPoint(const Point2DType &rPoint, OrientationType OrientationDest)
{
    // Take bool from Edge.
    double current_distance = 1e15;
    int edge_id = -1;
    Point2DType intersection_point{};
    bool is_positive = OrientationDest == Orientation::Positive;
    OrientationType orientation_origin = is_positive ? Orientation::Negative : Orientation::Positive;

    auto& r_edges_dest = GetEdges(OrientationDest);
    auto& r_edges_origin = GetEdges(orientation_origin);

    for (IndexType i = 0; i < r_edges_origin.size(); ++i) {
        const auto &v1 = V1byEdgeId(i, orientation_origin);
        const auto &v2 = V2byEdgeId(i, orientation_origin);

        double pos1 = rPoint[0];
        double pos2 = 0.0;
        // Keep looser tolerance here.
        if (pos1 > v1[0] - mSnapTolerance && pos1 < v2[0] + mSnapTolerance)
        {
            pos2 = v1[1] + (v2[1] - v1[1]) / (v2[0] - v1[0]) * (rPoint[0] - v1[0]);
            double distance = (!is_positive) ? (rPoint[1] - pos2) : (pos2 - rPoint[1]);
            if ((distance < current_distance) && distance > mSnapTolerance)
            {
                current_distance = distance;
            }
        }
    }

    for (IndexType i = 0; i < r_edges_dest.size(); ++i)
    {
        const auto &v1 = V1byEdgeId(i, OrientationDest);
        const auto &v2 = V2byEdgeId(i, OrientationDest);

        double pos1 = rPoint[0];
        double pos2 = 0.0;

        if (pos1 >= v1[0] + mSnapTolerance && pos1 <= v2[0] - mSnapTolerance)
        {
            pos2 = v1[1] + (v2[1] - v1[1]) / (v2[0] - v1[0]) * (rPoint[0] - v1[0]);
            double distance = (!is_positive) ? (rPoint[1] - pos2) : (pos2 - rPoint[1]);

            if ((distance < current_distance-ZEROTOL) && distance > mSnapTolerance)
            {
                current_distance = distance;
                intersection_point[0] = pos1;
                intersection_point[1] = pos2;
                edge_id = i;
            }
        }
    }
    if (edge_id > -1)
    {
        r_edges_dest[edge_id].AddSplitPoint(intersection_point);
    }
}

void TrimmedDomainOnPlane::SplitEdgesAtSplitPoint(OrientationType Orientation)
{
    auto& r_edges = GetEdges(Orientation);

    IndexType pos = 0;
    IndexType size = r_edges.size();

    while( pos < size ) {
        auto& edge = r_edges[pos];
        auto& split_points = edge.GetSplitPoints();

        auto normal = edge.Normal();
        const IndexType num_split_points = split_points.size();
        if (num_split_points > 0) {
            std::sort(split_points.begin(), split_points.end(), [](auto &point_a, auto &point_b) -> bool
                { return point_a[0] < point_b[0]; });

            // Insert egde (vertex V1 + first split point)
            IndexType index1 = edge.V1();
            IndexType index2 = InsertVertex(split_points[0], Orientation);
            r_edges.push_back(Edge2D(index1, index2, normal) );

            // Insert egdes (only split points.)
            for (IndexType j = 0; j < num_split_points - 1; ++j) {
                index1 = InsertVertex(split_points[j], Orientation);
                index2 = InsertVertex(split_points[j + 1], Orientation);
                r_edges.push_back(Edge2D(index1, index2, normal));
            }

            // Insert egde (last split point + vertex V2)
            index1 = InsertVertex(split_points[num_split_points - 1], Orientation);
            index2 = edge.V2();
            r_edges.push_back(Edge2D(index1, index2, normal));

            split_points.clear();

            // Remove original edge
            r_edges.erase( r_edges.begin() + pos );
            --size;
        }
        else {
            // Keep original edge as is and increment to next one.
            ++pos;
        }
    }
}

////////////////////////
/// Setter Functions ///
////////////////////////

bool TrimmedDomainOnPlane::InsertEdge(const Point2DType& rV1, const Point2DType& rV2, const Point2DType& rNormal, OrientationType Orientation ){
    // Get unique vertex indices.
    auto indices = GetUniqueVertexIDs(rV1, rV2, Orientation);
    // Only insert if indices are not the same.
    if (indices.first != indices.second) {
        InsertVertex(rV1, indices.first, Orientation);
        InsertVertex(rV2, indices.second, Orientation);
        auto& r_edges = GetEdges(Orientation);
        r_edges.push_back(Edge2D(indices.first, indices.second, rNormal));
        return true;
    }
    return false;
}

void TrimmedDomainOnPlane::InsertEdge(const Point3DType& rV1, const Point3DType& rV2, const Point3DType &rNormal)
{
    // Insantiate 2D points.
    Point2DType v1{};
    Point2DType v2{};

    // Make sure edge is oriented along DIRINDEX1 such that x2 > x1.
    if (rV2[DIRINDEX1] > rV1[DIRINDEX1]) {
        v1 = {rV1[DIRINDEX1], rV1[DIRINDEX2]};
        v2 = {rV2[DIRINDEX1], rV2[DIRINDEX2]};
    }
    else {
        v1 = {rV2[DIRINDEX1], rV2[DIRINDEX2]};
        v2 = {rV1[DIRINDEX1], rV1[DIRINDEX2]};
    }

    Point2DType normal;
    normal[0] = rNormal[DIRINDEX1];
    normal[1] = rNormal[DIRINDEX2];
    if( std::abs( v1[0] - v2[0] ) > mSnapTolerance ){ //|| std::abs( v1[1] - v2[1] ) > mSnapTolerance ){
        if ((rNormal[DIRINDEX2] > ZEROTOL) ) { // Positive oriented
            InsertEdge(v1, v2, normal, Orientation::Positive);
        }
        else if ((rNormal[DIRINDEX2] < -ZEROTOL) ) { // Negative oriented
            InsertEdge(v1, v2, normal, Orientation::Negative);
        }
        else {
            // Vertical oriented
            if( std::abs( v1[1] - v2[1] ) > mSnapTolerance ){
                InsertEdge(v1, v2, normal, Orientation::Vertical);
            }
        }
    } else {
        // Vertical oriented
        if( std::abs( v1[1] - v2[1] ) > mSnapTolerance ){
            InsertEdge(v1, v2, normal, Orientation::Vertical);
        }
    }
}

IndexType TrimmedDomainOnPlane::InsertVertex(const Point2DType &rPoint, OrientationType Orientation) {
    auto &r_vertices = GetVertices(Orientation);
    IndexType index = r_vertices.size();
    r_vertices.push_back(rPoint);
    return index;
}

void TrimmedDomainOnPlane::InsertVertex(const Point2DType &rPoint, IndexType NewIndex, OrientationType Orientation) {
    auto &r_vertices = GetVertices(Orientation);
    auto& r_vertices_set = GetVerticesSet(Orientation);
    IndexType index = r_vertices.size();
    if (NewIndex == index) {
        r_vertices.push_back(rPoint);
        r_vertices_set.insert(--r_vertices.end());
        //QuESo_ERROR_IF("TrimmedDomainOnPlane::InsertVertex", !res.second) << "Vetrex already exists.\n";
    }
    else if(NewIndex > index) {
        QuESo_ERROR << "Given index out of range.\n";
    }
}

void TrimmedDomainOnPlane::AddIntersectedVertexPositive(Edge2D* pEdge, std::vector<std::pair<double, bool>>& rVertices){
    auto status = pEdge->IsVertexOnUpperBoundary();
    if( status.first ){
        auto normal = pEdge->Normal();
        auto v1 = mVerticesPositive[pEdge->V1()];
        if (normal[0] > -ZEROTOL) { // Just additional check. Actually not required.
            rVertices.push_back( std::make_pair(v1[0], false)  );
        }
    }
    else if( status.second ){
        auto normal = pEdge->Normal();
        auto v2 = mVerticesPositive[pEdge->V2()];
        if ( normal[0] < ZEROTOL ){ // Just additional check. Actually not required.
            rVertices.push_back( std::make_pair(v2[0], true)  );
        }
    }
}

void TrimmedDomainOnPlane::AddIntersectedVertexNegative(Edge2D* pEdge, std::vector<std::pair<double, bool>>& rVertices){
    auto status = pEdge->IsVertexOnUpperBoundary();
    if( status.first ){
        auto normal = pEdge->Normal();
        auto v1 = mVerticesNegative[pEdge->V1()];
        if( normal[0] < ZEROTOL ) { // Just additional check. Actually not required.
            rVertices.push_back( std::make_pair(v1[0], true)  );
        }
    }
    else if( status.second ){
        auto normal = pEdge->Normal();
        auto v2 = mVerticesNegative[pEdge->V2()];
        if(normal[0] > -ZEROTOL) { // Just additional check. Actually not required.
            rVertices.push_back( std::make_pair(v2[0], false)  );
        }
    }
}

void TrimmedDomainOnPlane::AddIntersectedVertexVertical(Edge2D* pEdge, std::vector<std::pair<double, bool>>& rVertices){
    auto status = pEdge->IsVertexOnUpperBoundary();
    if( status.first || status.second ){
        auto normal = pEdge->Normal();
        auto v1 = mVerticesVertical[pEdge->V1()];
        if( normal[0] < ZEROTOL ) {
            rVertices.push_back( std::make_pair(v1[0], true)  );
        } else {
            rVertices.push_back( std::make_pair(v1[0], false)  );
        }
    }
}

void TrimmedDomainOnPlane::RemoveDublicateVerticalEdges(
        std::vector<Edge2D>& rEdges, std::vector<Point2DType>& rVertices){

    const double tolerance = mSnapTolerance;
    // 1. Remove all dublicate verticals, that have equal normals.
    rEdges.erase(std::unique(rEdges.begin(), rEdges.end(), [&rVertices, tolerance](const auto &rValue1, const auto &rValue2){
        const bool normal_1_is_neg = rValue1.Normal()[0] < 0.0;
        const bool normal_2_is_neg = rValue2.Normal()[0] < 0.0;
        return normal_1_is_neg == normal_2_is_neg
            && std::abs(rVertices[rValue1.V1()][0] - (rVertices[rValue2.V1()][0])) < tolerance;
        }), rEdges.end() );

    // 2. Only keep normals that are unique (Only contained once!).
    std::vector<Edge2D> r_edges_copy( rEdges.begin(), rEdges.end() );
    rEdges.clear();
    for( auto r_edge : r_edges_copy){
        IndexType count = std::count_if(r_edges_copy.begin(), r_edges_copy.end(), [&r_edge, &rVertices, tolerance](const auto& rValue)
            { return std::abs(rVertices[r_edge.V1()][0] - (rVertices[rValue.V1()][0])) < tolerance; });

        if(count == 1 )
            rEdges.push_back(r_edge);
    }
}

void TrimmedDomainOnPlane::RemoveDublicateEdges(
        std::vector<Edge2D>& rEdges, std::vector<Point2DType>& rVertices){

    const double tolerance = mSnapTolerance;
    rEdges.erase(std::unique(rEdges.begin(), rEdges.end(), [&rVertices, tolerance](const auto &rValue1, const auto &rValue2){
        return std::abs(rVertices[rValue1.V1()][0] - (rVertices[rValue2.V1()][0])) < tolerance
            && std::abs(rVertices[rValue1.V2()][0] - (rVertices[rValue2.V2()][0])) < tolerance; }), rEdges.end() );
}

////////////////////////
/// Getter Functions ///
////////////////////////

std::pair<IndexType, IndexType> TrimmedDomainOnPlane::GetUniqueVertexIDs(const Point2DType &rV1, const Point2DType &rV2, OrientationType Orientation) const
{
    const auto& r_vertices = GetVertices(Orientation);
    const auto& r_vertices_set = GetVerticesSet(Orientation);

    // Instaniate a tmp_vector, as the following functio call require iterators.
    std::vector<Point2DType> tmp_vector = {rV1, rV2};
    // If points are the same.
    if( !PointComparison(mSnapTolerance)(tmp_vector.begin(), ++tmp_vector.begin())  &&
            !PointComparison(mSnapTolerance)(++tmp_vector.begin(), tmp_vector.begin()) ){
        return std::make_pair<IndexType, IndexType>(0UL, 0UL);
    }

    // r_vertices_set is only used to have a fast search here.
    const auto v1_res = r_vertices_set.find(tmp_vector.begin());
    const auto v2_res = r_vertices_set.find(++tmp_vector.begin());

    // If same points are found and both are not r_vertices_set.end().
    if( v1_res == v2_res && v1_res != r_vertices_set.end() ){
        return std::make_pair<IndexType, IndexType>(0UL, 0UL);
    }

    IndexType index_1 = 0UL;
    IndexType index_2 = 0UL;
    IndexType v1_is_new = 0UL;
    if (v1_res != r_vertices_set.end()) { // Vertex 1 already exists
        index_1 = std::distance<std::vector<Point2DType>::const_iterator>(r_vertices.cbegin(), (*v1_res));
    }
    else { // Add new vertex 1
        index_1 = r_vertices.size();
        v1_is_new = 1UL;
    }

    if (v2_res != r_vertices_set.end()) { // Vertex 2 already exists
        index_2 = std::distance<std::vector<Point2DType>::const_iterator>(r_vertices.cbegin(), (*v2_res));
    }
    else { // Add new vertex 2
        index_2 = r_vertices.size() + v1_is_new;
    }

    return std::make_pair(index_1, index_2);
}

const Point2DType& TrimmedDomainOnPlane::V1byEdgeId(IndexType EdgeId, OrientationType Orientation) const {
    switch (Orientation)
    {
    case Orientation::Positive:
        return mVerticesPositive[mEdgesPositiveOriented[EdgeId].V1()];
    case Orientation::Negative:
        return mVerticesNegative[mEdgesNegativeOriented[EdgeId].V1()];
    case Orientation::Vertical:
        return mVerticesVertical[mEdgesVertical[EdgeId].V1()];
    default:
        QuESo_ERROR << "Given Orientation not available.\n";
    }
}

const Point2DType& TrimmedDomainOnPlane::V2byEdgeId(IndexType EdgeId, OrientationType Orientation) const {
    switch (Orientation)
    {
    case Orientation::Positive:
        return mVerticesPositive[mEdgesPositiveOriented[EdgeId].V2()];
    case Orientation::Negative:
        return mVerticesNegative[mEdgesNegativeOriented[EdgeId].V2()];
    case Orientation::Vertical:
        return mVerticesVertical[mEdgesVertical[EdgeId].V2()];
    default:
        QuESo_ERROR << "Given Orientation not available.\n";
    }
}

const std::vector<Edge2D>& TrimmedDomainOnPlane::GetEdges(OrientationType Orientation) const {
    switch (Orientation)
    {
    case Orientation::Positive:
        return mEdgesPositiveOriented;
    case Orientation::Negative:
        return mEdgesNegativeOriented;
    case Orientation::Vertical:
        return mEdgesVertical;
    default:
        QuESo_ERROR << "Given Orientation not available.\n";
    }
}

std::vector<Edge2D>& TrimmedDomainOnPlane::GetEdges(OrientationType Orientation) {
    switch (Orientation)
    {
    case Orientation::Positive:
        return mEdgesPositiveOriented;
    case Orientation::Negative:
        return mEdgesNegativeOriented;
    case Orientation::Vertical:
        return mEdgesVertical;
    default:
        QuESo_ERROR << "Given Orientation not available.\n";
    }
}

const std::vector<Point2DType>& TrimmedDomainOnPlane::GetVertices(OrientationType Orientation) const
{
    switch (Orientation)
    {
    case Orientation::Positive:
        return mVerticesPositive;
    case Orientation::Negative:
        return mVerticesNegative;
    case Orientation::Vertical:
        return mVerticesVertical;
    default:
        QuESo_ERROR << "Given Orientation not available.\n";
    }
}

std::vector<Point2DType>& TrimmedDomainOnPlane::GetVertices(OrientationType Orientation)
{
    switch (Orientation)
    {
    case Orientation::Positive:
        return mVerticesPositive;
    case Orientation::Negative:
        return mVerticesNegative;
    case Orientation::Vertical:
        return mVerticesVertical;
    default:
        QuESo_ERROR << "Given Orientation not available.\n";
    }
}

Point2DSetType& TrimmedDomainOnPlane::GetVerticesSet(OrientationType Orientation) {
    switch (Orientation)
    {
    case Orientation::Positive:
        return *mVerticesSetPositive.get();
    case Orientation::Negative:
        return *mVerticesSetNegative.get();
    case Orientation::Vertical:
        return *mVerticesSetVertical.get();
    default:
        QuESo_ERROR << "Given Orientation not available.\n";
        break;
    }
}

const Point2DSetType& TrimmedDomainOnPlane::GetVerticesSet(OrientationType Orientation) const {
    switch (Orientation)
    {
    case Orientation::Positive:
        return *mVerticesSetPositive.get();
    case Orientation::Negative:
        return *mVerticesSetNegative.get();
    case Orientation::Vertical:
        return *mVerticesSetVertical.get();
    default:
        QuESo_ERROR << "Given Orientation not available.\n";
        break;
    }
}

IndexType TrimmedDomainOnPlane::GetNumberEdges(OrientationType Orientation) const {
    switch (Orientation)
    {
    case Orientation::Positive:
        return mEdgesPositiveOriented.size();
    case Orientation::Negative:
        return mEdgesNegativeOriented.size();
    case Orientation::Vertical:
        return mEdgesVertical.size();
    default:
        QuESo_ERROR << "Given Orientation not available.\n";
        break;
    }
}

double TrimmedDomainOnPlane::GetPlanePosition() const {
    return mUpperBoundary ? mUpperBound[DIRINDEX3] : mLowerBound[DIRINDEX3];
}

double TrimmedDomainOnPlane::GetMaxDIRINDEX2(std::vector<Point2DType>& rPoints) const {
    double max_value = LOWESTD;
    for( auto& point : rPoints ){
        if( point[1] > max_value ){
            max_value = point[1];
        }
    }
    return max_value;
}

bool TrimmedDomainOnPlane::PointExists(const Point2DType& rPoint , const Point2DSetType& rPointSet) const{
    std::vector<Point2DType> tmp_vector = {rPoint};
    auto it = rPointSet.find(tmp_vector.begin());
    return it != rPointSet.end();
}


} // End namespace queso