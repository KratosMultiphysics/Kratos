// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <map>
#include <numeric>
#include <algorithm>

/// Project includes
#include "utilities/mesh_utilities.h"
#include "utilities/math_utilities.hpp"

namespace queso {

typedef MeshUtilities::TriangleMeshPtrType TriangleMeshPtrType;

void MeshUtilities::Refine(TriangleMeshInterface& rTriangleMesh, IndexType MinNumberOfTriangles){
    IndexType original_size = rTriangleMesh.NumOfTriangles();
    IndexType size = original_size;

    double max_area = -1.0;
    for( IndexType pos = 0; pos < size; ++pos){
        max_area = std::max<double>( max_area, rTriangleMesh.Area(pos));
    }

    const auto& r_vertices = rTriangleMesh.GetVertices();
    rTriangleMesh.Reserve(4*size);
    IndexType pos = 0;
    while( pos < size ){
        // Make sure this is a copy!!
        auto vertex_ids = rTriangleMesh.VertexIds(pos);
        const Vector3d p1 = r_vertices[vertex_ids[0]];
        const Vector3d p2 = r_vertices[vertex_ids[1]];
        const Vector3d p3 = r_vertices[vertex_ids[2]];

        const double area = rTriangleMesh.Area(pos);
        if( area > 0.5*max_area ){
            IndexType e1 = rTriangleMesh.AddVertex( Math::AddAndMult(0.5, p1, p2) );
            IndexType e2 = rTriangleMesh.AddVertex( Math::AddAndMult(0.5, p2, p3) );
            IndexType e3 = rTriangleMesh.AddVertex( Math::AddAndMult(0.5, p3, p1) );

            const auto normal = rTriangleMesh.Normal(pos);
            rTriangleMesh.AddTriangle( {vertex_ids[0], e1, e3} );
            rTriangleMesh.AddTriangle( {e1, vertex_ids[1], e2} );
            rTriangleMesh.AddTriangle( {e2, vertex_ids[2], e3} );
            rTriangleMesh.AddTriangle( {e1, e2, e3} );
            rTriangleMesh.AddNormal( normal );
            rTriangleMesh.AddNormal( normal );
            rTriangleMesh.AddNormal( normal );
            rTriangleMesh.AddNormal( normal );
            size += 4;

            rTriangleMesh.RemoveTriangle(pos);
            rTriangleMesh.RemoveNormal(pos);
            --size;
        } else {
            ++pos;
        }
        if( pos > original_size ){
            rTriangleMesh.Reserve(4*original_size);
            original_size = size;
        }
    }
    if( rTriangleMesh.NumOfTriangles() < MinNumberOfTriangles ){
        Refine(rTriangleMesh, MinNumberOfTriangles );
    }
}

void MeshUtilities::Append(TriangleMeshInterface& rTriangleMesh, const TriangleMeshInterface& rNewMesh){
    std::vector<IndexType> indices(rNewMesh.NumOfTriangles());
    /// Fill vector with number of increasing values: 0,1,2..
    std::iota( indices.begin(), indices.end(), 0 );
    Append(rTriangleMesh, rNewMesh, indices);
}

void MeshUtilities::Append(TriangleMeshInterface& rTriangleMesh, const TriangleMeshInterface& rNewMesh, const std::vector<IndexType>& rIndices){

    const IndexType initial_number_triangles = rTriangleMesh.NumOfTriangles();
    const IndexType initial_number_vertices = rTriangleMesh.NumOfVertices();
    IndexType vertex_count = initial_number_vertices;
    std::map<IndexType, IndexType> index_map_vertices{};

    for( auto triangle : rIndices){
        const auto& tmp_indices = rNewMesh.VertexIds(triangle);
        Vector3i new_triangle{};
        IndexType ii = 0;
        for( auto index : tmp_indices ){
            // Insert index into index_map_vertices if map does not contain index.
            auto ret = index_map_vertices.insert( std::pair<IndexType,IndexType>(index, vertex_count) );
            if (ret.second==true) {
                new_triangle[ii] = vertex_count;
                vertex_count++;
            } else {
                new_triangle[ii] = index_map_vertices[index];
            }
            ++ii;
        }
        // Copy triangles and normals.
        rTriangleMesh.AddTriangle(new_triangle);
        rTriangleMesh.AddNormal( rNewMesh.Normal(triangle) );
    }

    auto& r_vertices = rTriangleMesh.GetVertices();
    r_vertices.resize(vertex_count);
    const auto& r_new_vertices = rNewMesh.GetVertices();

    // Copy vertices.
    for( auto index : index_map_vertices ){
        r_vertices[ index.second ] = r_new_vertices[ index.first ];
    }

    // Copy edges.
    const auto& new_edges_on_plane = rNewMesh.GetEdgesOnPlanes();
    for( IndexType plane_index = 0; plane_index < 6; ++plane_index){
        const auto& edges = new_edges_on_plane[plane_index];
        for( auto& edge : edges ){
            rTriangleMesh.AddEdgeOnPlane(plane_index, std::get<0>(edge)+initial_number_vertices,
                                                      std::get<1>(edge)+initial_number_vertices,
                                                      std::get<2>(edge)+initial_number_triangles );
        }
    }
}

Unique<TriangleMeshInterface> MeshUtilities::pGetCuboid(const PointType& rLowerPoint, const PointType& rUpperPoint){
    //
    //     2_______3                 y
    //     /      /|                ´|`
    //   6/_____7/ |                 |-->x
    //    | 0   |  /1               /
    //    |     | /                Z
    //   4|____5|/
    //
    auto p_new_triangle_mesh = MakeUnique<TriangleMesh>();

    p_new_triangle_mesh->AddVertex( {rLowerPoint[0], rLowerPoint[1], rLowerPoint[2]} ); //0
    p_new_triangle_mesh->AddVertex( {rUpperPoint[0], rLowerPoint[1], rLowerPoint[2]} ); //1
    p_new_triangle_mesh->AddVertex( {rLowerPoint[0], rUpperPoint[1], rLowerPoint[2]} ); //2
    p_new_triangle_mesh->AddVertex( {rUpperPoint[0], rUpperPoint[1], rLowerPoint[2]} ); //3
    p_new_triangle_mesh->AddVertex( {rLowerPoint[0], rLowerPoint[1], rUpperPoint[2]} ); //4
    p_new_triangle_mesh->AddVertex( {rUpperPoint[0], rLowerPoint[1], rUpperPoint[2]} ); //5
    p_new_triangle_mesh->AddVertex( {rLowerPoint[0], rUpperPoint[1], rUpperPoint[2]} ); //6
    p_new_triangle_mesh->AddVertex( {rUpperPoint[0], rUpperPoint[1], rUpperPoint[2]} ); //7

    // negative x
    p_new_triangle_mesh->AddTriangle({0, 6, 2});
    p_new_triangle_mesh->AddNormal({-1.0, 0.0, 0.0});
    p_new_triangle_mesh->AddTriangle({0, 4, 6});
    p_new_triangle_mesh->AddNormal({-1.0, 0.0, 0.0});

    // postive x
    p_new_triangle_mesh->AddTriangle({1, 7, 5});
    p_new_triangle_mesh->AddNormal({1.0, 0.0, 0.0});
    p_new_triangle_mesh->AddTriangle({1, 3, 7});
    p_new_triangle_mesh->AddNormal({1.0, 0.0, 0.0});

    // negative y
    p_new_triangle_mesh->AddTriangle({4, 1, 5});
    p_new_triangle_mesh->AddNormal({0.0, -1.0, 0.0});
    p_new_triangle_mesh->AddTriangle({4, 0, 1});
    p_new_triangle_mesh->AddNormal({0.0, -1.0, 0.0});

    // postive y
    p_new_triangle_mesh->AddTriangle({6, 7, 3});
    p_new_triangle_mesh->AddNormal({0.0, 1.0, 0.0});
    p_new_triangle_mesh->AddTriangle({6, 3, 2});
    p_new_triangle_mesh->AddNormal({0.0, 1.0, 0.0});

    // negative z
    p_new_triangle_mesh->AddTriangle({1, 0, 3});
    p_new_triangle_mesh->AddNormal({0.0, 0.0, -1.0});
    p_new_triangle_mesh->AddTriangle({0, 2, 3});
    p_new_triangle_mesh->AddNormal({0.0, 0.0, -1.0});

    // positive z
    p_new_triangle_mesh->AddTriangle({4, 5, 7});
    p_new_triangle_mesh->AddNormal({0.0, 0.0, 1.0});
    p_new_triangle_mesh->AddTriangle({4, 7, 6});
    p_new_triangle_mesh->AddNormal({0.0, 0.0, 1.0});

    return p_new_triangle_mesh;
}


double MeshUtilities::Area(const TriangleMeshInterface& rTriangleMesh){
    double area = 0.0;
    // Loop over all triangles
    for( IndexType i = 0; i < rTriangleMesh.NumOfTriangles(); ++i ){
        area += rTriangleMesh.Area(i);
    }
    return area;
}

double MeshUtilities::Volume(const TriangleMeshInterface& rTriangleMesh){
    double volume = 0.0;
    const IndexType num_triangles = rTriangleMesh.NumOfTriangles();
    // Loop over all triangles
    for( IndexType i = 0; i < num_triangles; ++i ){
        const auto p_points = rTriangleMesh.pGetIPsGlobal<BoundaryIntegrationPoint>(i, 0);
        const auto& r_points = *p_points;
        // Loop over all points.
        for( const auto& point : r_points ){
            const auto& normal = point.Normal();
            const double integrand = Math::Dot(normal, point.data() );
            const double integral = integrand * point.Weight();
            if( std::abs(integral) > 0.0 ) { // This skips possible NaN-values.
                volume += integral;
            }
        }
    }
    return std::abs(1.0/3.0*volume);
}

double MeshUtilities::VolumeOMP(const TriangleMeshInterface& rTriangleMesh){
    double volume = 0.0;
    const IndexType num_triangles = rTriangleMesh.NumOfTriangles();
    // Loop over all triangles in omp parallel.
    #pragma omp parallel for reduction(+ : volume)
    for( int i = 0; i < static_cast<int>(num_triangles); ++i ){
        const auto p_points = rTriangleMesh.pGetIPsGlobal<BoundaryIntegrationPoint>(i, 0);
        const auto& r_points = *p_points;
        // Loop over all points.
        for( const auto& point : r_points ){
            const auto& normal = point.Normal();
            const double integrand = Math::Dot(normal, point.data() );
            const double integral = integrand * point.Weight();
            if( std::abs(integral) > 0.0 ) { // This skips possible NaN-values.
                volume += integral;
            }
        }
    }
    return std::abs(1.0/3.0*volume);
}

double MeshUtilities::Volume(const TriangleMeshInterface& rTriangleMesh, IndexType Dir){
    double volume = 0.0;
    const IndexType num_triangles = rTriangleMesh.NumOfTriangles();

    QuESo_ERROR_IF(Dir < 0 || Dir > 2 ) << " Directional Index is out-of-range.\n";

    // Loop over all triangles
    for( IndexType i = 0; i < num_triangles; ++i ){
        const auto p_points = rTriangleMesh.pGetIPsGlobal<BoundaryIntegrationPoint>(i, 0);
        const auto& r_points = *p_points;
        // Loop over all points.
        for( const auto& point : r_points ){
            const auto& normal = point.Normal();
            const double integrand = normal[Dir]*point[Dir];
            const double integral = integrand * point.Weight();
            if( std::abs(integral) > 0.0 ) { // This skips possible NaN-values.
                volume += integral;
            }
        }
    }
    return std::abs(volume);
}

double MeshUtilities::MaxAspectRatio(const TriangleMeshInterface& rTriangleMesh){
    double max_aspect_ratio = MIND;
    for( IndexType i = 0; i < rTriangleMesh.NumOfTriangles(); ++i){
        const double aspect_ratio = rTriangleMesh.AspectRatio(i);
        if( aspect_ratio > max_aspect_ratio ){
            max_aspect_ratio = aspect_ratio;
        }
    }
    return max_aspect_ratio;
}

double MeshUtilities::AverageAspectRatio(const TriangleMeshInterface& rTriangleMesh){
    double average_aspect_ratio = 0.0;
    for( IndexType i = 0; i < rTriangleMesh.NumOfTriangles(); ++i){
        average_aspect_ratio += rTriangleMesh.AspectRatio(i);
    }
    return average_aspect_ratio/rTriangleMesh.NumOfTriangles();
}

double MeshUtilities::EstimateQuality(const TriangleMeshInterface& rTriangleMesh ){
    const IndexType num_triangles = rTriangleMesh.NumOfTriangles();
    double total_volume_1 = 0.0;
    double total_volume_2 = 0.0;
    double total_volume_3 = 0.0;
    double total_area = 0.0;
    PointType directional_areas = {0.0, 0.0, 0.0};
    for( IndexType i = 0; i < num_triangles; ++i ){
        const double area = rTriangleMesh.Area(i);
        const auto normal = rTriangleMesh.Normal(i);
        total_area += area;
        Math::AddSelf(directional_areas, Math::Mult(area, normal) );

        // Get integration points
        const auto p_points = rTriangleMesh.pGetIPsGlobal<BoundaryIntegrationPoint>(i, 0);
        const auto& r_points = *p_points;
        // Loop over all points.
        for( const auto& point : r_points ){
            total_volume_1 += normal[0]*point[0] * point.Weight();
            total_volume_2 += normal[1]*point[1] * point.Weight();
            total_volume_3 += normal[2]*point[2] * point.Weight();
        }
    }
    const double total_volume = std::abs(1.0/3.0*(total_volume_1 + total_volume_2 + total_volume_3));
    const double error_v1 = std::abs(total_volume_1 - total_volume) / total_volume;
    const double error_v2 = std::abs(total_volume_2 - total_volume) / total_volume;
    const double error_v3 = std::abs(total_volume_3 - total_volume) / total_volume;

    const double error_area = Math::Norm(directional_areas) / std::abs(total_area);

    const double max_error = std::max(std::max(std::max( error_v1, error_v2), error_v3), error_area );
    return max_error;
}


std::pair<PointType, PointType> MeshUtilities::BoundingBox(const TriangleMeshInterface& rTriangleMesh) {
    PointType lower_bound = {MAXD, MAXD, MAXD};
    PointType upper_bound = {LOWESTD, LOWESTD, LOWESTD};
    for( IndexType i = 0; i < rTriangleMesh.NumOfTriangles(); ++i){
        const auto& p1 = rTriangleMesh.P1(i);
        const auto& p2 = rTriangleMesh.P2(i);
        const auto& p3 = rTriangleMesh.P3(i);

        const PointType x_values = {p1[0], p2[0], p3[0]};
        const PointType y_values = {p1[1], p2[1], p3[1]};
        const PointType z_values = {p1[2], p2[2], p3[2]};

        auto x_min_max = std::minmax_element(x_values.begin(), x_values.end());
        auto y_min_max = std::minmax_element(y_values.begin(), y_values.end());
        auto z_min_max = std::minmax_element(z_values.begin(), z_values.end());

        lower_bound[0] = std::min<double>(*x_min_max.first, lower_bound[0]);
        upper_bound[0] = std::max<double>(*x_min_max.second, upper_bound[0]);

        lower_bound[1] = std::min<double>(*y_min_max.first, lower_bound[1]);
        upper_bound[1] = std::max<double>(*y_min_max.second, upper_bound[1]);

        lower_bound[2] = std::min<double>(*z_min_max.first, lower_bound[2]);
        upper_bound[2] = std::max<double>(*z_min_max.second, upper_bound[2]);
    }

    return std::make_pair(lower_bound, upper_bound);
}

} // End namespace queso