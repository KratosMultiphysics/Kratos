// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// Project includes
#include "embedding/polygon.h"
#include "utilities/math_utilities.hpp"

namespace queso {

// Function Definitions Polygon
template<IndexType SIZE>
IndexType Polygon<SIZE>::AddVertex(const PointType& rPoint ) {
    QuESo_ERROR_IF(mNumVertices >= SIZE) << "Size of Polygon is exceeded.\n";

    mVertices[mNumVertices].first = rPoint;
    mVertices[mNumVertices].second = {false};
    return mNumVertices++;
}

template<IndexType SIZE>
IndexType Polygon<SIZE>::AddVertex(const PointType& rPoint, const std::array<bool, 6>& rStatus) {
    QuESo_ERROR_IF(mNumVertices >= SIZE) << "Size of Polygon is exceeded.\n";

    mVertices[mNumVertices].first = rPoint;
    mVertices[mNumVertices].second = rStatus;
    return mNumVertices++;
}

template<IndexType SIZE>
IndexType Polygon<SIZE>::AddVertex(const VertexType& rPoint) {
    QuESo_ERROR_IF(mNumVertices >= SIZE) << "Size of Polygon is exceeded.\n";

    mVertices[mNumVertices] = rPoint;
    return mNumVertices++;
}

template<IndexType SIZE>
IndexType Polygon<SIZE>::NumVertices() const {
    return mNumVertices;
}

template<IndexType SIZE>
const typename Polygon<SIZE>::VertexType& Polygon<SIZE>::GetVertex(IndexType i) const {
    QuESo_ERROR_IF(i >= mNumVertices) << "Size of Polygon is exceeded.\n";

    return mVertices[i];
}

template<IndexType SIZE>
const typename Polygon<SIZE>::VertexType& Polygon<SIZE>::operator[] (IndexType i) const {
    QuESo_ERROR_IF(i >= mNumVertices) << "Size of Polygon is exceeded.\n";

    return mVertices[i];
}

template<IndexType SIZE>
const typename Polygon<SIZE>::VertexType& Polygon<SIZE>::GetLastVertex() const{
    return mVertices[mNumVertices-1];
}

template<IndexType SIZE>
void Polygon<SIZE>::Clear(){
    const std::array<bool,6> points_on_plane = {false}; // Inits all to false;
    std::fill(mVertices.begin(), mVertices.end(), std::make_pair(PointType{}, points_on_plane) );
    mNumVertices = 0;
}

template<IndexType SIZE>
void Polygon<SIZE>::AddToTriangleMesh(TriangleMeshInterface& rTriangleMesh) const {
    if(mNumVertices < 3){
        return;
    }

    const IndexType num_v = rTriangleMesh.NumOfVertices();
    const IndexType num_t = rTriangleMesh.NumOfTriangles();

    if( mNumVertices == 3 ){
        rTriangleMesh.AddVertex( mVertices[0].first );
        rTriangleMesh.AddVertex( mVertices[1].first );
        rTriangleMesh.AddVertex( mVertices[2].first );

        rTriangleMesh.AddTriangle( {num_v+0, num_v+1, num_v+2} );
        rTriangleMesh.AddNormal( mNormal );

        // Add edges, that are located on a plane, to the mesh.
        // Planes: (-x, +x, -y, y, -z, z)
        for( IndexType plane_index = 0; plane_index < 6; ++plane_index ){
            const bool v1_on_plane = mVertices[0].second[plane_index];
            const bool v2_on_plane = mVertices[1].second[plane_index];
            const bool v3_on_plane = mVertices[2].second[plane_index];
            if( (v1_on_plane + v2_on_plane + v3_on_plane) == 3 ) {
                QuESo_ERROR << "All vertices are set on plane.\n";
            }
            if( v1_on_plane && v2_on_plane ){
                rTriangleMesh.AddEdgeOnPlane(plane_index, num_v+0, num_v+1, num_t+0);
            }
            else if( v2_on_plane && v3_on_plane ){
                rTriangleMesh.AddEdgeOnPlane(plane_index, num_v+1, num_v+2, num_t+0);
            }
            else if( v3_on_plane && v1_on_plane ){
                rTriangleMesh.AddEdgeOnPlane(plane_index, num_v+2, num_v+0, num_t+0);
            }
        }
        return;
    }

    // Compute mean of vertices
    PointType centroid = {0.0, 0.0, 0.0};
    for( IndexType i = 0 ; i < mNumVertices; ++i){
        Math::AddSelf(centroid, mVertices[i].first);
    }
    Math::DivideSelf(centroid, static_cast<double>(mNumVertices));

    IndexType vertex_count = num_v;
    IndexType triangle_count = num_t;
    for( IndexType i = 0 ; i < mNumVertices-1; ++i){
        rTriangleMesh.AddVertex( mVertices[i].first );
        rTriangleMesh.AddVertex( mVertices[i+1].first );
        rTriangleMesh.AddVertex( centroid );

        rTriangleMesh.AddTriangle( {vertex_count+0, vertex_count+1, vertex_count+2} );
        rTriangleMesh.AddNormal( mNormal );

        for( IndexType plane_index = 0; plane_index < 6; ++plane_index ){
            const bool v1_on_plane = mVertices[i].second[plane_index];
            const bool v2_on_plane = mVertices[i+1].second[plane_index];
            if( v1_on_plane && v2_on_plane ){
                rTriangleMesh.AddEdgeOnPlane(plane_index, vertex_count+0, vertex_count+1, triangle_count);
            }
        }
        ++triangle_count;
        vertex_count += 3;
    }

    rTriangleMesh.AddVertex( mVertices[mNumVertices-1].first );
    rTriangleMesh.AddVertex( mVertices[0].first );
    rTriangleMesh.AddVertex( centroid );

    rTriangleMesh.AddTriangle( {vertex_count+0, vertex_count+1, vertex_count+2 } );
    rTriangleMesh.AddNormal( mNormal );

    for( IndexType plane_index = 0; plane_index < 6; ++plane_index ){
        const bool v1_on_plane = mVertices[mNumVertices-1].second[plane_index];
        const bool v2_on_plane = mVertices[0].second[plane_index];
        if( v1_on_plane && v2_on_plane ){
            rTriangleMesh.AddEdgeOnPlane(plane_index, vertex_count+0, vertex_count+1, triangle_count);
        }
    }

    return;
}

// Explicit instantiation Polygon
template class Polygon<9>;
template class Polygon<4>;

} // End namespace queso