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

#ifndef TRIANGLE_MESH_INCLUDE_H
#define TRIANGLE_MESH_INCLUDE_H

//// STL includes
#include <vector>
#include <array>
#include <cmath>
#include <memory>
#include <iostream>

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/containers/triangle_gauss_legendre_integration_points.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/triangle_mesh_interface.hpp"

namespace queso {

///@name QuESo Classes
///@{
/**
 * @class  TriangleMesh
 * @author Manuel Messmer
 * @brief  Simple implementation of a triangular surface mesh. Derives from TriangleMeshInterface.
 * @todo Refactor AddTriangle() and RemoveTriangle() such that normals are always passed.
 *       Put dot product etc. into math function.
*/
class TriangleMesh : public TriangleMeshInterface
{
public:
    ///@name Type Definitions
    ///@{

    typedef TriangleMeshInterface BaseType;
    typedef std::vector<std::vector<std::tuple<IndexType, IndexType, IndexType>>> EdgesOnPlanesVectorType;

    ///@}
    ///@name Life cycle
    ///@{

    TriangleMesh() = default;
    ~TriangleMesh() = default;
    TriangleMesh(const TriangleMesh& rVal) = default;
    TriangleMesh& operator = (const TriangleMesh& rOther) = default;

    ///@}
    ///@name Pure virtual Operations
    ///@{

    ///@brief Get number of triangles in mesh.
    IndexType NumOfTriangles() const override {
        return mTriangles.size();
    }

    ///@brief Get number of vertices in mesh.
    IndexType NumOfVertices() const override{
        return mVertices.size();
    }

    ///@brief Get triangle vertex 1
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3d& P1(IndexType TriangleId) const override {
        return mVertices[mTriangles[TriangleId][0]];
    }

    ///@brief Get triangle vertex 2
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3d& P2(IndexType TriangleId) const override {
        return mVertices[mTriangles[TriangleId][1]];
    }

    ///@brief Get triangle vertex 3
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3d& P3(IndexType TriangleId) const override {
        return mVertices[mTriangles[TriangleId][2]];
    }

    ///@}
    ///@name Virtual Operations
    ///@{

    /// @brief Clone object
    Unique<TriangleMeshInterface> Clone() override{
        return MakeUnique<TriangleMesh>(*this);
    }

    ///@brief Clear all containers.
    void Clear() override {
        BaseType::Clear();
        mTriangles.clear();
        mVertices.clear();
    }

    ///@brief Reserve all containers for normal vertices and triangles.
    ///       Note, call ReserveEdgesOnPlane() to reserve edges containers.
    ///@param Size
    void Reserve(IndexType Size) override {
        BaseType::Reserve(Size);
        mVertices.reserve(Size);
        mTriangles.reserve(Size);
    }

    ///@brief Add vertex to mesh.
    ///@param NewVertex
    IndexType AddVertex(const Vector3d& NewVertex) override {
        mVertices.push_back(NewVertex);
        return mVertices.size()-1;
    }

    ///@brief Add triangle to mesh.
    ///@param NewTriangle
    void AddTriangle(const Vector3i& NewTriangle) override {
        mTriangles.push_back(NewTriangle);
    }

    /// @brief Remove triangle by index.
    /// @param Index
    void RemoveTriangle(IndexType Index ) override {
        mTriangles.erase( mTriangles.begin() + Index );
    }

    ///@brief Get vertices from mesh. (const version)
    ///@return const std::vector<Vector3d>&
    const std::vector<Vector3d>& GetVertices() const override {
        return mVertices;
    }

    /// @brief Get vertices from mesh. (non-const version)
    /// @return std::vector<Vector3d>&
    std::vector<Vector3d>& GetVertices() override {
        return mVertices;
    }

    ///@brief Get triangles from mesh. (const version)
    ///@return const std::vector<Vector3i>&
    const std::vector<Vector3i>& GetTriangles() const override {
        return mTriangles;
    }

    ///@brief Get triangle vertex 3
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3i& VertexIds(IndexType TriangleId) const override {
        return mTriangles[TriangleId];
    }

    ///@brief Basic check of this TriangleMesh instance.
    bool Check() const override {
        // Check if mTriangles and mNormals are of the same size.
        if( mTriangles.size() != BaseType::NumOfNormals() ){
            std::cerr << "TriangleMesh :: Number of Triangles and Normals in mesh do not match.\n";
            return false;
        }
        // Check if all vertex ids exist.
        for( IndexType i = 0; i < mTriangles.size(); ++i ){
            for(IndexType j = 0; j < 3; ++j){
                if( mTriangles[i][j] >= mVertices.size() ){
                    std::cerr << "TriangleMesh :: Triangle/Vertex mismatch.\n";
                    return false;
                }
            }
        }
        return true;
    }
    ///@}

private:

    ///@}
    ///@name Private Member Variables
    ///@{
    std::vector<Vector3d> mVertices;
    std::vector<Vector3i> mTriangles;
    ///@}

}; // End of class TriangleMesh
///@} // End QuESo classes

} // End namespace queso

#endif // TRIANGLE_MESH_INCLUDE_H