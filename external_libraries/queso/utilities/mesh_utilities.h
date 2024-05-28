// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MESH_UTILITIES_INCLUDE_H
#define MESH_UTILITIES_INCLUDE_H

//// Project includes
#include "containers/triangle_mesh.hpp"

namespace queso {

///@name QuESo Classes
///@{

///
/**
 * @class  MeshUtilities
 * @author Manuel Messmer
 * @brief  Provides operators for the TriangleMesh.
*/
class MeshUtilities {
public:
    ///@name Type Definitions
    ///@{
    typedef Unique<TriangleMeshInterface> TriangleMeshPtrType;

    ///@}
    ///@name Public Operations
    ///@{

    /// @brief Refines triangle mesh. Always conducts one refinement loop, such that area < 0.5*area_max.
    /// @param rTriangleMesh Triangle mesh to refine.
    /// @param MinNumberOfTriangles Minimum number of triangles in final mesh.
    /// @todo May needs improvement, more parameters etc.
    static void Refine(TriangleMeshInterface& rTriangleMesh, IndexType MinNumberOfTriangles);

    /// @brief Appends rTriangleMesh by rNewMesh.
    /// @param rTriangleMesh Mesh to be appended.
    /// @param rNewMesh New mesh to be inserted in rTriangleMesh.
    static void Append(TriangleMeshInterface& rTriangleMesh, const TriangleMeshInterface& rNewMesh);

    /// @brief Append rTriangleMesh with some triangles in rNewMesh (given by Indices).
    /// @param rTriangleMesh Mesh to be appended.
    /// @param rNewMesh New mesh to be inserted in rTriangleMesh.
    /// @param rIndices Indices of triangles to be copied.
    static void Append(TriangleMeshInterface& rTriangleMesh, const TriangleMeshInterface& rNewMesh, const std::vector<IndexType>& rIndices);

    ///@brief Return meshed cuboid.
    ///@param rLowerPoint
    ///@param rUpperPoint
    ///@return Unique<TriangleMeshInterface>
    static TriangleMeshPtrType pGetCuboid(const PointType& rLowerPoint, const PointType& rUpperPoint);

    /// @brief Returns surface area of triangle mesh.
    /// @param rTriangleMesh
    /// @return double.
    static double Area(const TriangleMeshInterface& rTriangleMesh);

    ///@brief Returns enclosed volume by triangle mesh. Uses divergence theorem to compute volume.
    ///@param rTriangleMesh
    ///@return double
    static double Volume(const TriangleMeshInterface& rTriangleMesh);

    /// @brief Returns enclosed volume by triangle mesh. Uses divergence theorem to compute volume.
    /// Only applies divergence theorem in direction Dir: 0-x, 1-y, 2-z.
    /// @param rTriangleMesh
    /// @param Dir
    /// @return double
    static double Volume(const TriangleMeshInterface& rTriangleMesh, IndexType Dir);

    ///@brief Returns enclosed volume by triangle mesh (OMP-version). Uses divergence theorem to compute volume.
    ///@param rTriangleMesh
    ///@return double
    static double VolumeOMP(const TriangleMeshInterface& rTriangleMesh);

    ///@brief Returns a quaility measure of triangle mesh. This function computes the volume of the mesh using
    ///       the divergence theorem. It applies the divergence theorem individually in each space direction (volume_1, volume_2, volume_3).
    ///       Theoretically, all directional volumes must be equal. Additionally, this function computes the directional
    ///       areas: Sum of (normal * area). This value must be zero if the mesh is closed.
    ///       The maximum relative errors of volume_1, volume_2, volume_3 and directional_area is returned.
    ///@param rTriangleMesh
    ///@return double maximum error.
    static double EstimateQuality(const TriangleMeshInterface& rTriangleMesh);

    ///@brief Returns maximum aspect ratio of all triangles in triangle mesh.
    ///@param rTriangleMesh
    ///@return double
    static double MaxAspectRatio(const TriangleMeshInterface& rTriangleMesh);

    ///@brief Returns average aspect ratio of all triangles in triangle mesh.
    ///@param rTriangleMesh
    ///@return double
    static double AverageAspectRatio(const TriangleMeshInterface& rTriangleMesh);


    ///@brief Returns axis-aligned bounding box of triangle mesh.
    ///@param rTriangleMesh
    ///@return std::pair<PointType, PointType>
    static std::pair<PointType, PointType> BoundingBox(const TriangleMeshInterface& rTriangleMesh);

    ///@}
}; // End class MeshUtilities
///@} // End QuESo Classes
} // End namespace queso

#endif // MESH_UTILITIES_INCLUDE_H