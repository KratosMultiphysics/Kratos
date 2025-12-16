//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortes
//
//

#if !defined(KRATOS_VOXEL_UTILITIES)
#define  KRATOS_VOXEL_UTILITIES

// System includes

// External includes

// Project includes
#include "includes/geometrical_object.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "geometries/hexahedra_3d_8.h"
#include "intersection_utilities.h"

namespace Kratos
{ 
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class VoxelUtilities
 * @ingroup KratosCore
 * @brief Utilities to compute the real volume inside a voxel
 * @details This class provides static methods to compute (using different approximations) the portion of a 
 * voxel that is actually filled with volume, according to the known triangle elements that intersect the voxel
 * @author Ariadna Cortes
 */
class VoxelUtilities
{
public:

    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef Node<3>::Pointer NodePtrType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::GeometriesArrayType GeometryArrayType;
    typedef GeometryType::PointsArrayType PointsArrayType;

    /// Pointer definition of VoxelInsideVolume
    KRATOS_CLASS_POINTER_DEFINITION( VoxelUtilities );

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     */
    VoxelUtilities(){}

    /// Destructor
    virtual ~VoxelUtilities(){}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Aproximates the actual volume inside the voxel 
     * @param rVoxel references to the voxel whose actual volume will be approximated
     * @return Approximated volume 
     * @note This approximation assigns a fraction of volume (1/8) to each node of the
     * voxel, and counts it as volume if the node is inside the object (NodeDistance > 0)
     * This operation is VERY cheap
     */  
    static double NodesApproximation(const GeometryType& rVoxel);
    
    /**
     * @brief Aproximates the actual volume inside the voxel 
     * @param rVoxel references to the voxel whose actual volume will be approximated
     * @param rTriangles references to the triangles which intersect the voxel at some edge.
     * @return Approximated volume 
     * @note This approximation finds the portion of each edge that is part of the volume (using
     * intersection points with triangles of the mesh).
     */  
    static double EdgesPortionApproximation(const GeometryType& rVoxel, const GeometryArrayType& rTriangles);

    /**
     * @brief Aproximates the actual area inside a quadrilateral 
     * @param rFace references to the quadrilateral3D4 whose actual area will be approximated
     * @param rTriangles references to the triangles which may intersect the quadrilateral at some edge.
     * @return Approximated area 
     */  
    static double FaceArea(const GeometryType& rFace,const GeometryArrayType& rTriangles);

    /**
     * @brief Aproximates the portion of the edge (Line3D2) that represents volume
     * @param rDistances references to a sorted vector containing the relatve distances of each intersecting point 
     * with the first node of the edge 
     * @param rEnds references to the nodes at both sides of the edge
     * @return Approximated Portion of edge that contains volume 
     */  
    static double EdgeFilledPortion(std::vector<double>& Distances, const PointsArrayType& rEnds);

private:

    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    static constexpr std::size_t mNeighbours[4][2] = {{3,1},{0,2},{1,3},{2,0}}; /// Neighbour list

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Returns the area enclosed by a set of 4 points
     * @param rPoints the array of points
     * @return Area inside this points
     * @note This code is an approximation and won't work perfectly for "parabolic" quadrilaterals
     */  
    static double PointsArea(const PointsArrayType& rPoints);

    static double GetFactor(const PointsArrayType& rNodes, const int NodeIndex);

    static int GetCase(const PointsArrayType& rNodes, const int NodeIndex);

}; /* Class VoxelInsideVolumeUtility */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_VOXEL_UTILITIES  defined */