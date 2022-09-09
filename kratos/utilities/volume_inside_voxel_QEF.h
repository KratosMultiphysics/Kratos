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

#pragma once

// System includes

// External includes

// Project includes
#include "qef_utility.h"
#include "voxel_utilities.h"

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
 * @class VolumeInsideVoxelQEF
 * @ingroup KratosCore
 * @brief Utilities to compute the real volume inside a voxel using quadratic error functions
 * @details This class provides static methods to compute (using different approximations) the portion of a 
 * voxel that is actually filled with volume.
 * @author Ariadna Cortes
 */
class VolumeInsideVoxelQEF : public VoxelUtilities
{
public:

    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::GeometriesArrayType GeometryArrayType;
    typedef GeometryType::PointsArrayType PointsArrayType;


    /// Pointer definition of VoxelInsideVolumeQEF
    KRATOS_CLASS_POINTER_DEFINITION( VolumeInsideVoxelQEF );

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     */
    VolumeInsideVoxelQEF(){}

    /// Destructor
    virtual ~VolumeInsideVoxelQEF(){}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Aproximates the actual volume inside the voxel 
     * @param rVoxel references to the voxel whose actual volume will be approximated
     * @param rTriangles references to the triangles which intersect the voxel at some edge.
     * @return Approximated volume 
     * @note This is cheaper version of VolumeQEFApproximation. 
     */ 
    static double SimpleNodesQEFApproximation(
        const GeometryType& rVoxel, 
        const GeometryArrayType& rTriangles);

    /**
     * @brief Aproximates the actual volume inside the voxel based on the triangle intersections at the edges  
     * @param rVoxel references to the voxel whose actual volume will be approximated
     * @param rTriangles references to the triangles which intersect the voxel at some edge.
     * @return Approximated volume 
     */
    static double VolumeQEFApproximation(
        const GeometryType& rVoxel,
        const GeometryArrayType& rTriangles);

private:

    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * @brief Approximates the portion of a face that actually corresponds to area (assigning each node 
     * 1/numberOfNodes portion if it is inside the volume)
     * @param Nodes The nodes of the geometry
     * @param NodesInFaces matrix containing the index of the nodes of the geometry that belong to each faces
     * @return 
     * */
    static double GetPortion(
        PointsArrayType& rNodes,
        const DenseMatrix<unsigned int>& rNodesInFaces, 
        int Face);

    /**
     * @brief Calculates the distance from a face (given by its nodes) and the given point, normalized 
     * to the size of the face side.
     * @param Nodes The nodes of the geometry
     * @param NodesInFaces matrix containing the index of the nodes of the geometry that belong to each faces
     * @return Distance from the face to the specified point
     * */
    static double NormalizedDistanceToQEF(
        PointsArrayType& rNodes,
        const DenseMatrix<unsigned int>& rNodesInFaces, 
        const array_1d<double,3>& rPoint, 
        int Face);

    /**
     * @brief Calculates the distance from a face (given by its nodes) and the given point, normalized 
     * to the size of the face side.
     * @param Nodes The nodes of the geometry
     * @param NodesInFaces matrix containing the index of the nodes of the geometry that belong to each faces
     * @return Distance from the face to the specified point
     * */
    static double NormalizedDistanceToQEF(
        GeometryType& rFace, 
        const array_1d<double,3>& rPoint, 
        int Face);

}; /* Class VoxelInsideVolumeQEF */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/
