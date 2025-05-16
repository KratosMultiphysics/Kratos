//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortes
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "geometries/geometry.h"
#include "geometries/bounding_box.h"

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

// Forward Declaration class GeometricalObject (using pointers)
class GeometricalObject;

/**
 * @brief Utilities to compute the minimum error point in a 3D voxel intersected by a triangle mesh
 * This implementation is based on the algorithm explained here: https://www.mattkeeter.com/projects/qef/
 * @author Ariadna Cortes
 * @ingroup KratosCore
 */
class KRATOS_API(KRATOS_CORE) QuadraticErrorFunction
{
public:

    ///@name Type Definitions
    ///@{

    using NodeType = Node;
    using NodePtrType = Node::Pointer;
    using GeometryType = Geometry<NodeType>;
    using GeometryPtrType = GeometryType::Pointer;
    using GeometryArrayType = GeometryType::GeometriesArrayType;
    using PointsArrayType = GeometryType::PointsArrayType;

    /// Pointer definition of VoxelInsideVolume
    KRATOS_CLASS_POINTER_DEFINITION( QuadraticErrorFunction );

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     */
    QuadraticErrorFunction(){}

    /// Destructor
    virtual ~QuadraticErrorFunction(){}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Finds the QuadraticErrorFunction point of a voxel
     * @param rVoxel references to the voxel whose x-point will be calculated
     * @param rTriangles references to the triangles which intersect the voxel at some edge.
     * @return The QuadraticErrorFunction point (x,y,z)
     */
    static array_1d<double,3> QuadraticErrorFunctionPoint (
        const GeometryType& rVoxel,
        const GeometryArrayType& rTriangles
        );

    /**
     * @brief Finds the QuadraticErrorFunction point of a voxel
     * @param rVoxel references to the voxel whose x-point will be calculated
     * @param rTriangles references to the triangles which intersect the voxel at some edge.
     * @return The QuadraticErrorFunction point (x,y,z)
     */
    static array_1d<double,3> QuadraticErrorFunctionPoint (
        const BoundingBox<Point>& rBox,
        const std::vector<GeometricalObject*>& rTriangles
        );

    /**
     * @brief Calculates the normal vector to the surface of a 3D triangle
     * @param rTriangle reference to the triangle
     * @return Normal vector (x,y,z)
     */
    static array_1d<double,3> CalculateNormal(const GeometryType& rTriangle);

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

    static Point FirstEnd(int i, const BoundingBox<Point>& rBox);
    static Point SecondEnd(int i, const BoundingBox<Point>& rBox);

    ///@}
}; /* Class VoxelInsideVolumeUtility */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/