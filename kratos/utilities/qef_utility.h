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
//                   Vicente Mataix Ferrandiz
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
 * @class QuadraticErrorFunction (Quadratic Error Function)
 * @ingroup KratosCore
 * @brief Utilities to compute the minimum error point in a 3D voxel intersected by a triangle mesh
 * @details This implementation is based on the algorithm explained here: https://www.mattkeeter.com/projects/qef/
 * The QuadraticErrorFunction class contains methods to find the quadratic error function point of a voxel
 * and to calculate the normal vector to the surface of a 3D triangle. It uses Kratos geometry types and
 * provides static methods for these calculations.
 * @see Geometry
 * @see BoundingBox
 * @author Ariadna Cortes
 */
class KRATOS_API(KRATOS_CORE) QuadraticErrorFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// The node type
    using NodeType = Node;

    // The node pointer type
    using NodePtrType = Node::Pointer;

    /// The geometry type
    using GeometryType = Geometry<NodeType>;

    /// The geometry pointer type
    using GeometryPtrType = GeometryType::Pointer;

    /// The geometry array type
    using GeometryArrayType = GeometryType::GeometriesArrayType;

    /// The points array type
    using PointsArrayType = GeometryType::PointsArrayType;

    /**
     * @brief Auxiliary classes to store the matrices and vectors needed to compute the quadratic error function point
     */
    struct AuxiliaryClasses {
        BoundedMatrix<double, 3, 1> MatCenter;
        BoundedMatrix<double, 3, 3> ATA = ZeroMatrix(3, 3);
        BoundedMatrix<double, 3, 1> ATB = ZeroMatrix(3, 1);

        array_1d<double, 3> Normal;
        array_1d<double, 3> Intersection;
        BoundedMatrix<double, 3, 1> MatNormal;
        BoundedMatrix<double, 1, 3> MatNormalTrans;
        BoundedMatrix<double, 1, 1> MatAux;
        GeometryType::CoordinatesArrayType AuxCenter = ZeroVector(3);
    };

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
    QuadraticErrorFunction() = default;

    /// Destructor
    virtual ~QuadraticErrorFunction() = default;

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

    ///@}
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
     * @brief Update ATA and ATB matrices
     * @param rAuxiliaryClasses The auxiliary classes containing the matrices and vectors needed to compute the quadratic error function point
     */
    static void UpdateATAATB(AuxiliaryClasses& rAuxiliaryClasses);

    /**
     * @brief Computes the first endpoint of an edge in a bounding box.
     * @param rPoint The first endpoint of the edge.
     * @param i The index of the edge.
     * @param rBox The bounding box from which the endpoint is computed.
     * @return The first endpoint of the specified edge.
     */
    static void FirstEnd(
        Point& rPoint,
        const int i,
        const BoundingBox<Point>& rBox
        );

    /**
     * @brief Computes the second endpoint of an edge in a bounding box.
     * @param rPoint The point to be computed.
     * @param i The index of the edge.
     * @param rBox The bounding box from which the endpoint is computed.
     * @return The second endpoint of the specified edge.
     */
    static void SecondEnd(
        Point& rPoint,
        const int i,
        const BoundingBox<Point>& rBox
        );

    /**
     * @brief Computes the point that minimizes the quadratic error function.
     * @details This function computes the solution for the quadratic error function by using the
     * pseudo-inverse of the ATA matrix. The solution is given by:
     * \f[
     * x = center + ATA^{-1} \cdot (ATB - ATA \cdot center)
     * \f]
     * where ATA^{-1} is computed via an SVD decomposition of the ATA matrix.
     * @param rATA A symmetric 3x3 matrix representing the accumulated outer product of normals (i.e., AᵀA).
     * @param rATB A 3x1 vector representing the accumulated product of normals and their corresponding scalar offsets.
     * @param rMatCenter A 3x1 vector representing the center point.
     * @return array_1d<double,3> The computed point that minimizes the quadratic error function.
     */
    static array_1d<double,3> ComputeQuadraticErrorFunctionPoint(
        const BoundedMatrix<double,3,3>& rATA,
        const BoundedMatrix<double,3,1>& rATB,
        const BoundedMatrix<double,3,1>& rMatCenter
        );

    ///@}
}; /* Class VoxelInsideVolumeUtility */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/