//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes

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
// forward declaring ModelPart and Point to be avoid including heavy header here
class ModelPart;
class Point;

/**
 * @namespace DelaunatorUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities using the library triangle
 * @author Vicente Mataix Ferrandiz
 */
namespace DelaunatorUtilities
{
    /**
     * @brief This method creates a triangle mesh from a model part of nodes
     * @param rModelPart The model of the problem to mesh
     */
    void KRATOS_API(KRATOS_CORE) CreateTriangleMeshFromNodes(ModelPart& rModelPart);

    /**
     * @brief This method returns the triangles connectivity from a list of coordinates (using triangle library)
     * @param rCoordinates The list of coordinates, first X, then Y, for each point of the point cloud
     * @return The connectivity vector
     */
    std::vector<std::size_t> KRATOS_API(KRATOS_CORE) ComputeTrianglesConnectivity(const std::vector<double>& rCoordinates);

    /**
     * @brief This method returns the triangles connectivity from a list of coordinates (using triangle library)
     * @param rPoints The list of points
     * @return The connectivity vector
     */
    std::vector<std::size_t> KRATOS_API(KRATOS_CORE) ComputeTrianglesConnectivity(const std::vector<Point>& rPoints);

    /**
     * @brief This methods does a Constrained Delaunay Triangularization from a list of coordinates and segments (using triangle library)
     * @param rCoordinates The list of coordinates, first X, then Y, for each point of the point cloud
     * @param rSegments The list of segments, each segment is determined in one std::array with from its i and j nodal ids
     * @param AreaConstraint If provided, imposes that value as a constraint on the maximum area
     * @return A pair containing in first position a list with the triangles connectivities and in second position a list with the x and y nodal coordinates 
     */
    std::pair<std::vector<std::size_t>, std::vector<double>> KRATOS_API(KRATOS_CORE) ComputeTrianglesConnectivity(
        const std::vector<double>& rCoordinates,
        const std::vector<std::array<double,2>>& rSegments,
        const double AreaConstraint = 0);

}; // namespace DelaunatorUtilities
}  // namespace Kratos
