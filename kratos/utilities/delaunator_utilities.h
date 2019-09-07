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

#if !defined(KRATOS_DELAUNATOR_UTILITIES)
#define KRATOS_DELAUNATOR_UTILITIES

// System includes

// External includes
#define REAL double

#if !defined(KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED)
#define  KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED
#include "triangle.h"
#endif

// Project includes
#include "includes/model_part.h"

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
 * @namespace DelaunatorUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities using the library delaunator-cpp
 * @author Vicente Mataix Ferrandiz
 */
namespace DelaunatorUtilities
{
    extern "C" {
        void triangulate(char *, struct triangulateio *, struct triangulateio *,struct triangulateio *);
    }

    /**
     * @brief This method creates a triangle mesh from a model part of nodes
     * @param rModelPart The model of the problem to mesh
     */
    void KRATOS_API(KRATOS_CORE) CreateTriangleMeshFromNodes(ModelPart& rModelPart);

    /**
     * @brief This method returns the triangles connectivity from a list of coordinates
     * @param rCoordinates The list of coordinates, first X, then Y, for each point of the point cloud
     * @return The connectivity vector
     */
    std::vector<std::size_t> KRATOS_API(KRATOS_CORE) ComputeTrianglesConnectivity(const std::vector<double>& rCoordinates);

    /**
     * @brief This method returns the triangles connectivity from a list of coordinates (using triangle library)
     * @param rCoordinates The list of coordinates, first X, then Y, for each point of the point cloud
     * @return The connectivity vector
     */
    std::vector<std::size_t> KRATOS_API(KRATOS_CORE) ComputeTrianglesConnectivityWithTriangle(const std::vector<double>& rCoordinates);

}; // namespace DelaunatorUtilities
}  // namespace Kratos
#endif /* KRATOS_DELAUNATOR_UTILITIES defined */
