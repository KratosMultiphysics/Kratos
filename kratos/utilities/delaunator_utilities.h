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
    /**
     * @brief This method creates a triangle mesh from a model part of nodes
     * @param rModelPart The model of the problem to mesh
     */
    void KRATOS_API(KRATOS_CORE) CreateTriangleMeshFromNodes(ModelPart& rModelPart);

}; // namespace DelaunatorUtilities
}  // namespace Kratos
#endif /* KRATOS_DELAUNATOR_UTILITIES defined */
