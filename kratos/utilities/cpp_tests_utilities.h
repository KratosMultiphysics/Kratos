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

#if !defined(KRATOS_CPP_TESTS_UTILITIES)
#define KRATOS_CPP_TESTS_UTILITIES

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
 * @namespace CppTestsUtilities
 * @ingroup KratosCore
 * @brief This namespace includes utilities for simplifying the deploy of C++ tests
 * @details The following method are implemented :
 * - Create2DGeometry: Creates a simple mesh of triangles 
 * - Create3DGeometry: Creates a simple mesh of tetrahedra 
 * @author Vicente Mataix Ferrandiz
 */
namespace CppTestsUtilities
{
    /**
     * @brief This method creates a simple geometry in 2D (triangles) 
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rElementName The element name considered
     */
    void KRATOS_API(KRATOS_CORE) Create2DGeometry(
        ModelPart& rModelPart, 
        const std::string& rElementName
        );

    /**
     * @brief This method creates a simple geometry in 3D (tetrahedra) 
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rElementName The element name considered
     */
    void KRATOS_API(KRATOS_CORE) Create3DGeometry(
        ModelPart& rModelPart, 
        const std::string& rElementName
        );

}; // namespace CppTestsUtilities
}  // namespace Kratos
#endif /* KRATOS_CPP_TESTS_UTILITIES defined */
