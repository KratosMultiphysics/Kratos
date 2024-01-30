//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "testing/testing.h"

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
// forward declaring ModelPart and Model to be avoid including heavy header here
class ModelPart;
class Model;

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
     * @param rEntityName The entity name considered
     * @param Initialize If initialize the entities
     * @param Elements If create elements or conditions
     */
    void KRATOS_API(KRATOS_CORE) Create2DGeometry(
        ModelPart& rModelPart,
        const std::string& rEntityName = "Element2D3N",
        const bool Initialize = true,
        const bool Elements = true
        );

    /**
     * @brief This method creates a pure (Element) simple geometry in 2D (triangles)
     * @param rModelPart Reference to the ModelPart containing the problem
     */
    void KRATOS_API(KRATOS_CORE) CreateTestModelPartTriangle2D3N(ModelPart& rModelPart);

    /**
     * @brief This method creates a simple geometry in 2D (quadrilaterals)
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rEntityName The entity name considered
     * @param Initialize If initialize the entities
     * @param Elements If create elements or conditions
     */
    void KRATOS_API(KRATOS_CORE) Create2DQuadrilateralsGeometry(
        ModelPart& rModelPart, 
        const std::string& rEntityName = "Element2D4N",
        const bool Initialize = true,
        const bool Elements = true
        );

    /**
     * @brief This method creates a simple geometry in 3D (tetrahedra)
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rElementName The element name considered
     * @param Initialize If initialize the elements
     */
    void KRATOS_API(KRATOS_CORE) Create3DGeometry(
        ModelPart& rModelPart,
        const std::string& rElementName = "Element3D4N",
        const bool Initialize = true
        );

    /**
     * @brief This method creates a pure (Element) simple geometry in 3D (tetrahedra)
     * @param rModelPart Reference to the ModelPart containing the problem
     */
    void KRATOS_API(KRATOS_CORE) CreateTestModelPartTetrahedra3D4N(ModelPart& rModelPart);

    /**
     * @brief This method creates a simple geometry in 3D (hexahedra)
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rElementName The element name considered
     * @param Initialize If initialize the elements
     */
    void KRATOS_API(KRATOS_CORE) Create3DHexahedraGeometry(
        ModelPart& rModelPart,
        const std::string& rElementName = "Element3D8N",
        const bool Initialize = true
        );

    /**
     * @brief This method creates a simple geometry in 3D (tetrahedra quadratic)
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rElementName The element name considered
     * @param Initialize If initialize the elements
     */
    void KRATOS_API(KRATOS_CORE) Create3DQuadraticGeometry(
        ModelPart& rModelPart,
        const std::string& rElementName = "Element3D10N",
        const bool Initialize = true
        );

    /**
     * @brief This method creates a simple geometry sphere of triangles
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param rConditionName The condition name considered
     * @param Radius The radius of the sphere
     * @param rCenter The center of the sphere
     */
    void KRATOS_API(KRATOS_CORE) CreateSphereTriangularMesh(
        ModelPart& rModelPart,
        const std::string& rConditionName = "SurfaceCondition3D3N",
        const double Radius = 0.25,
        const std::array<double, 3>& rCenter = {0.0, 0.0, 0.0}
        );

    /**
     * @brief Create a cube skin model part.
     * @param rCurrentModel The current model.
     * @param HalfX The half-length of the cube in the X-direction.
     * @param HalfY The half-length of the cube in the Y-direction.
     * @param HalfZ The half-length of the cube in the Z-direction.
     * @param rDataCommunicator The data communicator.
     * @return ModelPart& The created cube skin model part.
     */
    KRATOS_API(KRATOS_CORE) ModelPart& CreateCubeSkinModelPart(
        Model& rCurrentModel,
        const double HalfX = 0.6,
        const double HalfY = 0.9,
        const double HalfZ = 0.3,
        const DataCommunicator& rDataCommunicator = Testing::GetDefaultDataCommunicator()
        );

    /**
     * @brief Create a cube model part.
     * @param rCurrentModel The current model.
     * @param rDataCommunicator The data communicator.
     * @return The created cube model part.
     */
    KRATOS_API(KRATOS_CORE) ModelPart& CreateCubeModelPart(
        Model& rCurrentModel,
        const DataCommunicator& rDataCommunicator = Testing::GetDefaultDataCommunicator()
        );

}; // namespace CppTestsUtilities
}  // namespace Kratos
