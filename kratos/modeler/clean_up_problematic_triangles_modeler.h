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
#include "modeler/modeler.h"
#include "includes/define_registry.h"
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
 * @class CleanUpProblematicTrianglesModeler
 * @ingroup KratosCore
 * @brief This class is responsible for cleaning up problematic geometries in a mesh, specifically null area geometries.
 * @details The class identifies and removes degenerated triangles where two nodes are in the same position. It ensures that the mesh remains water-tight by reconnecting the sides and renumbering the nodes and entities.
 * @note This is a solution to clean up problematic geometries in the mesh that for example can be found in OBJ files.
 * @see Modeler
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) CleanUpProblematicTrianglesModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the base type
    using BaseType = Modeler;

    // Define GeometryType
    using GeometryType = Geometry<Node>;

    /// Pointer definition of CleanUpProblematicTrianglesModeler
    KRATOS_CLASS_POINTER_DEFINITION(CleanUpProblematicTrianglesModeler);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    CleanUpProblematicTrianglesModeler() : Modeler() {};

    /**
     * @brief Constructor with Model
     * @param Model The model to be used.
     * @param ModelerParameters The parameters for the modeler.
     */
    CleanUpProblematicTrianglesModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters())
        : BaseType(rModel, ModelerParameters),
          mpModelPart(&rModel.GetModelPart(ModelerParameters["model_part_name"].GetString()))
    {
        mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    }

    /**
     * @brief Destructor.
     */
    ~CleanUpProblematicTrianglesModeler() override
    {
    }

    /**
     * @brief Creates the Modeler Pointer
     * @param Model The model to be used.
     * @param Parameters The parameters for the modeler.
     */
    Modeler::Pointer Create(
        Model& rModel,
        const Parameters ModelParameters
        ) const override
    {
        return Kratos::make_shared<CleanUpProblematicTrianglesModeler>(rModel, ModelParameters);
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method clean up the problematic geometries (null area geometries) in the mesh.
     * @details For the moment it will be assumed that the problematic geometries are triangles with two nodes in the same position. This means several this:
     * - The triangle is degenerated, one of the nodes is redundant and should be cleaned up.
     * - The triangle must be removed as well. The sides must be reconnected, so the mesh is water-tight.
     * - The nodes and entities must be renumbered
     * @param rThisModelPart Reference to the model part to clean up into.
     * @param rEntityType The entity type to create in the model part. Can be "element", "condition" or "geometry".
     * @param FirstNodeId The first node ID to create in the model part.
     * @param FirstELementId The first element ID to create in the model part.
     * @param FirstConditionId The first condition ID to create in the model part.
     * @param AreaTolerance The tolerance to consider a geometry as problematic.
     */
    static void CleanUpProblematicGeometriesInMesh(
        ModelPart& rThisModelPart,
        const std::string& rEntityType = "element",
        const IndexType FirstNodeId = 1,
        const IndexType FirstElementId = 1,
        const IndexType FirstConditionId = 1,
        const double AreaTolerance = 1.0e-4
        );

    ///@}
    ///@name Modeler Stages at Initialize
    ///@{

    /**
     * @brief Convert the geometry model or import analysis suitable models.
     */
    void SetupModelPart() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"({
            "model_part_name"    : "PLEASE_DEFINE_A_NAME",
            "entity_type"        : "element",
            "first_node_id"      : 1,
            "first_element_id"   : 1,
            "first_condition_id" : 1,
            "area_tolerance"     : 1.0e-4,
            "echo_level"         : 0
        })");
        return default_parameters;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CleanUpProblematicTrianglesModeler";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    KRATOS_REGISTRY_ADD_PROTOTYPE("Modelers.KratosMultiphysics", Modeler, CleanUpProblematicTrianglesModeler)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Modelers.All", Modeler, CleanUpProblematicTrianglesModeler)

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart* mpModelPart = nullptr;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Get the nodes of the model part.
     * @param rThisModelPart Reference to the model part.
     * @return The nodes of the model part.
     * @tparam TEntityType The entity type to get from the model part. Can be Element,  Condition or Geometry.
     */
    template <typename TEntityType>
    static void RemoveEntitiesAndNodes(ModelPart& rThisModelPart);

    /**
     * @brief Clean up the problematic geometries (null area geometries) in the mesh.
     * @details For the moment it will be assumed that the problematic geometries are triangles with two nodes in the same position. This means several this:
     * - The triangle is degenerated, one of the nodes is redundant and should be cleaned up.
     * - The triangle must be removed as well. The sides must be reconnected, so the mesh is water-tight.
     * - The nodes and entities must be renumbered
     * @param rThisModelPart Reference to the model part to clean up into.
     * @param rEntityType The entity type to create in the model part. Can be "element", "condition" or "geometry".
     * @param FirstNodeId The first node ID to create in the model part.
     * @param FirstEntityId The first entity ID to create in the model part.
     * @param AreaTolerance The tolerance to consider a geometry as problematic.
     * @tparam TEntityType The entity type to create in the model part. Can be Element, Condition or Geometry.
     */
    template <typename TEntityType>
    static void CleanUpProblematicGeometries(
        ModelPart& rThisModelPart,
        const IndexType FirstNodeId = 1,
        const IndexType FirstEntityId = 1,
        const double AreaTolerance = 1.0e-4
        );

    /**
    * @brief Computes the distance between two nodes.
    * @details This function calculates the Euclidean distance between two nodes using their coordinates. The distance is computed as the norm of the difference between the two position vectors.
    * @param rNode1 The first node.
    * @param rNode2 The second node.
    * @return The distance between the two nodes.
    */
    static double ComputeDistance(
        const Node& rNode1,
        const Node& rNode2
        );

    /**
    * @brief Computes the squared area of a triangle.
    * @details This function calculates the squared area of a triangle defined by three points in a 3D space. The area is calculated using the semi-perimeter formula.
    * @param rGeometry The geometry of the triangle containing three points.
    * @return The squared area of the triangle.
    */
    static double ComputeSquaredArea(const GeometryType& rGeometry);

    /**
    * @brief Computes the number of triangles with null area in a model part.
    * @details This function iterates over the entities in the given model part and counts the number of triangles that have a squared area less than the specified reference area.
    * @param rThisModelPart The model part containing the geometrical entities.
    * @param RefArea The reference area to compare against.
    * @param AreaTolerance The tolerance for determining if a triangle's area is considered null.
    * @return The number of triangles with null area.
    */
    template <typename TEntityType>
    static std::size_t ComputeNullAreaTriangles(
        ModelPart& rThisModelPart,
        const double RefArea,
        const double AreaTolerance
        );

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class CleanUpProblematicTrianglesModeler

// Helper function definitions specialization
template <>
void CleanUpProblematicTrianglesModeler::RemoveEntitiesAndNodes<Element>(ModelPart& rThisModelPart);

template <>
void CleanUpProblematicTrianglesModeler::RemoveEntitiesAndNodes<Condition>(ModelPart& rThisModelPart);

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >>(
    std::istream& rIStream,
    CleanUpProblematicTrianglesModeler& rThis);

/// output stream function

inline std::ostream & operator <<(
    std::ostream& rOStream,
    const CleanUpProblematicTrianglesModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


} // namespace Kratos.


