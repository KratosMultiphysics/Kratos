//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Suneth Warnakulasuriya
//

#if !defined(KRATOS_ROM_AUXILIARY_UTILITIES_H)
#define KRATOS_ROM_AUXILIARY_UTILITIES_H

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/key_hash.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "modified_shape_functions/modified_shape_functions.h"

// Application includes
#include "rom_application_variables.h"


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(ROM_APPLICATION) RomAuxiliaryUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using SizeType = std::size_t;

    using IndexType = std::size_t;

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using GeometryPointerType = Geometry<NodeType>::Pointer;

    using NodesPointerSetType = ModelPart::NodesContainerType;

    using ElementFacesMapType = std::unordered_map<
        std::vector<IndexType>,
        std::pair<bool, GeometryPointerType>,
        KeyHasherRange<std::vector<IndexType>>,
        KeyComparorRange<std::vector<IndexType>>>;

    ///@}
    ///@name Static Operations
    ///@{

    /**
     * @brief Sets the HROM model part from the nodal weights
     * Provided an origin model part and a parameters object containing the HROM weights,
     * this method stores the elements and conditions required by the HROM in the destination
     * model part. The new model part features the same submodelpart hierarchy that the origin has.
     * @param HRomWeights Parameters object containing the HROM elemental and condition weights
     * @param rOriginModelPart Origin model part (this is likely the computing model part)
     * @param rHRomComputingModelPart Destination model part to store the HROM mesh
     */
    static void SetHRomComputingModelPart(
        const Parameters HRomWeights,
        const ModelPart& rOriginModelPart,
        ModelPart& rHRomComputingModelPart);

    /**
     * @brief Sets the HROM skin visualization model part for a volumetric body
     * This function detects the skin of the origin modelpart and creates the corresponding skin
     * entities (conditions) which are stored in an auxiliary modelpart to visualize the HROM solution
     * It is important to mention that this function only works with volumetric (a.k.a. solid) bodies
     * @param HRomWeights
     * @param rOriginModelPart
     * @param rHRomVisualizationModelPart
     */
    static void SetHRomVolumetricVisualizationModelPart(
        const ModelPart& rOriginModelPart,
        ModelPart& rHRomVisualizationModelPart);

    /**
     * @brief Return the missing HROM condition parents element ids
     * This function loops the HROM conditions in the HROM weights and searches for their parents in the
     * provided model part. Once these are found, it is checked that they are not already added and returned.
     * Note that this functions requires the NEIGHBOUR_ELEMENTS to be already computed.
     * @param rModelPart Complete model part (all elements and conditions)
     * @param rHRomWeights Map containing the original HROM conditions and elements weights
     * @return std::vector<IndexType> List containing the ids of the missing parent elements
     */
    static std::vector<IndexType> GetHRomConditionParentsIds(
        const ModelPart& rModelPart,
        const std::map<std::string, std::map<IndexType, double>>& rHRomWeights);

    /**
     * @brief Identifies condition IDs from a given ModelPart that are not in the HROM weights
     * This function iterates through the conditions in the provided ModelPart, checks if their IDs exist in the HROM weights,
     * and includes them in a list if they are missing. The IDs of the absent conditions are returned.
     * @param rModelPart Complete hrom model part (all elements and conditions)
     * @param rModelPartWithConditionsToInclude ModelPart with conditions that should be checked against HROM weights
     * @param rHRomWeights Map containing the original HROM conditions and elements weights
     * @return std::vector<IndexType> List containing the IDs of the conditions absent in HROM weights
     */
    static std::vector<IndexType> GetConditionIdsNotInHRomModelPart(
        const ModelPart& rModelPart,
        const ModelPart& rModelPartWithConditionsToInclude,
        std::map<std::string, std::map<IndexType, double>>& rHRomWeights);

    /**
     * @brief Identifies element IDs from a given ModelPart that are not in the HROM weights
     * This function iterates through the elements in the provided ModelPart, checks if their IDs exist in the HROM weights,
     * and includes them in a list if they are missing. The IDs of the absent elements are returned.
     * @param rModelPart Complete hrom model part (all elements and conditions)
     * @param rModelPartWithElementsToInclude ModelPart with elements that should be checked against HROM weights
     * @param rHRomWeights Map containing the original HROM elements and conditions weights
     * @return std::vector<IndexType> List containing the IDs of the elements absent in HROM weights
     */
    static std::vector<IndexType> GetElementIdsNotInHRomModelPart(
        const ModelPart& rModelPart,
        const ModelPart& rModelPartWithElementsToInclude,
        std::map<std::string, std::map<IndexType, double>>& rHRomWeights);

    /**
     * @brief Return the a minimum condition for each HROM submodelpart
     * This function loops the HROM mesh submodelparts and checks if there is at least a minimum condition
     * for each submodelpart in the HROM condition weights. If there are no conditions it adds a null weight
     * condition to the HROM weights. This might be required for the HROM submodelparts visualization or BC imposition
     * @param rModelPart Complete model part (all elements and conditions)
     * @param rHRomConditionWeights Map containing the original HROM conditions and elements weights
     * @return std::vector<IndexType> List containing the ids of the missing conditions to be added to the weights
     */
    static std::vector<IndexType> GetHRomMinimumConditionsIds(
        const ModelPart& rModelPart,
        const std::map<IndexType, double>& rHRomConditionWeights);

    /**
     * @brief Project the ROM solution increment onto the nodal basis
     * For a given model part this function takes the ROM_SOLUTION_INCREMENT, which is assumed to be
     * stored in the root model part, and projects it with the nodal ROM_BASIS to calculate the solution increment
     * @param rRomVariableNames Vector containing the names of the ROM variables
     * @param rModelPart Model part onto which the projection is to be performed
     */
    static void ProjectRomSolutionIncrementToNodes(
        const std::vector<std::string> &rRomVariableNames,
        ModelPart &rModelPart);

    /**
     * @brief Obtain the elemental basis matrix for a particular element.
     * @param rPhiElemental The matrix to store the result in. Must have the appropiate size already.
     * @param rDofs The set of dofs of the element.
     * @param rGeom The geometry of the element.
     * @rVarToRowMapping A map from each variables's key to its row in the basis matrix.
     */
    static void GetPhiElemental(
        Matrix &rPhiElemental,
        const Element::DofsVectorType& rDofs,
        const Element::GeometryType& rGeom,
        const std::unordered_map<Kratos::VariableData::KeyType, Matrix::size_type>& rVarToRowMapping);


    /**
     * @brief Obtain the left elemental basis (Psi) matrix for a particular element.
     * @param rPsiElemental The matrix to store the result in. Must have the appropiate size already.
     * @param rDofs The set of dofs of the element.
     * @param rGeom The geometry of the element.
     * @param rVarToRowMapping A map from each variables's key to its row in the basis matrix.
     */
    static void GetPsiElemental(
        Matrix &rPsiElemental,
        const Element::DofsVectorType& rDofs,
        const Element::GeometryType& rGeom,
        const std::unordered_map<Kratos::VariableData::KeyType, Matrix::size_type>& rVarToRowMapping);

    ///@}

    private:

    ///@name Private Static Operations
    ///@{

    static void RecursiveHRomModelPartCreation(
        const NodesPointerSetType& rNodesSet,
        const std::vector<Element::Pointer>& rElementsVector,
        const std::vector<Condition::Pointer>& rConditionsVector,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart);

    static void RecursiveHRomMinimumConditionIds(
        const ModelPart& rModelPart,
        const std::map<IndexType, double>& rHRomConditionWeights,
        std::vector<IndexType>& rMinimumConditionsIds);

    static void RecursiveVisualizationSubModelPartCreation(
        const ModelPart& rOriginSubModelPart,
        ModelPart& rDestinationModelPart);

    ///@}

};

///@}

} // namespace Kratos

#endif // KRATOS_ROM_AUXILIARY_UTILITIES_H
