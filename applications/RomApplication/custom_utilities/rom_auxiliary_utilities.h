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
#include "processes/find_nodal_neighbours_process.h"

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
     * @brief Sets the HROM model part including neighboring entities based on the nodal weights
     *
     * Provided an origin model part and a parameters object containing the HROM weights,
     * this method stores not only the elements and conditions required by the HROM in the destination
     * model part, but also includes neighboring elements and conditions of the nodes involved.
     * The resulting model part features the same submodelpart hierarchy that the origin has, but is
     * augmented with the neighboring elements and conditions.
     *
     * @param HRomWeights Parameters object containing the HROM elemental and condition weights
     * @param rOriginModelPart Origin model part (this is likely the computing model part)
     * @param rHRomComputingModelPart Destination model part to store the HROM mesh augmented with neighbours
     */
    static void SetHRomComputingModelPartWithNeighbours(
        const Parameters HRomWeights,
        ModelPart& rOriginModelPart,
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
     * @brief Return the decremented element ids(-1 to account for numpy indexing) of the missing HROM condition parents
     * This function loops the HROM conditions in the HROM weights and searches for their parents in the
     * provided model part. Once these are found, it is checked that they are not already added and returned.
     * Note that this functions requires the NEIGHBOUR_ELEMENTS to be already computed.
     * @param rModelPart Complete model part (all elements and conditions)
     * @param rHRomWeights Map containing the original HROM conditions and elements weights
     * @return std::vector<IndexType> List containing the decremented element IDs (-1 to account for numpy indexing) of the missing parent elements
     */
    static std::vector<IndexType> GetHRomConditionParentsIds(
        const ModelPart& rModelPart,
        const std::map<std::string, std::map<IndexType, double>>& rHRomWeights);

    /**
     * @brief Retrieve the decremented (-1 to account for numpy indexing) IDs of elements neighboring nodes in a given sub-model part but not present in HRom weights.
     *
     * This function iterates over all the nodes in the provided sub-model part (`rGivenModelPart`) and collects
     * the decremented IDs (-1 to account for numpy indexing) of the elements neighboring each node. The neighboring elements are determined using the
     * 'NEIGHBOUR_ELEMENTS' value attached to each node. The function then checks if these elements are already
     * present in `rHRomWeights`. If not, their decremented (-1 to account for numpy indexing) IDs are added to the return vector. It's important to note that
     * this function assumes that the 'NEIGHBOUR_ELEMENTS' values are already computed for the nodes in the model part.
     *
     * @param rModelPart The complete model part which houses all the elements.
     * @param rGivenModelPart The sub-model part with nodes for which neighboring elements should be fetched.
     * @param rHRomWeights Map containing the original HROM conditions and elements weights.
     * @return std::vector<IndexType> A list of unique neighboring decremented element IDs (-1 to account for numpy indexing) not already present in HRom weights.
     */
    static std::vector<IndexType> GetNodalNeighbouringElementIdsNotInHRom(
        ModelPart& rModelPart,
        ModelPart& rGivenModelPart,
        const std::map<std::string, std::map<IndexType, double>>& rHRomWeights);

    /**
     * @brief Retrieve the decremented (-1 to account for numpy indexing) IDs of elements neighboring nodes in a given sub-model part.
     *
     * This function iterates over all the nodes in the provided sub-model part (`rGivenModelPart`) and collects
     * the decremented (-1 to account for numpy indexing) IDs of the elements neighboring each node. The neighboring elements are determined using the
     * 'NEIGHBOUR_ELEMENTS' value attached to each node. The function automatically ensures that each element ID
     * is unique, avoiding duplicates in the return vector. It's important to note that this function assumes that
     * the 'NEIGHBOUR_ELEMENTS' values are already computed for the nodes in the model part.
     *
     * @param rModelPart The complete model part which houses all the elements.
     * @param rGivenModelPart The sub-model part with nodes for which neighboring elements should be fetched.
     * @return std::vector<IndexType> A list of unique neighboring decremented element IDs (-1 to account for numpy indexing).
     */
    static std::vector<IndexType> GetNodalNeighbouringElementIds(
        ModelPart& rModelPart,
        ModelPart& rGivenModelPart);

    /**
     * @brief Identifies condition decremented (-1 to account for numpy indexing) IDs from a given ModelPart that are not in the HROM weights
     * This function iterates through the conditions in the provided ModelPart, checks if their IDs exist in the HROM weights,
     * and includes them in a list if they are missing. The decremented (-1 to account for numpy indexing) IDs of the absent conditions are returned.
     * @param rModelPart Complete hrom model part (all elements and conditions)
     * @param rModelPartWithConditionsToInclude ModelPart with conditions that should be checked against HROM weights
     * @return std::vector<IndexType> List containing the decremented condition IDs (-1 to account for numpy indexing) of the conditions absent in HROM weights
     */
    static std::vector<IndexType> GetConditionIdsNotInHRomModelPart(
        const ModelPart& rModelPartWithConditionsToInclude,
        std::map<std::string, std::map<IndexType, double>>& rHRomWeights);

    /**
     * @brief Identifies element decremented (-1 to account for numpy indexing) IDs from a given ModelPart that are not in the HROM weights
     * This function iterates through the elements in the provided ModelPart, checks if their decremented (-1 to account for numpy indexing) IDs exist in the HROM weights,
     * and includes them in a list if they are missing. The decremented (-1 to account for numpy indexing) IDs of the absent elements are returned.
     * @param rModelPartWithElementsToInclude ModelPart with elements that should be checked against HROM weights
     * @param rHRomWeights Map containing the original HROM elements and conditions weights
     * @return std::vector<IndexType> List containing decremented element IDs (-1 to account for numpy indexing) of the elements absent in HROM weights
     */
    static std::vector<IndexType> GetElementIdsNotInHRomModelPart(
        const ModelPart& rModelPartWithElementsToInclude,
        std::map<std::string, std::map<IndexType, double>>& rHRomWeights);

    /**
     * @brief Return the a minimum condition for each HROM submodelpart
     * This function loops the HROM mesh submodelparts and checks if there is at least a minimum condition
     * for each submodelpart in the HROM condition weights. If there are no conditions it adds a null weight
     * condition to the HROM weights. This might be required for the HROM submodelparts visualization or BC imposition
     * @param rHRomConditionWeights Map containing the original HROM conditions and elements weights
     * @return std::vector<IndexType> List containing the decremented condition IDs (-1 to account for numpy indexing) of the missing conditions to be added to the weights
     */
    static std::vector<IndexType> GetHRomMinimumConditionsIds(
        const ModelPart& rModelPart,
        const std::map<IndexType, double>& rHRomConditionWeights);


    /**
    * @brief Retrieve the decremented IDs (-1 to account for numpy indexing) of elements present in the provided model part.
    *
    * This function iterates over all the elements in the given model part (`rModelPart`) and collects
    * their IDs after decrementing each by one. The decremented IDs are then added to a vector which
    * is returned to the caller.
    *
    * @param rModelPart The model part containing the elements whose IDs are to be fetched.
    * @return std::vector<IndexType> A list of decremented element IDs (-1 to account for numpy indexing) from the provided model part.
    */
    static std::vector<IndexType> GetElementIdsInModelPart(
        const ModelPart& rModelPart
    );

    /**
    * @brief Retrieve the decremented IDs (-1 to account for numpy indexing) of conditions present in the provided model part.
    *
    * This function iterates over all the conditions in the given model part (`rModelPart`) and collects
    * their IDs after decrementing each by one. The decremented IDs are then added to a vector which
    * is returned to the caller.
    *
    * @param rModelPart The model part containing the conditions whose IDs are to be fetched.
    * @return std::vector<IndexType> A list of decremented condition IDs (-1 to account for numpy indexing) from the provided model part.
    */
    static std::vector<IndexType> GetConditionIdsInModelPart(
        const ModelPart& rModelPart
    );

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

    /**
     * @brief Obtain the JPhi elemental matrix for a particular element. 
     * JPhi represents the projection of the Jacobian onto the ROM_BASIS.
     * @param rJPhiElemental The matrix to store the result in. Must have the appropriate size already.
     * @param rDofs The set of degrees of freedom (DoFs) of the element.
     * @param rJPhi The JPhi matrix, from which rows are extracted according to the equation ID of each DoF.
     * 
     * This function loops over all the DoFs for the given element. For each DoF, it uses its equation ID to extract a 
     * corresponding row from the rJPhi matrix, which is then stored in the corresponding row of rJPhiElemental.
     */
    static void GetJPhiElemental(
        Matrix &rJPhiElemental,
        const Element::DofsVectorType& rDofs,
        const Matrix &rJPhi);
        
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
