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

    using NodeType = Node<3>;

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

    ///@}

    static void AppendConditionParentsToHRomWeights(
        const ModelPart& rModelPart,
        std::map<std::string, std::map<IndexType, double>>& rHromWeights);

    private:

    ///@name Private Static Operations
    ///@{

    static void RecursiveHRomModelPartCreation(
        const NodesPointerSetType& rNodesSet,
        const std::vector<Element::Pointer>& rElementsVector,
        const std::vector<Condition::Pointer>& rConditionsVector,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart);

    ///@}

};

///@}

} // namespace Kratos

#endif // KRATOS_ROM_AUXILIARY_UTILITIES_H
