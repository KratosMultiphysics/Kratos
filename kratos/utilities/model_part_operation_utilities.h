//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/key_hash.h"
#include "includes/model_part.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) ModelPartOperationUtilities {
public:
    ///@name Type definitions
    ///@{

    using CNodePointersType = std::vector<ModelPart::NodeType const*>;

    template<class KeyType, class ValueType>
    using RangedKeyMapType = std::unordered_map<KeyType, ValueType, KeyHasherRange<KeyType>, KeyComparorRange<KeyType>>;

    ///@}
    ///@name Static operations
    ///@{

    /**
     * @brief Checks the validity of main model part and check model parts to be used within operation utilities.
     *
     * This method checks whether all the nodes, and geometries in @ref rCheckModelParts are present in the
     * @ref rMainModelPart.
     *
     * If @ref ThrowError is shown, a meaning full error is shown.
     *
     * @param rMainModelPart        Main model part.
     * @param rCheckModelParts      List of check model parts.
     * @param ThrowError            To throw an error or not.
     * @return true                 If the given model parts are valid.
     * @return false                If the given model parts are invalid.
     */
    static bool CheckValidityOfModelPartsForOperations(
        const ModelPart& rMainModelPart,
        const std::vector<ModelPart const*>& rCheckModelParts,
        const bool ThrowError = false);

    /**
     * @brief Merge given union model parts to one model part.
     *
     * This method merges all the nodes, geometries in @ref rUnionModelParts in to one model part
     * given by @ref rOutputSubModelPartName. The output model part is a sub-model part of
     * @ref rMainModelPart.
     *
     * The merging is carried out by finding all the nodes, geometries in @ref rUnionModelParts
     * within the @ref rMainModelPart, and the output model part is populated with the found
     * entities in @ref rMainModelPart.
     *
     * If @ref AddNeighbourEntities is true, then all the neighbours of the merged model part
     * is found and added as well.
     *
     * Please make sure to check validity of model parts using @see CheckValidityOfModelPartsForOperations.
     *
     * @param rOutputSubModelPartName       Output sub model part name.
     * @param rMainModelPart                Main Model part.
     * @param rUnionModelParts              Union model parts list to be merged.
     * @param AddNeighbourEntities          To add or not the neighbours.
     * @return ModelPart&                   Output sub model part.
     */
    static ModelPart& Merge(
        const std::string& rOutputSubModelPartName,
        ModelPart& rMainModelPart,
        const std::vector<ModelPart const*>& rUnionModelParts,
        const bool AddNeighbourEntities);

    /**
     * @brief Substract given substraction model parts from main model part.
     *
     * This method substract all the nodes, geometries in @ref rSubstractionModelParts from model part
     * given by @ref rMainModelPart and output to sub model part given by name @ref rOutputSubModelPartName.
     * The output model part is a sub-model part of @ref rMainModelPart.
     *
     * The substraction is carried out by finding all the nodes, geometries in @ref rSubstractionModelParts
     * within the @ref rMainModelPart, and the output model part is populated with the not found
     * entities in @ref rMainModelPart.
     *
     * If @ref AddNeighbourEntities is true, then all the neighbours of the substracted model part
     * is found and added as well.
     *
     * Please make sure to check validity of model parts using @see CheckValidityOfModelPartsForOperations.
     *
     * @param rOutputSubModelPartName       Output sub model part name.
     * @param rMainModelPart                Main Model part.
     * @param rSubstractionModelParts       Model parts to substract from @ref rMainModelPart.
     * @param AddNeighbourEntities          To add or not the neighbours.
     * @return ModelPart&                   Output sub model part.
     */
    static ModelPart& Substract(
        const std::string& rOutputSubModelPartName,
        ModelPart& rMainModelPart,
        const std::vector<ModelPart const*>& rSubstractionModelParts,
        const bool AddNeighbourEntities);

    /**
     * @brief Intersect given model parts to one model part.
     *
     * This method finds the intersection of all the nodes, geometries in @ref rIntersectionModelParts in to one model part
     * given by @ref rOutputSubModelPartName. The output model part is a sub-model part of
     * @ref rMainModelPart.
     *
     * The intersection is carried out by finding all the nodes, geometries in @ref rIntersectionModelParts
     * within the @ref rMainModelPart which are intersecting, and the output model part is populated with the found
     * entities in @ref rMainModelPart.
     *
     * If @ref AddNeighbourEntities is true, then all the neighbours of the merged model part
     * is found and added as well.
     *
     * Please make sure to check validity of model parts using @see CheckValidityOfModelPartsForOperations.
     *
     * @param rOutputSubModelPartName       Output sub model part name.
     * @param rMainModelPart                Main Model part.
     * @param rIntersectionModelParts       Model parts to find intersection.
     * @param AddNeighbourEntities          To add or not the neighbours.
     * @return ModelPart&                   Output sub model part.
     */
    static ModelPart& Intersect(
        const std::string& rOutputSubModelPartName,
        ModelPart& rMainModelPart,
        const std::vector<ModelPart const*>& rIntersectionModelParts,
        const bool AddNeighbourEntities);

    /**
     * @brief Checks whether given model parts list has an intersection.
     *
     * This method checks whether @ref rIntersectionModelParts has intersections.
     *
     * The intersection is carried out by finding all the nodes, geometries in @ref rIntersectionModelParts
     * within the @ref rMainModelPart which are intersecting.
     *
     * @param rMainModelPart                Main Model part.
     * @param rIntersectionModelParts       Model parts to find intersection.
     * @return true                         If there are entities intersecting.
     * @return false                        If there aren't any entities intersecting.
     */
    static bool HasIntersection(
        ModelPart& rMainModelPart,
        const std::vector<ModelPart const*>& rIntersectionModelParts);

    ///@}
};

///@}

} // namespace Kratos