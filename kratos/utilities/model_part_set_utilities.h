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

class KRATOS_API(KRATOS_CORE) ModelPartSetUtilities {
public:
    ///@name Type definitions
    ///@{

    using CNodePointersType = std::vector<ModelPart::NodeType const*>;

    template<class KeyType, class ValueType>
    using RangedKeyMapType = std::unordered_map<KeyType, ValueType, KeyHasherRange<KeyType>, KeyComparorRange<KeyType>>;

    ///@}
    ///@name Static operations
    ///@{

    static bool CheckValidityOfModelPartsForOperations(
        const ModelPart& rMainModelPart,
        const std::vector<ModelPart const*>& rCheckModelParts,
        const bool ThrowError = false);

    static ModelPart& Union(
        const std::string& rOutputSubModelPartName,
        ModelPart& rMainModelPart,
        const std::vector<ModelPart const*>& rUnionModelParts,
        const bool AddNeighbourEntities);

    static ModelPart& Substraction(
        const std::string& rOutputSubModelPartName,
        ModelPart& rMainModelPart,
        const std::vector<ModelPart const*>& rSubstractionModelParts,
        const bool AddNeighbourEntities);

    static ModelPart& Intersection(
        const std::string& rOutputSubModelPartName,
        ModelPart& rMainModelPart,
        const std::vector<ModelPart const*>& rIntersectionModelParts,
        const bool AddNeighbourEntities);

    ///@}
private:
    ///@name Private static operations
    ///@{

    ///@}
};

///@}

} // namespace Kratos