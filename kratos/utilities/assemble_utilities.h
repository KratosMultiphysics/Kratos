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

#if !defined(KRATOS_ASSEMBLE_UTILITIES)
#define  KRATOS_ASSEMBLE_UTILITIES

// System includes
#include <unordered_map>
#include <vector>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "containers/variable.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) AssembleUtilities
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AssembleUtilities);

    template<class TDataType>
    using TMap = std::unordered_map<int, TDataType>;
    using NodeType = ModelPart::NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    virtual ~AssembleUtilities() = default;

    ///@}
    ///@name Operations
    ///@{

    virtual void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TMap<int>& rNodalValuesMap) const;

    virtual void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const;

    virtual void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const;

    virtual void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TMap<int>& rNodalValuesMap) const;

    virtual void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const;

    virtual void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const;

    virtual void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TMap<int>& rNodalValuesMap) const;

    virtual void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const;

    virtual void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const;

    virtual void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TMap<int>& rNodalValuesMap) const;

    virtual void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const;

    virtual void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const;

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    template<class TDataType>
    static void CheckHistoricalVariable(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!rModelPart.HasNodalSolutionStepVariable(rVariable))
            << rVariable.Name() << " not found in solution step variables list of "
            << rModelPart.Name() << ".\n";

        KRATOS_CATCH("");
    }

    template<class TDataType>
    void static UpdateHistoricalNodalValue(
        NodeType& rNode,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rNode.FastGetSolutionStepValue(rVariable) += rValue;
    }

    template<class TDataType, class TEntityType>
    void static UpdateNonHistoricalValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rEntity.GetValue(rVariable) += rValue;
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    template<class TContainerType>
    static TContainerType& GetContainer(ModelPart& rModelPart);

    /**
     * @brief Assembles entity values given in the map
     *
     * This method can assemble entity(nodal/elemental/condition) values according to given values map.
     * No clearing of entity values are done, therefore, assemble will add
     * values to existing values.
     *
     * @tparam TDataType
     * @tparam TUpdateFunction
     * @param rModelPart            Model part where nodes need to update
     * @param rVariable             Variable to store assembled values
     * @param rValuesMap            Values map with entity_id and value
     * @param rUpdateFunction       Update function
     */
    template<class TDataType, class TContainerType, class TUpdateFunction>
    void static AssembleDataWithEntityValuesMap(
        TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const TMap<TDataType>& rValuesMap,
        const TUpdateFunction&& rUpdateFunction)
    {
        KRATOS_TRY

        std::vector<int> entity_ids;
        entity_ids.reserve(rValuesMap.size());
        for (const auto& r_item : rValuesMap) {
            entity_ids.push_back(r_item.first);
        }

        IndexPartition<int>(entity_ids.size()).for_each([&](const int Index) {
            const int entity_id = entity_ids[Index];
            auto p_entity = rContainer.find(entity_id);

            KRATOS_ERROR_IF(p_entity == rContainer.end())
                << "Entity with id " << entity_id << " not found.\n";

            rUpdateFunction(*p_entity, rVariable, rValuesMap.find(entity_id)->second);
        });

        KRATOS_CATCH("");
    }

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_ASSEMBLE_UTILITIES