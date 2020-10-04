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
#include "containers/global_pointers_unordered_map.h"
#include "containers/variable.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/global_pointer_utilities.h"
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

    template<class TEntityType, class TDataType>
    using TGPMap = GlobalPointersUnorderedMap<TEntityType, TDataType>;

    using NodeType = ModelPart::NodeType;

    using NodesContainerType = ModelPart::NodesContainerType;

    using ElementType = ModelPart::ElementType;

    using ElementsContainerType = ModelPart::ElementsContainerType;

    using ConditionType = ModelPart::ConditionType;

    using ConditionsContainerType = ModelPart::ConditionsContainerType;

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
        const TMap<int>& rNodalValuesMap) const
    {
        this->AssembleCurrentDataWithValuesMap(
            rModelPart, rVariable,
            AssembleUtilities::GetGlobalPointersMap<NodesContainerType, int>(
                rModelPart, rNodalValuesMap));
    }

    virtual void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const
    {
        this->AssembleCurrentDataWithValuesMap(
            rModelPart, rVariable,
            AssembleUtilities::GetGlobalPointersMap<NodesContainerType, double>(
                rModelPart, rNodalValuesMap));
    }

    virtual void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const
    {
        this->AssembleCurrentDataWithValuesMap(
            rModelPart, rVariable,
            AssembleUtilities::GetGlobalPointersMap<NodesContainerType, array_1d<double, 3>>(
                rModelPart, rNodalValuesMap));
    }

    virtual void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TGPMap<NodeType, int>& rNodalValuesMap) const;

    virtual void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TGPMap<NodeType, double>& rNodalValuesMap) const;

    virtual void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TGPMap<NodeType, array_1d<double, 3>>& rNodalValuesMap) const;

    virtual void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TMap<int>& rNodalValuesMap) const
    {
        this->AssembleNonHistoricalNodalDataWithValuesMap(
            rModelPart, rVariable,
            AssembleUtilities::GetGlobalPointersMap<NodesContainerType, int>(
                rModelPart, rNodalValuesMap));
    }

    virtual void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const
    {
        this->AssembleNonHistoricalNodalDataWithValuesMap(
            rModelPart, rVariable,
            AssembleUtilities::GetGlobalPointersMap<NodesContainerType, double>(
                rModelPart, rNodalValuesMap));
    }

    virtual void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const
    {
        this->AssembleNonHistoricalNodalDataWithValuesMap(
            rModelPart, rVariable,
            AssembleUtilities::GetGlobalPointersMap<NodesContainerType, array_1d<double, 3>>(
                rModelPart, rNodalValuesMap));
    }

    virtual void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TGPMap<NodeType, int>& rNodalValuesMap) const;

    virtual void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TGPMap<NodeType, double>& rNodalValuesMap) const;

    virtual void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TGPMap<NodeType, array_1d<double, 3>>& rNodalValuesMap) const;

    virtual void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TMap<int>& rElementValuesMap) const
    {
        this->AssembleElementDataWithValuesMap(
            rModelPart, rVariable,
            AssembleUtilities::GetGlobalPointersMap<ElementsContainerType, int>(
                rModelPart, rElementValuesMap));
    }

    virtual void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rElementValuesMap) const
    {
        this->AssembleElementDataWithValuesMap(
            rModelPart, rVariable,
            AssembleUtilities::GetGlobalPointersMap<ElementsContainerType, double>(
                rModelPart, rElementValuesMap));
    }

    virtual void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rElementValuesMap) const
    {
        this->AssembleElementDataWithValuesMap(
            rModelPart, rVariable,
            AssembleUtilities::GetGlobalPointersMap<ElementsContainerType, array_1d<double, 3>>(
                rModelPart, rElementValuesMap));
    }

    virtual void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TGPMap<ElementType, int>& rNodalValuesMap) const;

    virtual void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TGPMap<ElementType, double>& rNodalValuesMap) const;

    virtual void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TGPMap<ElementType, array_1d<double, 3>>& rNodalValuesMap) const;

    virtual void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TMap<int>& rConditionValuesMap) const
    {
        this->AssembleConditionDataWithValuesMap(
            rModelPart, rVariable,
            AssembleUtilities::GetGlobalPointersMap<ConditionsContainerType, int>(
                rModelPart, rConditionValuesMap));
    }

    virtual void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rConditionValuesMap) const
    {
        this->AssembleConditionDataWithValuesMap(
            rModelPart, rVariable,
            AssembleUtilities::GetGlobalPointersMap<ConditionsContainerType, double>(
                rModelPart, rConditionValuesMap));
    }

    virtual void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rConditionValuesMap) const
    {
        this->AssembleConditionDataWithValuesMap(
            rModelPart, rVariable,
            AssembleUtilities::GetGlobalPointersMap<ConditionsContainerType, array_1d<double, 3>>(
                rModelPart, rConditionValuesMap));
    }

    virtual void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TGPMap<ConditionType, int>& rNodalValuesMap) const;

    virtual void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TGPMap<ConditionType, double>& rNodalValuesMap) const;

    virtual void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TGPMap<ConditionType, array_1d<double, 3>>& rNodalValuesMap) const;

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    template<class TContainerType>
    static TContainerType& GetContainer(ModelPart& rModelPart);

    template <class TContainerType, class TDataType>
    static TGPMap<typename TContainerType::value_type, TDataType> GetGlobalPointersMap(
        ModelPart& rModelPart,
        const TMap<TDataType>& rValuesMap)
    {
        KRATOS_TRY

        const auto& entity_id_list = GetKeys(rValuesMap);
        const auto& id_gp_map = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(
            GetContainer<TContainerType>(rModelPart), entity_id_list,
            rModelPart.GetCommunicator().GetDataCommunicator());

        TGPMap<typename TContainerType::value_type, TDataType> gp_map;
        for (const auto& r_item : id_gp_map) {
            gp_map[r_item.second] = rValuesMap.find(r_item.first)->second;
        }

        return gp_map;

        KRATOS_CATCH("");
    }

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
    static void UpdateHistoricalNodalValue(
        NodeType& rNode,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rNode.FastGetSolutionStepValue(rVariable) += rValue;
    }

    template<class TEntityType, class TDataType>
    static void UpdateNonHistoricalValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rEntity.GetValue(rVariable) += rValue;
    }

    template<class TKey, class TValue, class T1, class T2>
    static std::vector<TKey> GetKeys(const std::unordered_map<TKey, TValue, T1, T2>& rMap) {
        std::vector<TKey> key_list;
        key_list.reserve(rMap.size());
        for (const auto& r_item : rMap) {
            key_list.push_back(r_item.first);
        }

        return key_list;
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    /**
     * @brief Assembles entity values given in the map
     *
     * This method can assemble entity(nodal/elemental/condition) values according to given values map.
     * No clearing of entity values are done, therefore, assemble will add
     * values to existing values.
     *
     * @tparam TContainerType
     * @tparam TDataType
     * @tparam TUpdateFunction
     * @param rContainer            Entity container
     * @param rVariable             Variable to store assembled values
     * @param rValuesMap            Values map with entity global pointer and value
     * @param rUpdateFunction       Update function
     */
    template<class TContainerType, class TDataType, class TUpdateFunction>
    static void AssembleDataWithEntityValuesMap(
        TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const TGPMap<typename TContainerType::value_type, TDataType>& rValuesMap,
        const TUpdateFunction&& rUpdateFunction)
    {
        KRATOS_TRY

        const auto& entity_gps = GetKeys(rValuesMap);

        IndexPartition<int>(entity_gps.size()).for_each([&](const int Index) {
            auto entity_gp = entity_gps[Index];
            rUpdateFunction(*entity_gp, rVariable, rValuesMap.find(entity_gp)->second);
        });

        KRATOS_CATCH("");
    }

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_ASSEMBLE_UTILITIES