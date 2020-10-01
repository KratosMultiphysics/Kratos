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
        const TMap<int>& rNodalValuesMap) const
    {
        AssembleUtilities::CheckHistoricalVariable(rModelPart, rVariable);
        AssembleUtilities::AssembleCurrentDataWithNodalValuesMap(
            rModelPart, rVariable, rNodalValuesMap,
            AssembleUtilities::UpdateHistoricalNodalValue<int>);
    }

    virtual void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const
    {
        AssembleUtilities::CheckHistoricalVariable(rModelPart, rVariable);
        AssembleUtilities::AssembleCurrentDataWithNodalValuesMap(
            rModelPart, rVariable, rNodalValuesMap,
            AssembleUtilities::UpdateHistoricalNodalValue<double>);
    }

    virtual void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const
    {
        AssembleUtilities::CheckHistoricalVariable(rModelPart, rVariable);
        AssembleUtilities::AssembleCurrentDataWithNodalValuesMap(
            rModelPart, rVariable, rNodalValuesMap,
            AssembleUtilities::UpdateHistoricalNodalValue<array_1d<double, 3>>);
    }

    virtual void AssembleNonHistoricalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TMap<int>& rNodalValuesMap) const
    {
        AssembleUtilities::AssembleCurrentDataWithNodalValuesMap(
            rModelPart, rVariable, rNodalValuesMap,
            AssembleUtilities::UpdateNonHistoricalNodalValue<int>);
    }

    virtual void AssembleNonHistoricalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const
    {
        AssembleUtilities::AssembleCurrentDataWithNodalValuesMap(
            rModelPart, rVariable, rNodalValuesMap,
            AssembleUtilities::UpdateNonHistoricalNodalValue<double>);
    }

    virtual void AssembleNonHistoricalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const
    {
        AssembleUtilities::AssembleCurrentDataWithNodalValuesMap(
            rModelPart, rVariable, rNodalValuesMap,
            AssembleUtilities::UpdateNonHistoricalNodalValue<array_1d<double, 3>>);
    }

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

    template<class TDataType>
    void static UpdateNonHistoricalNodalValue(
        NodeType& rNode,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rNode.GetValue(rVariable) += rValue;
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    /**
     * @brief Assembles nodal values given in the map
     *
     * This method can assemble nodal values according to given nodal map.
     * No clearing of nodal values are done, therefore, assemble will add
     * values to existing nodal values.
     *
     * @tparam TDataType
     * @tparam TUpdateFunction
     * @param rModelPart            Model part where nodes need to update
     * @param rVariable             Variable to store assembled values
     * @param rNodalValuesMap       Nodal values map with node_id and value
     * @param rUpdateFunction       Update function
     */
    template<class TDataType, class TUpdateFunction>
    void static AssembleCurrentDataWithNodalValuesMap(
        ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const TMap<TDataType>& rNodalValuesMap,
        const TUpdateFunction&& rUpdateFunction)
    {
        KRATOS_TRY

        std::vector<int> node_ids;
        node_ids.reserve(rNodalValuesMap.size());
        for (const auto& r_item : rNodalValuesMap) {
            node_ids.push_back(r_item.first);
        }

        IndexPartition<int>(node_ids.size()).for_each([&](const int Index) {
            const int node_id = node_ids[Index];
            auto& r_node = rModelPart.GetNode(node_id);
            rUpdateFunction(r_node, rVariable, rNodalValuesMap.find(node_id)->second);
        });

        KRATOS_CATCH("");
    }

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_ASSEMBLE_UTILITIES