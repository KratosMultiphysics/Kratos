//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Ruben Zorrilla
//                   Vicente Mataix Ferrandiz
//
//

/* System includes */
#include <functional>

/* External includes */

/* Project includes */
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

ModelPart::NodesContainerType VariableUtils::SelectNodeList(
    const DoubleVarType& Variable,
    const double Value,
    const NodesContainerType& rOriginNodes
    )
{
    KRATOS_TRY

    NodesContainerType selected_nodes;
    for (auto it_node = rOriginNodes.begin(); it_node != rOriginNodes.end(); ++it_node) {
        if (std::abs(it_node->FastGetSolutionStepValue(Variable) - Value) <
            std::numeric_limits<double>::epsilon()) {
            selected_nodes.push_back(*(it_node.base()));
        }
    }

    return selected_nodes;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> VariableUtils::SumNonHistoricalNodeVectorVariable(
    const ArrayVarType& rVar,
    const ModelPart& rModelPart
    )
{
    KRATOS_TRY

    array_1d<double, 3> sum_value = ZeroVector(3);
    auto& r_comm = rModelPart.GetCommunicator();

    sum_value = block_for_each<SumReduction<array_1d<double,3>>>(r_comm.LocalMesh().Nodes(),[&](NodeType& rNode){
        return rNode.GetValue(rVar);
    });

    return r_comm.GetDataCommunicator().SumAll(sum_value);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> VariableUtils::SumConditionVectorVariable(
    const ArrayVarType& rVar,
    const ModelPart& rModelPart
    )
{
    KRATOS_TRY

    array_1d<double, 3> sum_value = ZeroVector(3);
    auto& r_comm = rModelPart.GetCommunicator();

    sum_value = block_for_each<SumReduction<array_1d<double,3>>>(r_comm.LocalMesh().Conditions(),[&](ConditionType& rCond){
        return rCond.GetValue(rVar);
    });

    return r_comm.GetDataCommunicator().SumAll(sum_value);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> VariableUtils::SumElementVectorVariable(
    const ArrayVarType& rVar,
    const ModelPart& rModelPart
    )
{
    KRATOS_TRY

    array_1d<double, 3> sum_value = ZeroVector(3);
    auto& r_comm = rModelPart.GetCommunicator();

    sum_value = block_for_each<SumReduction<array_1d<double,3>>>(r_comm.LocalMesh().Elements(),[&](ElementType& rElem){
        return rElem.GetValue(rVar);
    });

    return r_comm.GetDataCommunicator().SumAll(sum_value);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

bool VariableUtils::CheckVariableKeys()
{
    KRATOS_TRY

    CheckVariableKeysHelper< Variable<double> >();
    CheckVariableKeysHelper< Variable<array_1d<double,3> > >();
    CheckVariableKeysHelper< Variable<array_1d<double,4> > >();
    CheckVariableKeysHelper< Variable<array_1d<double,6> > >();
    CheckVariableKeysHelper< Variable<array_1d<double,9> > >();
    CheckVariableKeysHelper< Variable<bool> >();
    CheckVariableKeysHelper< Variable<int> >();
    CheckVariableKeysHelper< Variable<unsigned int> >();
    CheckVariableKeysHelper< Variable<Vector> >();
    CheckVariableKeysHelper< Variable<Matrix> >();

    return true;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::UpdateCurrentToInitialConfiguration(const ModelPart::NodesContainerType& rNodes)
{
    KRATOS_TRY;

    block_for_each(rNodes, [&](Node<3>& rNode){
        noalias(rNode.Coordinates()) = rNode.GetInitialPosition();
    });

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::UpdateInitialToCurrentConfiguration(const ModelPart::NodesContainerType& rNodes)
{
    KRATOS_TRY;

    block_for_each(rNodes, [&](Node<3>& rNode){
        noalias(rNode.GetInitialPosition().Coordinates()) = rNode.Coordinates();
    });

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::UpdateCurrentPosition(
    const ModelPart::NodesContainerType& rNodes,
    const ArrayVarType& rUpdateVariable,
    const IndexType BufferPosition
    )
{
    KRATOS_TRY;

    block_for_each(rNodes, [&](Node<3>& rNode){
        const auto& r_update_coords = rNode.FastGetSolutionStepValue(rUpdateVariable, BufferPosition);
        noalias(rNode.Coordinates()) = (rNode.GetInitialPosition()).Coordinates() + r_update_coords;
    });

    KRATOS_CATCH("");
}

template <>
KRATOS_API(KRATOS_CORE) ModelPart::NodesContainerType& VariableUtils::GetContainer<ModelPart::NodesContainerType>(ModelPart& rModelPart)
{
    return rModelPart.Nodes();
}


template <>
KRATOS_API(KRATOS_CORE) ModelPart::ElementsContainerType& VariableUtils::GetContainer<ModelPart::ElementsContainerType>(ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template <>
KRATOS_API(KRATOS_CORE) ModelPart::ConditionsContainerType& VariableUtils::GetContainer<ModelPart::ConditionsContainerType>(ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

template <>
KRATOS_API(KRATOS_CORE) const ModelPart::NodesContainerType& VariableUtils::GetContainer<ModelPart::NodesContainerType>(const ModelPart& rModelPart)
{
    return rModelPart.Nodes();
}

template <>
KRATOS_API(KRATOS_CORE) const ModelPart::ElementsContainerType& VariableUtils::GetContainer<ModelPart::ElementsContainerType>(const ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template <>
KRATOS_API(KRATOS_CORE) const ModelPart::ConditionsContainerType& VariableUtils::GetContainer<ModelPart::ConditionsContainerType>(const ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

template <class TDataType, class TContainerType, class TWeightDataType>
void VariableUtils::WeightedAccumulateVariableOnNodes(
    ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const Variable<TWeightDataType>& rWeightVariable,
    const bool IsInverseWeightProvided)
{
    KRATOS_TRY

    SetNonHistoricalVariableToZero(rVariable, rModelPart.Nodes());

    auto& r_entities = GetContainer<TContainerType>(rModelPart);
    const int n_entities = r_entities.size();

    const std::function<double(const Node<3>&)>& r_weight_method =
        (IsInverseWeightProvided) ?
        static_cast<std::function<double(const Node<3>&)>>([rWeightVariable](const Node<3>& rNode) -> double {return 1.0 / rNode.GetValue(rWeightVariable);}) :
        static_cast<std::function<double(const Node<3>&)>>([rWeightVariable](const Node<3>& rNode) -> double {return rNode.GetValue(rWeightVariable);});

#pragma omp parallel for
    for (int i_entity = 0; i_entity < n_entities; ++i_entity)
    {
        auto it_entity = r_entities.begin() + i_entity;
        auto& r_geometry = it_entity->GetGeometry();

        const auto& r_value = it_entity->GetValue(rVariable);
        for (int i_node = 0; i_node < static_cast<int>(r_geometry.PointsNumber()); ++i_node)
        {
            auto& r_node = r_geometry[i_node];

            KRATOS_DEBUG_ERROR_IF(!r_node.Has(rWeightVariable))
                << "Non-historical nodal " << rWeightVariable.Name() << " at "
                << r_node << " is not initialized in " << rModelPart.Name()
                << ". Please initialize it first.";

            const double weight = r_weight_method(r_node);

            r_node.SetLock();
            r_node.GetValue(rVariable) += r_value * weight;
            r_node.UnSetLock();
        }
    }

    rModelPart.GetCommunicator().AssembleNonHistoricalData(rVariable);

    KRATOS_CATCH("");
}

// template instantiations
template KRATOS_API(KRATOS_CORE) void VariableUtils::WeightedAccumulateVariableOnNodes<double, ModelPart::ConditionsContainerType, int>(
    ModelPart&, const Variable<double>&, const Variable<int>&, const bool);
template KRATOS_API(KRATOS_CORE) void VariableUtils::WeightedAccumulateVariableOnNodes<array_1d<double, 3>, ModelPart::ConditionsContainerType, int>(
    ModelPart&, const Variable<array_1d<double, 3>>&, const Variable<int>&, const bool);

template KRATOS_API(KRATOS_CORE) void VariableUtils::WeightedAccumulateVariableOnNodes<double, ModelPart::ElementsContainerType, int>(
    ModelPart&, const Variable<double>&, const Variable<int>&, const bool);
template KRATOS_API(KRATOS_CORE) void VariableUtils::WeightedAccumulateVariableOnNodes<array_1d<double, 3>, ModelPart::ElementsContainerType, int>(
    ModelPart&, const Variable<array_1d<double, 3>>&, const Variable<int>&, const bool);

template KRATOS_API(KRATOS_CORE) void VariableUtils::WeightedAccumulateVariableOnNodes<double, ModelPart::ConditionsContainerType, double>(
    ModelPart&, const Variable<double>&, const Variable<double>&, const bool);
template KRATOS_API(KRATOS_CORE) void VariableUtils::WeightedAccumulateVariableOnNodes<array_1d<double, 3>, ModelPart::ConditionsContainerType, double>(
    ModelPart&, const Variable<array_1d<double, 3>>&, const Variable<double>&, const bool);

template KRATOS_API(KRATOS_CORE) void VariableUtils::WeightedAccumulateVariableOnNodes<double, ModelPart::ElementsContainerType, double>(
    ModelPart&, const Variable<double>&, const Variable<double>&, const bool);
template KRATOS_API(KRATOS_CORE) void VariableUtils::WeightedAccumulateVariableOnNodes<array_1d<double, 3>, ModelPart::ElementsContainerType, double>(
    ModelPart&, const Variable<array_1d<double, 3>>&, const Variable<double>&, const bool);

} /* namespace Kratos.*/
