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

void VariableUtils::SetVectorVar(
    const ArrayVarType& rVariable,
    const array_1d<double, 3 >& Value,
    NodesContainerType& rNodes
    )
{
    KRATOS_TRY

    #pragma omp parallel for
    for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
        NodesContainerType::iterator it_node = rNodes.begin() + k;
        noalias(it_node->FastGetSolutionStepValue(rVariable)) = Value;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::SetVectorVarForFlag(
    const ArrayVarType& rVariable,
    const array_1d<double, 3 >& Value,
    NodesContainerType& rNodes,
    const Flags Flag,
    const bool Check
    )
{
    KRATOS_TRY

    #pragma omp parallel for
    for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
        NodesContainerType::iterator it_node = rNodes.begin() + k;
        if (it_node->Is(Flag) == Check) noalias(it_node->FastGetSolutionStepValue(rVariable)) = Value;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::SetNonHistoricalVectorVar(
    const ArrayVarType& rVariable,
    const array_1d<double, 3 >& Value,
    NodesContainerType& rNodes
    )
{
    KRATOS_TRY

    #pragma omp parallel for
    for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
        NodesContainerType::iterator it_node = rNodes.begin() + k;
        it_node->SetValue(rVariable, Value);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::SaveVectorVar(
    const ArrayVarType& OriginVariable,
    const ArrayVarType& SavedVariable,
    NodesContainerType& rNodes
    )
{
    KRATOS_TRY

    #pragma omp parallel for
    for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
        NodesContainerType::iterator it_node = rNodes.begin() + k;
        it_node->SetValue(SavedVariable, it_node->FastGetSolutionStepValue(OriginVariable));
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::SaveScalarVar(
    const DoubleVarType& OriginVariable,
    const DoubleVarType& SavedVariable,
    NodesContainerType& rNodes
    )
{
    KRATOS_TRY

    #pragma omp parallel for
    for (int k = 0; k < static_cast<int> (rNodes.size()); ++k) {
        NodesContainerType::iterator it_node = rNodes.begin() + k;
        it_node->SetValue(SavedVariable,it_node->FastGetSolutionStepValue(OriginVariable));
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::SaveVectorNonHistoricalVar(
    const ArrayVarType& OriginVariable,
    const ArrayVarType& SavedVariable,
    NodesContainerType& rNodes
    )
{
    KRATOS_TRY

    #pragma omp parallel for
    for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
        NodesContainerType::iterator it_node = rNodes.begin() + k;
        it_node->SetValue(SavedVariable, it_node->GetValue(OriginVariable));
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::SaveScalarNonHistoricalVar(
    const DoubleVarType& OriginVariable,
    const DoubleVarType& SavedVariable,
    NodesContainerType& rNodes
    )
{
    KRATOS_TRY

    #pragma omp parallel for
    for (int k = 0; k < static_cast<int> (rNodes.size()); ++k) {
        NodesContainerType::iterator it_node = rNodes.begin() + k;
        it_node->SetValue(SavedVariable,it_node->GetValue(OriginVariable));
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::CopyVectorVar(
    const ArrayVarType& OriginVariable,
    const ArrayVarType& DestinationVariable,
    NodesContainerType& rNodes
    )
{
    KRATOS_TRY

    #pragma omp parallel for
    for (int k = 0; k < static_cast<int> (rNodes.size()); ++k) {
        NodesContainerType::iterator it_node = rNodes.begin() + k;
        noalias(it_node->FastGetSolutionStepValue(DestinationVariable)) = it_node->FastGetSolutionStepValue(OriginVariable);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::CopyComponentVar(
    const ComponentVarType& OriginVariable,
    const ComponentVarType& DestinationVariable,
    NodesContainerType& rNodes
    )
{
    KRATOS_TRY

    #pragma omp parallel for
    for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
        NodesContainerType::iterator it_node = rNodes.begin() + k;
        it_node->FastGetSolutionStepValue(DestinationVariable) = it_node->FastGetSolutionStepValue(OriginVariable);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::CopyScalarVar(
    const DoubleVarType& OriginVariable,
    const DoubleVarType& DestinationVariable,
    NodesContainerType& rNodes
    )
{
    KRATOS_TRY

    #pragma omp parallel for
    for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
        NodesContainerType::iterator it_node = rNodes.begin() + k;
        it_node->FastGetSolutionStepValue(DestinationVariable) = it_node->FastGetSolutionStepValue(OriginVariable);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

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

    #pragma omp parallel
    {
        array_1d<double, 3> private_sum_value = ZeroVector(3);

        #pragma omp for
        for (int k = 0; k < static_cast<int>(r_comm.LocalMesh().NumberOfNodes()); ++k) {
            const auto it_node = r_comm.LocalMesh().NodesBegin() + k;
            private_sum_value += it_node->GetValue(rVar);
        }

        for (int j = 0; j < static_cast<int>(sum_value.size()); ++j) {
            #pragma omp atomic
            sum_value[j] += private_sum_value[j];
        }
    }

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
    const auto& r_comm = rModelPart.GetCommunicator();

    #pragma omp parallel
    {
        array_1d<double, 3> private_sum_value = ZeroVector(3);

        #pragma omp for
        for (int k = 0; k < static_cast<int>(r_comm.LocalMesh().NumberOfConditions()); ++k) {
            const auto it_cond = r_comm.LocalMesh().ConditionsBegin() + k;
            private_sum_value += it_cond->GetValue(rVar);
        }

        for (int j = 0; j < static_cast<int>(sum_value.size()); ++j) {
            #pragma omp atomic
            sum_value[j] += private_sum_value[j];
        }
    }

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

    #pragma omp parallel
    {
        array_1d<double, 3> private_sum_value = ZeroVector(3);

        #pragma omp for
        for (int k = 0; k < static_cast<int>(r_comm.LocalMesh().NumberOfElements()); ++k) {
            const auto it_elem = r_comm.LocalMesh().ElementsBegin() + k;
            private_sum_value += it_elem->GetValue(rVar);
        }

        for (int j = 0; j < static_cast<int>(sum_value.size()); ++j) {
            #pragma omp atomic
            sum_value[j] += private_sum_value[j];
        }
    }

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
    CheckVariableKeysHelper< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >();
    CheckVariableKeysHelper< VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >();
    CheckVariableKeysHelper< VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >();
    CheckVariableKeysHelper< VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >();

    return true;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

bool VariableUtils::CheckDofs(ModelPart& rModelPart)
{
    KRATOS_TRY

    for(auto& node : rModelPart.Nodes()) {
        for (Kratos::unique_ptr<Dof<double>>& p_dof : node.GetDofs()) {
//                 KRATOS_ERROR_IF_NOT(node.SolutionStepsDataHas(dof.GetVariable())) << "Node : " << node << " does not have allocated space for the variable " << dof << std::endl;
            KRATOS_CHECK_VARIABLE_KEY(p_dof->GetVariable());

        }
    }

    return true;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::UpdateCurrentToInitialConfiguration(const ModelPart::NodesContainerType& rNodes) {
    KRATOS_TRY;

    const int num_nodes = static_cast<int>(rNodes.size());
    const auto it_node_begin = rNodes.begin();

    #pragma omp parallel for
    for (int i=0; i<num_nodes; i++) {
        const auto it_node  = it_node_begin + i;
        noalias(it_node->Coordinates()) = it_node->GetInitialPosition();
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::UpdateInitialToCurrentConfiguration(const ModelPart::NodesContainerType& rNodes) {
    KRATOS_TRY;

    const int num_nodes = static_cast<int>(rNodes.size());
    const auto it_node_begin = rNodes.begin();

    #pragma omp parallel for
    for (int i=0; i<num_nodes; i++) {
        const auto it_node  = it_node_begin + i;
        noalias(it_node->GetInitialPosition().Coordinates()) = it_node->Coordinates();
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::UpdateCurrentPosition(
    const ModelPart::NodesContainerType &rNodes,
    const ArrayVarType &rUpdateVariable,
    const IndexType BufferPosition
    )
{
    KRATOS_TRY;

    const int num_nodes = static_cast<int>(rNodes.size());
    const auto it_node_begin = rNodes.begin();

    #pragma omp parallel for
    for (int i_node = 0; i_node < num_nodes; ++i_node) {
        const auto it_node  = it_node_begin + i_node;
        const auto& r_update_coords = it_node->FastGetSolutionStepValue(rUpdateVariable, BufferPosition);
        noalias(it_node->Coordinates()) = (it_node->GetInitialPosition()).Coordinates() + r_update_coords;
    }

    KRATOS_CATCH("");
}

void VariableUtils::AuxiliaryInitializeValue(double &rValue)
{
    rValue = 0.0;
}

void VariableUtils::AuxiliaryInitializeValue(array_1d<double, 3> &rValue)
{
    rValue = ZeroVector(3);
}

void VariableUtils::AuxiliaryAtomicAdd(
    const double &rPrivateValue,
    double &rSumValue)
{
#pragma omp atomic
        rSumValue += rPrivateValue;
}

void VariableUtils::AuxiliaryAtomicAdd(
    const array_1d<double, 3> &rPrivateValue,
    array_1d<double, 3> &rSumValue)
{
#pragma omp atomic
        rSumValue[0] += rPrivateValue[0];
#pragma omp atomic
        rSumValue[1] += rPrivateValue[1];
#pragma omp atomic
        rSumValue[2] += rPrivateValue[2];
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
