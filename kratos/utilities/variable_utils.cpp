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

} /* namespace Kratos.*/
