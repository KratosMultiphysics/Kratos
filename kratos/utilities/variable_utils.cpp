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
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"

/* External includes */

/* Project includes */

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
        NodesContainerType::iterator i = rNodes.begin() + k;
        noalias(i->FastGetSolutionStepValue(rVariable)) = Value;
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
    DoubleVarType& SavedVariable,
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

void VariableUtils::CopyVectorVar(
    const ArrayVarType& OriginVariable,
    ArrayVarType& DestinationVariable,
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

void VariableUtils::CopyScalarVar(
    const DoubleVarType& OriginVariable,
    DoubleVarType& DestinationVariable,
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

void VariableUtils::SetToZero_VectorVar(
    const ArrayVarType& Variable, 
    NodesContainerType& rNodes
    )
{
    KRATOS_TRY
    
    #pragma omp parallel for
    for (int k = 0; k < static_cast<int> (rNodes.size()); ++k) {
        NodesContainerType::iterator i = rNodes.begin() + k;
        noalias(i->FastGetSolutionStepValue(Variable)) = ZeroVector(3);
    }
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::SetToZero_ScalarVar(
    const DoubleVarType& Variable, 
    NodesContainerType& rNodes
    )
{
    KRATOS_TRY
    
    #pragma omp parallel for
    for (int k = 0; k < static_cast<int> (rNodes.size()); ++k) {
        NodesContainerType::iterator i = rNodes.begin() + k;
        i->FastGetSolutionStepValue(Variable) = 0.0;
    }
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart::NodesContainerType VariableUtils::SelectNodeList(
    const DoubleVarType& Variable,
    const double& Value,
    NodesContainerType& rOriginNodes
    )
{
    KRATOS_TRY
    
    NodesContainerType selected_nodes;
    for (NodesContainerType::iterator it_node = rOriginNodes.begin(); it_node != rOriginNodes.end(); ++it_node) {
        if (it_node->FastGetSolutionStepValue(Variable) == Value)
            selected_nodes.push_back(*(it_node.base()));
    }

    return selected_nodes;
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> VariableUtils::SumNonHistoricalNodeVectorVariable(
    const Variable<array_1d<double, 3> >& rVar,
    ModelPart& rModelPart
    )
{
    KRATOS_TRY

    array_1d<double, 3> sum_value = ZeroVector(3);

    #pragma omp parallel
    {
        array_1d<double, 3> private_sum_value = ZeroVector(3);

        #pragma omp for
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); ++k) {
            NodesContainerType::iterator it_node = rModelPart.NodesBegin() + k;
            private_sum_value += it_node->GetValue(rVar);
        }

        for (int j = 0; j < static_cast<int>(sum_value.size()); ++j) {
            #pragma omp atomic
            sum_value[j] += private_sum_value[j];
        }
    }

    rModelPart.GetCommunicator().SumAll(sum_value);

    return sum_value;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> VariableUtils::SumHistoricalNodeVectorVariable(
    const Variable<array_1d<double, 3> >& rVar,
    ModelPart& rModelPart,
    const unsigned int& rBuffStep
    )
{
    KRATOS_TRY

    array_1d<double, 3> sum_value = ZeroVector(3);

    #pragma omp parallel
    {
        array_1d<double, 3> private_sum_value = ZeroVector(3);

        #pragma omp for
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); ++k) {
            NodesContainerType::iterator it_node = rModelPart.NodesBegin() + k;
            private_sum_value += it_node->GetSolutionStepValue(rVar, rBuffStep);
        }

        for (int j = 0; j < static_cast<int>(sum_value.size()); ++j) {
            #pragma omp atomic
            sum_value[j] += private_sum_value[j];
        }
    }

    rModelPart.GetCommunicator().SumAll(sum_value);

    return sum_value;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> VariableUtils::SumConditionVectorVariable(
    const Variable<array_1d<double, 3> >& rVar,
    ModelPart& rModelPart
    )
{
    KRATOS_TRY

    array_1d<double, 3> sum_value = ZeroVector(3);

    #pragma omp parallel
    {
        array_1d<double, 3> private_sum_value = ZeroVector(3);

        #pragma omp for
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfConditions()); ++k) {
            ConditionsContainerType::iterator it_cond = rModelPart.ConditionsBegin() + k;
            private_sum_value += it_cond->GetValue(rVar);
        }

        for (int j = 0; j < static_cast<int>(sum_value.size()); ++j) {
            #pragma omp atomic
            sum_value[j] += private_sum_value[j];
        }
    }

    rModelPart.GetCommunicator().SumAll(sum_value);

    return sum_value;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> VariableUtils::SumElementVectorVariable(
    const Variable<array_1d<double, 3> >& rVar,
    ModelPart& rModelPart
    )
{
    KRATOS_TRY

    array_1d<double, 3> sum_value = ZeroVector(3);

    #pragma omp parallel
    {
        array_1d<double, 3> private_sum_value = ZeroVector(3);

        #pragma omp for
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfElements()); ++k) {
            ElementsContainerType::iterator it_elem = rModelPart.ElementsBegin() + k;
            private_sum_value += it_elem->GetValue(rVar);
        }

        for (int j = 0; j < static_cast<int>(sum_value.size()); ++j) {
            #pragma omp atomic
            sum_value[j] += private_sum_value[j];
        }
    }

    rModelPart.GetCommunicator().SumAll(sum_value);

    return sum_value;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

bool VariableUtils::CheckVariableKeys()
{
    KRATOS_TRY

    CheckVariableKeysHelper< Variable<double> >();
    CheckVariableKeysHelper< Variable<array_1d<double,3> > >();
    CheckVariableKeysHelper< Variable<bool> >();
    CheckVariableKeysHelper< Variable<int> >();
    CheckVariableKeysHelper< Variable<unsigned int> >();
    CheckVariableKeysHelper< Variable<Vector> >();
    CheckVariableKeysHelper< Variable<Matrix> >();
    CheckVariableKeysHelper< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >();
        
    return true;
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

bool VariableUtils::CheckDofs(ModelPart& rModelPart)
{
    KRATOS_TRY

    for(auto& node : rModelPart.Nodes()) {
        for (auto& dof : node.GetDofs()) {
//                 KRATOS_ERROR_IF_NOT(node.SolutionStepsDataHas(dof.GetVariable())) << "Node : " << node << " does not have allocated space for the variable " << dof << std::endl;
            KRATOS_ERROR_IF(dof.GetVariable().Key() == 0) << "Found a zero key on a dof of node " << node << std::endl;

        }
    }
    
    return true;
    
    KRATOS_CATCH("")
}
} /* namespace Kratos.*/
