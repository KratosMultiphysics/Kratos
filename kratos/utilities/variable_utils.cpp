//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Ruben Zorrilla
//                   Vicente Mataix Ferrandiz
//

// System includes
#include <algorithm>
#include <functional>
#include <unordered_set>

// External includes

// Project includes
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

void VariableUtils::AddDofsList(
    const std::vector<std::string>& rDofsVarNamesList,
    ModelPart& rModelPart)
{
    // Create a set with the variables to be added as DOFs
    std::unordered_set<const Variable<double>*> dofs_vars_set;
    const IndexType n_dofs = rDofsVarNamesList.size();
    for (IndexType i = 0; i < n_dofs; ++i) {
        const std::string var_name = rDofsVarNamesList[i];
        KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(var_name)) << "Provided variable \'" << var_name << "\' is not in KratosComponents." << std::endl;
        const auto& r_dof_var = KratosComponents<Variable<double>>::Get(var_name);
        const Variable<double>* p_dof_var = &r_dof_var;
        dofs_vars_set.insert(p_dof_var);
        rModelPart.GetNodalSolutionStepVariablesList().AddDof(&r_dof_var);
    }

    // Add the DOFs to the model part nodes
    block_for_each(rModelPart.Nodes(), [&dofs_vars_set](Node& rNode){
        for (auto p_var : dofs_vars_set) {
            rNode.AddDof(*p_var);
        }
    });
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::AddDofsWithReactionsList(
    const std::vector<std::array<std::string,2>>& rDofsAndReactionsNamesList,
    ModelPart& rModelPart)
{
    // Create auxiliary hasher and comparors for the variable array type
    struct VariableNamesArrayHasher
    {
        std::size_t operator()(const std::array<const Variable<double>*, 2>& rDofAndReactArray) const
        {
            std::size_t seed = 0;
            HashCombine(seed, rDofAndReactArray[0]);
            HashCombine(seed, rDofAndReactArray[1]);
            return seed;
        }
    };

    struct VariableNamesArrayComparor
    {
        bool operator()(
            const std::array<const Variable<double>*, 2>& rDofAndReactArray1,
            const std::array<const Variable<double>*, 2>& rDofAndReactArray2) const
        {
            return ((rDofAndReactArray1[0] == rDofAndReactArray2[0]) && (rDofAndReactArray1[1] == rDofAndReactArray2[1]));
        }
    };

    // Create a set with the variables to be added as DOFs and reactions
    std::vector<std::string> aux_dof_var_names;
    std::vector<std::string> aux_react_var_names;
    std::array<const Variable<double>*,2> aux_dof_vars;
    const IndexType n_dofs = rDofsAndReactionsNamesList.size();
    std::unordered_set<std::array<const Variable<double>*, 2>, VariableNamesArrayHasher, VariableNamesArrayComparor> dofs_and_react_vars_set;
    for (IndexType i = 0; i < n_dofs; ++i) {
        // Get current DOF data pair
        const auto& r_dof_data = rDofsAndReactionsNamesList[i];
        const std::string var_name = r_dof_data[0];
        const std::string react_name = r_dof_data[1];

        // Check current pair data values
        KRATOS_ERROR_IF(var_name == react_name) << "DOF and reaction variable name are the same." << std::endl;
        KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(var_name)) << "Provided variable \'" << var_name << "\' is not in KratosComponents." << std::endl;
        KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(react_name)) << "Provided reaction \'" << react_name << "\' is not in KratosComponents." << std::endl;

        // Get variables from KratosComponents and insert in the auxiliary set
        // Note that using a set ensures that the insertion is done once if repeated data is provided
        const auto& r_dof_var = KratosComponents<Variable<double>>::Get(var_name);
        const auto& r_react_var = KratosComponents<Variable<double>>::Get(react_name);
        aux_dof_vars[0] = &r_dof_var;
        aux_dof_vars[1] = &r_react_var;
        auto insert_result = dofs_and_react_vars_set.insert(aux_dof_vars);

        // If the pair has been inserted, check that neither the DOF nor the reaction variables have been used before in a different pair
        if (std::get<1>(insert_result)) {
            KRATOS_ERROR_IF(std::any_of(aux_dof_var_names.begin(), aux_dof_var_names.end(), [&var_name](std::string& r_val){return r_val == var_name;}))
                << var_name << " has been already added as DOF. Use a different variable." << std::endl;
            KRATOS_ERROR_IF(std::any_of(aux_react_var_names.begin(), aux_react_var_names.end(), [&react_name](std::string& r_val){return r_val == react_name;}))
                << react_name << " has been already added as reaction to DOF. Use a different variable." << std::endl;
            aux_dof_var_names.push_back(var_name);
            aux_react_var_names.push_back(react_name);
        }

        // Add the current pair as a DOF and reaction to the nodal historical variables list
        rModelPart.GetNodalSolutionStepVariablesList().AddDof(&r_dof_var, &r_react_var);
    }

    // Add the DOFs and reactions to the model part nodes
    block_for_each(rModelPart.Nodes(), [&dofs_and_react_vars_set](Node& rNode){
        for (auto& r_dof_data : dofs_and_react_vars_set) {
            rNode.AddDof(*(r_dof_data[0]), *(r_dof_data[1]));
        }
    });
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

    block_for_each(rNodes, [&](Node& rNode){
        noalias(rNode.Coordinates()) = rNode.GetInitialPosition();
    });

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void VariableUtils::UpdateInitialToCurrentConfiguration(const ModelPart::NodesContainerType& rNodes)
{
    KRATOS_TRY;

    block_for_each(rNodes, [&](Node& rNode){
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

    block_for_each(rNodes, [&](Node& rNode){
        noalias(rNode.Coordinates()) = (rNode.GetInitialPosition()).Coordinates() + rNode.FastGetSolutionStepValue(rUpdateVariable, BufferPosition);
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

template<class TDataType>
void VariableUtils::AuxiliaryHistoricalValueSetter(
    const Variable<TDataType>& rVariable,
    const TDataType& rValue,
    NodeType& rNode)
{
    rNode.FastGetSolutionStepValue(rVariable) = rValue;
}

template<>
KRATOS_API(KRATOS_CORE) void VariableUtils::AuxiliaryHistoricalValueSetter(
    const Variable<array_1d<double,3>>& rVariable,
    const array_1d<double,3>& rValue,
    NodeType& rNode)
{
    noalias(rNode.FastGetSolutionStepValue(rVariable)) = rValue;
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

    const std::function<double(const Node&)>& r_weight_method =
        (IsInverseWeightProvided) ?
        static_cast<std::function<double(const Node&)>>([&rWeightVariable](const Node& rNode) -> double {return 1.0 / rNode.GetValue(rWeightVariable);}) :
        static_cast<std::function<double(const Node&)>>([&rWeightVariable](const Node& rNode) -> double {return rNode.GetValue(rWeightVariable);});

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

template<class TVectorType>
TVectorType VariableUtils::GetCurrentPositionsVector(
    const ModelPart::NodesContainerType& rNodes,
    const unsigned int Dimension)
{
    KRATOS_ERROR_IF(Dimension>3) << "Only Dimension<=3 is admitted by the function" << std::endl;
    TVectorType pos(rNodes.size()*Dimension);

    IndexPartition<unsigned int>(rNodes.size()).for_each(
    [&](unsigned int i){
        auto& coords = (rNodes.begin()+i)->Coordinates();
        for(unsigned int k=0; k<Dimension; k++)
            pos[i*Dimension+k] = coords[k];
        }
    );
    return pos;
}

template<class TVectorType>
TVectorType VariableUtils::GetInitialPositionsVector(
    const ModelPart::NodesContainerType& rNodes,
    const unsigned int Dimension)
{
    KRATOS_ERROR_IF(Dimension>3) << "Only Dimension<=3 is admitted by the function" << std::endl;
    TVectorType pos(rNodes.size()*Dimension);

    IndexPartition<unsigned int>(rNodes.size()).for_each(
    [&](unsigned int i){
        auto& coords = (rNodes.begin()+i)->GetInitialPosition();
        for(unsigned int k=0; k<Dimension; k++)
            pos[i*Dimension+k] = coords[k];
        }
    );
    return pos;
}

KRATOS_API(KRATOS_CORE) void VariableUtils::SetCurrentPositionsVector(
    ModelPart::NodesContainerType& rNodes,
    const Vector& rPositions)
{
    KRATOS_ERROR_IF(rPositions.size()%rNodes.size()!=0) << "Incompatible number of nodes and position data" << std::endl;

    unsigned int Dimension = rPositions.size()/rNodes.size();

    IndexPartition<unsigned int>(rNodes.size()).for_each(
    [&](unsigned int i){
        auto& coords = (rNodes.begin()+i)->Coordinates();
        for(unsigned int k=0; k<Dimension; k++)
            coords[k] = rPositions[i*Dimension+k];
        }
    );
}

KRATOS_API(KRATOS_CORE) void VariableUtils::SetInitialPositionsVector(
    ModelPart::NodesContainerType& rNodes,
    const Vector& rPositions)
{
    KRATOS_ERROR_IF(rPositions.size()%rNodes.size()!=0) << "Incompatible number of nodes and position data" << std::endl;

    unsigned int Dimension = rPositions.size()/rNodes.size();

    IndexPartition<unsigned int>(rNodes.size()).for_each(
    [&](unsigned int i){
        auto& coords = (rNodes.begin()+i)->GetInitialPosition();
        for(unsigned int k=0; k<Dimension; k++)
            coords[k] = rPositions[i*Dimension+k];
        }
    );
}

KRATOS_API(KRATOS_CORE) Vector VariableUtils::GetSolutionStepValuesVector(
                            const ModelPart::NodesContainerType& rNodes,
                            const Variable<array_1d<double,3>>& rVar,
                            const unsigned int Step,
                            const unsigned int Dimension
                            )
{
    Vector out(rNodes.size()*Dimension);

    IndexPartition<unsigned int>(rNodes.size()).for_each(
        [&](unsigned int i){
            auto& v = (rNodes.begin()+i)->FastGetSolutionStepValue(rVar,Step);
            for(unsigned int k=0; k<Dimension; k++)
                out[i*Dimension+k] = v[k];
            }
        );
    return out;
}

KRATOS_API(KRATOS_CORE) Vector VariableUtils::GetSolutionStepValuesVector(
                            const ModelPart::NodesContainerType& rNodes,
                            const Variable<double>& rVar,
                            const unsigned int Step
                            )
{
    Vector out(rNodes.size());

    IndexPartition<unsigned int>(rNodes.size()).for_each(
        [&](unsigned int i){
            out[i] = (rNodes.begin()+i)->FastGetSolutionStepValue(rVar,Step);
            }
        );
    return out;
}

KRATOS_API(KRATOS_CORE) void VariableUtils::SetSolutionStepValuesVector(
                            ModelPart::NodesContainerType& rNodes,
                            const Variable<array_1d<double,3>>& rVar,
                            const Vector& rData,
                            const unsigned int Step
                            )
{
    KRATOS_ERROR_IF(rData.size()%rNodes.size()!=0) << "Incompatible number of nodes and position data" << std::endl;

    const unsigned int Dimension = rData.size()/rNodes.size();

    IndexPartition<unsigned int>(rNodes.size()).for_each(
        [&](unsigned int i){
            auto& v = (rNodes.begin()+i)->FastGetSolutionStepValue(rVar,Step);
            for(unsigned int k=0; k<Dimension; k++)
                v[k] = rData[i*Dimension+k];
            }
        );
}

KRATOS_API(KRATOS_CORE) void VariableUtils::SetSolutionStepValuesVector(
                            ModelPart::NodesContainerType& rNodes,
                            const Variable<double>& rVar,
                            const Vector& rData,
                            const unsigned int Step
                            )
{
    KRATOS_ERROR_IF(rData.size()%rNodes.size()!=0) << "Incompatible number of nodes and position data" << std::endl;

    IndexPartition<unsigned int>(rNodes.size()).for_each(
        [&](unsigned int i){
            auto& v = (rNodes.begin()+i)->FastGetSolutionStepValue(rVar,Step);
            v = rData[i];
            }
        );
}

KRATOS_API(KRATOS_CORE) Vector VariableUtils::GetValuesVector(
    const ModelPart::NodesContainerType& rNodes,
    const Variable<array_1d<double,3>>& rVariable,
    const unsigned int Dimension
    )
{
    Vector out(rNodes.size()*Dimension);

    IndexPartition<unsigned int>(rNodes.size()).for_each(
        [&](unsigned int i){
            const auto& r_v = (rNodes.begin()+i)->GetValue(rVariable);
            for(unsigned int k=0; k<Dimension; k++)
                out[i*Dimension+k] = r_v[k];
            }
        );
    return out;
}

KRATOS_API(KRATOS_CORE) Vector VariableUtils::GetValuesVector(
    const ModelPart::NodesContainerType& rNodes,
    const Variable<double>& rVar
    )
{
    Vector out(rNodes.size());

    IndexPartition<unsigned int>(rNodes.size()).for_each(
        [&](unsigned int i){
            out[i] = (rNodes.begin()+i)->GetValue(rVar);
            }
        );
    return out;
}

KRATOS_API(KRATOS_CORE) void VariableUtils::SetValuesVector(
    ModelPart::NodesContainerType& rNodes,
    const Variable<array_1d<double,3>>& rVar,
    const Vector& rData
    )
{
    KRATOS_ERROR_IF(rData.size()%rNodes.size()!=0) << "Incompatible number of nodes and position data" << std::endl;

    const unsigned int Dimension = rData.size()/rNodes.size();

    IndexPartition<unsigned int>(rNodes.size()).for_each(
        [&](unsigned int i){
            auto& r_v = (rNodes.begin()+i)->GetValue(rVar);
            for(unsigned int k=0; k<Dimension; k++)
                r_v[k] = rData[i*Dimension+k];
            }
        );
}

KRATOS_API(KRATOS_CORE) void VariableUtils::SetValuesVector(
    ModelPart::NodesContainerType& rNodes,
    const Variable<double>& rVar,
    const Vector& rData
    )
{
    KRATOS_ERROR_IF(rData.size()%rNodes.size()!=0) << "Incompatible number of nodes and position data" << std::endl;

    IndexPartition<unsigned int>(rNodes.size()).for_each(
        [&](unsigned int i){
            auto& r_v = (rNodes.begin()+i)->GetValue(rVar);
            r_v = rData[i];
            }
        );
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

template KRATOS_API(KRATOS_CORE) Vector VariableUtils::GetInitialPositionsVector<Vector>(const ModelPart::NodesContainerType&, const unsigned int Dimension);
template KRATOS_API(KRATOS_CORE) std::vector<double> VariableUtils::GetInitialPositionsVector<std::vector<double>>(const ModelPart::NodesContainerType&, const unsigned int Dimension);

template KRATOS_API(KRATOS_CORE) Vector VariableUtils::GetCurrentPositionsVector<Vector>(const ModelPart::NodesContainerType&, const unsigned int Dimension);
template KRATOS_API(KRATOS_CORE) std::vector<double> VariableUtils::GetCurrentPositionsVector<std::vector<double>>(const ModelPart::NodesContainerType&, const unsigned int Dimension);

template KRATOS_API(KRATOS_CORE) void VariableUtils::AuxiliaryHistoricalValueSetter<int>(const Variable<int>&, const int&, NodeType&);
template KRATOS_API(KRATOS_CORE) void VariableUtils::AuxiliaryHistoricalValueSetter<double>(const Variable<double>&, const double&, NodeType&);
template KRATOS_API(KRATOS_CORE) void VariableUtils::AuxiliaryHistoricalValueSetter<Vector>(const Variable<Vector>&, const Vector&, NodeType&);
template KRATOS_API(KRATOS_CORE) void VariableUtils::AuxiliaryHistoricalValueSetter<Matrix>(const Variable<Matrix>&, const Matrix&, NodeType&);

} /* namespace Kratos.*/
