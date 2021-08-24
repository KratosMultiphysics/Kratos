//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
//                  Rahul Kikkeri Nagaraja
//

// Project includes
#include "custom_processes/apply_rans_chimera_process_monolithic.h"

namespace Kratos {

template <int TDim>
    ApplyRANSChimeraProcessMonolithic<TDim>::ApplyRANSChimeraProcessMonolithic(
        ModelPart& rMainModelPart, Parameters iParameters, SolvingVariablesVectorType rSolvingVariablesVector)
        : BaseType(rMainModelPart, iParameters), mSolvingVariableNamesVector(rSolvingVariablesVector)
    {
        BaseType::mNumberofSolvingVariables = mSolvingVariableNamesVector.size(); // updating the value set by BaseType constructor        
    }
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

template <int TDim>
    std::string ApplyRANSChimeraProcessMonolithic<TDim>::Info() const
    {
        return "ApplyRANSChimeraProcessMonolithic";
    }

template <int TDim>
    void ApplyRANSChimeraProcessMonolithic<TDim>::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyRANSChimeraProcessMonolithic" << std::endl;
    }

template <int TDim>
    void ApplyRANSChimeraProcessMonolithic<TDim>::PrintData(std::ostream& rOStream) const
    {
        KRATOS_INFO("ApplyRANSChimeraProcessMonolithic") << std::endl;
    }

template <int TDim>
    void ApplyRANSChimeraProcessMonolithic<TDim>::ApplyContinuityWithMpcs(ModelPart& rBoundaryModelPart, PointLocatorType& rBinLocator)
    {
        // loop over nodes and find the triangle in which it falls, then do
        // interpolation
        MasterSlaveConstraintContainerType master_slave_constraint_container;
        ApplyRANSChimeraProcessMonolithic::ReserveMemoryForConstraintContainers(
            rBoundaryModelPart, master_slave_constraint_container);

        const int n_boundary_nodes = static_cast<int>(rBoundaryModelPart.Nodes().size());
        for (int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn) {
            ModelPart::NodesContainerType::iterator i_boundary_node =
                rBoundaryModelPart.NodesBegin() + i_bn;
            BaseType::mNodeIdToConstraintIdsMap[i_boundary_node->Id()].reserve(150);
        }

        ApplyRANSChimeraProcessMonolithic::FormulateConstraints(rBoundaryModelPart, rBinLocator,
                                       master_slave_constraint_container);
        BuiltinTimer mpc_add_time;
        ApplyRANSChimeraProcessMonolithic::AddConstraintsToModelpart(BaseType::mrMainModelPart,
                                            master_slave_constraint_container);
        KRATOS_INFO_IF(
            "Adding of MPCs from containers to modelpart took         : ",
            BaseType::mEchoLevel > 1)
            << mpc_add_time.ElapsedSeconds() << " seconds" << std::endl;
    }
    
template <int TDim>
void ApplyRANSChimeraProcessMonolithic<TDim>::ReserveMemoryForConstraintContainers(
    ModelPart& rModelPart, MasterSlaveConstraintContainerType& rMsConstraintContainer)
{
    const IndexType num_of_constraints = rModelPart.NumberOfNodes() * 6; 
    // 6 Max value--> considering 3 velocities, 1 Pressure and 2 for RANS variables
    rMsConstraintContainer.reserve(num_of_constraints);
}

template <int TDim>
void ApplyRANSChimeraProcessMonolithic<TDim>::FormulateConstraints(
    ModelPart& rBoundaryModelPart,
    PointLocatorType& rBinLocator,
    MasterSlaveConstraintContainerType& rMsConstraintContainer)
{
    std::vector<int> vector_of_non_found_nodes;
    const int n_boundary_nodes = static_cast<int>(rBoundaryModelPart.Nodes().size());
    std::vector<int> constraints_id_vector;
    int num_constraints_required = (BaseType::mNumberofSolvingVariables) * (rBoundaryModelPart.Nodes().size());        
    BaseType::CreateConstraintIds(constraints_id_vector, num_constraints_required);

    IndexType found_counter = 0;
    IndexType removed_counter = 0;

    BuiltinTimer loop_over_b_nodes;
#pragma omp parallel for shared(constraints_id_vector,  \
                                rMsConstraintContainer, \
                                rBinLocator) reduction(+ : found_counter)
    for (int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn) {
        ModelPart::NodesContainerType::iterator i_boundary_node =
            rBoundaryModelPart.NodesBegin() + i_bn;
        NodeType& r_boundary_node = *(*(i_boundary_node.base()));
        Element::Pointer r_host_element;
        Vector weights;
        bool is_found = BaseType::SearchNode(rBinLocator, r_boundary_node, r_host_element, weights);
        if (is_found) {
            unsigned int start_constraint_id = i_bn * (BaseType::mNumberofSolvingVariables) * (r_host_element->GetGeometry().size());
            removed_counter += BaseType::RemoveExistingConstraintsForNode(r_boundary_node);
            ApplyRANSChimeraProcessMonolithic::MakeConstraints(r_boundary_node, r_host_element, 
                            weights, rMsConstraintContainer,
                            constraints_id_vector, start_constraint_id);
            found_counter += 1;
        }
    }

    double loop_time = loop_over_b_nodes.ElapsedSeconds();
    KRATOS_INFO_IF(
        "ApplyChimera : Loop over boundary nodes took             : ", BaseType::mEchoLevel > 0)
        << loop_time << " seconds" << std::endl;
    KRATOS_INFO_IF(
        "ApplyChimera : Number of Boundary nodes                  : ", BaseType::mEchoLevel > 1)
        << n_boundary_nodes << std::endl;
    KRATOS_INFO_IF(
        "ApplyChimera : Number of Boundary nodes found            : ", BaseType::mEchoLevel > 1)
        << found_counter << std::endl;
    KRATOS_INFO_IF(
        "ApplyChimera : Number of Boundary nodes not found        : ", BaseType::mEchoLevel > 1)
        << n_boundary_nodes - found_counter << std::endl;
    KRATOS_INFO_IF(
        "ApplyChimera : Number of constraints made                : ", BaseType::mEchoLevel > 1)
        << found_counter * (BaseType::mNumberofSolvingVariables) * 3 << std::endl;
    KRATOS_INFO_IF(
        "ApplyChimera : Number of constraints removed             : ", BaseType::mEchoLevel > 1)
        << removed_counter << std::endl;
}

template <int TDim>
void ApplyRANSChimeraProcessMonolithic<TDim>::MakeConstraints(
    NodeType& rNodeToFind,
    Element::Pointer& rHostElement,
    Vector& rWeights,
    MasterSlaveConstraintContainerType& rMsConstraintContainer,
    std::vector<int>& rConstraintIdVector,
    const IndexType StartConstraintId)
{
    Geometry<NodeType>& r_geom = rHostElement->GetGeometry();
    int init_index = 0;
    for (const auto& variable_name : mSolvingVariableNamesVector)
    {
        const Variable<double>& variable = KratosComponents<Variable<double>>::Get(variable_name);
        ApplyRANSChimeraProcessMonolithic::ApplyContinuityWithElement(r_geom, rNodeToFind, rWeights, variable,
                                             StartConstraintId + init_index,
                                             rConstraintIdVector, rMsConstraintContainer);
        init_index += (BaseType::mNumberofSolvingVariables);
    }
}

template <int TDim>
template <typename TVariableType>
void ApplyRANSChimeraProcessMonolithic<TDim>::AddMasterSlaveRelation(
    MasterSlaveConstraintContainerType& rMasterSlaveContainer,
    const LinearMasterSlaveConstraint& rCloneConstraint,
    unsigned int ConstraintId,
    NodeType& rMasterNode,
    TVariableType& rMasterVariable,
    NodeType& rSlaveNode,
    TVariableType& rSlaveVariable,
    const double Weight,
    const double Constant)
{
    rSlaveNode.Set(SLAVE);
    ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint =
        rCloneConstraint.Create(ConstraintId, rMasterNode, rMasterVariable,
                                rSlaveNode, rSlaveVariable, Weight, Constant);
    p_new_constraint->Set(TO_ERASE);
    BaseType::mNodeIdToConstraintIdsMap[rSlaveNode.Id()].push_back(ConstraintId);
    rMasterSlaveContainer.push_back(p_new_constraint);
}

template <int TDim>
template <typename TVariableType>
void ApplyRANSChimeraProcessMonolithic<TDim>::ApplyContinuityWithElement(
    Geometry<NodeType>& rGeometry,
    NodeType& rBoundaryNode,
    Vector& rShapeFuncWeights,
    TVariableType& rVariable,
    unsigned int StartIndex,
    std::vector<int>& rConstraintIdVector,
    MasterSlaveConstraintContainerType& rMsContainer)
{
    const auto& r_clone_constraint = new LinearMasterSlaveConstraint();
    // Initialise the boundary nodes dofs to 0 at ever time steps
    rBoundaryNode.FastGetSolutionStepValue(rVariable, 0) = 0.0;

    for (IndexType i = 0; i < rGeometry.size(); i++) {
        // Interpolation of rVariable
        rBoundaryNode.FastGetSolutionStepValue(rVariable, 0) +=
            rGeometry[i].GetDof(rVariable).GetSolutionStepValue(0) * rShapeFuncWeights[i];
        ApplyRANSChimeraProcessMonolithic::AddMasterSlaveRelation(rMsContainer, *r_clone_constraint,
                               rConstraintIdVector[StartIndex++], rGeometry[i], rVariable,
                               rBoundaryNode, rVariable, rShapeFuncWeights[i]);
    } // end of loop over host element nodes

    // Setting the buffer 1 same buffer 0
    rBoundaryNode.FastGetSolutionStepValue(rVariable, 1) =
        rBoundaryNode.FastGetSolutionStepValue(rVariable, 0);
}

template <int TDim>
void ApplyRANSChimeraProcessMonolithic<TDim>::AddConstraintsToModelpart(ModelPart& rModelPart,
                                                   MasterSlaveConstraintContainerType& rMsConstraintContainer)
{
    IndexType n_total_constraints = static_cast<int>(rMsConstraintContainer.size());

    auto& constraints = rModelPart.MasterSlaveConstraints();
    constraints.reserve(n_total_constraints);

    auto& constraints_data = constraints.GetContainer();
    constraints_data.insert(constraints_data.end(), rMsConstraintContainer.ptr_begin(), rMsConstraintContainer.ptr_end());

    constraints.Sort();
}

// Template declarations
template class ApplyRANSChimeraProcessMonolithic<2>;
// template class ApplyRANSChimeraProcessMonolithic<3>;

} // namespace Kratos.

