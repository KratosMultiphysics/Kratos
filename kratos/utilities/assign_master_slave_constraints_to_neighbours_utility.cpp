//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Sebastian Ares de Parga Regalado
//

// System includes
#include <string>
#include <iostream>

// External includes
#include "concurrentqueue/concurrentqueue.h"

// Project includes
#include "utilities/assign_master_slave_constraints_to_neighbours_utility.h"

namespace Kratos 
{
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/**
 * @class AssignMasterSlaveConstraintsToNeighboursUtility
 * @ingroup Kratos Core
 * @brief Assing Master-Slave Constraints to Neighbouring Nodes
 * @details This class provides a method to search for neighbouring nodes of one node
 * within a given radius and assign a master-slave constraint using radial basis functions.
 * @author Sebastian Ares de Parga Regalado
 */

/// Default constructor.
AssignMasterSlaveConstraintsToNeighboursUtility::AssignMasterSlaveConstraintsToNeighboursUtility(NodesContainerType& rMasterStructureNodes){
    KRATOS_TRY;
    NodesContainerType::ContainerType& nodes_model_part = rMasterStructureNodes.GetContainer();
    mpBins = Kratos::make_unique<NodeBinsType>(nodes_model_part.begin(), nodes_model_part.end());
    mMaxNumberOfNodes = rMasterStructureNodes.size();
    KRATOS_CATCH(""); 
} 

//Search for the neighbouring nodes (in rMasterStructureNodes) of the given rNode
void AssignMasterSlaveConstraintsToNeighboursUtility::SearchNodesInRadiusForNode(
    NodeType::Pointer pSlaveNode,
    double const Radius,
    ResultNodesContainerType& rResults)
{

    KRATOS_TRY;

    ResultNodesContainerType local_results(mMaxNumberOfNodes);
    std::size_t num_of_results = 0;
    auto results_iterator = local_results.begin();
    num_of_results = mpBins->SearchInRadius(*pSlaveNode, Radius, results_iterator, mMaxNumberOfNodes);
    rResults.insert(rResults.begin(), local_results.begin(), local_results.begin() + num_of_results);

    KRATOS_CATCH("");
}

void AssignMasterSlaveConstraintsToNeighboursUtility::GetDofsAndCoordinatesForSlaveNode(
    NodeType::Pointer pSlaveNode,
    const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList,
    std::vector<DofPointerVectorType>& rSlaveDofs,
    array_1d<double, 3>& rSlaveCoordinates
)
{
    KRATOS_TRY;

    // Check if the node has the required DOFs for all variables in the list
    bool has_all_dofs = true;
    for (const auto& variable : rVariableList) {
        if (!pSlaveNode->HasDofFor(variable.get())) {
            has_all_dofs = false;
            break;
        }
    }

    if (has_all_dofs) {
        // Get Dofs for each variable in the list
        rSlaveDofs.resize(rVariableList.size());
        for (std::size_t i = 0; i < rVariableList.size(); ++i) {
            const auto& variable = rVariableList[i];
            rSlaveDofs[i].resize(1);
            rSlaveDofs[i][0] = pSlaveNode->pGetDof(variable.get());
        }
    } else {
        std::stringstream variables_str;
        variables_str << "[";
        for (const auto& variable : rVariableList) {
            variables_str << variable.get() << ", ";
        }
        variables_str.seekp(-2, std::ios_base::end); // Remove the last ", "
        variables_str << "]";
        KRATOS_ERROR << "The node with ID " << pSlaveNode->Id() << " does not have the required DOFs for one of the variables " << variables_str.str() << std::endl;
    }

    rSlaveCoordinates = pSlaveNode->Coordinates();

    KRATOS_CATCH("");
}

// Get Dofs and Coordinates arrays for a given variable list. (For nodes)
void AssignMasterSlaveConstraintsToNeighboursUtility::GetDofsAndCoordinatesForCloudOfNodes(
    const ResultNodesContainerType& CloudOfNodesArray,
    const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList,
    std::vector<DofPointerVectorType>& rCloudOfDofs,
    Matrix& rCloudOfNodesCoordinates
)
{
    KRATOS_TRY;
    rCloudOfDofs.resize(rVariableList.size());
    rCloudOfNodesCoordinates.resize(CloudOfNodesArray.size(), 3);

    // Resize rCloudOfDofs
    for (auto& dofs : rCloudOfDofs)
    {
        dofs.resize(CloudOfNodesArray.size());
    }

    for (int i = 0; i < static_cast<int>(CloudOfNodesArray.size()); i++)
    {
        NodeType::Pointer pSlaveNode = CloudOfNodesArray[i];

        for (int j = 0; j < static_cast<int>(rVariableList.size()); j++)
        {
            const auto& rVariable = rVariableList[j];

            if (pSlaveNode->HasDofFor(rVariable.get()))
            {
                rCloudOfDofs[j][i] = pSlaveNode->pGetDof(rVariable.get());
            }
            else
            {
                std::stringstream variable_str;
                variable_str << rVariable.get();
                KRATOS_ERROR << "The node with ID " << pSlaveNode->Id() << " does not have the required DOF for the variable " << variable_str.str() << std::endl;
            }
        }

        noalias(row(rCloudOfNodesCoordinates, i)) = pSlaveNode->Coordinates();
    }
    KRATOS_CATCH("");
}

void AssignMasterSlaveConstraintsToNeighboursUtility::AssignMasterSlaveConstraintsToNodes(
    NodesContainerType pSlaveNodes,
    double const Radius,
    ModelPart& rComputingModelPart,
    const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList,
    double const MinNumOfNeighNodes
    )
{
    KRATOS_TRY;

    BuiltinTimer build_and_assign_mscs;

    // Declare thread-local storage (TLS) variables
    thread_local std::vector<DofPointerVectorType> r_cloud_of_dofs;
    thread_local std::vector<DofPointerVectorType> r_slave_dofs;
    thread_local Matrix r_cloud_of_nodes_coordinates;
    thread_local array_1d<double, 3> r_slave_coordinates;
    thread_local Vector r_n_container;
    thread_local Vector constant_vector;
    thread_local Matrix shape_matrix;

    int prev_num_mscs = rComputingModelPart.NumberOfMasterSlaveConstraints();

    if (prev_num_mscs > 0) {
        KRATOS_WARNING("AssignMasterSlaveConstraintsToNeighboursUtility") << "Previous Master-Slave Constraints exist in the ModelPart. The new Constraints may interact with the existing ones." << std::endl;
    }

    auto& p_nodes_array = pSlaveNodes.GetContainer();

    const auto& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");

    // Create a concurrent queue for constraints
    using MasterSlaveConstarintQueue = moodycamel::ConcurrentQueue<ModelPart::MasterSlaveConstraintType::Pointer>;

    MasterSlaveConstarintQueue concurrent_constraints;

    // Declare a counter variable outside the loop as std::atomic<int>
    std::atomic<int> i(0);

    block_for_each(p_nodes_array, [&](Node<3>::Pointer& pSlaveNode) {

        double local_radius = Radius;

        ResultNodesContainerType r_cloud_of_nodes(0);

        while (r_cloud_of_nodes.size() < MinNumOfNeighNodes) {
            r_cloud_of_nodes.clear();
            SearchNodesInRadiusForNode(pSlaveNode, local_radius, r_cloud_of_nodes);
            local_radius += Radius;
        }

        // Get Dofs and Coordinates
        GetDofsAndCoordinatesForSlaveNode(pSlaveNode, rVariableList, r_slave_dofs, r_slave_coordinates);
        GetDofsAndCoordinatesForCloudOfNodes(r_cloud_of_nodes, rVariableList, r_cloud_of_dofs, r_cloud_of_nodes_coordinates);
        
        // Calculate shape functions
        RBFShapeFunctionsUtility::CalculateShapeFunctions(r_cloud_of_nodes_coordinates, r_slave_coordinates, r_n_container);

        // Create MasterSlaveConstraints
        shape_matrix.resize(1, r_n_container.size());
        noalias(row(shape_matrix, 0)) = r_n_container; // Shape functions matrix
        constant_vector.resize(r_n_container.size());
        IndexType it = i.fetch_add(1) * rVariableList.size(); // Atomically increment the counter and get the previous value

        for (std::size_t j = 0; j < rVariableList.size(); ++j)
        {
            concurrent_constraints.enqueue(r_clone_constraint.Create(prev_num_mscs + it + j + 1, const_cast<DofPointerVectorType&>(r_cloud_of_dofs[j]), r_slave_dofs[j], shape_matrix, constant_vector));
        }
    });

    // Create a temporary container to store the master-slave constraints
    ConstraintContainerType temp_constraints;
    ModelPart::MasterSlaveConstraintType::Pointer p_constraint;

    // Dequeue the constraints from the concurrent queue and add them to the temporary container
    while (concurrent_constraints.try_dequeue(p_constraint)) {
        temp_constraints.push_back(p_constraint);
    }

    // Add the constraints to the rComputingModelPart in a single call
    rComputingModelPart.AddMasterSlaveConstraints(temp_constraints.begin(), temp_constraints.end());

    KRATOS_INFO("AssignMasterSlaveConstraintsToNeighboursUtility") << "Build and Assign Master-Slave Constraints Time: " << build_and_assign_mscs.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("");
}
}  // namespace Kratos.

