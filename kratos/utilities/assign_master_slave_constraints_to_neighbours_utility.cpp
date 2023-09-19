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

// Destructor
AssignMasterSlaveConstraintsToNeighboursUtility::~AssignMasterSlaveConstraintsToNeighboursUtility() {}

// Search for neighboring nodes (in rMasterStructureNodes) of the given slave node
void AssignMasterSlaveConstraintsToNeighboursUtility::SearchNodesInRadiusForNode(
    NodeType::Pointer pSlaveNode,
    double const Radius,
    ResultNodesContainerType& rCloudOfNodes,
    ResultNodesContainerType& rLocalResults)
{

    KRATOS_TRY;

    // ResultNodesContainerType local_results(mMaxNumberOfNodes);
    std::size_t num_of_results = 0;
    auto results_iterator = rLocalResults.begin();
    num_of_results = mpBins->SearchInRadius(*pSlaveNode, Radius, results_iterator, mMaxNumberOfNodes);
    rCloudOfNodes.insert(rCloudOfNodes.begin(), rLocalResults.begin(), rLocalResults.begin() + num_of_results);

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
        // Get DOFs for each variable in the list
        rSlaveDofs.resize(rVariableList.size());
        for (std::size_t i = 0; i < rVariableList.size(); ++i) {
            const auto& variable = rVariableList[i];
            rSlaveDofs[i].resize(1);
            rSlaveDofs[i][0] = pSlaveNode->pGetDof(variable.get());
        }
    } else {
        std::stringstream variables_str;
        variables_str << "[";
        for (const auto& variable : rVariableList) 
            variables_str << variable.get() << ", ";
        variables_str.seekp(-2, std::ios_base::end); // Remove the last ", "
        variables_str << "]";
        KRATOS_ERROR << "The node with ID " << pSlaveNode->Id() << " does not have the required DOFs for one of the variables " << variables_str.str() << std::endl;
    }

    rSlaveCoordinates = pSlaveNode->Coordinates();

    KRATOS_CATCH("");
}

// Get DOFs and Coordinates arrays for a given variable list for a cloud of nodes
void AssignMasterSlaveConstraintsToNeighboursUtility::GetDofsAndCoordinatesForCloudOfNodes(
    const ResultNodesContainerType& CloudOfNodesArray,
    const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList,
    std::vector<DofPointerVectorType>& rCloudOfDofs,
    Matrix& rCloudOfNodesCoordinates
)
{
    KRATOS_TRY;

    // Resize rCloudOfDofs to match the size of the variable list
    rCloudOfDofs.resize(rVariableList.size());

    // Resize rCloudOfDofs vectors to match the size of the cloud of nodes
    for (auto& dofs : rCloudOfDofs)
        dofs.resize(CloudOfNodesArray.size());

    // Resize rCloudOfNodesCoordinates to match the size of the cloud of nodes
    rCloudOfNodesCoordinates.resize(CloudOfNodesArray.size(), 3);

    // Loop over the cloud of nodes
    for (std::size_t i = 0; i < CloudOfNodesArray.size(); i++) {
        auto p_slave_node = CloudOfNodesArray[i];

        // Loop over the variable list
        for (std::size_t j = 0; j < rVariableList.size(); j++) {
            const auto& r_variable = rVariableList[j];

            // Check if the slave node has the required DOF for the current variable
            if (p_slave_node->HasDofFor(r_variable.get())) {
                rCloudOfDofs[j][i] = p_slave_node->pGetDof(r_variable.get()); // Assign the DOF pointer to the corresponding position in rCloudOfDofs
            } else {
                std::stringstream variable_str;
                variable_str << r_variable.get();
                KRATOS_ERROR << "The node with ID " << p_slave_node->Id() << " does not have the required DOF for the variable " << variable_str.str() << std::endl;
            }
        }

        // Assign the coordinates of the slave node to the corresponding row in rCloudOfNodesCoordinates
        noalias(row(rCloudOfNodesCoordinates, i)) = p_slave_node->Coordinates();
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
    thread_local std::vector<DofPointerVectorType> r_cloud_of_dofs;           // The DOFs of the cloud of nodes
    thread_local std::vector<DofPointerVectorType> r_slave_dofs;              // The DOFs of the slave node
    thread_local Matrix r_cloud_of_nodes_coordinates;                          // The coordinates of the cloud of nodes
    thread_local array_1d<double, 3> r_slave_coordinates;                      // The coordinates of the slave node
    thread_local Vector r_n_container;                                         // Container for shape functions
    thread_local Vector constant_vector;                                       // Constant vector for shape functions
    thread_local Matrix shape_matrix;                                          // Shape functions matrix

    int prev_num_mscs = rComputingModelPart.NumberOfMasterSlaveConstraints();

    KRATOS_WARNING_IF("AssignMasterSlaveConstraintsToNeighboursUtility", prev_num_mscs > 0) << "Previous Master-Slave Constraints exist in the ModelPart. The new Constraints may interact with the existing ones." << std::endl;

    auto& p_nodes_array = pSlaveNodes.GetContainer();  // Get the container of slave nodes

    const auto& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");  // Get the clone constraint

    // Create a concurrent queue for constraints
    using MasterSlaveConstraintQueue = moodycamel::ConcurrentQueue<ModelPart::MasterSlaveConstraintType::Pointer>;
    MasterSlaveConstraintQueue concurrent_constraints;

    // Declare a counter variable outside the loop as std::atomic<int>
    std::atomic<int> i(0);

    block_for_each(p_nodes_array, [&](Node::Pointer& p_slave_node) {

        double local_radius = Radius;

        ResultNodesContainerType r_cloud_of_nodes(0);  // Container for the cloud of nodes
        ResultNodesContainerType r_local_results(mMaxNumberOfNodes);  // Container for node search results

        while (r_cloud_of_nodes.size() < MinNumOfNeighNodes) {
            r_cloud_of_nodes.clear();
            SearchNodesInRadiusForNode(p_slave_node, local_radius, r_cloud_of_nodes, r_local_results);  // Search for nodes in the radius
            local_radius += Radius;
        }

        // Get Dofs and Coordinates for the slave node and the cloud of nodes
        GetDofsAndCoordinatesForSlaveNode(p_slave_node, rVariableList, r_slave_dofs, r_slave_coordinates);
        GetDofsAndCoordinatesForCloudOfNodes(r_cloud_of_nodes, rVariableList, r_cloud_of_dofs, r_cloud_of_nodes_coordinates);
        
        // Calculate shape functions
        RBFShapeFunctionsUtility::CalculateShapeFunctions(r_cloud_of_nodes_coordinates, r_slave_coordinates, r_n_container);

        // Create MasterSlaveConstraints
        shape_matrix.resize(1, r_n_container.size());
        noalias(row(shape_matrix, 0)) = r_n_container; // Shape functions matrix
        constant_vector.resize(r_n_container.size());
        IndexType it = i.fetch_add(1) * rVariableList.size(); // Atomically increment the counter and get the previous value

        for (std::size_t j = 0; j < rVariableList.size(); ++j) {
            concurrent_constraints.enqueue(r_clone_constraint.Create(prev_num_mscs + it + j + 1, const_cast<DofPointerVectorType&>(r_cloud_of_dofs[j]), r_slave_dofs[j], shape_matrix, constant_vector));
        }
    });

    // Create a temporary container to store the master-slave constraints
    ConstraintContainerType temp_constraints;
    ModelPart::MasterSlaveConstraintType::Pointer p_constraint;

    // Dequeue the constraints from the concurrent queue and add them to the temporary container
    while (concurrent_constraints.try_dequeue(p_constraint))
        temp_constraints.push_back(p_constraint);

    // Add the constraints to the rComputingModelPart in a single call
    rComputingModelPart.AddMasterSlaveConstraints(temp_constraints.begin(), temp_constraints.end());

    KRATOS_INFO("AssignMasterSlaveConstraintsToNeighboursUtility") << "Build and Assign Master-Slave Constraints Time: " << build_and_assign_mscs.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("");
}
}  // namespace Kratos.

