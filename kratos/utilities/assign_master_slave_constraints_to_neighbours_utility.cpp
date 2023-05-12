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
AssignMasterSlaveConstraintsToNeighboursUtility::AssignMasterSlaveConstraintsToNeighboursUtility(NodesContainerType& rStructureNodes){
    KRATOS_TRY;
    NodesContainerType::ContainerType& nodes_model_part = rStructureNodes.GetContainer();
    mpBins = Kratos::make_unique<NodeBinsType>(nodes_model_part.begin(), nodes_model_part.end());
    mMaxNumberOfNodes = rStructureNodes.size();
    KRATOS_CATCH(""); 
} 

// Search for the neighbouring nodes (in rStructureNodes) of each rNode on a given array of rNodes
void AssignMasterSlaveConstraintsToNeighboursUtility::SearchNodesInRadiusForNodes(
    const NodesContainerType& pNodes,
    const double Radius,
    const double MinNumOfNeighNodes,
    VectorResultNodesContainerType& rResults)
{
    KRATOS_TRY;
    auto& p_nodes_array = pNodes.GetContainer();

    // Resize rResults vector
    rResults.resize(p_nodes_array.size());

    // Declare a counter variable outside the loop as std::atomic<int>
    std::atomic<int> i(0);

    block_for_each(p_nodes_array, [&](const Node<3>::Pointer& pNode)
    {
        double local_radius = Radius;

        IndexType it = i.fetch_add(1); // Atomically increment the counter and get the previous value

        while (rResults[it].size()<MinNumOfNeighNodes){
            rResults[it].clear();
            SearchNodesInRadiusForNode(pNode, local_radius, rResults[it]);
            local_radius += Radius;
        }
    });

    // Optimize memory usage
    rResults.shrink_to_fit();

    KRATOS_CATCH("");
}


//Search for the neighbouring nodes (in rStructureNodes) of the given rNode
void AssignMasterSlaveConstraintsToNeighboursUtility::SearchNodesInRadiusForNode(
    NodeType::Pointer pNode,
    double const Radius,
    ResultNodesContainerType& rResults)
{

    KRATOS_TRY;

    ResultNodesContainerType local_results(mMaxNumberOfNodes);
    std::size_t num_of_results = 0;
    auto results_iterator = local_results.begin();
    num_of_results = mpBins->SearchInRadius(*pNode, Radius, results_iterator, mMaxNumberOfNodes);
    rResults.insert(rResults.begin(), local_results.begin(), local_results.begin() + num_of_results);

    KRATOS_CATCH("");
}

const Variable<double>& AssignMasterSlaveConstraintsToNeighboursUtility::GetComponentVariable(
    const Variable<array_1d<double, 3>>& rVectorVariable, 
    const std::size_t ComponentIndex)
{
    static const std::array<std::string, 3> component_suffixes = {"_X", "_Y", "_Z"};

    if (ComponentIndex < 3) {
        std::string component_variable_name = rVectorVariable.Name() + component_suffixes[ComponentIndex];
        return KratosComponents<Variable<double>>::Get(component_variable_name);
    } else {
        KRATOS_ERROR << "Component variable for " << rVectorVariable.Name() << " at index " << ComponentIndex << " not found" << std::endl;
    }
}

void AssignMasterSlaveConstraintsToNeighboursUtility::GetDofsAndCoordinatesForNode(
    NodeType::Pointer pNode,
    const Variable<double>& rVariable,
    DofPointerVectorType& rCloudOfDofs,
    array_1d<double,3>& rSlaveCoordinates
    )
{
    KRATOS_TRY;
    rCloudOfDofs.push_back(pNode->pGetDof(rVariable));
    rSlaveCoordinates = pNode->Coordinates();
    KRATOS_CATCH("");
}

void AssignMasterSlaveConstraintsToNeighboursUtility::GetDofsAndCoordinatesForNode(
    NodeType::Pointer pNode,
    std::array<const Kratos::Variable<double>*, 3> ComponentVariables,
    DofPointerVectorType& rCloudOfDofsX,
    DofPointerVectorType& rCloudOfDofsY,
    DofPointerVectorType& rCloudOfDofsZ,
    array_1d<double,3>& rSlaveCoordinates
    )
{
    KRATOS_TRY;

    // Check if the node has the required DOFs
    if (pNode->HasDofFor(*ComponentVariables[0]) &&
        pNode->HasDofFor(*ComponentVariables[1]) &&
        pNode->HasDofFor(*ComponentVariables[2])) {

        rCloudOfDofsX.push_back(pNode->pGetDof(*ComponentVariables[0]));
        rCloudOfDofsY.push_back(pNode->pGetDof(*ComponentVariables[1]));
        rCloudOfDofsZ.push_back(pNode->pGetDof(*ComponentVariables[2]));

    } else {
            std::stringstream variables_str;
            variables_str << "[" << *ComponentVariables[0] << ", " << *ComponentVariables[1] << ", " << *ComponentVariables[2] << "]";
            KRATOS_ERROR << "The node with ID " << pNode->Id() << " does not have the required DOFs for one of the variables " << variables_str.str() << std::endl;
    }
    rSlaveCoordinates = pNode->Coordinates();
    KRATOS_CATCH("");
}

// Get Dofs and Coordinates arrays for a given variable double. (For nodes)
void AssignMasterSlaveConstraintsToNeighboursUtility::GetDofsAndCoordinatesForNodes(
    ResultNodesContainerType NodesArray,
    const Variable<double>& rVariable,
    DofPointerVectorType& rCloudOfDofs,
    Matrix& rCloudOfNodesCoordinates
    )
{
    KRATOS_TRY;
    rCloudOfDofs.resize(NodesArray.size());
    rCloudOfNodesCoordinates.resize(NodesArray.size(),3);
    for(unsigned int i = 0; i < NodesArray.size(); i++){
        rCloudOfDofs[i] = NodesArray[i]->pGetDof(rVariable);
        noalias(row(rCloudOfNodesCoordinates,i)) = NodesArray[i]->Coordinates();
    }
    KRATOS_CATCH("");
}

void AssignMasterSlaveConstraintsToNeighboursUtility::GetDofsAndCoordinatesForNodes(
    ResultNodesContainerType NodesArray,
    std::array<const Kratos::Variable<double>*, 3> ComponentVariables,
    DofPointerVectorType& rCloudOfDofsX,
    DofPointerVectorType& rCloudOfDofsY,
    DofPointerVectorType& rCloudOfDofsZ,
    Matrix& rCloudOfNodesCoordinates
)
{
    KRATOS_TRY;
    rCloudOfDofsX.resize(NodesArray.size());
    rCloudOfDofsY.resize(NodesArray.size());
    rCloudOfDofsZ.resize(NodesArray.size());
    rCloudOfNodesCoordinates.resize(NodesArray.size(), 3);

    for(int i = 0; i < static_cast<int>(NodesArray.size()); i++)
    {
        NodeType::Pointer pNode = NodesArray[i];

        // Check if the node has the required DOFs
        if (pNode->HasDofFor(*ComponentVariables[0]) &&
            pNode->HasDofFor(*ComponentVariables[1]) &&
            pNode->HasDofFor(*ComponentVariables[2])){
          rCloudOfDofsX[i] = pNode->pGetDof(*ComponentVariables[0]);
          rCloudOfDofsY[i] = pNode->pGetDof(*ComponentVariables[1]);
          rCloudOfDofsZ[i] = pNode->pGetDof(*ComponentVariables[2]);

        } else {
            std::stringstream variables_str;
            variables_str << "[" << *ComponentVariables[0] << ", " << *ComponentVariables[1] << ", " << *ComponentVariables[2] << "]";
            KRATOS_ERROR << "The node with ID " << pNode->Id() << " does not have the required DOFs for one of the variables " << variables_str.str() << std::endl;
        }
        noalias(row(rCloudOfNodesCoordinates, i)) = pNode->Coordinates();
    }
    KRATOS_CATCH("");
}

void AssignMasterSlaveConstraintsToNeighboursUtility::AssignMPCsToNodes(
    NodesContainerType pNodes,
    double const Radius,
    ModelPart& rComputingModelPart,
    const Variable<double>& rVariable,
    double const MinNumOfNeighNodes
    )
{
    KRATOS_TRY;

    BuiltinTimer build_and_assign_mpcs;

    int prev_num_mpcs = rComputingModelPart.NumberOfMasterSlaveConstraints();

    auto& p_nodes_array = pNodes.GetContainer();

    ModelPart::MasterSlaveConstraintType const& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");

    // Create a concurrent queue for constraints
    using MasterSlaveConstarintQueue = moodycamel::ConcurrentQueue<ModelPart::MasterSlaveConstraintType::Pointer>;

    MasterSlaveConstarintQueue concurrent_constraints;

    // Declare a counter variable outside the loop as std::atomic<int>
    std::atomic<int> i(0);

    block_for_each(p_nodes_array, [&](Node<3>::Pointer& pNode) {

        double local_radius = Radius;

        ResultNodesContainerType r_result(0);

        while (r_result.size() < MinNumOfNeighNodes) {
            r_result.clear();
            SearchNodesInRadiusForNode(pNode, local_radius, r_result);
            local_radius += Radius;
        }

        // Get Dofs and Coordinates
        DofPointerVectorType r_cloud_of_dofs, r_slave_dof;
        Matrix r_cloud_of_nodes_coordinates;
        array_1d<double, 3> r_slave_coordinates;
        GetDofsAndCoordinatesForNode(pNode, rVariable, r_slave_dof, r_slave_coordinates);
        GetDofsAndCoordinatesForNodes(r_result, rVariable, r_cloud_of_dofs, r_cloud_of_nodes_coordinates);

        // Calculate shape functions
        Vector r_n_container;
        RBFShapeFunctionsUtility::CalculateShapeFunctions(r_cloud_of_nodes_coordinates, r_slave_coordinates, r_n_container);

        // Create MPCs
        Matrix shape_matrix(1, r_n_container.size());
        noalias(row(shape_matrix, 0)) = r_n_container; // Shape functions matrix
        const Vector constant_vector = ZeroVector(r_n_container.size());
        IndexType it = i.fetch_add(1); // Atomically increment the counter and get the previous value
        concurrent_constraints.enqueue(r_clone_constraint.Create(prev_num_mpcs + it + 1, r_cloud_of_dofs, r_slave_dof, shape_matrix, constant_vector));
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

    KRATOS_INFO("AssignMasterSlaveConstraintsToNeighboursUtility") << "Build and Assign MPCs Time: " << build_and_assign_mpcs.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("");
}

void AssignMasterSlaveConstraintsToNeighboursUtility::AssignMPCsToNodes(
    NodesContainerType pNodes,
    double const Radius,
    ModelPart& rComputingModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    double const MinNumOfNeighNodes
    )
{
    KRATOS_TRY;

    BuiltinTimer build_and_assign_mpcs;
    
    // Create array of component variables
    std::array<const Kratos::Variable<double>*, 3> component_variables = {
        &GetComponentVariable(rVariable, 0),
        &GetComponentVariable(rVariable, 1),
        &GetComponentVariable(rVariable, 2)
    };

    int prev_num_mpcs = rComputingModelPart.NumberOfMasterSlaveConstraints();

    auto& p_nodes_array = pNodes.GetContainer();

    ModelPart::MasterSlaveConstraintType const& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");

    // Create a concurrent queue for constraints
    using MasterSlaveConstarintQueue = moodycamel::ConcurrentQueue<ModelPart::MasterSlaveConstraintType::Pointer>;

    MasterSlaveConstarintQueue concurrent_constraints;

    // Declare a counter variable outside the loop as std::atomic<int>
    std::atomic<int> i(0);

    block_for_each(p_nodes_array, [&](Node<3>::Pointer& pNode) {

        double local_radius = Radius;

        ResultNodesContainerType r_result(0);

        while (r_result.size() < MinNumOfNeighNodes) {
            r_result.clear();
            SearchNodesInRadiusForNode(pNode, local_radius, r_result);
            local_radius += Radius;
        }

        // Get Dofs and Coordinates
        DofPointerVectorType r_cloud_of_dofs_x,r_cloud_of_dofs_y,r_cloud_of_dofs_z,r_slave_dof_x,r_slave_dof_y,r_slave_dof_z; //Struct local storage
        Matrix r_cloud_of_nodes_coordinates; 
        array_1d<double, 3> r_slave_coordinates; //Struct local storage
        GetDofsAndCoordinatesForNode(pNode, component_variables, r_slave_dof_x,r_slave_dof_y,r_slave_dof_z, r_slave_coordinates);
        GetDofsAndCoordinatesForNodes(r_result, component_variables, r_cloud_of_dofs_x, r_cloud_of_dofs_y, r_cloud_of_dofs_z, r_cloud_of_nodes_coordinates);
        
        // Calculate shape functions
        Vector r_n_container;//Struct local storage
        RBFShapeFunctionsUtility::CalculateShapeFunctions(r_cloud_of_nodes_coordinates, r_slave_coordinates, r_n_container);

        // Create MPCs
        Matrix shape_matrix(1, r_n_container.size());
        noalias(row(shape_matrix, 0)) = r_n_container; // Shape functions matrix
        const Vector constant_vector = ZeroVector(r_n_container.size());
        IndexType it = i.fetch_add(1)*3; // Atomically increment the counter and get the previous value
        concurrent_constraints.enqueue(r_clone_constraint.Create(prev_num_mpcs + it + 1, r_cloud_of_dofs_x, r_slave_dof_x, shape_matrix, constant_vector));
        concurrent_constraints.enqueue(r_clone_constraint.Create(prev_num_mpcs + it + 2, r_cloud_of_dofs_y, r_slave_dof_y, shape_matrix, constant_vector));
        concurrent_constraints.enqueue(r_clone_constraint.Create(prev_num_mpcs + it + 3, r_cloud_of_dofs_z, r_slave_dof_z, shape_matrix, constant_vector));
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

    KRATOS_INFO("AssignMasterSlaveConstraintsToNeighboursUtility") << "Build and Assign MPCs Time: " << build_and_assign_mpcs.ElapsedSeconds() << std::endl;
    KRATOS_CATCH("");
}
}  // namespace Kratos.