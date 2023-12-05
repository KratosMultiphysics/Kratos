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
        auto p_slave_node = CloudOfNodesArray[i];//////CHANGE SLAVE NODE NAME

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

void AssignMasterSlaveConstraintsToNeighboursUtility::GetDofsAndCoordinatesForCloudOfNodesFlattened(
    const ResultNodesContainerType& CloudOfNodesArray,
    const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList,
    DofPointerVectorType& rCloudOfDofs,
    Matrix& rCloudOfNodesCoordinates
)
{
    KRATOS_TRY;

    // Resize rCloudOfDofs to accommodate all DOFs for all variables
    rCloudOfDofs.resize(rVariableList.size() * CloudOfNodesArray.size());

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
                std::size_t counter = i * rVariableList.size() + j;
                rCloudOfDofs[counter] = p_slave_node->pGetDof(r_variable.get());
            } else {
                KRATOS_ERROR << "The node with ID " << p_slave_node->Id() << " does not have the required DOF for the variable " << r_variable.get().Name() << std::endl;
            }
        }

        // Assign the coordinates of the slave node to the corresponding row in rCloudOfNodesCoordinates
        noalias(row(rCloudOfNodesCoordinates, i)) = p_slave_node->Coordinates();
    }

    KRATOS_CATCH("");
}

struct CustomComparator {
    static constexpr double TOLERANCE = 1e-6;

    bool operator()(const std::tuple<double, double, double>& a, const std::tuple<double, double, double>& b) const {
        return std::fabs(std::get<0>(a) - std::get<0>(b)) > TOLERANCE ||
               std::fabs(std::get<1>(a) - std::get<1>(b)) > TOLERANCE ||
               std::fabs(std::get<2>(a) - std::get<2>(b)) > TOLERANCE;
    }
};

void AssignMasterSlaveConstraintsToNeighboursUtility::FilterUniqueNodesForRBF(
    ResultNodesContainerType& CloudOfNodesArray,
    Matrix& rCloudOfNodesCoordinates
) {
    KRATOS_TRY;

    // Use a set with a custom comparator to keep track of unique coordinates
    std::set<std::tuple<double, double, double>, CustomComparator> unique_coords;

    // Temporary vector to store unique nodes
    std::vector<NodeType::Pointer> temp_nodes;

    // Loop over the cloud of nodes
    for (auto& p_node : CloudOfNodesArray) {
        auto coords = p_node->Coordinates();
        std::tuple<double, double, double> current_coords(coords[0], coords[1], coords[2]);

        // Check if the coordinates are unique
        if (unique_coords.find(current_coords) == unique_coords.end()) {
            unique_coords.insert(current_coords);
            temp_nodes.push_back(p_node);
        }
    }

    // Clear the original CloudOfNodesArray and repopulate with unique nodes
    CloudOfNodesArray.clear();
    for (auto& p_node : temp_nodes) {
        CloudOfNodesArray.push_back(p_node);
    }

    // Resize rCloudOfNodesCoordinates to match the size of unique nodes
    rCloudOfNodesCoordinates.resize(temp_nodes.size(), 3);

    // Fill rCloudOfNodesCoordinates with unique coordinates
    for (std::size_t i = 0; i < temp_nodes.size(); i++) {
        auto coords = temp_nodes[i]->Coordinates();
        for (std::size_t j = 0; j < 3; j++) {
            rCloudOfNodesCoordinates(i, j) = coords[j];
        }
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

        if (p_slave_node->Id() == 929){
            KRATOS_WATCH(r_cloud_of_nodes)
            KRATOS_WATCH(r_n_container)
            KRATOS_WATCH(shape_matrix)
        }

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

void AssignMasterSlaveConstraintsToNeighboursUtility::AssignRotationalMasterSlaveConstraintsToNodes(
    NodesContainerType pSlaveNodes,
    double const Radius,
    ModelPart& rComputingModelPart,
    const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList,
    double const MinNumOfNeighNodes,
    const std::unordered_map<IndexType, std::pair<IndexType, double>>& GhostToOriginalMapping
)
{
    KRATOS_TRY;

    BuiltinTimer build_and_assign_mscs;

    thread_local std::vector<DofPointerVectorType> r_cloud_of_dofs;
    thread_local DofPointerVectorType r_cloud_of_master_dofs;
    thread_local std::vector<DofPointerVectorType> r_slave_dofs;
    thread_local Matrix r_cloud_of_nodes_coordinates;
    thread_local Matrix r_cloud_of_master_nodes_coordinates;
    thread_local array_1d<double, 3> r_slave_coordinates;
    thread_local Vector r_n_container;
    thread_local Vector constant_vector;
    thread_local Matrix shape_matrix_x, shape_matrix_y, shape_matrix_z;

    int prev_num_mscs = rComputingModelPart.NumberOfMasterSlaveConstraints();

    KRATOS_WARNING_IF("AssignRotationalMasterSlaveConstraintsToNodes", prev_num_mscs > 0) << "Previous Master-Slave Constraints exist in the ModelPart. The new Constraints may interact with the existing ones." << std::endl;

    auto& p_nodes_array = pSlaveNodes.GetContainer();

    const auto& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");

    using MasterSlaveConstraintQueue = moodycamel::ConcurrentQueue<ModelPart::MasterSlaveConstraintType::Pointer>;
    MasterSlaveConstraintQueue concurrent_constraints;

    std::atomic<int> i(0);

    block_for_each(p_nodes_array, [&](Node::Pointer& p_slave_node) {

        // Skip nodes with IDs 1174 or 1015 950 1111
        if (p_slave_node->Id() == 1174 || p_slave_node->Id() == 1015) {
            return; // Skip the current iteration of the loop
        }   

        double local_radius = Radius;
        ResultNodesContainerType r_cloud_of_nodes(0);
        ResultNodesContainerType r_local_results(mMaxNumberOfNodes);

        while (r_cloud_of_nodes.size() < MinNumOfNeighNodes) {
            r_cloud_of_nodes.clear();
            SearchNodesInRadiusForNode(p_slave_node, local_radius, r_cloud_of_nodes, r_local_results);
            local_radius += Radius;
        }

        GetDofsAndCoordinatesForSlaveNode(p_slave_node, rVariableList, r_slave_dofs, r_slave_coordinates);
        FilterUniqueNodesForRBF(r_cloud_of_nodes, r_cloud_of_nodes_coordinates);

        ResultNodesContainerType r_cloud_of_master_nodes;
        r_cloud_of_master_nodes.reserve(r_cloud_of_nodes.size());
        for (auto& ghost_node : r_cloud_of_nodes) {
            auto it = GhostToOriginalMapping.find(ghost_node->Id());
            if (it != GhostToOriginalMapping.end()) {
                IndexType master_node_id = it->second.first; // Extracting the original/master node ID from the mapping
                Node::Pointer p_master_node = rComputingModelPart.pGetNode(master_node_id);
                r_cloud_of_master_nodes.push_back(p_master_node);
            } else {
                KRATOS_ERROR << "Ghost Node ID " << ghost_node->Id() << " not found in GhostToOriginalMapping!" << std::endl;
            }
        }

        GetDofsAndCoordinatesForCloudOfNodesFlattened(r_cloud_of_master_nodes, rVariableList, r_cloud_of_master_dofs, r_cloud_of_master_nodes_coordinates);

        RBFShapeFunctionsUtility::CalculateShapeFunctions(r_cloud_of_nodes_coordinates, r_slave_coordinates, r_n_container);

        Vector ListOfAngles;
        ExtractAnglesFromMapping(r_cloud_of_nodes, GhostToOriginalMapping, ListOfAngles);

        ObtainShapeMatrix(shape_matrix_x, shape_matrix_y, shape_matrix_z, r_n_container, ListOfAngles);
        constant_vector.resize(shape_matrix_x.size2());
        if (p_slave_node->Id() == 929){
            KRATOS_WATCH(r_cloud_of_nodes)
            KRATOS_WATCH(r_cloud_of_nodes_coordinates)
            KRATOS_WATCH(r_slave_coordinates)
            KRATOS_WATCH(r_cloud_of_master_nodes)
            KRATOS_WATCH(r_cloud_of_master_dofs)
            KRATOS_WATCH(r_n_container)
            KRATOS_WATCH(ListOfAngles)
            KRATOS_WATCH(shape_matrix_x)
            KRATOS_WATCH(shape_matrix_y)
            KRATOS_WATCH(shape_matrix_z)
        }

        IndexType it = i.fetch_add(1) * rVariableList.size();

        std::vector<Matrix> shape_matrices = {shape_matrix_x, shape_matrix_y, shape_matrix_z};

        for (std::size_t j = 0; j < rVariableList.size(); ++j) {
            if (j >= shape_matrices.size()) {
                KRATOS_ERROR << "Unexpected index in rVariableList. Only 3 DoFs (x, y, z) are supported." << std::endl;
            }
            concurrent_constraints.enqueue(r_clone_constraint.Create(prev_num_mscs + it + j + 1, r_cloud_of_master_dofs, r_slave_dofs[j], shape_matrices[j], constant_vector));
        }


    });

    ConstraintContainerType temp_constraints;
    ModelPart::MasterSlaveConstraintType::Pointer p_constraint;

    while (concurrent_constraints.try_dequeue(p_constraint))
        temp_constraints.push_back(p_constraint);

    rComputingModelPart.AddMasterSlaveConstraints(temp_constraints.begin(), temp_constraints.end());

    KRATOS_INFO("AssignRotationalMasterSlaveConstraintsToNodes") << "Build and Assign Rotational Master-Slave Constraints Time: " << build_and_assign_mscs.ElapsedSeconds() << std::endl;

    KRATOS_CATCH("");
}

void AssignMasterSlaveConstraintsToNeighboursUtility::ObtainShapeMatrix(
    Matrix& shape_matrix_x,
    Matrix& shape_matrix_y,
    Matrix& shape_matrix_z,
    Vector& rN_container,
    Vector& ListOfAngles
)
{
    std::vector<double> rotation_x_global;
    std::vector<double> rotation_y_global;
    std::vector<double> rotation_z_global;
    std::vector<double> rotation_x(3);
    std::vector<double> rotation_y(3);
    std::vector<double> rotation_z(3);
    for (int i=0; i < static_cast<int>(rN_container.size()); i++){
        //TODO: Any kind of rotation. x-axis, y-axis, z-axis. Make it a function.
        rotation_x[0] = rN_container[i]*cos(ListOfAngles[i]);
        rotation_x[1] = -rN_container[i]*sin(ListOfAngles[i]);
        rotation_x[2] = rN_container[i]*0.0;
        rotation_y[0] = rN_container[i]*sin(ListOfAngles[i]);
        rotation_y[1] = rN_container[i]*cos(ListOfAngles[i]);
        rotation_y[2] = rN_container[i]*0.0;
        rotation_z[0] = rN_container[i]*0.0;
        rotation_z[1] = rN_container[i]*0.0;
        rotation_z[2] = rN_container[i]*1.0;
        // noalias(row(shape_matrix,i)) = aux_vector_x;
        rotation_x_global.insert(rotation_x_global.end(),rotation_x.begin(),rotation_x.end());
        rotation_y_global.insert(rotation_y_global.end(),rotation_y.begin(),rotation_y.end());
        rotation_z_global.insert(rotation_z_global.end(),rotation_z.begin(),rotation_z.end());
    }
    shape_matrix_x.resize(1, rotation_x_global.size());
    shape_matrix_y.resize(1, rotation_y_global.size());
    shape_matrix_z.resize(1, rotation_z_global.size());
    for (int i=0; i<static_cast<int>(rotation_x_global.size());i++){
        shape_matrix_x(0,i) = rotation_x_global[i];
        shape_matrix_y(0,i) = rotation_y_global[i];
        shape_matrix_z(0,i) = rotation_z_global[i];
    }
}

// void AssignMasterSlaveConstraintsToNeighboursUtility::AssignRotationalMasterSlaveConstraintsToNodes(
//     NodesContainerType pSlaveNodes,
//     double const Radius,
//     ModelPart& rComputingModelPart,
//     const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList,
//     double const MinNumOfNeighNodes,
//     const std::unordered_map<IndexType, std::pair<IndexType, double>>& GhostToOriginalMapping
//     )
// {
//     KRATOS_TRY;

//     BuiltinTimer build_and_assign_mscs;

//     thread_local std::vector<DofPointerVectorType> r_cloud_of_dofs;
//     thread_local std::vector<DofPointerVectorType> r_slave_dofs;
//     thread_local Matrix r_cloud_of_nodes_coordinates;
//     thread_local array_1d<double, 3> r_slave_coordinates;
//     thread_local Vector r_n_container;
//     thread_local Vector constant_vector;
//     thread_local Matrix shape_matrix;

//     int prev_num_mscs = rComputingModelPart.NumberOfMasterSlaveConstraints();

//     KRATOS_WARNING_IF("AssignRotationalMasterSlaveConstraintsToNodes", prev_num_mscs > 0) << "Previous Master-Slave Constraints exist in the ModelPart. The new Constraints may interact with the existing ones." << std::endl;

//     auto& p_nodes_array = pSlaveNodes.GetContainer();

//     const auto& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");

//     using MasterSlaveConstraintQueue = moodycamel::ConcurrentQueue<ModelPart::MasterSlaveConstraintType::Pointer>;
//     MasterSlaveConstraintQueue concurrent_constraints;

//     std::atomic<int> i(0);

//     block_for_each(p_nodes_array, [&](Node::Pointer& p_slave_node) {

//         double local_radius = Radius;
//         ResultNodesContainerType r_cloud_of_nodes(0);
//         ResultNodesContainerType r_local_results(mMaxNumberOfNodes);

//         while (r_cloud_of_nodes.size() < MinNumOfNeighNodes) {
//             r_cloud_of_nodes.clear();
//             SearchNodesInRadiusForNode(p_slave_node, local_radius, r_cloud_of_nodes, r_local_results);
//             local_radius += Radius;
//         }

//         GetDofsAndCoordinatesForSlaveNode(p_slave_node, rVariableList, r_slave_dofs, r_slave_coordinates);
//         GetDofsAndCoordinatesForCloudOfNodes(r_cloud_of_nodes, rVariableList, r_cloud_of_dofs, r_cloud_of_nodes_coordinates);

//         RBFShapeFunctionsUtility::CalculateShapeFunctions(r_cloud_of_nodes_coordinates, r_slave_coordinates, r_n_container);

//         shape_matrix.resize(1, r_n_container.size());
//         noalias(row(shape_matrix, 0)) = r_n_container;
//         constant_vector.resize(r_n_container.size());

//         // Check if the current slave node is a ghost node and adjust for rotation if necessary
//         if (GhostToOriginalMapping.find(p_slave_node->Id()) != GhostToOriginalMapping.end()) {
//             auto mapping_data = GhostToOriginalMapping.at(p_slave_node->Id());
//             double rotation_angle = mapping_data.second;
//             AdjustConstraintForRotation(rotation_angle, shape_matrix);
//         }

//         IndexType it = i.fetch_add(1) * rVariableList.size();

//         for (std::size_t j = 0; j < rVariableList.size(); ++j) {
//             concurrent_constraints.enqueue(r_clone_constraint.Create(prev_num_mscs + it + j + 1, const_cast<DofPointerVectorType&>(r_cloud_of_dofs[j]), r_slave_dofs[j], shape_matrix, constant_vector));
//         }

//     });

//     ConstraintContainerType temp_constraints;
//     ModelPart::MasterSlaveConstraintType::Pointer p_constraint;

//     while (concurrent_constraints.try_dequeue(p_constraint))
//         temp_constraints.push_back(p_constraint);

//     rComputingModelPart.AddMasterSlaveConstraints(temp_constraints.begin(), temp_constraints.end());

//     KRATOS_INFO("AssignRotationalMasterSlaveConstraintsToNodes") << "Build and Assign Rotational Master-Slave Constraints Time: " << build_and_assign_mscs.ElapsedSeconds() << std::endl;

//     KRATOS_CATCH("");
// }

void AssignMasterSlaveConstraintsToNeighboursUtility::AdjustConstraintForRotation(double rotation_angle, Matrix& shape_matrix)
{
    KRATOS_TRY;

    double cos_theta = std::cos(rotation_angle * M_PI / 180.0); // Convert degrees to radians
    double sin_theta = std::sin(rotation_angle * M_PI / 180.0);

    for (std::size_t i = 0; i < shape_matrix.size1(); ++i) {
        double original_w1 = shape_matrix(i, 0);
        double original_w2 = shape_matrix(i, 1);

        shape_matrix(i, 0) = original_w1 * cos_theta - original_w2 * sin_theta;
        shape_matrix(i, 1) = original_w1 * sin_theta + original_w2 * cos_theta;
    }

    KRATOS_CATCH("");
}

void AssignMasterSlaveConstraintsToNeighboursUtility::ExtractAnglesFromMapping(
    const ResultNodesContainerType& Results,
    const std::unordered_map<IndexType, std::pair<IndexType, double>>& GhostToOriginalMapping,
    Vector& ListOfAngles
)
{
    ListOfAngles.resize(Results.size(), false);
    for(int i = 0; i < static_cast<int>(Results.size()); i++)
    {
        auto it = GhostToOriginalMapping.find(Results[i]->Id());
        if(it != GhostToOriginalMapping.end())
        {
            ListOfAngles[i] = it->second.second; // Extracting the angle from the mapping
        }
        else
        {
            KRATOS_ERROR << "Node ID " << Results[i]->Id() << " not found in GhostToOriginalMapping!" << std::endl;
        }
    }
}

}  // namespace Kratos.

void AssignMasterSlaveConstraintsToNeighboursUtility::InterpolateVelocityToNodes(
    NodesContainerType pDestinationNodes,
    double const Radius,
    const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList,
    double const MinNumOfNeighNodes
    )
{
    KRATOS_TRY;

    BuiltinTimer build_and_assign_mscs;

    // Declare thread-local storage (TLS) variables
    thread_local std::vector<DofPointerVectorType> r_cloud_of_dofs;           // The DOFs of the cloud of nodes
    thread_local std::vector<DofPointerVectorType> r_destination_dofs;              // The DOFs of the destination node
    thread_local Matrix r_cloud_of_nodes_coordinates;                          // The coordinates of the cloud of nodes
    thread_local array_1d<double, 3> r_destination_coordinates;                      // The coordinates of the destination node
    thread_local Vector r_n_container;                                         // Container for shape functions

    auto& p_nodes_array = pDestinationNodes.GetContainer();  // Get the container of destination nodes

    // Declare a counter variable outside the loop as std::atomic<int>
    std::atomic<int> i(0);

    block_for_each(p_nodes_array, [&](Node::Pointer& p_destination_node) {

        double local_radius = Radius;

        ResultNodesContainerType r_cloud_of_nodes(0);  // Container for the cloud of nodes
        ResultNodesContainerType r_local_results(mMaxNumberOfNodes);  // Container for node search results

        while (r_cloud_of_nodes.size() < MinNumOfNeighNodes) {
            r_cloud_of_nodes.clear();
            SearchNodesInRadiusForNode(p_destination_node, local_radius, r_cloud_of_nodes, r_local_results);  // Search for nodes in the radius
            local_radius += Radius;
        }

        // Get Dofs and Coordinates for the slave node and the cloud of nodes
        GetDofsAndCoordinatesForSlaveNode(p_destination_node, rVariableList, r_destination_dofs, r_destination_coordinates);
        GetDofsAndCoordinatesForCloudOfNodes(r_cloud_of_nodes, rVariableList, r_cloud_of_dofs, r_cloud_of_nodes_coordinates);

        // Calculate shape functions
        RBFShapeFunctionsUtility::CalculateShapeFunctions(r_cloud_of_nodes_coordinates, r_destination_coordinates, r_n_container);

        // IndexType it = i.fetch_add(1) * rVariableList.size(); // Atomically increment the counter and get the previous value

        for (std::size_t j = 0; j < rVariableList.size(); ++j) {
            double interpolated_variable_value = 0.0;
            for (std::size_t k = 0; k < r_cloud_of_nodes.size(); ++k){
                auto r_current_node = r_cloud_of_nodes[k];
                interpolated_variable_value += r_current_node->FastGetSolutionStepValue(rVariableList[j].get()) * r_n_container[k];
            }
            p_destination_node->FastGetSolutionStepValue(rVariableList[j].get(), interpolated_variable_value);
        }
    });
    KRATOS_CATCH("");
}  // namespace Kratos.
}