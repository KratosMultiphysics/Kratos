//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Sebastian Ares de Parga Regalado
//
// This file is partly copied from
// "DEMApplication/custom_utilities/omp_dem_search.h"
// and modified to create only once the binary tree

// System includes
#include <string>
#include <iostream>

// Project includes
#include "utilities/assign_mpcs_to_neighbours_utility.h"

// Configures

// External includes

namespace Kratos {
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/**
 * @class AssignMPCsToNeighboursUtility
 * @ingroup Kratos Core
 * @brief Assing MPCs to Neighbouring Nodes
 * @details This class provides a method to search for neighbouring nodes of one node
 * within a given radius and assign a multipoint constraint (MPC) using radial basis functions.
 * @author Sebastian Ares de Parga Regalado
 */

/// Default constructor.
AssignMPCsToNeighboursUtility::AssignMPCsToNeighboursUtility(NodesContainerType& rStructureNodes){
    KRATOS_TRY;
    NodesContainerType::ContainerType& nodes_model_part = rStructureNodes.GetContainer();
    mpBins = Kratos::make_unique<NodeBinsType>(nodes_model_part.begin(), nodes_model_part.end());
    mMaxNumberOfNodes = rStructureNodes.size();
    KRATOS_CATCH(""); 
} 

// Search for the neighbouring nodes (in rStructureNodes) of each rNode on a given array of rNodes
void AssignMPCsToNeighboursUtility::SearchNodesInRadiusForNodes(
    NodesContainerType const& rNodes,
    const double Radius,
    VectorResultNodesContainerType& rResults)
{
    KRATOS_TRY;
    NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(rNodes.GetContainer());

    rResults.resize(nodes_array.size());

    #pragma omp parallel
    {
        ResultNodesContainerType  localResults(mMaxNumberOfNodes);
        std::size_t               NumberOfResults = 0;

        #pragma omp for
        for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
        {
            ResultNodesContainerType::iterator ResultsPointer    = localResults.begin();

            NumberOfResults = mpBins->SearchInRadius(*(nodes_array[i]),Radius,ResultsPointer,mMaxNumberOfNodes);

            rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);
        }
    }
    KRATOS_CATCH("");
}


//Search for the neighbouring nodes (in rStructureNodes) of the given rNode
void AssignMPCsToNeighboursUtility::SearchNodesInRadiusForNode(
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


void AssignMPCsToNeighboursUtility::GetDofsAndCoordinatesForNodes(
NodeType::Pointer pNode,
const Variable<array_1d<double, 3>>& rVariable,
DofPointerVectorType& rCloud_of_dofs_x,
DofPointerVectorType& rCloud_of_dofs_y,
DofPointerVectorType& rCloud_of_dofs_z,
array_1d<double,3>& rSlave_coordinates
){
KRATOS_TRY;
if (rVariable==VELOCITY){
    rCloud_of_dofs_x.push_back(pNode->pGetDof(VELOCITY_X));
    rCloud_of_dofs_y.push_back(pNode->pGetDof(VELOCITY_Y));
    rCloud_of_dofs_z.push_back(pNode->pGetDof(VELOCITY_Z));
}
else if (rVariable==DISPLACEMENT){
    // NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(pNodes.GetContainer());
    rCloud_of_dofs_x.push_back(pNode->pGetDof(DISPLACEMENT_X));
    rCloud_of_dofs_y.push_back(pNode->pGetDof(DISPLACEMENT_Y));
    rCloud_of_dofs_z.push_back(pNode->pGetDof(DISPLACEMENT_Z));
}
else{
    KRATOS_ERROR << "This function is not yet implemented for variable " << rVariable << std::endl;
}
rSlave_coordinates = pNode->Coordinates();
KRATOS_CATCH("");
}

// Get Dofs and coordinates arrays for a given variable array. (For nodes)
void AssignMPCsToNeighboursUtility::GetDofsAndCoordinatesForNodes(
ResultNodesContainerType nodes_array,
const Variable<array_1d<double, 3>>& rVariable,
DofPointerVectorType& rCloud_of_dofs_x,
DofPointerVectorType& rCloud_of_dofs_y,
DofPointerVectorType& rCloud_of_dofs_z,
Matrix& rCloud_of_nodes_coordinates
)
{
KRATOS_TRY;
rCloud_of_dofs_x.resize(nodes_array.size());
rCloud_of_dofs_y.resize(nodes_array.size());
rCloud_of_dofs_z.resize(nodes_array.size());
rCloud_of_nodes_coordinates.resize(nodes_array.size(),3);
if (rVariable==VELOCITY){
    for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
    {
    rCloud_of_dofs_x[i] = nodes_array[i]->pGetDof(VELOCITY_X);
    rCloud_of_dofs_y[i] = nodes_array[i]->pGetDof(VELOCITY_Y);
    rCloud_of_dofs_z[i] = nodes_array[i]->pGetDof(VELOCITY_Z);
    noalias(row(rCloud_of_nodes_coordinates,i)) = nodes_array[i]->Coordinates();
    }
}
else if (rVariable==DISPLACEMENT){
    for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
    {
    rCloud_of_dofs_x[i] = nodes_array[i]->pGetDof(DISPLACEMENT_X);
    rCloud_of_dofs_y[i] = nodes_array[i]->pGetDof(DISPLACEMENT_Y);
    rCloud_of_dofs_z[i] = nodes_array[i]->pGetDof(DISPLACEMENT_Z);
    noalias(row(rCloud_of_nodes_coordinates,i)) = nodes_array[i]->Coordinates();
    }
}
else{
    KRATOS_ERROR << "This function is not yet implemented for variable " << rVariable << std::endl;
}
KRATOS_CATCH("");
}

void AssignMPCsToNeighboursUtility::GetDofsAndCoordinatesForNodes(
NodeType::Pointer pNode,
const Variable<double>& rVariable,
DofPointerVectorType& rCloud_of_dofs,
array_1d<double,3>& rSlave_coordinates
){
KRATOS_TRY;
rCloud_of_dofs.push_back(pNode->pGetDof(rVariable));
rSlave_coordinates = pNode->Coordinates();
KRATOS_CATCH("");
}

// Get Dofs and Coordinates arrays for a given variable double. (For nodes)
void AssignMPCsToNeighboursUtility::GetDofsAndCoordinatesForNodes(
ResultNodesContainerType nodes_array,
const Variable<double>& rVariable,
DofPointerVectorType& rCloud_of_dofs,
Matrix& rCloud_of_nodes_coordinates
)
{
KRATOS_TRY;
rCloud_of_dofs.resize(nodes_array.size());
rCloud_of_nodes_coordinates.resize(nodes_array.size(),3);
for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
{
    rCloud_of_dofs[i] = nodes_array[i]->pGetDof(rVariable);
    noalias(row(rCloud_of_nodes_coordinates,i)) = nodes_array[i]->Coordinates();
}
KRATOS_CATCH("");
}

void AssignMPCsToNeighboursUtility::AssignRotationToNodes(
NodesContainerType pNodes,
Matrix RotationMatrix
)
{
//TODO: Do it in parallel
KRATOS_TRY;
NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(pNodes.GetContainer());

// Assign rotation to given nodes
for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
{
    auto coordinate_vector = nodes_array[i]->GetInitialPosition().Coordinates();
    Vector rotated_position(coordinate_vector.size());
    noalias(rotated_position) = prod(RotationMatrix,coordinate_vector);
    nodes_array[i]->X() = rotated_position[0];
    nodes_array[i]->Y() = rotated_position[1];
    nodes_array[i]->Z() = rotated_position[2];
}
KRATOS_CATCH("");
}


void AssignMPCsToNeighboursUtility::AssignMPCsToNodes(
    NodesContainerType pNodes,
    double const Radius,
    ModelPart& rComputingModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    double const MinNumOfNeighNodes
    )
    //TODO: Do it in parallel and allow rVariable to be a double variable i.e. TEMPERATURE (Make a template)
{
    KRATOS_TRY;
    // #pragma omp parallel
    {
    BuiltinTimer build_and_assign_mpcs;

    NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(pNodes.GetContainer());

    ModelPart::MasterSlaveConstraintType const& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");
    ConstraintContainerType all_constraints;
    all_constraints.reserve(nodes_array.size()*3);
    int prev_num_mpcs = rComputingModelPart.NumberOfMasterSlaveConstraints();

    // #pragma omp for
        for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
        {
        //Set slave nodes' slip condition to false
        // nodes_array[i]->Set(SLIP,false);

        // Search cloud of nodes for slave node
        ResultNodesContainerType  Results;
        double localRadius = Radius;
        while (Results.size()<MinNumOfNeighNodes){
            Results.clear();
            SearchNodesInRadiusForNode(nodes_array[i], localRadius, Results);
            localRadius += Radius;
        }

        // Get Dofs and Coordinates
        DofPointerVectorType rCloud_of_dofs_x,rCloud_of_dofs_y,rCloud_of_dofs_z,rSlave_dof_x,rSlave_dof_y,rSlave_dof_z;
        Matrix rCloud_of_nodes_coordinates;
        array_1d<double,3> rSlave_coordinates;
        GetDofsAndCoordinatesForNodes(nodes_array[i], rVariable, rSlave_dof_x,rSlave_dof_y,rSlave_dof_z, rSlave_coordinates);
        GetDofsAndCoordinatesForNodes(Results, rVariable, rCloud_of_dofs_x, rCloud_of_dofs_y, rCloud_of_dofs_z, rCloud_of_nodes_coordinates);

        // Calculate shape functions
        Vector rN_container;
        RBFShapeFunctionsUtility::CalculateShapeFunctions(rCloud_of_nodes_coordinates,rSlave_coordinates,rN_container);

        //Create MPCs
        Matrix shape_matrix(1,rN_container.size());
        noalias(row(shape_matrix,0)) = rN_container;// Shape functions matrix
        const Vector constant_vector = ZeroVector(rN_container.size());
        IndexType it = i*3;
        all_constraints.push_back(r_clone_constraint.Create(prev_num_mpcs + it + 1,rCloud_of_dofs_x,rSlave_dof_x,shape_matrix,constant_vector));
        all_constraints.push_back(r_clone_constraint.Create(prev_num_mpcs + it + 2,rCloud_of_dofs_y,rSlave_dof_y,shape_matrix,constant_vector));
        all_constraints.push_back(r_clone_constraint.Create(prev_num_mpcs + it + 3,rCloud_of_dofs_z,rSlave_dof_z,shape_matrix,constant_vector));
        }
        rComputingModelPart.AddMasterSlaveConstraints(all_constraints.begin(),all_constraints.end());
        KRATOS_INFO("AssignMPCsToNeighboursUtility") << "Build and Assign MPCs Time: " << build_and_assign_mpcs.ElapsedSeconds() << std::endl;
    }
    KRATOS_CATCH("");
}

void AssignMPCsToNeighboursUtility::AssignMPCsToNodes(
    NodesContainerType pNodes,
    double const Radius,
    ModelPart& rComputingModelPart,
    const Variable<double>& rVariable,
    double const MinNumOfNeighNodes
    )
    //TODO: Do it in parallel and allow rVariable to be a double variable i.e. TEMPERATURE (Make a template)
{
    KRATOS_TRY;
    // #pragma omp parallel
    {
    BuiltinTimer build_and_assign_mpcs;

    NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(pNodes.GetContainer());

    ModelPart::MasterSlaveConstraintType const& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");
    ConstraintContainerType all_constraints;
    all_constraints.reserve(nodes_array.size()*MinNumOfNeighNodes);
    int prev_num_mpcs = rComputingModelPart.NumberOfMasterSlaveConstraints();

    // #pragma omp for
        for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
        {
        //Set slave nodes' slip condition to false
        // nodes_array[i]->Set(SLIP,false);

        // Search cloud of nodes for slave node
        ResultNodesContainerType  Results;
        double localRadius = Radius;
        while (Results.size()<MinNumOfNeighNodes){
            Results.clear();
            SearchNodesInRadiusForNode(nodes_array[i], localRadius, Results);
            localRadius += Radius;
        }

        // Get Dofs and Coordinates
        DofPointerVectorType rCloud_of_dofs,rSlave_dof;
        Matrix rCloud_of_nodes_coordinates;
        array_1d<double,3> rSlave_coordinates;
        GetDofsAndCoordinatesForNodes(nodes_array[i], rVariable, rSlave_dof, rSlave_coordinates);
        GetDofsAndCoordinatesForNodes(Results, rVariable, rCloud_of_dofs, rCloud_of_nodes_coordinates);

        // Calculate shape functions
        Vector rN_container;
        RBFShapeFunctionsUtility::CalculateShapeFunctions(rCloud_of_nodes_coordinates,rSlave_coordinates,rN_container);

        //Create MPCs
        Matrix shape_matrix(1,rN_container.size());
        noalias(row(shape_matrix,0)) = rN_container;// Shape functions matrix
        const Vector constant_vector = ZeroVector(rN_container.size());
        IndexType it = i;
        all_constraints.push_back(r_clone_constraint.Create(prev_num_mpcs + it + 1,rCloud_of_dofs,rSlave_dof,shape_matrix,constant_vector));
        }
        rComputingModelPart.AddMasterSlaveConstraints(all_constraints.begin(),all_constraints.end());
        KRATOS_INFO("AssignMPCsToNeighboursUtility") << "Build and Assign MPCs Time: " << build_and_assign_mpcs.ElapsedSeconds() << std::endl;
    }
    KRATOS_CATCH("");
}

}  // namespace Kratos.