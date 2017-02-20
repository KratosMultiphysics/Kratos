//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Original author:  Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "partitioned_fsi_utilities.hpp"

namespace Kratos
{

// Fluid interface domain velocity imposition
// template<class TSpace, unsigned int TDim>
// void PartitionedFSIUtilities<TSpace, TDim>::SetFluidInterfaceVelocity(VectorType& rFluidInterfaceVelocity)
// {
//     unsigned int i = 0;
//
//     if (TDim == 2)
//     {
//         for (ModelPart::NodesContainerType::iterator it_node = mrFluidInterfaceModelPart.NodesBegin(); it_node<mrFluidInterfaceModelPart.NodesEnd(); it_node++)
//         {
//             array_1d<double,3> aux_velocity;
//             aux_velocity[0] = rFluidInterfaceVelocity[i];
//             aux_velocity[1] = rFluidInterfaceVelocity[i+1];
//             aux_velocity[2] = 0.0;
//
//             it_node->Fix(VELOCITY_X);
//             it_node->Fix(VELOCITY_Y);
//             it_node->Fix(VELOCITY_Z);
//
//             it_node->FastGetSolutionStepValue(VELOCITY) = aux_velocity;
//             i += 2;
//         }
//     }
//     else
//     {
//         for (ModelPart::NodesContainerType::iterator it_node = mrFluidInterfaceModelPart.NodesBegin(); it_node<mrFluidInterfaceModelPart.NodesEnd(); it_node++)
//         {
//             array_1d<double,3> aux_velocity;
//             aux_velocity[0] = rFluidInterfaceVelocity[i];
//             aux_velocity[1] = rFluidInterfaceVelocity[i+1];
//             aux_velocity[2] = rFluidInterfaceVelocity[i+2];
//
//             it_node->Fix(VELOCITY_X);
//             it_node->Fix(VELOCITY_Y);
//             it_node->Fix(VELOCITY_Z);
//
//             it_node->FastGetSolutionStepValue(VELOCITY) = aux_velocity;
//             i += 3;
//         }
//     }
// }
//
// // Fluid interface domain nodal force imposition
// template<class TSpace, unsigned int TDim>
// void PartitionedFSIUtilities<TSpace, TDim>::SetFluidInterfaceNodalForce(VectorType& rFluidInterfaceNodalForce)
// {
//     unsigned int i = 0;
//
//     if (TDim == 2)
//     {
//         for (ModelPart::NodesContainerType::iterator it_node = mrFluidInterfaceModelPart.NodesBegin(); it_node<mrFluidInterfaceModelPart.NodesEnd(); it_node++)
//         {
//             array_1d<double,3> aux_force;
//             aux_force[0] = rFluidInterfaceNodalForce[i];
//             aux_force[1] = rFluidInterfaceNodalForce[i+1];
//             aux_force[2] = 0.0;
//
//             it_node->Fix(FORCE_X);
//             it_node->Fix(FORCE_Y);
//             it_node->Fix(FORCE_Z);
//
//             it_node->FastGetSolutionStepValue(FORCE) = aux_force;
//             i += 2;
//         }
//     }
//     else
//     {
//         for (ModelPart::NodesContainerType::iterator it_node = mrFluidInterfaceModelPart.NodesBegin(); it_node<mrFluidInterfaceModelPart.NodesEnd(); it_node++)
//         {
//             array_1d<double,3> aux_force;
//             aux_force[0] = rFluidInterfaceNodalForce[i];
//             aux_force[1] = rFluidInterfaceNodalForce[i+1];
//             aux_force[2] = rFluidInterfaceNodalForce[i+2];
//
//             it_node->Fix(FORCE_X);
//             it_node->Fix(FORCE_Y);
//             it_node->Fix(FORCE_Z);
//
//             it_node->FastGetSolutionStepValue(FORCE) = aux_force;
//             i += 3;
//         }
//     }
// }
//
// // Fluid interface domain velocity residual computation
// template<class TSpace, unsigned int TDim>
// typename PartitionedFSIUtilities<TSpace, TDim>::VectorType PartitionedFSIUtilities<TSpace, TDim>::ComputeFluidInterfaceVelocityResidual()
// {
//
//     // Compute the residual
//     unsigned int residual_size = (this->GetFluidInterfaceProblemSize)*TDim;
//     double err_j = 0.0;
//     double total_weight = 0.0;
//     double fluid_interface_residual_norm = 0.0;
//     VectorType fluid_interface_residual(residual_size);
//
//     unsigned int i = 0;
//
//     for (ModelPart::NodesContainerType::iterator it_node = mrFluidInterfaceModelPart.NodesBegin(); it_node<mrFluidInterfaceModelPart.NodesEnd(); it_node++)
//     {
//         double err_node = 0.0;
//         array_1d<double,3> velocity_fluid = it_node->FastGetSolutionStepValue(VELOCITY);
//         array_1d<double,3> velocity_fluid_projected = it_node->FastGetSolutionStepValue(VECTOR_PROJECTED);
//
//         for (int j=0; j<TDim; j++)
//         {
//             err_j = velocity_fluid[j] - velocity_fluid_projected[j];
//             fluid_interface_residual[i+j] = err_j;
//             err_node += std::pow(err_j,2);
//         }
//
//         double weight = it_node->GetValue(NODAL_AREA);
//         total_weight += weight;
//         fluid_interface_residual_norm += weight*err_j;
//
//         i += TDim;
//     }
//
//     mrFluidInterfaceModelPart.GetCommunicator().SumAll(total_weight);
//     mrFluidInterfaceModelPart.GetCommunicator().SumAll(fluid_interface_residual_norm);
//
//     // Store its weighted L2 norm in the fluid process info
//     mrFluidInterfaceModelPart.GetProcessInfo().GetValue(FSI_INTERFACE_RESIDUAL_NORM) = fluid_interface_residual_norm/total_weight;
//
//     return fluid_interface_residual;
//
// }

} // Namespace Kratos.
