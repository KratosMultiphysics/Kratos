//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"

// Application includes
#include "two_fluid_navier_stokes_wall_condition.h"


namespace Kratos
{

template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
int TwoFluidNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = BaseType::Check(rCurrentProcessInfo);
    if (check != 0) {
        return check;
    } else {
        // Checks on nodes
        const auto& r_geom = BaseType::GetGeometry();
        for (const auto& r_node : r_geom) {
            if(r_node.SolutionStepsDataHas(DYNAMIC_VISCOSITY) == false){
                KRATOS_ERROR << "missing DYNAMIC_VISCOSITY variable on solution step data for node " << r_node.Id();
            }
        }

        return check;
    }

    KRATOS_CATCH("");
}

// template<unsigned int TDim, unsigned int TNumNodes, class TWallModel>
// void TwoFluidNavierStokesWallCondition<TDim,TNumNodes,TWallModel>::ComputeGaussPointNavierSlipRHSContribution(
//     array_1d<double,LocalSize>& rRightHandSideVector,
//     const ConditionDataStruct& rDataStruct )
// {
//     KRATOS_TRY

//     const auto& r_geom = this->GetGeometry();
//     const double zero_tol = 1.0e-12;
//     array_1d<double,3> nodal_normal;
//     BoundedMatrix<double, TNumNodes, TNumNodes> nodal_projection_matrix;
//     for (unsigned int i_node = 0; i_node < TNumNodes; i_node++){
//         // Finding the nodal tangential projection matrix (I - n x n)
//         nodal_normal = r_geom[i_node].FastGetSolutionStepValue(NORMAL);
//         const double norm = norm_2(nodal_normal);
//         if (norm > zero_tol) {
//             nodal_normal /= norm;
//         }
//         FluidElementUtilities<3>::SetTangentialProjectionMatrix(nodal_normal, nodal_projection_matrix);

//         // Finding the coefficient to relate velocity to drag
//         const double viscosity = r_geom[i_node].GetSolutionStepValue(DYNAMIC_VISCOSITY);
//         const double navier_slip_length = r_geom[i_node].GetValue(SLIP_LENGTH);
//         KRATOS_ERROR_IF_NOT( navier_slip_length > 0.0 ) << "Negative or zero slip length was defined" << std::endl;
//         const double nodal_beta = viscosity / navier_slip_length;

//         const array_1d<double, TNumNodes> N = rDataStruct.N;
//         const double wGauss = rDataStruct.wGauss;
//         Vector interpolated_velocity = ZeroVector(TNumNodes);
//         for( unsigned int comp = 0; comp < TNumNodes; comp++){
//             for (unsigned int i = 0; i < TNumNodes; i++){
//                 // necessary because VELOCITY with 3 entries even in 2D case
//                 interpolated_velocity[i] -= N[comp] * r_geom[comp].FastGetSolutionStepValue(VELOCITY)[i];
//             }
//         }
//         // Application of the nodal projection matrix
//         const array_1d<double,TNumNodes> nodal_entry_rhs = prod( nodal_projection_matrix, (wGauss * N[i_node] * nodal_beta * interpolated_velocity) );
//         for (unsigned int entry = 0; entry < TNumNodes; entry++){
//             rRightHandSideVector( i_node*(TNumNodes+1) + entry ) += nodal_entry_rhs[entry];
//         }
//     }

//     KRATOS_CATCH("")
// }


// template<unsigned int TDim, unsigned int TNumNodes, class TWallModel>
// void TwoFluidNavierStokesWallCondition<TDim,TNumNodes,TWallModel>::ComputeGaussPointNavierSlipLHSContribution(
//     BoundedMatrix<double,LocalSize,LocalSize>& rLeftHandSideMatrix,
//     const ConditionDataStruct& rDataStruct )
// {
//     KRATOS_TRY

//     const GeometryType& r_geom = this->GetGeometry();

//     array_1d<double, TNumNodes> N = rDataStruct.N;
//     const double wGauss = rDataStruct.wGauss;

//     for(unsigned int inode = 0; inode < TNumNodes; inode++){

//         // finding the nodal projection matrix nodal_projection_matrix = ( [I] - (na)(na) )
//         BoundedMatrix<double, TNumNodes, TNumNodes> nodal_projection_matrix;
//         array_1d<double,3> nodal_normal = r_geom[inode].FastGetSolutionStepValue(NORMAL);
//         double sum_of_squares = 0.0;
//         for (unsigned int j = 0; j < 3; j++){
//             sum_of_squares += nodal_normal[j] * nodal_normal[j];
//         }
//         nodal_normal /= sqrt(sum_of_squares);
//         FluidElementUtilities<3>::SetTangentialProjectionMatrix( nodal_normal, nodal_projection_matrix );

//         // finding the coefficient to relate velocity to drag
//         const double viscosity = r_geom[inode].GetSolutionStepValue(DYNAMIC_VISCOSITY);
//         const double navier_slip_length = r_geom[inode].GetValue(SLIP_LENGTH);
//         KRATOS_ERROR_IF_NOT( navier_slip_length > 0.0 ) << "Negative or zero slip length was defined" << std::endl;
//         const double nodal_beta = viscosity / navier_slip_length;

//         for(unsigned int jnode = 0; jnode < TNumNodes; jnode++){

//             const BoundedMatrix<double, TNumNodes, TNumNodes> nodal_lhs_contribution = wGauss * nodal_beta * N[inode] * N[jnode] * nodal_projection_matrix;

//             for( unsigned int i = 0; i < TNumNodes; i++){
//                 for( unsigned int j = 0; j < TNumNodes; j++){

//                     const unsigned int istart = inode * (TNumNodes+1);
//                     const unsigned int jstart = jnode * (TNumNodes+1);
//                     rLeftHandSideMatrix(istart + i, jstart + j) += nodal_lhs_contribution(i,j);
//                 }
//             }
//         }
//     }

//     KRATOS_CATCH("")
// }


template class TwoFluidNavierStokesWallCondition<2,2>;
template class TwoFluidNavierStokesWallCondition<3,3>;

} // namespace Kratos
