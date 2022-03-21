//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//


// System includes


// External includes


// Project includes
#include "includes/define.h"


// Application includes
#include "compressible_navier_stokes_explicit_condition.h"


namespace Kratos {

/**
 * Returns the integration method for computation of midpoint magnitudes.
 * Computation of RHS integration method is chosen in the symbolic generator.
 */
template<>
GeometryData::IntegrationMethod CompressibleNavierStokesExplicitCondition<2,2>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

template<>
BoundedVector<double, 8> CompressibleNavierStokesExplicitCondition<2,2>::CalculateRightHandSideInternal(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BoundedVector<double, BlockSize*NumNodes> rRightHandSideBoundedVector = ZeroVector(BlockSize*NumNodes);

    const auto data = ConditionData(rCurrentProcessInfo);
    [[maybe_unused]] const auto& DN_DX = data.DN_DX;

    //substitute_rhs_2D_fluxes

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(NumNodes);

    return rRightHandSideBoundedVector;
    KRATOS_CATCH("")
}


template class CompressibleNavierStokesExplicitCondition<2,2>;
using CompressibleNavierStokesExplicitCondition2D2N = CompressibleNavierStokesExplicitCondition<2,2>;

}