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


// Application includes
#include "compressible_navier_stokes_explicit_condition.h"


namespace Kratos {

/**
 * Returns the integration method for computation of midpoint magnitudes.
 * Computation of RHS integration method is chosen in the symbolic generator.
 */
template<>
GeometryData::IntegrationMethod CompressibleNavierStokesExplicitCondition<2,3>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

template<>
void CompressibleNavierStokesExplicitCondition<2,3>::CalculateRightHandSideInternal(
    BoundedVector<double, BlockSize * NumNodes>& rRightHandSideBoundedVector,
    const ProcessInfo& rCurrentProcessInfo)
{

}


template class CompressibleNavierStokesExplicitCondition<2,3>;
using CompressibleNavierStokesExplicitCondition2D2N = CompressibleNavierStokesExplicitCondition<2,3>;

}