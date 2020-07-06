//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Elisa Magliozzi
//

#include "custom_elements/compressible_biphase_navier_stokes_explicit.h"

namespace Kratos {

template<>
void CompressibleBiphaseNavierStokesExplicit<3>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    
    KRATOS_CATCH("")
}

template<>
void CompressibleBiphaseNavierStokesExplicit<3>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    
    KRATOS_CATCH("");
}

template<>
void CompressibleBiphaseNavierStokesExplicit<3>::ComputeGaussPointLHSContribution(BoundedMatrix<double,24,24>& lhs, const ElementDataStruct& data)
{
    
}

template<>
void CompressibleBiphaseNavierStokesExplicit<3>::ComputeGaussPointRHSContribution(array_1d<double,24>& rhs, const ElementDataStruct& data)
{


}

}