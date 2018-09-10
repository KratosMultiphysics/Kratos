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

#include "embedded_ausas_navier_stokes_wall_condition.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see EmbeddedAusasNavierStokesWallCondition::EquationIdVector
 */
template <>
void EmbeddedAusasNavierStokesWallCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,
                                                                   ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes = 2;
    const unsigned int LocalSize = 6;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

/**
 * @see EmbeddedAusasNavierStokesWallCondition::EquationIdVector
 */
template <>
void EmbeddedAusasNavierStokesWallCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,
                                                                   ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

/**
 * @see EmbeddedAusasNavierStokesWallCondition::GetDofList
 */
template <>
void EmbeddedAusasNavierStokesWallCondition<2,2>::GetDofList(DofsVectorType& rElementalDofList,
                                                             ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 2;
    const SizeType LocalSize = 6;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

/**
 * @see EmbeddedAusasNavierStokesWallCondition::GetDofList
 */
template <>
void EmbeddedAusasNavierStokesWallCondition<3,3>::GetDofList(DofsVectorType& rElementalDofList,
                                                             ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

/// Computes the Gauss pt. LHS contribution
/**
* @param lhs_gauss reference to the local LHS matrix
* @param data Gauss pt. data structure
*/
template<unsigned int TDim, unsigned int TNumNodes>
void EmbeddedAusasNavierStokesWallCondition<TDim,TNumNodes>::ComputeGaussPointLHSContribution(BoundedMatrix<double, TNumNodes*(TDim+1), TNumNodes*(TDim+1)>& lhs_gauss,
                                                                                              const ConditionDataStruct& rData)
{
    const unsigned int LocalSize = TDim+1;
    noalias(lhs_gauss) = ZeroMatrix(TNumNodes*LocalSize, TNumNodes*LocalSize);

    // LHS boundary term coming from the integration by parts of the mass conservation equation
    for (unsigned int i=0; i<TNumNodes; ++i)
    {
        const unsigned int row = i*LocalSize + TDim;

        for (unsigned int j=0; j<TNumNodes; ++j)
        {
            for (unsigned int dim=0; dim<TDim; ++dim)
            {
                const unsigned int col = j*LocalSize + dim;
                lhs_gauss(row, col) = rData.wGauss * rData.N[i] * rData.N[j] * rData.Normal[dim];
            }
        }
    }

}

/// Computes the Gauss pt. RHS contribution
/**
* @param rhs_gauss reference to the local RHS vector
* @param data Gauss pt. data structure
*/
template<unsigned int TDim, unsigned int TNumNodes>
void EmbeddedAusasNavierStokesWallCondition<TDim,TNumNodes>::ComputeGaussPointRHSContribution(array_1d<double, TNumNodes*(TDim+1)>& rhs_gauss,
                                                                                              const ConditionDataStruct& rData)
{
    // Initialize the local RHS
    const unsigned int LocalSize = TDim+1;
    noalias(rhs_gauss) = ZeroVector(TNumNodes*LocalSize);

    // Gauss pt. Neumann BC contribution
    this->ComputeRHSNeumannContribution(rhs_gauss, rData);

    // Gauss pt. outlet inflow prevention contribution
    if (this->Is(OUTLET))
    {
        this->ComputeRHSOutletInflowContribution(rhs_gauss, rData);
    }

    // RHS boundary term coming from the integration by parts of the mass conservation equation
    for (unsigned int i=0; i<TNumNodes; ++i)
    {
        const unsigned int row = i*LocalSize + TDim;

        for (unsigned int j=0; j<TNumNodes; ++j)
        {
            for (unsigned int dim=0; dim<TDim; ++dim)
            {
                rhs_gauss[row] -= rData.wGauss * rData.N[i] * rData.N[j] * rData.Normal[dim] * rData.v(j,dim);
            }
        }
    }
}

/// Computes the condition RHS Neumann BC contribution
/**
* @param rhs_gauss reference to the local RHS vector
* @param data Gauss pt. data structure
*/
template<unsigned int TDim, unsigned int TNumNodes>
void EmbeddedAusasNavierStokesWallCondition<TDim,TNumNodes>::ComputeRHSNeumannContribution(array_1d<double, TNumNodes*(TDim+1)>& rhs_gauss,
                                                                                           const ConditionDataStruct& rData)
{
    const unsigned int LocalSize = TDim+1;
    const GeometryType& rGeom = this->GetGeometry();

    // Add Neumann BC contribution
    for (unsigned int i=0; i<TNumNodes; ++i)
    {
        const double pext = rGeom[i].FastGetSolutionStepValue(EXTERNAL_PRESSURE);

        for (unsigned int j=0; j<TNumNodes; ++j)
        {
            unsigned int row = j*LocalSize;
            for (unsigned int d=0; d<TDim; ++d)
            {
                rhs_gauss[row + d] -= rData.wGauss * rData.N[j] * rData.N[i] * pext * rData.Normal[d];
            }
        }
    }
}

/// Computes the condition RHS outlet inflow prevention contribution
/**
* @param rhs_gauss reference to the local RHS vector
* @param rData Gauss pt. data structure
*/
template<unsigned int TDim, unsigned int TNumNodes>
void EmbeddedAusasNavierStokesWallCondition<TDim,TNumNodes>::ComputeRHSOutletInflowContribution(
    array_1d<double, TNumNodes*(TDim+1)>& rhs_gauss,
    const ConditionDataStruct& rData)
{
    const unsigned int LocalSize = TDim+1;
    const GeometryType& r_geom = this->GetGeometry();

    // Compute Gauss pt. density, velocity norm and velocity projection
    double rho_gauss = 0.0;
    array_1d<double, 3> v_gauss = ZeroVector(3);
    for (unsigned int i=0; i<TNumNodes; ++i) {
        const double rho = r_geom[i].FastGetSolutionStepValue(DENSITY);
        const array_1d<double, 3>& r_vel = r_geom[i].FastGetSolutionStepValue(VELOCITY);
        rho_gauss += rData.N[i]*rho;
        v_gauss += rData.N[i]*r_vel;
    }

    const double v_gauss_proj = inner_prod(v_gauss, rData.Normal);
    const double v_gauss_squared_norm = std::pow(v_gauss[0],2) + std::pow(v_gauss[1],2) + std::pow(v_gauss[2],2);

    // Add outlet inflow prevention contribution
    const double delta = rData.delta;
    const double U_0 = rData.charVel;
    const double S_0 = 0.5*(1-tanh(v_gauss_proj/(U_0*delta)));

    for (unsigned int i=0; i<TNumNodes; ++i) {
        unsigned int row = i*LocalSize;
        for (unsigned int d=0; d<TDim; ++d) {
            rhs_gauss[row+d] += rData.wGauss*rData.N[i]*0.5*rho_gauss*v_gauss_squared_norm*S_0*rData.Normal[d];
        }
    }
}

template class EmbeddedAusasNavierStokesWallCondition<2,2>;
template class EmbeddedAusasNavierStokesWallCondition<3,3>;

} // namespace Kratos
