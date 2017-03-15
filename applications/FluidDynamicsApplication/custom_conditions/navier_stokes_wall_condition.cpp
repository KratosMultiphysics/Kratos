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

#include "navier_stokes_wall_condition.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see NavierStokesWallCondition::EquationIdVector
 */
template <>
void NavierStokesWallCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,
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
 * @see NavierStokesWallCondition::EquationIdVector
 */
template <>
void NavierStokesWallCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,
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
 * @see NavierStokesWallCondition::GetDofList
 */
template <>
void NavierStokesWallCondition<2,2>::GetDofList(DofsVectorType& rElementalDofList,
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
 * @see NavierStokesWallCondition::GetDofList
 */
template <>
void NavierStokesWallCondition<3,3>::GetDofList(DofsVectorType& rElementalDofList,
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




// template<unsigned int TDim, unsigned int TNumNodes>
// void NavierStokesWallCondition<TDim,TNumNodes>::CalculateLocalVelocityContribution(MatrixType &rDampMatrix,
//                                                                                    VectorType &rRightHandSideVector,
//                                                                                    ProcessInfo &rCurrentProcessInfo)
// {
//     // Initialize local contributions
//     const SizeType LocalSize = (TDim + 1) * TNumNodes;
//
//     if (rDampMatrix.size1() != LocalSize)
//         rDampMatrix.resize(LocalSize,LocalSize);
//     if (rRightHandSideVector.size() != LocalSize)
//         rRightHandSideVector.resize(LocalSize);
//
//     noalias(rDampMatrix) = ZeroMatrix(LocalSize,LocalSize);
//     noalias(rRightHandSideVector) = ZeroVector(LocalSize);
//
//     this->ApplyNeumannCondition(rDampMatrix,rRightHandSideVector);
//
//     this->ApplyWallLaw(rDampMatrix,rRightHandSideVector,rCurrentProcessInfo);
// }

/// Returns a list of the element's Dofs
/**
* @param ElementalDofList the list of DOFs
* @param Ngauss Gauss point shape function values
*/
template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeGaussPointLHSContribution(bounded_matrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& lhs_gauss,
                                                                                 const ConditionDataStruct& data)
{
    const unsigned int LocalSize = TDim+1;
    noalias(lhs_gauss) = ZeroMatrix(TNumNodes*LocalSize,TNumNodes*LocalSize);
}

/// Returns a list of the element's Dofs
/**
* @param ElementalDofList the list of DOFs
* @param Ngauss  Gauss point shape function values
*/
template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs_gauss,
                                                                                 const ConditionDataStruct& data)
{
    const unsigned int LocalSize = TDim+1;
    const GeometryType& rGeom = this->GetGeometry();

    // Initialize the local RHS
    rhs_gauss = ZeroVector(TNumNodes*LocalSize);

    // Gauss pt. Neumann BC contribution
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        //unsigned int row = i*LocalSize;
        const NodeType& rConstNode = rGeom[i];
        const double pext = rConstNode.FastGetSolutionStepValue(EXTERNAL_PRESSURE);

        for (unsigned int j = 0; j < TNumNodes; j++)
        {
            unsigned int row = j*LocalSize;
            for (unsigned int d = 0; d < TDim;d++)
            {
                rhs_gauss[row+d] -= data.wGauss*data.N[j]*data.N[i]*pext*data.Normal[d];
            }
        }
    }
}

// protected funcions

template <>
void NavierStokesWallCondition<2,2>::CalculateNormal(array_1d<double,3>& An)
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = - (pGeometry[1].X() - pGeometry[0].X());
    An[2] =    0.00;

}

template <>
void NavierStokesWallCondition<3,3>::CalculateNormal(array_1d<double,3>& An )
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    array_1d<double,3> v1,v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);
    An *= 0.5;
}

// template<unsigned int TDim, unsigned int TNumNodes>
// void NavierStokesWallCondition<TDim,TNumNodes>::ApplyNeumannCondition(MatrixType &rLocalMatrix, VectorType &rLocalVector)
// {
//     const NavierStokesWallCondition<TDim,TNumNodes>& rConstThis = *this;
//     if (rConstThis.GetValue(IS_STRUCTURE) == 0.0)
//     {
//         const unsigned int LocalSize = TDim+1;
//         const GeometryType& rGeom = this->GetGeometry();
//         const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
//         const unsigned int NumGauss = IntegrationPoints.size();
//
//         MatrixType NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
//
//         array_1d<double,3> Normal;
//         this->CalculateNormal(Normal); //this already contains the area
//         double A = std::sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);
//         Normal /= A;
//
//         // CAUTION: "Jacobian" is 2.0*A for triangles but 0.5*A for lines
//         double J = (TDim == 2) ? 0.5*A : 2.0*A;
//
//         for (unsigned int g = 0; g < NumGauss; g++)
//         {
//             Vector N = row(NContainer,g);
//             double Weight = J * IntegrationPoints[g].Weight();
//
//             // Neumann boundary condition
//             for (unsigned int i = 0; i < TNumNodes; i++)
//             {
//                 //unsigned int row = i*LocalSize;
//                 const NodeType& rConstNode = this->GetGeometry()[i];
//                 const double pext = rConstNode.FastGetSolutionStepValue(EXTERNAL_PRESSURE);
//
//                 for (unsigned int j = 0; j < TNumNodes; j++)
//                 {
//                     unsigned int row = j*LocalSize;
//                     for (unsigned int d = 0; d < TDim;d++)
//                         rLocalVector[row+d] -= Weight*N[j]*N[i]*pext*Normal[d];
//                 }
//             }
//         }
//     }
// }

template class NavierStokesWallCondition<2,2>;
template class NavierStokesWallCondition<3,3>;

} // namespace Kratos
