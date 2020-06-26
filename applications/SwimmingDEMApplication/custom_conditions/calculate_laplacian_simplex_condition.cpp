//  Main authors:    Guillermo Casas gcasas@gmail.com

#include "calculate_laplacian_simplex_condition.h"

namespace Kratos
{

///@name
///@{

/**
 * @see ComputeLaplacianSimplexCondition::EquationIdVector
 */
template< unsigned int TDim, unsigned int TNumNodes >
void ComputeLaplacianSimplexCondition<TDim,TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
{
    const SizeType LocalSize = TDim * TNumNodes;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_Y).EquationId();
        if (TDim == 3){
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_Z).EquationId();
        }
    }
}

/**
 * @see ComputeLaplacianSimplexCondition::GetDofList
 */
template< unsigned int TDim, unsigned int TNumNodes >
void ComputeLaplacianSimplexCondition<TDim,TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
{
    const SizeType LocalSize = TDim * TNumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_Y);
        if (TDim == 3){
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_Z);
        }
    }
}

// protected funcions

template <>
void ComputeLaplacianSimplexCondition<2,2>::CalculateNormal(array_1d<double,3>& An)
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = - (pGeometry[1].X() - pGeometry[0].X());
    An[2] =    0.00;

}

template <>
void ComputeLaplacianSimplexCondition<3,3>::CalculateNormal(array_1d<double,3>& An )
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

// Explicit instantiations
template class ComputeLaplacianSimplexCondition<2, 2>;
template class ComputeLaplacianSimplexCondition<3, 3>;
} // namespace Kratos
