#include "wall_condition.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see WallCondition::EquationIdVector
 */
template <>
void WallCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,
                                                    ProcessInfo& rCurrentProcessInfo)
{
    int step = rCurrentProcessInfo[FRACTIONAL_STEP];
    if ( step == 1 )
    {
        const unsigned int NumNodes = 2;
        const unsigned int LocalSize = 4;
        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        }
    }
    else
    {
        rResult.resize(0,false);
    }
}

/**
 * @see WallCondition::EquationIdVector
 */
template <>
void WallCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,
                                                    ProcessInfo& rCurrentProcessInfo)
{
    int step = rCurrentProcessInfo[FRACTIONAL_STEP];
    if ( step == 1 )
    {
        const SizeType NumNodes = 3;
        const SizeType LocalSize = 9;
        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        }
    }
    else
    {
        rResult.resize(0,false);
    }
}

/**
 * @see WallCondition::GetDofList
 */
template <>
void WallCondition<2,2>::GetDofList(DofsVectorType& rElementalDofList,
                                              ProcessInfo& rCurrentProcessInfo)
{
 	int step = rCurrentProcessInfo[FRACTIONAL_STEP];
    if ( step == 1 )
    {
        const SizeType NumNodes = 2;
        const SizeType LocalSize = 4;

        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        }
    }
    else
    {
        rElementalDofList.resize(0);
    }
}

/**
 * @see WallCondition::GetDofList
 */
template <>
void WallCondition<3,3>::GetDofList(DofsVectorType& rElementalDofList,
                                    ProcessInfo& rCurrentProcessInfo)
{
	int step = rCurrentProcessInfo[FRACTIONAL_STEP];
    if ( step == 1 )
    {
        const SizeType NumNodes = 3;
        const SizeType LocalSize = 9;

        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        }
    }
    else
    {
        rElementalDofList.resize(0);
    }
}

template <>
 void WallCondition<2,2>::CalculateNormal(array_1d<double,3>& An)
    {
        Geometry<Node<3> >& pGeometry = this->GetGeometry();

        An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
        An[1] = - (pGeometry[1].X() - pGeometry[0].X());
        An[2] =    0.00;

    }

template <>
void WallCondition<3,3>::CalculateNormal(array_1d<double,3>& An )
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
        An *= -0.5;

        // 				noalias((it)->GetValue(NORMAL)) = An;
    }

} // namespace Kratos
