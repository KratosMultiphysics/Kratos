#include "wall_condition_discontinuous.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see WallConditionDiscontinuous::EquationIdVector
 */
template <>
void WallConditionDiscontinuous<2,2>::EquationIdVector(EquationIdVectorType& rResult,
                                                    ProcessInfo& rCurrentProcessInfo)
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    unsigned int step = r_process_info[FRACTIONAL_STEP];
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
    else if(step == 5)
    {
        const SizeType NumNodes = 2;
        const SizeType LocalSize = 2;
        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
        }
    }
    else
    {
        rResult.resize(0,false);
    }
}

/**
 * @see WallConditionDiscontinuous::EquationIdVector
 */
template <>
void WallConditionDiscontinuous<3,3>::EquationIdVector(EquationIdVectorType& rResult,
                                                    ProcessInfo& rCurrentProcessInfo)
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    unsigned int step = r_process_info[FRACTIONAL_STEP];
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
    else if(step == 5)
    {
        const SizeType NumNodes = 3;
        const SizeType LocalSize = 3;
        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
        }
    }
    else
    {
        rResult.resize(0,false);
    }
}

/**
 * @see WallConditionDiscontinuous::GetDofList
 */
template <>
void WallConditionDiscontinuous<2,2>::GetDofList(DofsVectorType& rElementalDofList,
                                              ProcessInfo& rCurrentProcessInfo)
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    unsigned int step = r_process_info[FRACTIONAL_STEP];
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
    else if ( step == 5 )
    {
        const SizeType NumNodes = 2;
        const SizeType LocalSize = 2;

        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
        }
    }

    else
    {
        rElementalDofList.resize(0);
    }
}

/**
 * @see WallConditionDiscontinuous::GetDofList
 */
template <>
void WallConditionDiscontinuous<3,3>::GetDofList(DofsVectorType& rElementalDofList,
                                    ProcessInfo& rCurrentProcessInfo)
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    unsigned int step = r_process_info[FRACTIONAL_STEP];
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
    else if ( step == 5 )
    {
        const SizeType NumNodes = 3;
        const SizeType LocalSize = 3;

        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
        }
    }
    else
    {
        rElementalDofList.resize(0);
    }
}



} // namespace Kratos
