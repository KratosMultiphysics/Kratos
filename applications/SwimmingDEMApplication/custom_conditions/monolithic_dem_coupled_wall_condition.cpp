//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Guillermo Casas gcasas@gmail.com
//
#include "monolithic_dem_coupled_wall_condition.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see MonolithicDEMCoupledWallCondition::EquationIdVector
 */
template <>
void MonolithicDEMCoupledWallCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) {
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

    else {
        const unsigned int NumNodes = 2;
        const unsigned int LocalSize = 4;
        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_X).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_Y).EquationId();
        }
    }
}

/**
 * @see MonolithicDEMCoupledWallCondition::EquationIdVector
 */
template <>
void MonolithicDEMCoupledWallCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) {
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

    else {
        const SizeType NumNodes = 3;
        const SizeType LocalSize = 9;
        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_X).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_Y).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_Z).EquationId();
        }
    }
}

/**
 * @see MonolithicDEMCoupledWallCondition::GetDofList
 */
template <>
void MonolithicDEMCoupledWallCondition<2,2>::GetDofList(DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) {
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

    else {
        const SizeType NumNodes = 2;
        const SizeType LocalSize = 4;

        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_X);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_Y);
        }
    }
}

/**
 * @see MonolithicDEMCoupledWallCondition::GetDofList
 */
template <>
void MonolithicDEMCoupledWallCondition<3,3>::GetDofList(DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) {
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

    else {
        const SizeType NumNodes = 3;
        const SizeType LocalSize = 9;

        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_X);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_Y);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_Z);
        }
    }
}

// protected funcions

template class MonolithicDEMCoupledWallCondition<2,2>;
template class MonolithicDEMCoupledWallCondition<3,3>;

} // namespace Kratos
