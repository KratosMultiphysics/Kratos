//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: Kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "shallow_water_application_variables.h"
#include "custom_conditions/nothing_condition.hpp"

namespace Kratos
{

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void NothingCondition<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        unsigned int local_size = TNumNodes*1;
        if(rResult.size() != local_size)
            rResult.resize(local_size,false);                           // False says not to preserve existing storage!!
        
        GeometryType& rGeom = GetGeometry();
        int counter=0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[counter++] = rGeom[i].GetDof(HEIGHT).EquationId();
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void NothingCondition<TNumNodes>::GetDofList(DofsVectorType& rConditionDofList,ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const unsigned int local_size = TNumNodes*1;
        if(rConditionDofList.size() != local_size)
            rConditionDofList.resize(local_size);

        GeometryType& rGeom = GetGeometry();
        int counter=0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rConditionDofList[counter++] = rGeom[i].pGetDof(HEIGHT);
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void NothingCondition<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Resize the Left and Right Hand side
        unsigned int local_size = TNumNodes*1;
        if(rLeftHandSideMatrix.size1() != local_size)
            rLeftHandSideMatrix.resize(local_size,local_size,false);    // False says not to preserve existing storage!!

        if(rRightHandSideVector.size() != local_size)
            rRightHandSideVector.resize(local_size,false);              // False says not to preserve existing storage!!

        noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size,local_size);
        
        noalias(rRightHandSideVector) = ZeroVector(local_size);

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void NothingCondition<TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Resize the Right Hand Side
        unsigned int local_size = TNumNodes*1;
        if(rRightHandSideVector.size() != local_size)
            rRightHandSideVector.resize(local_size,false);              // False says not to preserve existing storage!!

        noalias(rRightHandSideVector) = ZeroVector(local_size);

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

template class NothingCondition<2>;

} // namespace Kratos
