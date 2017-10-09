/*
==============================================================================
KratosShallowWaterApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//
//   Project Name:        Kratos
//   Last modified by:    Miguel Masó Sotomayor
//   Date:                September 27th 2017
//   Revision:            1.0
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/nothing_condition.hpp"
#include "shallow_water_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void NothingCondition<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        unsigned int local_size = TNumNodes*3;
        if(rResult.size() != local_size)
            rResult.resize(local_size,false);                           // False says not to preserve existing storage!!
        
        GeometryType& rGeom = GetGeometry();
        int counter=0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[counter++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
            rResult[counter++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
            rResult[counter++] = rGeom[i].GetDof(HEIGHT).EquationId();
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void NothingCondition<TNumNodes>::GetDofList(DofsVectorType& rConditionDofList,ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const unsigned int local_size = TNumNodes*3;
        if(rConditionDofList.size() != local_size)
            rConditionDofList.resize(local_size);

        GeometryType& rGeom = GetGeometry();
        int counter=0;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rConditionDofList[counter++] = rGeom[i].pGetDof(VELOCITY_X);
            rConditionDofList[counter++] = rGeom[i].pGetDof(VELOCITY_Y);
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
        unsigned int local_size = TNumNodes*3;
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
        unsigned int local_size = TNumNodes*3;
        if(rRightHandSideVector.size() != local_size)
            rRightHandSideVector.resize(local_size,false);              // False says not to preserve existing storage!!

        noalias(rRightHandSideVector) = ZeroVector(local_size);

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

template class NothingCondition<2>;

} // namespace Kratos
