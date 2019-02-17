/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2008-07-24 16:46:55 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/pointforce2D.h"
#include "structural_application.h"
#include "utilities/math_utils.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
PointForce2D::PointForce2D(IndexType NewId, GeometryType::Pointer
                           pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
PointForce2D::PointForce2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer PointForce2D::Create(IndexType NewId, NodesArrayType
                                        const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PointForce2D(NewId,
                              GetGeometry().Create(ThisNodes), pProperties));
}

PointForce2D::~PointForce2D()
{
}


//************************************************************************************
//************************************************************************************
void PointForce2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    if(rRightHandSideVector.size() != 2)
        rRightHandSideVector.resize(2,false);

    array_1d<double,3>& force = GetGeometry()[0].GetSolutionStepValue(FORCE);
    rRightHandSideVector[0] = force[0];
    rRightHandSideVector[1] = force[1];

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void PointForce2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rLeftHandSideMatrix.size1() != 2)
        rLeftHandSideMatrix.resize(2,2,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(2,2);

    if(rRightHandSideVector.size() != 2)
        rRightHandSideVector.resize(2,false);

    array_1d<double,3>& force = GetGeometry()[0].GetSolutionStepValue(FORCE);
    rRightHandSideVector[0] = force[0];
    rRightHandSideVector[1] = force[1];

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************
void PointForce2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int index;
    unsigned int dim = 2;
    rResult.resize(number_of_nodes*dim);
    for (int i=0; i<number_of_nodes; i++)
    {
        index = i*dim;
        rResult[index] = (GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId());
        rResult[index+1] = (GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId());

    }
}

//************************************************************************************
//************************************************************************************
void PointForce2D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
    unsigned int dim = 2;
    ConditionalDofList.resize(GetGeometry().size()*dim);
    unsigned int index;
    for (unsigned int i=0; i<GetGeometry().size(); i++)
    {

        index = i*dim;
        ConditionalDofList[index] = (GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        ConditionalDofList[index+1] = (GetGeometry()[i].pGetDof(DISPLACEMENT_Y));

    }
}
} // Namespace Kratos



