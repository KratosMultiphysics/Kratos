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
//   Last modified by:    $Author: janosch $
//   Date:                $Date: 2008-10-23 12:26:35 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_conditions/node_tying_lagrange.h"
#include "structural_application.h"
#include "utilities/math_utils.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	NodeTyingLagrange::NodeTyingLagrange(IndexType NewId, GeometryType::Pointer 
pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	NodeTyingLagrange::NodeTyingLagrange(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}

	Condition::Pointer NodeTyingLagrange::Create(IndexType NewId, NodesArrayType 
const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new NodeTyingLagrange(NewId, 
GetGeometry().Create(ThisNodes), pProperties));
	}

	NodeTyingLagrange::~NodeTyingLagrange()
	{
	}
    

	//************************************************************************************
	//************************************************************************************
	void NodeTyingLagrange::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		if(rRightHandSideVector.size() != 9)
			rRightHandSideVector.resize(9,false);
        noalias(rRightHandSideVector) = ZeroVector(9);

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void NodeTyingLagrange::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		if(rLeftHandSideMatrix.size1() != 9)
			rLeftHandSideMatrix.resize(9,9,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(9,9);

		if(rRightHandSideVector.size() != 9)
			rRightHandSideVector.resize(9,false);
        noalias(rRightHandSideVector) = ZeroVector(9);
        
        //adding lagrange multipliers to LHS
        for( int j=0; j<3; j++ )
        {
            rLeftHandSideMatrix(j,j+6) = 1.0;
            rLeftHandSideMatrix(j+6,j) = 1.0;
            rLeftHandSideMatrix(j+3,j+6) = -1.0;
            rLeftHandSideMatrix(j+6,j+3) = -1.0;
        }
        KRATOS_WATCH(rLeftHandSideMatrix);
		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void NodeTyingLagrange::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
        if(rResult.size() != 9)
            rResult.resize(9,false);
        rResult[0] = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
        rResult[1] = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[2] = GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[3] = GetGeometry()[1].GetDof(DISPLACEMENT_X).EquationId();
        rResult[4] = GetGeometry()[1].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[5] = GetGeometry()[1].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[6] = GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_X).EquationId();
        rResult[7] = GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_Y).EquationId();
        rResult[8] = GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_Z).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	  void NodeTyingLagrange::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
	{
        std::cout << "Initializing Lagrangian Node Tying Condition" << std::endl;
        ConditionalDofList.resize(9);
        ConditionalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_X));
        ConditionalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_Y));
        ConditionalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_Z));
        ConditionalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_X));
        ConditionalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_Y));
        ConditionalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_Z));
        ConditionalDofList.push_back(GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_X));
        ConditionalDofList.push_back(GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_Y));
        ConditionalDofList.push_back(GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_Z));
	}
} // Namespace Kratos

 

