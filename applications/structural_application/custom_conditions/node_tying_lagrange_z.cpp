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
#include "custom_conditions/node_tying_lagrange_z.h"
#include "structural_application.h"
#include "utilities/math_utils.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	NodeTyingLagrangeZ::NodeTyingLagrangeZ(IndexType NewId, GeometryType::Pointer 
pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	NodeTyingLagrangeZ::NodeTyingLagrangeZ(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}

	Condition::Pointer NodeTyingLagrangeZ::Create(IndexType NewId, NodesArrayType 
const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new NodeTyingLagrangeZ(NewId, 
GetGeometry().Create(ThisNodes), pProperties));
	}

	NodeTyingLagrangeZ::~NodeTyingLagrangeZ()
	{
	}
    

	//************************************************************************************
	//************************************************************************************
	void NodeTyingLagrangeZ::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		if(rRightHandSideVector.size() != 3)
			rRightHandSideVector.resize(3,false);
        noalias(rRightHandSideVector) = ZeroVector(3);

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void NodeTyingLagrangeZ::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
        
		if(rLeftHandSideMatrix.size1() != 3)
			rLeftHandSideMatrix.resize(3,3,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(3,3);

        if(rRightHandSideVector.size() != 3)
			rRightHandSideVector.resize(3,false);
        noalias(rRightHandSideVector) = ZeroVector(3);
        
        //adding lagrange multipliers to LHS
        rLeftHandSideMatrix(0,2) = 1.0;
        rLeftHandSideMatrix(2,0) = 1.0;
        rLeftHandSideMatrix(1,2) = -1.0;
        rLeftHandSideMatrix(2,1) = -1.0;
        
//         KRATOS_WATCH(rLeftHandSideMatrix);
		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void NodeTyingLagrangeZ::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
        if(rResult.size() != 3)
            rResult.resize(3,false);
        rResult[0] = GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[1] = GetGeometry()[1].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[2] = GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_Z).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	  void NodeTyingLagrangeZ::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
	{
        std::cout << "Initializing Lagrangian Node Tying Condition" << std::endl;
        ConditionalDofList.resize(3);
        ConditionalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_Z));
        ConditionalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_Z));
        ConditionalDofList.push_back(GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_Z));
	}
} // Namespace Kratos

 

