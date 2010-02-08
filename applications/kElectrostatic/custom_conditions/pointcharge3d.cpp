/*
==============================================================================
KratosR1ElectrostaticApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2010
Pooyan Dadvand, Riccardo Rossi, Javier Mora 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
mora@cimne.upc.edu
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
//   Last modified by:    $Author: jmora $
//   Date:                $Date: 2010-02-02 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"

#include "custom_conditions/pointcharge3D.h"
#include "kElectrostatic.h"

#include "utilities/math_utils.h"



namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	PointCharge3D::PointCharge3D(IndexType NewId, GeometryType::Pointer 
pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	PointCharge3D::PointCharge3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}

	Condition::Pointer PointCharge3D::Create(IndexType NewId, NodesArrayType 
const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new PointCharge3D(NewId, 
GetGeometry().Create(ThisNodes), pProperties));
	}

	PointCharge3D::~PointCharge3D()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void PointCharge3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		if(rRightHandSideVector.size() != 1)
			rRightHandSideVector.resize(1,false);

		double& point_source = (this)->GetValue(ELECTROSTATIC_POINT_CHARGE);
		//double point_source = GetGeometry()[0].FastGetSolutionStepValue(ELECTROSTATIC_POINT_CHARGE);

		rRightHandSideVector[0] = point_source;
		
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PointCharge3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		if(rLeftHandSideMatrix.size1() != 1)
			rLeftHandSideMatrix.resize(1,1,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(1,1);

		if(rRightHandSideVector.size() != 1)
			rRightHandSideVector.resize(1,false);

		double& point_source = (this)->GetValue(ELECTROSTATIC_POINT_CHARGE);
		//double point_source = GetGeometry()[0].FastGetSolutionStepValue(ELECTROSTATIC_POINT_CHARGE);
		rRightHandSideVector[0] = point_source;
		
		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void PointCharge3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int index;
		unsigned int dim = 1;
		rResult.resize(number_of_nodes*dim);
		for (int i=0;i<number_of_nodes;i++)
		{
			index = i*dim;
			rResult[index] = (GetGeometry()[i].GetDof(ELECTROSTATIC_POTENTIAL).EquationId());
			
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void PointCharge3D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int dim = 1;
		ConditionalDofList.resize(GetGeometry().size()*dim);
		unsigned int index;
		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			
			index = i*dim;
			ConditionalDofList[index] = (GetGeometry()[i].pGetDof(ELECTROSTATIC_POTENTIAL));
			
		}
	}
} // Namespace Kratos

 

