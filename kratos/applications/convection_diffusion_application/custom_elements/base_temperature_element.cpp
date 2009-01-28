/*
==============================================================================
KratosConvectionDiffusionApplication 
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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes 
 

// External includes 


// Project includes 
#include "includes/define.h"
#include "new_elements/base_temperature_element.h"

#include "utilities/math_utils.h"

namespace Kratos
{
	//initialization of static members
	Element::GeometryType::JacobiansType BaseTemperatureElement::msJ0(1);
	Matrix BaseTemperatureElement::msDN_DX(0,0);
	Matrix BaseTemperatureElement::msInvJ0(0,0);
	Vector BaseTemperatureElement::mstemp(0);

	//************************************************************************************
	//************************************************************************************
	BaseTemperatureElement::BaseTemperatureElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!

	}

	//************************************************************************************
	//************************************************************************************
	BaseTemperatureElement::BaseTemperatureElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
		//setting up the nodal degrees of freedom
		for(int i = 0 ; i != GetGeometry().size() ; ++i)
		{
			(GetGeometry()[i].pAddDof(TEMPERATURE));
		}

		//resizing static members
		unsigned int dim = GetGeometry().WorkingSpaceDimension();
		//msJ0.resize(GetGeometry().IntegrationPointsNumber());
		msDN_DX.resize(GetGeometry().size(),dim);
		msInvJ0.resize(dim,dim);
		mstemp.resize(GetGeometry().size());
	}

	Element::Pointer BaseTemperatureElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new BaseTemperatureElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	BaseTemperatureElement::~BaseTemperatureElement()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void BaseTemperatureElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			unsigned int number_of_points = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();

		//resizing as needed the LHS
		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points,number_of_points); //resetting LHS


		//resizing as needed the RHS
		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points);
		rRightHandSideVector = ZeroVector(number_of_points); //resetting RHS

		//reading integration points and local gradients
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients();

		GetGeometry().Jacobian(msJ0);  

		double DetJ0;
		for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
		{
			//calculating inverse jacobian and jacobian determinant			
			MathUtils<double>::InvertMatrix(msJ0[PointNumber],msInvJ0,DetJ0);

			//Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
			//opb_prod(DN_De[PointNumber],InvJ0,DN_DX,true); 
			noalias(msDN_DX) = prod(DN_De[PointNumber],msInvJ0);
			//noalias(msDN_DX) = prod(DN_De[PointNumber],InvJ0);

			double IntToReferenceWeight = integration_points[PointNumber].Weight() * DetJ0;
			noalias(rLeftHandSideMatrix) += IntToReferenceWeight * prod(msDN_DX, trans(msDN_DX)); //
			//noalias(rLeftHandSideMatrix) += IntToReferenceWeight * prod(msDN_DX, trans(msDN_DX)); //
			
					}

		//calculating external forces
		noalias(rRightHandSideVector) = ZeroVector(number_of_points); //case of zero ext forces

		//if residua is needed, subtracting the internal forces
		// RHS = ExtForces - K*temp;
//		if(rCurrentProcessInfo[RESIDUA_IS_NEEDED] == true)
		{
			for (unsigned int i=0;i<number_of_points;i++)
				mstemp[i] = - GetGeometry()[i].GetSolutionStepValue(TEMPERATURE) ; //this includes the - sign
			
			//axpy_prod(rLeftHandSideMatrix, temp, rRightHandSideVector, false);  //RHS -= K*temp 
			noalias(rRightHandSideVector) = prod(rLeftHandSideMatrix,mstemp);
		}

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void BaseTemperatureElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		MatrixType temp(0,0);
		CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
	}






	//************************************************************************************
	//************************************************************************************
	void BaseTemperatureElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes);
		for (int i=0;i<number_of_nodes;i++)
		{
			rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void BaseTemperatureElement::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);
		for (int i=0;i<number_of_nodes;i++)
		{
			ElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);
		}		
		
		//ElementalDofList.resize(0);
		//for (unsigned int i=0;i<GetGeometry().size();i++)
		//{
		//	ElementalDofList.push_back(GetGeometry()[i].pGetDof(TEMPERATURE));
		//}

	}
} // Namespace Kratos


