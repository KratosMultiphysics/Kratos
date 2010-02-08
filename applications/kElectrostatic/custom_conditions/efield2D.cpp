 /*
==============================================================================
KratosR1ElectrostaticApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi,Javier Mora
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
//   Date:                $Date: 2009-09-06 10:30:32 $
//   Revision:            $Revision: 1.4 $
//
// 


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_conditions/efield2D.h"
#include "utilities/math_utils.h"
#include "kElectrostatic.h"


namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	Efield2D::Efield2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Efield2D::Efield2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}

	Condition::Pointer Efield2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new Efield2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	Efield2D::~Efield2D()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void Efield2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = false;
		bool CalculateResidualVectorFlag = true;
		MatrixType temp = Matrix();

		CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
	}

	//************************************************************************************
	//************************************************************************************
	void Efield2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = true;
		
		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
	}

	//************************************************************************************
	//************************************************************************************
	void Efield2D::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										ProcessInfo& rCurrentProcessInfo,
										bool CalculateStiffnessMatrixFlag,
										bool CalculateResidualVectorFlag)
	{
		KRATOS_TRY

		unsigned int number_of_nodes = GetGeometry().size();

		//resizing as needed the LHS
		unsigned int MatSize=number_of_nodes;

		//calculate lenght
		double x21 = GetGeometry()[1].X() - GetGeometry()[0].X();
		double y21 = GetGeometry()[1].Y() - GetGeometry()[0].Y();
		double length = x21*x21 + y21*y21;
		length = sqrt(length);

		double xm = (GetGeometry()[1].X() + GetGeometry()[0].X())/2;
		double ym = (GetGeometry()[1].Y() + GetGeometry()[0].Y())/2;
		double distance = sqrt(xm*xm+ym*ym);

		//double& infinit_coefficient = (this)->GetValue(INFINIT_COEFFICIENT);
		double& infinit_coefficient = (this)->GetValue(INFINIT_COEFFICIENT);
		double convection_coefficient = pow(distance,infinit_coefficient);
		if(infinit_coefficient==0.0)
			convection_coefficient=0.0;
		else
			convection_coefficient = 1/convection_coefficient;

		std::cout << "infinit coefficient=" << infinit_coefficient << std::endl;
		std::cout << "convection_coefficient=" << convection_coefficient << std::endl;

		const double& V0 = GetGeometry()[0].FastGetSolutionStepValue(ELECTROSTATIC_POTENTIAL);
		const double& V1 = GetGeometry()[1].FastGetSolutionStepValue(ELECTROSTATIC_POTENTIAL);

		std::cout << "V0=" << V0 << std::endl;
		std::cout << "V1=" << V1 << std::endl;

		const array_1d<double,3>& ConditionalField = (this)->GetValue(ELECTRIC_DISPLACEMENT_FIELD);
//		const array_1d<double,3>& ConditionalField = GetGeometry()[0].FastGetSolutionStepValue(ELECTRIC_DISPLACEMENT_FIELD);

		const double& c0 = ConditionalField[0];//*x21/length;
		const double& c1 = ConditionalField[1];//*y21/length;
		const double& c2 = ConditionalField[2]*0.0;

		std::cout << "c0=" << c0 << std::endl;
		std::cout << "c1=" << c1 << std::endl;
		std::cout << "c2=" << c2 << std::endl;

		if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
		{
			if(rLeftHandSideMatrix.size1() != MatSize )
				rLeftHandSideMatrix.resize(MatSize,MatSize,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize);	
			
			rLeftHandSideMatrix(0,0) = convection_coefficient * 0.5 * length;
			rLeftHandSideMatrix(1,1) = convection_coefficient * 0.5 * length; 
			//rLeftHandSideMatrix(0,0) = convection_coefficient * length/3.0;
			//rLeftHandSideMatrix(1,1) = convection_coefficient * length/3.0; 
			//rLeftHandSideMatrix(1,0) = convection_coefficient * length/6.0;
			//rLeftHandSideMatrix(0,1) = convection_coefficient * length/6.0; 

		}


		//resizing as needed the RHS
		if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
		{
			if(rRightHandSideVector.size() != MatSize )
				rRightHandSideVector.resize(MatSize,false);
			
			//rRightHandSideVector[0] =  c0;//  -  convection_coefficient *  V0 ;
			//rRightHandSideVector[1] =  c1;//  -  convection_coefficient *  V1 ;


			rRightHandSideVector[0] =  c0 -  convection_coefficient *  V0 ;
			rRightHandSideVector[1] =  c1 -  convection_coefficient *  V1 ;


			rRightHandSideVector *= 0.5*length; //integramos
			
		}

		KRATOS_CATCH("")
	}



	//************************************************************************************
	//************************************************************************************
	void Efield2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			rResult[i] = (GetGeometry()[i].GetDof(ELECTROSTATIC_POTENTIAL)).EquationId();
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void Efield2D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
	{
		ConditionalDofList.resize(GetGeometry().size());
		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			ConditionalDofList[i] = (GetGeometry()[i].pGetDof(ELECTROSTATIC_POTENTIAL));
		}
	}
	  
} // Namespace Kratos


