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
//   Revision:            $Revision: 1.4 $
//
// 


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_conditions/thermal_face3D.h"
#include "utilities/math_utils.h"
#include "convection_diffusion_application.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	ThermalFace3D::ThermalFace3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	ThermalFace3D::ThermalFace3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}

	Condition::Pointer ThermalFace3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new ThermalFace3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	ThermalFace3D::~ThermalFace3D()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void ThermalFace3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = false;
		bool CalculateResidualVectorFlag = true;
		MatrixType temp = Matrix();

		CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
	}

	//************************************************************************************
	//************************************************************************************
	void ThermalFace3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = true;
		
		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
	}

	//************************************************************************************
	//************************************************************************************
	void ThermalFace3D::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										ProcessInfo& rCurrentProcessInfo,
										bool CalculateStiffnessMatrixFlag,
										bool CalculateResidualVectorFlag)
	{
		KRATOS_TRY

		unsigned int number_of_nodes = GetGeometry().size();

		//resizing as needed the LHS
		unsigned int MatSize=number_of_nodes;

		//calculate area
		double area = GetGeometry().Area();
		
		ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
		const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
		const Variable<double>& rSurfaceSourceVar = my_settings->GetSurfaceSourceVariable();

		const double& ambient_temperature = GetProperties()[AMBIENT_TEMPERATURE];
		double StefenBoltzmann = 5.67e-8;
		double emissivity = GetProperties()[EMISSIVITY];
		double convection_coefficient = GetProperties()[CONVECTION_COEFFICIENT];

		const double& T0 = GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);
		const double& T1 = GetGeometry()[1].FastGetSolutionStepValue(rUnknownVar);
		const double& T2 = GetGeometry()[2].FastGetSolutionStepValue(rUnknownVar);

		const double& q0 = GetGeometry()[0].FastGetSolutionStepValue(rSurfaceSourceVar);
		const double& q1 = GetGeometry()[1].FastGetSolutionStepValue(rSurfaceSourceVar);
		const double& q2 = GetGeometry()[2].FastGetSolutionStepValue(rSurfaceSourceVar);

		if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
		{
			if(rLeftHandSideMatrix.size1() != MatSize )
				rLeftHandSideMatrix.resize(MatSize,MatSize,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize);
			rLeftHandSideMatrix(0,0) = ( convection_coefficient + emissivity*StefenBoltzmann*4.0*pow(T0,3)  )* 0.333333333333 * area;
			rLeftHandSideMatrix(1,1) = ( convection_coefficient + emissivity*StefenBoltzmann*4.0*pow(T1,3)  )* 0.333333333333 * area; 
			rLeftHandSideMatrix(2,2) = ( convection_coefficient + emissivity*StefenBoltzmann*4.0*pow(T2,3)  )* 0.333333333333 * area; 
		}

		//resizing as needed the RHS
		double aux = pow(ambient_temperature,4);
		if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
		{
			if(rRightHandSideVector.size() != MatSize )
				rRightHandSideVector.resize(MatSize,false);
			
			rRightHandSideVector[0] =  emissivity*q0 - emissivity*StefenBoltzmann*(pow(T0,4) - aux)   
									-  convection_coefficient * ( T0 - ambient_temperature);

			rRightHandSideVector[1] =  emissivity*q1  - emissivity*StefenBoltzmann*(pow(T1,4) - aux)  
									-  convection_coefficient * ( T1 - ambient_temperature);

			rRightHandSideVector[2] =  emissivity*q2  - emissivity*StefenBoltzmann*(pow(T2,4) - aux)  
									-  convection_coefficient * ( T2 - ambient_temperature);
			
			rRightHandSideVector *= 0.3333333333333*area;
			
		}
		KRATOS_CATCH("")
	}



	//************************************************************************************
	//************************************************************************************
	void ThermalFace3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();

		ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        	const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();


		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			rResult[i] = (GetGeometry()[i].GetDof(rUnknownVar)).EquationId();
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void ThermalFace3D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
	{
		ConditionalDofList.resize(GetGeometry().size());

		ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        	const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			ConditionalDofList[i] = (GetGeometry()[i].pGetDof(rUnknownVar));
		}
	}
	  
} // Namespace Kratos


