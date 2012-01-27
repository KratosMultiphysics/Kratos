/*
==============================================================================
KratosPFEMApplication 
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
//   Date:                $Date: 2008-02-14 09:41:09 $
//   Revision:            $Revision: 1.2 $
//
// 


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/environment_contact.h"
#include "utilities/math_utils.h"
#include "thermo_mechanical_application.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	EnvironmentContact::EnvironmentContact(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	EnvironmentContact::EnvironmentContact(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{

	}

	Condition::Pointer EnvironmentContact::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new EnvironmentContact(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	EnvironmentContact::~EnvironmentContact()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void EnvironmentContact::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		MatrixType temp = Matrix();
// 		CalculateLocalSystem(temp, rRightHandSideVector,  rCurrentProcessInfo);

	}

	//************************************************************************************
	//************************************************************************************
	void EnvironmentContact::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
	  
	  int nodes_number = GetGeometry().size();
	  
	  if(rLeftHandSideMatrix.size1() != 1)
	  {
			  rLeftHandSideMatrix.resize(1,1,false);
			  rRightHandSideVector.resize(1,false);
			  
	  }

		    rRightHandSideVector = ZeroVector(1);
		    array_1d<double,2> length_normal;
		    length_normal = GetGeometry()[0].FastGetSolutionStepValue(NORMAL);		
		    double length = norm_2(length_normal);
		    
		    //take thermal properties
		  ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
		  const Variable<double>& rTransferCoefficientVar = my_settings->GetTransferCoefficientVariable();
		  
	  //         double HTC_Alpha = GetProperties()[rDiffusionVar]; 
		  double HTC_Alpha = GetGeometry()[0].FastGetSolutionStepValue(rTransferCoefficientVar);
		  HTC_Alpha *= length;

	  if(length == 0.0)
	    KRATOS_WATCH("NORMAL is ZERO")
		  rLeftHandSideMatrix(0,0) = HTC_Alpha;
			  
		  //Residual
		  const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();	
		  array_1d<double, 1 > unknown_vec;
		  unknown_vec[0] =  GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);	
		  
		  double amb_T = rCurrentProcessInfo[AMBIENT_TEMPERATURE];
		  unknown_vec[0] -= amb_T;

		  noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, unknown_vec);	
	  
	}

        //************************************************************************************
	//************************************************************************************

	void EnvironmentContact::CalculateLocalVelocityContribution(MatrixType& rDampMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

// 	int nodes_number = 4;
// 	int dim = 2;
// 	unsigned int matsize = nodes_number*(dim);
// 
// 	if(rDampMatrix.size1() != matsize)
// 			rDampMatrix.resize(matsize,matsize,false); //false says not to preserve existing storage!!




		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void EnvironmentContact::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										ProcessInfo& rCurrentProcessInfo,
										bool CalculateStiffnessMatrixFlag,
										bool CalculateResidualVectorFlag)
	{
		KRATOS_TRY


		KRATOS_CATCH("")
	}



	//************************************************************************************
	//************************************************************************************
	void EnvironmentContact::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);	
		const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
		
			if(rResult.size() != number_of_nodes)
				rResult.resize(number_of_nodes,false);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
			}
		KRATOS_CATCH("")
		
	}

	//************************************************************************************
	//************************************************************************************
	  void EnvironmentContact::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);	
		const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

			if(ElementalDofList.size() != number_of_nodes)
				ElementalDofList.resize(number_of_nodes);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
                                ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

			}
		KRATOS_CATCH("");
	}



	  
} // Namespace Kratos


