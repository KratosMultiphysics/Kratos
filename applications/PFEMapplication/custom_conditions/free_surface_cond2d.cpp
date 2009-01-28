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
//   Date:                $Date: 2007-03-28 17:37:14 $
//   Revision:            $Revision: 1.3 $
//
// 


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_conditions/free_surface_cond2d.h"
#include "utilities/math_utils.h"
#include "PFEM_application.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	FreeSurfaceCond2d::FreeSurfaceCond2d(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	FreeSurfaceCond2d::FreeSurfaceCond2d(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}

	Condition::Pointer FreeSurfaceCond2d::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new FreeSurfaceCond2d(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	FreeSurfaceCond2d::~FreeSurfaceCond2d()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void FreeSurfaceCond2d::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		unsigned int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];
		
		if(FractionalStepNumber == 4) //pressure step
		{
			bool CalculateStiffnessMatrixFlag = false;
			bool CalculateResidualVectorFlag = true;
			MatrixType temp = Matrix();
			CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
		}
		else
		{
			if(rRightHandSideVector.size() != 0)
						rRightHandSideVector.resize(0,false);
		}
	}

	//************************************************************************************
	//************************************************************************************
	void FreeSurfaceCond2d::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = true;

		unsigned int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

//		KRATOS_WATCH(FractionalStepNumber);

		if(FractionalStepNumber == 4) //pressure step
		{
			CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
		}
		else
		{
			if(rLeftHandSideMatrix.size1() != 0)
			{
				rLeftHandSideMatrix.resize(0,0,false);
				rRightHandSideVector.resize(0,false);
			}
		}
	}

	//************************************************************************************
	//************************************************************************************
	void FreeSurfaceCond2d::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										ProcessInfo& rCurrentProcessInfo,
										bool CalculateStiffnessMatrixFlag,
										bool CalculateResidualVectorFlag)
	{
		KRATOS_TRY

		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int MatSize=number_of_nodes;

		//calculate lenght
		double x21 = GetGeometry()[1].X() - GetGeometry()[0].X();
		double y21 = GetGeometry()[1].Y() - GetGeometry()[0].Y();
		double lenght = x21*x21 + y21*y21;
		lenght = sqrt(lenght);

		double dt = rCurrentProcessInfo[DELTA_TIME];

		//proposal by riccardo
		double density = 0.5*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY) + GetGeometry()[1].FastGetSolutionStepValue(DENSITY));

		double sound_speed = 10.0;
//		double K = density * sound_speed; //K = density * sound_speed*sound_speed

		//chapuza to make this of the same order of the other terms in the system independently on the deltatime
//		double K = density * 0.0001 / (dt*dt);  //.this would give a c = 1 for dt = 0.01...which we tested to be good

//		K = 1000;

	


		if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
		{
			if(rLeftHandSideMatrix.size1() != MatSize )
				rLeftHandSideMatrix.resize(MatSize,MatSize,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize);	
//			if(GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1)
//				rLeftHandSideMatrix(0,0)  = 0.5 * lenght*lenght / (  K *dt);
//			if(GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1)
//				rLeftHandSideMatrix(1,1) = 0.5 * lenght*lenght / (  K * dt);
			if(GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1)
				rLeftHandSideMatrix(0,0)  = 0.5 * lenght  / ( density * sound_speed * dt );
			if(GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1)
				rLeftHandSideMatrix(1,1) = 0.5 * lenght / (  density * sound_speed *dt );
		}

		if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
		{
			if(rRightHandSideVector.size() != MatSize )
				rRightHandSideVector.resize(MatSize,false);

//			noalias(rRightHandSideVector) = ZeroVector(2);

			array_1d<double,2> temp;
//			temp[0] = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1) - GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
//			temp[1] = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1) - GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
			temp[0] = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT) - GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
			temp[1] = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT) - GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
			noalias(rRightHandSideVector) =  prod( rLeftHandSideMatrix , temp);
		}
//		KRATOS_WATCH(rLeftHandSideMatrix);
//		KRATOS_WATCH(rRightHandSideVector);

		KRATOS_CATCH("")
	}



	//************************************************************************************
	//************************************************************************************
	void FreeSurfaceCond2d::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();

		unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber == 4) //pressure step
		{
			if(rResult.size() != number_of_nodes)
				rResult.resize(number_of_nodes,false);
			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				rResult[i] = (GetGeometry()[i].GetDof(PRESSURE)).EquationId();
			}
		}
		else
			if(rResult.size() != 0)
				rResult.resize(0,false);

	}

	//************************************************************************************
	//************************************************************************************
	  void FreeSurfaceCond2d::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber == 4) //pressure step
		{
			if(ConditionalDofList.size() != GetGeometry().size())
				ConditionalDofList.resize(GetGeometry().size());
			for (unsigned int i=0;i<GetGeometry().size();i++)
			{
				ConditionalDofList[i] = (GetGeometry()[i].pGetDof(PRESSURE));
			}
		}
		else
			if(ConditionalDofList.size() != 0)
				ConditionalDofList.resize(0);
	}
	  
} // Namespace Kratos


