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
#include "custom_conditions/monolithic2d_neumann.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	Monolithic2DNeumann::Monolithic2DNeumann(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Monolithic2DNeumann::Monolithic2DNeumann(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{

	}

	Condition::Pointer Monolithic2DNeumann::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new Monolithic2DNeumann(NewId, GetGeometry().Create(ThisNodes), pProperties));
KRATOS_WATCH("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
	}

	Monolithic2DNeumann::~Monolithic2DNeumann()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void Monolithic2DNeumann::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		MatrixType temp = Matrix();
		CalculateLocalSystem(temp, rRightHandSideVector,  rCurrentProcessInfo);

	}

	//************************************************************************************
	//************************************************************************************
	void Monolithic2DNeumann::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
	  //		bool CalculateStiffnessMatrixFlag = true;
	  //bool CalculateResidualVectorFlag = true;
		//KRATOS_WATCH("@@@@@@@@ INSIDE CONDITION@@@@@@@@@@@@@@@@@@@");

			if(rLeftHandSideMatrix.size1() != 4)
			{
					rLeftHandSideMatrix.resize(4,4,false);
					rRightHandSideVector.resize(4,false);
					
			}

			noalias(rLeftHandSideMatrix) = ZeroMatrix(4,4);

			//calculate normal to element.(normal follows the cross rule)
			array_1d<double,2> An,edge;
			edge[0] = GetGeometry()[1].X() - GetGeometry()[0].X();
			edge[1] = GetGeometry()[1].Y() - GetGeometry()[0].Y();
						
			
			double norm = edge[0]*edge[0] + edge[1]*edge[1];
			norm = pow(norm,0.5);

			An[0] = -1.0*edge[1];
			An[1] = edge[0];
			//An /= norm; this is then simplified by length of element in integration so is not divided.

			double mean_ex_p = 0.0;

			for(unsigned int i = 0; i<2 ; i++)
				mean_ex_p += 0.5*GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_PRESSURE); 

			double p0 = GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_PRESSURE); 
			rRightHandSideVector[0] = -0.5*An[0]*p0;
			rRightHandSideVector[1] = -0.5*An[1]*p0;

			double p1 = GetGeometry()[1].FastGetSolutionStepValue(EXTERNAL_PRESSURE); 
			rRightHandSideVector[2] = -0.5*An[0]*p1;
			rRightHandSideVector[3] = -0.5*An[1]*p1;

			//	if(mean_ex_p !=GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_PRESSURE))
			//		mean_ex_p = 0.0;
			//KRATOS_WATCH(mean_ex_p);

/*			for(unsigned int ii = 0; ii< 2; ++ii)
				{

					int id = (2 + 1)*(ii);
					rRightHandSideVector[id] = mean_ex_p * An[0]* 0.5;
					rRightHandSideVector[id + 1] = mean_ex_p * An[1]* 0.5;
					rRightHandSideVector[id + 2] = 0.0;
				//KRATOS_WATCH(An);
				}*/
// 			KRATOS_WATCH(An);
//KRATOS_WATCH(p0);
//KRATOS_WATCH(p1);
//KRATOS_WATCH("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
 		//	KRATOS_WATCH(rRightHandSideVector);
		


		
	}

        //************************************************************************************
	//************************************************************************************

	void Monolithic2DNeumann::CalculateLocalVelocityContribution(MatrixType& rDampMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

	int nodes_number = 2;
	int dim = 2;
	unsigned int matsize = nodes_number*(dim);

	if(rDampMatrix.size1() != matsize)
			rDampMatrix.resize(matsize,matsize,false); //false says not to preserve existing storage!!


	noalias(rDampMatrix) = ZeroMatrix(matsize,matsize); 

		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Monolithic2DNeumann::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										ProcessInfo& rCurrentProcessInfo,
										bool CalculateStiffnessMatrixFlag,
										bool CalculateResidualVectorFlag)
	{
		KRATOS_TRY

/*		unsigned int number_of_nodes = GetGeometry().size();
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
*/
		KRATOS_CATCH("")
	}



	//************************************************************************************
	//************************************************************************************
	void Monolithic2DNeumann::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;
		unsigned int node_size = dim;

		
			if(rResult.size() != number_of_nodes*node_size)
				rResult.resize(number_of_nodes*node_size,false);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				rResult[i*node_size] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
				rResult[i*node_size+1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();

			}
		KRATOS_CATCH("")
		
	}

	//************************************************************************************
	//************************************************************************************
	  void Monolithic2DNeumann::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;
		unsigned int node_size = dim;

			if(ElementalDofList.size() != number_of_nodes*node_size)
				ElementalDofList.resize(number_of_nodes*node_size);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				ElementalDofList[i*node_size] = GetGeometry()[i].pGetDof(VELOCITY_X);
				ElementalDofList[i*node_size+1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
			}
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//*************************************************************************************

	  void Monolithic2DNeumann::GetFirstDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes * (dim );
		if(values.size() != MatSize)   values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i * (dim );
			values[index] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y,Step);
			//values[index + 2] = GetGeometry()[i].GetSolutionStepValue(PRESSURE,Step);

		}
	}
	//************************************************************************************
	//************************************************************************************
	  void Monolithic2DNeumann::GetSecondDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes * (dim);
		if(values.size() != MatSize) values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i * (dim );
			values[index] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y,Step);
			//values[index + 2] = 0.0;
		}
	}
	//************************************************************************************
	//************************************************************************************

	  
} // Namespace Kratos


