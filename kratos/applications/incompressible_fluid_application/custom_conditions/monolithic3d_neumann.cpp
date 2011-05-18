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
//   Last modified by:    $Author: antonia $
//   Date:                $Date: 2011-03-31 09:41:09 $
//   Revision:            $Revision: 1.2 $
//
// 


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_conditions/monolithic3d_neumann.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	Monolithic3DNeumann::Monolithic3DNeumann(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Monolithic3DNeumann::Monolithic3DNeumann(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{

	}

	Condition::Pointer Monolithic3DNeumann::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new Monolithic3DNeumann(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	Monolithic3DNeumann::~Monolithic3DNeumann()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void Monolithic3DNeumann::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		MatrixType temp = Matrix();
		CalculateLocalSystem(temp, rRightHandSideVector,  rCurrentProcessInfo);

	}

	//************************************************************************************
	//************************************************************************************
	void Monolithic3DNeumann::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
	  //		bool CalculateStiffnessMatrixFlag = true;
	  //bool CalculateResidualVectorFlag = true;
		//KRATOS_WATCH("@@@@@@@@ INSIDE CONDITION@@@@@@@@@@@@@@@@@@@");
// 			if(rRightHandSideVector.size() != 9)
// 					rRightHandSideVector.resize(9,false);
			if(rLeftHandSideMatrix.size1() != 9)
			{
					rLeftHandSideMatrix.resize(9,9,false);
					rRightHandSideVector.resize(9,false);
					
			}
			noalias(rLeftHandSideMatrix) = ZeroMatrix(9,9);
			
			//calculate normal to element
			array_1d<double,3> An,v1,v2;
			v1[0] = GetGeometry()[1].X() - GetGeometry()[0].X();
			v1[1] = GetGeometry()[1].Y() - GetGeometry()[0].Y();
			v1[2] = GetGeometry()[1].Z() - GetGeometry()[0].Z();
						
			v2[0] = GetGeometry()[2].X() - GetGeometry()[0].X();
			v2[1] = GetGeometry()[2].Y() - GetGeometry()[0].Y();
			v2[2] = GetGeometry()[2].Z() - GetGeometry()[0].Z();

			MathUtils<double>::CrossProduct(An,v1,v2);
			An *= 0.5 ;
			
			unsigned int is_structure = this->GetValue(IS_STRUCTURE);

			if(is_structure != 1.0)
			{
			    for(unsigned int i = 0; i<GetGeometry().size(); i++)
				    rRightHandSideVector[i] = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_PRESSURE) * An[i] * 0.3333333333333333; 
			}
			else
			{
			    double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
			    double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
			    double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
			    rRightHandSideVector[0] = -(2.0*p0+p1+p2) * An[0] / 12.0; 
			    rRightHandSideVector[1] = -(p0+2.0*p1+p2) * An[1] / 12.0; 
			    rRightHandSideVector[2] = -(p0+p1+2.0*p2) * An[2] / 12.0;
			    
			    
// 			if(rLeftHandSideMatrix.size1() != 4)
// 			{
// 					rLeftHandSideMatrix.resize(4,4,false);
// 					rRightHandSideVector.resize(4,false);
// 					
// 			}
// 
// 			noalias(rLeftHandSideMatrix) = ZeroMatrix(4,4);
// 
// 			//calculate normal to element.(normal follows the cross rule)
// 			array_1d<double,2> An,edge;
// 			edge[0] = GetGeometry()[1].X() - GetGeometry()[0].X();
// 			edge[1] = GetGeometry()[1].Y() - GetGeometry()[0].Y();
// 						
// 			
// 			double norm = edge[0]*edge[0] + edge[1]*edge[1];
// 			norm = pow(norm,0.5);
// 
// 			An[0] = -1.0*edge[1];
// 			An[1] = edge[0];
// 			//An /= norm; this is then simplified by length of element in integration so is not divided.
// 
// 			double mean_ex_p = 0.0;
// 
// 			for(unsigned int i = 0; i<2 ; i++)
// 				mean_ex_p += 0.5*GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_PRESSURE); 
// 
// 			double p0 = GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_PRESSURE); 
// 			rRightHandSideVector[0] = -0.5*An[0]*p0;
// 			rRightHandSideVector[1] = -0.5*An[1]*p0;
// 
// 			double p1 = GetGeometry()[1].FastGetSolutionStepValue(EXTERNAL_PRESSURE); 
// 			rRightHandSideVector[2] = -0.5*An[0]*p1;
// 			rRightHandSideVector[3] = -0.5*An[1]*p1;

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
//  			KRATOS_WATCH(rRightHandSideVector);
    
// //------------------------------For debugging------------------------------------------    
// 			double& AI0 = GetGeometry()[0].FastGetSolutionStepValue(AUX_INDEX);
// 			double& AI1 = GetGeometry()[1].FastGetSolutionStepValue(AUX_INDEX);
// 			double& AI2 = GetGeometry()[2].FastGetSolutionStepValue(AUX_INDEX);
// 			AI0 = 5.0;
// 			AI1 = 5.0;
// 			AI2 = 5.0;


		
	}

        //************************************************************************************
	//************************************************************************************

	void Monolithic3DNeumann::CalculateLocalVelocityContribution(MatrixType& rDampMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

	int nodes_number = 3;
	int dim = 3;
	unsigned int matsize = nodes_number*(dim);

	if(rDampMatrix.size1() != matsize)
			rDampMatrix.resize(matsize,matsize,false); //false says not to preserve existing storage!!


	noalias(rDampMatrix) = ZeroMatrix(matsize,matsize); 

		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Monolithic3DNeumann::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										ProcessInfo& rCurrentProcessInfo,
										bool CalculateStiffnessMatrixFlag,
										bool CalculateResidualVectorFlag)
	{
		KRATOS_TRY

		KRATOS_CATCH("")
	}



	//************************************************************************************
	//************************************************************************************
	void Monolithic3DNeumann::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
	  void Monolithic3DNeumann::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
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

	  void Monolithic3DNeumann::GetFirstDerivativesVector(Vector& values, int Step)
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
			values[index + 2] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Z,Step);

			//values[index + 2] = GetGeometry()[i].GetSolutionStepValue(PRESSURE,Step);

		}
	}
	//************************************************************************************
	//************************************************************************************
	  void Monolithic3DNeumann::GetSecondDerivativesVector(Vector& values, int Step)
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
			values[index + 2] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Z,Step);

			//values[index + 2] = 0.0;
		}
	}
	//************************************************************************************
	//************************************************************************************

	  
} // Namespace Kratos


