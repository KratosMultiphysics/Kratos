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
#include "custom_conditions/Point_Neumann_Monolithic2D.h"
#include "utilities/math_utils.h"
#include "ULF_application.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	PointNeumannMonolithic2D::PointNeumannMonolithic2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	PointNeumannMonolithic2D::PointNeumannMonolithic2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
		//KRATOS_WATCH("CREATING ============= PointNeumannMonolithic2D  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
	}

	Condition::Pointer PointNeumannMonolithic2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new PointNeumannMonolithic2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_WATCH("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
	}

	PointNeumannMonolithic2D::~PointNeumannMonolithic2D()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void PointNeumannMonolithic2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags

		MatrixType temp = Matrix();
		CalculateLocalSystem(temp, rRightHandSideVector,  rCurrentProcessInfo);

	}

	//************************************************************************************
	//************************************************************************************
	void PointNeumannMonolithic2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
			if(rLeftHandSideMatrix.size1() != 2)
			{
					rLeftHandSideMatrix.resize(2,2,false);
					rRightHandSideVector.resize(2,false);
					
			}

			noalias(rLeftHandSideMatrix) = ZeroMatrix(2,2);
			
			// Pressure component of the stress
			//calculate normal to element.(normal follows the cross rule)
			array_1d<double,2> An,edge,vel;
			An = GetGeometry()[0].FastGetSolutionStepValue(NORMAL);
			double ext_pr = -1.0*GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_PRESSURE);

			//add to RHS
			rRightHandSideVector[0] = An[0]*ext_pr;
			rRightHandSideVector[1] = An[1]*ext_pr;

			// Viscous component of the stress
			edge[0] = -An[1];
			edge[1] = An[0];
			vel = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
			double temp = vel[0]*edge[0] + vel[1]*edge[1];
			if(temp < 0.0)
			{
			  edge[0] *= -1.0;
			  edge[1] *= -1.0;
			}
			
// 			double viscst_xx = GetGeometry()[0].FastGetSolutionStepValue(VISCOUS_STRESSX_X);
// 			double viscst_xy = GetGeometry()[0].FastGetSolutionStepValue(VISCOUS_STRESSX_Y);
// 			double viscst_yx = GetGeometry()[0].FastGetSolutionStepValue(VISCOUS_STRESSY_X);
// 			double viscst_yy = GetGeometry()[0].FastGetSolutionStepValue(VISCOUS_STRESSY_Y);
// 			rRightHandSideVector[0] += edge[0]*viscst_xx + edge[1]*viscst_xy;
// 			rRightHandSideVector[1] += edge[0]*viscst_yx + edge[1]*viscst_yy;
			
			
	}
	//************************************************************************************
	//************************************************************************************
	void PointNeumannMonolithic2D::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										ProcessInfo& rCurrentProcessInfo,
										bool CalculateStiffnessMatrixFlag,
										bool CalculateResidualVectorFlag)
	{
		KRATOS_TRY

KRATOS_THROW_ERROR(std::logic_error,"Method not implemented!!!!","");

		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void PointNeumannMonolithic2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;
		unsigned int node_size = dim;
// 		const double dt = CurrentProcessInfo.GetValue(DELTA_TIME);
			if(rResult.size() != number_of_nodes*node_size)
				rResult.resize(number_of_nodes*node_size,false);
			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				/*
				rResult[i*node_size] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
				rResult[i*node_size+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
				//rResult[i*node_size+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
				*/
				rResult[i*node_size] = (GetGeometry()[i].GetDof(VELOCITY_X).EquationId());
				rResult[i*node_size+1] = (GetGeometry()[i].GetDof(VELOCITY_Y).EquationId());
// 				rResult[i*node_size+2] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
			}
		KRATOS_CATCH("")
		
	}

	//************************************************************************************
	//************************************************************************************
	  void PointNeumannMonolithic2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = 1;
		unsigned int dim = 2;
		unsigned int node_size = dim;

			if(ElementalDofList.size() != number_of_nodes*node_size)
				ElementalDofList.resize(number_of_nodes*node_size);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				/*
				ElementalDofList[i*node_size] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
				ElementalDofList[i*node_size+1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
				//ElementalDofList[i*node_size+2] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
				*/
				ElementalDofList[i*node_size] = GetGeometry()[i].pGetDof(VELOCITY_X);
				ElementalDofList[i*node_size+1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
// 				ElementalDofList[i*node_size+2] = GetGeometry()[i].pGetDof(PRESSURE);
			}
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//*************************************************************************************

	  void PointNeumannMonolithic2D::GetFirstDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = 2;
		unsigned int MatSize = number_of_nodes * (dim );

		if(values.size() != MatSize)   values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i * (dim );
			values[index] = 0.0;
			values[index + 1] = 0.0;
// 			values[index + 2] = 0.0;
// 			values[index + 2] = GetGeometry()[i].GetSolutionStepValue(PRESSURE,Step);

		}
	}
	//************************************************************************************
	//************************************************************************************
	  void PointNeumannMonolithic2D::GetSecondDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
// 		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		const unsigned int dim = 2;
		unsigned int MatSize = number_of_nodes * (dim);
		if(values.size() != MatSize) values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i * (dim );
			values[index] = 0.0;
			values[index + 1] = 0.0;
// 			values[index + 2] = 0.0;
			//values[index + 2] = 0.0;
		}
	}
	//************************************************************************************
	//************************************************************************************

	  
} // Namespace Kratos


