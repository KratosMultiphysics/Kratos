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
#include "custom_conditions/Point_Neumann3D.h"
#include "utilities/math_utils.h"
#include "ULF_application.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	PointNeumann3D::PointNeumann3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	PointNeumann3D::PointNeumann3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
		//KRATOS_WATCH("CREATING ============= POINTNEUMANN3D  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
	}

	Condition::Pointer PointNeumann3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new PointNeumann3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_WATCH("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
	}

	PointNeumann3D::~PointNeumann3D()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void PointNeumann3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags

		MatrixType temp = Matrix();
		CalculateLocalSystem(temp, rRightHandSideVector,  rCurrentProcessInfo);

	}

	//************************************************************************************
	//************************************************************************************
	void PointNeumann3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{

	//KRATOS_WATCH("POINT NEUMANNNN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11")
		//calculation flags
	  //		bool CalculateStiffnessMatrixFlag = true;
	  //bool CalculateResidualVectorFlag = true;
		//KRATOS_WATCH("@@@@@@@@ INSIDE CONDITION@@@@@@@@@@@@@@@@@@@");

			if(rLeftHandSideMatrix.size1() != 3)
			{
					rLeftHandSideMatrix.resize(3,3,false);
					rRightHandSideVector.resize(3,false);
					
			}

			noalias(rLeftHandSideMatrix) = ZeroMatrix(3,3);

	int nodes_number = 1;
	int dim = 3;
	unsigned int matsize = nodes_number*(dim);

	//if(rDampingMatrix.size1() != matsize)
	//		rDampingMatrix.resize(matsize,matsize,false); //false says not to preserve existing storage!!


	//noalias(rDampingMatrix) = ZeroMatrix(matsize,matsize); 



			//calculate normal to element.(normal follows the cross rule)
			array_1d<double,3> An,edge;
			An = GetGeometry()[0].FastGetSolutionStepValue(NORMAL);
			//KRATOS_WATCH(An)
			double ext_pr = -1.0*GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_PRESSURE);


			double temp=An[0]*ext_pr;
			double temp1=An[1]*ext_pr;
			double temp2=An[2]*ext_pr;

//	KRATOS_WATCH("¿¿¿¿¿¿¿¿¿")
			//KRATOS_WATCH(temp)
			//KRATOS_WATCH(temp1)
			//KRATOS_WATCH(temp2)
			//KRATOS_WATCH("¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿/()(/)(/()/()/)(/)(/()¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿")
			//add to RHS
			rRightHandSideVector[0] = An[0]*ext_pr;
			rRightHandSideVector[1] = An[1]*ext_pr;
			rRightHandSideVector[2] = An[2]*ext_pr;

			//KRATOS_WATCH("external_pressure contrib =============================================================================================")
			//KRATOS_WATCH(rRightHandSideVector)
			
			//AND NOW WE ADD THE VISCOUS STRESS CONTRIBUTION
			/*
			const array_1d<double, 3 > &  s_x=GetGeometry()[0].FastGetSolutionStepValue(STRESSX);											   
			const array_1d<double, 3 > &  s_y=GetGeometry()[0].FastGetSolutionStepValue(STRESSY);
			const array_1d<double, 3 > &  s_z=GetGeometry()[0].FastGetSolutionStepValue(STRESSZ);
			
			rRightHandSideVector[0] += s_x[0]*An[0]+s_x[1]*An[1]+s_x[0]*An[2];
			rRightHandSideVector[1] += s_y[0]*An[0]+s_y[1]*An[1]+s_y[0]*An[2];
			rRightHandSideVector[2] += s_z[0]*An[0]+s_z[1]*An[1]+s_z[0]*An[2];
			*/
			//JUST TO CHECK THE DIRECTION
			//GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY_X)=s_x[0]*An[0]+s_x[1]*An[1]+s_x[0]*An[2];
			//GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY_Y)=s_y[0]*An[0]+s_y[1]*An[1]+s_y[0]*An[2];
			//GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY_Z)=s_z[0]*An[0]+s_z[1]*An[1]+s_z[0]*An[2];





		
	}

        //************************************************************************************
	//************************************************************************************
	/*
	void PointNeumann3D::CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
	KRATOS_WATCH("")
	int nodes_number = 1;
	int dim = 3;
	unsigned int matsize = nodes_number*(dim);

	if(rDampingMatrix.size1() != matsize)
			rDampingMatrix.resize(matsize,matsize,false); //false says not to preserve existing storage!!


	noalias(rDampingMatrix) = ZeroMatrix(matsize,matsize); 



			//calculate normal to element.(normal follows the cross rule)
			array_1d<double,3> An,edge;
			An = GetGeometry()[0].FastGetSolutionStepValue(NORMAL);

			double ext_pr = -1.0*GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_PRESSURE);

			//add to RHS
			rRightHandSideVector[0] = An[0]*ext_pr;
			rRightHandSideVector[1] = An[1]*ext_pr;
			rRightHandSideVector[2] = An[2]*ext_pr;

		KRATOS_CATCH("")
	}
	*/
	//************************************************************************************
	//************************************************************************************
	void PointNeumann3D::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
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
	void PointNeumann3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 3;
		unsigned int node_size = dim;

		
			if(rResult.size() != number_of_nodes*node_size)
				rResult.resize(number_of_nodes*node_size,false);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				rResult[i*node_size] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
				rResult[i*node_size+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
				rResult[i*node_size+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();

			}
		KRATOS_CATCH("")
		
	}

	//************************************************************************************
	//************************************************************************************
	  void PointNeumann3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = 1;
		unsigned int dim = 3;
		unsigned int node_size = dim;

			if(ElementalDofList.size() != number_of_nodes*node_size)
				ElementalDofList.resize(number_of_nodes*node_size);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				ElementalDofList[i*node_size] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
				ElementalDofList[i*node_size+1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
				ElementalDofList[i*node_size+2] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
			}
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//*************************************************************************************

	  void PointNeumann3D::GetFirstDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = 3;
		unsigned int MatSize = number_of_nodes * (dim );

		if(values.size() != MatSize)   values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i * (dim );
			values[index] = 0.0;
			values[index + 1] = 0.0;
			values[index + 2] = 0.0;
			//values[index + 2] = GetGeometry()[i].GetSolutionStepValue(PRESSURE,Step);

		}
	}
	//************************************************************************************
	//************************************************************************************
	  void PointNeumann3D::GetSecondDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes * (dim);
		if(values.size() != MatSize) values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i * (dim );
			values[index] = 0.0;
			values[index + 1] = 0.0;
			values[index + 2] = 0.0;
			//values[index + 2] = 0.0;
		}
	}
	//************************************************************************************
	//************************************************************************************

	  
} // Namespace Kratos


