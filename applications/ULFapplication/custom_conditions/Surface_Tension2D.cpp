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
#include "custom_conditions/Surface_Tension2D.h"
#include "utilities/math_utils.h"
#include "ULF_application.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	SurfaceTension2D::SurfaceTension2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	SurfaceTension2D::SurfaceTension2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
		//KRATOS_WATCH("CREATING ============= TEMSIOMCOND  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
	}

	Condition::Pointer SurfaceTension2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new SurfaceTension2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		//KRATOS_WATCH("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
	}

	SurfaceTension2D::~SurfaceTension2D()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void SurfaceTension2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags

		MatrixType temp = Matrix(2,2);
		//MatrixType temp = Matrix();
		CalculateLocalSystem(temp, rRightHandSideVector,  rCurrentProcessInfo);

	}

	//************************************************************************************
	//************************************************************************************
	void SurfaceTension2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
	  rRightHandSideVector.resize(2,false);
	  noalias(rRightHandSideVector) = ZeroVector(2);
	  rLeftHandSideMatrix.resize(2,2,false);
	  noalias(rLeftHandSideMatrix) = ZeroMatrix(2,2);
	  
	  // 20140917: Remember to remove the FractionalStepNumber variable and the if condition if the scheme is not the fractional step!!
	  int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];
	  
	  if(FractionalStepNumber == 1)
	  {
	    array_1d<double,2> An;
	    An = GetGeometry()[0].FastGetSolutionStepValue(NORMAL);
	    double curv = GetGeometry()[0].FastGetSolutionStepValue(CURVATURE);
	    double gamma = 0.072; // surface tension coefficient between air and water, N m-1
	  
	    rRightHandSideVector[0] = -curv*gamma*An[0];
	    rRightHandSideVector[1] = -curv*gamma*An[1];
	    GetGeometry()[0].FastGetSolutionStepValue(FORCE_X) = -curv*gamma*An[0];
	    GetGeometry()[0].FastGetSolutionStepValue(FORCE_Y) = -curv*gamma*An[1];
	  }
	  
	  /*
	  boost::numeric::ublas::bounded_matrix<double,2,2> Dn = ZeroMatrix(2,2);
	  boost::numeric::ublas::bounded_matrix<double,2,2> Dnt = ZeroMatrix(2,2);
	  Dn(0,0) = GetGeometry()[0].FastGetSolutionStepValue(VISCOUS_STRESSX_X);
	  Dn(0,1) = GetGeometry()[0].FastGetSolutionStepValue(VISCOUS_STRESSX_Y);
	  Dn(1,0) = GetGeometry()[0].FastGetSolutionStepValue(VISCOUS_STRESSY_X);
	  Dn(1,1) = GetGeometry()[0].FastGetSolutionStepValue(VISCOUS_STRESSY_Y);
	  Dnt(0,1) = Dn(1,0);
	  Dnt(1,0) = Dn(0,1);
	  array_1d<double,2> temp_vec_nv = ZeroVector(2);
	  const array_1d<double,3>& v0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
	  */
	  //double dt = rCurrentProcessInfo[DELTA_TIME];
	  //double Area;
	  //boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
	  //array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
	  //GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);
	  
	  
	  //KRATOS_WATCH( GetGeometry()[0].FastGetSolutionStepValue(FORCE_X))
	  //im->FastGetSolutionStepValue(VISCOUS_STRESSX_X) = dn_dx(0,0);
	  //im->FastGetSolutionStepValue(VISCOUS_STRESSX_Y) = dn_dx(0,1);
	  //im->FastGetSolutionStepValue(VISCOUS_STRESSY_X) = dn_dx(1,0);
	  //im->FastGetSolutionStepValue(VISCOUS_STRESSY_Y) = dn_dx(1,1);
	  /*
	  double avg_curv = GetGeometry()[0].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
	  rLeftHandSideMatrix(0,0) = avg_curv*gamma*Dn(0,0); // -(-curv*gamma*Dn(0,0)) == curv*...
	  rLeftHandSideMatrix(0,1) = avg_curv*gamma*Dn(0,1);
	  rLeftHandSideMatrix(1,0) = avg_curv*gamma*Dn(1,0);
	  rLeftHandSideMatrix(1,1) = avg_curv*gamma*Dn(1,1);
	  */
	  
	  //NOW RHS-= L u_old
	  //temp_vec_nv[0] = v0[0];
	  //temp_vec_nv[1] = v0[1];
	  //noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,temp_vec_nv);
	  //KRATOS_WATCH(prod(rLeftHandSideMatrix,temp_vec_nv))
	  //rRightHandSideVector[0] -= (rLeftHandSideMatrix(0,0)*temp_vec_nv[0] + rLeftHandSideMatrix(0,1)*temp_vec_nv[1]);
	  //rRightHandSideVector[1] -= (rLeftHandSideMatrix(1,0)*temp_vec_nv[0] + rLeftHandSideMatrix(1,1)*temp_vec_nv[1]);
	  
	  //noalias(rRightHandSideVector) += Area*dt* (prod(msWorkMatrix,temp_vec_nv)) ;
	}

        //************************************************************************************
	//************************************************************************************
	/*
	void SurfaceTension2D::CalculateLocalVelocityContribution(MatrixType& rDampMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
	KRATOS_WATCH("")
	int nodes_number = 1;
	int dim = 3;
	unsigned int matsize = nodes_number*(dim);

	if(rDampMatrix.size1() != matsize)
			rDampMatrix.resize(matsize,matsize,false); //false says not to preserve existing storage!!


	noalias(rDampMatrix) = ZeroMatrix(matsize,matsize); 



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
	void SurfaceTension2D::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
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
	void SurfaceTension2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		//noalias(rResult) = ZeroVector(2);
		rResult.resize(2,false);
		rResult[0] = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
		rResult[1] = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
// 		unsigned int number_of_nodes = GetGeometry().PointsNumber();
// 		unsigned int dim = 2;
// 		unsigned int node_size = dim;
// 
// 		
// 			if(rResult.size() != number_of_nodes*node_size)
// 				rResult.resize(number_of_nodes*node_size,false);	
// 
// 			for (unsigned int i=0;i<number_of_nodes;i++)
// 			{
// 				rResult[i*node_size] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
// 				rResult[i*node_size+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
// 
// 			}
		KRATOS_CATCH("")
		
	}

	//************************************************************************************
	//************************************************************************************
	  void SurfaceTension2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		//noalias(ElementalDofList) = ZeroVector(2);
		ElementalDofList.resize(2);
		ElementalDofList[0] = GetGeometry()[0].pGetDof(DISPLACEMENT_X);
		ElementalDofList[1] = GetGeometry()[0].pGetDof(DISPLACEMENT_Y);
// 		unsigned int number_of_nodes = 1;
// 		unsigned int dim = 2;
// 		unsigned int node_size = dim;
// 
// 			if(ElementalDofList.size() != number_of_nodes*node_size)
// 				ElementalDofList.resize(number_of_nodes*node_size);	
// 
// 			for (unsigned int i=0;i<number_of_nodes;i++)
// 			{
// 				ElementalDofList[i*node_size] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
// 				ElementalDofList[i*node_size+1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
// 			}

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//*************************************************************************************

	  void SurfaceTension2D::GetFirstDerivativesVector(Vector& values, int Step)
	{
		values.resize(2,false);
		values[0] = 0.0;
		values[1] = 0.0;
// 		const unsigned int number_of_nodes = GetGeometry().size();
// 		const unsigned int dim = 2;
// 		unsigned int MatSize = number_of_nodes * (dim );
// 
// 		if(values.size() != MatSize)   values.resize(MatSize,false);
// 		for (unsigned int i=0;i<number_of_nodes;i++)
// 		{
// 			unsigned int index = i * (dim );
// 			values[index] = 0.0;
// 			values[index + 1] = 0.0;
// 			//values[index + 2] = GetGeometry()[i].GetSolutionStepValue(PRESSURE,Step);
// 
// 		}
	}
	//************************************************************************************
	//************************************************************************************
	  void SurfaceTension2D::GetSecondDerivativesVector(Vector& values, int Step)
	{
		values.resize(2,false);
		values[0] = 0.0;
		values[1] = 0.0;
// 		const unsigned int number_of_nodes = GetGeometry().size();
// 		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
// 		unsigned int MatSize = number_of_nodes * (dim);
// 		if(values.size() != MatSize) values.resize(MatSize,false);
// 		for (unsigned int i=0;i<number_of_nodes;i++)
// 		{
// 			unsigned int index = i * (dim );
// 			values[index] = 0.0;
// 			values[index + 1] = 0.0;
// 			//values[index + 2] = 0.0;
// 		}
	}
	//************************************************************************************
	//************************************************************************************

	  
} // Namespace Kratos


