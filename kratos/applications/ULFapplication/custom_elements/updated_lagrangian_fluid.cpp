/*
==============================================================================
KratosULFApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
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
//   Date:                $Date: 2007-05-16 13:59:01 $
//   Revision:            $Revision: 1.6 $
//
//
 
//#define GRADPN_FORM
//#define STOKES

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_fluid.h"
#include "utilities/math_utils.h"
#include "ULF_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
	//static variables
	

	
	//************************************************************************************
	//************************************************************************************
	UpdatedLagrangianFluid::UpdatedLagrangianFluid(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	UpdatedLagrangianFluid::UpdatedLagrangianFluid(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer UpdatedLagrangianFluid::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		KRATOS_TRY

		return Element::Pointer(new UpdatedLagrangianFluid(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	UpdatedLagrangianFluid::~UpdatedLagrangianFluid()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		boost::numeric::ublas::bounded_matrix<double,3,6> msB;
		boost::numeric::ublas::bounded_matrix<double,3,3> ms_constitutive_matrix;
		boost::numeric::ublas::bounded_matrix<double,3,6> ms_temp;

		boost::numeric::ublas::bounded_matrix<double,3,2> UpdatedLagrangianFluid::msDN_Dx;
	  	array_1d<double,3> msN; //dimension = number of nodes
	
		static array_1d<double,6> ms_temp_vec;
		const double& density = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
							GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
							GetGeometry()[2].FastGetSolutionStepValue(DENSITY));
		double K = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS)+
							GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS) +
							GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS));
//		double K = GetProperties()[BULK_MODULUS];
		K *= density;

		unsigned int dim = 2;
		unsigned int number_of_nodes = 3;

		if(rLeftHandSideMatrix.size1() != 6)
			rLeftHandSideMatrix.resize(6,6,false);

		if(rRightHandSideVector.size() != 6)
			rRightHandSideVector.resize(6,false);

		//calculate current area
		double current_area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_Dx, msN, current_area);


		//writing the body force
		const array_1d<double,3>& body_force = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+
							GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +
							GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE));

		for(unsigned int i = 0; i<number_of_nodes; i++)
		{
			rRightHandSideVector[i*2] = body_force[0]* density * mA0 * 0.3333333333333;
			rRightHandSideVector[i*2+1] = body_force[1] * density * mA0 * 0.3333333333333;
//			rRightHandSideVector[i*2] = body_force[0]* density * current_area * 0.3333333333333;
			//rRightHandSideVector[i*2+1] = body_force[1] * density * current_area * 0.3333333333333;
		}

		//VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
		// rLeftHandSideMatrix += Laplacian * nu;
		//filling matrix B
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = dim*i;
			msB(0,index+0)=msDN_Dx(i,0);					msB(0,index+1)= 0.0;
			msB(1,index+0)=0.0;								msB(1,index+1)= msDN_Dx(i,1);
			msB(2,index+0)= msDN_Dx(i,1);					msB(2,index+1)= msDN_Dx(i,0); 
		}			
				
		//constitutive tensor
		ms_constitutive_matrix(0,0) = K;	ms_constitutive_matrix(0,1) = K ;	ms_constitutive_matrix(0,2) = 0.0;
		ms_constitutive_matrix(1,0) = K; 	ms_constitutive_matrix(1,1) = K;	ms_constitutive_matrix(1,2) = 0.0;
		ms_constitutive_matrix(2,0) = 0.0;	ms_constitutive_matrix(2,1) = 0.0;	ms_constitutive_matrix(2,2) = 0.0;
			
		//calculating viscous contributions
		ms_temp = prod( ms_constitutive_matrix , msB);
		noalias(rLeftHandSideMatrix) = prod( trans(msB) , ms_temp);
		rLeftHandSideMatrix *= -current_area;

		//get the old pressure
		double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
		double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
		double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);

		//adding pressure gradient 
		double pavg = p0old + p1old + p2old; //calculate old pressure over the element
		pavg *= 0.33333333333333333333333 * current_area;

		//add to pavg the increment of pressure inside the time step
		double dp = K * (current_area - mA0);
		pavg += dp;

		rRightHandSideVector[0] += msDN_Dx(0,0)*pavg;
		rRightHandSideVector[1] += msDN_Dx(0,1)*pavg;
		rRightHandSideVector[2] += msDN_Dx(1,0)*pavg;
		rRightHandSideVector[3] += msDN_Dx(1,1)*pavg;
		rRightHandSideVector[4] += msDN_Dx(2,0)*pavg;
		rRightHandSideVector[5] += msDN_Dx(2,1)*pavg;
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{

		boost::numeric::ublas::bounded_matrix<double,3,6> msB;
		boost::numeric::ublas::bounded_matrix<double,3,3> ms_constitutive_matrix;
		boost::numeric::ublas::bounded_matrix<double,3,6> ms_temp;

		boost::numeric::ublas::bounded_matrix<double,3,2> UpdatedLagrangianFluid::msDN_Dx;
	  	array_1d<double,3> msN; //dimension = number of nodes		
		const double& density = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
							GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
							GetGeometry()[2].FastGetSolutionStepValue(DENSITY));		
		//const double& density = GetProperties()[DENSITY];
		double K = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS)+
							GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS) +
							GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS));
//		double K = GetProperties()[BULK_MODULUS];
		K *= density;

	//	KRATOS_ERROR(std::logic_error,"not goooood","");
		if(rRightHandSideVector.size() != 6)
			rRightHandSideVector.resize(6,false);
		unsigned int number_of_nodes = GetGeometry().size();

		//calculate current area
		double current_area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_Dx, msN, current_area);

		//writing the body force
		const array_1d<double,3>& body_force = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+
							GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +
							GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE));
		//const array_1d<double,3>& body_force = GetProperties()[BODY_FORCE];
		for(unsigned int i = 0; i<number_of_nodes; i++)
		{
			rRightHandSideVector[i*2] = body_force[0]* density * mA0 * 0.3333333333333;
			rRightHandSideVector[i*2+1] = body_force[1] * density * mA0 * 0.3333333333333;
			//rRightHandSideVector[i*2] = body_force[0]* density * current_area * 0.3333333333333;
			//rRightHandSideVector[i*2+1] = body_force[1] * density * current_area * 0.3333333333333;
		}





		//get the old pressure
		double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
		double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
		double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);

		//adding pressure gradient 
		double pavg = p0old + p1old + p2old; //calculate old pressure over the element
		pavg *= 0.33333333333333333333333 * current_area;

		//add to pavg the increment of pressure inside the time step
		double dp = K * (current_area - mA0);
		pavg += dp;

		rRightHandSideVector[0] += msDN_Dx(0,0)*pavg;
		rRightHandSideVector[1] += msDN_Dx(0,1)*pavg;
		rRightHandSideVector[2] += msDN_Dx(1,0)*pavg;
		rRightHandSideVector[3] += msDN_Dx(1,1)*pavg;
		rRightHandSideVector[4] += msDN_Dx(2,0)*pavg;
		rRightHandSideVector[5] += msDN_Dx(2,1)*pavg;
  	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		
		const double& density = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
							GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
							GetGeometry()[2].FastGetSolutionStepValue(DENSITY));	
		//lumped
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int NumberOfNodes = GetGeometry().size();
		unsigned int MatSize = dimension * NumberOfNodes;
		
		if(rMassMatrix.size1() != MatSize)
			rMassMatrix.resize(MatSize,MatSize,false);

		noalias(rMassMatrix) = ZeroMatrix(MatSize,MatSize);

		double nodal_mass = mA0 * density * 0.333333333333333333;

		for(unsigned int i=0; i<NumberOfNodes; i++)
		{
			for(unsigned int j=0; j<dimension; j++)
			{
				unsigned int index = i*dimension + j;
				rMassMatrix(index,index) = nodal_mass;
			}
		}

KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		boost::numeric::ublas::bounded_matrix<double,3,6> msB;
		boost::numeric::ublas::bounded_matrix<double,3,3> ms_constitutive_matrix;
		boost::numeric::ublas::bounded_matrix<double,3,6> ms_temp;

		boost::numeric::ublas::bounded_matrix<double,3,2> UpdatedLagrangianFluid::msDN_Dx;
	  	array_1d<double,3> msN; //dimension = number of nodes

		unsigned int NumberOfNodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();
		
		if(rDampMatrix.size1() != 6)
			rDampMatrix.resize(6,6,false);

		
		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_Dx, msN, Area);

		//getting properties
		const double& nu = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY)+
							GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY) +
							GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY));
		//double nu = GetProperties()[VISCOSITY];
		//double density = GetProperties()[DENSITY];
		const double& density = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
							GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
							GetGeometry()[2].FastGetSolutionStepValue(DENSITY));
		//VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
		// rLeftHandSideMatrix += Laplacian * nu;
		//filling matrix B
		for (unsigned int i=0;i<NumberOfNodes;i++)
		{
			unsigned int index = dim*i;
			msB(0,index+0)=msDN_Dx(i,0);							msB(0,index+1)= 0.0;
			msB(1,index+0)=0.0;								msB(1,index+1)= msDN_Dx(i,1);
			msB(2,index+0)= msDN_Dx(i,1);							msB(2,index+1)= msDN_Dx(i,0); 
		}			
				
		//constitutive tensor
		ms_constitutive_matrix(0,0) = (4.0/3.0)*nu*density;	ms_constitutive_matrix(0,1) = -2.0/3.0*nu*density;		ms_constitutive_matrix(0,2) = 0.0;
		ms_constitutive_matrix(1,0) = -2.0/3.0*nu*density; 	ms_constitutive_matrix(1,1) = 4.0/3.0*nu*density;		ms_constitutive_matrix(1,2) = 0.0;
		ms_constitutive_matrix(2,0) = 0.0;					ms_constitutive_matrix(2,1) = 0.0;						ms_constitutive_matrix(2,2) = nu*density;
			
		//calculating viscous contributions
		ms_temp = prod( ms_constitutive_matrix , msB);
		noalias(rDampMatrix) = prod( trans(msB) , ms_temp);

		rDampMatrix *= Area;

		KRATOS_CATCH("")
	}	
	  
	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void UpdatedLagrangianFluid::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		//save original Area
		mA0 = GeometryUtils::CalculateVolume2D(GetGeometry());

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;

		if(rResult.size() != number_of_nodes*dim)
			rResult.resize(number_of_nodes*dim,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			rResult[i*dim] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			rResult[i*dim+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
		}
		
	}

	//************************************************************************************
	//************************************************************************************
	  void UpdatedLagrangianFluid::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;

		if(ElementalDofList.size() != number_of_nodes*dim)
			ElementalDofList.resize(number_of_nodes*dim);	

		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			ElementalDofList[i*dim] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
			ElementalDofList[i*dim+1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
		}
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::GetValuesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes*dim;
		if(values.size() != MatSize)	values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i*dim;
			const array_1d<double,3>& disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,Step);
			values[index] = disp[0];
			values[index + 1] = disp[1];
		}
	}
	
	
	//************************************************************************************
	//************************************************************************************
	  void UpdatedLagrangianFluid::GetFirstDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes*dim;
		if(values.size() != MatSize)   values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i*dim;
			const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,Step);
			values[index] = vel[0];
			values[index + 1] = vel[1];
		}
	}
	//************************************************************************************
	//************************************************************************************
	  void UpdatedLagrangianFluid::GetSecondDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes*dim;
		if(values.size() != MatSize) values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i*dim;
			const array_1d<double,3>& acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION,Step);
			values[index] = acc[0];
			values[index + 1] = acc[1];
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void  UpdatedLagrangianFluid::Calculate(const Variable<double >& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
	  {
		if(rVariable == PRESSURE)
		{
			const double density = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
								GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
								GetGeometry()[2].FastGetSolutionStepValue(DENSITY));

			//const double& density = GetProperties()[DENSITY];
			double K = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS)+
								GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS) +
								GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS));
	//		double K = GetProperties()[BULK_MODULUS];
			K *= density;
	
			double current_area = GeometryUtils::CalculateVolume2D(GetGeometry());

			//save in msN the old pressure 
			msN[0] = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
			msN[1] = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
			msN[2] = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
			double diag_term = current_area/6.0; 
			double out_term = current_area/12.0;

			//add to pavg the increment of pressure inside the time step
			double dp_area = K * 0.3333333333333333333*(current_area - mA0);
	
			//P = Mconsistent*pold + dp
	//		GetGeometry()[0].FastGetSolutionStepValue(PRESSURE) += dp_area;
	//		GetGeometry()[1].FastGetSolutionStepValue(PRESSURE) += dp_area;
	//		GetGeometry()[2].FastGetSolutionStepValue(PRESSURE) += dp_area;
			GetGeometry()[0].FastGetSolutionStepValue(PRESSURE) += dp_area + diag_term*msN[0] + out_term*msN[1]  + out_term*msN[2] ;
			GetGeometry()[1].FastGetSolutionStepValue(PRESSURE) += dp_area + out_term*msN[0] + diag_term*msN[1]  + out_term*msN[2] ;
			GetGeometry()[2].FastGetSolutionStepValue(PRESSURE) += dp_area + out_term*msN[0] + out_term*msN[1]  + diag_term*msN[2] ;
		}
		else if (rVariable == IS_FLUID)
		{
			GetGeometry()[0].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;
			GetGeometry()[1].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;
			GetGeometry()[2].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;

		}


	}

} // Namespace Kratos


