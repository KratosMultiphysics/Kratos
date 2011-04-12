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
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2008-04-30 14:34:17 $
//   Revision:            $Revision: 1.4 $
//
//
 
//#define GRADPN_FORM
//#define STOKES

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_fluid3D.h"
#include "utilities/math_utils.h"
#include "ULF_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{


	//************************************************************************************
	//************************************************************************************
	UpdatedLagrangianFluid3D::UpdatedLagrangianFluid3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	UpdatedLagrangianFluid3D::UpdatedLagrangianFluid3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer UpdatedLagrangianFluid3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		KRATOS_TRY

		//KRATOS_WATCH("3D CONSTRUCTOR PAOLO IS THE BEST")
		return Element::Pointer(new UpdatedLagrangianFluid3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	UpdatedLagrangianFluid3D::~UpdatedLagrangianFluid3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		//aux variables
		boost::numeric::ublas::bounded_matrix<double,6,12> msB = ZeroMatrix(6,12);
		boost::numeric::ublas::bounded_matrix<double,6,6> ms_constitutive_matrix;
		boost::numeric::ublas::bounded_matrix<double,6,12> ms_temp;
	
		//the one below - not changed yet (coz I dont know what is stored in this vector)
		array_1d<double,6> ms_temp_vec;
	
		boost::numeric::ublas::bounded_matrix<double,4,3> msDN_Dx;
	  	array_1d<double,4> msN; //dimension = number of nodes

		const double density = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
							GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
							GetGeometry()[2].FastGetSolutionStepValue(DENSITY) + 
							GetGeometry()[3].FastGetSolutionStepValue(DENSITY));
		double K = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS)+
							GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS) +
							GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS)+
							GetGeometry()[3].FastGetSolutionStepValue(BULK_MODULUS));
		K *= density;

		unsigned int dim = 3;
		unsigned int number_of_nodes = 4;

		if(rLeftHandSideMatrix.size1() != 12)
			rLeftHandSideMatrix.resize(12,12,false);

		if(rRightHandSideVector.size() != 12)
			rRightHandSideVector.resize(12,false);

		//calculate current volume
		double current_volume;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_Dx, msN, current_volume);


		//writing the body force
		const array_1d<double,3>& body_force = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+
							GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +
							GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE) +
							GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE));

		for(unsigned int i = 0; i<number_of_nodes; i++)
		{
			rRightHandSideVector[i*3] = body_force[0]* density * mA0 * 0.25;
			rRightHandSideVector[i*3+1] = body_force[1] * density * mA0 * 0.25;
			rRightHandSideVector[i*3+2] = body_force[2] * density * mA0 * 0.25;
			
		}

		//VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
		// rLeftHandSideMatrix += Laplacian * nu;
		//filling matrix B
		for(unsigned int i = 0; i<number_of_nodes; i++)
		{
			unsigned int start = dim*i;

			msB(0,start) =	msDN_Dx(i,0);
			msB(1,start+1)=	msDN_Dx(i,1);
			msB(2,start+2)= msDN_Dx(i,2);
			msB(3,start) =	msDN_Dx(i,1);		msB(3,start+1) = msDN_Dx(i,0);
			msB(4,start) =	msDN_Dx(i,2);		msB(4,start+2) = msDN_Dx(i,0);
			msB(5,start+1)= msDN_Dx(i,2);		msB(5,start+2) = msDN_Dx(i,1);
		}
				
				
		//constitutive tensor
		ms_constitutive_matrix(0,0) = K; ms_constitutive_matrix(0,1) = K ;	ms_constitutive_matrix(0,2) = K;
		ms_constitutive_matrix(0,3)= 0.0;	ms_constitutive_matrix(0,4) = 0.0 ;	ms_constitutive_matrix(0,5) = 0.0;
				
		ms_constitutive_matrix(1,0) = K; 	ms_constitutive_matrix(1,1) = K;	ms_constitutive_matrix(1,2) = K;
		ms_constitutive_matrix(1,3) = 0.0; 	ms_constitutive_matrix(1,4) = 0.0;	ms_constitutive_matrix(1,5) = 0.0;

		ms_constitutive_matrix(2,0) = K;	ms_constitutive_matrix(2,1) = K;	ms_constitutive_matrix(2,2) = K;
		ms_constitutive_matrix(2,3) = 0.0;	ms_constitutive_matrix(2,4) = 0.0;	ms_constitutive_matrix(2,5) = 0.0;
		
		ms_constitutive_matrix(3,0) = 0.0;	ms_constitutive_matrix(3,1) = 0.0;	ms_constitutive_matrix(3,2) = 0.0;
		ms_constitutive_matrix(3,3) = 0.0;	ms_constitutive_matrix(3,4) = 0.0;	ms_constitutive_matrix(3,5) = 0.0;

		ms_constitutive_matrix(4,0) = 0.0;	ms_constitutive_matrix(4,1) = 0.0;	ms_constitutive_matrix(4,2) = 0.0;
		ms_constitutive_matrix(4,3) = 0.0;	ms_constitutive_matrix(4,4) = 0.0;	ms_constitutive_matrix(4,5) = 0.0;
		
		ms_constitutive_matrix(5,0) = 0.0;	ms_constitutive_matrix(5,1) = 0.0;	ms_constitutive_matrix(5,2) = 0.0;
		ms_constitutive_matrix(5,3) = 0.0;	ms_constitutive_matrix(5,4) = 0.0;	ms_constitutive_matrix(5,5) = 0.0;
		
		
		//calculating viscous contributions
		ms_temp = prod( ms_constitutive_matrix , msB);
		noalias(rLeftHandSideMatrix) = prod( trans(msB) , ms_temp);
		rLeftHandSideMatrix *= -current_volume;

		//get the old pressure
		double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
		double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
		double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
		double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE,1);

		//adding pressure gradient 
		double pavg = p0old + p1old + p2old + p3old; //calculate old pressure over the element
		pavg *= 0.25 * current_volume;

		//add to pavg the increment of pressure inside the time step
		double dp = K * (current_volume - mA0);
		pavg += dp;

		rRightHandSideVector[0] += msDN_Dx(0,0)*pavg;
		rRightHandSideVector[1] += msDN_Dx(0,1)*pavg;
		rRightHandSideVector[2] += msDN_Dx(0,2)*pavg;
		
		rRightHandSideVector[3] += msDN_Dx(1,0)*pavg;
		rRightHandSideVector[4] += msDN_Dx(1,1)*pavg;
		rRightHandSideVector[5] += msDN_Dx(1,2)*pavg;

		rRightHandSideVector[6] += msDN_Dx(2,0)*pavg;
		rRightHandSideVector[7] += msDN_Dx(2,1)*pavg;
		rRightHandSideVector[8] += msDN_Dx(2,2)*pavg;

		rRightHandSideVector[9]	 += msDN_Dx(3,0)*pavg;
		rRightHandSideVector[10] += msDN_Dx(3,1)*pavg;
		rRightHandSideVector[11] += msDN_Dx(3,2)*pavg;

		
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		
		//aux variables
		boost::numeric::ublas::bounded_matrix<double,6,12> msB = ZeroMatrix(6,12);
		boost::numeric::ublas::bounded_matrix<double,6,6> ms_constitutive_matrix;
		boost::numeric::ublas::bounded_matrix<double,6,12> ms_temp;
	
		//the one below - not changed yet (coz I dont know what is stored in this vector)
		array_1d<double,6> ms_temp_vec;
	
		boost::numeric::ublas::bounded_matrix<double,4,3> msDN_Dx;
	  	array_1d<double,4> msN; //dimension = number of nodes
		
		const double& density = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
							GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
							GetGeometry()[2].FastGetSolutionStepValue(DENSITY)+
							GetGeometry()[3].FastGetSolutionStepValue(DENSITY));	
		double K = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS)+
							GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS) +
							GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS)+
							GetGeometry()[3].FastGetSolutionStepValue(BULK_MODULUS));
		K *= density;

	//	KRATOS_ERROR(std::logic_error,"not goooood","");
		if(rRightHandSideVector.size() != 12)
			rRightHandSideVector.resize(12,false);
		unsigned int number_of_nodes = GetGeometry().size();

		//calculate current area
		double current_volume;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_Dx, msN, current_volume);

		//writing the body force
		const array_1d<double,3>& body_force = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+
							GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +
							GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE)+
							GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE));
	

		for(unsigned int i = 0; i<number_of_nodes; i++)
		{
			rRightHandSideVector[i*3] = body_force[0]* density * mA0 * 0.25;
			rRightHandSideVector[i*3+1] = body_force[1] * density * mA0 * 0.25;
			rRightHandSideVector[i*3+2] = body_force[2] * density * mA0 * 0.25;
		}


		//get the old pressure
		double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
		double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
		double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
		double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE,1);

		//adding pressure gradient 
		double pavg = p0old + p1old + p2old + p3old; //calculate old pressure over the element
		pavg *= 0.25 * current_volume;

		//add to pavg the increment of pressure inside the time step
		double dp = K * (current_volume - mA0);
		pavg += dp;

		rRightHandSideVector[0] += msDN_Dx(0,0)*pavg;
		rRightHandSideVector[1] += msDN_Dx(0,1)*pavg;
		rRightHandSideVector[2] += msDN_Dx(0,2)*pavg;
		
		rRightHandSideVector[3] += msDN_Dx(1,0)*pavg;
		rRightHandSideVector[4] += msDN_Dx(1,1)*pavg;
		rRightHandSideVector[5] += msDN_Dx(1,2)*pavg;

		rRightHandSideVector[6] += msDN_Dx(2,0)*pavg;
		rRightHandSideVector[7] += msDN_Dx(2,1)*pavg;
		rRightHandSideVector[8] += msDN_Dx(2,2)*pavg;

		rRightHandSideVector[9]	 += msDN_Dx(3,0)*pavg;
		rRightHandSideVector[10] += msDN_Dx(3,1)*pavg;
		rRightHandSideVector[11] += msDN_Dx(3,2)*pavg;
  	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid3D::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		//aux variables
		boost::numeric::ublas::bounded_matrix<double,6,12> msB = ZeroMatrix(6,12);
		boost::numeric::ublas::bounded_matrix<double,6,6> ms_constitutive_matrix;
		boost::numeric::ublas::bounded_matrix<double,6,12> ms_temp;
	
		//the one below - not changed yet (coz I dont know what is stored in this vector)
		array_1d<double,6> ms_temp_vec;
	
		boost::numeric::ublas::bounded_matrix<double,4,3> msDN_Dx;
	  	array_1d<double,4> msN; //dimension = number of nodes
		
		const double& density = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
							GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
							GetGeometry()[2].FastGetSolutionStepValue(DENSITY)+	
							GetGeometry()[3].FastGetSolutionStepValue(DENSITY));	
		//lumped
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int NumberOfNodes = GetGeometry().size();
		unsigned int MatSize = dimension * NumberOfNodes;
		
		if(rMassMatrix.size1() != MatSize)
			rMassMatrix.resize(MatSize,MatSize,false);

		noalias(rMassMatrix) = ZeroMatrix(MatSize,MatSize);

		double nodal_mass = mA0 * density * 0.25;

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
	void UpdatedLagrangianFluid3D::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		//aux variables
		boost::numeric::ublas::bounded_matrix<double,6,12> msB = ZeroMatrix(6,12);
		boost::numeric::ublas::bounded_matrix<double,6,6> ms_constitutive_matrix;
		boost::numeric::ublas::bounded_matrix<double,6,12> ms_temp;
	
		//the one below - not changed yet (coz I dont know what is stored in this vector)
		array_1d<double,6> ms_temp_vec;
	
		boost::numeric::ublas::bounded_matrix<double,4,3> msDN_Dx;
	  	array_1d<double,4> msN; //dimension = number of nodes

		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();
		
		if(rDampMatrix.size1() != 12)
			rDampMatrix.resize(12,12,false);

		
		//getting data for the given geometry
		double current_volume;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_Dx, msN, current_volume);

		//getting properties
		const double& nu = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY)+
							GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY) +
							GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY)+
							GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY));
		
		//double nu = GetProperties()[VISCOSITY];
		//double density = GetProperties()[DENSITY];
		const double& density = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
							GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
							GetGeometry()[2].FastGetSolutionStepValue(DENSITY)+
							GetGeometry()[3].FastGetSolutionStepValue(DENSITY));
		//VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
		// rLeftHandSideMatrix += Laplacian * nu;
		//filling matrix B
		for(unsigned int i = 0; i<number_of_nodes; i++)
		{
			unsigned int start = dim*i;

			msB(0,start) =	msDN_Dx(i,0); 
			msB(1,start+1)=	msDN_Dx(i,1);
			msB(2,start+2)= msDN_Dx(i,2);
			msB(3,start) =	msDN_Dx(i,1);		msB(3,start+1) = msDN_Dx(i,0);
			msB(4,start) =	msDN_Dx(i,2);		msB(4,start+2) = msDN_Dx(i,0);
			msB(5,start+1)= msDN_Dx(i,2);		msB(5,start+2) = msDN_Dx(i,1);
		}			
		
		const double& a = nu*density;
		//constitutive tensor
		ms_constitutive_matrix(0,0) = (4.0/3.0)*a;	ms_constitutive_matrix(0,1) = -2.0/3.0*a;	ms_constitutive_matrix(0,2) = -2.0/3.0*a;
		ms_constitutive_matrix(0,3) = 0.0;			ms_constitutive_matrix(0,4) = 0.0;			ms_constitutive_matrix(0,5) = 0.0;
		
		ms_constitutive_matrix(1,0) = -2.0/3.0*a; 	ms_constitutive_matrix(1,1) = 4.0/3.0*a;	ms_constitutive_matrix(1,2) = -2.0/3.0*a;
		ms_constitutive_matrix(1,3) = 0.0;		 	ms_constitutive_matrix(1,4) = 0.0;			ms_constitutive_matrix(1,5) = 0.0;

		ms_constitutive_matrix(2,0) = -2.0/3.0*a;	ms_constitutive_matrix(2,1) = -2.0/3.0*a;	ms_constitutive_matrix(2,2) = 4.0/3.0*a;
		ms_constitutive_matrix(2,3) = 0.0;			ms_constitutive_matrix(2,4) = 0.0;			ms_constitutive_matrix(2,5) = 0.0;
		
		ms_constitutive_matrix(3,0) = 0.0;			ms_constitutive_matrix(3,1) = 0.0;			ms_constitutive_matrix(3,2) = 0.0;
		ms_constitutive_matrix(3,3) = a;			ms_constitutive_matrix(3,4) = 0.0;			ms_constitutive_matrix(3,5) = 0.0;
		
		ms_constitutive_matrix(4,0) = 0.0;			ms_constitutive_matrix(4,1) = 0.0;			ms_constitutive_matrix(4,2) = 0.0;
		ms_constitutive_matrix(4,3) = 0.0;			ms_constitutive_matrix(4,4) = a;			ms_constitutive_matrix(4,5) = 0.0;
		
		ms_constitutive_matrix(5,0) = 0.0;			ms_constitutive_matrix(5,1) = 0.0;			ms_constitutive_matrix(5,2) = 0.0;
		ms_constitutive_matrix(5,3) = 0.0;			ms_constitutive_matrix(5,4) = 0.0;			ms_constitutive_matrix(5,5) = a;
		
		//calculating viscous contributions
		ms_temp = prod( ms_constitutive_matrix , msB);
		noalias(rDampMatrix) = prod( trans(msB) , ms_temp);

		rDampMatrix *= current_volume;
		KRATOS_CATCH("")
	}	
	  
	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void UpdatedLagrangianFluid3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		//save original Volume
		mA0 = GeometryUtils::CalculateVolume3D(GetGeometry());

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 3;

		if(rResult.size() != number_of_nodes*dim)
			rResult.resize(number_of_nodes*dim,false);	

		int guess_pos = GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			rResult[i*dim] = GetGeometry()[i].GetDof(DISPLACEMENT_X,guess_pos).EquationId();
			rResult[i*dim+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,guess_pos+1).EquationId();
			rResult[i*dim+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z,guess_pos+2).EquationId();
		}

	}

	//************************************************************************************
	//************************************************************************************
	  void UpdatedLagrangianFluid3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 3;

		if(ElementalDofList.size() != number_of_nodes*dim)
			ElementalDofList.resize(number_of_nodes*dim);	

		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			ElementalDofList[i*dim] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
			ElementalDofList[i*dim+1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
			ElementalDofList[i*dim+2] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
		}
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid3D::GetValuesVector(Vector& values, int Step)
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
			values[index + 2] = disp[2];
		}
	}
	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid3D::GetFirstDerivativesVector(Vector& values, int Step)
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
			values[index + 2] = vel[2];
		}
	}
	//************************************************************************************
	//************************************************************************************
	  void UpdatedLagrangianFluid3D::GetSecondDerivativesVector(Vector& values, int Step)
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
			values[index + 2] = acc[2];
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void  UpdatedLagrangianFluid3D::Calculate(const Variable<double >& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
	  {
		//aux variables
		boost::numeric::ublas::bounded_matrix<double,6,12> msB = ZeroMatrix(6,12);
		boost::numeric::ublas::bounded_matrix<double,6,6> ms_constitutive_matrix;
		boost::numeric::ublas::bounded_matrix<double,6,12> ms_temp;
	
		//the one below - not changed yet (coz I dont know what is stored in this vector)
		array_1d<double,6> ms_temp_vec;
	
		boost::numeric::ublas::bounded_matrix<double,4,3> msDN_Dx;
	  	array_1d<double,4> msN; //dimension = number of nodes
		if(rVariable == PRESSURE)
		{
			const double density = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
								GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
								GetGeometry()[2].FastGetSolutionStepValue(DENSITY)+
								GetGeometry()[3].FastGetSolutionStepValue(DENSITY));
			//const double& density = GetProperties()[DENSITY];
			double K = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS)+
								GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS) +
								GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS)+
								GetGeometry()[3].FastGetSolutionStepValue(BULK_MODULUS));

			K *= density;
	
			double current_volume = GeometryUtils::CalculateVolume3D(GetGeometry());

			//save in msN the old pressure 
			msN[0] = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
			msN[1] = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
			msN[2] = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
			msN[3] = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE,1);
			double diag_term = 0.1*current_volume; 
			double out_term = 0.05*current_volume;
			
	
			//add to pavg the increment of pressure inside the time step
			double dp_volume = K * 0.25*(current_volume - mA0);
	
	//		GetGeometry()[0].FastGetSolutionStepValue(PRESSURE) += dp_volume;
	//		GetGeometry()[1].FastGetSolutionStepValue(PRESSURE) += dp_volume;
	//		GetGeometry()[2].FastGetSolutionStepValue(PRESSURE) += dp_volume;
	//		GetGeometry()[3].FastGetSolutionStepValue(PRESSURE) += dp_volume;
			//p = Mconsistent*pold + dp
			GetGeometry()[0].FastGetSolutionStepValue(PRESSURE) += dp_volume + diag_term*msN[0] + out_term*msN[1]  + out_term*msN[2]  + out_term*msN[3];
			GetGeometry()[1].FastGetSolutionStepValue(PRESSURE) += dp_volume + out_term*msN[0]  + diag_term*msN[1] + out_term*msN[2]  + out_term*msN[3] ;
			GetGeometry()[2].FastGetSolutionStepValue(PRESSURE) += dp_volume + out_term*msN[0]  + out_term*msN[1]  + diag_term*msN[2] + out_term*msN[3];
			GetGeometry()[3].FastGetSolutionStepValue(PRESSURE) += dp_volume + out_term*msN[0]  + out_term*msN[1]  + out_term*msN[2]  + diag_term*msN[3] ;
		}
		else if(rVariable == IS_FLUID)
		{
			GetGeometry()[0].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;
			GetGeometry()[1].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;
			GetGeometry()[2].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;
			GetGeometry()[3].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;

			double current_volume;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_Dx, msN, current_volume);
		//KRATOS_WATCH(current_volume)
		//KRATOS_WATCH(GetGeometry()[0].FastGetSolutionStepValue(PRESSURE))
		//KRATOS_WATCH(msDN_Dx)
		GetGeometry()[0].FastGetSolutionStepValue(FORCE_X)+=0.25*current_volume*msDN_Dx(0,0)*GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
		GetGeometry()[0].FastGetSolutionStepValue(FORCE_Y)+=0.25*current_volume*msDN_Dx(0,1)*GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
		GetGeometry()[0].FastGetSolutionStepValue(FORCE_Z)+=0.25*current_volume*msDN_Dx(0,2)*GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);

		GetGeometry()[1].FastGetSolutionStepValue(FORCE_X)+=0.25*current_volume*msDN_Dx(1,0)*GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
		GetGeometry()[1].FastGetSolutionStepValue(FORCE_Y)+=0.25*current_volume*msDN_Dx(1,1)*GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
		GetGeometry()[1].FastGetSolutionStepValue(FORCE_Z)+=0.25*current_volume*msDN_Dx(1,2)*GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);

		GetGeometry()[2].FastGetSolutionStepValue(FORCE_X)+=0.25*current_volume*msDN_Dx(2,0)*GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
		GetGeometry()[2].FastGetSolutionStepValue(FORCE_Y)+=0.25*current_volume*msDN_Dx(2,1)*GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
		GetGeometry()[2].FastGetSolutionStepValue(FORCE_Z)+=0.25*current_volume*msDN_Dx(2,2)*GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);

		GetGeometry()[3].FastGetSolutionStepValue(FORCE_X)+=0.25*current_volume*msDN_Dx(3,0)*GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
		GetGeometry()[3].FastGetSolutionStepValue(FORCE_Y)+=0.25*current_volume*msDN_Dx(3,1)*GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
		GetGeometry()[3].FastGetSolutionStepValue(FORCE_Z)+=0.25*current_volume*msDN_Dx(3,2)*GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);

		//KRATOS_WATCH(GetGeometry()[0].FastGetSolutionStepValue(FORCE_X))
		}

	}

} // Namespace Kratos


