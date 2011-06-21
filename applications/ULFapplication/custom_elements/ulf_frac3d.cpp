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
//   Date:                $Date: 2009-01-15 14:50:24 $
//   Revision:            $Revision: 1.3 $
//
//
 
//#define GRADPN_FORM
//#define STOKES

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/ulf_frac3d.h"
#include "utilities/math_utils.h"
#include "ULF_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
	//static variables
	
	namespace UlfFrac3D_auxiliaries
	{

	boost::numeric::ublas::bounded_matrix<double,6,12> msB; //
	boost::numeric::ublas::bounded_matrix<double,6,6> ms_constitutive_matrix; //
	boost::numeric::ublas::bounded_matrix<double,6,12> ms_temp; //
	array_1d<double,6> ms_temp_vec;
	
		boost::numeric::ublas::bounded_matrix<double,4,4> msWorkMatrix = ZeroMatrix(4,4); //
		boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX = ZeroMatrix(4,3); //

		array_1d<double,4> ms_temp_vec_np = ZeroVector(4); //dimension = number of nodes //
		array_1d<double,3> ms_aux0 = ZeroVector(3); //dimension = number of space dimensions //
		array_1d<double,4> ms_aux1 = ZeroVector(4); //dimension = number of nodes //
		array_1d<double,3> ms_vel_gauss = ZeroVector(3); //dimesion coincides with space dimension//
		array_1d<double,4> msN = ZeroVector(4);//
	} 
	using  namespace UlfFrac3D_auxiliaries;

	//************************************************************************************
	//************************************************************************************
	UlfFrac3D::UlfFrac3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	UlfFrac3D::UlfFrac3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer UlfFrac3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		KRATOS_TRY
		return Element::Pointer(new UlfFrac3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	UlfFrac3D::~UlfFrac3D()
	{
	}
	//************************************************************************************
	//************************************************************************************
	void UlfFrac3D::CalculateLumpedMass()		
	{
	//note that for the compressible case, rho will be also a variable
	
	const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
	const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
	const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
	const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

	//double Area;
	//GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);
	double Volume = GeometryUtils::CalculateVolume3D(GetGeometry());
	double lumped_mass_fac = Volume * 0.25;
	//filling in the diagonal of the lumped mass matrix,  (later I can change it to vector...)
	GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho0;
	GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho1;	
	GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho2;
	GetGeometry()[3].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho3;				
	}

	//************************************************************************************
	//************************************************************************************
			
	void UlfFrac3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		//KRATOS_WATCH("Current FRACTIONAL STEP ISSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS")
		
		int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

		//KRATOS_WATCH(FractionalStepNumber)

		if(FractionalStepNumber == 1) //first step of the fractional step solution
		{
			//int ComponentIndex = FractionalStepNumber - 1;
			//Stage1(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo, ComponentIndex);
			VelocityStep(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
		}
		else if (FractionalStepNumber == 2)//second step of the fractional step solution
		{
			PressureStep(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
		}


		
		KRATOS_CATCH("")
	}
	

	//************************************************************************************
	//************************************************************************************
	void UlfFrac3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
	KRATOS_TRY
	//KRATOS_WATCH("Current FRACTIONAL STEP IS (CalcRightHandSide)!!!!!!!!!!!!!!!!!!")
		
	int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

	//KRATOS_WATCH(FractionalStepNumber)

	if(FractionalStepNumber == 1) //first step of the fractional step solution
	{
		const double& density = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
							GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
							GetGeometry()[2].FastGetSolutionStepValue(DENSITY) +
							GetGeometry()[3].FastGetSolutionStepValue(DENSITY));		
		
		if(rRightHandSideVector.size() != 12)
			rRightHandSideVector.resize(12,false);
		unsigned int number_of_nodes = GetGeometry().size();

		//calculate current area
		double current_vol;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, current_vol);

		//writing the body force
		const array_1d<double,3>& body_force = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+
							GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +
							GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE) +
							GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE));

		//const array_1d<double,3>& body_force = GetProperties()[BODY_FORCE];
		for(unsigned int i = 0; i<number_of_nodes; i++)
		{
			rRightHandSideVector[i*3] = body_force[0]* density * mA0 * 0.25;
			rRightHandSideVector[i*3+1] = body_force[1] * density * mA0 * 0.25;
			rRightHandSideVector[i*3+2] = body_force[2] * density * mA0 * 0.25;		
		}
		//get the nodal pressure
		double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
		double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
		double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
		double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);

		//adding pressure gradient 
		double pavg = p0 + p1 + p2 + p3; //calculate old pressure over the element
		pavg *= 0.25 * current_vol;

		rRightHandSideVector[0] += msDN_DX(0,0)*pavg;
		rRightHandSideVector[1] += msDN_DX(0,1)*pavg;
		rRightHandSideVector[2] += msDN_DX(0,2)*pavg;
		
		rRightHandSideVector[3] += msDN_DX(1,0)*pavg;
		rRightHandSideVector[4] += msDN_DX(1,1)*pavg;
		rRightHandSideVector[5] += msDN_DX(1,2)*pavg;

		rRightHandSideVector[6] += msDN_DX(2,0)*pavg;
		rRightHandSideVector[7] += msDN_DX(2,1)*pavg;
		rRightHandSideVector[8] += msDN_DX(2,2)*pavg;

		rRightHandSideVector[9]	 += msDN_DX(3,0)*pavg;
		rRightHandSideVector[10] += msDN_DX(3,1)*pavg;
		rRightHandSideVector[11] += msDN_DX(3,2)*pavg;
	}
	KRATOS_CATCH("")
  	}

	//************************************************************************************
	//************************************************************************************
	void UlfFrac3D::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
				
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
	void UlfFrac3D::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();
		
		if(rDampMatrix.size1() != 12)
			rDampMatrix.resize(12,12,false);

		
		//getting data for the given geometry
		double current_volume;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, current_volume);

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

			msB(0,start) =	msDN_DX(i,0); 
			msB(1,start+1)=	msDN_DX(i,1);
			msB(2,start+2)= msDN_DX(i,2);
			msB(3,start) =	msDN_DX(i,1);		msB(3,start+1) = msDN_DX(i,0);
			msB(4,start) =	msDN_DX(i,2);		msB(4,start+2) = msDN_DX(i,0);
			msB(5,start+1)= msDN_DX(i,2);		msB(5,start+2) = msDN_DX(i,1);
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
	void UlfFrac3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		//save original Area
		mA0 = GeometryUtils::CalculateVolume3D(GetGeometry());		
		
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void UlfFrac3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		

		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 3;		

		unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if (FractionalStepNumber==1)
			{
			if(rResult.size() != number_of_nodes*dim)
				rResult.resize(number_of_nodes*dim,false);	
			for (unsigned int i=0;i<number_of_nodes;i++)
				{
				rResult[i*dim] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
				rResult[i*dim+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
				rResult[i*dim+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
				}
			}
		else if (FractionalStepNumber==2)
			{
			if(rResult.size() != number_of_nodes)
				rResult.resize(number_of_nodes,false);	

						
				for (unsigned int i=0;i<number_of_nodes;i++)
					rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
			}

	}

	//************************************************************************************
	//************************************************************************************
	void UlfFrac3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;

		
		
		unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if (FractionalStepNumber==1)
			{
			if(ElementalDofList.size() != number_of_nodes*dim)
				ElementalDofList.resize(number_of_nodes*dim);	
			for (unsigned int i=0;i<number_of_nodes;i++)
				{
				ElementalDofList[i*dim] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
				ElementalDofList[i*dim+1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
				ElementalDofList[i*dim+2] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
				}
			}

		else if(FractionalStepNumber == 2) // pressure correction step
			{
			if(ElementalDofList.size() != number_of_nodes)
				ElementalDofList.resize(number_of_nodes);	

								
			for (unsigned int i=0;i<number_of_nodes;i++)
				ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);

			}
	}

	//************************************************************************************
	//************************************************************************************
	void UlfFrac3D::GetValuesVector(Vector& values, int Step)
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
	  void UlfFrac3D::GetFirstDerivativesVector(Vector& values, int Step)
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
	void UlfFrac3D::GetSecondDerivativesVector(Vector& values, int Step)
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
	//this function is executed for assembly of the local system corresponding to the modified step1 of FRACSTEP 
	void UlfFrac3D::VelocityStep(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		//KRATOS_WATCH("Velocity step!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		const double& density = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
							GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
							GetGeometry()[2].FastGetSolutionStepValue(DENSITY)+
							GetGeometry()[3].FastGetSolutionStepValue(DENSITY));
		
		//unsigned int dim = 2;
		unsigned int number_of_nodes = 4;

		if(rLeftHandSideMatrix.size1() != 12)
			rLeftHandSideMatrix.resize(12,12,false);

		if(rRightHandSideVector.size() != 12)
			rRightHandSideVector.resize(12,false);

		//calculate current area
		double current_vol;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, current_vol);


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

		noalias(rLeftHandSideMatrix) = ZeroMatrix(12,12);
		//get the nodal pressure
		double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
		double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
		double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
		double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);

		//adding pressure gradient 
		double pavg = p0 + p1 + p2 + p3; //calculate old pressure over the element
		pavg *= 0.25 * current_vol;
 
		rRightHandSideVector[0] += msDN_DX(0,0)*pavg;
		rRightHandSideVector[1] += msDN_DX(0,1)*pavg;
		rRightHandSideVector[2] += msDN_DX(0,2)*pavg;
		
		rRightHandSideVector[3] += msDN_DX(1,0)*pavg;
		rRightHandSideVector[4] += msDN_DX(1,1)*pavg;
		rRightHandSideVector[5] += msDN_DX(1,2)*pavg;

		rRightHandSideVector[6] += msDN_DX(2,0)*pavg;
		rRightHandSideVector[7] += msDN_DX(2,1)*pavg;
		rRightHandSideVector[8] += msDN_DX(2,2)*pavg;

		rRightHandSideVector[9]	 += msDN_DX(3,0)*pavg;
		rRightHandSideVector[10] += msDN_DX(3,1)*pavg;
		rRightHandSideVector[11] += msDN_DX(3,2)*pavg;


		
		KRATOS_CATCH("")
	}
	////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////

	void UlfFrac3D::PressureStep(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		
		if(rRightHandSideVector.size() != 4)
		{
			rLeftHandSideMatrix.resize(4,4,false);
			rRightHandSideVector.resize(4,false);
		}

		double dt = rCurrentProcessInfo[DELTA_TIME];		
		

		//fract. vel, that is calculated in the first Fractional Step.. but is saved inside the "VELOCITY" VARIABLE
		//so, u_n os VELOCITY, 1 and u_n-1 VELOCITY,2 
		const array_1d<double,3>& d0_tilde = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
		const array_1d<double,3>& d0_n = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT,1);	
		const array_1d<double,3>& acc0_n = GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION,1);		
		const array_1d<double,3>& acc0 = GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION);	
		//const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
		//const array_1d<double,3>& fv0_old = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
		double p0_old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const array_1d<double,3>& ff0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);	
		
		const array_1d<double,3>& d1_tilde = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT);
		const array_1d<double,3>& d1_n = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT,1);	
		const array_1d<double,3>& acc1_n = GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION,1);
		const array_1d<double,3>& acc1 = GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION);		
		//const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
		//const array_1d<double,3>& fv1_old = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
		double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE); 
		double p1_old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT); 	 	
		const array_1d<double,3>& ff1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);	
		
		
		const array_1d<double,3>& d2_tilde = GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT);
		const array_1d<double,3>& d2_n = GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT,1);	
		const array_1d<double,3>& acc2_n = GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION,1);	
		const array_1d<double,3>& acc2 = GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION);	
		//const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);	
		//const array_1d<double,3>& fv2_old = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
		double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE); 
		double p2_old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT); 	 	
		const array_1d<double,3>& ff2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);

		const array_1d<double,3>& d3_tilde = GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT);
		const array_1d<double,3>& d3_n = GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT,1);	
		const array_1d<double,3>& acc3_n = GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION,1);	
		const array_1d<double,3>& acc3 = GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION);	
		//const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);	
		//const array_1d<double,3>& fv2_old = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
		const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
		const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);
		double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE); 
		double p3_old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE_OLD_IT); 	 	
		const array_1d<double,3>& ff3 = GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE);
		

		
			

		//getting data for the given geometry
		double Volume;
		GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Volume);

		//calculating average density and viscosity
		double nu = 0.25*(nu0 + nu1 + nu2 + nu3 );
 		double density = 0.25*(rho0 + rho1 + rho2 + rho3 );

		//calculating parameter tau (saved internally to each element)
		double h = pow(12.00*Volume,0.3333333);
		h*=(2.0/sqrt(3.0));

		double tau = 1.00 / ( 4.00*nu/(h*h) + 1.0/dt);
	
						
		//AND NOW WE ADD THE RESPECTIVE CONTRIBUTIONS TO THE RHS AND LHS of THE SECOND FRAC STEP
		//we use Backward Euler for this step, therefore stab. contribution no RHS +=Tau1*(gradQ, residual)
		//								   and LHS +=Tau1*(gradQ, gradP)
		//laplacian term	       L = Dt * gradN * trans(gradN);
		//stabilization term       Spp = tau * gradN * trans(gradN);
		//WATCH OUT for DIVISION with RHO - check if it changes or not in case of Momentum being the primary Variable
		//
		//	msWorkMatrix stores the element laplacian
		//
		double alpha_bossak=-0.3;		
		double c1=0.25*(1.0-alpha_bossak);

		noalias(msWorkMatrix)=prod(msDN_DX,trans(msDN_DX));

		noalias(rLeftHandSideMatrix) = (1.0/density)*(c1*dt + tau) * Volume*msWorkMatrix;

		//rhs consists of D*u_tilda (divergence of the Fractional velocity) and the term: Tau1*(nabla_q, residual_gausspoint)
		//fv is u_tilda
		
		//////////////////////////////////////////////////////////
		////////////		AND NOW RHS	//////////////////
		//////////////////////////////////////////////////////////	
	
		//Dirichlet contribution  (that is: LHS*p_new)
		ms_temp_vec_np[0] = p0; 
		ms_temp_vec_np[1] = p1; 
		ms_temp_vec_np[2] = p2; 
		ms_temp_vec_np[3] = p3; 

		//LHS is already multiplied by AREA
		noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,ms_temp_vec_np);
		
		//NOW RHS-=dt L p_old		
		//changing the meaning of temp_vec_np
		ms_temp_vec_np[0] = p0_old; 
		ms_temp_vec_np[1] = p1_old; 
		ms_temp_vec_np[2] = p2_old; 
		ms_temp_vec_np[3] = p3_old; 

		noalias(rRightHandSideVector) += (1.0/density)*Volume*c1*dt* (prod(msWorkMatrix,ms_temp_vec_np)) ;
	
		//***************************************************************************
		
		//here we have the Dd_tila term
		double Gaux;

		Gaux =  msDN_DX(0,0)*(d0_tilde[0]-d0_n[0]) + msDN_DX(0,1)*(d0_tilde[1]-d0_n[1]) + msDN_DX(0,2)*(d0_tilde[2]-d0_n[2]);
		Gaux += msDN_DX(1,0)*(d1_tilde[0]-d1_n[0]) + msDN_DX(1,1)*(d1_tilde[1]-d1_n[1]) + msDN_DX(1,2)*(d1_tilde[2]-d1_n[2]);
		Gaux += msDN_DX(2,0)*(d2_tilde[0]-d2_n[0]) + msDN_DX(2,1)*(d2_tilde[1]-d2_n[1]) + msDN_DX(2,2)*(d2_tilde[2]-d2_n[2]);
		Gaux += msDN_DX(3,0)*(d3_tilde[0]-d3_n[0]) + msDN_DX(3,1)*(d3_tilde[1]-d3_n[1]) + msDN_DX(3,2)*(d3_tilde[2]-d3_n[2]);

		Gaux/=dt;

		rRightHandSideVector[0] -= Volume*Gaux * msN[0]; 
		rRightHandSideVector[1] -= Volume*Gaux * msN[1]; 
		rRightHandSideVector[2] -= Volume*Gaux * msN[2]; 
		rRightHandSideVector[3] -= Volume*Gaux * msN[3];

		Gaux=0.0;
		
		Gaux =  msDN_DX(0,0)*acc0_n[0] + msDN_DX(0,1)*acc0_n[1] + msDN_DX(0,2)*acc0_n[2];
		Gaux += msDN_DX(1,0)*acc1_n[0] + msDN_DX(1,1)*acc1_n[1] + msDN_DX(1,2)*acc1_n[2];
		Gaux += msDN_DX(2,0)*acc2_n[0] + msDN_DX(2,1)*acc2_n[1] + msDN_DX(2,2)*acc2_n[2];
		Gaux += msDN_DX(3,0)*acc3_n[0] + msDN_DX(3,1)*acc3_n[1] + msDN_DX(3,2)*acc3_n[2];

		//double alpha_bossak=-0.3;
		double gamma= 0.5-alpha_bossak;
		double beta=0.25*pow((1.00-alpha_bossak),2.0);
		double newmark_coef=0.5*dt*((gamma/beta)-2.0);

		rRightHandSideVector[0] -= newmark_coef*Volume*Gaux * msN[0]; 
		rRightHandSideVector[1] -= newmark_coef*Volume*Gaux * msN[1]; 
		rRightHandSideVector[2] -= newmark_coef*Volume*Gaux * msN[2];
		rRightHandSideVector[3] -= newmark_coef*Volume*Gaux * msN[3];  
						
		//RHS = +tau*nablaN*f, we reuse aux
		//ms_aux0 stores ff_gauss; 
		
		ms_aux0=0.25*(ff0+ff1+ff2+ff3);
		//ms_aux1 - is the product of: (nabla q, f)
		ms_aux1[0]=msDN_DX(0,0)*ms_aux0[0]+msDN_DX(0,1)*ms_aux0[1] + msDN_DX(0,2)*ms_aux0[2];
		ms_aux1[1]=msDN_DX(1,0)*ms_aux0[0]+msDN_DX(1,1)*ms_aux0[1] + msDN_DX(1,2)*ms_aux0[2];
		ms_aux1[2]=msDN_DX(2,0)*ms_aux0[0]+msDN_DX(2,1)*ms_aux0[1] + msDN_DX(2,2)*ms_aux0[2];
		ms_aux1[3]=msDN_DX(3,0)*ms_aux0[0]+msDN_DX(3,1)*ms_aux0[1] + msDN_DX(3,2)*ms_aux0[2];
				
		rRightHandSideVector += tau*Volume*ms_aux1;
				
		//RHS += -tau*nablaN*Vel (interia term)
		//time dt, because we are working using displacements
		
		ms_vel_gauss[0]=0.25*(acc0[0]+acc1[0]+acc2[0]+acc3[0])*dt;
		ms_vel_gauss[1]=0.25*(acc0[1]+acc1[1]+acc2[1]+acc3[1])*dt;
		ms_vel_gauss[2]=0.25*(acc0[1]+acc1[1]+acc2[1]+acc3[2])*dt;
		
		//and now we reuse ms_aux1

		ms_aux1=prod(msDN_DX,ms_vel_gauss);
		
		noalias(rRightHandSideVector) -= tau*(1.0/c1)*Volume*ms_aux1; //tau*density*Area*ms_aux1;  
				
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// here for the FSI problems involving flexible structures only
		// we add a Laplacian term, taking into account the mass of the structure.. this
		// Laplacian is computed using 3 integration points
		
		//first we identify the interface elements (if such exist)
		
		//int n_interf0=GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
		//int n_interf1=GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
		//int n_interf2=GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
		//this is just to try
		//double dens_str=1000.0;
		//if (n_interf>=3)			
		//	KRATOS_ERROR(std::logic_error,  "Something is wrong: fluid element cannot have all 3 nodes at the FSI boundary " , "");	

		//3 is the total number of nodes
		//so we add the inverse of the fluid density for the non-interface nodes (fluid only)
		//and the inverse of the sum of the fluid and structure density of the interface nodes
		//1/3 - stands as we are integrating the Laplcian using 3 integration points located at the nodes
		
		//rLeftHandSideMatrix*=0.0;
		//adding: L+LT
		/*
		if (n_interf0==1)
			{
			KRATOS_WATCH("First IS INTERF")
			rLeftHandSideMatrix(0,0)+=0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(0,0);

			rLeftHandSideMatrix(1,0)+=0.5*0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(1,0);
			rLeftHandSideMatrix(0,1)+=0.5*0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(1,0);

			rLeftHandSideMatrix(2,0)+=0.5*0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(2,0);
			rLeftHandSideMatrix(0,2)+=0.5*0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(2,0);
			}
		if (n_interf1==1)
			{
			KRATOS_WATCH("Second IS INTERF")
			rLeftHandSideMatrix(0,1)+=0.5*0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(0,1);
			rLeftHandSideMatrix(1,0)+=0.5*0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(0,1);

			rLeftHandSideMatrix(1,1)+=0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(1,1);

			rLeftHandSideMatrix(2,1)+=0.5*0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(2,1);
			rLeftHandSideMatrix(1,2)+=0.5*0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(2,1);
			}

		if (n_interf2==1)
			{
			KRATOS_WATCH("Thirs IS INTERF")
			rLeftHandSideMatrix(0,2)+=0.5*0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(0,2);
			rLeftHandSideMatrix(2,0)+=0.5*0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(0,2);

			rLeftHandSideMatrix(1,2)+=0.5*0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(1,2);
			rLeftHandSideMatrix(2,1)+=0.5*0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(1,2);

			rLeftHandSideMatrix(2,2)+=0.3333333333333333*(1.0/(dens_str))*(c1*dt) * Area*msWorkMatrix(2,2);
			}

		*/	
			
		//noalias(rLeftHandSideMatrix)+=0.3333333333333333*(n_interf/(density+dens_str))*(c1*dt) * Area*msWorkMatrix;
		//0.3333333333333333*((3.0-n_interf)*1.0/density+n_interf/(density+dens_str))*(c1*dt + tau) * Area*msWorkMatrix;
		
		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	
	  void  UlfFrac3D::Calculate(const Variable<double >& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
	  {		
		// this function assembles the pressure equation exactly in the same way as in the classical second stage of fractional techniques
		if(rVariable == PRESSURE)
		{
		KRATOS_WATCH("EMPTY FUNCTION FOR THIS ELEMENT - PRESSUREUES ARE UPDATED INSIDE BUILDER AND SOLVER");
		}
		else if (rVariable == IS_FLUID)
		{
			
		GetGeometry()[0].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;
		GetGeometry()[1].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;
		GetGeometry()[2].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;
		GetGeometry()[3].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;

		}
		else if(rVariable == NODAL_MASS)
		{
		CalculateLumpedMass();	
		}
		
		else
			KRATOS_ERROR(std::logic_error,  "You are doing something wrong  FCT calculate... of nodal_mass with wring parameters.. " , "");

	}


} // Namespace Kratos


