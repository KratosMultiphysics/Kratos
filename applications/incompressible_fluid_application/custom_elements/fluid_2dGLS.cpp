/*
==============================================================================
KratosR1IncompressibleFluidApplication 
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
//   Last modified by:    $Author: jmarti $
//   Date:                $Date: 2008-11-26 08:17:41 $
//   Revision:            $Revision: 1.0 $
//
//
 
//#define GRADPN_FORM

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/fluid_2dGLS.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{

//THIS IS A COMPRESSIBLE FLUID ELEMENT, WITH GLS STABILIZATION, RUNGE-KUTTA Momentum Time integration, FRACTIONAL STEP
	//************************************************************************************
	//************************************************************************************
	Fluid2DGLS::Fluid2DGLS(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Fluid2DGLS::Fluid2DGLS(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{


	}
	void Fluid2DGLS::CalculateLumpedMass()		
	{
	//note that for the compressible case, rho will be also a variable
	
	const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
	const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
	const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);


	
	//double Area;
	//GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);
	double Area = GeometryUtils::CalculateVolume2D(GetGeometry());
	double lumped_mass_fac = Area * 0.33333333333333333;

	double & m0 = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS);
	double & m1 = GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS);
	double & m2 = GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS);

	#pragma omp atomic
	m0+=lumped_mass_fac*rho0;
	#pragma omp atomic
	m1+=lumped_mass_fac*rho1;	
	#pragma omp atomic	
	m2+=lumped_mass_fac*rho2;


	//filling in the diagonal of the lumped mass matrix,  (later I can change it to vector...)
/*	GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho0;
	GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho1;	
	GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho2;

*/

	}

	Element::Pointer Fluid2DGLS::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		KRATOS_TRY

		return Element::Pointer(new Fluid2DGLS(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	Fluid2DGLS::~Fluid2DGLS()
	{
	}

	//************************************************************************************
	//************************************************************************************
	

	//************************************************************************************
	//************************************************************************************
	void Fluid2DGLS::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}
	
	//************************************************************************************
	//************************************************************************************
	
	void Fluid2DGLS::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		//KRATOS_WATCH("Empty function for this element")		
		//KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}

void Fluid2DGLS::CalculateGalerkinMomentumResidual(VectorType& GalerkinRHS)
		{
		KRATOS_TRY
		///////////////////////NECESSARY LOCALS///////////////////////////////////////////
		boost::numeric::ublas::bounded_matrix<double,3,3> msWorkMatrix = ZeroMatrix(3,3);
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
		array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
		boost::numeric::ublas::bounded_matrix<double,2,6> msConvOp = ZeroMatrix(2,6);
	       	boost::numeric::ublas::bounded_matrix<double,6,6> msAuxMat = ZeroMatrix(6,6);
		array_1d<double,6> msAuxVec = ZeroVector(6); //dimension = number of nodes
		array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
		array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension
		///////////////////////////////////////////////////////////////////////////////////

	
		//first we compute  the force term and pressure gradient terms:
		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

		//getting the velocity on the nodes and other necessary variabless
		const array_1d<double,3> vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
		double p_n0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		
		const array_1d<double,3> vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
		double p_n1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
		
		const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
		double p_n2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
			
		//====================================================================
		//calculating viscosity and density
		double nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
		double density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );
		
		//VISCOUS CONTRIBUTION 
		// += Laplacian * nu; --> ONE GAUSS POINT
		//msWorkMatrix is used now to store the element laplacian 3x3
		noalias(msWorkMatrix) = Area*density*nu * prod(msDN_DX,trans(msDN_DX));
				
		//x comp
		GalerkinRHS[0]=-1.0*(msWorkMatrix(0,0)*vel0[0]+msWorkMatrix(0,1)*vel1[0]+msWorkMatrix(0,2)*vel2[0]);
		//y comp
		GalerkinRHS[1]=-1.0*(msWorkMatrix(0,0)*vel0[1]+msWorkMatrix(0,1)*vel1[1]+msWorkMatrix(0,2)*vel2[1]);

		//x comp
		GalerkinRHS[2]=-1.0*(msWorkMatrix(1,0)*vel0[0]+msWorkMatrix(1,1)*vel1[0]+msWorkMatrix(1,2)*vel2[0]);
		//y comp
		GalerkinRHS[3]=-1.0*(msWorkMatrix(1,0)*vel0[1]+msWorkMatrix(1,1)*vel1[1]+msWorkMatrix(1,2)*vel2[1]);

		//x comp
		GalerkinRHS[4]=-1.0*(msWorkMatrix(2,0)*vel0[0]+msWorkMatrix(2,1)*vel1[0]+msWorkMatrix(2,2)*vel2[0]);
		//y comp
		GalerkinRHS[5]=-1.0*(msWorkMatrix(2,0)*vel0[1]+msWorkMatrix(2,1)*vel1[1]+msWorkMatrix(2,2)*vel2[1]);

		//convective contribution. Note that N[0]=N[1]=N[2]=0.33333 for our 1-Gauss Point integration
		//WATCH OUT THAT I AM NOT PLUGGING THE NEW ADVECTIVE VELOCITY WHEN EXECUTING RUNGE KUTTA - I use u_n (and not the intermediate u_aux)
		//KRATOS_WATCH("Now lets see N - they should be all equal 0.3333")
		//KRATOS_WATCH(msN)
		//KRATOS_WATCH(msDN_DX)
		ms_adv_vel[0] = msN[0]*(vel0[0])+msN[1]*(vel1[0])+msN[2]*(vel2[0]);
		ms_adv_vel[1] = msN[0]*(vel0[1])+msN[1]*(vel1[1])+msN[2]*(vel2[1]);

		//calculate convective term	
		int nodes_number = 3;
		int dof = 2;
		
		
		//and now we add the pressure gradient and the force term
		//external forces (component)
		const array_1d<double,3> body_force = 0.333333333333333*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+
							GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +
							GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE));
		unsigned int number_of_nodes=3;
		for(unsigned int i = 0; i<number_of_nodes; i++)
		{
			//f=A*N_I*b, N_I=0.33333333 for 1 Gauss point
			GalerkinRHS[i*2] += body_force[0]* density * Area * 0.3333333333333;
			GalerkinRHS[i*2+1] += body_force[1] * density * Area * 0.3333333333333;
		}
		
		
		
		//Now we shall add the Gp term

		double p_avg=0.333333333333*(p_n0+p_n1+p_n2)*Area;
		GalerkinRHS[0]+=msDN_DX(0,0)*p_avg;
		GalerkinRHS[1]+=msDN_DX(0,1)*p_avg; 

		GalerkinRHS[2]+=msDN_DX(1,0)*p_avg; 
		GalerkinRHS[3]+=msDN_DX(1,1)*p_avg; 

		GalerkinRHS[4]+=msDN_DX(2,0)*p_avg; 
		GalerkinRHS[5]+=msDN_DX(2,1)*p_avg; 
		
		
		/*GetGeometry()[0].FastGetSolutionStepValue(FORCE_X) += GalerkinRHS[0];
		GetGeometry()[0].FastGetSolutionStepValue(FORCE_Y) += GalerkinRHS[1];
		GetGeometry()[1].FastGetSolutionStepValue(FORCE_X) += GalerkinRHS[2];
		GetGeometry()[1].FastGetSolutionStepValue(FORCE_Y) += GalerkinRHS[3];
		GetGeometry()[2].FastGetSolutionStepValue(FORCE_X) += GalerkinRHS[4];
		GetGeometry()[2].FastGetSolutionStepValue(FORCE_Y) += GalerkinRHS[5];*/

		GetGeometry()[0].SetLock();
		array_1d<double,3>& rhs0 = GetGeometry()[0].FastGetSolutionStepValue(FORCE);
		rhs0[0] += GalerkinRHS[0];
		rhs0[1] += GalerkinRHS[1];
		GetGeometry()[0].UnSetLock();

		GetGeometry()[1].SetLock();
		array_1d<double,3>& rhs1 = GetGeometry()[1].FastGetSolutionStepValue(FORCE);
		rhs1[0] += GalerkinRHS[2];
		rhs1[1] += GalerkinRHS[3];
		GetGeometry()[1].UnSetLock();
	
		GetGeometry()[2].SetLock();
		array_1d<double,3>& rhs2 = GetGeometry()[2].FastGetSolutionStepValue(FORCE);
		rhs2[0] += GalerkinRHS[4];
		rhs2[1] += GalerkinRHS[5];
		GetGeometry()[2].UnSetLock();



				
		KRATOS_CATCH("")

	}

	/*void Fluid2DGLS::CalculateRHSVector(VectorType& Galerkin_RHS, double& dt)
	{
		KRATOS_TRY
		
		
		KRATOS_CATCH("")

	}*/

	
	void Fluid2DGLS::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		
		///////////////////////NECESSARY LOCALS///////////////////////////////////////////
		boost::numeric::ublas::bounded_matrix<double,3,3> msWorkMatrix = ZeroMatrix(3,3);
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
		array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
		boost::numeric::ublas::bounded_matrix<double,2,6> msConvOp = ZeroMatrix(2,6);
	       	boost::numeric::ublas::bounded_matrix<double,6,6> msAuxMat = ZeroMatrix(6,6);
		array_1d<double,6> msAuxVec = ZeroVector(6); //dimension = number of nodes
		array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
		array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension
		array_1d<double,3> ms_temp_vec_np = ZeroVector(3); //dimension = number of nodes
		array_1d<double,3> ms_aux0 = ZeroVector(3); //dimension = number of nodes
		array_1d<double,3> ms_aux1 = ZeroVector(3); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,6,6> msAuxMat1 = ZeroMatrix(6,6);
    		boost::numeric::ublas::bounded_matrix<double,6,3> msAuxMat2 = ZeroMatrix(6,3);
		boost::numeric::ublas::bounded_matrix<double,2,2> msGrad_ug = ZeroMatrix(2,2);
		array_1d<double,6> msStabMomRes = ZeroVector(6); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,6,3> msGradOp = ZeroMatrix(6,3);
		///////////////////////////////////////////////////////////////////////////////////
		
		if(rRightHandSideVector.size() != 3)
		{
			rLeftHandSideMatrix.resize(3,3,false);
			rRightHandSideVector.resize(3,false);
		}

		double dt = rCurrentProcessInfo[DELTA_TIME];
		
		//fract. vel, that is calculated in the first Fractional Step.. but is saved inside the "VELOCITY" VARIABLE
		//so, u_n os VELOCITY, 1 and u_n-1 VELOCITY,2 
		const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& fv0_old = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
		const array_1d<double,3>& fv0_n_1 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,2);
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
		//double p0_old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		double p0_old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1); 	 	 	 	
		const array_1d<double,3>& ff0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);	

		const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& fv1_old = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
		const array_1d<double,3>& fv1_n_1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,2);
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
		double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE); 
		//double p1_old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT); 	 	
		double p1_old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1); 	 	 	 	
		const array_1d<double,3>& ff1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);	

		const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);	
		const array_1d<double,3>& fv2_old = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
		const array_1d<double,3>& fv2_n_1 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,2);
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
		double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE); 
		double p2_old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1); 	 	 	 	
		//old iteration can be used if we want to iterate between 1st and 2nd fractional steps
		//double p2_old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT); 	 	
		const array_1d<double,3>& ff2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);	


		//in msAuxVec we store the velocity, (not the frac step vel, but u_n, the one that enters the stabilization)
		msAuxVec[0]=fv0[0];
		msAuxVec[1]=fv0[1];
		msAuxVec[2]=fv1[0];
		msAuxVec[3]=fv1[1];
		msAuxVec[4]=fv2[0];
		msAuxVec[5]=fv2[1];

		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

		//calculating average density and viscosity
		double nu = 0.33333333333333*(nu0 + nu1 + nu2 );
 		double density = 0.33333333333333*(rho0 + rho1 + rho2 );

		//ms_vel_gauss[i] =  msN[0]*(fv0[i]) + msN[1]*(fv1[i]) +  msN[2]*(fv2[i]);
		//but with one integration N=0.333333333
		ms_vel_gauss[0] =  0.33333333333333*(fv0[0]+fv1[0]+fv2[0]);
		ms_vel_gauss[1] =  0.33333333333333*(fv0[1]+fv1[1]+fv2[1]);

		//calculating parameter tau (saved internally to each element)
		double h = sqrt(2.00*Area);
		double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
		norm_u = sqrt(norm_u);
		double tau = 1.00 / ( 4.00*nu/(h*h) /*+2.00*norm_u/h */+ 1.0/dt);
		
		//AND NOW WE ADD THE RESPECTIVE CONTRIBUTIONS TO THE RHS AND LHS of THE SECOND FRAC STEP
		//we use Backward Euler for this step, therefore stab. contribution no RHS +=Tau1*(gradQ, residual)
		//								   and LHS +=Tau1*(gradQ, gradP)
		//laplacian term	       L = Dt * gradN * trans(gradN);
		//stabilization term       Spp = tau * gradN * trans(gradN);
		//WATCH OUT for DIVISION with RHO - check if it changes or not in case of Momentum being the primary Variable
		//
		//	msWorkMatrix stores the element laplacian
		//
		noalias(msWorkMatrix)=prod(msDN_DX,trans(msDN_DX));
		noalias(rLeftHandSideMatrix) = (dt/2.0 + tau) * Area*msWorkMatrix;
		//rhs consists of D*u_tilda (divergence of the Fractional velocity) and the term: Tau1*(nabla_q, residual_gausspoint)
		//fv is u_tilda
		
		//////////////////////////////////////////////////////////
		////////////		AND NOW RHS	//////////////////
		//////////////////////////////////////////////////////////	
	
		//Dirichlet contribution  (that is: LHS*p_new)
		ms_temp_vec_np[0] = p0; 
		ms_temp_vec_np[1] = p1; 
		ms_temp_vec_np[2] = p2; 
		//LHS is already multiplied by AREA
		noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,ms_temp_vec_np);
		
		//NOW RHS-=dt L p_old		
		//changing the meaning of temp_vec_np
		ms_temp_vec_np[0] = p0_old; 
		ms_temp_vec_np[1] = p1_old; 
		ms_temp_vec_np[2] = p2_old; 

		noalias(rRightHandSideVector) += Area* dt/2.0 * (prod(msWorkMatrix,ms_temp_vec_np)) ;
		
		//***************************************************************************
		
		//here we have the Du_tila term
		double Gaux;
		Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1];
		Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1];
		Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1];
		rRightHandSideVector[0] -= density*Area*Gaux * msN[0]; 
		rRightHandSideVector[1] -= density*Area*Gaux * msN[1]; 
		rRightHandSideVector[2] -= density*Area*Gaux * msN[2]; 
		
		
		
		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	

	void Fluid2DGLS::FinalFractionalStep(const ProcessInfo& rCurrentProcessInfo)		
	{
			KRATOS_TRY
			KRATOS_ERROR(std::logic_error,  "METHOD NOT IMPL inside the element Final Fractional Step is done within the low_mach strategy.. " , "");
			
			KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	
			
	



	//************************************************************************************
	//************************************************************************************
	void Fluid2DGLS::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
	  {
		//if the VAr is NODAL_MASS, we calculate the lumped mass		
		if(rVariable == NODAL_MASS)
		{
		CalculateLumpedMass();	
		}	
		else
			KRATOS_ERROR(std::logic_error,  "You are doing something wrong  FCT calculate... of nodal_mass with wring parameters.. " , "");
	  }
	//************************************************************************************
	//************************************************************************************
	void Fluid2DGLS::Calculate(const Variable<array_1d<double,3> >& rVariable, array_1d<double,3>& Output, const ProcessInfo& rCurrentProcessInfo)
	  {
			
		if(rVariable == VELOCITY && Output[0]==1.0)
		//we use "Output" as a switch between 1st Frac Step and last Frac Step(Correction Step)
		{
		//here the residual will be temporarily written
		Vector TmpRhs(6);
		// first we write the Galerkin contributions to the momentum residual						
		CalculateGalerkinMomentumResidual(TmpRhs);
		//and now the stabilization terms added
		double dt = rCurrentProcessInfo[DELTA_TIME];
		//CalculateRHSVector(TmpRhs,  dt);

		}
		else if(rVariable == VELOCITY && Output[0]==2.0)
		{
		FinalFractionalStep(rCurrentProcessInfo);
		}
		else 
		{
			KRATOS_ERROR(std::logic_error,  "You are doing something wrong in ur fractional step.... " , "");
		}


	}
	//************************************************************************************
	//************************************************************************************
	//************************************************************************************
	//************************************************************************************
	void Fluid2DGLS::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		
		
			if(rResult.size() != number_of_nodes)
				rResult.resize(number_of_nodes,false);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
			}
		KRATOS_CATCH("")
		
	}

	//************************************************************************************
	//************************************************************************************
	  void Fluid2DGLS::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		
			if(ElementalDofList.size() != number_of_nodes)
				ElementalDofList.resize(number_of_nodes);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);
			}
		KRATOS_CATCH("");
	}



} // Namespace Kratos


