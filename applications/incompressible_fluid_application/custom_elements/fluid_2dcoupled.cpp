/*
==============================================================================
KratosIncompressibleFluidApplication 
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
//   Date:                $Date: 2009-01-15 11:11:35 $
//   Revision:            $Revision: 1.6 $
//
//
 
//#define GRADPN_FORM
//#define STOKES

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/fluid_2dcoupled.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
	//static variables
	boost::numeric::ublas::bounded_matrix<double,3,3> Fluid2DCoupled::msaux_matrix;
	boost::numeric::ublas::bounded_matrix<double,3,3> Fluid2DCoupled::msMassFactors;
	boost::numeric::ublas::bounded_matrix<double,3,2> Fluid2DCoupled::msDN_DX;
  	array_1d<double,3> Fluid2DCoupled::msN; //dimension = number of nodes
	array_1d<double,2> Fluid2DCoupled::ms_vel_gauss; //dimesion coincides with space dimension
  	array_1d<double,3> Fluid2DCoupled::ms_temp_vec_np; //dimension = number of nodes
  	array_1d<double,3> Fluid2DCoupled::ms_u_DN; //dimension = number of nodes


	//************************************************************************************
	//************************************************************************************
	Fluid2DCoupled::Fluid2DCoupled(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Fluid2DCoupled::Fluid2DCoupled(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
		msMassFactors(0,0) = 1.00/3.00; msMassFactors(0,1) = 0.00;		msMassFactors(0,2) = 0.00;
		msMassFactors(1,0) = 0.00;		msMassFactors(1,1) = 1.00/3.00; msMassFactors(1,2) = 0.00;
		msMassFactors(2,0) = 0.00;		msMassFactors(2,1) = 0.00;		msMassFactors(2,2) = 1.00/3.00;		
	}

	Element::Pointer Fluid2DCoupled::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		KRATOS_TRY

		return Element::Pointer(new Fluid2DCoupled(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	Fluid2DCoupled::~Fluid2DCoupled()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void Fluid2DCoupled::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber < 3) //first step of the fractional step solution
		{
			Stage1(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
		}
		else if (FractionalStepNumber == 4)//second step of the fractional step solution
		{
			Stage2(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
		}

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void Fluid2DCoupled::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}

	//************************************************************************************
	//************************************************************************************
	//calculation by component of the fractional step velocity corresponding to the first stage
	void Fluid2DCoupled::Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										   ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		unsigned int number_of_points = 3;
		unsigned int dim = 2;
		unsigned int matsize = number_of_points *dim;
		
		if(rLeftHandSideMatrix.size1() != matsize)
			rLeftHandSideMatrix.resize(matsize,matsize,false);

		if(rRightHandSideVector.size() != matsize)
			rRightHandSideVector.resize(matsize,false);

		//getting data for the given geometry
		double Area;
		//CalculateGeometryData(msDN_DX,msN,Area);
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);

		//getting the velocity vector on the nodes

		//getting the velocity on the nodes
		const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL,0);
		const array_1d<double,3>& w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj0 = GetGeometry()[0].FastGetSolutionStepValue(CONV_PROJ);
		double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

		const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj1 = GetGeometry()[1].FastGetSolutionStepValue(CONV_PROJ);
		double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

		const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj2 = GetGeometry()[2].FastGetSolutionStepValue(CONV_PROJ);
		double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

		// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
		ms_vel_gauss[0] =  (fv0[0]-w0[0]) +(fv1[0]-w1[0]) +  (fv2[0]-w2[0]);
		ms_vel_gauss[1] =  (fv0[1]-w0[1]) +(fv1[1]-w1[1]) +  (fv2[1]-w2[1]);
		ms_vel_gauss[0] *= 0.33333333333333333333;
		ms_vel_gauss[1] *= 0.33333333333333333333;

		//calculating viscosity
		double nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
 		double density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );

		//getting the BDF2 coefficients (not fixed to allow variable time step)
		//the coefficients INCLUDE the time step
		const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

		//calculating parameter tau (saved internally to each element)
		double c1 = 4.00;
		double c2 = 2.00;
		double h = sqrt(2.00*Area);
		double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
		norm_u = sqrt(norm_u);
		double tau = 1.00 / ( c1*nu/(h*h) + c2*norm_u/h );
	
		//VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
		// rLeftHandSideMatrix += Laplacian * nu;
/*		boost::numeric::ublas::bounded_matrix<double,3,6> B;
		boost::numeric::ublas::bounded_matrix<double,3,3> constitutive_matrix;

			//filling matrix B
			for (unsigned int i=0;i<number_of_points;i++)
			{
				unsigned int index = dim*i;

				B(0,index+0)=msDN_DX(i,0);							B(0,index+1)= 0.0;
				B(1,index+0)=0.0;									B(1,index+1)= msDN_DX(i,1);
				B(2,index+0)= msDN_DX(i,1);							B(2,index+1)= msDN_DX(i,0); 
			}			
				
			//constitutive tensor
			constitutive_matrix(0,0) = 2.0*nu;	 constitutive_matrix(0,1) = 0.0;		constitutive_matrix(0,2) = 0.0;
			constitutive_matrix(1,0) = 0.0;		 constitutive_matrix(1,1) = 2.0*nu;		constitutive_matrix(1,2) = 0.0;
			constitutive_matrix(2,0) = 0.0;		 constitutive_matrix(2,1) = 0.0;		constitutive_matrix(2,2) = nu;
		
			
			//calculating viscous contributions
			boost::numeric::ublas::bounded_matrix<double,3,6> temp = prod( constitutive_matrix , B);
			noalias(rLeftHandSideMatrix) = prod( trans(B) , temp);
*/	
		//performs Btrans * D * B without writing any of the matrices
		CalculateViscousMatrix(rLeftHandSideMatrix, msDN_DX, nu);

		//CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
		noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
		noalias(msaux_matrix) = outer_prod(msN,ms_u_DN);

		//CONVECTION STABILIZING CONTRIBUTION (Suu)
		noalias(msaux_matrix) += tau * outer_prod(ms_u_DN,ms_u_DN);

			//INERTIA CONTRIBUTION
		noalias(msaux_matrix) += BDFcoeffs[0] * msMassFactors;

		//adding all contributions to the stiffness matrix
		ExpandAndAddReducedMatrix(rLeftHandSideMatrix,msaux_matrix,dim);

		//multiplication by the area
		rLeftHandSideMatrix *= Area * density;


		// *****************************************
		//CALCULATION OF THE RHS
		const array_1d<double,3>& force0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
		const array_1d<double,3>& force1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
		const array_1d<double,3>& force2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);		
		
		array_1d<double,3> rhs_aux;
		for( unsigned int component_index = 0; component_index < dim; component_index++)
		{
			//external forces (component)
			double force_component = 0.3333333333333333*(force0[component_index] + force1[component_index] + force2[component_index]);
			noalias(rhs_aux) = (force_component )*msN;

			//adding pressure gradient (integrated by parts)
			double p_avg = p0old + p1old + p2old;
			p_avg *= 0.3333333333333333 / density;
			rhs_aux[0] += msDN_DX(0,component_index)*p_avg; 
			rhs_aux[1] += msDN_DX(1,component_index)*p_avg; 
			rhs_aux[2] += msDN_DX(2,component_index)*p_avg;

			//adding the inertia terms
			// RHS += M*vhistory 
			//calculating the historical velocity
			//noalias(ms_temp_vec_np) = ZeroVector(3);		
			for(unsigned int iii = 0; iii<number_of_points; iii++)
			{
				const array_1d<double,3>& v = (GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY,1) );
				ms_temp_vec_np[iii] = -BDFcoeffs[1]*v[component_index];
			}
			for(unsigned int step = 2; step<BDFcoeffs.size(); step++)
			{
				for(unsigned int iii = 0; iii<number_of_points; iii++)
				{
					const array_1d<double,3>& v = (GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY,step) );
					ms_temp_vec_np[iii] -= BDFcoeffs[step]*v[component_index];
				}
			}	
			noalias(rhs_aux) += prod(msMassFactors,ms_temp_vec_np) ;

			//RHS += Suy * proj[component] 
			double proj_component = proj0[component_index] 
								  + proj1[component_index] 
								  + proj2[component_index];
			proj_component *= 0.3333333333333333;
			noalias(rhs_aux) += (tau*proj_component)*ms_u_DN;   

			//writing the rhs_aux in its place
			for( unsigned int i = 0; i < number_of_points; i++)
			{
				rRightHandSideVector[i*dim + component_index] = rhs_aux[i];
			}
		}
		
		//multiplying by area
		rRightHandSideVector *= Area * density;

		//LHS stabilization contribution
		Vector fvvect(6);
		for( unsigned int component_index = 0; component_index < dim; component_index++)
		{
				fvvect[0 + component_index] = fv0[component_index];
				fvvect[2 + component_index] = fv1[component_index];
				fvvect[4 + component_index] = fv2[component_index];
		}
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,fvvect);
		

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	//calculation by component of the fractional step velocity corresponding to the first stage
	void Fluid2DCoupled::Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										   ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		unsigned int number_of_points = 3;

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);

		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);
		
		const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
		double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE); 
		double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT); 
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		
		const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
		double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE); 
		double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT); 
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
		
		const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);
		double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE); 
		double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT); 
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

		// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
		ms_vel_gauss[0] =  msN[0]*(fv0[0]-w0[0]) + msN[1]*(fv1[0]-w1[0]) +  msN[2]*(fv2[0]-w2[0]);
		ms_vel_gauss[1] =  msN[0]*(fv0[1]-w0[1]) + msN[1]*(fv1[1]-w1[1]) +  msN[2]*(fv2[1]-w2[1]);

		//calculating convective auxiliary vector

		noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
		
		//calculating average density and viscosity
		double nu = 0.33333333333333*(nu0 + nu1 + nu2 );
 		double density = 0.33333333333333*(rho0 + rho1 + rho2 );

		//calculating parameter tau (saved internally to each element)
		double h = sqrt(2.00*Area);
		double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
		norm_u = sqrt(norm_u);
		double tau = 1.00 / ( 4.00*nu/(h*h) + 2.00*norm_u/h );

		//getting the BDF2 coefficients (not fixed to allow variable time step)
		//the coefficients INCLUDE the time step
		const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

		//CALCULATION OF THE LEFT HAND SIDE
		//laplacian term	       L = Dt * gradN * trans(gradN);
		//stabilization term       Spp = tau * gradN * trans(gradN);
		noalias(rLeftHandSideMatrix) = ((1.00/BDFcoeffs[0] + tau)/density) * prod(msDN_DX,trans(msDN_DX));

		//calculation of the RHS
		// RHS = -G*vfrac
		double Gaux;
		Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1];
		Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1];
		Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1];
		rRightHandSideVector[0] = - Gaux * msN[0]; 
		rRightHandSideVector[1] = - Gaux * msN[1]; 
		rRightHandSideVector[2] = - Gaux * msN[2]; 

		//attention!! changing the meaning of ms_vel_gauss
		// RHS += Sz * proj 
		//contrib of Spy*proj
		ms_vel_gauss[0] = msN[0]*proj0[0] + msN[1]*proj1[0] + msN[2]*proj2[0];
		ms_vel_gauss[1] = msN[0]*proj0[1] + msN[1]*proj1[1] + msN[2]*proj2[1];
		ms_vel_gauss *= tau;
		noalias(rRightHandSideVector) += prod(msDN_DX , ms_vel_gauss);   

		//dirichlet contribution
		ms_temp_vec_np[0] = p0; 
		ms_temp_vec_np[1] = p1; 
		ms_temp_vec_np[2] = p2; 
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);

		// RHS += dt * L * pold 
		ms_temp_vec_np[0] = p0old; 
		ms_temp_vec_np[1] = p1old; 
		ms_temp_vec_np[2] = p2old; 
		noalias(ms_vel_gauss) = prod(trans(msDN_DX),ms_temp_vec_np);
		noalias(rRightHandSideVector) += (1.00/BDFcoeffs[0]/density) * prod(msDN_DX,ms_vel_gauss); 

		//multiplicating by the area
		rLeftHandSideMatrix *= Area;
		rRightHandSideVector *= Area;

		//adding contributions to nodal areas following the corresponding lumping term
		double nodal_contrib = 0.333333333333333333333333333 * Area*density;
		GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
		GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
		GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;

		KRATOS_CATCH("");
	}
	  
	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void Fluid2DCoupled::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
 
		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);		
		
		if(FractionalStepNumber  == 5) //calculation of stabilization terms
		{

			array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
			array_1d<double,3>& w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
			array_1d<double,3>& press_proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
			array_1d<double,3>& conv_proj0 = GetGeometry()[0].FastGetSolutionStepValue(CONV_PROJ);
			double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
			const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
			
			array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
			array_1d<double,3>& w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
			array_1d<double,3>& press_proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
			array_1d<double,3>& conv_proj1 = GetGeometry()[1].FastGetSolutionStepValue(CONV_PROJ);
			double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
			const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

			array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
			array_1d<double,3>& w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
			array_1d<double,3>& press_proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);
			array_1d<double,3>& conv_proj2 = GetGeometry()[2].FastGetSolutionStepValue(CONV_PROJ);
			double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
			const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

			double density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );

			//calculation of the pressure gradient (saved in ms_vel_gauss)
			//note that here we calculate it "strong"
			ms_vel_gauss[0] = msDN_DX(0,0)*(p0) + msDN_DX(1,0)*(p1) + msDN_DX(2,0)*(p2);
			ms_vel_gauss[1] = msDN_DX(0,1)*(p0) + msDN_DX(1,1)*(p1) + msDN_DX(2,1)*(p2);
			ms_vel_gauss *= Area;

			//press_proj += G*p
			press_proj0[0] += msN[0]*ms_vel_gauss[0]; 
			press_proj0[1] += msN[0]*ms_vel_gauss[1]; 

			press_proj1[0] += msN[1]*ms_vel_gauss[0]; 
			press_proj1[1] += msN[1]*ms_vel_gauss[1]; 

			press_proj2[0] += msN[2]*ms_vel_gauss[0]; 
			press_proj2[1] += msN[2]*ms_vel_gauss[1]; 

			// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
			ms_vel_gauss[0] =  msN[0]*(fv0[0]-w0[0]) + msN[1]*(fv1[0]-w1[0]) +  msN[2]*(fv2[0]-w2[0]);
			ms_vel_gauss[1] =  msN[0]*(fv0[1]-w0[1]) + msN[1]*(fv1[1]-w1[1]) +  msN[2]*(fv2[1]-w2[1]);

			//calculating convective auxiliary vector
			noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);

			//attention changing the meaning of ms_vel_gauss!!
			ms_vel_gauss[0] = ms_u_DN[0] * fv0[0] + ms_u_DN[1] * fv1[0] + ms_u_DN[2] * fv2[0];
			ms_vel_gauss[1] = ms_u_DN[0] * fv0[1] + ms_u_DN[1] * fv1[1] + ms_u_DN[2] * fv2[1];
			ms_vel_gauss *= Area * density;

			// conv_proj += C*u
			conv_proj0[0] += msN[0]*ms_vel_gauss[0]; 
			conv_proj0[1] += msN[0]*ms_vel_gauss[1];

			conv_proj1[0] += msN[1]*ms_vel_gauss[0]; 
			conv_proj1[1] += msN[1]*ms_vel_gauss[1];

			conv_proj2[0] += msN[2]*ms_vel_gauss[0];
			conv_proj2[1] += msN[2]*ms_vel_gauss[1];
		}		
		else if(FractionalStepNumber == 6) //calculation of velocities
		{
			//outside of the element it is performed a loop on the elements in which it is considered nodally
			//v = vfrac - Dt/Mnodal * G*(p-pold)
			// the term G*(p-pold) needs to be assembled by element. the contribution is thus directly added to the nodal contrib

			double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
			double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
			array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
			//const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

			double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
			double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
			array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
			//const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

			double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
			double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
			array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
			//const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

			//double density = 0.33333333333333333333333*(rho0 + rho1 + rho2);
			

#if defined(GRADPN_FORM)
			//adding pressure gradient integrated by parts
			// fv += G*p(n+1) - grad p(n)
			double p_avg =  msN[0]*p0 + msN[1]*p1 + msN[2]*p2;
			p_avg *= Area/density;
			fv0[0] += msDN_DX(0,0)*p_avg; fv0[1] += msDN_DX(0,1)*p_avg;
			fv1[0] += msDN_DX(1,0)*p_avg; fv1[1] += msDN_DX(1,1)*p_avg;
			fv2[0] += msDN_DX(2,0)*p_avg; fv2[1] += msDN_DX(2,1)*p_avg;
			//fv -= grad_pn
			ms_temp_vec_np[0] = p0old;
			ms_temp_vec_np[1] = p1old;
			ms_temp_vec_np[2] = p2old;
			noalias(ms_vel_gauss) = prod(trans(msDN_DX),ms_temp_vec_np);
			ms_vel_gauss *= (Area);
			fv0[0] += (msN[0])*ms_vel_gauss[0]; fv0[1] += (msN[0])*ms_vel_gauss[1];
			fv1[0] += (msN[1])*ms_vel_gauss[0]; fv1[1] += (msN[1])*ms_vel_gauss[1];
			fv2[0] += (msN[2])*ms_vel_gauss[0]; fv2[1] += (msN[2])*ms_vel_gauss[1];
#else //G pn form
			double p_avg =	msN[0]*(p0 - p0old) 
							+ msN[1]*(p1 - p1old)
							+ msN[2]*(p2 - p2old) ;
			p_avg *= Area;
			fv0[0] += msDN_DX(0,0)*p_avg; fv0[1] += msDN_DX(0,1)*p_avg;
			fv1[0] += msDN_DX(1,0)*p_avg; fv1[1] += msDN_DX(1,1)*p_avg;
			fv2[0] += msDN_DX(2,0)*p_avg; fv2[1] += msDN_DX(2,1)*p_avg;
#endif



		}


		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void Fluid2DCoupled::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;

		unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber == 1) //step 1
		{
			if(rResult.size() != number_of_nodes*dim)
				rResult.resize(number_of_nodes*dim,false);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				rResult[i*dim] = GetGeometry()[i].GetDof(FRACT_VEL_X).EquationId();
				rResult[i*dim+1] = GetGeometry()[i].GetDof(FRACT_VEL_Y).EquationId();
			}
		}
		else if(FractionalStepNumber == 4) // pressure correction step
		{
			if(rResult.size() != number_of_nodes)
				rResult.resize(number_of_nodes,false);	
			for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void Fluid2DCoupled::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;

		unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber == 1) //step 1
		{
			if(ElementalDofList.size() != number_of_nodes*dim)
				ElementalDofList.resize(number_of_nodes*dim);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				ElementalDofList[i*dim] = GetGeometry()[i].pGetDof(FRACT_VEL_X);
				ElementalDofList[i*dim+1] = GetGeometry()[i].pGetDof(FRACT_VEL_Y);
			}
		}
		else if(FractionalStepNumber == 4) // pressure correction step
		{
			if(ElementalDofList.size() != number_of_nodes)
				ElementalDofList.resize(number_of_nodes);	
			
			for (unsigned int i=0;i<number_of_nodes;i++)
				ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);
		}
	}


	//************************************************************************************
	//************************************************************************************
	void Fluid2DCoupled::CalculateViscousMatrix(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double& nu)
	{
		K(0,0) = 2.0 * pow(DN_DX(0,0), 2) * nu + nu * pow(DN_DX(0,1), 2);
		K(0,1) = DN_DX(0,1) * nu * DN_DX(0,0);
		K(0,2) = 2.0 * DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,1) * nu * DN_DX(1,1);
		K(0,3) = DN_DX(0,1) * nu * DN_DX(1,0);
		K(0,4) = 2.0 * DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,1) * nu * DN_DX(2,1);
		K(0,5) = DN_DX(0,1) * nu * DN_DX(2,0);
//		K(1,0) = DN_DX(0,1) * nu * DN_DX(0,0);
		K(1,1) = 2.0 * nu * pow(DN_DX(0,1), 2) + pow(DN_DX(0,0), 2) * nu;
		K(1,2) = DN_DX(0,0) * nu * DN_DX(1,1);
		K(1,3) = 2.0 * DN_DX(0,1) * nu * DN_DX(1,1) + DN_DX(0,0) * nu * DN_DX(1,0);
		K(1,4) = DN_DX(0,0) * nu * DN_DX(2,1);
		K(1,5) = 2.0 * DN_DX(0,1) * nu * DN_DX(2,1) + DN_DX(0,0) * nu * DN_DX(2,0);
//		K(2,0) = 2.0 * DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,1) * nu * DN_DX(1,1);
//		K(2,1) = DN_DX(0,0) * nu * DN_DX(1,1);
		K(2,2) = 2.0 * pow(DN_DX(1,0), 2) * nu + nu * pow(DN_DX(1,1), 2);
		K(2,3) = DN_DX(1,1) * nu * DN_DX(1,0);
		K(2,4) = 2.0 * DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,1) * nu * DN_DX(2,1);
		K(2,5) = DN_DX(1,1) * nu * DN_DX(2,0);
//		K(3,0) = DN_DX(0,1) * nu * DN_DX(1,0);
//		K(3,1) = 2.0 * DN_DX(0,1) * nu * DN_DX(1,1) + DN_DX(0,0) * nu * DN_DX(1,0);
//		K(3,2) = DN_DX(1,1) * nu * DN_DX(1,0);
		K(3,3) = 2.0 * nu * pow(DN_DX(1,1), 2) + pow(DN_DX(1,0), 2) * nu;
		K(3,4) = DN_DX(1,0) * nu * DN_DX(2,1);
		K(3,5) = 2.0 * DN_DX(1,1) * nu * DN_DX(2,1) + DN_DX(1,0) * nu * DN_DX(2,0);
//		K(4,0) = 2.0 * DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,1) * nu * DN_DX(2,1);
//		K(4,1) = DN_DX(0,0) * nu * DN_DX(2,1);
//		K(4,2) = 2.0 * DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,1) * nu * DN_DX(2,1);
//		K(4,3) = DN_DX(1,0) * nu * DN_DX(2,1);
		K(4,4) = 2.0 * pow(DN_DX(2,0), 2) * nu + nu * pow(DN_DX(2,1), 2);
		K(4,5) = DN_DX(2,1) * nu * DN_DX(2,0);
//		K(5,0) = DN_DX(0,1) * nu * DN_DX(2,0);
//		K(5,1) = 2.0 * DN_DX(0,1) * nu * DN_DX(2,1) + DN_DX(0,0) * nu * DN_DX(2,0);
//		K(5,2) = DN_DX(1,1) * nu * DN_DX(2,0);
//		K(5,3) = 2.0 * DN_DX(1,1) * nu * DN_DX(2,1) + DN_DX(1,0) * nu * DN_DX(2,0);
//		K(5,4) = DN_DX(2,1) * nu * DN_DX(2,0);
		K(5,5) = 2.0 * nu * pow(DN_DX(2,1), 2) + pow(DN_DX(2,0), 2) * nu;	

		//filling the symmetric part
		for(unsigned int i = 1; i<K.size1(); i++)
			for(unsigned int j = 0; j<i; j++)
				K(i,j) = K(j,i);
	}

		//***********************************************************************
		//***********************************************************************
		//performs the Kroneker product of the Reduced Matrix with the identity matrix of 
		//size "dimension" ADDING to the destination matrix
		inline void  Fluid2DCoupled::ExpandAndAddReducedMatrix(
			MatrixType& Destination,
			boost::numeric::ublas::bounded_matrix<double,3,3>& ReducedMatrix,
			const unsigned int dimension)
		{
			KRATOS_TRY
			for (unsigned int i=0;i<ReducedMatrix.size2();i++)
			{
				int rowindex = i*dimension;
				for (unsigned int j=0;j<ReducedMatrix.size2();j++)
				{
					unsigned int colindex = j*dimension;
					for(unsigned int ii=0;ii<dimension;ii++)
						Destination(rowindex+ii,colindex+ii)+=ReducedMatrix(i,j);
				}
			}
			KRATOS_CATCH("")
		}

} // Namespace Kratos


