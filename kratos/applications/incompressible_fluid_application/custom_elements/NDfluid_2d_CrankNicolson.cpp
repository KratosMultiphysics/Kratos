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
//   Last modified by:    $Author: antonia $
//   Date:                $Date: 2008-07-25 07:31:19 $
//   Revision:            $Revision: 1.2 $
//
//
 
//#define GRADPN_FORM
//#define STOKES

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/NDfluid_2d_CrankNicolson.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
    namespace NDFluid2DCrankNicolsonauxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,3,3> msMassFactors = 1.0/3.0*IdentityMatrix(3,3);

        boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX = ZeroMatrix(3,2);

        array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes

        array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension

	array_1d<double,2> ms_vel_gauss_old = ZeroVector(2); //dimesion coincides with space dimension

        array_1d<double,3> ms_temp_vec_np = ZeroVector(3); //dimension = number of nodes

        array_1d<double,3> ms_u_DN = ZeroVector(3); //dimension = number of nodes
    }
    using  namespace NDFluid2DCrankNicolsonauxiliaries;


	//************************************************************************************
	//************************************************************************************
	NDFluid2DCrankNicolson::NDFluid2DCrankNicolson(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	NDFluid2DCrankNicolson::NDFluid2DCrankNicolson(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
		//setting up the nodal degrees of freedom
//		for(unsigned int i = 0 ; i != GetGeometry().size() ; ++i)
//		{
//			(GetGeometry()[i].pAddDof(FRACT_VEL_X));
//			(GetGeometry()[i].pAddDof(FRACT_VEL_Y));
//			(GetGeometry()[i].pAddDof(PRESSURE));
//		}

		//filling the mass factors
		//msMassFactors(0,0) = 1.00/6.00;  msMassFactors(0,1) = 1.00/12.00; msMassFactors(0,2) = 1.00/12.00;
		//msMassFactors(1,0) = 1.00/12.00; msMassFactors(1,1) = 1.00/6.00;  msMassFactors(1,2) = 1.00/12.00;
		//msMassFactors(2,0) = 1.00/12.00; msMassFactors(2,1) = 1.00/12.00; msMassFactors(2,2) = 1.00/6.00;
		
//		msMassFactors(0,0) = 1.00/3.00; msMassFactors(0,1) = 0.00;		msMassFactors(0,2) = 0.00;
//		msMassFactors(1,0) = 0.00;		msMassFactors(1,1) = 1.00/3.00; msMassFactors(1,2) = 0.00;
//		msMassFactors(2,0) = 0.00;		msMassFactors(2,1) = 0.00;		msMassFactors(2,2) = 1.00/3.00;		
	}

	Element::Pointer NDFluid2DCrankNicolson::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		KRATOS_TRY

		return Element::Pointer(new NDFluid2DCrankNicolson(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	NDFluid2DCrankNicolson::~NDFluid2DCrankNicolson()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void NDFluid2DCrankNicolson::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber < 3) //first step of the fractional step solution
		{
			int ComponentIndex = FractionalStepNumber - 1;
			Stage1(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo, ComponentIndex);
		}
		else if (FractionalStepNumber == 4)//second step of the fractional step solution
		{
			Stage2(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
		}

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void NDFluid2DCrankNicolson::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}

	//************************************************************************************
	//************************************************************************************
	//calculation by component of the fractional step velocity corresponding to the first stage
	void NDFluid2DCrankNicolson::Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										   ProcessInfo& rCurrentProcessInfo, unsigned int ComponentIndex)
	{
		KRATOS_TRY;

		const unsigned int number_of_points = 3;

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false); //false says not to preserve existing storage!!

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false); //false says not to preserve existing storage!!

		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

		//getting the velocity vector on the nodes

		//getting the velocity on the nodes
		const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL,0);
		const array_1d<double,3>& fv0_old = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
		const array_1d<double,3>& w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& w0_old = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY,1);
		const array_1d<double,3>& proj0 = GetGeometry()[0].FastGetSolutionStepValue(CONV_PROJ);
		const array_1d<double,3>& proj0_old = GetGeometry()[0].FastGetSolutionStepValue(CONV_PROJ,1);
		double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		const double fcomp0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];
		const double eps0 = GetGeometry()[0].FastGetSolutionStepValue(POROSITY);
		const double dp0 = GetGeometry()[0].FastGetSolutionStepValue(DIAMETER);
		

		const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& fv1_old = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
		const array_1d<double,3>& w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& w1_old = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY,1);
		const array_1d<double,3>& proj1 = GetGeometry()[1].FastGetSolutionStepValue(CONV_PROJ);
		double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const array_1d<double,3>& proj1_old = GetGeometry()[1].FastGetSolutionStepValue(CONV_PROJ,1);
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
		const double fcomp1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];
		const double eps1 = GetGeometry()[1].FastGetSolutionStepValue(POROSITY);
		const double dp1 = GetGeometry()[1].FastGetSolutionStepValue(DIAMETER);
		

		const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& fv2_old = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
		const array_1d<double,3>& w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& w2_old = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY,1);
		const array_1d<double,3>& proj2 = GetGeometry()[2].FastGetSolutionStepValue(CONV_PROJ);
		const array_1d<double,3>& proj2_old = GetGeometry()[2].FastGetSolutionStepValue(CONV_PROJ,1);
		double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
		const double fcomp2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];
		const double eps2 = GetGeometry()[2].FastGetSolutionStepValue(POROSITY);
		const double dp2 = GetGeometry()[2].FastGetSolutionStepValue(DIAMETER);


		// 
		// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
		ms_vel_gauss[0] =  msN[0]*(fv0[0]-w0[0]) + msN[1]*(fv1[0]-w1[0]) +  msN[2]*(fv2[0]-w2[0]);
		ms_vel_gauss[1] =  msN[0]*(fv0[1]-w0[1]) + msN[1]*(fv1[1]-w1[1]) +  msN[2]*(fv2[1]-w2[1]);

		//vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
		ms_vel_gauss_old[0] =  msN[0]*(fv0_old[0]-w0_old[0]) + msN[1]*(fv1_old[0]-w1_old[0]) +  msN[2]*(fv2_old[0]-w2_old[0]);
		ms_vel_gauss_old[1] =  msN[0]*(fv0_old[1]-w0_old[1]) + msN[1]*(fv1_old[1]-w1_old[1]) +  msN[2]*(fv2_old[1]-w2_old[1]);		


		//ms_vel_gauss = v at (n+1)/2;
		ms_vel_gauss[0] += ms_vel_gauss_old[0];
		ms_vel_gauss[0] *= 0.5;
		ms_vel_gauss[1] += ms_vel_gauss_old[1];
		ms_vel_gauss[1] *= 0.5;

		//calculating viscosity
		double nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
 		double density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );

		//DIAMETER of the element
		double dp = 0.3333333333333333333333*(dp0 + dp1 + dp2);
		//POROSITY of the element: average value of the porosity
		double eps = 0.3333333333333333333333*(eps0 + eps1 + eps2 );

		//1/PERMEABILITY of the element: average value of the porosity
		double kinv = 0.0;

		//Calculate the elemental Kinv in function of the nodal K of each element.

		//Version 1: we can calculate the elemental kinv from the nodal kinv;
		//THERE IS AN ERROR: IN THE INTERPHASE ELEMENTS A WATER NODE HAS TO BE ''MORE IMPORTANT'' THAN A POROUS ONE!!!
// 		if(kinv0 != 0.0 || kinv1 != 0.0 || kinv2 != 0.0) //if it is a fluid element
// 		{	double k0 = 0.0;
// 			double k1 = 0.0;
// 			double k2 = 0.0;
// 			if(kinv0 != 0.0)	
// 				k0 = 1.0/kinv0;
// 			if(kinv1 != 0.0)
// 				k1 = 1.0/kinv1;
// 			if(kinv2 != 0.0)	
// 				k2 = 1.0/kinv2;
// 			kinv = 3.0/(k0 + k1 + k2 );
// 		}
// 		
		//Calculate kinv = 1/ k(eps_elem);
// if(rLeftHandSideMatrix.size1() != number_of_points)
// 			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
		
		//Version 2:we can calculate the elemental kinv from the already calculate elemental eps;
		if( (eps0 != 1) | (eps1 != 1) | (eps2 != 1) )
			kinv = 150*(1-eps)*(1-eps)/(eps*eps*eps*dp*dp);
		
		//getting the BDF2 coefficients (not fixed to allow variable time step)
		//the coefficients INCLUDE the time step
		//const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS]; 


		//Getting delta time value for the restriction over tau.
		double Dt = rCurrentProcessInfo[DELTA_TIME];
		
		array_1d<double,2> BDFcoeffs = ZeroVector(2);
		BDFcoeffs[0]= 1.0 / Dt; 		//coeff for step n+1;
		BDFcoeffs[1]= -1.0 / Dt;		//coeff for step n;




		//calculating parameter tau (saved internally to each element)
		double c1 = 4.00;
		double c2 = 2.00;
		double h = sqrt(2.00*Area);

		double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
		norm_u = sqrt(norm_u); //norm_u calculated at (n+1)/2;//
		double tau = 1.00 / ( c1*nu/(h*h) + c2*norm_u/h );


		

		// *****************************************
		//CALCULATION OF THE LHS

		//CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
		noalias(ms_u_DN) =  prod(msDN_DX , ms_vel_gauss);
		noalias(rLeftHandSideMatrix) = 0.5 * outer_prod(msN,ms_u_DN)/(eps*eps);

 		//CONVECTION STABILIZING CONTRIBUTION (Suu)
 		noalias(rLeftHandSideMatrix) += 0.5 * tau/(eps*eps) * outer_prod(ms_u_DN,ms_u_DN); 

		//VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
		// rLeftHandSideMatrix += Laplacian * nu;
		noalias(rLeftHandSideMatrix) += 0.5 * nu/eps * prod(msDN_DX,trans(msDN_DX));

		//INERTIA CONTRIBUTION
		//  rLeftHandSideMatrix += M/Dt
		noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * msMassFactors/eps;

 		//DARCY linear CONTRIBUTION
		//  rLeftHandSideMatrix -= nu/permeability (using the cinematic viscosity it is already divided for density: then we are going to multiplicate for density again);
 		noalias(rLeftHandSideMatrix) -= 0.5 * nu*msMassFactors*kinv;

 		//DARCY non linear CONTRIBUTION (brinkmann)
		//rLeftHandSideMatrix -= 1.75*|u(n+1/2)|/[(150*k)^0.5*eps^(3/2)]
 		noalias(rLeftHandSideMatrix) -= 0.5 * msMassFactors*norm_u*1.75*sqrt(kinv)/12.2474487/sqrt(eps*eps*eps);

		//multiplication by the area
		rLeftHandSideMatrix *= (Area * density);
		

		// *****************************************
		//CALCULATION OF THE RHS

		//external forces (component)
		double force_component = 0.3333333333333333*(fcomp0 + fcomp1 + fcomp2);

// 		KRATOS_WATCH(force_component);
// 		KRATOS_WATCH(p0old);
// 		KRATOS_WATCH(p1old);
// 		KRATOS_WATCH(p2old);

		//adding pressure gradient (integrated by parts)
		noalias(rRightHandSideVector) = (force_component )*msN;

//3PG-------------
//p_avg turn out to be p0_avg p1_avg and p2_avg
//3PG-------------

		double p_avg = msN[0]*p0old + msN[1]*p1old + msN[2]*p2old;
		p_avg /= density;

// 		KRATOS_WATCH(p_avg);

		rRightHandSideVector[0] += msDN_DX(0,ComponentIndex)*p_avg; 
		rRightHandSideVector[1] += msDN_DX(1,ComponentIndex)*p_avg; 
		rRightHandSideVector[2] += msDN_DX(2,ComponentIndex)*p_avg;

// 		KRATOS_WATCH(rRightHandSideVector);
		//adding the inertia terms
		// RHS += M*vhistory
		//calculating the historical velocity
		noalias(ms_temp_vec_np) = ZeroVector(3);		

		for(int iii = 0; iii<3; iii++)
		{
			const array_1d<double,3>& v = (GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY,1) );
			ms_temp_vec_np[iii] = BDFcoeffs[1]*v[ComponentIndex];
		}
	
		noalias(rRightHandSideVector) -= prod(msMassFactors,ms_temp_vec_np)/eps ;
// 		KRATOS_WATCH(prod(msMassFactors,ms_temp_vec_np)/eps);

//3PG-------------
//proj_component calculated on the 3 gauss points
//3PG-------------


		//RHS += Suy * proj[component] 
		double proj_component = msN[0]*proj0[ComponentIndex] 
							  + msN[1]*proj1[ComponentIndex] 
							  + msN[2]*proj2[ComponentIndex];
		double proj_old_component = msN[0]*proj0_old[ComponentIndex] 
							  + msN[1]*proj1_old[ComponentIndex] 
							  + msN[2]*proj2_old[ComponentIndex];
		proj_component += proj_old_component;
		proj_component *= 0.5;    //proj_component calculate in t_n+1/2;
 		noalias(rRightHandSideVector) += (tau*proj_component)/(eps*eps)*ms_u_DN;   

		//multiplying by area
		rRightHandSideVector  *= (Area * density);

		
		ms_temp_vec_np[0] = fv0_old[ComponentIndex]; 
		ms_temp_vec_np[1] = fv1_old[ComponentIndex]; 
		ms_temp_vec_np[2] = fv2_old[ComponentIndex]; 


		//there is a part of the lhs which is already included;
		for(int iii = 0; iii<3; iii++)
		{
			ms_temp_vec_np[iii] *= BDFcoeffs[0];
		}
		
		noalias(rRightHandSideVector) +=  (Area * density)*prod(msMassFactors,ms_temp_vec_np)/eps ;	
// 		KRATOS_WATCH((Area * density)*prod(msMassFactors,ms_temp_vec_np)/eps);

		//suubtracting the dirichlet term
		// RHS -= LHS*FracVel
		ms_temp_vec_np[0] = fv0[ComponentIndex] + fv0_old[ComponentIndex]; 
		ms_temp_vec_np[1] = fv1[ComponentIndex] + fv1_old[ComponentIndex];
		ms_temp_vec_np[2] = fv2[ComponentIndex] + fv2_old[ComponentIndex]; 

		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);
		
// 		KRATOS_WATCH(prod(rLeftHandSideMatrix,ms_temp_vec_np));
// 		
// 		KRATOS_WATCH(fv0); 
// 		KRATOS_WATCH(fv1); 
// 		KRATOS_WATCH(fv2); 
// 		KRATOS_WATCH(fv0_old); 
// 		KRATOS_WATCH(fv1_old); 
// 		KRATOS_WATCH(fv2_old); 	


// 		KRATOS_WATCH(rRightHandSideVector); 	

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	//calculation by component of the fractional step velocity corresponding to the first stage
	void NDFluid2DCrankNicolson::Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
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
		GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);


		const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
		double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE); 
		double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT); 
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		const double eps0 = GetGeometry()[0].FastGetSolutionStepValue(POROSITY);
	
		
		const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
		double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE); 
		double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT); 
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
		const double eps1 = GetGeometry()[1].FastGetSolutionStepValue(POROSITY);

		
		const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);
		double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE); 
		double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT); 
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
		const double eps2 = GetGeometry()[2].FastGetSolutionStepValue(POROSITY);

//3PG-------------
//ms_vel_gauss
//3PG-------------
		

		// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
		ms_vel_gauss[0] =  msN[0]*(fv0[0]-w0[0]) + msN[1]*(fv1[0]-w1[0]) +  msN[2]*(fv2[0]-w2[0]);
		ms_vel_gauss[1] =  msN[0]*(fv0[1]-w0[1]) + msN[1]*(fv1[1]-w1[1]) +  msN[2]*(fv2[1]-w2[1]);

		//calculating convective auxiliary vector
//3PG-------------
//ms_u_DN
//3PG-------------
		noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
		
		//calculating average density and viscosity
		double nu = 0.33333333333333*(nu0 + nu1 + nu2 );
 		double density = 0.33333333333333*(rho0 + rho1 + rho2 );

		//DIAMETER of the element
		//double dp = 0.33333333333333*(d0+d1+d2)
		//POROSITY of the element: average value of the porosity
		double eps = 0.3333333333333333333333*(eps0 + eps1 + eps2 );
		
		//Getting delta time value for the restriction over tau.
		double Dt = rCurrentProcessInfo[DELTA_TIME];
		
		//calculating parameter tau (saved internally to each element)
		double h = sqrt(2.00*Area);
		double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
		norm_u = sqrt(norm_u);
		double tau = 1.00 / ( 4.00*nu/(h*h) + 2.00*norm_u/h );

		//tau = min{tau; Dt}
		if(tau > Dt)
		{
			tau = Dt;	
		}

		//getting the BDF2 coefficients (not fixed to allow variable time step)
		//the coefficients INCLUDE the time step
		const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

		//CALCULATION OF THE LEFT HAND SIDE
		//laplacian term	       L = Dt * gradN  *eps * trans(gradN);
		//stabilization term       Spp = tau * gradN * eps * trans(gradN);
		noalias(rLeftHandSideMatrix) = ((1.00/BDFcoeffs[0] + tau)/density*eps) * prod(msDN_DX,trans(msDN_DX));

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
		//Inserting the influence of POROSITY!
		// RHS += Sz * proj*eps 
		//contrib of Spy*proj
		ms_vel_gauss[0] = msN[0]*proj0[0] + msN[1]*proj1[0] + msN[2]*proj2[0];
		ms_vel_gauss[1] = msN[0]*proj0[1] + msN[1]*proj1[1] + msN[2]*proj2[1];
		ms_vel_gauss *= tau*eps;
		noalias(rRightHandSideVector) += prod(msDN_DX , ms_vel_gauss);   

		//dirichlet contribution
		ms_temp_vec_np[0] = p0; 
		ms_temp_vec_np[1] = p1; 
		ms_temp_vec_np[2] = p2; 
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);

		// RHS += dt * L *eps* pold 
		ms_temp_vec_np[0] = p0old; 
		ms_temp_vec_np[1] = p1old; 
		ms_temp_vec_np[2] = p2old; 
		noalias(ms_vel_gauss) = prod(trans(msDN_DX),ms_temp_vec_np);
		noalias(rRightHandSideVector) += (1.00/BDFcoeffs[0]/density*eps) * prod(msDN_DX,ms_vel_gauss); 

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
	void NDFluid2DCrankNicolson::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
 
		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);
		
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
// 			const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
// 			const double eps0 = GetGeometry()[0].FastGetSolutionStepValue(POROSITY);

			double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
			double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
			array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
// 			const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
// 			const double eps1 = GetGeometry()[1].FastGetSolutionStepValue(POROSITY);


			double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
			double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
			array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
// 			const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
// 			const double eps2 = GetGeometry()[2].FastGetSolutionStepValue(POROSITY);

// 			double density = 0.33333333333333333333333*(rho0 + rho1 + rho2);
// 			double eps = 0.33333333333333333333333*(eps0 + eps1 + eps2);


			double p_avg =	msN[0]*(p0 - p0old) 
							+ msN[1]*(p1 - p1old)
							+ msN[2]*(p2 - p2old) ;
			p_avg *= Area;
// 			KRATOS_WATCH(p_avg);
// 			KRATOS_WATCH(eps);
			fv0[0] += msDN_DX(0,0)*p_avg; fv0[1] += msDN_DX(0,1)*p_avg;
			fv1[0] += msDN_DX(1,0)*p_avg; fv1[1] += msDN_DX(1,1)*p_avg;
			fv2[0] += msDN_DX(2,0)*p_avg; fv2[1] += msDN_DX(2,1)*p_avg;

			


		}


		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void NDFluid2DCrankNicolson::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber == 1) //step 1
			for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_X).EquationId();
		else if(FractionalStepNumber == 2) //step 2
			for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_Y).EquationId();

		else if(FractionalStepNumber == 4) // pressure correction step
			for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
		}

	//************************************************************************************
	//************************************************************************************
	  void NDFluid2DCrankNicolson::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber == 1) //step 1
			for (unsigned int i=0;i<number_of_nodes;i++)
				ElementalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_X);
		else if(FractionalStepNumber == 2) //step 2
			for (unsigned int i=0;i<number_of_nodes;i++)
				ElementalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_Y);
		else if(FractionalStepNumber == 4) // pressure correction step
			for (unsigned int i=0;i<number_of_nodes;i++)
				ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);

	}

	//************************************************************************************
	//************************************************************************************
/*	inline void NDFluid2DCrankNicolson::CalculateGeometryData(boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, array_1d<double,3>& N, double& Area)
	{
		double det_J = GetGeometry()[1].X() * GetGeometry()[2].Y()
	               - GetGeometry()[1].X() * GetGeometry()[0].Y()
	               - GetGeometry()[0].X() * GetGeometry()[2].Y()
 				   - GetGeometry()[1].Y() * GetGeometry()[2].X() 
 	               + GetGeometry()[1].Y() * GetGeometry()[0].X() 
	               + GetGeometry()[0].Y() * GetGeometry()[2].X();

		double temp = 1.00/det_J;

		Area = fabs(det_J) * 0.5;
		
		// X derivatives
		DN_DX(0, 0) = (GetGeometry()[1].Y() - GetGeometry()[2].Y()) *temp;	
		DN_DX(1, 0) = (GetGeometry()[2].Y() - GetGeometry()[0].Y()) *temp;
		DN_DX(2, 0) = (GetGeometry()[0].Y() - GetGeometry()[1].Y()) *temp;
		// Y derivatives
		DN_DX(0, 1) = (GetGeometry()[2].X() - GetGeometry()[1].X()) *temp;
		DN_DX(1, 1) = (GetGeometry()[0].X() - GetGeometry()[2].X()) *temp;
		DN_DX(2, 1) = (GetGeometry()[1].X() - GetGeometry()[0].X()) *temp;

		//shape functions
		N[0] = 0.333333333333333333333333333;
		N[1] = 0.333333333333333333333333333;
		N[2] = 0.333333333333333333333333333;
	}
*/
//#undef GRADPN_FORM
} // Namespace Kratos


