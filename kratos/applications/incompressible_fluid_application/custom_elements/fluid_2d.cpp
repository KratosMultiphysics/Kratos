/* b
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
//   Date:                $Date: 2009-01-13 15:39:56 $
//   Revision:            $Revision: 1.14 $
//
//
 
//#define GRADPN_FORM

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/fluid_2d.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{

	//************************************************************************************
	//************************************************************************************
	Fluid2D::Fluid2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Fluid2D::Fluid2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
		//setting up the nodal degrees of freedom
//		for(unsigned int i = 0 ; i != GetGeometry().size() ; ++i)
//		{
//			(GetGeometry()[i].pAddDof(FRACT_VEL_X));
//			(GetGeometry()[i].pAddDof(FRACT_VEL_Y));
//			(GetGeometry()[i].pAddDof(PRESSURE));
//		}


	}

	Element::Pointer Fluid2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		KRATOS_TRY

		return Element::Pointer(new Fluid2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	Fluid2D::~Fluid2D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void Fluid2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
	void Fluid2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}

	//************************************************************************************
	//************************************************************************************
	//calculation by component of the fractional step velocity corresponding to the first stage
	void Fluid2D::Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										   ProcessInfo& rCurrentProcessInfo, unsigned int ComponentIndex)
	{
		KRATOS_TRY;

		const unsigned int number_of_points = 3;

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false); //false says not to preserve existing storage!!

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false); //false says not to preserve existing storage!!
			
		boost::numeric::ublas::bounded_matrix<double,3,3> msMassFactors = 1.0/3.0*IdentityMatrix(3,3);
		boost::numeric::ublas::bounded_matrix<double,3,3> msWorkMatrix;
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
		array_1d<double,3> msN; 
		array_1d<double,2> ms_vel_gauss; 
		array_1d<double,3> ms_temp_vec_np; 
		array_1d<double,3> ms_u_DN; 
		array_1d<double,3> ms_aux0; 
		array_1d<double,3> ms_aux1; 
		array_1d<double,3> ms_aux2; 

		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

		//getting the velocity vector on the nodes

		//getting the velocity on the nodes
		const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL,0);
		const array_1d<double,3>& w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj0 = GetGeometry()[0].FastGetSolutionStepValue(CONV_PROJ);
		double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		const double fcomp0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];

		const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj1 = GetGeometry()[1].FastGetSolutionStepValue(CONV_PROJ);
		double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
		const double fcomp1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];

		const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj2 = GetGeometry()[2].FastGetSolutionStepValue(CONV_PROJ);
		double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
		const double fcomp2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];
		
		//====================================================================
		//calculating viscosity
		double nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
		double density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );

		//getting the BDF2 coefficients (not fixed to allow variable time step)
		//the coefficients INCLUDE the time step
		const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
		
		//VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
		// rLeftHandSideMatrix += Laplacian * nu; --> ONE GAUSS POINT
		noalias(rLeftHandSideMatrix) = nu * prod(msDN_DX,trans(msDN_DX));

		//INERTIA CONTRIBUTION
                		//filling the mass factors

		//  rLeftHandSideMatrix += M/Dt
		noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * msMassFactors;
		
		// *****************************************
		//CALCULATION OF THE RHS

		//external forces (component)
		double force_component = 0.3333333333333333*(fcomp0 + fcomp1 + fcomp2);

#if defined(GRADPN_FORM)
		//std::cout << "grad pn form " << std::endl;
		//calculating pressure grad (component)
		double p_grad   = msDN_DX(0,ComponentIndex)*p0old 
					+ msDN_DX(1,ComponentIndex)*p1old 
					+ msDN_DX(2,ComponentIndex)*p2old;
		p_grad /= density;
		// RHS = Fext - grad*pn
		noalias(rRightHandSideVector) = (force_component - p_grad)*msN;
#else 
		//adding pressure gradient (integrated by parts)
		noalias(rRightHandSideVector) = (force_component )*msN;
		double p_avg = msN[0]*p0old + msN[1]*p1old + msN[2]*p2old;
		p_avg /= density;
		rRightHandSideVector[0] += msDN_DX(0,ComponentIndex)*p_avg; 
		rRightHandSideVector[1] += msDN_DX(1,ComponentIndex)*p_avg; 
		rRightHandSideVector[2] += msDN_DX(2,ComponentIndex)*p_avg;
#endif

		//adding the inertia terms
		// RHS += M*vhistory 
		//calculating the historical velocity
		noalias(ms_temp_vec_np) = ZeroVector(3);		
		for(unsigned int step = 1; step<BDFcoeffs.size(); step++)
		{
			for(int iii = 0; iii<3; iii++)
			{
				const array_1d<double,3>& v = (GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY,step) );
				ms_temp_vec_np[iii] -= BDFcoeffs[step]*v[ComponentIndex];
			}
			
		}	
		noalias(rRightHandSideVector) += prod(msMassFactors,ms_temp_vec_np) ;
		
		//multiplying by area
		rLeftHandSideMatrix *= (Area * density);
		rRightHandSideVector *= (Area * density);
		
		
		
		//====================================================================
		//calculation of convective and stabilizing terms ... using 3 gauss points (on the sides)
//		double c1 = 4.00; double c2 = 2.00;
		double h = sqrt(2.00*Area/3.0);
		double norm_u = 0.0; double tau=0.0;
		double area_density_third = Area * density * 0.33333333333333333333;
		
		// ******************* GAUSS1 ***************************
		msN[0]=0.5; msN[1]=0.5; msN[2]=0.0;
		
		// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) 
		//(note that the fractional step vel is used)
		ms_vel_gauss[0] =  msN[0]*(fv0[0]-w0[0]) + msN[1]*(fv1[0]-w1[0]) +  msN[2]*(fv2[0]-w2[0]);
		ms_vel_gauss[1] =  msN[0]*(fv0[1]-w0[1]) + msN[1]*(fv1[1]-w1[1]) +  msN[2]*(fv2[1]-w2[1]);
// KRATOS_WATCH(ms_vel_gauss);

		//calculating parameter tau 		
		norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
		norm_u = sqrt(norm_u);

                tau = CalculateTau(h,nu,norm_u,rCurrentProcessInfo);

		//CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
		noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
		noalias(msWorkMatrix) = outer_prod(msN,ms_u_DN);

		//CONVECTION STABILIZING CONTRIBUTION (Suu)
		noalias(msWorkMatrix) += tau * outer_prod(ms_u_DN,ms_u_DN);
		
		//adding gauss point contribution
		noalias(rLeftHandSideMatrix) += area_density_third * msWorkMatrix;

		//RHS += Suy * proj[component] 
		double proj_component = msN[0]*proj0[ComponentIndex] 
					+ msN[1]*proj1[ComponentIndex] 
					+ msN[2]*proj2[ComponentIndex];
		noalias(rRightHandSideVector) += (area_density_third*tau*proj_component)*ms_u_DN;  
		 
		
		// ******************* GAUSS2 ***************************
		msN[0]=0.0; msN[1]=0.5; msN[2]=0.5;
		
		// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) 
		//(note that the fractional step vel is used)
		ms_vel_gauss[0] =  msN[0]*(fv0[0]-w0[0]) + msN[1]*(fv1[0]-w1[0]) +  msN[2]*(fv2[0]-w2[0]);
		ms_vel_gauss[1] =  msN[0]*(fv0[1]-w0[1]) + msN[1]*(fv1[1]-w1[1]) +  msN[2]*(fv2[1]-w2[1]);
// KRATOS_WATCH(ms_vel_gauss);

		//calculating parameter tau 		
		norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
		norm_u = sqrt(norm_u);

                tau = CalculateTau(h,nu,norm_u,rCurrentProcessInfo);

		//CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
		noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
		noalias(msWorkMatrix) = outer_prod(msN,ms_u_DN);

		//CONVECTION STABILIZING CONTRIBUTION (Suu)
		noalias(msWorkMatrix) += tau * outer_prod(ms_u_DN,ms_u_DN);
		
		//adding gauss point contribution
		noalias(rLeftHandSideMatrix) += area_density_third * msWorkMatrix;

		//RHS += Suy * proj[component] 
		proj_component = msN[0]*proj0[ComponentIndex] 
					+ msN[1]*proj1[ComponentIndex] 
					+ msN[2]*proj2[ComponentIndex];
		noalias(rRightHandSideVector) += (area_density_third*tau*proj_component)*ms_u_DN;  
		
		// ******************* GAUSS3 ***************************
		msN[0]=0.5; msN[1]=0.0; msN[2]=0.5;
		
		// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) 
		//(note that the fractional step vel is used)
		ms_vel_gauss[0] =  msN[0]*(fv0[0]-w0[0]) + msN[1]*(fv1[0]-w1[0]) +  msN[2]*(fv2[0]-w2[0]);
		ms_vel_gauss[1] =  msN[0]*(fv0[1]-w0[1]) + msN[1]*(fv1[1]-w1[1]) +  msN[2]*(fv2[1]-w2[1]);
// KRATOS_WATCH(ms_vel_gauss);

		//calculating parameter tau 		
		norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
		norm_u = sqrt(norm_u);


		tau = CalculateTau(h,nu,norm_u,rCurrentProcessInfo);
// KRATOS_WATCH(tau);

		//CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
		noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
		noalias(msWorkMatrix) = outer_prod(msN,ms_u_DN);

		//CONVECTION STABILIZING CONTRIBUTION (Suu)
		noalias(msWorkMatrix) += tau * outer_prod(ms_u_DN,ms_u_DN);
		
		//adding gauss point contribution
		noalias(rLeftHandSideMatrix) += area_density_third * msWorkMatrix;

		//RHS += Suy * proj[component] 
		proj_component = msN[0]*proj0[ComponentIndex] 
					+ msN[1]*proj1[ComponentIndex] 
					+ msN[2]*proj2[ComponentIndex];
		noalias(rRightHandSideVector) += (area_density_third*tau*proj_component)*ms_u_DN;  		
		
		

		//suubtracting the dirichlet term
		// RHS -= LHS*FracVel
		ms_temp_vec_np[0] = fv0[ComponentIndex]; 
		ms_temp_vec_np[1] = fv1[ComponentIndex]; 
		ms_temp_vec_np[2] = fv2[ComponentIndex]; 
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);
// rLeftHandSideMatrix *= density;
// rRightHandSideVector *= density;
		
// 		KRATOS_WATCH(rLeftHandSideMatrix);
// 		KRATOS_WATCH(rRightHandSideVector);

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	//calculation by component of the fractional step velocity corresponding to the first stage
	void Fluid2D::Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										   ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		unsigned int number_of_points = 3;

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);


		//boost::numeric::ublas::bounded_matrix<double,3,3> msMassFactors = 1.0/3.0*IdentityMatrix(3,3);
		boost::numeric::ublas::bounded_matrix<double,3,3> msWorkMatrix;
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
		array_1d<double,3> msN; 
		array_1d<double,2> ms_vel_gauss; 
		array_1d<double,3> ms_temp_vec_np; 
		array_1d<double,3> ms_u_DN; 
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

		
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
		double tau = CalculateTau(h,nu,norm_u,rCurrentProcessInfo);

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
//		GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
//		GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
//		GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;

                GetGeometry()[0].SetLock();
                GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
                GetGeometry()[0].UnSetLock();

                GetGeometry()[1].SetLock();
                GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
                GetGeometry()[1].UnSetLock();

                GetGeometry()[2].SetLock();
                GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
                GetGeometry()[2].UnSetLock();

		KRATOS_CATCH("");
	}
	  
	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void Fluid2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
 		boost::numeric::ublas::bounded_matrix<double,3,3> msMassFactors = 1.0/3.0*IdentityMatrix(3,3);
		boost::numeric::ublas::bounded_matrix<double,3,3> msWorkMatrix;
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
		array_1d<double,3> msN; 
		array_1d<double,2> ms_vel_gauss; 
		array_1d<double,3> ms_temp_vec_np; 
		array_1d<double,3> ms_u_DN; 
		array_1d<double,3> ms_aux0; 
		array_1d<double,3> ms_aux1; 
		array_1d<double,3> ms_aux2; 

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
			//note that here we calculate it "strong"
			ms_vel_gauss[0] = msDN_DX(0,0)*(p0) + msDN_DX(1,0)*(p1) + msDN_DX(2,0)*(p2);
			ms_vel_gauss[1] = msDN_DX(0,1)*(p0) + msDN_DX(1,1)*(p1) + msDN_DX(2,1)*(p2);
			ms_vel_gauss *= Area;

			//press_proj += G*p
                        GetGeometry()[0].SetLock();
			press_proj0[0] += msN[0]*ms_vel_gauss[0]; 
			press_proj0[1] += msN[0]*ms_vel_gauss[1]; 
                        GetGeometry()[0].UnSetLock();

                        GetGeometry()[1].SetLock();
			press_proj1[0] += msN[1]*ms_vel_gauss[0]; 
			press_proj1[1] += msN[1]*ms_vel_gauss[1]; 
                        GetGeometry()[1].UnSetLock();
                        
                        GetGeometry()[2].SetLock();
			press_proj2[0] += msN[2]*ms_vel_gauss[0]; 
			press_proj2[1] += msN[2]*ms_vel_gauss[1]; 
                        GetGeometry()[2].UnSetLock();
			
			//CONVECTIVE PROJECTION			
			
			// ******************* GAUSS1 ***************************
			msN[0]=0.5; msN[1]=0.5; msN[2]=0.0;
		
			// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) 
			//(note that the fractional step vel is used)
			ms_vel_gauss[0] =  msN[0]*(fv0[0]-w0[0]) + msN[1]*(fv1[0]-w1[0]) +  msN[2]*(fv2[0]-w2[0]);
			ms_vel_gauss[1] =  msN[0]*(fv0[1]-w0[1]) + msN[1]*(fv1[1]-w1[1]) +  msN[2]*(fv2[1]-w2[1]);

			//calculating convective auxiliary vector
			noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);

			//attention changing the meaning of ms_vel_gauss!!
			ms_vel_gauss[0] = ms_u_DN[0] * fv0[0] + ms_u_DN[1] * fv1[0] + ms_u_DN[2] * fv2[0];
			ms_vel_gauss[1] = ms_u_DN[0] * fv0[1] + ms_u_DN[1] * fv1[1] + ms_u_DN[2] * fv2[1];
			
			// conv_proj += C*u
			ms_aux0[0] = msN[0]*ms_vel_gauss[0]; 
			ms_aux0[1] = msN[0]*ms_vel_gauss[1];

			ms_aux1[0] = msN[1]*ms_vel_gauss[0]; 
			ms_aux1[1] = msN[1]*ms_vel_gauss[1];

			ms_aux2[0] = msN[2]*ms_vel_gauss[0];
			ms_aux2[1] = msN[2]*ms_vel_gauss[1];
			
			
			// ******************* GAUSS2 ***************************
			msN[0]=0.0; msN[1]=0.5; msN[2]=0.5;
		
			// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) 
			//(note that the fractional step vel is used)
			ms_vel_gauss[0] =  msN[0]*(fv0[0]-w0[0]) + msN[1]*(fv1[0]-w1[0]) +  msN[2]*(fv2[0]-w2[0]);
			ms_vel_gauss[1] =  msN[0]*(fv0[1]-w0[1]) + msN[1]*(fv1[1]-w1[1]) +  msN[2]*(fv2[1]-w2[1]);

			//calculating convective auxiliary vector
			noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);

			//attention changing the meaning of ms_vel_gauss!!
			ms_vel_gauss[0] = ms_u_DN[0] * fv0[0] + ms_u_DN[1] * fv1[0] + ms_u_DN[2] * fv2[0];
			ms_vel_gauss[1] = ms_u_DN[0] * fv0[1] + ms_u_DN[1] * fv1[1] + ms_u_DN[2] * fv2[1];
			
			// conv_proj += C*u
			ms_aux0[0] += msN[0]*ms_vel_gauss[0]; 
			ms_aux0[1] += msN[0]*ms_vel_gauss[1];

			ms_aux1[0] += msN[1]*ms_vel_gauss[0]; 
			ms_aux1[1] += msN[1]*ms_vel_gauss[1];

			ms_aux2[0] += msN[2]*ms_vel_gauss[0];
			ms_aux2[1] += msN[2]*ms_vel_gauss[1];
						
			// ******************* GAUSS3 ***************************
			msN[0]=0.5; msN[1]=0.0; msN[2]=0.5;
		
			// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) 
			//(note that the fractional step vel is used)
			ms_vel_gauss[0] =  msN[0]*(fv0[0]-w0[0]) + msN[1]*(fv1[0]-w1[0]) +  msN[2]*(fv2[0]-w2[0]);
			ms_vel_gauss[1] =  msN[0]*(fv0[1]-w0[1]) + msN[1]*(fv1[1]-w1[1]) +  msN[2]*(fv2[1]-w2[1]);

			//calculating convective auxiliary vector
			noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);

			//attention changing the meaning of ms_vel_gauss!!
			ms_vel_gauss[0] = ms_u_DN[0] * fv0[0] + ms_u_DN[1] * fv1[0] + ms_u_DN[2] * fv2[0];
			ms_vel_gauss[1] = ms_u_DN[0] * fv0[1] + ms_u_DN[1] * fv1[1] + ms_u_DN[2] * fv2[1];
			
			// conv_proj += C*u
			ms_aux0[0] += msN[0]*ms_vel_gauss[0]; 
			ms_aux0[1] += msN[0]*ms_vel_gauss[1];

			ms_aux1[0] += msN[1]*ms_vel_gauss[0]; 
			ms_aux1[1] += msN[1]*ms_vel_gauss[1];

			ms_aux2[0] += msN[2]*ms_vel_gauss[0];
			ms_aux2[1] += msN[2]*ms_vel_gauss[1];
			
			//access to database

			// conv_proj += C*u
			double temp = Area*density*0.33333333333333333333333333;

                        GetGeometry()[0].SetLock();
			conv_proj0[0] += temp*ms_aux0[0]; 
			conv_proj0[1] +=  temp*ms_aux0[1];
                        GetGeometry()[0].UnSetLock();

                        GetGeometry()[1].SetLock();
			conv_proj1[0] +=  temp*ms_aux1[0]; 
			conv_proj1[1] +=  temp*ms_aux1[1];
                        GetGeometry()[1].UnSetLock();

                        GetGeometry()[2].SetLock();
			conv_proj2[0] +=  temp*ms_aux2[0];
			conv_proj2[1] +=  temp*ms_aux2[1];
                        GetGeometry()[2].UnSetLock();
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
	void Fluid2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
	  void Fluid2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
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
	inline double Fluid2D::CalculateTau(const double h, const double nu, const double norm_u, ProcessInfo& CurrentProcessInfo)
	{
              const double c1 = 4.00;
              const double c2 = 2.00;

              const double dyn_st_beta = CurrentProcessInfo[DYNAMIC_TAU];

              const double inv_dt_coeff = CurrentProcessInfo[BDF_COEFFICIENTS][0];
              double tau = 1.00 / (dyn_st_beta*inv_dt_coeff +  c1*nu/(h*h) + c2*norm_u/h );

//                if (dyn_st_switch)
//                {
//                    const double inv_dt_coeff = CurrentProcessInfo[BDF_COEFFICIENTS][0];
//                    tau = 1.00 / (inv_dt_coeff +  c1*nu/(h*h) + c2*norm_u/h );
////                    KRATOS_WATCH(tau);
//                }
//                else
//                {
//                    tau = 1.00 / (c1*nu/(h*h) + c2*norm_u/h );
//                }

                return tau;

	}


//#undef GRADPN_FORM
} // Namespace Kratos


