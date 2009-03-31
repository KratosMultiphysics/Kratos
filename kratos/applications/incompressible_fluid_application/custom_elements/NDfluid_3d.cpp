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
//   Date:                $Date: 2008-03-11 18:20:02 $
//   Revision:            $Revision: 1.3 $
//
//


// System includes 
//#define GRADPN_FORM //the grad(pn) is used instead of the G(pn) in doing the splitting
 
// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/NDfluid_3d.h"
#include "incompressible_fluid_application.h"
#include "utilities/math_utils.h" 
#include "utilities/geometry_utilities.h" 

namespace Kratos 
{
    //space defined to allow parallelism
    namespace NDFluid3DAuxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,4,4> msMassFactors = 0.25*IdentityMatrix(4,4);
//         #pragma omp threadprivate(msMassFactors)

        boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX = ZeroMatrix(4,3);
//         #pragma omp threadprivate(msDN_DX)

        array_1d<double,4> msN = ZeroVector(4); //dimension = number of nodes
//         #pragma omp threadprivate(msN)

        array_1d<double,3> ms_aux = ZeroVector(3); //dimension coincides with space dimension
//         #pragma omp threadprivate(ms_aux)

        array_1d<double,3> ms_vel_gauss = ZeroVector(3); //dimesion coincides with space dimension
//         #pragma omp threadprivate(ms_vel_gauss)

        array_1d<double,4> ms_temp_vec_np = ZeroVector(4); //dimension = number of nodes
//         #pragma omp threadprivate(ms_temp_vec_np)

        array_1d<double,4> ms_u_DN = ZeroVector(4); //dimension = number of nodes
//         #pragma omp threadprivate(ms_u_DN)
 
    }
	using namespace NDFluid3DAuxiliaries;
    

	
		


	//************************************************************************************
	//************************************************************************************
	NDFluid3D::NDFluid3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	NDFluid3D::NDFluid3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer NDFluid3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		KRATOS_TRY
		return Element::Pointer(new NDFluid3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	NDFluid3D::~NDFluid3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void NDFluid3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber <= 3) //first step of the fractional step solution
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
	void NDFluid3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}

	//************************************************************************************
	//************************************************************************************
	//calculation by component of the fractional step velocity corresponding to the first stage
	void NDFluid3D::Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										   ProcessInfo& rCurrentProcessInfo, unsigned int ComponentIndex)
	{
		KRATOS_TRY;

		unsigned int number_of_points = 4;

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);

		//getting data for the given geometry
		double Volume;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Volume);
		//CalculateGeometryData(msDN_DX,msN,Volume);

		//getting the velocity vector on the nodes

		//getting the velocity on the nodes
		const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL,0);
		const array_1d<double,3>& w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj0 = GetGeometry()[0].FastGetSolutionStepValue(CONV_PROJ);
		const double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		const double fcomp0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];
		const double eps0 = GetGeometry()[0].FastGetSolutionStepValue(POROSITY);
// 		const double kinv0 = GetGeometry()[0].FastGetSolutionStepValue(PERMEABILITY_INV);
		const double dp0 = GetGeometry()[0].FastGetSolutionStepValue(DIAMETER);

		const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj1 = GetGeometry()[1].FastGetSolutionStepValue(CONV_PROJ);
		const double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
		const double fcomp1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];
		const double eps1 = GetGeometry()[1].FastGetSolutionStepValue(POROSITY);
// 		const double kinv1 = GetGeometry()[1].FastGetSolutionStepValue(PERMEABILITY_INV);
		const double dp1 = GetGeometry()[1].FastGetSolutionStepValue(DIAMETER);

		const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj2 = GetGeometry()[2].FastGetSolutionStepValue(CONV_PROJ);
		const double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
		const double fcomp2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];
		const double eps2 = GetGeometry()[2].FastGetSolutionStepValue(POROSITY);
// 		const double kin2 = GetGeometry()[2].FastGetSolutionStepValue(PERMEABILITY_INV);
		const double dp2 = GetGeometry()[2].FastGetSolutionStepValue(DIAMETER);

		const array_1d<double,3>& fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj3 = GetGeometry()[3].FastGetSolutionStepValue(CONV_PROJ);
		const double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
		const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);
		const double fcomp3 = GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];
		const double eps3 = GetGeometry()[3].FastGetSolutionStepValue(POROSITY);
// 		const double kinv3 = GetGeometry()[3].FastGetSolutionStepValue(PERMEABILITY_INV);
		const double dp3 = GetGeometry()[3].FastGetSolutionStepValue(DIAMETER);

		// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
		noalias(ms_aux) = fv0;	noalias(ms_aux) -= w0;	
		noalias(ms_vel_gauss) = msN[0]*ms_aux;

		noalias(ms_aux) = fv1;	noalias(ms_aux) -= w1;	
		noalias(ms_vel_gauss) += msN[1]*ms_aux;

		noalias(ms_aux) = fv2;	noalias(ms_aux) -= w2;	
		noalias(ms_vel_gauss) += msN[2]*ms_aux;

		noalias(ms_aux) = fv3;	noalias(ms_aux) -= w3;	
		noalias(ms_vel_gauss) += msN[3]*ms_aux;

		//calculating viscosity
		double nu = 0.25*(nu0 + nu1 + nu2 + nu3);
 		double density = 0.25*(rho0 + rho1 + rho2 + rho3);

		//DIAMETER of the element
		double dp = 0.25*(dp0 + dp1 + dp2 + dp3);
// 		//POROSITY of the element: average value of the porosity of each node;
		double eps = 0.25*(eps0 + eps1 + eps2 + eps3);

		//1/PERMEABILITY of the element: average value of the porosity
		double kinv = 0.0;
		//we can calculate the elemental kinv from the already calculate elemental eps;
		//If the nodes of the element are fluid one then kinv= 0, otherwise kinv = 150(1-eps)²/eps³*dp²;
		if((eps0 != 1) || (eps1 != 1) || (eps2 != 1) || (eps3 != 1) )
			kinv = 150*(1-eps)*(1-eps)/(eps*eps*eps*dp*dp);

		//getting the BDF2 coefficients (not fixed to allow variable time step)
		//the coefficients INCLUDE the time step
		const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

		//calculating parameter tau (saved internally to each element)
		double c1 = 4.00;
		double c2 = 2.00;
		double h = CalculateH(Volume);
		//double h = pow(6.00*Volume,0.3333333);
		double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1] + ms_vel_gauss[2]*ms_vel_gauss[2];
		norm_u = sqrt(norm_u);
		double tau = 1.00 / ( c1*nu/(h*h) + c2*norm_u/h );

		//adjusting the stablization by a constant factor
		//double stab_factor = GetProperties()[STABILIZATION_FACTOR];
		//if(stab_factor != 0.0)
		//	tau *= GetProperties()[STABILIZATION_FACTOR];

		// *****************************************
		//CALCULATION OF THE LHS

		//CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
		noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
		noalias(rLeftHandSideMatrix) = outer_prod(msN,ms_u_DN)/(eps*eps);

		//CONVECTION STABILIZING CONTRIBUTION (Suu) 
		noalias(rLeftHandSideMatrix) += tau/(eps*eps) * outer_prod(ms_u_DN,ms_u_DN);
 
		//VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
		// rLeftHandSideMatrix += Laplacian * nu;
		noalias(rLeftHandSideMatrix) += nu/eps * prod(msDN_DX,trans(msDN_DX));

		//INERTIA CONTRIBUTION
		//  rLeftHandSideMatrix += M*BDFcoeffs[0]
		noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * msMassFactors/eps;

		//DARCY linear CONTRIBUTION
		//  rLeftHandSideMatrix -= nu/permeability (using the cinematic viscosity it is already divided for density: then we are going to multiplicate for density again);
		noalias(rLeftHandSideMatrix) += nu*msMassFactors*kinv;

// 		//DARCY non linear CONTRIBUTION (brinkmann)
// 		//rLeftHandSideMatrix -= 1.75*|u|/[(150*k)^0.5*eps^(3/2)]
// 		noalias(rLeftHandSideMatrix) += msMassFactors*norm_u*1.75*sqrt(kinv)/12.2474487/sqrt(eps*eps*eps)/density;

		//multiplication by the Volume
		rLeftHandSideMatrix *= (Volume * density);

		// *****************************************
		//CALCULATION OF THE RHS

		//external forces (component)
		double force_component = 0.25*(fcomp0 + fcomp1 + fcomp2 + fcomp3);

		//adding pressure gradient (integrated by parts)
		noalias(rRightHandSideVector) = (force_component )*msN;
		double p_avg = p0old + p1old + p2old + p3old;
		p_avg *= 0.25 / density;
		rRightHandSideVector[0] += msDN_DX(0,ComponentIndex)*p_avg; 
		rRightHandSideVector[1] += msDN_DX(1,ComponentIndex)*p_avg; 
		rRightHandSideVector[2] += msDN_DX(2,ComponentIndex)*p_avg;
		rRightHandSideVector[3] += msDN_DX(3,ComponentIndex)*p_avg;


		//adding the inertia terms
		// RHS += M*vhistory 
		//calculating the historical velocity
		noalias(ms_temp_vec_np) = ZeroVector(4);		
		for(unsigned int step = 1; step<BDFcoeffs.size(); step++)
		{ 
			for(int iii = 0; iii<4; iii++)
			{
				const array_1d<double,3>& v = (GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY,step) );
				ms_temp_vec_np[iii] -= BDFcoeffs[step]*v[ComponentIndex];
			}
			
		}	
		noalias(rRightHandSideVector) += prod(msMassFactors,ms_temp_vec_np)/eps ;

		//RHS += Suy * proj[component] 
		double proj_component = msN[0]*proj0[ComponentIndex] 
							  + msN[1]*proj1[ComponentIndex] 
							  + msN[2]*proj2[ComponentIndex]
							  + msN[3]*proj3[ComponentIndex];
		noalias(rRightHandSideVector) += (tau*proj_component)/(eps*eps)*ms_u_DN;   

		//multiplying by Volume
		rRightHandSideVector *= (Volume * density);

		//suubtracting the dirichlet term
		// RHS -= LHS*FracVel
		ms_temp_vec_np[0] = fv0[ComponentIndex]; 
		ms_temp_vec_np[1] = fv1[ComponentIndex]; 
		ms_temp_vec_np[2] = fv2[ComponentIndex]; 
		ms_temp_vec_np[3] = fv3[ComponentIndex]; 
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);


		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	//calculation by component of the fractional step velocity corresponding to the first stage
	void NDFluid3D::Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										   ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		unsigned int number_of_points = 4;

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points);

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points);

		//getting data for the given geometry
		double Volume;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Volume);
		//CalculateGeometryData(msDN_DX,msN,Volume);

		
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

		const array_1d<double,3>& fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj3 = GetGeometry()[3].FastGetSolutionStepValue(PRESS_PROJ);
		double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE); 
		double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE_OLD_IT); 
		const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
		const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);
		const double eps3 = GetGeometry()[3].FastGetSolutionStepValue(POROSITY);

		// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
		noalias(ms_aux) = fv0;		
		ms_aux[0] -= w0[0]; ms_aux[1] -= w0[1]; ms_aux[2] -= w0[2];
		noalias(ms_vel_gauss) = msN[0]*ms_aux;

		noalias(ms_aux) = fv1;	
		ms_aux[0] -= w1[0]; ms_aux[1] -= w1[1]; ms_aux[2] -= w1[2];
		noalias(ms_vel_gauss) += msN[1]*ms_aux;

		noalias(ms_aux) = fv2;		
		ms_aux[0] -= w2[0]; ms_aux[1] -= w2[1]; ms_aux[2] -= w2[2];
		noalias(ms_vel_gauss) += msN[2]*ms_aux;

		noalias(ms_aux) = fv3;	
		ms_aux[0] -= w3[0]; ms_aux[1] -= w3[1]; ms_aux[2] -= w3[2];
		noalias(ms_vel_gauss) += msN[3]*ms_aux;

		//calculating avergage density and viscosity
		double nu = 0.25*(nu0 + nu1 + nu2 + nu3);
 		double density = 0.25*(rho0 + rho1 + rho2 + rho3);

		//DIAMETER of the element
		//double dp = 0.33333333333333*(d0+d1+d2)
		//POROSITY of the element: average value of the porosity
		double eps = 0.25*(eps0 + eps1 + eps2 + eps3);
		
		//Getting delta time value for the restriction over tau.
		double Dt = rCurrentProcessInfo[DELTA_TIME];

		//calculating parameter tau (saved internally to each element)
		double c1 = 4.00;
		double c2 = 2.00;
		//double h = pow(6.00*Volume,0.3333333333);
		double h = CalculateH(Volume);
		double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1] + ms_vel_gauss[2]*ms_vel_gauss[2];
		norm_u = sqrt(norm_u);
		double tau = 1.00 / ( c1*nu/(h*h) + c2*norm_u/h );
		

		//tau = min{tau; Dt}
		if(tau > Dt)
		{
			tau = Dt;	
		}

		//getting the BDF2 coefficients (not fixed to allow variable time step)
		//the coefficients INCLUDE the time step
		const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

		//CALCULATION OF THE LEFT HAND SIDE
		//laplacian term	       L = Dt * gradN * trans(gradN);
		//stabilization term       Spp = tau * gradN * trans(gradN);
		noalias(rLeftHandSideMatrix) = ((1.00/BDFcoeffs[0] + tau)/density*eps) * prod(msDN_DX,trans(msDN_DX));

		//calculation of the RHS
		// RHS = -G*vfrac
		double Gaux;
		Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1] + msDN_DX(0,2)*fv0[2];
		Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1] + msDN_DX(1,2)*fv1[2];
		Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1] + msDN_DX(2,2)*fv2[2];
		Gaux += msDN_DX(3,0)*fv3[0] + msDN_DX(3,1)*fv3[1] + msDN_DX(3,2)*fv3[2];

		rRightHandSideVector[0] = - Gaux * msN[0]; 
		rRightHandSideVector[1] = - Gaux * msN[1]; 
		rRightHandSideVector[2] = - Gaux * msN[2]; 
		rRightHandSideVector[3] = - Gaux * msN[3]; 
		//std::cout << Id() << " Gtrans fv " << rRightHandSideVector << std::endl;

		//attention!! changing the meaning of ms_vel_gauss
		// RHS += Sz * proj*eps 
		//contrib of Spy*proj
		noalias(ms_vel_gauss) = msN[0]*proj0;
		noalias(ms_vel_gauss) += msN[1]*proj1;
		noalias(ms_vel_gauss) += msN[2]*proj2;
		noalias(ms_vel_gauss) += msN[3]*proj3;
		ms_vel_gauss *= tau * eps;
		noalias(rRightHandSideVector) += prod(msDN_DX , ms_vel_gauss);   
		//std::cout << Id() << " proj " << prod(msDN_DX , ms_vel_gauss) << std::endl;

		//dirichlet contribution
		ms_temp_vec_np[0] = p0; 
		ms_temp_vec_np[1] = p1; 
		ms_temp_vec_np[2] = p2; 
		ms_temp_vec_np[3] = p3; 
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);

		// RHS += dt * L *eps* pold  
		ms_temp_vec_np[0] = p0old; 
		ms_temp_vec_np[1] = p1old; 
		ms_temp_vec_np[2] = p2old; 
		ms_temp_vec_np[3] = p3old; 
		noalias(ms_vel_gauss) = prod(trans(msDN_DX),ms_temp_vec_np);
		noalias(rRightHandSideVector) += (1.00/(BDFcoeffs[0]*density)*eps ) * prod(msDN_DX,ms_vel_gauss); 

		//multiplicating by the Volume
		rLeftHandSideMatrix *= Volume;
		rRightHandSideVector *= Volume;

		//adding contributions to nodal Volumes following the corresponding lumping term
		double nodal_contrib = 0.25 * Volume * density;
		GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
		GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
		GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
		GetGeometry()[3].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;

		KRATOS_CATCH("");
	}
	  
	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void NDFluid3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
 
		//getting data for the given geometry
		double Volume;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Volume);
		//CalculateGeometryData(msDN_DX,msN,Volume);
		
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

			array_1d<double,3>& fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);
			array_1d<double,3>& w3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);
			array_1d<double,3>& press_proj3 = GetGeometry()[3].FastGetSolutionStepValue(PRESS_PROJ);
			array_1d<double,3>& conv_proj3 = GetGeometry()[3].FastGetSolutionStepValue(CONV_PROJ);
			double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
			const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

			double density = 0.25*(rho0 + rho1 + rho2 + rho3);

			//calculation of the pressure gradient (saved in ms_vel_gauss)
			ms_temp_vec_np[0] = p0;
			ms_temp_vec_np[1] = p1;
			ms_temp_vec_np[2] = p2;
			ms_temp_vec_np[3] = p3;
			noalias(ms_vel_gauss) = prod(trans(msDN_DX),ms_temp_vec_np);
// 			ms_vel_gauss *= Volume/density;
			ms_vel_gauss *= Volume;

			//press_proj += G*p
			noalias(press_proj0) += msN[0]*ms_vel_gauss;
			noalias(press_proj1) += msN[1]*ms_vel_gauss;
			noalias(press_proj2) += msN[2]*ms_vel_gauss;
			noalias(press_proj3) += msN[3]*ms_vel_gauss;

			// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
			noalias(ms_aux) = fv0;		
			ms_aux[0] -= w0[0]; ms_aux[1] -= w0[1]; ms_aux[2] -= w0[2];
			noalias(ms_vel_gauss) = msN[0]*ms_aux;

			noalias(ms_aux) = fv1;	
			ms_aux[0] -= w1[0]; ms_aux[1] -= w1[1]; ms_aux[2] -= w1[2];
			noalias(ms_vel_gauss) += msN[1]*ms_aux;

			noalias(ms_aux) = fv2;		
			ms_aux[0] -= w2[0]; ms_aux[1] -= w2[1]; ms_aux[2] -= w2[2];
			noalias(ms_vel_gauss) += msN[2]*ms_aux;

			noalias(ms_aux) = fv3;	
			ms_aux[0] -= w3[0]; ms_aux[1] -= w3[1]; ms_aux[2] -= w3[2]; 
			noalias(ms_vel_gauss) += msN[3]*ms_aux;


			//calculating convective auxiliary vector
			noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);

			//attention changing the meaning of ms_vel_gauss!!
			noalias(ms_vel_gauss) = ms_u_DN[0] * fv0;
			noalias(ms_vel_gauss) += ms_u_DN[1] * fv1;
			noalias(ms_vel_gauss) += ms_u_DN[2] * fv2;
			noalias(ms_vel_gauss) += ms_u_DN[3] * fv3;
// 			ms_vel_gauss *= Volume;
			ms_vel_gauss *= Volume * density;

			// conv_proj += C*u
			noalias(conv_proj0) += msN[0]*ms_vel_gauss;
			noalias(conv_proj1) += msN[1]*ms_vel_gauss;
			noalias(conv_proj2) += msN[2]*ms_vel_gauss;
			noalias(conv_proj3) += msN[3]*ms_vel_gauss;
 
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

			double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
			double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE_OLD_IT);
			array_1d<double,3>& fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);
// 			const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);
// 			const double eps3 = GetGeometry()[3].FastGetSolutionStepValue(POROSITY);

// 			double density = 0.25*(rho0 + rho1 + rho2 + rho3);
// 			double eps = 0.25*(eps0 + eps1 + eps2 + eps3);

 //G pn form
	
			double p_avg =	msN[0]*(p0 - p0old) 
					+ msN[1]*(p1 - p1old)
					+ msN[2]*(p2 - p2old) 
					+ msN[3]*(p3 - p3old);
// 			p_avg *= Volume/density;
			p_avg *= Volume;
			noalias(fv0)  += p_avg * row(msDN_DX,0);  
			noalias(fv1)  += p_avg * row(msDN_DX,1); 
			noalias(fv2)  += p_avg * row(msDN_DX,2); 
			noalias(fv3)  += p_avg * row(msDN_DX,3); 

		}


		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void NDFluid3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes);	

		int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber == 1) //step 1
			for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_X).EquationId();
		else if(FractionalStepNumber == 2) //step 2
			for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_Y).EquationId();
		else if(FractionalStepNumber == 3) //step 3
			for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_Z).EquationId();

		else if(FractionalStepNumber == 4) // pressure correction step
			for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
		}

	//************************************************************************************
	//************************************************************************************
	  void NDFluid3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber == 1) //step 1
			for (unsigned int i=0;i<number_of_nodes;i++)
				ElementalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_X);
		else if(FractionalStepNumber == 2) //step 2
			for (unsigned int i=0;i<number_of_nodes;i++)
				ElementalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_Y);
		else if(FractionalStepNumber == 3) //step 2
			for (unsigned int i=0;i<number_of_nodes;i++)
				ElementalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_Z);
		else if(FractionalStepNumber == 4) // pressure correction step
			for (unsigned int i=0;i<number_of_nodes;i++)
				ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);

	}

	//************************************************************************************
	//************************************************************************************
/*	inline void NDFluid3D::CalculateGeometryData(boost::numeric::ublas::bounded_matrix<double,4,3>& DN_DX, array_1d<double,4>& N, double& Volume)
	{
		double x0 = GetGeometry()[0].X();
		double x1 = GetGeometry()[1].X();
		double x2 = GetGeometry()[2].X();
		double x3 = GetGeometry()[3].X();

		double y0 = GetGeometry()[0].Y();
		double y1 = GetGeometry()[1].Y();
		double y2 = GetGeometry()[2].Y();
		double y3 = GetGeometry()[3].Y();

		double z0 = GetGeometry()[0].Z();
		double z1 = GetGeometry()[1].Z();
		double z2 = GetGeometry()[2].Z();
		double z3 = GetGeometry()[3].Z();

		//local derivatives
		msDN_De(0,0) = -1; msDN_De(0,1) = -1; msDN_De(0,2) = -1; //(0,0,0)
		msDN_De(1,0) = 1;  msDN_De(1,1) = 0;  msDN_De(1,2) = 0; //(1,0,0)
		msDN_De(2,0) = 0;  msDN_De(2,1) = 1;  msDN_De(2,2) = 0; //(0,1,0)
		msDN_De(3,0) = 0;  msDN_De(3,1) = 0;  msDN_De(3,2) = 1; //(0,0,1)
 
		//calculation of the jacobian
		msJ(0,0) = x1-x0; msJ(0,1) = x2-x0; msJ(0,2) = x3-x0;
		msJ(1,0) = y1-y0; msJ(1,1) = y2-y0; msJ(1,2) = y3-y0;
		msJ(2,0) = z1-z0; msJ(2,1) = z2-z0; msJ(2,2) = z3-z0;

		//inverse of the jacobian
		//first column
		msJinv(0,0) = msJ(1,1)*msJ(2,2) - msJ(1,2)*msJ(2,1);
		msJinv(1,0) = -msJ(1,0)*msJ(2,2) + msJ(1,2)*msJ(2,0);
		msJinv(2,0) = msJ(1,0)*msJ(2,1) - msJ(1,1)*msJ(2,0);		
		//second column
		msJinv(0,1) = -msJ(0,1)*msJ(2,2) + msJ(0,2)*msJ(2,1);
		msJinv(1,1) = msJ(0,0)*msJ(2,2) - msJ(0,2)*msJ(2,0);
		msJinv(2,1) = -msJ(0,0)*msJ(2,1) + msJ(0,1)*msJ(2,0);
		//third column
		msJinv(0,2) = msJ(0,1)*msJ(1,2) - msJ(0,2)*msJ(1,1);
		msJinv(1,2) = -msJ(0,0)*msJ(1,2) + msJ(0,2)*msJ(1,0);
		msJinv(2,2) = msJ(0,0)*msJ(1,1) - msJ(0,1)*msJ(1,0);
		//calculation of determinant (of the input matrix)
 
		double detJ = msJ(0,0)*msJinv(0,0) + msJ(0,1)*msJinv(1,0) + msJ(0,2)*msJinv(2,0);	
		//finalizing the calculation of the inverted matrix
		msJinv /= detJ;

		Volume = detJ*0.1666666666666666666667;

		//cartesian derivatives
		noalias(msDN_DX) = prod(msDN_De,msJinv);

		//shape function values
		N[0] = 0.25; //(0,0,0)
		N[1] = 0.25;  //(1,0,0)
		N[2] = 0.25;   //(0,1,0)
		N[3] = 0.25;  //(0,0,1)


	}
*/
	
	inline double NDFluid3D::CalculateH(double Volume)
	{
		double h = pow(6.00*Volume,0.3333333);

		//double h = 0.0;
		//for(int j=0; j<3; j++)
		//{
		//	for(int i = j+1; i<4; i++)
		//	{ 
		//		double l;
		//		l =  pow(GetGeometry()[i].X() - GetGeometry()[j].X()   ,2);
		//		l += pow(GetGeometry()[i].Y() - GetGeometry()[j].Y()   ,2);
		//		l += pow(GetGeometry()[i].Z() - GetGeometry()[j].Z()   ,2);

		//		if(l > h) h = l;
		//	}
		//}
		//h = sqrt(h);
		
		return h;
	}

} // Namespace Kratos


