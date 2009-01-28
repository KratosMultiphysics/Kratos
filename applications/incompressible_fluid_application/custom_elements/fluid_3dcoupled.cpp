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
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.5 $
//
//


// System includes 
//#define GRADPN_FORM //the grad(pn) is used instead of the G(pn) in doing the splitting
 
// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/fluid_3dcoupled.h"
#include "incompressible_fluid_application.h"
#include "utilities/math_utils.h" 
#include "utilities/geometry_utilities.h" 

namespace Kratos 
{
	//static variables
	boost::numeric::ublas::bounded_matrix<double,4,4> Fluid3DCoupled::msaux_matrix;
	boost::numeric::ublas::bounded_matrix<double,4,4> Fluid3DCoupled::msMassFactors;
	boost::numeric::ublas::bounded_matrix<double,4,3> Fluid3DCoupled::msDN_DX;
  	array_1d<double,4> Fluid3DCoupled::msN; //dimension = number of nodes
	array_1d<double,3> Fluid3DCoupled::ms_aux; //dimension coincides with space dimension
	array_1d<double,3> Fluid3DCoupled::ms_vel_gauss; //dimesion coincides with space dimension
  	array_1d<double,4> Fluid3DCoupled::ms_temp_vec_np; //dimension = number of nodes
  	array_1d<double,4> Fluid3DCoupled::ms_u_DN; //dimension = number of nodes
	
	
		


	//************************************************************************************
	//************************************************************************************
	Fluid3DCoupled::Fluid3DCoupled(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Fluid3DCoupled::Fluid3DCoupled(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
		InitializeAuxiliaries();
	}

	Element::Pointer Fluid3DCoupled::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		KRATOS_TRY
		return Element::Pointer(new Fluid3DCoupled(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	Fluid3DCoupled::~Fluid3DCoupled()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void Fluid3DCoupled::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
	void Fluid3DCoupled::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}

	//************************************************************************************
	//************************************************************************************
	//calculation by component of the fractional step velocity corresponding to the first stage
	void Fluid3DCoupled::Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										   ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		const unsigned int number_of_points = 4;
		const unsigned int dim = 3;
		unsigned int matsize = number_of_points *dim;
		
		if(rLeftHandSideMatrix.size1() != matsize)
			rLeftHandSideMatrix.resize(matsize,matsize,false);

		if(rRightHandSideVector.size() != matsize)
			rRightHandSideVector.resize(matsize,false);
		
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

		const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj1 = GetGeometry()[1].FastGetSolutionStepValue(CONV_PROJ);
		const double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

		const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj2 = GetGeometry()[2].FastGetSolutionStepValue(CONV_PROJ);
		const double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

		const array_1d<double,3>& fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj3 = GetGeometry()[3].FastGetSolutionStepValue(CONV_PROJ);
		const double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
		const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

		// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
		noalias(ms_aux) =  fv0;	noalias(ms_aux) -= w0;	
		noalias(ms_aux) += fv1;	noalias(ms_aux) -= w1;	
		noalias(ms_aux) += fv2;	noalias(ms_aux) -= w2;	
		noalias(ms_aux) += fv3;	noalias(ms_aux) -= w3;	
		noalias(ms_vel_gauss) = 0.25*ms_aux;

		//calculating viscosity
		double nu = 0.25*(nu0 + nu1 + nu2 + nu3);
 		double density = 0.25*(rho0 + rho1 + rho2 + rho3);

		//getting the BDF2 coefficients (not fixed to allow variable time step)
		//the coefficients INCLUDE the time step
		const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

		//calculating parameter tau (saved internally to each element)
		double c1 = 4.00;
		double c2 = 2.00;
		double h = CalculateH(Volume);
		double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1] + ms_vel_gauss[2]*ms_vel_gauss[2];
		norm_u = sqrt(norm_u);
		double tau = 1.00 / ( c1*nu/(h*h) + c2*norm_u/h );

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
		rLeftHandSideMatrix *= Volume;
		

		// *****************************************
		//CALCULATION OF THE RHS
		const array_1d<double,3>& force0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
		const array_1d<double,3>& force1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
		const array_1d<double,3>& force2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);		
		const array_1d<double,3>& force3 = GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE);		

		array_1d<double,4> rhs_aux;
		for( unsigned int component_index = 0; component_index < dim; component_index++)
		{
			//external forces (component)
			double force_component = 0.25*(force0[component_index] +
										force1[component_index] + 
										force2[component_index] + 
										force3[component_index]);
			noalias(rhs_aux) = (force_component )*msN;

			//adding pressure gradient (integrated by parts)
			double p_avg = p0old + p1old + p2old + p3old;
			p_avg *= 0.25 / density;
			rhs_aux[0] += msDN_DX(0,component_index)*p_avg; 
			rhs_aux[1] += msDN_DX(1,component_index)*p_avg; 
			rhs_aux[2] += msDN_DX(2,component_index)*p_avg;
			rhs_aux[3] += msDN_DX(3,component_index)*p_avg;

			//adding the inertia terms
			// RHS += M*vhistory 
			//calculating the historical velocity

			noalias(ms_temp_vec_np) = ZeroVector(4);
			for(unsigned int step = 1; step<BDFcoeffs.size(); step++)
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
								  + proj2[component_index]
								  + proj3[component_index];
			proj_component *= 0.25;
			noalias(rhs_aux) += (tau*proj_component)*ms_u_DN;   

			//writing the rhs_aux in its place
			for( unsigned int i = 0; i < number_of_points; i++)
			{
				rRightHandSideVector[i*dim + component_index] = rhs_aux[i];
			}
		}
		
		//multiplying by area
		rRightHandSideVector *= Volume;

		//LHS dirichlet contribution
		Vector fvvect(12);
		for( unsigned int component_index = 0; component_index < dim; component_index++)
		{
				fvvect[0 + component_index] = fv0[component_index];
				fvvect[3 + component_index] = fv1[component_index];
				fvvect[6 + component_index] = fv2[component_index];
				fvvect[9 + component_index] = fv3[component_index];
		}
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,fvvect);
//KRATOS_WATCH(rRightHandSideVector)
//KRATOS_WATCH(rLeftHandSideMatrix)

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	//calculation by component of the fractional step velocity corresponding to the first stage
	void Fluid3DCoupled::Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
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

		const array_1d<double,3>& fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& w3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);
		const array_1d<double,3>& proj3 = GetGeometry()[3].FastGetSolutionStepValue(PRESS_PROJ);
		double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE); 
		double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE_OLD_IT); 
		const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
		const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

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

		//calculating parameter tau (saved internally to each element)
		double c1 = 4.00;
		double c2 = 2.00;
		//double h = pow(6.00*Volume,0.3333333333);
		double h = CalculateH(Volume);
		double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1] + ms_vel_gauss[2]*ms_vel_gauss[2];
		norm_u = sqrt(norm_u);
		double tau = 1.00 / ( c1*nu/(h*h) + c2*norm_u/h );
		
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
		// RHS += Sz * proj 
		//contrib of Spy*proj
		noalias(ms_vel_gauss) = msN[0]*proj0;
		noalias(ms_vel_gauss) += msN[1]*proj1;
		noalias(ms_vel_gauss) += msN[2]*proj2;
		noalias(ms_vel_gauss) += msN[3]*proj3;
		ms_vel_gauss *= tau;
		noalias(rRightHandSideVector) += prod(msDN_DX , ms_vel_gauss);   
		//std::cout << Id() << " proj " << prod(msDN_DX , ms_vel_gauss) << std::endl;

		//dirichlet contribution
		ms_temp_vec_np[0] = p0; 
		ms_temp_vec_np[1] = p1; 
		ms_temp_vec_np[2] = p2; 
		ms_temp_vec_np[3] = p3; 
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);

		// RHS += dt * L * pold 
		ms_temp_vec_np[0] = p0old; 
		ms_temp_vec_np[1] = p1old; 
		ms_temp_vec_np[2] = p2old; 
		ms_temp_vec_np[3] = p3old; 
		noalias(ms_vel_gauss) = prod(trans(msDN_DX),ms_temp_vec_np);
		noalias(rRightHandSideVector) += (1.00/(BDFcoeffs[0]*density) ) * prod(msDN_DX,ms_vel_gauss); 

		//multiplicating by the Volume
		rLeftHandSideMatrix *= Volume;
		rRightHandSideVector *= Volume;

		//adding contributions to nodal Volumes following the corresponding lumping term
		double nodal_contrib = 0.25 * Volume;
		GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib;
		GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib;
		GetGeometry()[2].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib;
		GetGeometry()[3].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib;

		KRATOS_CATCH("");
	}
	  
	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void Fluid3DCoupled::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
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
			ms_vel_gauss *= Volume/density;

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
			ms_vel_gauss *= Volume;

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
			const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

			double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
			double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
			array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
			const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

			double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
			double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
			array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
			const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

			double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
			double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE_OLD_IT);
			array_1d<double,3>& fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);
			const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

			double density = 0.25*(rho0 + rho1 + rho2 + rho3);

			double p_avg =	msN[0]*(p0 - p0old) 
					+ msN[1]*(p1 - p1old)
					+ msN[2]*(p2 - p2old) 
					+ msN[3]*(p3 - p3old);
			p_avg *= Volume/density;
			noalias(fv0)  += p_avg * row(msDN_DX,0);  
			noalias(fv1)  += p_avg * row(msDN_DX,1); 
			noalias(fv2)  += p_avg * row(msDN_DX,2); 
			noalias(fv3)  += p_avg * row(msDN_DX,3); 
		}


		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void Fluid3DCoupled::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		const unsigned int dim = 3;

		unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber == 1) //step 1
		{
			if(rResult.size() != number_of_nodes*dim)
				rResult.resize(number_of_nodes*dim,false);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				rResult[i*dim] = GetGeometry()[i].GetDof(FRACT_VEL_X).EquationId();
				rResult[i*dim+1] = GetGeometry()[i].GetDof(FRACT_VEL_Y).EquationId();
				rResult[i*dim+2] = GetGeometry()[i].GetDof(FRACT_VEL_Z).EquationId();
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
	  void Fluid3DCoupled::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 3;

		unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber == 1) //step 1
		{
			if(ElementalDofList.size() != number_of_nodes*dim)
				ElementalDofList.resize(number_of_nodes*dim);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				ElementalDofList[i*dim] = GetGeometry()[i].pGetDof(FRACT_VEL_X);
				ElementalDofList[i*dim+1] = GetGeometry()[i].pGetDof(FRACT_VEL_Y);
				ElementalDofList[i*dim+2] = GetGeometry()[i].pGetDof(FRACT_VEL_Z);
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
	void Fluid3DCoupled::CalculateViscousMatrix(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double,4,3>& DN_DX, const double& nu)
	{
		//horrorful ... just to make sure it works
		Matrix B(6,12);
		noalias(B) = ZeroMatrix(6,12);
		unsigned int start;
		for(unsigned int i = 0; i<4; i++)
		{
			start=i*3;
			B(0,start) = DN_DX(i,0);
			B(1,start+1) = DN_DX(i,1);
			B(2,start+2) = DN_DX(i,2);
			B(3,start) = DN_DX(i,1);	B(3,start+1) = DN_DX(i,0);
			B(4,start) = DN_DX(i,2);	B(4,start+2) = DN_DX(i,0);
			B(5,start+1) = DN_DX(i,2);	B(5,start+2) = DN_DX(i,1);
		}

		Matrix msCapx(6,6);
		noalias(msCapx) = ZeroMatrix(6,6);
		msCapx(0,0)=2.0*nu; msCapx(1,1)=2.0*nu; msCapx(2,2)=2.0*nu;
		msCapx(3,3)=nu;     msCapx(4,4)=nu;     msCapx(5,5)=nu;

		Matrix temp(6,12);
		noalias(temp) = prod(msCapx,B);
		
		noalias(K) = prod(trans(B),temp);


/*		K(0,1) = DN_DX(0,1) * nu * DN_DX(0,0);
		K(0,2) = DN_DX(0,2) * nu * DN_DX(0,0);
		K(0,3) = 2.00 * DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,1) * nu * DN_DX(1,1) + DN_DX(0,2) * nu * DN_DX(1,2);
		K(0,4) = DN_DX(0,1) * nu * DN_DX(1,0);
		K(0,5) = DN_DX(0,2) * nu * DN_DX(1,0);
		K(0,6) = 2.00 * DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,1) * nu * DN_DX(2,1) + DN_DX(0,2) * nu * DN_DX(2,2);
		K(0,7) = DN_DX(0,1) * nu * DN_DX(2,0);
		K(0,8) = DN_DX(0,2) * nu * DN_DX(2,0);
		K(0,9) = 2.00 * DN_DX(0,0) * nu * DN_DX(3,0) + DN_DX(0,1) * nu * DN_DX(3,1) + DN_DX(0,2) * nu * DN_DX(3,2);
		K(0,10) = DN_DX(0,1) * nu * DN_DX(3,0);
		K(0,11) = DN_DX(0,2) * nu * DN_DX(3,0);
		K(1,0) = DN_DX(0,1) * nu * DN_DX(0,0);
		K(1,1) = 2.00 * nu * pow(DN_DX(0,1), 2) + pow(DN_DX(0,0), 2) * nu + nu * pow(DN_DX(0,2), 2);
		K(1,2) = DN_DX(0,2) * nu * DN_DX(0,1);
		K(1,3) = DN_DX(0,0) * nu * DN_DX(1,1);
		K(1,4) = 2.00 * DN_DX(0,1) * nu * DN_DX(1,1) + DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,2) * nu * DN_DX(1,2);
		K(1,5) = DN_DX(0,2) * nu * DN_DX(1,1);
		K(1,6) = DN_DX(0,0) * nu * DN_DX(2,1);
		K(1,7) = 2.00 * DN_DX(0,1) * nu * DN_DX(2,1) + DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,2) * nu * DN_DX(2,2);
		K(1,8) = DN_DX(0,2) * nu * DN_DX(2,1);
		K(1,9) = DN_DX(0,0) * nu * DN_DX(3,1);
		K(1,10) = 2.00 * DN_DX(0,1) * nu * DN_DX(3,1) + DN_DX(0,0) * nu * DN_DX(3,0) + DN_DX(0,2) * nu * DN_DX(3,2);
		K(1,11) = DN_DX(0,2) * nu * DN_DX(3,1);
		K(2,0) = DN_DX(0,2) * nu * DN_DX(0,0);
		K(2,1) = DN_DX(0,2) * nu * DN_DX(0,1);
		K(2,2) = 2.00 * nu * pow(DN_DX(0,2), 2) + pow(DN_DX(0,0), 2) * nu + nu * pow(DN_DX(0,1), 2);
		K(2,3) = DN_DX(0,0) * nu * DN_DX(1,2);
		K(2,4) = DN_DX(0,1) * nu * DN_DX(1,2);
		K(2,5) = 2.00 * DN_DX(0,2) * nu * DN_DX(1,2) + DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,1) * nu * DN_DX(1,1);
		K(2,6) = DN_DX(0,0) * nu * DN_DX(2,2);
		K(2,7) = DN_DX(0,1) * nu * DN_DX(2,2);
		K(2,8) = 2.00 * DN_DX(0,2) * nu * DN_DX(2,2) + DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,1) * nu * DN_DX(2,1);
		K(2,9) = DN_DX(0,0) * nu * DN_DX(3,2);
		K(2,10) = DN_DX(0,1) * nu * DN_DX(3,2);
		K(2,11) = 2.00 * DN_DX(0,2) * nu * DN_DX(3,2) + DN_DX(0,0) * nu * DN_DX(3,0) + DN_DX(0,1) * nu * DN_DX(3,1);
		K(3,0) = 2.00 * DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,1) * nu * DN_DX(1,1) + DN_DX(0,2) * nu * DN_DX(1,2);
		K(3,1) = DN_DX(0,0) * nu * DN_DX(1,1);
		K(3,2) = DN_DX(0,0) * nu * DN_DX(1,2);
		K(3,3) = 2.00 * pow(DN_DX(1,0), 2) * nu + nu * pow(DN_DX(1,1), 2) + nu * pow(DN_DX(1,2), 2);
		K(3,4) = DN_DX(1,1) * nu * DN_DX(1,0);
		K(3,5) = DN_DX(1,2) * nu * DN_DX(1,0);
		K(3,6) = 2.00 * DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,1) * nu * DN_DX(2,1) + DN_DX(1,2) * nu * DN_DX(2,2);
		K(3,7) = DN_DX(1,1) * nu * DN_DX(2,0);
		K(3,8) = DN_DX(1,2) * nu * DN_DX(2,0);
		K(3,9) = 2.00 * DN_DX(1,0) * nu * DN_DX(3,0) + DN_DX(1,1) * nu * DN_DX(3,1) + DN_DX(1,2) * nu * DN_DX(3,2);
		K(3,10) = DN_DX(1,1) * nu * DN_DX(3,0);
		K(3,11) = DN_DX(1,2) * nu * DN_DX(3,0);
		K(4,0) = DN_DX(0,1) * nu * DN_DX(1,0);
		K(4,1) = 2.00 * DN_DX(0,1) * nu * DN_DX(1,1) + DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,2) * nu * DN_DX(1,2);
		K(4,2) = DN_DX(0,1) * nu * DN_DX(1,2);
		K(4,3) = DN_DX(1,1) * nu * DN_DX(1,0);
		K(4,4) = 2.00 * nu * pow(DN_DX(1,1), 2) + pow(DN_DX(1,0), 2) * nu + nu * pow(DN_DX(1,2), 2);
		K(4,5) = DN_DX(1,2) * nu * DN_DX(1,1);
		K(4,6) = DN_DX(1,0) * nu * DN_DX(2,1);
		K(4,7) = 2.00 * DN_DX(1,1) * nu * DN_DX(2,1) + DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,2) * nu * DN_DX(2,2);
		K(4,8) = DN_DX(1,2) * nu * DN_DX(2,1);
		K(4,9) = DN_DX(1,0) * nu * DN_DX(3,1);
		K(4,10) = 2.00 * DN_DX(1,1) * nu * DN_DX(3,1) + DN_DX(1,0) * nu * DN_DX(3,0) + DN_DX(1,2) * nu * DN_DX(3,2);
		K(4,11) = DN_DX(1,2) * nu * DN_DX(3,1);
		K(5,0) = DN_DX(0,2) * nu * DN_DX(1,0);
		K(5,1) = DN_DX(0,2) * nu * DN_DX(1,1);
		K(5,2) = 2.00 * DN_DX(0,2) * nu * DN_DX(1,2) + DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,1) * nu * DN_DX(1,1);
		K(5,3) = DN_DX(1,2) * nu * DN_DX(1,0);
		K(5,4) = DN_DX(1,2) * nu * DN_DX(1,1);
		K(5,5) = 2.00 * nu * pow(DN_DX(1,2), 2) + pow(DN_DX(1,0), 2) * nu + nu * pow(DN_DX(1,1), 2);
		K(5,6) = DN_DX(1,0) * nu * DN_DX(2,2);
		K(5,7) = DN_DX(1,1) * nu * DN_DX(2,2);
		K(5,8) = 2.00 * DN_DX(1,2) * nu * DN_DX(2,2) + DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,1) * nu * DN_DX(2,1);
		K(5,9) = DN_DX(1,0) * nu * DN_DX(3,2);
		K(5,10) = DN_DX(1,1) * nu * DN_DX(3,2);
		K(5,11) = 2.00 * DN_DX(1,2) * nu * DN_DX(3,2) + DN_DX(1,0) * nu * DN_DX(3,0) + DN_DX(1,1) * nu * DN_DX(3,1);
		K(6,0) = 2.00 * DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,1) * nu * DN_DX(2,1) + DN_DX(0,2) * nu * DN_DX(2,2);
		K(6,1) = DN_DX(0,0) * nu * DN_DX(2,1);
		K(6,2) = DN_DX(0,0) * nu * DN_DX(2,2);
		K(6,3) = 2.00 * DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,1) * nu * DN_DX(2,1) + DN_DX(1,2) * nu * DN_DX(2,2);
		K(6,4) = DN_DX(1,0) * nu * DN_DX(2,1);
		K(6,5) = DN_DX(1,0) * nu * DN_DX(2,2);
		K(6,6) = 2.00 * pow(DN_DX(2,0), 2) * nu + nu * pow(DN_DX(2,1), 2) + nu * pow(DN_DX(2,2), 2);
		K(6,7) = DN_DX(2,1) * nu * DN_DX(2,0);
		K(6,8) = DN_DX(2,2) * nu * DN_DX(2,0);
		K(6,9) = 2.00 * DN_DX(2,0) * nu * DN_DX(3,0) + DN_DX(2,1) * nu * DN_DX(3,1) + DN_DX(2,2) * nu * DN_DX(3,2);
		K(6,10) = DN_DX(2,1) * nu * DN_DX(3,0);
		K(6,11) = DN_DX(2,2) * nu * DN_DX(3,0);
		K(7,0) = DN_DX(0,1) * nu * DN_DX(2,0);
		K(7,1) = 2.00 * DN_DX(0,1) * nu * DN_DX(2,1) + DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,2) * nu * DN_DX(2,2);
		K(7,2) = DN_DX(0,1) * nu * DN_DX(2,2);
		K(7,3) = DN_DX(1,1) * nu * DN_DX(2,0);
		K(7,4) = 2.00 * DN_DX(1,1) * nu * DN_DX(2,1) + DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,2) * nu * DN_DX(2,2);
		K(7,5) = DN_DX(1,1) * nu * DN_DX(2,2);
		K(7,6) = DN_DX(2,1) * nu * DN_DX(2,0);
		K(7,7) = 2.00 * nu * pow(DN_DX(2,1), 2) + pow(DN_DX(2,0), 2) * nu + nu * pow(DN_DX(2,2), 2);
		K(7,8) = DN_DX(2,2) * nu * DN_DX(2,1);
		K(7,9) = DN_DX(2,0) * nu * DN_DX(3,1);
		K(7,10) = 2.00 * DN_DX(2,1) * nu * DN_DX(3,1) + DN_DX(2,0) * nu * DN_DX(3,0) + DN_DX(2,2) * nu * DN_DX(3,2);
		K(7,11) = DN_DX(2,2) * nu * DN_DX(3,1);
		K(8,0) = DN_DX(0,2) * nu * DN_DX(2,0);
		K(8,1) = DN_DX(0,2) * nu * DN_DX(2,1);
		K(8,2) = 2.00 * DN_DX(0,2) * nu * DN_DX(2,2) + DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,1) * nu * DN_DX(2,1);
		K(8,3) = DN_DX(1,2) * nu * DN_DX(2,0);
		K(8,4) = DN_DX(1,2) * nu * DN_DX(2,1);
		K(8,5) = 2.00 * DN_DX(1,2) * nu * DN_DX(2,2) + DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,1) * nu * DN_DX(2,1);
		K(8,6) = DN_DX(2,2) * nu * DN_DX(2,0);
		K(8,7) = DN_DX(2,2) * nu * DN_DX(2,1);
		K(8,8) = 2.00 * nu * pow(DN_DX(2,2), 2) + pow(DN_DX(2,0), 2) * nu + nu * pow(DN_DX(2,1), 2);
		K(8,9) = DN_DX(2,0) * nu * DN_DX(3,2);
		K(8,10) = DN_DX(2,1) * nu * DN_DX(3,2);
		K(8,11) = 2.00 * DN_DX(2,2) * nu * DN_DX(3,2) + DN_DX(2,0) * nu * DN_DX(3,0) + DN_DX(2,1) * nu * DN_DX(3,1);
		K(9,0) = 2.00 * DN_DX(0,0) * nu * DN_DX(3,0) + DN_DX(0,1) * nu * DN_DX(3,1) + DN_DX(0,2) * nu * DN_DX(3,2);
		K(9,1) = DN_DX(0,0) * nu * DN_DX(3,1);
		K(9,2) = DN_DX(0,0) * nu * DN_DX(3,2);
		K(9,3) = 2.00 * DN_DX(1,0) * nu * DN_DX(3,0) + DN_DX(1,1) * nu * DN_DX(3,1) + DN_DX(1,2) * nu * DN_DX(3,2);
		K(9,4) = DN_DX(1,0) * nu * DN_DX(3,1);
		K(9,5) = DN_DX(1,0) * nu * DN_DX(3,2);
		K(9,6) = 2.00 * DN_DX(2,0) * nu * DN_DX(3,0) + DN_DX(2,1) * nu * DN_DX(3,1) + DN_DX(2,2) * nu * DN_DX(3,2);
		K(9,7) = DN_DX(2,0) * nu * DN_DX(3,1);
		K(9,8) = DN_DX(2,0) * nu * DN_DX(3,2);
		K(9,9) = 2.00 * pow(DN_DX(3,0), 2) * nu + nu * pow(DN_DX(3,1), 2) + nu * pow(DN_DX(3,2), 2);
		K(9,10) = DN_DX(3,1) * nu * DN_DX(3,0);
		K(9,11) = DN_DX(3,2) * nu * DN_DX(3,0);
		K(10,0) = DN_DX(0,1) * nu * DN_DX(3,0);
		K(10,1) = 2.00 * DN_DX(0,1) * nu * DN_DX(3,1) + DN_DX(0,0) * nu * DN_DX(3,0) + DN_DX(0,2) * nu * DN_DX(3,2);
		K(10,2) = DN_DX(0,1) * nu * DN_DX(3,2);
		K(10,3) = DN_DX(1,1) * nu * DN_DX(3,0);
		K(10,4) = 2.00 * DN_DX(1,1) * nu * DN_DX(3,1) + DN_DX(1,0) * nu * DN_DX(3,0) + DN_DX(1,2) * nu * DN_DX(3,2);
		K(10,5) = DN_DX(1,1) * nu * DN_DX(3,2);
		K(10,6) = DN_DX(2,1) * nu * DN_DX(3,0);
		K(10,7) = 2.00 * DN_DX(2,1) * nu * DN_DX(3,1) + DN_DX(2,0) * nu * DN_DX(3,0) + DN_DX(2,2) * nu * DN_DX(3,2);
		K(10,8) = DN_DX(2,1) * nu * DN_DX(3,2);
		K(10,9) = DN_DX(3,1) * nu * DN_DX(3,0);
		K(10,10) = 2.00 * nu * pow(DN_DX(3,1), 2) + pow(DN_DX(3,0), 2) * nu + nu * pow(DN_DX(3,2), 2);
		K(10,11) = DN_DX(3,2) * nu * DN_DX(3,1);
		K(11,0) = DN_DX(0,2) * nu * DN_DX(3,0);
		K(11,1) = DN_DX(0,2) * nu * DN_DX(3,1);
		K(11,2) = 2.00 * DN_DX(0,2) * nu * DN_DX(3,2) + DN_DX(0,0) * nu * DN_DX(3,0) + DN_DX(0,1) * nu * DN_DX(3,1);
		K(11,3) = DN_DX(1,2) * nu * DN_DX(3,0);
		K(11,4) = DN_DX(1,2) * nu * DN_DX(3,1);
		K(11,5) = 2.00 * DN_DX(1,2) * nu * DN_DX(3,2) + DN_DX(1,0) * nu * DN_DX(3,0) + DN_DX(1,1) * nu * DN_DX(3,1);
		K(11,6) = DN_DX(2,2) * nu * DN_DX(3,0);
		K(11,7) = DN_DX(2,2) * nu * DN_DX(3,1);
		K(11,8) = 2.00 * DN_DX(2,2) * nu * DN_DX(3,2) + DN_DX(2,0) * nu * DN_DX(3,0) + DN_DX(2,1) * nu * DN_DX(3,1);
		K(11,9) = DN_DX(3,2) * nu * DN_DX(3,0);
		K(11,10) = DN_DX(3,2) * nu * DN_DX(3,1);
		K(11,11) = 2.00 * nu * pow(DN_DX(3,2), 2) + pow(DN_DX(3,0), 2) * nu + nu * pow(DN_DX(3,1), 2);

		//filling the symmetric part
		for(unsigned int i = 1; i<K.size1(); i++)
			for(unsigned int j = 0; j<i; j++)
				K(i,j) = K(j,i);
	*/
	}

		//***********************************************************************
		//***********************************************************************
		//performs the Kroneker product of the Reduced Matrix with the identity matrix of 
		//size "dimension" ADDING to the destination matrix
		inline void  Fluid3DCoupled::ExpandAndAddReducedMatrix(
			MatrixType& Destination,
			boost::numeric::ublas::bounded_matrix<double,4,4>& ReducedMatrix,
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

		//***********************************************************************
		//***********************************************************************

	inline double Fluid3DCoupled::CalculateH(double Volume)
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


