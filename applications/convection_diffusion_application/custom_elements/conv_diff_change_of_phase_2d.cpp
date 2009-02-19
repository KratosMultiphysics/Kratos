/*
==============================================================================
KratosConvectionDiffusionApplication 
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
//   Date:                $Date: 2009-01-13 15:39:55 $
//   Revision:            $Revision: 1.2 $
//
//
 

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/conv_diff_change_of_phase_2d.h"
#include "convection_diffusion_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
	//static variables
	boost::numeric::ublas::bounded_matrix<double,3,3> ConvDiffChangeOfPhase2D::msMassFactors;
	boost::numeric::ublas::bounded_matrix<double,3,2> ConvDiffChangeOfPhase2D::msDN_DX;
  	array_1d<double,3> ConvDiffChangeOfPhase2D::msN; //dimension = number of nodes
	array_1d<double,2> ConvDiffChangeOfPhase2D::ms_vel_gauss; //dimesion coincides with space dimension
  	array_1d<double,3> ConvDiffChangeOfPhase2D::ms_temp_vec_np; //dimension = number of nodes
  	array_1d<double,3> ConvDiffChangeOfPhase2D::ms_u_DN; //dimension = number of nodes

	//************************************************************************************
	//************************************************************************************
	ConvDiffChangeOfPhase2D::ConvDiffChangeOfPhase2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	ConvDiffChangeOfPhase2D::ConvDiffChangeOfPhase2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
		//filling the mass factors
//	msMassFactors(0,0) = 1.00/6.00;  msMassFactors(0,1) = 1.00/12.00; msMassFactors(0,2) = 1.00/12.00;
//		msMassFactors(1,0) = 1.00/12.00; msMassFactors(1,1) = 1.00/6.00;  msMassFactors(1,2) = 1.00/12.00;
//		msMassFactors(2,0) = 1.00/12.00; msMassFactors(2,1) = 1.00/12.00; msMassFactors(2,2) = 1.00/6.00;
//
		msMassFactors(0,0) = 1.00/3.00; msMassFactors(0,1) = 0.00;		msMassFactors(0,2) = 0.00;
		msMassFactors(1,0) = 0.00;		msMassFactors(1,1) = 1.00/3.00; msMassFactors(1,2) = 0.00;
		msMassFactors(2,0) = 0.00;		msMassFactors(2,1) = 0.00;		msMassFactors(2,2) = 1.00/3.00;	
		
	}

	Element::Pointer ConvDiffChangeOfPhase2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new ConvDiffChangeOfPhase2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	ConvDiffChangeOfPhase2D::~ConvDiffChangeOfPhase2D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void ConvDiffChangeOfPhase2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int number_of_points = GetGeometry().size();
		const double lumping_factor = 1.00/double(number_of_points);
unsigned int TDim = 2;

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);

		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);

		//calculating viscosity
		double conductivity = GetGeometry()[0].FastGetSolutionStepValue(CONDUCTIVITY);
		double density = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		double specific_heat = GetGeometry()[0].FastGetSolutionStepValue(SPECIFIC_HEAT);
		double heat_flux = GetGeometry()[0].FastGetSolutionStepValue(HEAT_FLUX);
		double proj = GetGeometry()[0].FastGetSolutionStepValue(TEMP_CONV_PROJ);
		
		const array_1d<double,3>& v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& w = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
		for(unsigned int j = 0; j<TDim; j++)
			ms_vel_gauss[j] = v[j] - w[j];
		
		for(unsigned int i = 1; i<number_of_points; i++)
		{
			conductivity += GetGeometry()[i].FastGetSolutionStepValue(CONDUCTIVITY);
			density += GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
			specific_heat += GetGeometry()[i].FastGetSolutionStepValue(SPECIFIC_HEAT);
			heat_flux += GetGeometry()[i].FastGetSolutionStepValue(HEAT_FLUX);
			proj += GetGeometry()[i].FastGetSolutionStepValue(TEMP_CONV_PROJ);
			
			const array_1d<double,3>& v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
			const array_1d<double,3>& w = GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);
			for(unsigned int j = 0; j<TDim; j++)
				ms_vel_gauss[j] += v[j] - w[j];
			
		}
		conductivity *= lumping_factor;
		density *= lumping_factor;
		specific_heat *= lumping_factor;
		heat_flux *= lumping_factor;
		proj *= lumping_factor;
		ms_vel_gauss *= lumping_factor;
		
		double alpha = conductivity/(density*specific_heat);

		//calculating parameter tau 
		double c1 = 4.00;
		double c2 = 2.00;
		double h = sqrt(2.00*Area);
		double norm_u =norm_2(ms_vel_gauss);
		double tau = 1.00 / ( c1*alpha/(h*h) + c2*norm_u/h );

		//getting the BDF2 coefficients (not fixed to allow variable time step)
		//the coefficients INCLUDE the time step
		const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
		
		//CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
		noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
		noalias(rLeftHandSideMatrix) = outer_prod(msN,ms_u_DN);

		//CONVECTION STABILIZING CONTRIBUTION (Suu)
		noalias(rLeftHandSideMatrix) += (tau) * outer_prod(ms_u_DN,ms_u_DN);

		//VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
		noalias(rLeftHandSideMatrix) += (alpha ) * prod(msDN_DX,trans(msDN_DX));

		//INERTIA CONTRIBUTION
		noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * msMassFactors;

		// RHS = Fext (heat per unit mass)
		noalias(rRightHandSideVector) = (heat_flux/specific_heat)*msN;

		//RHS += Suy * proj[component] 
		noalias(rRightHandSideVector) += (tau*proj)*ms_u_DN;  

		//adding the inertia terms
		// RHS += M*vhistory 
		//calculating the historical velocity
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			ms_temp_vec_np[iii] =  BDFcoeffs[1]*GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE,1);
		for(unsigned int step = 2; step<BDFcoeffs.size(); step++)
		{
			for(unsigned int iii = 0; iii<number_of_points; iii++)
				ms_temp_vec_np[iii] += BDFcoeffs[step]*GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE,step);
		}	
		noalias(rRightHandSideVector) -= prod(msMassFactors,ms_temp_vec_np) ;

		//subtracting the dirichlet term
		// RHS -= LHS*temperatures
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);
		
		//multiplying by area, rho and density
		rRightHandSideVector *= (Area * density*specific_heat);
		rLeftHandSideMatrix *= (Area * density*specific_heat);

		
		//Phase change contribution
		if(GetGeometry()[0].FastGetSolutionStepValue(LATENT_HEAT)>0.0)
		{
	
		double Tnew0 = GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
		double Told0 = GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE,1);
		//double Tvold0 = GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE,2);
		double Tnew1 = GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE);
		double Told1 = GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE,1);
		//double Tvold1 = GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE,2);
	        double Tnew2 = GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE);
		double Told2 = GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE,1);
		//double Tvold2 = GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE,2);
		
		double rho = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

		
		double x1 = GetGeometry()[0].X();
		double y1 = GetGeometry()[0].Y();
		double x2 = GetGeometry()[1].X();
		double y2 = GetGeometry()[1].Y();
		double x3 = GetGeometry()[2].X();
		double y3 = GetGeometry()[2].Y();

		double Tm = (GetGeometry()[0].FastGetSolutionStepValue(MELT_TEMPERATURE_1)+GetGeometry()[0].FastGetSolutionStepValue(MELT_TEMPERATURE_2))/2.0;


		double anew;
		if((Tnew0-Tnew2)==0.0)
		{
			double a = 0.0;
			anew = k(a);
			
		}
		else 
		{
			double a = (Tnew0-Tm)/(Tnew0-Tnew2);
			anew = k(a);
			
		}
		
		

		double bnew;
		if((Tnew0-Tnew1)==0.0)
		{
			double b = 0.0;
			bnew = k(b);
		}
		else
		{
			double b = (Tnew0-Tm)/(Tnew0-Tnew1);
			bnew = k(b);
		}
		
		
		double cnew;
		if((Tnew2-Tnew1)==0.0)
		{
			double c = 0.0;
			cnew = k(c);
		}
		else
		{
			double c = (Tnew2-Tm)/(Tnew2-Tnew1);
			cnew = k(c);
		}
 	
	

			
		double xx = (1.0-anew)*x1+anew*x3;
		double yx = (1.0-anew)*y1+anew*y3;
		double xy = (1.0-bnew)*x1+bnew*x2;
		double yy = (1.0-bnew)*y1+bnew*y2;
		double xz = cnew*x2+(1.0-cnew)*x3;
		double yz = cnew*y2+(1.0-cnew)*y3;
		
		
		double W1=25.0/48.0;
		double W2=-27.0/48.0;
		
		Matrix temp_LHS(3,3);
		noalias(temp_LHS)=ZeroMatrix(3,3);
		Vector temp_RHS(3);
		noalias(temp_RHS)=ZeroVector(3);

		/*
		//BDF1 + Gausspoints in nodes
 		//First node
		//if((tnew0 > T1 and tnew0 <= T2) or (told0 > T1 and told0 <= T2)
			DD(0.333333333333*Area, rho, 1.0, 0.0, 0.0, Tnew0, Told0, Tvold0, temp_LHS, temp_RHS,  rCurrentProcessInfo);
		
		//Second node
		//if((tnew1 > T1 and tnew1 <= T2) or (told1 > T1 and told1 <= T2))
			DD(0.333333333333*Area, rho, 0.0, 1.0, 0.0, Tnew1, Told1, Tvold1, temp_LHS, temp_RHS,  rCurrentProcessInfo);
		
		//Third node
		//if((tnew2 > T1 and tnew2 <= T2) or (told2 > T1 and told2 <= T2))
			DD(0.333333333333*Area, rho, 0.0, 0.0, 1.0, Tnew2, Told2, Tvold2, temp_LHS, temp_RHS,  rCurrentProcessInfo);
		*/

		
		//Element Part I
		Matrix E1(3,3);
		E1(0,0)=1.0;      E1(0,1)=0.0;  E1(0,2)=0.0;
		E1(1,0)=1.0-bnew; E1(1,1)=bnew; E1(1,2)=0.0;
		E1(2,0)=1.0-anew; E1(2,1)=0.0;  E1(2,2)=anew;

		

		double AE1=AA(x1,y1,xy,yy,xx,yx);
		
		//Integration point I of Element Part I
		Vector NI(5);
		noalias(NI)=ZeroVector(5);
		NI=Z(E1, 0.6, 0.2, 0.2, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W1*AE1, rho, NI[0], NI[1], NI[2], NI[3], NI[4], temp_LHS, temp_RHS, rCurrentProcessInfo);
		
		//Integration point II of Element Part I
		Vector NII(5);
		noalias(NII)=ZeroVector(5);
		NII=Z(E1, 0.2, 0.6, 0.2, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W1*AE1, rho, NII[0], NII[1], NII[2], NII[3], NII[4], temp_LHS, temp_RHS, rCurrentProcessInfo);
		
		//Integration point 3 of Element Part I
		Vector N3(5);
		noalias(N3)=ZeroVector(5);
		N3=Z(E1, 0.33333333, 0.33333333, 0.33333333, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W2*AE1, rho, N3[0], N3[1], N3[2], N3[3], N3[4], temp_LHS, temp_RHS, rCurrentProcessInfo);
		
		//Integration point 4 of Element Part I
		Vector N4(5);
		noalias(N4)=ZeroVector(5);
		N4=Z(E1, 0.2, 0.2, 0.6, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W1*AE1, rho, N4[0], N4[1], N4[2], N4[3], N4[4], temp_LHS, temp_RHS, rCurrentProcessInfo);
	

		//Element Part II
		Matrix E2(3,3);
		E2(0,0)=1.0-anew;  E2(0,1)=0.0;   E2(0,2)=anew;
		E2(1,0)=1.0-bnew;  E2(1,1)=bnew;  E2(1,2)=0.0;
		E2(2,0)=0.0;       E2(2,1)=cnew;  E2(2,2)=1.0-cnew;

		double AE2=AA(xx,yx,xy,yy,xz,yz);

		//Integration point I of Element Part II
		noalias(NI)=ZeroVector(5);
		NI=Z(E2, 0.6, 0.2, 0.2, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W1*AE2, rho, NI[0], NI[1], NI[2], NI[3], NI[4], temp_LHS, temp_RHS, rCurrentProcessInfo);

		//Integration point II of Element Part II
		noalias(NII)=ZeroVector(5);
		NII=Z(E2, 0.2, 0.6, 0.2, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W1*AE2, rho, NII[0], NII[1], NII[2], NII[3], NII[4], temp_LHS, temp_RHS, rCurrentProcessInfo);
	
		//Integration point 3 of Element Part II
		noalias(N3)=ZeroVector(5);
		N3=Z(E2, 0.33333333, 0.33333333, 0.33333333, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W2*AE2, rho, N3[0], N3[1], N3[2], N3[3], N3[4], temp_LHS, temp_RHS, rCurrentProcessInfo);
		
		//Integration point 4 of Element Part II
		noalias(N4)=ZeroVector(5);
		N4=Z(E2, 0.2, 0.2, 0.6, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W1*AE2, rho, N4[0], N4[1], N4[2], N4[3], N4[4], temp_LHS, temp_RHS, rCurrentProcessInfo);


		
		//Element Part 3
		Matrix E3(3,3);
		E3(0,0)=1.0-bnew;  E3(0,1)=bnew;   E3(0,2)=0.0;
		E3(1,0)=0.0;       E3(1,1)=1.0;    E3(1,2)=0.0;
		E3(2,0)=0.0;       E3(2,1)=cnew;   E3(2,2)=1.0-cnew;
		

		double AE3=AA(x2,y2,xz,yz,xy,yy);
	
		//Integration point I of Element Part 3
		noalias(NI)=ZeroVector(5);
		NI=Z(E3, 0.6, 0.2, 0.2, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W1*AE3, rho, NI[0], NI[1], NI[2], NI[3], NI[4], temp_LHS, temp_RHS, rCurrentProcessInfo);

		//Integration point II of Element Part 3
		noalias(NII)=ZeroVector(5);
		NII=Z(E3, 0.2, 0.6, 0.2, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W1*AE3, rho, NII[0], NII[1], NII[2], NII[3], NII[4], temp_LHS, temp_RHS, rCurrentProcessInfo);

		//Integration point 3 of Element Part 3
		noalias(N3)=ZeroVector(5);
		N3=Z(E3, 0.33333333, 0.33333333, 0.33333333, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W2*AE3, rho, N3[0], N3[1], N3[2], N3[3], N3[4], temp_LHS, temp_RHS, rCurrentProcessInfo);

		//Integration point 4 of Element Part 3
		noalias(N4)=ZeroVector(5);
		N4=Z(E3, 0.2, 0.2, 0.6, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W1*AE3, rho, N4[0], N4[1], N4[2], N4[3], N4[4], temp_LHS, temp_RHS, rCurrentProcessInfo);



		//Element Part 4
		Matrix E4(3,3);
		E4(0,0)=1.0-anew;  E4(0,1)=0.0;    E4(0,2)=anew;
		E4(1,0)=0.0;       E4(1,1)=cnew;   E4(1,2)=1.0-cnew;
		E4(2,0)=0.0;       E4(2,1)=0.0;    E4(2,2)=1.0;


		double AE4=AA(xx,yx,xz,yz,x3,y3);

		//Integration point I of Element Part 4
		noalias(NI)=ZeroVector(5);
		NI=Z(E4, 0.6, 0.2, 0.2, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W1*AE4, rho, NI[0], NI[1], NI[2], NI[3], NI[4], temp_LHS, temp_RHS, rCurrentProcessInfo);

		//Integration point II of Element Part 4
		noalias(NII)=ZeroVector(5);
		NII=Z(E4, 0.2, 0.6, 0.2, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W1*AE4, rho, NII[0], NII[1], NII[2], NII[3], NII[4], temp_LHS, temp_RHS, rCurrentProcessInfo);

		//Integration point 3 of Element Part 4
		noalias(N3)=ZeroVector(5);
		N3=Z(E4, 0.33333333, 0.33333333, 0.33333333, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W2*AE4, rho, N3[0], N3[1], N3[2], N3[3], N3[4], temp_LHS, temp_RHS, rCurrentProcessInfo);

		//Integration point 4 of Element Part 4
		noalias(N4)=ZeroVector(5);
		N4=Z(E4, 0.2, 0.2, 0.6, Tnew0, Tnew1, Tnew2, Told0, Told1, Told2);
		CC(W1*AE4, rho, N4[0], N4[1], N4[2], N4[3], N4[4], temp_LHS, temp_RHS, rCurrentProcessInfo);
			
	

		rRightHandSideVector[0] += temp_RHS[0];
		rRightHandSideVector[1] += temp_RHS[1];
		rRightHandSideVector[2] += temp_RHS[2];
		noalias(rLeftHandSideMatrix) += temp_LHS;
		}
 	
		KRATOS_CATCH("");
	}


	double ConvDiffChangeOfPhase2D::k(double a)
	{
		double anew = 0.5;
		if(a <= 0.05 || a >=0.95)
			return anew;
		else 
			return a;
	}

	//Non-isothermal phase-change function
	double ConvDiffChangeOfPhase2D::f(double T)
	{
		double Tm1 = GetGeometry()[0].FastGetSolutionStepValue(MELT_TEMPERATURE_1);
		double Tm2 = GetGeometry()[0].FastGetSolutionStepValue(MELT_TEMPERATURE_2);
		double L = GetGeometry()[0].FastGetSolutionStepValue(LATENT_HEAT);
		if(T <= Tm1)
			return 0.0;
		else if((Tm1 < T) && (T <= Tm2))
			return (T-Tm1)/(Tm2-Tm1)*L;
		else
			return L;
	}

	//Isothermal phase-change function
	double ConvDiffChangeOfPhase2D::g(double T)
	{
		double Tm = 0.5*(GetGeometry()[0].FastGetSolutionStepValue(MELT_TEMPERATURE_1)+GetGeometry()[0].FastGetSolutionStepValue(MELT_TEMPERATURE_2));
		double L =  GetGeometry()[0].FastGetSolutionStepValue(LATENT_HEAT);
		if(T <= Tm)
			return 0.0;
		else
			return L;
	}
	
	//Area calculation fo current "phase-change" element
	double ConvDiffChangeOfPhase2D::AA(double x1, double y1, double x2, double y2, double x3, double y3)
	{
		
		double A =0.5 * ((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
		return A;
	}
	
	Vector ConvDiffChangeOfPhase2D::Z(Matrix M, double N1, double N2, double N3,double Tnew0, double Tnew1,double Tnew2,double Told0, double Told1, double Told2)
	{
		
		Vector temp(5);
		noalias(temp)=ZeroVector(5);
		temp[0]=N1*M(0,0)+N2*M(1,0)+N3*M(2,0);
		temp[1]=N1*M(0,1)+N2*M(1,1)+N3*M(2,1);
		temp[2]=N1*M(0,2)+N2*M(1,2)+N3*M(2,2);
		temp[3]=Tnew0*temp[0]+Tnew1*temp[1]+Tnew2*temp[2];
		temp[4]=Told0*temp[0]+Told1*temp[1]+Told2*temp[2];
		return temp;
		
	}
		
	
	//Calculation of Contribution of Latent heat for: Gauss Integration  + BDF1 
	void ConvDiffChangeOfPhase2D::CC(double int_weight,double rho, double N1, double N2, double N3, double Tn2, double Tn1, MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
	const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
	double temp = int_weight * rho * (BDFcoeffs[0] * f(Tn2) + BDFcoeffs[1] * f(Tn1));
	
	
	rRightHandSideVector[0]-= N1 * temp;
	rRightHandSideVector[1]-= N2 * temp;
	rRightHandSideVector[2]-= N3 * temp;

	if (fabs(Tn1-Tn2)> 1e-9)
		{
		rLeftHandSideMatrix(0,0) += (N1 * N1)*(temp / (Tn2 - Tn1)); rLeftHandSideMatrix(0,1) += (N1 * N2)*(temp / (Tn2 - Tn1)); rLeftHandSideMatrix(0,2) += (N1 * N3)*(temp / (Tn2 - Tn1));  
		rLeftHandSideMatrix(1,0) += (N2 * N1)*(temp / (Tn2 - Tn1)); rLeftHandSideMatrix(1,1) += (N2 * N2)*(temp / (Tn2 - Tn1)); rLeftHandSideMatrix(1,2) += (N2 * N3)*(temp / (Tn2 - Tn1));  
		rLeftHandSideMatrix(2,0) += (N3 * N1)*(temp / (Tn2 - Tn1)); rLeftHandSideMatrix(2,1) += (N3 * N2)*(temp / (Tn2 - Tn1)); rLeftHandSideMatrix(2,2) += (N3 * N3)*(temp / (Tn2 - Tn1)); 

	
		}
	}
	//Calculation of Contribution of Latent heat for: Gauss Integration  + BDF2 
	void ConvDiffChangeOfPhase2D::DD(double int_weight,double rho, double N1, double N2, double N3,double Tn2,double Tn1,double Tn0, MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
	const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
	double temp = int_weight * rho * (BDFcoeffs[0] * f(Tn2) + BDFcoeffs[1] * f(Tn1));

	for(unsigned int step = 2; step<BDFcoeffs.size(); step++) //for second order accuracy
			temp += int_weight * rho * BDFcoeffs[step]*f(Tn0);

	rRightHandSideVector[0]-= N1 * temp;
	rRightHandSideVector[1]-= N2 * temp;
	rRightHandSideVector[2]-= N3 * temp;

	if (fabs(Tn2-Tn1)> 1e-9)
		{
		rLeftHandSideMatrix(0,0) += (N1 * N1 + N1 * N2 + N1 * N3 ) * temp / (Tn2 - Tn1);   
		rLeftHandSideMatrix(1,1) += (N2 * N1 + N2 * N2 + N2 * N3 ) * temp / (Tn2 - Tn1);  
		rLeftHandSideMatrix(2,2) += (N3 * N1 + N3 * N2 + N3 * N3 ) * temp / (Tn2 - Tn1); 
		}
	}
	
	
	


	//************************************************************************************
	//************************************************************************************
	void ConvDiffChangeOfPhase2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}	 

	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void ConvDiffChangeOfPhase2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
 
		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);		
		
		if(FractionalStepNumber  == 2) //calculation of temperature convective projection
		{
			const unsigned int number_of_points = GetGeometry().size();
			const double lumping_factor = 1.00/double(number_of_points);
			unsigned int TDim = 2;

			//calculating viscosity
			ms_temp_vec_np[0] = GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);			
			const array_1d<double,3>& v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
			const array_1d<double,3>& w = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
			for(unsigned int j = 0; j<TDim; j++)
				ms_vel_gauss[j] = v[j] - w[j];
			
			for(unsigned int i = 1; i<number_of_points; i++)
			{
				ms_temp_vec_np[i] = GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE);				
				const array_1d<double,3>& v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
				const array_1d<double,3>& w = GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);
				for(unsigned int j = 0; j<TDim; j++)
					ms_vel_gauss[j] += v[j] - w[j];
				
			}
			ms_vel_gauss *= lumping_factor;

			//calculating convective auxiliary vector
			noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
			double temp_conv = inner_prod( ms_u_DN , ms_temp_vec_np);
			temp_conv *= Area;

			for(unsigned int i = 0; i<number_of_points; i++)
			{
				GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += lumping_factor*Area;
				GetGeometry()[i].FastGetSolutionStepValue(TEMP_CONV_PROJ) += lumping_factor*temp_conv; ;
			}
		}
		KRATOS_CATCH("");
	}


	//************************************************************************************
	//************************************************************************************
	void ConvDiffChangeOfPhase2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	  void ConvDiffChangeOfPhase2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);

	}

} // Namespace Kratos


