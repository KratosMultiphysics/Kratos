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
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.4 $
//
//
 

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/conv_diff_3d.h"
#include "convection_diffusion_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
// 	//space defined to allow parallelism
// 	 namespace ConvDiff3DAuxiliaries
// 	 {
// 	     boost::numeric::ublas::bounded_matrix<double,4,4> msMassFactors = 0.25*IdentityMatrix(4,4);
// 	     #pragma omp threadprivate(msMassFactors)
// 
// 	     boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX = ZeroMatrix(4,3);
// 	     #pragma omp threadprivate(msDN_DX)
// 
// 	     array_1d<double,4> msN = ZeroVector(4); //dimension = number of nodes
// 	     #pragma omp threadprivate(msN)
// 
// 	     array_1d<double,3> ms_vel_gauss = ZeroVector(3); //dimesion coincides with space dimension
// 	     #pragma omp threadprivate(ms_vel_gauss)
// 
// 	     array_1d<double,4> ms_temp_vec_np = ZeroVector(4); //dimension = number of nodes
// 	     #pragma omp threadprivate(ms_temp_vec_np)
// 
// 	     array_1d<double,4> ms_u_DN = ZeroVector(4); //dimension = number of nodes
// 	     #pragma omp threadprivate(ms_u_DN)
// 
// 	 }
// 	 using namespace ConvDiff3DAuxiliaries;

	//************************************************************************************
	//************************************************************************************
	ConvDiff3D::ConvDiff3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	ConvDiff3D::ConvDiff3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer ConvDiff3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new ConvDiff3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	ConvDiff3D::~ConvDiff3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void ConvDiff3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int number_of_points = GetGeometry().size();
		const double lumping_factor = 1.00/double(number_of_points);
		unsigned int TDim = 3;

		const boost::numeric::ublas::bounded_matrix<double,4,4> msMassFactors = 0.25*IdentityMatrix(4,4);
		boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX;
		array_1d<double,4> msN;
		array_1d<double,3> ms_vel_gauss;
		array_1d<double,4> ms_temp_vec_np;
		array_1d<double,4> ms_u_DN ;
	     
		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);

		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);
		ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

		//calculating viscosity
		const Variable<double>&  rDensityVar = my_settings->GetDensityVariable();
		const Variable<array_1d<double,3> >& rConvectionVar = my_settings->GetConvectionVariable();
		const Variable<double>& rDiffusionVar =my_settings->GetDiffusionVariable();
		const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
        	const Variable<double>& rSourceVar =my_settings->GetVolumeSourceVariable();
        	const Variable<double>& rSurfaceSourceVar =my_settings->GetSurfaceSourceVariable();
		const Variable<array_1d<double,3> >& rMeshVelocityVar =my_settings->GetMeshVelocityVariable();





		double conductivity = GetGeometry()[0].FastGetSolutionStepValue(rDiffusionVar);
		double density = GetGeometry()[0].FastGetSolutionStepValue(rDensityVar);
		double specific_heat = GetGeometry()[0].FastGetSolutionStepValue(SPECIFIC_HEAT);
		double heat_flux = GetGeometry()[0].FastGetSolutionStepValue(rSourceVar);
		double proj = GetGeometry()[0].FastGetSolutionStepValue(TEMP_CONV_PROJ);
		
		const array_1d<double,3>& v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar);
		for(unsigned int j = 0; j<TDim; j++)
			ms_vel_gauss[j] = v[j] - w[j];
		
		for(unsigned int i = 1; i<number_of_points; i++)
		{
			conductivity += GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
			density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
			specific_heat += GetGeometry()[i].FastGetSolutionStepValue(SPECIFIC_HEAT);
			heat_flux += GetGeometry()[i].FastGetSolutionStepValue(rSourceVar);
			proj += GetGeometry()[i].FastGetSolutionStepValue(TEMP_CONV_PROJ);
			
			const array_1d<double,3>& v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
			const array_1d<double,3>& w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
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
		double h =  pow(6.00*Area,0.3333333);
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

		// RHS = Fext
		noalias(rRightHandSideVector) = (heat_flux)*msN;

		//RHS += Suy * proj[component] 
		noalias(rRightHandSideVector) += (tau*proj)*ms_u_DN;  

		//adding the inertia terms
		// RHS += M*vhistory 
		//calculating the historical velocity
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			ms_temp_vec_np[iii] =  BDFcoeffs[1]*GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar,1);
		for(unsigned int step = 2; step<BDFcoeffs.size(); step++)
		{
			for(unsigned int iii = 0; iii<number_of_points; iii++)
				ms_temp_vec_np[iii] += BDFcoeffs[step]*GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar,step);
		}	
		noalias(rRightHandSideVector) -= prod(msMassFactors,ms_temp_vec_np) ;

		//subtracting the dirichlet term
		// RHS -= LHS*temperatures
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);
		
		//multiplying by area, rho and density
		rRightHandSideVector *= (Area * density*specific_heat);
		rLeftHandSideMatrix *= (Area * density*specific_heat);
//KRATOS_WATCH(rLeftHandSideMatrix)
//KRATOS_WATCH(rRightHandSideVector)
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void ConvDiff3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}	 

	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void ConvDiff3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
		
// 		const boost::numeric::ublas::bounded_matrix<double,4,4> msMassFactors = 0.25*IdentityMatrix(4,4);
		boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX;
		array_1d<double,4> msN;
		array_1d<double,3> ms_vel_gauss;
		array_1d<double,4> ms_temp_vec_np;
		array_1d<double,4> ms_u_DN ;
 
		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);		
		ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
		const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
		const Variable<array_1d<double,3> >& rConvectionVar = my_settings->GetConvectionVariable();
		const Variable<array_1d<double,3> >& rMeshVelocityVar =my_settings->GetMeshVelocityVariable();


		
		if(FractionalStepNumber  == 2) //calculation of temperature convective projection
		{
			const unsigned int number_of_points = GetGeometry().size();
			const double lumping_factor = 1.00/double(number_of_points);
			unsigned int TDim = 3;

			//calculating viscosity
			ms_temp_vec_np[0] = GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);			
			const array_1d<double,3>& v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
			const array_1d<double,3>& w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar);
			for(unsigned int j = 0; j<TDim; j++)
				ms_vel_gauss[j] = v[j] - w[j];
			
			for(unsigned int i = 1; i<number_of_points; i++)
			{
				ms_temp_vec_np[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);				
				const array_1d<double,3>& v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
				const array_1d<double,3>& w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
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
	void ConvDiff3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{

		ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
		const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	  void ConvDiff3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
		const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();

		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

	}

} // Namespace Kratos


