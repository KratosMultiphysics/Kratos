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
#include "custom_elements/conv_diff_2d.h"
#include "convection_diffusion_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
    namespace ConvDiff2Dauxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,3,3> msMassFactors = 1.0/3.0*IdentityMatrix(3,3);
        #pragma omp threadprivate(msMassFactors)

        boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
        #pragma omp threadprivate(msDN_DX)

        array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
        #pragma omp threadprivate(msN)

        array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension
        #pragma omp threadprivate(ms_vel_gauss)

        array_1d<double,3> ms_temp_vec_np = ZeroVector(3); //dimension = number of nodes
        #pragma omp threadprivate(ms_temp_vec_np)

        array_1d<double,3> ms_u_DN = ZeroVector(3); //dimension = number of nodes
        #pragma omp threadprivate(ms_u_DN)

    }
    using  namespace ConvDiff2Dauxiliaries;

	//************************************************************************************
	//************************************************************************************
	ConvDiff2D::ConvDiff2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	ConvDiff2D::ConvDiff2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
		
	}

	Element::Pointer ConvDiff2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new ConvDiff2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	ConvDiff2D::~ConvDiff2D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void ConvDiff2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
		//double tau = 1.00 / ( c1*alpha/(h*h) + c2*norm_u/h );
		double tau1=( h*h )/(c1 * conductivity + c2 * density * specific_heat * norm_u * h);
		//getting the BDF2 coefficients (not fixed to allow variable time step)
		//the coefficients INCLUDE the time step
		const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
		
		//CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
		noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
		noalias(rLeftHandSideMatrix) = (density*specific_heat) * outer_prod(msN,ms_u_DN);

		//CONVECTION STABILIZING CONTRIBUTION (Suu)
		noalias(rLeftHandSideMatrix) += density * specific_heat * density * specific_heat * tau1 * outer_prod(ms_u_DN,ms_u_DN);

		//VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
		noalias(rLeftHandSideMatrix) += (conductivity )* prod(msDN_DX,trans(msDN_DX));

                		//filling the mass factors
//	msMassFactors(0,0) = 1.00/6.00;  msMassFactors(0,1) = 1.00/12.00; msMassFactors(0,2) = 1.00/12.00;
//		msMassFactors(1,0) = 1.00/12.00; msMassFactors(1,1) = 1.00/6.00;  msMassFactors(1,2) = 1.00/12.00;
//		msMassFactors(2,0) = 1.00/12.00; msMassFactors(2,1) = 1.00/12.00; msMassFactors(2,2) = 1.00/6.00;
		msMassFactors(0,0) = 1.00/3.00; msMassFactors(0,1) = 0.00;		msMassFactors(0,2) = 0.00;
		msMassFactors(1,0) = 0.00;		msMassFactors(1,1) = 1.00/3.00; msMassFactors(1,2) = 0.00;
		msMassFactors(2,0) = 0.00;		msMassFactors(2,1) = 0.00;		msMassFactors(2,2) = 1.00/3.00;

		//INERTIA CONTRIBUTION
		noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * (density*specific_heat) * msMassFactors;

		// RHS = Fext (heat per unit mass)
		noalias(rRightHandSideVector) = (heat_flux * density )*msN;

		//RHS += Suy * proj[component] 
		noalias(rRightHandSideVector) += density*specific_heat * density*specific_heat * (tau1*proj)*ms_u_DN;  

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
		noalias(rRightHandSideVector) -= prod(msMassFactors,ms_temp_vec_np*density*specific_heat) ;

		//subtracting the dirichlet term
		// RHS -= LHS*temperatures
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);
		
		//multiplying by area, rho and density
		rRightHandSideVector *= Area;
		rLeftHandSideMatrix *= Area;
		//multiplying by area
		//rRightHandSideVector *= Area;
		//rLeftHandSideMatrix *= Area;
//KRATOS_WATCH(rLeftHandSideMatrix)
//KRATOS_WATCH(rRightHandSideVector)
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void ConvDiff2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}	 

	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void ConvDiff2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
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
	void ConvDiff2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	  void ConvDiff2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);

	}

} // Namespace Kratos


