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
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2009-01-21 14:15:02 $
//   Revision:            $Revision: 1.6 $
//
//
 
//#define GRADPN_FORM
//#define STOKES

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/asgs_2d.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{



	//************************************************************************************
	//************************************************************************************
	ASGSPRDC::ASGSPRDC(IndexType NewId, GeometryType::Pointer pGeometry)
		: ASGS2D(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	ASGSPRDC::ASGSPRDC(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: ASGS2D(NewId, pGeometry, pProperties)
	{
		
	}

	Element::Pointer ASGSPRDC::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		
		KRATOS_TRY
		return Element::Pointer(new ASGSPRDC(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	ASGSPRDC::~ASGSPRDC()
	{
	}

	//************************************************************************************
	//************************************************************************************


	void ASGSPRDC::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

	int nodes_number = 3;
	int dim = 2;
	unsigned int matsize = nodes_number*(dim+1);

	if(rLeftHandSideMatrix.size1() != matsize)
			rLeftHandSideMatrix.resize(matsize,matsize); //false says not to preserve existing storage!!

	if(rRightHandSideVector.size() != matsize)
			rRightHandSideVector.resize(matsize); //false says not to preserve existing storage!!

	
	noalias(rLeftHandSideMatrix) = ZeroMatrix(matsize,matsize); 
	noalias(rRightHandSideVector) = ZeroVector(matsize); 
	
	double delta_t= rCurrentProcessInfo[DELTA_TIME];
	
	boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
        array_1d<double,3> N; 
        array_1d<double,2> ms_adv_vel;	
	
	//getting data for the given geometry
	double Area;
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);
	

	//mass & viscous term	
	//double timearea = Area * delta_t;
	CalculateMassContribution(rLeftHandSideMatrix,delta_t,Area); 
	CalculateViscousTerm(rLeftHandSideMatrix, DN_DX, Area);
	
	//Advective term
	double tauone;
	double tautwo;
	CalculateTau(tauone, tautwo, delta_t, Area, rCurrentProcessInfo);

	CalculateAdvectiveTerm(rLeftHandSideMatrix, DN_DX, tauone, tautwo, delta_t, Area);
		
	//calculate pressure term
	CalculatePressureTerm(rLeftHandSideMatrix, DN_DX, N, delta_t,Area);

	//compute projections
	
	//stabilization terms
	CalculateDivStblTerm(rLeftHandSideMatrix, DN_DX, tautwo, Area);
	CalculateAdvStblAllTerms(rLeftHandSideMatrix,rRightHandSideVector, DN_DX, N, tauone,delta_t, Area);
	CalculateGradStblAllTerms(rLeftHandSideMatrix,rRightHandSideVector,DN_DX, delta_t, tauone, Area);
	//KRATOS_WATCH(rRightHandSideVector);
	
	//add body force and momentum
	//AddBodyForceAndMomentum(rRightHandSideVector, N, delta_t, Area);
	AddBodyForceAndMomentum(rRightHandSideVector, N, delta_t, Area, tauone,tautwo);

	//add projections
	//AddProjectionForces(rRightHandSideVector,DN_DX,Area);	
	
	//calculate residual
	CalculateResidual(rLeftHandSideMatrix, rRightHandSideVector);


	KRATOS_CATCH("")

	}

	//************************************************************************************
	//************************************************************************************
	void ASGSPRDC::CalculateResidual(const MatrixType& K, VectorType& F)
	{
	KRATOS_TRY

	int nodes_number = 3;
	int dof = 2;


	array_1d<double,9> UP = ZeroVector(9);
	for ( int ii = 0; ii < nodes_number; ++ii)
	    {
		int index = ii * (dof + 1);
		UP[index] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY,0)[0];
		UP[index + 1] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY,0)[1];
		UP[index + 2] = GetGeometry()[ii].FastGetSolutionStepValue(AIR_PRESSURE,0);
	    }
	F -= prod(K,UP);

	KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void ASGSPRDC::ComputeProjections(array_1d<double,6>& adv_proj , array_1d<double,3>& div_proj, const 			boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const double tauone,const double tautwo,const array_1d<double,3>& N,const double area, const double time)
	{
	unsigned int number_of_nodes = GetGeometry().PointsNumber();
	unsigned int dim = 2;
	
	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	const array_1d<double,3>& adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double,3>& adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double,3>& adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);

	array_1d<double,2> mean_adv_vel; 

	mean_adv_vel[0] = N[0]*(adv_vel0[0]-mesh_vel0[0])+N[1]*(adv_vel1[0]-mesh_vel1[0])+N[2]*(adv_vel2[0]-mesh_vel2[0]);
	mean_adv_vel[1] = N[0]*(adv_vel0[1]-mesh_vel0[1])+N[1]*(adv_vel1[1]-mesh_vel1[1])+N[2]*(adv_vel2[1]-mesh_vel2[1]);


	double const_adv_proj_X =0.0;
	double const_adv_proj_Y =0.0;
	double mean_div_proj =0.0;
	array_1d<double,3> mean_new_vel = ZeroVector(3);
	array_1d<double,3> mean_old_vel = ZeroVector(3);
	array_1d<double,3> mean_bdf = ZeroVector(3);


	for (unsigned int i=0;i<number_of_nodes;i++)
	 {
		
		//int index = i*dim;
		double pr = GetGeometry()[i].FastGetSolutionStepValue(AIR_PRESSURE);
		const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& old_vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
		const array_1d<double,3>& bdf = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);

		//const array_1d<double,2>& bdf = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
		// to consider the jump gradp/ro is calculated
			//pr = pr/density;

		//adv_proj = PI(ro*dv/dt + ro*a.gradU + gradP - f) considering lumped mass PI() = ()
		//calculate constant part of RES ->ro*a.gradU + gradP
		const_adv_proj_X += ( pr * DN_DX(i,0) + density*(mean_adv_vel[0]*DN_DX(i,0) + mean_adv_vel[1]*DN_DX(i,1))*vel[0] );
		const_adv_proj_Y +=  (pr * DN_DX(i,1) + density*(mean_adv_vel[0]*DN_DX(i,0) + mean_adv_vel[1]*DN_DX(i,1))*vel[1] );

		//div_proj = PI(ro*divU)
		mean_div_proj += density*(DN_DX(i,0)*vel[0] + DN_DX(i,1)*vel[1]);

		//calcuale mean velocity and body force
		mean_new_vel += 0.3333333333333333333333333333*vel;
		mean_old_vel += 0.3333333333333333333333333333*old_vel;
		mean_bdf += 0.3333333333333333333333333333*bdf;
	}



	for (unsigned int i=0;i<number_of_nodes;i++)
	 {
		int index = i*dim;

		adv_proj[index] = area*N[i]*(density*(mean_new_vel[0]-mean_old_vel[0])/time + const_adv_proj_X - density*mean_bdf[0]);
		adv_proj[index +1] = area*N[i]*(density*(mean_new_vel[1]-mean_old_vel[1])/time + const_adv_proj_Y - density*mean_bdf[1]);

		div_proj[i] = area*N[i]*density*mean_div_proj;

		//update projections
		array_1d<double,3>& advtermproj = GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
		advtermproj[0]+= adv_proj[index];
		advtermproj[1]+= adv_proj[index +1];

		double& divtermproj = GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ);
		divtermproj += div_proj[i] ;

		

		//calculate nodal area
		
		GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += 0.333333333333333333*area;
	}
	/*for (unsigned int i=0;i<number_of_nodes;i++)
	 {
		int index = i*dim;
		const double pr = GetGeometry()[i].FastGetSolutionStepValue(AIR_PRESSURE);
		const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);

		const array_1d<double,3>& mesh_vel = GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);
		array_1d<double,2> adv_vel = ZeroVector(2);
		adv_vel[0] = vel[0]-mesh_vel[0];
		adv_vel[1] = vel[1]-mesh_vel[1];

		//adv_proj = PI(ro*a.gradU + gradP) considering lumped mass PI() = ()
		adv_proj[index] = area*N[i]*(pr * DN_DX(i,0) +density*(adv_vel[0]*DN_DX(i,0) + adv_vel[1]*DN_DX(i,1))*vel[0]);
		adv_proj[index +1] =  area*N[i]*(pr * DN_DX(i,1) + density*(adv_vel[0]*DN_DX(i,0) + adv_vel[1]*DN_DX(i,1))*vel[1]);

		//div_proj = PI(ro*divU)
		div_proj[i] = area*N[i]*density*(DN_DX(i,0)*vel[0] + DN_DX(i,1)*vel[1]);

		//update projections
		array_1d<double,3>& advtermproj = GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
		advtermproj[0]+= tauone*adv_proj[index];
		advtermproj[1]+= tauone*adv_proj[index +1];

		double& divtermproj = GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ);
		divtermproj += tautwo*div_proj[i] ;
		
		//calculate nodal area
		GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += 0.333333333333333333*area;
	 }*/

	
	}

	//************************************************************************************
	//************************************************************************************
	void ASGSPRDC::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;
		unsigned int node_size = dim+1;

		
			if(rResult.size() != number_of_nodes*node_size)
				rResult.resize(number_of_nodes*node_size,false);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				rResult[i*node_size] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
				rResult[i*node_size+1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
				rResult[i*node_size+2] = GetGeometry()[i].GetDof(AIR_PRESSURE).EquationId();
			}
		KRATOS_CATCH("")
		
	}

	//************************************************************************************
	//************************************************************************************
	  void ASGSPRDC::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;
		unsigned int node_size = dim+1;

			if(ElementalDofList.size() != number_of_nodes*node_size)
				ElementalDofList.resize(number_of_nodes*node_size);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				ElementalDofList[i*node_size] = GetGeometry()[i].pGetDof(VELOCITY_X);
				ElementalDofList[i*node_size+1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
				ElementalDofList[i*node_size+2] = GetGeometry()[i].pGetDof(AIR_PRESSURE);
			}
		KRATOS_CATCH("");
	}

	//*************************************************************************************
	//*************************************************************************************
	void ASGSPRDC::calculatedensity(Geometry< Node<3> > geom, double& density, double& viscosity)
	{

	/*double kk = 0.0;
	for(int ii=0;ii<3;++ii)
		if(geom[ii].GetSolutionStepValue(IS_STRUCTURE) != 1.0)
			{
				kk++;
				density +=geom[ii].FastGetSolutionStepValue(DENSITY);
			}

	density/=kk;*/
	/*
		density = ZeroVector(3);
	for(int ii=0;ii<3;++ii)
		density[ii] = geom[ii].FastGetSolutionStepValue(DENSITY);*/
	
	const double rho0 = geom[0].FastGetSolutionStepValue(DENSITY_AIR);
	const double rho1 = geom[1].FastGetSolutionStepValue(DENSITY_AIR);
	const double rho2 = geom[2].FastGetSolutionStepValue(DENSITY_AIR);
	 density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );

	const double visc0 = geom[0].FastGetSolutionStepValue(VISCOSITY_AIR);
	const double visc1 = geom[1].FastGetSolutionStepValue(VISCOSITY_AIR);
	const double visc2 = geom[2].FastGetSolutionStepValue(VISCOSITY_AIR);
	 viscosity = 0.3333333333333333333333*(visc0 + visc1 + visc2 );

/*
	double first = geom[0].FastGetSolutionStepValue(IS_WATER);
	double second = geom[1].FastGetSolutionStepValue(IS_WATER);
	double third = geom[2].FastGetSolutionStepValue(IS_WATER);

	density = 0.0;

	if(first == second && second==third)
	  {
		//for inside the domain totally inside one fluid
		density = geom[0].FastGetSolutionStepValue(DENSITY);	
		viscosity = geom[0].FastGetSolutionStepValue(VISCOSITY);	
	  }
	else
	  {
		//for element having common node between two fluids or boundary element with IS_WATER==1 inside the domain
		for(int ii=0;ii<3;++ii)
			{
			  if(geom[ii].GetSolutionStepValue(IS_WATER) == 1.0 && geom[ii].GetSolutionStepValue(IS_STRUCTURE) != 1.0)
				{
			  	 density = geom[ii].FastGetSolutionStepValue(DENSITY_WATER);
				 viscosity = geom[ii].FastGetSolutionStepValue(VISCOSITY_WATER);

				}
			}
		//for boundary element with IS_WATER==1 on the boundary
		if(density == 0.0)
			for(int ii=0;ii<3;++ii)	
			 {
			  if(geom[ii].GetSolutionStepValue(IS_WATER) == 0.0)
				{
				 density = geom[ii].FastGetSolutionStepValue(DENSITY_AIR);
				 viscosity = geom[ii].FastGetSolutionStepValue(VISCOSITY_AIR);

			
				}
			 }	
			
	  }
*/

	//Here we calculate Dynamic viscosity from Kinemeatic viscosity
	viscosity *= density;

	}
	//*************************************************************************************
	//*************************************************************************************

} // Namespace Kratos


