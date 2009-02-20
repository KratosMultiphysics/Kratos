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
#include "custom_elements/fluid_2d_asgs.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
        namespace Fluid2DASGSauxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
        #pragma omp threadprivate(DN_DX)

        array_1d<double,3> N = ZeroVector(3); //dimension = number of nodes
        #pragma omp threadprivate(N)

        array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
        #pragma omp threadprivate(ms_adv_vel)

    }
    using  namespace Fluid2DASGSauxiliaries;


	//************************************************************************************
	//************************************************************************************
	Fluid2DASGS::Fluid2DASGS(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Fluid2DASGS::Fluid2DASGS(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
		
	}

	Element::Pointer Fluid2DASGS::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		
		KRATOS_TRY
		return Element::Pointer(new Fluid2DASGS(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	Fluid2DASGS::~Fluid2DASGS()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void Fluid2DASGS::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
	
	
	
	//getting data for the given geometry
	double Area;
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);
	

	//mass & viscous term	
	//double timearea = Area * delta_t;
	CalculateMassContribution(rLeftHandSideMatrix,delta_t,Area); 
	CalculateViscousTerm(rLeftHandSideMatrix, DN_DX, Area);
	
	//Advective term
	double thawone;
	double thawtwo;
	CalculateThaw(thawone, thawtwo, delta_t, Area, rCurrentProcessInfo);

	CalculateAdvectiveTerm(rLeftHandSideMatrix, DN_DX, thawone, thawtwo, delta_t, Area);
		
	//calculate pressure term
	CalculatePressureTerm(rLeftHandSideMatrix, DN_DX, N, delta_t,Area);

	//compute projections
	
	//stabilization terms
	CalculateDivStblTerm(rLeftHandSideMatrix, DN_DX, thawtwo, Area);
	CalculateAdvStblAllTerms(rLeftHandSideMatrix,rRightHandSideVector, DN_DX, N, thawone,delta_t, Area);
	CalculateGradStblAllTerms(rLeftHandSideMatrix,rRightHandSideVector,DN_DX, delta_t, thawone, Area);
	//KRATOS_WATCH(rRightHandSideVector);
	
	//add body force and momentum
	AddBodyForceAndMomentum(rRightHandSideVector, N, delta_t, Area);

	//add projections
	//AddProjectionForces(rRightHandSideVector,DN_DX,Area);	
	
	//calculate residual
	CalculateResidual(rLeftHandSideMatrix, rRightHandSideVector);

	

		KRATOS_CATCH("")
	}
	//***********************************************************************************++
	//**************************************************************************************++
	void Fluid2DASGS::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		MatrixType temp = Matrix();
		//CalculateLocalSystem(temp, rRightHandSideVector,  rCurrentProcessInfo);

		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DASGS::CalculateMassContribution(MatrixType& K,const double time,const double area)
	{
		KRATOS_TRY
	double lump_mass_fac = area * 0.333333333333333333333333;
	double density;
	calculatedensity(GetGeometry(), density);

	int nodes_number = 3;
	int dof = 2;
	for ( int nd = 0; nd< nodes_number; nd++)
	    {
		int row = nd*(dof + 1);
		for( int jj=0; jj< dof; jj++)
			K(row + jj, row + jj) += density/time*lump_mass_fac;
	    }
	
		KRATOS_CATCH("")
	
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DASGS::CalculateViscousTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double area)
	{
		KRATOS_TRY
	double mu;
	const double mu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
	const double mu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
	const double mu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
	mu = 0.333333333333333333333333*(mu0 + mu1 + mu2);
	
	double density;
	calculatedensity(GetGeometry(), density);


	//nu = nu/density;	

	int nodes_number = 3;
	int dof = 2;

	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);
		K(row,column) += mu*1*area*(DN_DX(ii,0)*DN_DX(jj,0) + DN_DX(ii,1)*DN_DX(jj,1));
		K(row + 1,column + 1) += mu*1*area*(DN_DX(ii,0)*DN_DX(jj,0) + DN_DX(ii,1)*DN_DX(jj,1));
		   }	
	    }
					

		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DASGS::CalculateAdvectiveTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double thawone, const double thawtwo, const double time,const double area)
	{
		KRATOS_TRY
	//calculate mean advective velocity and thaws
	const array_1d<double,3>& adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double,3>& adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double,3>& adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);

	double density;
	calculatedensity(GetGeometry(), density);

	ms_adv_vel[0] = N[0]*(adv_vel0[0]-mesh_vel0[0])+N[1]*(adv_vel1[0]-mesh_vel1[0])+N[2]*(adv_vel2[0]-mesh_vel2[0]);
	ms_adv_vel[1] = N[0]*(adv_vel0[1]-mesh_vel0[1])+N[1]*(adv_vel1[1]-mesh_vel1[1])+N[2]*(adv_vel2[1]-mesh_vel2[1]);

	//ms_adv_vel[0] = 0.0;
	//ms_adv_vel[1] = 0.0;
	                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

	//calculate convective term	
	int nodes_number = 3;
	int dof = 2;
	int matsize = dof*nodes_number;

	boost::numeric::ublas::bounded_matrix<double,2,6> conv_opr = ZeroMatrix(dof,matsize);
	boost::numeric::ublas::bounded_matrix<double,6,2> shape_func = ZeroMatrix(matsize, dof);

	for (int ii = 0; ii< nodes_number; ii++)
	    {
		int column = ii*dof;
		conv_opr(0,column) = DN_DX(ii,0)*ms_adv_vel[0] + DN_DX(ii,1)*ms_adv_vel[1];
		conv_opr(1,column + 1) = conv_opr(0,column);

		shape_func(column,0) = N[ii];
		shape_func(column + 1, 1) = shape_func(column,0);
	    }
	boost::numeric::ublas::bounded_matrix<double,6,6> temp_convterm = ZeroMatrix(matsize,matsize);
	temp_convterm = prod(shape_func, conv_opr);

	//double fac = thawone/time;
	//fac = 0.0; 
	//temp_convterm *= ((1 + fac*density)*density); // For the simplicity of implementation the stabilization term thaw1*ro/deltat*(a.gradV U(n+1,t)) is added 	
	
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		int loc_row = ii*dof;
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);
			int loc_column = jj*dof;

			K(row,column) += area*density*temp_convterm(loc_row,loc_column);
			K(row + 1,column + 1) += area*density*temp_convterm(loc_row + 1,loc_column + 1);
		   }
	    }

		KRATOS_CATCH("")

	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DASGS::CalculatePressureTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>&  N, const double time,const double area)
	{
		KRATOS_TRY
	int nodes_number = 3;
	int dof = 2;

	double density;
	calculatedensity(GetGeometry(), density);


	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;

			K(row,column) += -1*area * N(jj) * DN_DX(ii,0);
			K(column,row) += 1*area *density* N(jj) * DN_DX(ii,0);

			K(row + 1,column) += -1*area * N(jj) * DN_DX(ii,1);
			K(column,row + 1) += 1*area *density* N(jj) * DN_DX(ii,1);
		   }
	    }

	
		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DASGS::CalculateDivStblTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double thawtwo,const double area)
	{
		KRATOS_TRY
	int nodes_number = 3;
	int dof = 2;
	int matsize = dof*nodes_number;

	boost::numeric::ublas::bounded_matrix<double,1,6> div_opr = ZeroMatrix(1,matsize);
	for(int ii=0; ii<nodes_number; ii++)
	  {
		int index = dof*ii;
		div_opr(0,index) = DN_DX(ii,0);
		div_opr(0,index + 1) = DN_DX(ii,1);
	  }


	double density;
	calculatedensity(GetGeometry(), density);


	boost::numeric::ublas::bounded_matrix<double,6,6> temp_div = ZeroMatrix(matsize,matsize);
	temp_div = thawtwo * prod(trans(div_opr),div_opr);

	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		int loc_row = ii*dof;
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);
			int loc_column = jj*dof;

			K(row,column) += 1*area*density*temp_div(loc_row,loc_column);
			K(row,column + 1) += 1*area*density*temp_div(loc_row,loc_column + 1);
			K(row + 1,column) += 1*area*density*temp_div(loc_row + 1,loc_column);
			K(row + 1,column + 1) += 1*area*density*temp_div(loc_row + 1,loc_column + 1);
		   }
	    }

		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DASGS::CalculateAdvStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N, const double thawone,const double time,const double area)
	{
		KRATOS_TRY
	const array_1d<double,3>& adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double,3>& adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double,3>& adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);

	ms_adv_vel[0] = N[0]*(adv_vel0[0]-mesh_vel0[0])+N[1]*(adv_vel1[0]-mesh_vel1[0])+N[2]*(adv_vel2[0]-mesh_vel2[0]);
	ms_adv_vel[1] = N[0]*(adv_vel0[1]-mesh_vel0[1])+N[1]*(adv_vel1[1]-mesh_vel1[1])+N[2]*(adv_vel2[1]-mesh_vel2[1]);

		//ms_adv_vel[0] = 0.0;
		//ms_adv_vel[1] = 0.0;

	//calculate convective term	
	int nodes_number = 3;
	int dof = 2;
	int matsize = dof*nodes_number;

	boost::numeric::ublas::bounded_matrix<double,2,6> conv_opr = ZeroMatrix(dof,matsize);
	boost::numeric::ublas::bounded_matrix<double,6,2> shape_func = ZeroMatrix(matsize, dof);

	for (int ii = 0; ii< nodes_number; ii++)
	    {
		int column = ii*dof;
		conv_opr(0,column) = DN_DX(ii,0)*ms_adv_vel[0] + DN_DX(ii,1)*ms_adv_vel[1];
		conv_opr(1,column + 1) = conv_opr(0,column);
		
		shape_func(column,0) = N[ii];
		shape_func(column + 1, 1) = shape_func(column,0);
	    }

	//build (a.grad V)(ro*a.grad U) stabilization term & assemble
	boost::numeric::ublas::bounded_matrix<double,6,6> adv_stblterm = ZeroMatrix(matsize,matsize);		
	adv_stblterm = thawone * prod(trans(conv_opr),conv_opr);


	double density;
	calculatedensity(GetGeometry(), density);

	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		int loc_row = ii*dof;
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);
			int loc_column = jj*dof;
			
			K(row,column) += 1*area*density*adv_stblterm(loc_row,loc_column);
			K(row,column + 1) += 1*area*density*adv_stblterm(loc_row,loc_column + 1);
			K(row + 1,column) += 1*area*density*adv_stblterm(loc_row + 1,loc_column);
			K(row + 1,column + 1) += 1*area*density*adv_stblterm(loc_row + 1,loc_column + 1);
		   }
	    }
		
	//build 1*thaw1*(a.grad V)(grad P) & 1*thaw1*(grad q)(ro*a.grad U) stabilization terms & assemble
	boost::numeric::ublas::bounded_matrix<double,6,3> grad_stblterm = ZeroMatrix(matsize,nodes_number);
	grad_stblterm = thawone * prod(trans(conv_opr),trans(DN_DX));

	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		int loc_row = ii*dof;
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;

			K(row,column) += 1.0*area *1.0* grad_stblterm(loc_row,jj);
			K(row + 1,column) += 1.0*area *1.0* grad_stblterm(loc_row + 1, jj);

			K(column, row) += 1.0*area * density*grad_stblterm(loc_row, jj);
			K(column, row + 1) += 1.0*area *density* grad_stblterm(loc_row + 1, jj);
		   }
	    }

	//thaw1*ro/dt*U(n+1,i+1).(1.0*a.grad V)
	boost::numeric::ublas::bounded_matrix<double,6,6> temp_convterm = ZeroMatrix(matsize,matsize);
	temp_convterm = prod(trans(conv_opr),trans(shape_func));

	double fac = thawone/time*density;	
	
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		int loc_row = ii*dof;
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);
			int loc_column = jj*dof;

			K(row,column) += area*fac*temp_convterm(loc_row,loc_column);
			K(row + 1,column + 1) += area*fac*temp_convterm(loc_row + 1,loc_column + 1);
		   }
	    }


	//build (1.0*a.grad V) (Fbody + ro/dt * U(n)) stabilization term & assemble
	array_1d<double,2> bdf = ZeroVector(2);
	const array_1d<double,2> bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double,2> bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double,2> bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);

	const array_1d<double,3>& vel0_n = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
	const array_1d<double,3>& vel1_n = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
	const array_1d<double,3>& vel2_n = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);

	bdf[0] = N[0]*(density*bdf0[0] + density/time * vel0_n[0] ) +  N[1]*(density*bdf1[0] + density/time * vel1_n[0]) + N[2]*(density*bdf2[0] + density/time * vel2_n[0]);
	bdf[1] =  N[0]*(density*bdf0[1] + density/time * vel0_n[1] ) +  N[1]*(density*bdf1[1] + density/time * vel1_n[1]) + N[2]*(density*bdf2[1] + density/time * vel2_n[1]);
	

	array_1d<double,6> fbd_stblterm = ZeroVector(matsize);
	fbd_stblterm = thawone *1.0* prod(trans(conv_opr),bdf);

	for(int ii = 0; ii< nodes_number; ++ii)
	  {
		int index = ii*(dof + 1);
		int loc_index = ii*dof;
		F[index] += 1*area*fbd_stblterm[loc_index];
		F[index + 1] += 1*area*fbd_stblterm[loc_index + 1];
	  }

	KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DASGS::CalculateGradStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double time,const double thawone,const double area)
	{
		KRATOS_TRY
	int nodes_number = 3;
	int dof = 2;
	
	double density;
	calculatedensity(GetGeometry(), density);

	//build 1*(grad q . grad p) stabilization term & assemble
	boost::numeric::ublas::bounded_matrix<double,3,3> gard_opr = ZeroMatrix(nodes_number,nodes_number);
	gard_opr = 1.0*thawone * prod(DN_DX,trans(DN_DX)); 

	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1) + dof;
		
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;

		K(row,column) += area *1*gard_opr(ii,jj);

		   }
	    }

	//build 1*thaw1*ro/deltat U grad q)
	double fac = thawone*density/time;
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;

			//K(row,column) += -1*area * fac* N(ii) * DN_DX(jj,0);
			K(column,row) += 1*area * fac* N(ii) * DN_DX(jj,0);

			//K(row + 1,column) += -1*area * fac* N(ii) * DN_DX(jj,1);
			K(column,row + 1) += 1*area * fac* N(ii) * DN_DX(jj,1);
		   }
	    }


	//build 1*(grad q) (Fbody + ro/dt * U(n+1,i)) stabilization term & assemble
	array_1d<double,2> bdf = ZeroVector(2);
	const array_1d<double,2> bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double,2> bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double,2> bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);


	const array_1d<double,3>& vel0_n = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
	const array_1d<double,3>& vel1_n = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
	const array_1d<double,3>& vel2_n = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);


	bdf[0] = N[0]*(density*bdf0[0] + density/time * vel0_n[0] ) +  N[1]*(density*bdf1[0] + density/time * vel1_n[0]) + N[2]*(density*bdf2[0] + density/time * vel2_n[0]);
	bdf[1] =  N[0]*(density*bdf0[1] + density/time * vel0_n[1] ) +  N[1]*(density*bdf1[1] + density/time * vel1_n[1]) + N[2]*(density*bdf2[1] + density/time * vel2_n[1]);

	array_1d<double,3> fbd_stblterm = ZeroVector(nodes_number);
	fbd_stblterm = thawone * prod(DN_DX,bdf);
	

	for(int ii = 0; ii< nodes_number; ++ii)
	  {
		int index = ii*(dof + 1) + dof;
		F[index] += 1*area*fbd_stblterm[ii];
	  }

	KRATOS_CATCH("")

	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DASGS::AddBodyForceAndMomentum(VectorType& F,const array_1d<double,3>& N, const double time,const double area)
	{
		KRATOS_TRY
	int nodes_number = 3;
	int dof = 2;
	

	double lump_mass_fac = area * 0.333333333333333333333333;

	double density;
	calculatedensity(GetGeometry(), density);

	//body  & momentum term force
	for ( int ii = 0; ii < nodes_number; ii++)
	   {
		int index = ii*(dof + 1) ;
		const array_1d<double,2> bdf = GetGeometry()[ii].FastGetSolutionStepValue(BODY_FORCE);
		const array_1d<double,2> ndvel = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY,1);

		
		
		F[index] += (area*N[ii]*density*bdf[0] + density/time*lump_mass_fac*ndvel[0]);
		F[index + 1] += (area*N[ii]*density*bdf[1] + density/time*lump_mass_fac*ndvel[1]);
	   }
	

	KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DASGS::CalculateResidual(const MatrixType& K, VectorType& F)
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
		UP[index + 2] = GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE,0);
	    }

	F -= prod(K,UP);

	KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DASGS::ComputeProjections(array_1d<double,6>& adv_proj , array_1d<double,3>& div_proj, const 			boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const double thawone,const double thawtwo,const array_1d<double,3>& N,const double area, const double time)
	{
	unsigned int number_of_nodes = GetGeometry().PointsNumber();
	unsigned int dim = 2;
	
	double density;
	calculatedensity(GetGeometry(), density);

	const array_1d<double,3>& adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double,3>& adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double,3>& adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);



	ms_adv_vel[0] = N[0]*(adv_vel0[0]-mesh_vel0[0])+N[1]*(adv_vel1[0]-mesh_vel1[0])+N[2]*(adv_vel2[0]-mesh_vel2[0]);
	ms_adv_vel[1] = N[0]*(adv_vel0[1]-mesh_vel0[1])+N[1]*(adv_vel1[1]-mesh_vel1[1])+N[2]*(adv_vel2[1]-mesh_vel2[1]);


	double const_adv_proj_X =0.0;
	double const_adv_proj_Y =0.0;
	double mean_div_proj =0.0;
	array_1d<double,3> mean_new_vel = ZeroVector(3);
	array_1d<double,3> mean_old_vel = ZeroVector(3);
	array_1d<double,3> mean_bdf = ZeroVector(3);


	for (unsigned int i=0;i<number_of_nodes;i++)
	 {
		
		//int index = i*dim;
		double pr = GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
		const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& old_vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
		const array_1d<double,3>& bdf = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);

		//const array_1d<double,2>& bdf = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
		// to consider the jump gradp/ro is calculated
			//pr = pr/density;

		//adv_proj = PI(ro*dv/dt + ro*a.gradU + gradP - f) considering lumped mass PI() = ()
		//calculate constant part of RES ->ro*a.gradU + gradP
		const_adv_proj_X += ( pr * DN_DX(i,0) + density*(ms_adv_vel[0]*DN_DX(i,0) + ms_adv_vel[1]*DN_DX(i,1))*vel[0] );
		const_adv_proj_Y +=  (pr * DN_DX(i,1) + density*(ms_adv_vel[0]*DN_DX(i,0) + ms_adv_vel[1]*DN_DX(i,1))*vel[1] );

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
		const double pr = GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
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
		advtermproj[0]+= thawone*adv_proj[index];
		advtermproj[1]+= thawone*adv_proj[index +1];

		double& divtermproj = GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ);
		divtermproj += thawtwo*div_proj[i] ;
		
		//calculate nodal area
		GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += 0.333333333333333333*area;
	 }*/

	
	}

	//************************************************************************************
	//************************************************************************************

	void Fluid2DASGS::AddProjectionForces(VectorType& F, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double area)
	{
	unsigned int number_of_nodes = GetGeometry().PointsNumber();
	unsigned int dim = 2;

	double density;
	calculatedensity(GetGeometry(), density);

	 const array_1d<double,3> advproj_0 = GetGeometry()[0].FastGetSolutionStepValue(ADVPROJ);
	 const array_1d<double,3> advproj_1 = GetGeometry()[1].FastGetSolutionStepValue(ADVPROJ);
	 const array_1d<double,3> advproj_2 = GetGeometry()[2].FastGetSolutionStepValue(ADVPROJ);

	 const double div_proj_0 = GetGeometry()[0].FastGetSolutionStepValue(DIVPROJ);
	 const double div_proj_1 = GetGeometry()[1].FastGetSolutionStepValue(DIVPROJ);
	 const double div_proj_2 = GetGeometry()[2].FastGetSolutionStepValue(DIVPROJ);


	
	//mean values
	double mean_x_adv = 0.3333333333333333*(advproj_0[0] +advproj_1[0] + advproj_2[0]); 
	double mean_y_adv = 0.3333333333333333*(advproj_0[1] +advproj_1[1] + advproj_2[1]); 
	
	double mean_div = 0.3333333333333333*(div_proj_0 +div_proj_1 + div_proj_2); 
	
	for (unsigned int ii=0;ii<number_of_nodes;ii++)
	 {
		int index = ii*(dim + 1) ;
		const array_1d<double,3>& vel = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY);

		const array_1d<double,3>& mesh_vel = GetGeometry()[ii].FastGetSolutionStepValue(MESH_VELOCITY);
		array_1d<double,2> adv_vel = ZeroVector(2);
		adv_vel[0] = vel[0]-mesh_vel[0];
		adv_vel[1] = vel[1]-mesh_vel[1];
	
		//thawone*ro*(xi,a.gradv)
		double proj;

		proj = mean_x_adv*(adv_vel[0]*DN_DX(ii,0) + adv_vel[1]*DN_DX(ii,1));
		F[index] += (m_thawone*1.0*area*proj);

		proj = mean_y_adv*(adv_vel[0]*DN_DX(ii,0) + adv_vel[1]*DN_DX(ii,1));
		F[index +1 ] += (m_thawone*1.0*area*proj);
	

		//thawone*(xi,gradq)
		proj = (mean_x_adv*DN_DX(ii,0) + mean_y_adv*DN_DX(ii,1));
		F[index + 2] += (m_thawone*area*proj);

		//thawtwo*(divv)

		F[index] += (m_thawtwo*area*mean_div*DN_DX(ii,0));
		F[index +1] += (m_thawtwo*area*mean_div*DN_DX(ii,1));
		
	}

	


	}

	//************************************************************************************
	//************************************************************************************
	void Fluid2DASGS::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
				rResult[i*node_size+2] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
			}
		KRATOS_CATCH("")
		
	}

	//************************************************************************************
	//************************************************************************************
	  void Fluid2DASGS::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
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
				ElementalDofList[i*node_size+2] = GetGeometry()[i].pGetDof(PRESSURE);
			}
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************

   void Fluid2DASGS::Calculate( const Variable<array_1d<double,3> >& rVariable, 
                                    array_1d<double,3>& Output, 
                                    const ProcessInfo& rCurrentProcessInfo)
	{

	array_1d<double,6> adv_proj = ZeroVector(6);
	array_1d<double,3> div_proj = ZeroVector(3);

	double delta_t= rCurrentProcessInfo[DELTA_TIME];

	//getting data for the given geometry
	double Area;
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);

	ComputeProjections(adv_proj, div_proj, DN_DX,m_thawone,m_thawtwo,N,Area, delta_t);

	}

	//************************************************************************************
	//************************************************************************************

	void Fluid2DASGS::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
	{
		
            if(rVariable==THAWONE)
            {
                for( unsigned int PointNumber = 0; 
                     PointNumber < 1; 
                     PointNumber++ )
                {
                    rValues[PointNumber] = m_thawone;
                }
            }
            if(rVariable==THAWTWO)
            {
                for( unsigned int PointNumber = 0; 
                     PointNumber < 1; PointNumber++ )
                {
		
                    rValues[PointNumber] = m_thawtwo;
                }
            }

        }
	//*************************************************************************************
	//*************************************************************************************
	void Fluid2DASGS::calculatedensity(Geometry< Node<3> > geom, double& density)
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
	
	const double rho0 = geom[0].FastGetSolutionStepValue(DENSITY);
	const double rho1 = geom[1].FastGetSolutionStepValue(DENSITY);
	const double rho2 = geom[2].FastGetSolutionStepValue(DENSITY);
	 density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );




	}
	//*************************************************************************************
	//*************************************************************************************
		void Fluid2DASGS::CalculateThaw(double& thawone, double& thawtwo, const double time,const double area,const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
	//calculate mean advective velocity and thaws
	const array_1d<double,3>& adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double,3>& adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double,3>& adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);



	ms_adv_vel[0] = N[0]*(adv_vel0[0]-mesh_vel0[0])+N[1]*(adv_vel1[0]-mesh_vel1[0])+N[2]*(adv_vel2[0]-mesh_vel2[0]);
	ms_adv_vel[1] = N[0]*(adv_vel0[1]-mesh_vel0[1])+N[1]*(adv_vel1[1]-mesh_vel1[1])+N[2]*(adv_vel2[1]-mesh_vel2[1]);

	//ms_adv_vel[0] = 0.0;
	//ms_adv_vel[1] = 0.0;
	                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

	double advvel_norm = ms_adv_vel[0]*ms_adv_vel[0]+ms_adv_vel[1]*ms_adv_vel[1];
	advvel_norm = sqrt(advvel_norm);
	
	double ele_length = 2.0*sqrt(area/3.00);
	
	double mu;
	const double mu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
	const double mu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
	const double mu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
	mu = 0.333333333333333333333333*(mu0 + mu1 + mu2);

	double density;
	calculatedensity(GetGeometry(), density);

	int dyn_st_switch = rCurrentProcessInfo[FRACTIONAL_STEP];
	
		if(dyn_st_switch)
		  {
		
			thawone = 1.0/(1.0/time + 4.0*mu/(ele_length*ele_length*density)+2.0*advvel_norm*1.0/ele_length);
		  }
		else
		 {
			
			thawone = 1.0/(0.0+ 4.0*mu/(ele_length*ele_length*density)+2.0*advvel_norm*1.0/ele_length);
		  }
		
	thawtwo = mu/density + 1.0*ele_length*advvel_norm/2.0;

	m_thawone = thawone;
	m_thawtwo = thawtwo;

	KRATOS_CATCH("")

	}
	//*************************************************************************************
	//*************************************************************************************
} // Namespace Kratos


