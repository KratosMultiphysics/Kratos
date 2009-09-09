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
#include "custom_elements/fluid_2d_split.h"
#include "utilities/math_utils.h"
#include "utilities/divide_elem_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
        namespace Fluid2DSplitauxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
        #pragma omp threadprivate(DN_DX)

        array_1d<double,3> N = ZeroVector(3); //dimension = number of nodes
        #pragma omp threadprivate(N)

        array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
        #pragma omp threadprivate(ms_adv_vel)

    }
    using  namespace Fluid2DSplitauxiliaries;


	//************************************************************************************
	//************************************************************************************
	Fluid2DSplit::Fluid2DSplit(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Fluid2DSplit::Fluid2DSplit(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
		
	}

	Element::Pointer Fluid2DSplit::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		
		KRATOS_TRY
		return Element::Pointer(new Fluid2DSplit(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	Fluid2DSplit::~Fluid2DSplit()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
	//The shape functions are calculated on the unique gauss point (their value is then 1/3. 
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);
	 

	 if(GetValue(IS_DIVIDED) == 1.0)
	 {
		boost::numeric::ublas::bounded_matrix<double,4,2> aux_gp = ZeroMatrix(4,2);
		array_1d<double,4> A_on_agp = ZeroVector(4);
		boost::numeric::ublas::bounded_matrix<double,4,3> N_on_agp = ZeroMatrix(4,3);
		array_1d<double,4> dist_on_agp = ZeroVector(4);

		DivideElemUtils::DivideElement_2D(GetGeometry(),  aux_gp, A_on_agp, N_on_agp, dist_on_agp);
 
		for(unsigned int i = 0 ; i< aux_gp.size1() ; i++)
		{
		  if (dist_on_agp(i) <= 0.0)
		  {
// KRATOS_WATCH(this->Id())
		        double Area = A_on_agp(i);
// KRATOS_WATCH(Area);

		      
		        for (unsigned int j = 0; j < N.size(); j++)
			     N[j] = N_on_agp(i,j);
// KRATOS_WATCH(N)
// 		        double tauone = 0.0;//delta_t;
// 		        double tautwo = 0.0;//delta_t;*/
		        double tauone, tautwo;
 		        CalculateTau(tauone, tautwo, delta_t, Area , rCurrentProcessInfo);


		        //add body force and momentum
		        AddBodyForceAndMomentum(rRightHandSideVector, N, delta_t, Area, tauone,tautwo);

		        //add projections
		        if(rCurrentProcessInfo[OSS_SWITCH] == 1.0)
			   AddProjectionForces(rRightHandSideVector,DN_DX,Area,tauone, tautwo);	
		  }

		}

	 }
	 else if(GetValue(IS_DIVIDED) == -1.0)
	   {

// 		double tauone = 0.0;//delta_t;
// 		double tautwo = 0.0;//delta_t;
		double tauone, tautwo;
		CalculateTau(tauone, tautwo, delta_t, Area, rCurrentProcessInfo);


		//add body force and momentum
		AddBodyForceAndMomentum(rRightHandSideVector, N, delta_t, Area, tauone,tautwo);

		//add projections
		if(rCurrentProcessInfo[OSS_SWITCH] == 1.0)
		    AddProjectionForces(rRightHandSideVector,DN_DX,Area,tauone, tautwo);	
	   }
	 

		KRATOS_CATCH("")
	}
	//***********************************************************************************++
	//**************************************************************************************++
	void Fluid2DSplit::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		MatrixType temp = Matrix();
		CalculateLocalSystem(temp, rRightHandSideVector,  rCurrentProcessInfo);

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateMassContribution(MatrixType& K,const double time,const double area)
	{
	KRATOS_TRY
	double lump_mass_fac = area * 0.333333333333333333333333;
	double density;
	double mu;
	double eps;
	CalculateDensity(GetGeometry(), density, mu, eps);

// KRATOS_WATCH(mu)
// KRATOS_WATCH(eps)


	int nodes_number = 3;
	int dof = 2;
	for ( int nd = 0; nd< nodes_number; nd++)
	    {
		int row = nd*(dof + 1);
		for( int jj=0; jj< dof; jj++)
			K(row + jj, row + jj) += density * lump_mass_fac;
	    }
	
		KRATOS_CATCH("")
	
	}

	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		//lumped
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int NumberOfNodes = GetGeometry().size();
		unsigned int MatSize = (dimension + 1) * NumberOfNodes;
		if(rMassMatrix.size1() != MatSize)
			rMassMatrix.resize(MatSize,MatSize,false);

		rMassMatrix = ZeroMatrix(MatSize,MatSize);
		double delta_t= rCurrentProcessInfo[DELTA_TIME];
		

		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);
	
	       //In the case of the mass matrix of a divided element I have to calculate the contribution of the fluid part to the shape functions of the nodes of the elements.
	       //And finally we lumpe the resulting matrix.
		if(GetValue(IS_DIVIDED) == 1.0)
		{	
		        
		        int nodes_number = 3;
		        int dof = 2;
		        int matsize = dof*nodes_number;
		        double density;
		        double mu;
		        double eps;
		        CalculateDensity(GetGeometry(), density, mu, eps);

		        boost::numeric::ublas::bounded_matrix<double,4,2> aux_gp = ZeroMatrix(4,2);
		        array_1d<double,4> A_on_agp = ZeroVector(4);
		        boost::numeric::ublas::bounded_matrix<double,4,3> N_on_agp = ZeroMatrix(4,3);
		        array_1d<double,4> dist_on_agp = ZeroVector(4);

		        DivideElemUtils::DivideElement_2D(GetGeometry(), aux_gp, A_on_agp, N_on_agp, dist_on_agp);
		        
		        boost::numeric::ublas::bounded_matrix<double,6,6> temp_MassMatr = ZeroMatrix(matsize,matsize);

		        for(unsigned int i = 0 ; i< aux_gp.size1() ; i++)
		        {
			 //if the virtual sub element is a fluid element
			 if (dist_on_agp(i) <= 0.0)
			 {
			       double Area = A_on_agp(i);
			       //shape functions of the element calculated on the auxiliary gauss points (gp of the sub elements)
			       for (unsigned int j = 0; j < N.size(); j++)
				    N[j] = N_on_agp(i,j);

			       boost::numeric::ublas::bounded_matrix<double,2,6> shape_func = ZeroMatrix(dof, matsize);

			       for (unsigned int ii = 0; ii< N.size(); ii++)
				  {
				        int column = ii*2;
				        shape_func(0,column) = N[ii]; 
				        shape_func(1,column + 1) = shape_func(0, column);
				  }
		      // |  N0(agp_i)* N0(agp_i)		0	N0(agp_i)* N1(agp_i)	0	N0(agp_i)* N2(agp_i)	0	  |	
		      // |  0			N0(agp_i)* N0(agp_i)	0	N0(agp_i)* N1(agp_i)	0	N0(agp_i)* N2(agp_i)|	
        //temp_MassMatrix=	|N1(agp_i)* N0(agp_i)		0	N1(agp_i)* N1(agp_i)	0	N1(agp_i)* N2(agp_i)	0	  | * density * Area	(sum on agp_i)
		      // |  0			N1(agp_i)* N0(agp_i)	0	N1(agp_i)* N1(agp_i)	0	N1(agp_i)* N2(agp_i)|
		      // |  N2(agp_i)* N0(agp_i)		0	N2(agp_i)* N1(agp_i)	0	N2(agp_i)* N2(agp_i)	0	  |	
		      // |  0			N2(agp_i)* N0(agp_i)	0	N2(agp_i)* N1(agp_i)	0	N2(agp_i)* N2(agp_i)|	
			       

			       noalias(temp_MassMatr) += Area * density * prod(trans(shape_func),shape_func);	    
 
			 }
		        }

		      
		        //construing the lumped mass matrix from the 6*6 tempMassMatrix;
		        for ( int ii = 0; ii < nodes_number; ii++)
		        {
			     int row = ii*(dof+1);
			     int loc_row = ii*dof;
			     for( int jj=0; jj < nodes_number; jj++)
			       {
// 				    int column = jj*(dof+1);
				    int loc_column = jj*dof;

				    rMassMatrix(row,row) += temp_MassMatr(loc_row,loc_column);
				    rMassMatrix(row + 1,row + 1) += temp_MassMatr(loc_row + 1,loc_column + 1);
			       }
		        }

		        for(unsigned int i = 0 ; i< aux_gp.size1() ; i++)
		        {
			   //if the virtual sub element is a fluid element
			   if (dist_on_agp(i) <= 0.0)
			   {
				double Area = A_on_agp(i);
				//shape functions of the element calculated on the auxiliary gauss points (gp of the sub elements)
				for (unsigned int j = 0; j < N.size(); j++)
				      N[j] = N_on_agp(i,j);

	 // 		        double tauone = 0.0;//delta_t;
	 // 		        double tautwo = 0.0;//delta_t;*/
				double tauone, tautwo;
//comment***************************************************************
				CalculateTau(tauone, tautwo, delta_t, Area , rCurrentProcessInfo);
				//add stablilization terms due to advective term (a)grad(V) * rho*Acc
				CalculateAdvMassStblTerms(rMassMatrix, DN_DX, N,tauone,Area);
				//add stablilization terms due to grad term grad(q) * rho*Acc
				CalculateGradMassStblTerms(rMassMatrix, DN_DX, tauone,Area);
// 				//add stabilization due to darcy term (rho*Acc, (Av+B|v|v)/rho)
// // 				/*CalculateDarcyMassStblTerms(rMassMatrix, DN_DX, tauone,Area);*/
//comment***************************************************************
			   }	
		        }
		}
		else if(GetValue(IS_DIVIDED) == -1.0)
		  {
		        //Calculate tau
        /*		double tauone = 0.0;//delta_t;
		        double tautwo = 0.0;// delta_t;*/
		        double tauone, tautwo;
		        CalculateTau(tauone, tautwo, delta_t, Area, rCurrentProcessInfo);
//comment***************************************************************
		        CalculateMassContribution(rMassMatrix,delta_t,Area); 
		        //add stablilization terms due to advective term (a)grad(V) * ro*Acce
		        CalculateAdvMassStblTerms(rMassMatrix, DN_DX, N,tauone,Area);
		        //add stablilization terms due to grad term grad(q) * ro*Acce
		        CalculateGradMassStblTerms(rMassMatrix, DN_DX, tauone,Area);
// 		        //add stabilization due to darcy term (rho*Acc, (Av+B|v|v)/rho)
// // 		        /*CalculateDarcyMassStblTerms(rMassMatrix, DN_DX, tauone,Area);*/
//comment***************************************************************

		  }
	 
//***end*/
		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateLocalVelocityContribution(MatrixType& rDampMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
	int nodes_number = 3;
	int dim = 2;
	unsigned int matsize = nodes_number*(dim+1);

	if(rDampMatrix.size1() != matsize)
			rDampMatrix.resize(matsize,matsize,false); //false says not to preserve existing storage!!


	noalias(rDampMatrix) = ZeroMatrix(matsize,matsize); 
	
	double delta_t= rCurrentProcessInfo[DELTA_TIME];
	
	
	
	//getting data for the given geometry
	double Area;
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);
	

	 if(GetValue(IS_DIVIDED) == 1.0)
	 {
		boost::numeric::ublas::bounded_matrix<double,4,2> aux_gp = ZeroMatrix(4,2);
		array_1d<double,4> A_on_agp = ZeroVector(4);
		boost::numeric::ublas::bounded_matrix<double,4,3> N_on_agp = ZeroMatrix(4,3);
		array_1d<double,4> dist_on_agp = ZeroVector(4);

		DivideElemUtils::DivideElement_2D(GetGeometry(), aux_gp, A_on_agp, N_on_agp, dist_on_agp);
		for(unsigned int i = 0 ; i< aux_gp.size1() ; i++)
		{
		  if (dist_on_agp(i) <= 0.0)
		  {

		        double Area = A_on_agp(i);
		      
		        for (unsigned int j = 0; j < N.size(); j++)
			     N[j] = N_on_agp(i,j);

		        //viscous term	
//comment *****************************************************************************
		        CalculateViscousTerm(rDampMatrix, DN_DX, Area);
//comment *****************************************************************************
		        //Advective term
// 		        double tauone = 0.0;//delta_t;
// 		        double tautwo = 0.0;//delta_t;
		        double tauone, tautwo;
  		        CalculateTau(tauone, tautwo, delta_t, Area, rCurrentProcessInfo);
//comment *****************************************************************************
		        CalculateAdvectiveTerm(rDampMatrix, DN_DX, tauone, tautwo, delta_t, Area);
//comment *****************************************************************************				
// 		        //calculate pressure term
		        CalculatePressureTerm(rDampMatrix, DN_DX, N, delta_t,Area);
//decomment*********************************************************************

		        //calculate Darcy term
		        CalculateDarcyTerm_SubElem(rDampMatrix, N, Area);
//decomment*********************************************************************

		        //compute projections

		        //stabilization terms
		        CalculateDivStblTerm(rDampMatrix, DN_DX, tautwo, Area);
		        CalculateAdvStblAllTerms(rDampMatrix,rRightHandSideVector, DN_DX, N, tauone,delta_t, Area);
		        CalculateGradStblAllTerms(rDampMatrix,rRightHandSideVector,DN_DX, delta_t, tauone, Area);
//  /*		        CalculateDarcyStblAllTerms(rDampMatrix,rRightHandSideVector,DN_DX, delta_t, tauone, Area);*/
// 		        //KRATOS_WATCH(rRightHandSideVector);

		  }

		}
	 }
	 else if(GetValue(IS_DIVIDED) == -1.0)
	 {
/*KRATOS_WATCH("K mass_matrix")
KRATOS_WATCH(rDampMatrix)	*/		  //viscous term	
//comment *****************************************************************************
 		  CalculateViscousTerm(rDampMatrix, DN_DX, Area);
//comment *****************************************************************************
/*KRATOS_WATCH("K viscous")
KRATOS_WATCH(rDampMatrix)	*/		  
		  //Advective term
		  //		  double tauone =0.0; //delta_t;
		  //		  double tautwo =0.0;// delta_t;
		  double tauone, tautwo;
 		  CalculateTau(tauone, tautwo, delta_t, Area, rCurrentProcessInfo);
//comment *****************************************************************************
 		  CalculateAdvectiveTerm(rDampMatrix, DN_DX, tauone, tautwo, delta_t, Area);
// //comment *****************************************************************************
/*KRATOS_WATCH("K advection")
KRATOS_WATCH(rDampMatrix)	*/		   

		  //calculate pressure term
		  CalculatePressureTerm(rDampMatrix, DN_DX, N, delta_t,Area);
// KRATOS_WATCH("K pressure")
// KRATOS_WATCH(rDampMatrix)
//decomment*********************************************************************
		  //calculate Darcy term
		  CalculateDarcyTerm(rDampMatrix, Area);
//decomment*********************************************************************
// KRATOS_WATCH("K darcy****************************************************")
// KRATOS_WATCH(rDampMatrix)
		  //compute projections

		  //stabilization terms
//comment *******************************************************************************
		  CalculateDivStblTerm(rDampMatrix, DN_DX, tautwo, Area);
		  CalculateAdvStblAllTerms(rDampMatrix,rRightHandSideVector, DN_DX, N, tauone,delta_t, Area);
		  CalculateGradStblAllTerms(rDampMatrix,rRightHandSideVector,DN_DX, delta_t, tauone, Area);
// 		  /*CalculateDarcyStblAllTerms(rDampMatrix,rRightHandSideVector,DN_DX, delta_t, tauone, Area);*/
//comment *******************************************************************************
		  //KRATOS_WATCH(rRightHandSideVector);
	
	 }
//delete *******************
//calculating 



      


//delete *******************


		KRATOS_CATCH("")
	}
        


	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateViscousTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double area)
	{
		KRATOS_TRY
	double mu;
	/*const double mu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
	const double mu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
	const double mu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
	mu = 0.333333333333333333333333*(mu0 + mu1 + mu2);*/
	
	double density;
	double eps;
	CalculateDensity(GetGeometry(), density, mu, eps);
//KRATOS_WATCH(mu);

	//nu = mu/density;	

	int nodes_number = 3;
	int dof = 2;

	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);
		K(row,column) += mu * area*(DN_DX(ii,0)*DN_DX(jj,0) + DN_DX(ii,1)*DN_DX(jj,1));
		K(row + 1,column + 1) += mu * area*(DN_DX(ii,0)*DN_DX(jj,0) + DN_DX(ii,1)*DN_DX(jj,1));
		   }	
	    }
					

		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateAdvectiveTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double tauone, const double tautwo, const double time,const double area)
	{
		KRATOS_TRY
	//calculate mean advective velocity and taus
	const array_1d<double,3>& adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double,3>& adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double,3>& adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,0);
	const array_1d<double,3>& mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);

	double density;
	double mu;
	double eps;
	CalculateDensity(GetGeometry(), density, mu, eps);

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

	//double fac = tauone/time;
	//fac = 0.0; 
	//temp_convterm *= ((1 + fac*density)*density); // For the simplicity of implementation the stabilization term tau1*rho/deltat*(a.gradV U(n+1,t)) is added 	
	
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		int loc_row = ii*dof;
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);
			int loc_column = jj*dof;

			K(row,column) += area * density * temp_convterm(loc_row,loc_column);
			K(row + 1,column + 1) += area * density * temp_convterm(loc_row + 1,loc_column + 1);
		   }
	    }

		KRATOS_CATCH("")

	}
	//************************************************************************************
	//************************************************************************************
	 void Fluid2DSplit::CalculateDarcyTerm_SubElem(MatrixType& K, const array_1d<double,3>&  N, const double area)
	 {
		KRATOS_TRY

		double density;
		double mu;
		double eps;
		CalculateDensity(GetGeometry(), density, mu, eps);
// KRATOS_WATCH(mu);
		double dp = 0.01; //diameter of the particle	
		double kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*dp*dp);
// KRATOS_WATCH(kinv);
		const array_1d<double,3>& vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,0);
		const array_1d<double,3>& vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,0);
		const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,0);
// KRATOS_WATCH(vel0);	
// KRATOS_WATCH(vel1);
// KRATOS_WATCH(vel2);
		int nodes_number = 3;
		int dof = 2;
		int matsize = dof*nodes_number;

		array_1d<double,3> norm_vel_2 = ZeroVector(nodes_number);
		//	vector with the norm^2 of the velocity of the three nodes at the previous iteration;
		for (int ii = 0; ii < nodes_number; ii++)
		{
			 norm_vel_2[0] += vel0[ii]*vel0[ii];
			 norm_vel_2[1] += vel1[ii]*vel1[ii];
			 norm_vel_2[2] += vel2[ii]*vel2[ii];
		}
/*
		      // |  N0(agp_i)* N0(agp_i)		0	N0(agp_i)* N1(agp_i)	0	N0(agp_i)* N2(agp_i)	0	  |	
		      // |  0			N0(agp_i)* N0(agp_i)	0	N0(agp_i)* N1(agp_i)	0	N0(agp_i)* N2(agp_i)|	
        //temp_sfprod=	|N1(agp_i)* N0(agp_i)		0	N1(agp_i)* N1(agp_i)	0	N1(agp_i)* N2(agp_i)	0	  | * Darcy term	
		      // |  0			N1(agp_i)* N0(agp_i)	0	N1(agp_i)* N1(agp_i)	0	N1(agp_i)* N2(agp_i)|
		      // |  N2(agp_i)* N0(agp_i)		0	N2(agp_i)* N1(agp_i)	0	N2(agp_i)* N2(agp_i)	0	  |	
		      // |  0			N2(agp_i)* N0(agp_i)	0	N2(agp_i)* N1(agp_i)	0	N2(agp_i)* N2(agp_i)|	
			       
 
*/

		boost::numeric::ublas::bounded_matrix<double,2,6> shape_func = ZeroMatrix(dof, matsize);

		for (int ii = 0; ii< nodes_number; ii++)
		    {
			 int column = ii*dof;
			 shape_func(0,column) = N[ii];
			 shape_func(1,column + 1) = shape_func(0,column);
		    }
// KRATOS_WATCH(shape_func)
		boost::numeric::ublas::bounded_matrix<double,6,6> temp_sfprod = ZeroMatrix(matsize,matsize);
		noalias(temp_sfprod) = area * prod(trans(shape_func),shape_func);	    
// 		KRATOS_WATCH(temp_sfprod)

// Consistent form
// 		for ( int ii = 0; ii < nodes_number; ii++)
// 		    {
// 			 int row = ii*(dof+1);
// 			 int loc_row = ii*dof;
// 			 for( int jj=0; jj < nodes_number; jj++)
// 			   {
// 				int column = jj*(dof+1);
// 				int loc_column = jj*dof;
// 				//DARCY TERM linear part
// 				K(row,column) += (kinv * mu) * temp_sfprod(loc_row,loc_column);
// /*double aux = (kinv * mu) * temp_sfprod(loc_row,loc_column);
// KRATOS_WATCH(aux)*/
// 				K(row + 1,column + 1) += (kinv * mu) * temp_sfprod(loc_row + 1,loc_column + 1);
// /*double aux2 = (kinv * mu) * temp_sfprod(loc_row + 1,loc_column + 1);
// KRATOS_WATCH(aux2)*/
// 				//DARCY TERM nonlinear part
// 				K(row,column) += (1.75 * density /eps * sqrt(norm_vel_2[jj] *  kinv / (eps * 150.0))) * temp_sfprod(loc_row,loc_column);
// 				K(row + 1,column + 1) += (1.75 * density /eps * sqrt(norm_vel_2[jj] *  kinv / (eps * 150.0))) * temp_sfprod(loc_row + 1,loc_column + 1);
// // KRATOS_WATCH(K(row,column))
// 
// 			   }
// 		    }
// KRATOS_WATCH("after")
// 	KRATOS_WATCH(mu)
// Lumped form
		for ( int ii = 0; ii < nodes_number; ii++)
		    {
			 int row = ii*(dof+1);
			 int loc_row = ii*dof;
			 for( int jj=0; jj < nodes_number; jj++)
			   {
/*				int column = jj*(dof+1);*/
				int loc_column = jj*dof;
				//DARCY TERM linear part
				K(row,row) += (kinv * mu) * temp_sfprod(loc_row,loc_column);
/*double aux = (kinv * mu) * temp_sfprod(loc_row,loc_column);
KRATOS_WATCH(aux)*/
				K(row + 1,row + 1) += (kinv * mu) * temp_sfprod(loc_row + 1,loc_column + 1);
/*double aux2 = (kinv * mu) * temp_sfprod(loc_row + 1,loc_column + 1);
KRATOS_WATCH(aux2)*/
				//DARCY TERM nonlinear part
				K(row,row) += (1.75 * density /eps * sqrt(norm_vel_2[jj] *  kinv / (eps * 150.0))) * temp_sfprod(loc_row,loc_column);
				K(row + 1,row + 1) += (1.75 * density /eps * sqrt(norm_vel_2[jj] *  kinv / (eps * 150.0))) * temp_sfprod(loc_row + 1,loc_column + 1);
			   }
// KRATOS_WATCH(K(row,column))
		    }

		KRATOS_CATCH("")
	  }

	 void Fluid2DSplit::CalculateDarcyTerm(MatrixType& K, const double area)
	 {
		KRATOS_TRY
		double lump_mass_fac = area * 0.333333333333333333333333;
		double density;
		double mu;
		double eps;
		CalculateDensity(GetGeometry(), density, mu, eps);
		
		double dp = 0.01; //diameter of the particle	
		double kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*dp*dp);
// KRATOS_WATCH(kinv);
// KRATOS_WATCH(mu);
// KRATOS_WATCH(lump_mass_fac);
// KRATOS_WATCH(density);
// KRATOS_WATCH(eps);

		const array_1d<double,3>& vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,0);
		const array_1d<double,3>& vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,0);
		const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,0);
		
		int nodes_number = 3;
		int dof = 2;
// 		int matsize = dof*nodes_number;

		array_1d<double,3> norm_vel_2 = ZeroVector(nodes_number);
		//	vector with the norm^2 of the velocity of the three nodes at the previous iteration;
		for (int ii = 0; ii < nodes_number; ii++)
		{
			 norm_vel_2[0] += vel0[ii]*vel0[ii];
			 norm_vel_2[1] += vel1[ii]*vel1[ii];
			 norm_vel_2[2] += vel2[ii]*vel2[ii];
		}
//delete *******************************************************************
// 	boost::numeric::ublas::bounded_matrix<double,2,2> temporary_DARCYterm = ZeroMatrix(2,2);
// // 	array_1d<double,2> AUX_VEL = ZeroVector(2);

//delete *******************************************************************
		for ( int nd = 0; nd< nodes_number; nd++)
		{
		      int row = nd*(dof + 1);
		      for( int jj=0; jj< dof; jj++)
		      {	   //DARCY TERM linear part:
			   K(row + jj, row + jj) +=  (kinv * mu) * lump_mass_fac;
			   //DARCY TERM nonlinear part:
			   K(row + jj, row + jj) += (1.75 * density /eps * sqrt(norm_vel_2[nd] *  kinv / (eps * 150.0)))* lump_mass_fac;
		      
//delete *******************************************************************
// 			   temporary_DARCYterm(jj, jj) =  (kinv * mu) * lump_mass_fac;
// // 			   //DARCY TERM nonlinear part:
// 			   temporary_DARCYterm(jj, jj) += (1.75 * density /eps * sqrt(norm_vel_2[nd] *  kinv / (eps * 150.0)))* lump_mass_fac;
//delete *******************************************************************
		      }
//delete *******************************************************************
// AUX_VEL(0)=GetGeometry()[nd].FastGetSolutionStepValue(VELOCITY_X,0);
// AUX_VEL(1)=GetGeometry()[nd].FastGetSolutionStepValue(VELOCITY_Y,0);
// // 
// KRATOS_WATCH(AUX_VEL);
// // KRATOS_WATCH(nd);
// KRATOS_WATCH(prod(temporary_DARCYterm, AUX_VEL));
/*			   AUX_VEL(row) = GetGeometry()[nd].FastGetSolutionStepValue(VELOCITY_X,0);
			   AUX_VEL(row + 1) = GetGeometry()[nd].FastGetSolutionStepValue(VELOCITY_Y,0);*/
//delete *******************************************************************
		}
//delete *******************************************************************
// KRATOS_WATCH(AUX_VEL);
// KRATOS_WATCH(temporary_DARCYterm);
//delete *******************************************************************

// KRATOS_WATCH(area)
// KRATOS_WATCH(kinv)
// KRATOS_WATCH(lump_mass_fac)
// 		KRATOS_WATCH((kinv * mu) * lump_mass_fac)
//delete *******************************************************************
// 		KRATOS_WATCH(prod(temp_DARCYterm,AUX_VEL));
//delete *******************************************************************
		KRATOS_CATCH("")
	 }
	//************************************************************************************
	//************************************************************************************
	//Calculate the divergence and the gradient operators
	void Fluid2DSplit::CalculatePressureTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>&  N, const double time,const double area)
	{
		KRATOS_TRY
	int nodes_number = 3;
	int dof = 2;

	double density;
	double mu;
	double eps;
	CalculateDensity(GetGeometry(), density, mu, eps);
//KRATOS_WATCH(this->Id())
//KRATOS_WATCH(density)

        //Pj mean point of the edge ik
        //N_on_meanp(i,j) = shape function od the node i calculated on the point Pj
	boost::numeric::ublas::bounded_matrix<double,3,3> N_on_meanp = ZeroMatrix(3,3);
	N_on_meanp(0,0) = 0.0; N_on_meanp(0,1) = 0.5; N_on_meanp(0,2) = 0.5;
	N_on_meanp(1,0) = 0.5; N_on_meanp(1,1) = 0.0; N_on_meanp(1,2) = 0.5;
	N_on_meanp(2,0) = 0.5; N_on_meanp(2,1) = 0.5; N_on_meanp(2,2) = 0.0;


	// the matrix l_n contain the product of edge lengh l_jk*n_i (1st column x-comp, 2nd column y-comp, i = row index).
	boost::numeric::ublas::bounded_matrix<double,3,2> l_n = ZeroMatrix(3,2);
	l_n(0,0) = -2.0 * area * DN_DX(0,0); l_n(0,1) = -2.0 * area * DN_DX(0,1);
	l_n(1,0) = -2.0 * area * DN_DX(1,0); l_n(1,1) = -2.0 * area * DN_DX(1,1); 
	l_n(2,0) = -2.0 * area * DN_DX(2,0); l_n(2,1) = -2.0 * area * DN_DX(2,1);   

// 	  const double& rho0 = 	GetGeometry()[0].FastGetSolutionStepValue(DENSITY,0);
// 	  const double& rho1 =	GetGeometry()[1].FastGetSolutionStepValue(DENSITY,0);
// 	  const double& rho2 =	GetGeometry()[2].FastGetSolutionStepValue(DENSITY,0);
// 	  const double& eps0 = 	GetGeometry()[0].FastGetSolutionStepValue(POROSITY,0);
// 	  const double& eps1 =	GetGeometry()[1].FastGetSolutionStepValue(POROSITY,0);
// 	  const double& eps2 =	GetGeometry()[2].FastGetSolutionStepValue(POROSITY,0);
// 	  array_1d<double,3> rho = ZeroVector(3);
// 	 rho[0]= rho0 ; //* eps0;
// 	 rho[1] = rho1 ; //* eps1;
// 	 rho[2] = rho2 ; //* eps2;



// KRATOS_WATCH(l_n);
//KRATOS_WATCH(rho);
//delete *******************************************************************
// 	boost::numeric::ublas::bounded_matrix<double,9,9> temp_PRESSterm = ZeroMatrix(9,9);
// 	array_1d<double,9> AUX_PRES = ZeroVector(9);
// array_1d<double,3> AUX_PRES = ZeroVector(3);
//delete *******************************************************************


	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;
			//**************************************************************
 			//Elemental gradient of pressure term (momentum equation)
			K(row,column) -= area * N(jj) * DN_DX(ii,0);
			K(row + 1,column) -= area * N(jj) * DN_DX(ii,1);
// 			 double prova = 0.0;
// 			 double prova2 = 0.0;

			//**************************************************************
// 			//Elemental divergence terms (continuity equation)

			 //Formulation n2 calculate int_Area(q * Div( rho * u ) ) = - int_Area(Grad q * (rho*u)  + int_bound(rho * v dot n *q)
			 //in order to take into account the jump in the density value. NOT EFFECTIVE ----> exactly equal to Formulation n1
// 			for(int kk = 0; kk < nodes_number; kk++)
// 			{
// 				//	//int_bound(rho * v dot n *q)
// 				K(column,row) +=  density*N_on_meanp(jj,kk)*l_n(kk,0)*N_on_meanp(ii,kk);
// 				K(column,row + 1) +=  density*N_on_meanp(jj,kk)*l_n(kk,1)*N_on_meanp(ii,kk);
// // // 				prova += density*N_on_meanp(jj,kk)*l_n(kk,0)*N_on_meanp(ii,kk);
// // // 				prova2  +=  density*N_on_meanp(jj,kk)*l_n(kk,1)*N_on_meanp(ii,kk);
// 			}
// KRATOS_WATCH("K(column, row)")	KRATOS_WATCH(column)  	 KRATOS_WATCH(row)
// 				KRATOS_WATCH(prova);
// KRATOS_WATCH(row+1)
// KRATOS_WATCH(prova2)
//			//- int_Area(Grad q * (rho*u) 
//			K(column,row) -= area*density*DN_DX(jj,0)*N(ii);
//			K(column,row + 1) -= area*density*DN_DX(jj,1)*N(ii);
// KRATOS_WATCH("prova_2parte")
// KRATOS_WATCH(row)
// KRATOS_WATCH(area*density*(jj,0)*N(ii));
// KRATOS_WATCH(row+1)
// KRATOS_WATCH(area*density*DN_DX(jj,1)*N(ii));
// KRATOS_WATCH("old version ****************************************************************************")
// KRATOS_WATCH(area *density* N(jj) * DN_DX(ii,0));
// KRATOS_WATCH(row+1)
// KRATOS_WATCH(area *density* N(jj) * DN_DX(ii,1));
// 



// 	 		// Fomulation n1  int( q * rho * Div( u ))
			K(column,row) += area *density * N(jj) * DN_DX(ii,0);
			K(column,row + 1) += area *density * N(jj) * DN_DX(ii,1);

//delete *******************************************************************
// 			temp_PRESSterm(row,column) -= area * N(jj) * DN_DX(ii,0);
// 			temp_PRESSterm(row+1,column) -= area * N(jj) * DN_DX(ii,1);
// 			AUX_PRES(column) = GetGeometry()[jj].FastGetSolutionStepValue(PRESSURE,0);
// 			AUX_PRES(jj) = GetGeometry()[jj].FastGetSolutionStepValue(PRESSURE,0);

//delete *******************************************************************


// 			// Fomulation n3  Ni * DN_DX(j,k)* (rho_i * u_i)  ******* int(q * Div( rho * u ) ) (RHO NODAL!!!!!!!!!!)
// 			K(column,row) += area * N(jj) * DN_DX(ii,0) * rho(ii);
// 			K(column,row + 1) += area * N(jj) * DN_DX(ii,1) * rho (ii);
 



		   }
	    }
//delete *******************************************************************
// KRATOS_WATCH(temp_PRESSterm);
//KRATOS_WATCH(AUX_PRES);
// KRATOS_WATCH(prod(trans(DN_DX),AUX_PRES))
//delete *******************************************************************

// KRATOS_WATCH(K);
//delete *******************************************************************
// 			KRATOS_WATCH(prod(temp_PRESSterm,AUX_PRES));
//delete *******************************************************************

	    KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateDivStblTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double tautwo,const double area)
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
	double mu;
	double eps;
	CalculateDensity(GetGeometry(), density, mu, eps);


	boost::numeric::ublas::bounded_matrix<double,6,6> temp_div = ZeroMatrix(matsize,matsize);
	temp_div = tautwo * prod(trans(div_opr),div_opr);
//delete *******************************************************************
// 	boost::numeric::ublas::bounded_matrix<double,9,9> temp_divSTBLterm = ZeroMatrix(9,9);
//delete *******************************************************************
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		int loc_row = ii*dof;
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);
			int loc_column = jj*dof;

			K(row,column) += area*density*temp_div(loc_row,loc_column);
			K(row,column + 1) += area*density*temp_div(loc_row,loc_column + 1);
			K(row + 1,column) += area*density*temp_div(loc_row + 1,loc_column);
			K(row + 1,column + 1) += area*density*temp_div(loc_row + 1,loc_column + 1);
//delete *******************************************************************
// 			temp_divSTBLterm(row,column) +=area*density*temp_div(loc_row,loc_column);
// 			temp_divSTBLterm(row,column + 1) += area*density*temp_div(loc_row,loc_column + 1);
// 			temp_divSTBLterm(row + 1,column) += area*density*temp_div(loc_row + 1,loc_column);
// 			temp_divSTBLterm(row + 1,column + 1) += area*density*temp_div(loc_row + 1,loc_column + 1);
//delete *******************************************************************
		   }
	    }
//delete *******************************************************************
// 			KRATOS_WATCH(temp_divSTBLterm);
//delete *******************************************************************
		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateAdvStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N, const double tauone,const double time,const double area)
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
	adv_stblterm = tauone * prod(trans(conv_opr),conv_opr);


	double density;
	double mu;
	double eps;
	CalculateDensity(GetGeometry(), density, mu, eps);

	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		int loc_row = ii*dof;
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);
			int loc_column = jj*dof;
			
			K(row,column) += area*density*adv_stblterm(loc_row,loc_column);
			K(row,column + 1) += area*density*adv_stblterm(loc_row,loc_column + 1);
			K(row + 1,column) += area*density*adv_stblterm(loc_row + 1,loc_column);
			K(row + 1,column + 1) += area*density*adv_stblterm(loc_row + 1,loc_column + 1);
		   }
	    }
//decomment*********************************************************************
	
	 //build Darcy Au+B|u|u * (a.grad v) stabilization term & assemble + build (rho a.grad u)* (Av+B|v|v)/rho   stabilization term & assemble
         double dp = 0.01; //diameter of the particle	
         double kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*dp*dp);

// 	 const array_1d<double,3>& vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,0);
// 	 const array_1d<double,3>& vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,0);
// 	 const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,0);

	 array_1d<double,3> norm_vel_2 = ZeroVector(nodes_number);
	 //	vector with the norm^2 of the velocity of the three nodes at the previous iteration;
	 for (int ii = 0; ii < nodes_number; ii++)
	 {
		  norm_vel_2[0] += adv_vel0[ii]*adv_vel0[ii];
		  norm_vel_2[1] += adv_vel1[ii]*adv_vel1[ii];
		  norm_vel_2[2] += adv_vel2[ii]*adv_vel2[ii];
	 }

	boost::numeric::ublas::bounded_matrix<double,2,6> darcy_opr = ZeroMatrix(dof,matsize);

	for (int ii = 0; ii< nodes_number; ii++)
	    {
		int column = ii*dof;
					
		//DARCY TERM linear part /rho
		darcy_opr(0,column) = N[ii] * kinv * mu / density;
		darcy_opr(1,column + 1) = darcy_opr(0,column);
		//DARCY TERM nonlinear part /rho
		darcy_opr(0,column) += N[ii] * 1.75  /eps * sqrt(norm_vel_2[ii] *  kinv / (eps * 150.0));
		darcy_opr(1,column + 1) += N[ii] * 1.75  /eps * sqrt(norm_vel_2[ii] *  kinv / (eps * 150.0));

	    }
//delete *******************************************************************
// 	boost::numeric::ublas::bounded_matrix<double,9,9> temp_pressSTBLterm = ZeroMatrix(9,9);
// 	boost::numeric::ublas::bounded_matrix<double,9,9> temp_darcySTBLterm = ZeroMatrix(9,9);
// 	array_1d<double,9> AUX_PRES = ZeroVector(9);
// 	array_1d<double,9> AUX_VEL = ZeroVector(9);
//delete *******************************************************************


	boost::numeric::ublas::bounded_matrix<double,6,6> darcy_stblterm = ZeroMatrix(matsize,matsize);		
	darcy_stblterm = tauone * prod(trans(conv_opr),darcy_opr);

	 for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		int loc_row = ii*dof;
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);
			int loc_column = jj*dof;
			//Au+B|u|u * (a.grad v)
			K(row,column) += area*density *darcy_stblterm(loc_row,loc_column);
			K(row,column + 1) += area*density *darcy_stblterm(loc_row,loc_column + 1);
			K(row + 1,column) += area*density *darcy_stblterm(loc_row + 1,loc_column);
			K(row + 1,column + 1) += area*density *darcy_stblterm(loc_row + 1,loc_column + 1);

//delete *******************************************************************
// 			temp_darcySTBLterm(row,column) += area*density *darcy_stblterm(loc_row,loc_column);
// 			temp_darcySTBLterm(row,column + 1) += area*density *darcy_stblterm(loc_row,loc_column + 1);
// 			temp_darcySTBLterm(row + 1,column) += area*density *darcy_stblterm(loc_row + 1,loc_column);
// 			temp_darcySTBLterm(row + 1,column + 1) += area*density *darcy_stblterm(loc_row + 1,loc_column + 1);
// 			AUX_VEL(column) = GetGeometry()[jj].FastGetSolutionStepValue(VELOCITY_X,0);
// 			AUX_VEL(column+1) = GetGeometry()[jj].FastGetSolutionStepValue(VELOCITY_Y,0);
//delete *******************************************************************


// 			//(rho a.grad u)* (Av+B|v|v)/rho
// 			K(column, row) += area*density *darcy_stblterm(loc_row,loc_column);
// 			K(column + 1, row) += area*density *darcy_stblterm(loc_row,loc_column + 1);
// 			K(column, row + 1) += area*density *darcy_stblterm(loc_row + 1,loc_column);
// 			K(column + 1, row + 1) += area*density *darcy_stblterm(loc_row + 1,loc_column + 1);

		   }
	    }
//decomment*********************************************************************


	//build 1*tau1*(a.grad V)(grad P) & 1*tau1*(grad q)(ro*a.grad U) stabilization terms & assemble
	boost::numeric::ublas::bounded_matrix<double,6,3> grad_stblterm = ZeroMatrix(matsize,nodes_number);
	grad_stblterm = tauone * prod(trans(conv_opr),trans(DN_DX));

	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		int loc_row = ii*dof;
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;
		        //1*tau1*(a.grad V)(grad P)
			K(row,column) += area * grad_stblterm(loc_row,jj);
			K(row + 1,column) += area * grad_stblterm(loc_row + 1, jj);

//delete *******************************************************************
/*
			temp_pressSTBLterm(row,column) += area * grad_stblterm(loc_row,jj);
			temp_pressSTBLterm(row + 1,column) += area * grad_stblterm(loc_row + 1, jj);
			AUX_PRES(column) = GetGeometry()[jj].FastGetSolutionStepValue(PRESSURE,0);*/
//delete *******************************************************************

			//1*tau1*(grad q)(ro*a.grad U) 
			K(column, row) += area * density*grad_stblterm(loc_row, jj);
			K(column, row + 1) += area *density* grad_stblterm(loc_row + 1, jj);
		   }
	    }

//delete *******************************************************************
// KRATOS_WATCH(this->Id())
// KRATOS_WATCH(prod(temp_darcySTBLterm,AUX_VEL));
// KRATOS_WATCH(prod(temp_pressSTBLterm, AUX_PRES));
//delete ******************************************************************* 

	//build (1.0*a.grad V) (Fbody) stabilization term & assemble
	array_1d<double,2> bdf = ZeroVector(2);
	const array_1d<double,2> bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double,2> bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double,2> bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);


	bdf[0] = N[0]*(density*bdf0[0] ) +  N[1]*(density*bdf1[0]) + N[2]*(density*bdf2[0]);
	bdf[1] =  N[0]*(density*bdf0[1]) +  N[1]*(density*bdf1[1] ) + N[2]*(density*bdf2[1] );
	

	array_1d<double,6> fbd_stblterm = ZeroVector(matsize);
	fbd_stblterm = tauone * prod(trans(conv_opr),bdf);

	for(int ii = 0; ii< nodes_number; ++ii)
	  {
		int index = ii*(dof + 1);
		int loc_index = ii*dof;
		F[index] += area*fbd_stblterm[loc_index];
		F[index + 1] += area*fbd_stblterm[loc_index + 1];
	  }
	KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateAdvMassStblTerms(MatrixType& M,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N, const double tauone,const double area)
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

	//calculate density
	double density;
	double mu;
	double eps;
	CalculateDensity(GetGeometry(), density, mu, eps);

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


	//tau1*rho*Nacc.(1.0*a.grad V)
	boost::numeric::ublas::bounded_matrix<double,6,6> temp_convterm = ZeroMatrix(matsize,matsize);
	temp_convterm = prod(trans(conv_opr),trans(shape_func));

	double fac = tauone*density;	
	
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		int loc_row = ii*dof;
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);
			int loc_column = jj*dof;

			M(row,column) += area*fac*temp_convterm(loc_row,loc_column);
			M(row + 1,column + 1) += area*fac*temp_convterm(loc_row + 1,loc_column + 1);
		   }
	    }


	KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateGradStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double time,const double tauone,const double area)
	{
		KRATOS_TRY
	int nodes_number = 3;
	int dof = 2;
	
	double density;
	double mu;
	double eps;
	CalculateDensity(GetGeometry(), density, mu, eps);

	//build 1*(grad q . grad p) stabilization term & assemble
	boost::numeric::ublas::bounded_matrix<double,3,3> gard_opr = ZeroMatrix(nodes_number,nodes_number);
	gard_opr = tauone * prod(DN_DX,trans(DN_DX)); 


	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1) + dof;
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;
			K(row,column) += area *gard_opr(ii,jj);
		   }
	    }
//decomment*********************************************************************
	//build Darcy Au+B|u|u * grad q stabilization term and assemble and grad p * Av+B|v|v and assemble
         double dp = 0.01; //diameter of the particle	
         double kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*dp*dp);

	 const array_1d<double,3>& vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,0);
	 const array_1d<double,3>& vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,0);
	 const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,0);

	 array_1d<double,3> norm_vel_2 = ZeroVector(nodes_number);
	 //	vector with the norm^2 of the velocity of the three nodes at the previous iteration;
	 for (int ii = 0; ii < nodes_number; ii++)
	 {
		  norm_vel_2[0] += vel0[ii]*vel0[ii];
		  norm_vel_2[1] += vel1[ii]*vel1[ii];
		  norm_vel_2[2] += vel2[ii]*vel2[ii];
	 }

	double fac_1 = tauone * kinv * mu / density;

	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		double fac_2 = tauone * 1.75 /eps * sqrt(norm_vel_2[ii] *  kinv / (eps * 150.0));

		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;
			//Au+B|u|u * grad q
			//DARCY TERM linear part /rho
			K(column,row) += area * density * fac_1* N[ii] * DN_DX(jj,0);
			K(column,row + 1) += area * density * fac_1* N[ii] * DN_DX(jj,1);
			//DARCY TERM nonlinear part /rho
			K(column,row) += area * density * fac_2* N[ii] * DN_DX(jj,0);
			K(column,row + 1) += area * density * fac_2* N[ii] * DN_DX(jj,1);


		   }
	    }
//decomment*********************************************************************

//delete *******************************************************************
// KRATOS_WATCH(prod(temp_darcySTBLterm,AUX_VEL));
// KRATOS_WATCH(prod(temp_pressSTBLterm, AUX_PRES));
//delete ******************************************************************* 
//COMMENT*********************************************************************************     

	//build 1*(grad q) (Fbody ) stabilization term & assemble
	array_1d<double,2> bdf = ZeroVector(2);
	const array_1d<double,2> bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double,2> bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double,2> bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);


	bdf[0] = N[0]*(density*bdf0[0] ) +  N[1]*(density*bdf1[0] ) + N[2]*(density*bdf2[0]);
	bdf[1] =  N[0]*(density*bdf0[1] ) +  N[1]*(density*bdf1[1] ) + N[2]*(density*bdf2[1]);

	array_1d<double,3> fbd_stblterm = ZeroVector(nodes_number);
	fbd_stblterm = tauone * prod(DN_DX,bdf);
	

	for(int ii = 0; ii< nodes_number; ++ii)
	  {
		int index = ii*(dof + 1) + dof;
		F[index] += area*fbd_stblterm[ii];
	  }
//COMMENT*********************************************************************************



	KRATOS_CATCH("")

	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateGradMassStblTerms(MatrixType& M,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const double tauone,const double area)
	{
		KRATOS_TRY
	int nodes_number = 3;
	int dof = 2;
	
	double density;
	double mu;
	double eps;
	CalculateDensity(GetGeometry(), density, mu, eps);

	//build 1*tau1*ro Nacc grad q)
	double fac = tauone*density;
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;

			//K(row,column) += -1*area * fac* N(ii) * DN_DX(jj,0);
			M(column,row) += area * fac* N[ii] * DN_DX(jj,0);

			//K(row + 1,column) += -1*area * fac* N(ii) * DN_DX(jj,1);
			M(column,row + 1) += area * fac* N[ii] * DN_DX(jj,1);
		   }
	    }

	KRATOS_CATCH("")

	}

	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateDarcyStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double time,const double tauone,const double area)
	 {
	 KRATOS_TRY

		int nodes_number = 3;
		int dof = 2;
		int matsize = nodes_number * dof;

		double density;
		double mu;
		double eps;
		CalculateDensity(GetGeometry(), density, mu, eps);

		//build grad p * Av+B|v|v and assemble
		//see function  CalculateGradStblAllTerms (line 1156)

		//build (rho a.grad u)* (Av+B|v|v)/rho and assemble
		//see function  CalculateAdvStblAllTerms (line 975)


		//build tau1(Au+B|u|u,(Av+B|v|v)/rho) stabilization term & assemble
		double dp = 0.01; //diameter of the particle	
		double kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*dp*dp);

		  const array_1d<double,3>& vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,0);
		  const array_1d<double,3>& vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,0);
		  const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,0);

		  array_1d<double,3> norm_vel_2 = ZeroVector(nodes_number);
		  //	vector with the norm^2 of the velocity of the three nodes at the previous iteration;
		  for (int ii = 0; ii < nodes_number; ii++)
		  {
			   norm_vel_2[0] += vel0[ii]*vel0[ii];
			   norm_vel_2[1] += vel1[ii]*vel1[ii];
			   norm_vel_2[2] += vel2[ii]*vel2[ii];
		  }

		  boost::numeric::ublas::bounded_matrix<double,2,6> darcy_opr = ZeroMatrix(dof,matsize);
		  for (int ii = 0; ii< nodes_number; ii++)
		  {
		        int column = ii*dof;
					       
		        //DARCY TERM linear part / rho
		        darcy_opr(0,column) = N[ii] *  kinv * mu / density;
		        darcy_opr(1,column + 1) = darcy_opr(0,column);
		        //DARCY TERM nonlinear part / rho
		        darcy_opr(0,column) += N[ii] * 1.75  /eps * sqrt(norm_vel_2[ii] *  kinv / (eps * 150.0));
		        darcy_opr(1,column + 1) += N[ii] * 1.75  /eps * sqrt(norm_vel_2[ii] *  kinv / (eps * 150.0));

		  }
		  boost::numeric::ublas::bounded_matrix<double,6,6> darcy_stblterm = ZeroMatrix(matsize,matsize);		
		  darcy_stblterm = tauone * prod(trans(darcy_opr),darcy_opr);

		  for ( int ii = 0; ii < nodes_number; ii++)
		    {
			 int row = ii*(dof+1);
			 int loc_row = ii*dof;
			 for( int jj=0; jj < nodes_number; jj++)
			   {
				int column = jj*(dof+1);
				int loc_column = jj*dof;
				
				K(row,column) += area*density *darcy_stblterm(loc_row,loc_column);
				K(row,column + 1) += area*density *darcy_stblterm(loc_row,loc_column + 1);
				K(row + 1,column) += area*density *darcy_stblterm(loc_row + 1,loc_column);
				K(row + 1,column + 1) += area*density *darcy_stblterm(loc_row + 1,loc_column + 1);

			   }
		    }

		//build tau1(Fbody ,(Av+B|v|v)/rho) stabilization term & assemble
		array_1d<double,2> bdf = ZeroVector(2);
		const array_1d<double,2> bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
		const array_1d<double,2> bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
		const array_1d<double,2> bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);


		bdf[0] = N[0]*(density*bdf0[0] ) +  N[1]*(density*bdf1[0]) + N[2]*(density*bdf2[0]);
		bdf[1] =  N[0]*(density*bdf0[1]) +  N[1]*(density*bdf1[1] ) + N[2]*(density*bdf2[1] );
		

		array_1d<double,6> fbd_stblterm = ZeroVector(matsize);
		fbd_stblterm = tauone * prod(trans(darcy_opr),bdf);

		for(int ii = 0; ii< nodes_number; ++ii)
		  {
			 int index = ii*(dof + 1);
			 int loc_index = ii*dof;
			 F[index] += area*fbd_stblterm[loc_index];
			 F[index + 1] += area*fbd_stblterm[loc_index + 1];
		  }


	 KRATOS_CATCH("")
	 }

	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateDarcyMassStblTerms(MatrixType& M,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const double tauone,const double area)
	 {
	 KRATOS_TRY

		int nodes_number = 3;
		int dof = 2;
		int matsize = dof*nodes_number;

		double density;
		double mu;
		double eps;
		CalculateDensity(GetGeometry(), density, mu, eps);


		//build tau1(rho Nacc,(Av+B|v|v)/rho) stabilization term & assemble
		double dp = 0.01; //diameter of the particle	
		double kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*dp*dp);

		  const array_1d<double,3>& vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,0);
		  const array_1d<double,3>& vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,0);
		  const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,0);

		  array_1d<double,3> norm_vel_2 = ZeroVector(nodes_number);
		  //	vector with the norm^2 of the velocity of the three nodes at the previous iteration;
		  for (int ii = 0; ii < nodes_number; ii++)
		  {
			   norm_vel_2[0] += vel0[ii]*vel0[ii];
			   norm_vel_2[1] += vel1[ii]*vel1[ii];
			   norm_vel_2[2] += vel2[ii]*vel2[ii];
		  }

		  boost::numeric::ublas::bounded_matrix<double,2,6> darcy_opr = ZeroMatrix(dof,matsize);
		  boost::numeric::ublas::bounded_matrix<double,6,2> shape_func = ZeroMatrix(matsize, dof);

		  for (int ii = 0; ii< nodes_number; ii++)
		  {
		        int column = ii*dof;
					       
		        //DARCY TERM linear part /rho
		        darcy_opr(0,column) = N[ii] *  kinv * mu / density;
		        darcy_opr(1,column + 1) = darcy_opr(0,column);
// KRATOS_WATCH("LIN PART")
// KRATOS_WATCH(darcy_opr)

		        //DARCY TERM nonlinear part /rho
		        darcy_opr(0,column) += N[ii] * 1.75  /eps * sqrt(norm_vel_2[ii] *  kinv / (eps * 150.0));
		        darcy_opr(1,column + 1) += N[ii] * 1.75  /eps * sqrt(norm_vel_2[ii] *  kinv / (eps * 150.0));
// KRATOS_WATCH("NON LIN PART")
// KRATOS_WATCH(darcy_opr)
		        shape_func(column,0) = N[ii];
		        shape_func(column + 1, 1) = shape_func(column,0);

		  }

	 	//tau1*rho*Nacc.(1.0*(Av+B|v|v))
		boost::numeric::ublas::bounded_matrix<double,6,6> temp_darcyterm = ZeroMatrix(matsize,matsize);
		temp_darcyterm = prod(trans(darcy_opr),trans(shape_func));
/*KRATOS_WATCH("Massmatr before")
KRATOS_WATCH(M)	*/	
		double fac = tauone*density;	
		for ( int ii = 0; ii < nodes_number; ii++)
		    {
			 int row = ii*(dof+1);
			 int loc_row = ii*dof;
			 for( int jj=0; jj < nodes_number; jj++)
			   {
				int column = jj*(dof+1);
				int loc_column = jj*dof;

				M(row,column) += area*fac*temp_darcyterm(loc_row,loc_column);
				M(row + 1,column + 1) += area*fac*temp_darcyterm(loc_row + 1,loc_column + 1);
			   }
		    }
// KRATOS_WATCH("Massmatr after darcy")
// KRATOS_WATCH(M)

	 KRATOS_CATCH("")
	 }

	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::AddBodyForceAndMomentum(VectorType& F,const array_1d<double,3>& N, const double time,const double area,const double tauone,const double tautwo)
	{
		KRATOS_TRY
	int nodes_number = 3;
	int dof = 2;
	

	//double lump_mass_fac = area * 0.333333333333333333333333;

	double density;
	double mu;
	double eps;
	CalculateDensity(GetGeometry(), density, mu, eps);

	//body  & momentum term force
	for ( int ii = 0; ii < nodes_number; ii++)
	   {
		int index = ii*(dof + 1) ;
		const array_1d<double,2> bdf = GetGeometry()[ii].FastGetSolutionStepValue(BODY_FORCE);

	
		F[index] += area*N[ii]*density*bdf[0] ;
		F[index + 1] += area*N[ii]*density*bdf[1];

	   }
	

	KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateResidual(const MatrixType& K, VectorType& F)
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
	void Fluid2DSplit::ComputeProjections(array_1d<double,6>& adv_proj , array_1d<double,3>& div_proj, const 			boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const double tauone,const double tautwo,const array_1d<double,3>& N,const double area, const double time)
	{
	unsigned int number_of_nodes = GetGeometry().PointsNumber();
	unsigned int dim = 2;

	double density;
	double mu;
	double eps;
	CalculateDensity(GetGeometry(), density, mu, eps);

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

	void Fluid2DSplit::AddProjectionForces(VectorType& F, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double area,const double tauone,const double tautwo)
	{
	unsigned int number_of_nodes = GetGeometry().PointsNumber();
	unsigned int dim = 2;

	double density;
	double mu;
	double eps;
	CalculateDensity(GetGeometry(), density, mu, eps);

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
	
		//tauone*ro*(xi,a.gradv)
		double proj;

		proj = mean_x_adv*(adv_vel[0]*DN_DX(ii,0) + adv_vel[1]*DN_DX(ii,1));
		F[index] += (tauone*1.0*area*proj);

		proj = mean_y_adv*(adv_vel[0]*DN_DX(ii,0) + adv_vel[1]*DN_DX(ii,1));
		F[index +1 ] += (tauone*1.0*area*proj);
	

		//tauone*(xi,gradq)
		proj = (mean_x_adv*DN_DX(ii,0) + mean_y_adv*DN_DX(ii,1));
		F[index + 2] += (tauone*area*proj);

		//tautwo*(divv)

		F[index] += (tautwo*area*mean_div*DN_DX(ii,0));
		F[index +1] += (tautwo*area*mean_div*DN_DX(ii,1));
		
	}

	


	}

// 	//************************************************************************************
// 	//************************************************************************************
// 	void Fluid2DSplit::GetElementalDensity(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
// 	{
// 
// 	//************************************************************************************
// 	//************************************************************************************
// 	void Fluid2DSplit::GetElementalPorosity(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
// 	{


	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;
		unsigned int node_size = dim+1;

		
			if(rResult.size() != number_of_nodes*node_size)
				rResult.resize(number_of_nodes*node_size,false);	
			
			if( this->GetValue(IS_DIVIDED) == 0.0 ) 
			{
				rResult.resize(0);
			}
			else
			{ 
				for (unsigned int i=0;i<number_of_nodes;i++)
				{
					 rResult[i*node_size] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
					 rResult[i*node_size+1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
					 rResult[i*node_size+2] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
				}
			}
		KRATOS_CATCH("")
		
	}

	//************************************************************************************
	//************************************************************************************
	  void Fluid2DSplit::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;
		unsigned int node_size = dim+1;


			if(ElementalDofList.size() != number_of_nodes*node_size)
				ElementalDofList.resize(number_of_nodes*node_size);	

			if( this->GetValue(IS_DIVIDED) == 0.0 ) 
			{
				ElementalDofList.resize(0);
			}
			else
			{ 
				for (unsigned int i=0;i<number_of_nodes;i++)
				{
					 ElementalDofList[i*node_size] = GetGeometry()[i].pGetDof(VELOCITY_X);
					 ElementalDofList[i*node_size+1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
					 ElementalDofList[i*node_size+2] = GetGeometry()[i].pGetDof(PRESSURE);
				}
			}
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************

   void Fluid2DSplit::Calculate( const Variable<array_1d<double,3> >& rVariable, 
                                    array_1d<double,3>& Output, 
                                    const ProcessInfo& rCurrentProcessInfo)
	{

	array_1d<double,6> adv_proj = ZeroVector(6);
	array_1d<double,3> div_proj = ZeroVector(3);

	double delta_t= rCurrentProcessInfo[DELTA_TIME];

	//getting data for the given geometry
	double Area;
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);

// 	double tauone = 0.0;//delta_t;
// 	double tautwo = 0.0;//delta_t;*/
        double tauone, tautwo;
 	CalculateTau(tauone, tautwo, delta_t, Area, rCurrentProcessInfo);

	ComputeProjections(adv_proj, div_proj, DN_DX,tauone,tautwo,N,Area, delta_t);


	}

	//************************************************************************************
	//************************************************************************************

	void Fluid2DSplit::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
	{

	double delta_t= rCurrentProcessInfo[DELTA_TIME];
		
	//getting data for the given geometry
	double Area;
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);
// /*	double tauone =0.0;// delta_t;
// 	double tautwo = 0.0;//delta_t;*/

       double tauone, tautwo;
	CalculateTau(tauone, tautwo, delta_t, Area, rCurrentProcessInfo);
            if(rVariable==THAWONE)
            {
                for( unsigned int PointNumber = 0; 
                     PointNumber < 1; 
                     PointNumber++ )
                {
                    rValues[PointNumber] = tauone;
                }
            }
            if(rVariable==THAWTWO)
            {
                for( unsigned int PointNumber = 0; 
                     PointNumber < 1; PointNumber++ )
                {
		
                    rValues[PointNumber] = tautwo;
                }
            }

        }
	//*************************************************************************************
	//*************************************************************************************
	void Fluid2DSplit::CalculateDensity(Geometry< Node<3> > geom, double& elemental_density, double& elemental_viscosity, double& elemental_porosity)
 	{//PAY ATTENTION: CALCULATION OF ELEMENTAL DENSITY WITHOUT THE EFFECT OF POROSITY : RHO === RHO_WATER
// 
// 	/*double kk = 0.0;
// 	for(int ii=0;ii<3;++ii)
// 		if(geom[ii].GetSolutionStepValue(IS_STRUCTURE) != 1.0)
// 			{
// 				kk++;
// 				density +=geom[ii].FastGetSolutionStepValue(DENSITY);
// 			}
// 
// 	density/=kk;*/
// 	/*
// 		density = ZeroVector(3);
// 	for(int ii=0;ii<3;++ii)
// 		density[ii] = geom[ii].FastGetSolutionStepValue(DENSITY);*/
// 	
// 	/*const double rho0 = geom[0].FastGetSolutionStepValue(DENSITY);
// 	const double rho1 = geom[1].FastGetSolutionStepValue(DENSITY);
// 	const double rho2 = geom[2].FastGetSolutionStepValue(DENSITY);
// 	 density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );*/
// 
// 

	////Check if some of the elements don't have a porosity assigned
	double eps0 = geom[0].FastGetSolutionStepValue(POROSITY);
	if(eps0 == 0.0)
	 {
	     eps0 = 1.0;
	 }
	double eps1 = geom[1].FastGetSolutionStepValue(POROSITY);
	if(eps1 == 0.0)
	 {
	     eps1 = 1.0;
	 }
	double eps2 = geom[2].FastGetSolutionStepValue(POROSITY);
	if(eps2 == 0.0)
	 {
	     eps2 = 1.0;
	 }

	elemental_density = 0.0;
	elemental_porosity = 1.0;

	if(eps0 == eps1 && eps1 == eps2)
	  {
		//for inside the domain totally inside one fluid
		elemental_porosity = eps0;
		elemental_density = geom[0].FastGetSolutionStepValue(DENSITY); //* eps0;	
		elemental_viscosity = geom[0].FastGetSolutionStepValue(VISCOSITY) * elemental_density;	//mu = nu * density //we assigne nu=1E-6 from Gid 
// 		KRATOS_WATCH("fluid nodes")
// 		KRATOS_WATCH(geom[0].Id();)
// 		KRATOS_WATCH(geom[1].Id();)
// 		KRATOS_WATCH(geom[2].Id();)
	  }
	else if(eps0 == eps1)
	 {
	   	elemental_porosity = eps0;
		elemental_density = geom[0].FastGetSolutionStepValue(DENSITY); // * eps0;	
		elemental_viscosity = geom[0].FastGetSolutionStepValue(VISCOSITY) * elemental_density;	//mu = nu * density 
	 }
	else if(eps1 == eps2)
	 {
		elemental_porosity = eps1;
		elemental_density = geom[1].FastGetSolutionStepValue(DENSITY); // * eps1;	
		elemental_viscosity = geom[1].FastGetSolutionStepValue(VISCOSITY) * elemental_density;  //mu = nu * density 
	 }
	else if(eps2 == eps0)
	 {
		elemental_porosity = eps2;
		elemental_density = geom[2].FastGetSolutionStepValue(DENSITY); // * eps2;	
		elemental_viscosity = geom[2].FastGetSolutionStepValue(VISCOSITY)* elemental_density;  //mu = nu * density 
	 }
	else { KRATOS_WATCH("ERROR!!! three different values of densities");}
/*
	KRATOS_WATCH("nodes of the element")
	KRATOS_WATCH(geom[0].Id())
        KRATOS_WATCH(geom[1].Id())
	KRATOS_WATCH(geom[2].Id())
	KRATOS_WATCH(elemental_porosity)
	KRATOS_WATCH(elemental_density)
	KRATOS_WATCH(elemental_viscosity)
      */

// KRATOS_WATCH(elemental_porosity);
// KRATOS_WATCH(elemental_density);
// KRATOS_WATCH(elemental_viscosity);

// 	  {
// 	 
// 		//for element having common node between two fluids or boundary element with POROSITY==1 inside the domain
// 		for(int ii = 0 ; ii < geom.size();++ii)
// 			{
// 			  if(geom[ii].GetSolutionStepValue(IS_POROUS) == 1.0 && geom[ii].GetSolutionStepValue(IS_STRUCTURE) != 1.0)
// 				{
// 			  	 elemental_density = geom[ii].FastGetSolutionStepValue(DENSITY);
// 				 elemental_viscosity = geom[ii].FastGetSolutionStepValue(VISCOSITY);
// 
// 				}
// 			}
// 		//for boundary element with IS_POROUS==1 on the boundary
// 		if(elemental_density == 0.0)
// 			 for(int ii = 0 ; ii < geom.size();++ii)
// 			 {
// 			  if(geom[ii].GetSolutionStepValue(IS_POROUS) == 0.0)
// 				{
// 				 elemental_density = geom[ii].FastGetSolutionStepValue(DENSITY);
// 				 elemental_viscosity = geom[ii].FastGetSolutionStepValue(VISCOSITY);
// 
// 			
// 				}
// 			 }	
// 			
 	}
	//*************************************************************************************
	//*************************************************************************************

		void Fluid2DSplit::CalculateTau(double& tauone, double& tautwo, const double time,const double area,const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
	//calculate mean advective velocity and taus
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
	//const double mu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
	//const double mu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
	//const double mu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
	//mu = 0.333333333333333333333333*(mu0 + mu1 + mu2);

	double density;
	double eps;
	CalculateDensity(GetGeometry(), density, mu, eps);

	int dyn_st_switch = rCurrentProcessInfo[DYNAMIC_TAU];
	
		if(dyn_st_switch)
		  {
		
			tauone = 1.0/(1.0/time + 4.0*mu/(ele_length*ele_length*density)+2.0*advvel_norm*1.0/ele_length);
		  }
		else
		 {
			
			tauone = 1.0/(0.0+ 4.0*mu/(ele_length*ele_length*density)+2.0*advvel_norm*1.0/ele_length);
		  }
		
	tautwo = mu/density + 1.0*ele_length*advvel_norm/2.0;


	KRATOS_CATCH("")

	}


	//*************************************************************************************
	//*************************************************************************************

	  void Fluid2DSplit::GetFirstDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes * (dim + 1);
		if(values.size() != MatSize)   values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i * (dim + 1);
			values[index] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y,Step);
			values[index + 2] = GetGeometry()[i].GetSolutionStepValue(PRESSURE,Step);

		}
	}
	//************************************************************************************
	//************************************************************************************
	  void Fluid2DSplit::GetSecondDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes * (dim + 1);
		if(values.size() != MatSize) values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i * (dim + 1);
			values[index] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y,Step);
			values[index + 2] = 0.0;
		}
	}
	//************************************************************************************
	//************************************************************************************
} // Namespace Kratos


