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

        array_1d<double,3> conv_opr = ZeroVector(3); //
        #pragma omp threadprivate(conv_opr)

        array_1d<double,3> darcy_opr = ZeroVector(3); //dimension = number of nodes
        #pragma omp threadprivate(darcy_opr)

        array_1d<double,6> div_opr = ZeroVector(6);
        #pragma omp threadprivate(div_opr)

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
	//***************************************boost::numeric::ublas::bounded_matrix<double,4,2> *********************************************
	void Fluid2DSplit::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

	int nodes_number = 3;
	int dim = 2;
	unsigned int matsize = nodes_number*(dim+1);

	if(rLeftHandSideMatrix.size1() != matsize)
			rLeftHandSideMatrix.resize(matsize,matsize, false); //false says not to preserve existing storage!!

	if(rRightHandSideVector.size() != matsize)
			rRightHandSideVector.resize(matsize, false); //false says not to preserve existing storage!!

	
	noalias(rLeftHandSideMatrix) = ZeroMatrix(matsize,matsize); 
	noalias(rRightHandSideVector) = ZeroVector(matsize); 
	
	double delta_t= rCurrentProcessInfo[DELTA_TIME];
	
	
	//getting data for the given geometry
	double Area;
	//The shape functions are calculated on the unique gauss point (their value is then 1/3). 
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
		  if (dist_on_agp[i] < 0.0)
		  {
		        double Area_se = A_on_agp(i);
		      
		        for (unsigned int j = 0; j < N.size(); j++)
			     N[j] = N_on_agp(i,j);
// 		        double tauone = 0.0;//delta_t;
// 		        double tautwo = 0.0;//delta_t;*/
		        double tauone, tautwo;
 		        CalculateTau(tauone, tautwo, delta_t, Area ,  rCurrentProcessInfo);


		        //add body force and momentum
		        AddBodyForceAndMomentum(rRightHandSideVector, N, delta_t, Area_se, tauone,tautwo);
		        
		        //add volume correction - q * div vn
// 		        AddVolumeCorrection(rRightHandSideVector, N, delta_t, Area_se) ;

		  }

		}

	 }
	 else if(GetValue(IS_DIVIDED) == -1.0)
	   {

// 		double tauone = 0.0;//delta_t;
// 		double tautwo = 0.0;//delta_t;
		double tauone, tautwo;
		CalculateTau(tauone, tautwo, delta_t, Area,  rCurrentProcessInfo);


		//add body force and momentum
		AddBodyForceAndMomentum(rRightHandSideVector, N, delta_t, Area, tauone,tautwo);
		       
		//add volume correction - q * div vn
// 		AddVolumeCorrection(rRightHandSideVector, N, delta_t, Area) ;

	   }
	 

		KRATOS_CATCH("")
	}
	//***********************************************************************************
	//**************************************************************************************
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
	double density;
	double mu;
	double eps;
	double dp;
	CalculateDensity(GetGeometry(), density, mu, eps, dp);
	double lump_mass_fac = density * area * 0.333333333333333333333333;

	int nodes_number = 3;
	int dof = 2;
	for ( int nd = 0; nd< nodes_number; nd++)
	    {
		int row = nd*(dof + 1);
		for( int jj=0; jj< dof; jj++)
			K(row + jj, row + jj) +=  lump_mass_fac;
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
// 		        int matsize = dof*nodes_number;
		        double density;
		        double mu;
		        double eps;
		        double dp;
		        CalculateDensity(GetGeometry(), density, mu, eps, dp);

		        boost::numeric::ublas::bounded_matrix<double,4,2> aux_gp = ZeroMatrix(4,2);
		        array_1d<double,4> A_on_agp = ZeroVector(4);
		        boost::numeric::ublas::bounded_matrix<double,4,3> N_on_agp = ZeroMatrix(4,3);
		        array_1d<double,4> dist_on_agp = ZeroVector(4);

		        DivideElemUtils::DivideElement_2D(GetGeometry(), aux_gp, A_on_agp, N_on_agp, dist_on_agp);
		        
		        boost::numeric::ublas::bounded_matrix<double,3,3> temp_MassMatr = ZeroMatrix(dof + 1, dof +1);

		        for(unsigned int i = 0 ; i< aux_gp.size1() ; i++)
		        {
			 //if the virtual sub element is a fluid element
			 if (dist_on_agp[i] < 0.0)
			 {
			       double Area = A_on_agp(i);
			       double fac = Area * density;
			       //shape functions of the element calculated on the auxiliary gauss points (gp of the sub elements)
			       for (unsigned int j = 0; j < N.size(); j++)
				    N[j] = N_on_agp(i,j);

		      // |  N0(agp_i)* N0(agp_i)		N0(agp_i)* N1(agp_i)	N0(agp_i)* N2(agp_i)	 |	
        //temp_MassMatrix=	|N1(agp_i)* N0(agp_i)		N1(agp_i)* N1(agp_i)	N1(agp_i)* N2(agp_i)	  | * density * Area	(sum on agp_i)
		      // |  N2(agp_i)* N0(agp_i)		N2(agp_i)* N1(agp_i)	N2(agp_i)* N2(agp_i)	  |	
			       
			       noalias(temp_MassMatr) +=  fac * outer_prod(N,N);	    

			 }
		        }

		        for ( int ii = 0; ii < nodes_number; ii++)
		        {
			     int row = ii*(dof+1);
			     for( int jj=0; jj < nodes_number; jj++)
			       {
				    rMassMatrix(row,row) += temp_MassMatr(ii,jj);
				    rMassMatrix(row + 1,row + 1) += temp_MassMatr(ii,jj);

			       }
		        }


		        for(unsigned int i = 0 ; i< aux_gp.size1() ; i++)
		        {
			   //if the virtual sub element is a fluid element
			   if (dist_on_agp[i] < 0.0)
			   {
				double Area_se = A_on_agp[i];
				//shape functions of the element calculated on the auxiliary gauss points (gp of the sub elements)
				for (unsigned int j = 0; j < N.size(); j++)
				      N[j] = N_on_agp(i,j);

	 // 		        double tauone = 0.0;//delta_t;
	 // 		        double tautwo = 0.0;//delta_t;*/
				double tauone, tautwo;
				CalculateTau(tauone, tautwo, delta_t, Area ,  rCurrentProcessInfo);
				//add stablilization terms due to advective term (a)grad(V) * rho*Acc
				CalculateAdvMassStblTerms(rMassMatrix, DN_DX, N,tauone,Area_se);
				//add stablilization terms due to grad term grad(q) * rho*Acc
				CalculateGradMassStblTerms(rMassMatrix, DN_DX, tauone,Area_se);
			   }	
		        }
		}
		else if(GetValue(IS_DIVIDED) == -1.0)
		  {
		        //Calculate tau
        /*		double tauone = 0.0;//delta_t;
		        double tautwo = 0.0;// delta_t;*/
		        double tauone, tautwo;
		        CalculateTau(tauone, tautwo, delta_t, Area,  rCurrentProcessInfo);
		        CalculateMassContribution(rMassMatrix,delta_t,Area); 
		        //add stablilization terms due to advective term (a)grad(V) * ro*Acce
		        CalculateAdvMassStblTerms(rMassMatrix, DN_DX, N,tauone,Area);
		        //add stablilization terms due to grad term grad(q) * ro*Acce
		        CalculateGradMassStblTerms(rMassMatrix, DN_DX, tauone,Area);

		  }
	 
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
			 if (dist_on_agp[i] < 0.0)
			 {
			       double Area_se = A_on_agp[i];
			     
			       for (unsigned int j = 0; j < N.size(); j++)
				    N[j] = N_on_agp(i,j);

			       //viscous term	
			       CalculateViscousTerm(rDampMatrix, DN_DX, Area);
			       //Advective term
        // 		        double tauone = 0.0;//delta_t;
        // 		        double tautwo = 0.0;//delta_t;
			       double tauone, tautwo;
			       CalculateTau(tauone, tautwo, delta_t, Area,  rCurrentProcessInfo);
			       CalculateAdvectiveTerm(rDampMatrix, DN_DX, tauone, tautwo, delta_t, Area_se);
        // 		        //calculate pressure term
			       CalculatePressureTerm(rDampMatrix, DN_DX, N, delta_t,Area_se);
			       //calculate Darcy term
			       CalculateDarcyTerm_SubElem(rDampMatrix, N, Area_se);

			       //compute projections

			       //stabilization terms
		//	       CalculateDivStblTerm(rDampMatrix, DN_DX, tautwo, Area_se);
			       CalculateAdvStblAllTerms(rDampMatrix,rRightHandSideVector, DN_DX, N, tauone,delta_t, Area_se);
			       CalculateGradStblAllTerms(rDampMatrix,rRightHandSideVector,DN_DX, delta_t, tauone, Area_se);

			 }

		        }
		}
		else if(GetValue(IS_DIVIDED) == -1.0)
		{

			 CalculateViscousTerm(rDampMatrix, DN_DX, Area);
		
			 //Advective term
			 //		  double tauone =0.0; //delta_t;
			 //		  double tautwo =0.0;// delta_t;
			 double tauone, tautwo;
			 CalculateTau(tauone, tautwo, delta_t, Area,  rCurrentProcessInfo);
			 CalculateAdvectiveTerm(rDampMatrix, DN_DX, tauone, tautwo, delta_t, Area);		   

			 //calculate pressure term
			 CalculatePressureTerm(rDampMatrix, DN_DX, N, delta_t,Area);

			 //calculate Darcy term
			 CalculateDarcyTerm(rDampMatrix, Area);

			 //compute projections

			 //stabilization terms
		//	 CalculateDivStblTerm(rDampMatrix, DN_DX, tautwo, Area);
			 CalculateAdvStblAllTerms(rDampMatrix,rRightHandSideVector, DN_DX, N, tauone,delta_t, Area);
			 CalculateGradStblAllTerms(rDampMatrix,rRightHandSideVector,DN_DX, delta_t, tauone, Area);
	       
		}

		KRATOS_CATCH("")
	}
        


	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateViscousTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double area)
	{
		KRATOS_TRY
	double mu;
	double density;
	double eps;
	double dp;
	CalculateDensity(GetGeometry(), density, mu, eps, dp);

	//nu = mu/density;	

	int nodes_number = 3;
	int dof = 2;
	 double fac = mu * area;
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		for( int jj=0; jj < nodes_number; jj++)
		   {
			 int column = jj*(dof+1);
			 K(row,column) += fac *(DN_DX(ii,0)*DN_DX(jj,0) + DN_DX(ii,1)*DN_DX(jj,1));
			 K(row + 1,column + 1) += fac *(DN_DX(ii,0)*DN_DX(jj,0) + DN_DX(ii,1)*DN_DX(jj,1));
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
	const array_1d<double,3>& adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
	const array_1d<double,3>& adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
	const array_1d<double,3>& adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);

	double density;
	double mu;
	double eps;
	double dp;
	CalculateDensity(GetGeometry(), density, mu, eps, dp);

	ms_adv_vel[0] = N[0]*adv_vel0[0]+N[1]*adv_vel1[0]+N[2]*adv_vel2[0];
	ms_adv_vel[1] = N[0]*adv_vel0[1]+N[1]*adv_vel1[1]+N[2]*adv_vel2[1];
//*****Convective vel == real fluid velocity
ms_adv_vel[0] /=eps;
ms_adv_vel[1] /=eps;
//****
	//ms_adv_vel[0] = 0.0;
	//ms_adv_vel[1] = 0.0;
	                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
	//calculate convective term	
	int nodes_number = 3;
	int dof = 2;

	for (int ii = 0; ii< nodes_number; ii++)
	    {
		conv_opr[ii] = DN_DX(ii,0)*ms_adv_vel[0] + DN_DX(ii,1)*ms_adv_vel[1];
	    }
	boost::numeric::ublas::bounded_matrix<double,3,3> temp_convterm;

	noalias(temp_convterm) = outer_prod(N, conv_opr);

	double fac = area * density;
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);

			K(row,column) += fac * temp_convterm(ii,jj);
			K(row + 1,column + 1) += fac * temp_convterm(ii,jj);
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
		double dp;
		CalculateDensity(GetGeometry(), density, mu, eps, dp);
			array_1d<double,3> nodal_eps = ZeroVector(3);
	for (int i = 0; i < 3; i++)
	 {
	  nodal_eps[i] = GetGeometry()[i].FastGetSolutionStepValue(POROSITY);
	 }
// KRATOS_WATCH("line 535")
// KRATOS_WATCH(eps)
// KRATOS_WATCH(dp)
// 		double dp = 0.01; //diameter of the particle	
		double kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*dp*dp);
		const array_1d<double,3>& vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);

		int nodes_number = 3;
		int dof = 2;

		array_1d<double,3> norm_vel_2 = ZeroVector(nodes_number);
		array_1d<double,3> vel_gp = ZeroVector(nodes_number);
		//	vector with the norm^2 of the velocity of the three nodes at the previous iteration;
		for (int ii = 0; ii < nodes_number; ii++)
		{
			 norm_vel_2[0] += vel0[ii]*vel0[ii];
			 norm_vel_2[1] += vel1[ii]*vel1[ii];
			 norm_vel_2[2] += vel2[ii]*vel2[ii];
			 vel_gp[ii] = N[0]*vel0[ii] + N[1]*vel1[ii] + N[2]*vel2[ii];
		}

		double norm_vel_gp = 0.0;
		for (int ii = 0; ii < nodes_number; ii++)
		{
			 norm_vel_gp +=  vel_gp[ii]*vel_gp[ii];
		}
/*
		      // |  N0(agp_i)* N0(agp_i)	N0(agp_i)* N1(agp_i)	N0(agp_i)* N2(agp_i)  |	
        //temp_sfprod=	| N1(agp_i)* N0(agp_i)	N1(agp_i)* N1(agp_i)	N1(agp_i)* N2(agp_i)  | * Darcy term	
		      // |  N2(agp_i)* N0(agp_i)	N2(agp_i)* N1(agp_i)	N2(agp_i)* N2(agp_i)  |	
			       
*/
	    
		boost::numeric::ublas::bounded_matrix<double,3,3> temp_sfprod = ZeroMatrix(dof +1,dof +1);
		noalias(temp_sfprod) = area * outer_prod(N,N);

		// Lumped form
		double fac_linear = kinv * mu;	
		double fac_nonlinear = (1.75 * density /eps * sqrt(norm_vel_gp *  kinv / (eps * 150.0))) ;

		for ( int ii = 0; ii < nodes_number; ii++)
		    {
			 int row = ii*(dof+1);
			 for( int jj=0; jj < nodes_number; jj++)
			   {
				//DARCY TERM linear part
				K(row,row) += nodal_eps[ii] * fac_linear * temp_sfprod(ii,jj);
				K(row + 1,row + 1) += nodal_eps[ii] * fac_linear * temp_sfprod(ii,jj);

				//DARCY TERM nonlinear part
				K(row    ,row    ) += nodal_eps[ii] * fac_nonlinear * temp_sfprod(ii,jj);
				K(row + 1,row + 1) += nodal_eps[ii] * fac_nonlinear * temp_sfprod(ii,jj);
			   }
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
		double dp;
		CalculateDensity(GetGeometry(), density, mu, eps, dp);
		array_1d<double,3> nodal_eps = ZeroVector(3);
	for (int i = 0; i < 3; i++)
	 {
	  nodal_eps[i] = GetGeometry()[i].FastGetSolutionStepValue(POROSITY);
	 }
/*KRATOS_WATCH("line 604")
KRATOS_WATCH(eps)
KRATOS_WATCH(dp)	*/	
// 		double dp = 0.01; //diameter of the particle	
		double kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*dp*dp);

		const array_1d<double,3>& vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
		
		int nodes_number = 3;
		int dof = 2;

		//	vector with the norm^2 of the velocity of the three nodes at the previous iteration;
		array_1d<double,3> norm_vel_2 = ZeroVector(nodes_number);
		array_1d<double,3> vel_gp = ZeroVector(nodes_number);
		//	vector with the norm^2 of the velocity of the three nodes at the previous iteration;
		for (int ii = 0; ii < nodes_number; ii++)
		{
			 norm_vel_2[0] += vel0[ii]*vel0[ii];
			 norm_vel_2[1] += vel1[ii]*vel1[ii];
			 norm_vel_2[2] += vel2[ii]*vel2[ii];
			 vel_gp[ii] = 0.33333333333 * (vel0[ii] + vel1[ii] + vel2[ii]);
		}
		double norm_vel_gp = 0.0;
		for (int ii = 0; ii < nodes_number; ii++)
		{
			 norm_vel_gp +=  vel_gp[ii]*vel_gp[ii];
		}		

		double fac_linear = kinv * mu;	
		double fac_nonlinear = (1.75 * density /eps * sqrt(norm_vel_gp *  kinv / (eps * 150.0))) ;

		for ( int nd = 0; nd< nodes_number; nd++)
		{
		      int row = nd*(dof + 1);
		      for( int jj=0; jj< dof; jj++)
		      {	   //DARCY TERM linear part:
			   K(row + jj, row + jj) +=  nodal_eps[nd] * fac_linear * lump_mass_fac;
			   //DARCY TERM nonlinear part:
			   K(row + jj, row + jj) += nodal_eps[nd] * fac_nonlinear * lump_mass_fac;
		      }
		}

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
	double dp;
	CalculateDensity(GetGeometry(), density, mu, eps, dp);
	
	array_1d<double,3> nodal_eps = ZeroVector(3);
	for (int i = 0; i < nodes_number; i++)
	 {
	  nodal_eps[i] = GetGeometry()[i].FastGetSolutionStepValue(POROSITY);
	 }
// KRATOS_WATCH(nodal_eps)
        //Pj mean point of the edge ik
        //N_on_meanp(i,j) = shape function od the node i calculated on the point Pj
	boost::numeric::ublas::bounded_matrix<double,3,3> N_on_meanp = ZeroMatrix(3,3);
	N_on_meanp(0,0) = 0.0; N_on_meanp(0,1) = 0.5; N_on_meanp(0,2) = 0.5;
	N_on_meanp(1,0) = 0.5; N_on_meanp(1,1) = 0.0; N_on_meanp(1,2) = 0.5;
	N_on_meanp(2,0) = 0.5; N_on_meanp(2,1) = 0.5; N_on_meanp(2,2) = 0.0;


	// the matrix l_n contain the product of edge lengh l_jk*n_i (1st column x-comp, 2nd column y-comp, i = row index).
	boost::numeric::ublas::bounded_matrix<double,3,2> l_n = ZeroMatrix(3,2); 
	noalias(l_n) = -2.0 * area *DN_DX;


	 double fac = area *density;
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;
			//**************************************************************
 			//Elemental gradient of pressure term (momentum equation)
			K(row,column) -= area * nodal_eps[ii] * N[jj] * DN_DX(ii,0);
			K(row + 1,column) -= area * nodal_eps[ii] * N[jj] * DN_DX(ii,1);

			//**************************************************************
// 			//Elemental divergence terms (continuity equation)
// 	 		// Fomulation n1  int( q * rho * Div( u ))
			K(column,row) += fac * N[jj] * DN_DX(ii,0);
			K(column,row + 1) += fac * N[jj] * DN_DX(ii,1);

		   }
	    }


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

	for(int ii=0; ii<nodes_number; ii++)
	  {
		int index = dof*ii;
		div_opr[index] = DN_DX(ii,0);
		div_opr[index + 1] = DN_DX(ii,1);
	  }


	double density;
	double mu;
	double eps;
	double dp;
	CalculateDensity(GetGeometry(), density, mu, eps, dp);
	array_1d<double,3> nodal_eps = ZeroVector(3);
	for (int i = 0; i < nodes_number; i++)
	 {
		  nodal_eps[i] = GetGeometry()[i].FastGetSolutionStepValue(POROSITY);
	}

	boost::numeric::ublas::bounded_matrix<double,6,6> temp_div = ZeroMatrix(matsize,matsize);
	noalias(temp_div) = tautwo * outer_prod(div_opr,div_opr);

	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		int loc_row = ii*dof;
/**/
		double fac = area*density /nodal_eps[ii];
/**/
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);
			int loc_column = jj*dof;

			K(row,column) += fac*temp_div(loc_row,loc_column);
			K(row,column + 1) += fac*temp_div(loc_row,loc_column + 1);
			K(row + 1,column) += fac*temp_div(loc_row + 1,loc_column);
			K(row + 1,column + 1) += fac*temp_div(loc_row + 1,loc_column + 1);

		   }
	    }

		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateAdvStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N, const double tauone,const double time,const double area)
	{

	KRATOS_TRY
		const array_1d<double,3>& adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);

		ms_adv_vel[0] = N[0]*adv_vel0[0]+N[1]*adv_vel1[0]+N[2]*adv_vel2[0];
		ms_adv_vel[1] = N[0]*adv_vel0[1]+N[1]*adv_vel1[1]+N[2]*adv_vel2[1];

		//calculate convective term	
		int nodes_number = 3;
		int dof = 2;
		int matsize = dof*nodes_number;

		double density;
		double mu;
		double eps;
		double dp;
		CalculateDensity(GetGeometry(), density, mu, eps, dp);

//*****Convective vel == real fluid velocity
ms_adv_vel[0] /=eps;
ms_adv_vel[1] /=eps;
//****


		for (int ii = 0; ii< nodes_number; ii++)
		    {
			 conv_opr[ii] = DN_DX(ii,0)*ms_adv_vel[0] + DN_DX(ii,1)*ms_adv_vel[1];
			 
		    }

		//build (a.grad V)(ro*a.grad U) stabilization term & assemble
		boost::numeric::ublas::bounded_matrix<double,3,3> adv_stblterm = ZeroMatrix(dof+1,dof +1);		
		adv_stblterm = tauone * outer_prod(conv_opr,conv_opr);
		array_1d<double,3> nodal_eps = ZeroVector(3);
		for (int i = 0; i < nodes_number; i++)
		  {
		  nodal_eps[i] = GetGeometry()[i].FastGetSolutionStepValue(POROSITY);
		  }

// KRATOS_WATCH("line 779")
// KRATOS_WATCH(eps)
// KRATOS_WATCH(dp)
		double fac = area * density;
		for ( int ii = 0; ii < nodes_number; ii++)
		    {
			 int row = ii*(dof+1);
			 for( int jj=0; jj < nodes_number; jj++)
			   {
				int column = jj*(dof+1);
			
				K(row,column) += fac * adv_stblterm(ii,jj);
				K(row + 1,column + 1) +=  fac * adv_stblterm(ii,jj);
			   }
		    }
	
	       //build Darcy Au+B|u|u * (a.grad v) stabilization term & assemble + build (rho a.grad u)* (Av+B|v|v)/rho   stabilization term & assemble
// 	       double dp = 0.01; //diameter of the particle	
	       double kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*dp*dp);

		array_1d<double,3> norm_vel_2 = ZeroVector(nodes_number);
		array_1d<double,3> vel_gp = ZeroVector(nodes_number);
		//	vector with the norm^2 of the velocity of the three nodes at the previous iteration;
		for (int ii = 0; ii < nodes_number; ii++)
		{
			 norm_vel_2[0] += adv_vel0[ii]*adv_vel0[ii];
			 norm_vel_2[1] += adv_vel1[ii]*adv_vel1[ii];
			 norm_vel_2[2] += adv_vel2[ii]*adv_vel2[ii];
			 vel_gp[ii] = N[0]*adv_vel0[ii] + N[1]*adv_vel1[ii] + N[2]*adv_vel2[ii];
		}
		double norm_vel_gp = 0.0;
		for (int ii = 0; ii < nodes_number; ii++)
		{
			 norm_vel_gp +=  vel_gp[ii]*vel_gp[ii];
		}

		double fac_linear = kinv * mu / density;
		double fac_nonlinear = 1.75  /eps * sqrt(norm_vel_gp *  kinv / (eps * 150.0));
		for (int ii = 0; ii< nodes_number; ii++)
		    {
			 //DARCY TERM linear part /rho
			 darcy_opr[ii] = N[ii] * fac_linear * nodal_eps[ii];
			 //DARCY TERM nonlinear part /rho
			 darcy_opr[ii] += N[ii] * fac_nonlinear * nodal_eps[ii];
		    }

		boost::numeric::ublas::bounded_matrix<double,3,3> darcy_stblterm = ZeroMatrix(dof +1,dof +1);
		
		darcy_stblterm = tauone * outer_prod(conv_opr,darcy_opr);

		  for ( int ii = 0; ii < nodes_number; ii++)
		    {
			 int row = ii*(dof+1);
			 for( int jj=0; jj < nodes_number; jj++)
			   {
				int column = jj*(dof+1);
// 				Au+B|u|u * (a.grad v)

				K(row,column) += fac *darcy_stblterm(ii,jj);
				K(row + 1,column + 1) += fac *darcy_stblterm(ii,jj);

			   }
		    }

// 		array_1d<double,3> nodal_eps = ZeroVector(3);
// 		for (int i = 0; i < nodes_number; i++)
// 		  {
// 		  nodal_eps[i] = GetGeometry()[i].FastGetSolutionStepValue(POROSITY);
// 		  }
// 		//build 1*tau1*(a.grad V)(grad P) & 1*tau1*(grad q)(ro*a.grad U) stabilization terms & assemble
		boost::numeric::ublas::bounded_matrix<double,6,3> grad_stblterm = ZeroMatrix(matsize,nodes_number);
		for ( int ii = 0; ii < nodes_number; ii++)
		{
			int row = ii*dof;
			for ( int jj = 0; jj < nodes_number; jj++)
		        {
				grad_stblterm(row,jj)   =  DN_DX(jj,0) * conv_opr(ii);
  				grad_stblterm(row+1,jj) =  DN_DX(jj,1) * conv_opr(ii);
			 }
		}

 		grad_stblterm *= tauone * area;

		for ( int ii = 0; ii < nodes_number; ii++)
		    {
			 int row = ii*(dof+1);
			 int loc_row = ii*dof;
			 for( int jj=0; jj < nodes_number; jj++)
			   {
				int column = jj*(dof+1) + dof;
/**/
				//1*tau1*(a.grad V)(grad P)
				K(row,column) += nodal_eps[ii] * grad_stblterm(loc_row,jj);
				K(row + 1,column) += nodal_eps[ii] * grad_stblterm(loc_row + 1, jj);
// 
// 				//1*tau1*(grad q)(ro*a.grad U) 
// 				K(column, row) += density * grad_stblterm(loc_row, jj);
// 				K(column, row + 1) += density * grad_stblterm(loc_row + 1, jj);

// 				//1*tau1*(a.grad V)(grad P)
// 				K(row,column) +=  grad_stblterm(loc_row,jj);
// 				K(row + 1,column) +=  grad_stblterm(loc_row + 1, jj);

				//1*tau1*(grad q)(ro*a.grad U) 
				K(column, row) += density * grad_stblterm(loc_row, jj)  /nodal_eps[jj];
				K(column, row + 1) += density * grad_stblterm(loc_row + 1, jj)  /nodal_eps[jj];
/**/			  
			  }
		    }


		//build (1.0*a.grad V) (Fbody) stabilization term & assemble
		array_1d<double,2> bdf = ZeroVector(2);

		const array_1d<double,2> bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
		const array_1d<double,2> bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
		const array_1d<double,2> bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);


		bdf[0] =  (N[0]*bdf0[0] +  N[1]*bdf1[0] + N[2]*bdf2[0]);
		bdf[1] =  (N[0]*bdf0[1] +  N[1]*bdf1[1] + N[2]*bdf2[1]);
		
		fac = tauone * area * density; 
		for(int ii = 0; ii< nodes_number; ++ii)
		  {
			 int index = ii*(dof + 1);
			 F[index]     += nodal_eps[ii] * fac * conv_opr[ii] * bdf[0];
			 F[index + 1] += nodal_eps[ii] * fac * conv_opr[ii] * bdf[1];

		  }

	KRATOS_CATCH("")

	}
	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::CalculateAdvMassStblTerms(MatrixType& M,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N, const double tauone,const double area)
	{
		KRATOS_TRY
	const array_1d<double,3>& adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
	const array_1d<double,3>& adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
	const array_1d<double,3>& adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);

	ms_adv_vel[0] = N[0]*adv_vel0[0]+N[1]*adv_vel1[0]+N[2]*adv_vel2[0];
	ms_adv_vel[1] = N[0]*adv_vel0[1]+N[1]*adv_vel1[1]+N[2]*adv_vel2[1];


	//calculate convective term	
	int nodes_number = 3;
	int dof = 2;
//	int matsize = dof*nodes_number;

	//calculate density
	double density;
	double mu;
	double eps;
	double dp;
	CalculateDensity(GetGeometry(), density, mu, eps, dp);
	array_1d<double,3> nodal_eps = ZeroVector(3);
	for (int i = 0; i < nodes_number; i++)
	  {
	  nodal_eps[i] = GetGeometry()[i].FastGetSolutionStepValue(POROSITY);
	  }
//*****Convective vel == real fluid velocity
ms_adv_vel[0] /=eps;
ms_adv_vel[1] /=eps;
//****

	for (int ii = 0; ii< nodes_number; ii++)
	    {
		conv_opr[ii] = DN_DX(ii,0)*ms_adv_vel[0] + DN_DX(ii,1)*ms_adv_vel[1];
	    }
	boost::numeric::ublas::bounded_matrix<double,3,3> temp_convterm;
	

	//tau1*rho*Nacc.(1.0*a.grad V)
	noalias(temp_convterm) = outer_prod(conv_opr,N);

	double fac = tauone*density*area;	
	
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1);

			M(row,column) += fac*temp_convterm(ii,jj) /*/ nodal_eps[ii]*/;
			M(row + 1,column + 1) += fac*temp_convterm(ii,jj) /*/ nodal_eps[ii]*/;
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
	double dp;
	CalculateDensity(GetGeometry(), density, mu, eps, dp);
	array_1d<double,3> nodal_eps = ZeroVector(3);
	for (int i = 0; i < nodes_number; i++)
	  {
	  nodal_eps[i] = GetGeometry()[i].FastGetSolutionStepValue(POROSITY);
	  }
	//build 1*(grad q . grad p) stabilization term & assemble
	boost::numeric::ublas::bounded_matrix<double,3,3> gard_opr = ZeroMatrix(nodes_number,nodes_number);
	gard_opr = tauone * area * prod(DN_DX,trans(DN_DX)); 
	
// 	 array_1d<double,3> nodal_eps = ZeroVector(3);
// 	for (int i = 0; i < nodes_number; i++)
// 	   {
// 	   nodal_eps[i] = GetGeometry()[i].FastGetSolutionStepValue(POROSITY);
// 	   }

	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1) + dof;
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;
			K(row,column) +=   gard_opr(ii,jj);
		   }
	    }
	//build Darcy Au+B|u|u * grad q stabilization term and assemble and assemble
//          double dp = 0.01; //diameter of the particle	
         double kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*dp*dp);

	 const array_1d<double,3>& vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
	 const array_1d<double,3>& vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
	 const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);

	 array_1d<double,3> norm_vel_2 = ZeroVector(nodes_number);
	 array_1d<double,3> vel_gp = ZeroVector(nodes_number);
	 //	vector with the norm^2 of the velocity of the three nodes at the previous iteration;
	 for (int ii = 0; ii < nodes_number; ii++)
	 {
		  norm_vel_2[0] += vel0[ii]*vel0[ii];
		  norm_vel_2[1] += vel1[ii]*vel1[ii];
		  norm_vel_2[2] += vel2[ii]*vel2[ii];
		  vel_gp[ii] = N[0]*vel0[ii] + N[1]*vel1[ii] + N[2]*vel2[ii];
	 }
	 double norm_vel_gp = 0.0;
	 for (int ii = 0; ii < nodes_number; ii++)
	 {
		  norm_vel_gp +=  vel_gp[ii]*vel_gp[ii];
	 }
	
	double fac_linear = tauone * area * kinv * mu ;

	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		double fac_nonlinear = tauone * area * 1.75  * density /eps * sqrt(norm_vel_gp *  kinv / (eps * 150.0));

		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;
			//Au+B|u|u * grad q
			//DARCY TERM linear part
			K(column,row) +=  fac_linear* N[ii] * DN_DX(jj,0)/*/nodal_eps[jj] *nodal_eps[jj]*/ ;
			K(column,row + 1) +=fac_linear* N[ii] * DN_DX(jj,1)/*/nodal_eps[jj] *nodal_eps[jj]*/ ;
			//DARCY TERM nonlinear part 
			K(column,row) +=fac_nonlinear* N[ii] * DN_DX(jj,0)/*/nodal_eps[jj] *nodal_eps[jj]*/ ;
			K(column,row + 1) +=fac_nonlinear* N[ii] * DN_DX(jj,1)/*/nodal_eps[jj] *nodal_eps[jj] */;


		   }
	    }

	//build 1*(grad q) (Fbody ) stabilization term & assemble
	array_1d<double,2> bdf = ZeroVector(2);
	const array_1d<double,2> bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double,2> bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double,2> bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);


	bdf[0] = N[0]*bdf0[0] +  N[1]*bdf1[0] + N[2]*bdf2[0];
	bdf[1] =  N[0]*bdf0[1] +  N[1]*bdf1[1] + N[2]*bdf2[1];

	array_1d<double,3> fbd_stblterm = ZeroVector(nodes_number);
	fbd_stblterm =  tauone * density * area * prod(DN_DX,bdf);

	for(int ii = 0; ii< nodes_number; ++ii)
	  {
		int index = ii*(dof + 1) + dof;
		F[index] += fbd_stblterm[ii] ;
	  }



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
	double dp;
	CalculateDensity(GetGeometry(), density, mu, eps, dp);
	array_1d<double,3> nodal_eps = ZeroVector(3);
	for (int i = 0; i < nodes_number; i++)
	  {
	  nodal_eps[i] = GetGeometry()[i].FastGetSolutionStepValue(POROSITY);
	  }	
// 	 array_1d<double,3> nodal_eps = ZeroVector(3);
// 	for (int i = 0; i < nodes_number; i++)
//     {
// 	   nodal_eps[i] = GetGeometry()[i].FastGetSolutionStepValue(POROSITY);
// 	   }

	//build 1*tau1*ro Nacc grad q)
	double fac = tauone * density * area;
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		int row = ii*(dof+1);
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;

			//K(row,column) += -1*area * fac* N(ii) * DN_DX(jj,0);
			M(column,row) +=  fac* N[ii] * DN_DX(jj,0)/*/nodal_eps[jj]*/;

			//K(row + 1,column) += -1*area * fac* N(ii) * DN_DX(jj,1);
			M(column,row + 1) +=  fac* N[ii] * DN_DX(jj,1)/*/nodal_eps[jj]*/;
		   }
	    }

	KRATOS_CATCH("")

	}

	//************************************************************************************
	//************************************************************************************
	void Fluid2DSplit::AddBodyForceAndMomentum(VectorType& F,const array_1d<double,3>& N, const double time,const double area,const double tauone,const double tautwo)
	{
		KRATOS_TRY
	int nodes_number = 3;
	int dof = 2;
	
	double density;
	double mu;
	double eps;
	double dp;
	CalculateDensity(GetGeometry(), density, mu, eps, dp);
	 
	array_1d<double,3> nodal_eps = ZeroVector(3);	  
	for (int i = 0; i < nodes_number; i++)
	 {
	   nodal_eps[i] = GetGeometry()[i].FastGetSolutionStepValue(POROSITY);
	 }
//  KRATOS_WATCH(nodal_eps)
	 double fac = area * density;
	//body  & momentum term force
	for ( int ii = 0; ii < nodes_number; ii++)
	   {
		int index = ii*(dof + 1) ;
		const array_1d<double,2> bdf = GetGeometry()[ii].FastGetSolutionStepValue(BODY_FORCE);
	
		F[index] += fac * nodal_eps[ii] * N[ii] * bdf[0] ;
		F[index + 1] += fac * nodal_eps[ii] * N[ii] * bdf[1];

	   }
	

	KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	 void Fluid2DSplit::AddVolumeCorrection(VectorType& F,const array_1d<double,3>& N, const double time,const double area)
	   {
		  KRATOS_TRY
		  int nodes_number = 3;
		  int dof = 2;
		  
		  double density;
		  double mu;
		  double eps;
		  double dp;
		  CalculateDensity(GetGeometry(), density, mu, eps, dp);   

		  const array_1d<double,3> vel0_old = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
		  const array_1d<double,3> vel1_old = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
		  const array_1d<double,3> vel2_old = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);

		  double fac = area *density;
		  for ( int ii = 0; ii < nodes_number; ii++)
		      {
			   int row = ii*(dof + 1) + dof ;
			   for( int jj=0; jj < dof; jj++)
			   {		        
				  F[row] -= fac * N(ii) * (DN_DX(0,jj) * vel0_old(jj) + DN_DX(1,jj) * vel1_old(jj)  + DN_DX(2,jj) * vel2_old(jj));
				  //F[row] -= area * density * DN_DX(ii) * (N[0] * vel0_old(jj) + N[1] * vel1_old(jj)  + N[2] * vel2_old(jj));
			   }
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
		UP[index] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY)[0];
		UP[index + 1] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY)[1];
		UP[index + 2] = GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE);
	    }

	F -= prod(K,UP);

	KRATOS_CATCH("")
	}

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
				rResult.resize(0, false);
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

	void Fluid2DSplit::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
	{

	double delta_t= rCurrentProcessInfo[DELTA_TIME];
		
	//getting data for the given geometry
	double Area;
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);
// /*	double tauone =0.0;// delta_t;
// 	double tautwo = 0.0;//delta_t;*/

       double tauone, tautwo;
	CalculateTau(tauone, tautwo, delta_t, Area,  rCurrentProcessInfo);
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
	void Fluid2DSplit::CalculateDensity(Geometry< Node<3> > geom, double& elemental_density, double& elemental_viscosity, double& elemental_porosity, double& elemental_diameter)
 	{//PAY ATTENTION: CALCULATION OF ELEMENTAL DENSITY WITHOUT THE EFFECT OF POROSITY : RHO === RHO_WATER

//*****begin provisionary_comment to main porosity*************************
// 	       //Check if some of the elements don't have a porosity assigned
// 	       double eps0 = geom[0].FastGetSolutionStepValue(POROSITY);
// 	       if(eps0 == 0.0)
// 		{   eps0 = 1.0;}
// 	       double eps1 = geom[1].FastGetSolutionStepValue(POROSITY);
// 	       if(eps1 == 0.0)
// 		{    eps1 = 1.0;}
// 	       double eps2 = geom[2].FastGetSolutionStepValue(POROSITY);
// 	       if(eps2 == 0.0)
// 		{    eps2 = 1.0;}
// 	       //Check if some of the elements don't have a DIAMETER assigned
// 	       double d0 = geom[0].FastGetSolutionStepValue(DIAMETER);
// 	       if(d0 == 0.0)
// 		{    d0 = 1.0;}
// 	       double d1 = geom[1].FastGetSolutionStepValue(DIAMETER);
// 	       if(d1 == 0.0)
// 		{    d1 = 1.0;}
// 	       double d2 = geom[2].FastGetSolutionStepValue(DIAMETER);
// 	       if(d2 == 0.0)
// 		{    d2 = 1.0;}
// 
// 	       elemental_density = 0.0;
// 	       elemental_porosity = 1.0;
// 	       elemental_diameter = 0.01;
// 
// 	       if(eps0 == eps1 && eps1 == eps2)
// 		{
// 		        //for inside the domain totally inside one fluid
// 		        elemental_porosity = eps0;
// 		        elemental_density = geom[0].FastGetSolutionStepValue(DENSITY); //* eps0;
// 		        elemental_diameter = geom[0].FastGetSolutionStepValue(DIAMETER); 
// 		        elemental_viscosity = geom[0].FastGetSolutionStepValue(VISCOSITY) * elemental_density;	//mu = nu * density //we assigne nu=1E-6 from Gid 
// 
// 		}
// 	       else if(eps0 == eps1)
// 		{
// 		        elemental_porosity = eps0;
// 		        elemental_density = geom[0].FastGetSolutionStepValue(DENSITY); // * eps0;	
// 		        elemental_diameter = geom[0].FastGetSolutionStepValue(DIAMETER); 
// 		        elemental_viscosity = geom[0].FastGetSolutionStepValue(VISCOSITY) * elemental_density;	//mu = nu * density 
// 		}
// 	       else if(eps1 == eps2)
// 		{
// 		        elemental_porosity = eps1;
// 		        elemental_density = geom[1].FastGetSolutionStepValue(DENSITY); // * eps1;	
// 		        elemental_diameter = geom[1].FastGetSolutionStepValue(DIAMETER); 
// 		        elemental_viscosity = geom[1].FastGetSolutionStepValue(VISCOSITY) * elemental_density;  //mu = nu * density 
// 		}
// 	       else if(eps2 == eps0)
// 		{
// 		        elemental_porosity = eps2;
// 		        elemental_density = geom[2].FastGetSolutionStepValue(DENSITY); // * eps2;
// 		        elemental_diameter = geom[2].FastGetSolutionStepValue(DIAMETER); 
// 		        elemental_viscosity = geom[2].FastGetSolutionStepValue(VISCOSITY)* elemental_density;  //mu = nu * density 
// 		}
// 	       else { KRATOS_WATCH("ERROR!!! three different values of densities");}
//*****end provisionary_comment to main porosity *************************
//*****begin provisionary_comment to average porosity *************************
	       elemental_porosity = 0.0;
	       elemental_density = 0.0;
	       elemental_diameter = 0.0;
	       elemental_viscosity = 0.0;
	       unsigned int count = 0;
	       for(unsigned int i=0; i< 3; i++)
	       {
		        count++;
// 		        elemental_porosity += geom[i].FastGetSolutionStepValue(DENSITY);
		        elemental_density += geom[i].FastGetSolutionStepValue(DENSITY); //* eps0;
		        elemental_diameter += geom[i].FastGetSolutionStepValue(DIAMETER); 
		        elemental_viscosity += geom[i].FastGetSolutionStepValue(VISCOSITY) * elemental_density;	//mu = nu * density //we assigne nu=1E-6 from Gid 
	       }
	       if (count == 0.0)
		    { KRATOS_WATCH("*********************************************************ERROR-> Calculate Density Fluid_2d_split *******************")}
	       else
 		    { 
// 		      elemental_porosity /= count;
		      elemental_density /= count;
		      elemental_diameter /= count;
		      elemental_viscosity /= count;
		    }
//*****end provisionary_comment to average porosity *************************
//*****begin provisionary_comment to dominant porosity *************************
		if(     geom[0].FastGetSolutionStepValue(POROSITY) == 1.0 )
		{    elemental_porosity = geom[0].FastGetSolutionStepValue(POROSITY);}
		else if(geom[1].FastGetSolutionStepValue(POROSITY) == 1.0)
		{    elemental_porosity = geom[1].FastGetSolutionStepValue(POROSITY);}
		else if(geom[2].FastGetSolutionStepValue(POROSITY) == 1.0)
		{    elemental_porosity = geom[2].FastGetSolutionStepValue(POROSITY);}
		else
		{ elemental_porosity = geom[0].FastGetSolutionStepValue(POROSITY);}

//*****end provisionary_comment to dominant porosity *************************
 	}

	//*************************************************************************************
	//*************************************************************************************

		void Fluid2DSplit::CalculateTau(double& tauone, double& tautwo, const double time,const double area, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
	//calculate mean advective velocity and taus
	const array_1d<double,3>& adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
	const array_1d<double,3>& adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
	const array_1d<double,3>& adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);




	ms_adv_vel[0] = N[0]*(adv_vel0[0])+N[1]*(adv_vel1[0])+N[2]*(adv_vel2[0]);
	ms_adv_vel[1] = N[0]*(adv_vel0[1])+N[1]*(adv_vel1[1])+N[2]*(adv_vel2[1]);
									                                                                                                                                                                                                  
	double advvel_norm = ms_adv_vel[0]*ms_adv_vel[0]+ms_adv_vel[1]*ms_adv_vel[1];
	advvel_norm = sqrt(advvel_norm);
	
	double ele_length = 2.0*sqrt(area/3.00);
	
	double mu;

	double density;
	double eps;
	double dp;
	CalculateDensity(GetGeometry(), density, mu, eps, dp);

	double advvel_norm_fluid = advvel_norm/eps;
//         	double dp = 0.01; //diameter of the particle	
	double kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*dp*dp);

	double fac_linear = eps * kinv * mu /density;	
	double fac_nonlinear = (1.75 * advvel_norm * sqrt(  kinv / (eps * 150.0))) ;     
/*	double fac_linear = kinv * mu /density;	
	double fac_nonlinear = (1.75 / eps * advvel_norm * sqrt(  kinv / (eps * 150.0))) ;  */   

	int dyn_st_switch = rCurrentProcessInfo[DYNAMIC_TAU];
	
		if(dyn_st_switch)
		  {
			tauone = 1.0/(1.0/time + 4.0*mu/(ele_length*ele_length*density)+2.0*advvel_norm_fluid*1.0/ele_length + fac_linear +  fac_nonlinear);
		  }
		else
		 {
			tauone = 1.0/(4.0*mu/(ele_length*ele_length*density)+2.0*advvel_norm_fluid*1.0/ele_length + fac_linear +  fac_nonlinear);
		  }
		
	tautwo = mu/density + 1.0*ele_length*advvel_norm_fluid/2.0;


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


