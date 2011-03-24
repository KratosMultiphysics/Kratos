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
//   Date:                $Date: 2009-01-21 14:15:02 $
//   Revision:            $Revision: 1.6 $
//
//

// #define NAVIERSTOKES //if not STOKES is solved
#define EXPONENCIAL_MODEL // if not -> BILINEAR_MODEL is calculated
// #define COMPRESSIBLE_MODEL_ML__PN1_PN__
// #define COMPRESSIBLE_MODEL_MLPN1_MCPN

// #define K0 //fixed tangent method
// #define Kvisc0 //fixed tangent method only for viscous terms

// System includes 

// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/nonewtonian_asgs_2d.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos {

    //************************************************************************************
    //************************************************************************************

    NoNewtonianASGS2D::NoNewtonianASGS2D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {
	//DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    NoNewtonianASGS2D::NoNewtonianASGS2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {

    }

    Element::Pointer NoNewtonianASGS2D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {

	KRATOS_TRY
	return Element::Pointer(new NoNewtonianASGS2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	KRATOS_CATCH("");
    }

    NoNewtonianASGS2D::~NoNewtonianASGS2D() {
    }

    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
	KRATOS_TRY
	
	boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
	array_1d<double, 3 > N;
	array_1d<double, 2 > ms_adv_vel;
	
	// KRATOS_WATCH("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ NoNewtonianASGS element ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	int nodes_number = 3;
	int dim = 2;
	unsigned int matsize = nodes_number * (dim + 1);

	if (rLeftHandSideMatrix.size1() != matsize)
	    rLeftHandSideMatrix.resize(matsize, matsize); //false says not to preserve existing storage!!

	if (rRightHandSideVector.size() != matsize)
	    rRightHandSideVector.resize(matsize); //false says not to preserve existing storage!!


	noalias(rLeftHandSideMatrix) = ZeroMatrix(matsize, matsize);
	noalias(rRightHandSideVector) = ZeroVector(matsize);

	double delta_t = rCurrentProcessInfo[DELTA_TIME];



	//getting data for the given geometry
	double Area;
	GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

	double tauone;
	double tautwo;
	CalculateTau(DN_DX,N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);


	//add body force and momentum
	AddBodyForceAndMomentum(rRightHandSideVector,DN_DX, N, delta_t, Area, tauone, tautwo);

// KRATOS_WATCH(rCurrentProcessInfo[OSS_SWITCH])

// 	      //add projections
 	      if (rCurrentProcessInfo[OSS_SWITCH] == 1)
		  AddProjectionForces(rRightHandSideVector, DN_DX, Area, tauone, tautwo);


	KRATOS_CATCH("")
    }
    //***********************************************************************************++
    //**************************************************************************************++

    void NoNewtonianASGS2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
	KRATOS_TRY

	MatrixType temp = Matrix();
	CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);

	KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateMassContribution(MatrixType& K, const double time, const double area) {
	KRATOS_TRY
		double lump_mass_fac = area * 0.333333333333333333333333;
	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	int nodes_number = 3;
	int dof = 2;
	for (int nd = 0; nd < nodes_number; nd++) {
	    int row = nd * (dof + 1);
	    for (int jj = 0; jj < dof; jj++)
		K(row + jj, row + jj) += density * lump_mass_fac;
	}

	KRATOS_CATCH("")

    }
    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {
	KRATOS_TRY

		//lumped
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
	unsigned int NumberOfNodes = GetGeometry().size();
	unsigned int MatSize = (dimension + 1) * NumberOfNodes;
	if (rMassMatrix.size1() != MatSize)
	    rMassMatrix.resize(MatSize, MatSize, false);

	rMassMatrix = ZeroMatrix(MatSize, MatSize);
	double delta_t = rCurrentProcessInfo[DELTA_TIME];

	boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
	array_1d<double, 3 > N;
	
	//getting data for the given geometry
	double Area;
	GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

	//Calculate tau
	double tauone;
	double tautwo;
	CalculateTau(DN_DX,N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);

	CalculateMassContribution(rMassMatrix, delta_t, Area);
/*20101216*/
	/*Stablization*/
	//add stablilization terms due to advective term (a)grad(V) * ro*Acce
// #ifdef NAVIERSTOKES
	CalculateAdvMassStblTerms(rMassMatrix, DN_DX, N, tauone, Area);
// #endif
	//add stablilization terms due to grad term grad(q) * ro*Acce
	CalculateGradMassStblTerms(rMassMatrix, DN_DX,N, tauone, Area);

	KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateLocalVelocityContribution(MatrixType& rDampMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
	KRATOS_TRY
	int nodes_number = 3;
	int dim = 2;
	unsigned int matsize = nodes_number * (dim + 1);

	if (rDampMatrix.size1() != matsize)
	    rDampMatrix.resize(matsize, matsize, false); //false says not to preserve existing storage!!
	    
	boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
	array_1d<double, 3 > N;
//         array_1d<double, 2 > ms_adv_vel;

	noalias(rDampMatrix) = ZeroMatrix(matsize, matsize);

	double delta_t = rCurrentProcessInfo[DELTA_TIME];
	unsigned int it_num= rCurrentProcessInfo[NL_ITERATION_NUMBER];
	const double m_coef = rCurrentProcessInfo[M];
	
	//getting data for the given geometry
	double Area;
	GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);


	/*         LHS           */
	/*Advective term*/
	double tauone;
	double tautwo;
	CalculateTau(DN_DX,N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);
// #ifdef NAVIERSTOKES
	CalculateAdvectiveTerm(rDampMatrix, DN_DX,N, tauone, tautwo, delta_t, Area);
// #endif
	/*Calculate Pressure term + divergence term of pressure equation*/
	CalculatePressureTerm(rDampMatrix, DN_DX, N, delta_t, Area);

	//compute projections
	/*Stablization*/
	//stabilization terms
/*20101216*/
	CalculateDivStblTerm(rDampMatrix, DN_DX, tautwo, Area);//tau 2
// #ifdef NAVIERSTOKES
	CalculateAdvStblAllTerms(rDampMatrix, rRightHandSideVector, DN_DX, N, tauone, delta_t, Area);
// #endif
	CalculateGradStblAllTerms(rDampMatrix, rRightHandSideVector, DN_DX,N, delta_t, tauone, Area);
	//KRATOS_WATCH(rRightHandSideVector);

	/*         RHS           */
	/*Internal Forces*/
	CalculateResidual(rDampMatrix, rRightHandSideVector,DN_DX, m_coef, Area);

	/*         LHS           */
	/*Viscous term*/

	CalculateViscousTerm(rDampMatrix,  DN_DX, it_num, m_coef, Area);

#ifdef K0 //Fixed tangent method 
KRATOS_WATCH("Fixed tangent method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	mlhs0.resize(matsize,matsize);
// 	unsigned int it_num= rCurrentProcessInfo[NL_ITERATION_NUMBER];
// KRATOS_WATCH(it_num);
// KRATOS_WATCH(rDampMatrix);

	if(it_num == 1)	
	  noalias(mlhs0) = rDampMatrix;
	else
	  rDampMatrix = mlhs0;
// KRATOS_WATCH(rDampMatrix);
// KRATOS_WATCH(mlhs0);
#endif	
	KRATOS_CATCH("")
    }


    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateViscousTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX,  const int it_num, const double m_coef, const double area) {
	KRATOS_TRY
	double mu;
	int nodes_number = 3;
	int dof = 2;
	double density;
	calculatedensity(GetGeometry(), density, mu);
	boost::numeric::ublas::bounded_matrix<double, 3, 6 > B = ZeroMatrix(3, 6);
	boost::numeric::ublas::bounded_matrix<double, 3, 3 > C = ZeroMatrix(3, 3);
	boost::numeric::ublas::bounded_matrix<double, 6, 6 > temp; 
	
// 	unsigned int it_num= rCurrentProcessInfo[NL_ITERATION_NUMBER];
// 	if(it_num == 1){
		temp = ZeroMatrix(6, 6);
		double app_mu;
		double app_mu_derivative;
		double gamma_dot;
		array_1d<double, 3 > grad_sym_vel = ZeroVector(3);
			
		//calculating operator B
		CalculateB(B, DN_DX);
		// KRATOS_WATCH(B)
		
		//Bingham Fluid
	      CalculateApparentViscosity(app_mu, app_mu_derivative, grad_sym_vel, gamma_dot, B, mu, m_coef);
// 	      CalculateNodalApparentViscosity(grad_sym_vel, gamma_dot, B, mu, m_coef);
		//Newtonian Fluid: leave decommented the CalculateApparentViscosity (we need grad_sym_vel) and decomment the following line
		// Remember to modify CalculateResidualand CalculateTau.
	// 	app_mu = mu;
		
		C(0, 0) = 2.0- 2.0/3.0;
		C(0, 1) = 0.0- 2.0/3.0;
		C(0, 2) = 0.0;
		C(1, 0) = 0.0- 2.0/3.0;
		C(1, 1) = 2.0- 2.0/3.0;
		C(1, 2) = 0.0;
		C(2, 0) = 0.0;
		C(2, 1) = 0.0;
		C(2, 2) = 1.0;
		
		C *= app_mu;

	// 	//PAY ATTENTION : 3rd component of B*u divided by two in order to recover eps12 or eps21 that is the same.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// 	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// 		grad_sym_vel[2] *= 0.5;
// 		double aux_coeff = 4.0 * app_mu_derivative / gamma_dot;
// 		for(unsigned int i=0 ; i < nodes_number; i++){
// 		  for(unsigned int j=0; j<nodes_number; j++){
// 		    C(i,j) += grad_sym_vel[i] * grad_sym_vel[j] * aux_coeff;
// 		  }
// 		}
		  

		
	//  	if(gamma_dot > 1e-3){
	// 	double aux_coeff = 4.0 * app_mu_derivative / gamma_dot;
	// 	C(0, 0) +=  grad_sym_vel[0] * grad_sym_vel[0] * aux_coeff;
	// 	C(0, 1) +=  grad_sym_vel[0] * grad_sym_vel[1] * aux_coeff;
	// 	C(0, 2) +=  grad_sym_vel[0] * grad_sym_vel[2] * aux_coeff * 0.5;
	// 	C(1, 0) +=  grad_sym_vel[1] * grad_sym_vel[0] * aux_coeff;
	// 	C(1, 1) +=  grad_sym_vel[1] * grad_sym_vel[1] * aux_coeff;
	// 	C(1, 2) +=  grad_sym_vel[1] * grad_sym_vel[2] * aux_coeff * 0.5;
	// 	C(2, 0) +=  grad_sym_vel[2] * grad_sym_vel[0] * aux_coeff * 0.5;
	// 	C(2, 1) +=  grad_sym_vel[2] * grad_sym_vel[1] * aux_coeff * 0.5;
	// 	C(2, 2) +=  grad_sym_vel[2] * grad_sym_vel[2] * aux_coeff * 0.25;
	// 		C(0, 0) +=  4.0 *  grad_sym_vel[0] * grad_sym_vel[0] *app_mu_derivative / gamma_dot;
	// 		C(0, 1) +=  4.0 *  grad_sym_vel[0] * grad_sym_vel[1] *app_mu_derivative / gamma_dot;
	// 		C(0, 2) +=  2.0 *  grad_sym_vel[0] * grad_sym_vel[2] *app_mu_derivative / gamma_dot;
	// 		C(1, 0) +=  4.0 *  grad_sym_vel[1] * grad_sym_vel[0] *app_mu_derivative / gamma_dot;
	// 		C(1, 1) +=  4.0 *  grad_sym_vel[1] * grad_sym_vel[1] *app_mu_derivative / gamma_dot;
	// 		C(1, 2) +=  2.0 *  grad_sym_vel[1] * grad_sym_vel[2] *app_mu_derivative / gamma_dot;
	// 		C(2, 0) +=  2.0 *  grad_sym_vel[2] * grad_sym_vel[0] *app_mu_derivative / gamma_dot;
	// 		C(2, 1) +=  2.0*  grad_sym_vel[2] * grad_sym_vel[1] *app_mu_derivative / gamma_dot;
	// 		C(2, 2) +=       grad_sym_vel[2] * grad_sym_vel[2] *app_mu_derivative / gamma_dot;
	//  	}

		// KRATOS_WATCH(C)
		//Calculating the viscous contribution to the LHS int(Btrans C B)dA
		//         temp = prod(trans(B),prod(C,B));
		for (unsigned int i = 0; i < B.size2(); i++) {
		    for (unsigned int j = 0; j < B.size2(); j++) {
			for (unsigned int l = 0; l < C.size1(); l++) {
			    for (unsigned int k = 0; k < C.size1(); k++) {
				temp(i, j) += B(l, i) * C(l, k) * B(k, j);
			    }
			}
		    }
		}
// 	}

#ifdef Kvisc0 //fixed tangent method only for viscous contribution
// 	unsigned int it_num= rCurrentProcessInfo[NL_ITERATION_NUMBER];
//KRATOS_WATCH("Fixed tangent method only for viscous term~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
// KRATOS_WATCH(it_num);
// KRATOS_WATCH(temp);
	if(it_num == 1)	
	  mKvisc0 = temp;
	else
	  temp = mKvisc0;
// KRATOS_WATCH(temp);
#endif
	// KRATOS_WATCH(temp)
	temp *= area;
	for (int ii = 0; ii < nodes_number; ii++) {
	    int row = ii * (dof + 1);
	    int loc_row = ii*dof;
	    for (int jj = 0; jj < nodes_number; jj++) {
		int column = jj * (dof + 1);
		int loc_column = jj*dof;

		K(row, column) += temp(loc_row, loc_column);
		K(row, column + 1) += temp(loc_row, loc_column + 1);
		K(row + 1, column + 1) += temp(loc_row + 1, loc_column + 1);
		K(row + 1, column) += temp(loc_row + 1, loc_column);
	    }
	}

	KRATOS_CATCH("")
    }
//     void NoNewtonianASGS2D::CalculateViscousTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX,  const int it_num, const double m_coef, const double area) {
// 	KRATOS_TRY
// 	double mu;
// 	int nodes_number = 3;
// 	int dof = 2;
// 	double density;
// 	calculatedensity(GetGeometry(), density, mu);
// 	boost::numeric::ublas::bounded_matrix<double, 3, 6 > B = ZeroMatrix(3, 6);
// 	boost::numeric::ublas::bounded_matrix<double, 3, 3 > C = ZeroMatrix(3, 3);
// 	boost::numeric::ublas::bounded_matrix<double, 6, 6 > temp; 
// 	
// 
// 		temp = ZeroMatrix(6, 6);
// 		array_1d<double, 3 >  app_mu = = ZeroVector(3);
// 		double app_mu_derivative;
// 		double gamma_dot;
// 		array_1d<double, 3 > grad_sym_vel = ZeroVector(3);
// 			
// 		//calculating operator B
// 		CalculateB(B, DN_DX);
// 		// KRATOS_WATCH(B)
// 		
// 		//Bingham Fluid
// 	      CalculateApparentViscosity(app_mu, app_mu_derivative, grad_sym_vel, gamma_dot, B, mu, m_coef);
// 		//Newtonian Fluid: leave decommented the CalculateApparentViscosity (we need grad_sym_vel) and decomment the following line
// 		// Remember to modify CalculateResidualand CalculateTau.
// 	// 	app_mu = mu;
// 		for(unsigned int ii = 0.0, ii < nodes_number, ii++){
// 		  C(0, 0) = 2.0- 2.0/3.0;
// 		  C(0, 1) = 0.0- 2.0/3.0;
// 		  C(0, 2) = 0.0;
// 		  C(1, 0) = 0.0- 2.0/3.0;
// 		  C(1, 1) = 2.0- 2.0/3.0;
// 		  C(1, 2) = 0.0;
// 		  C(2, 0) = 0.0;
// 		  C(2, 1) = 0.0;
// 		  C(2, 2) = 1.0;
// 		
// 		  C *= app_mu[ii];
// 
// 		  // KRATOS_WATCH(C)
// 		  //Calculating the viscous contribution to the LHS int(Btrans C B)dA
// 		  //         temp = prod(trans(B),prod(C,B));
// 		  for (unsigned int i = 0; i < B.size2(); i++) {
// 		      for (unsigned int j = 0; j < B.size2(); j++) {
// 			  for (unsigned int l = 0; l < C.size1(); l++) {
// 			      for (unsigned int k = 0; k < C.size1(); k++) {
// 				  temp(i, j) += B(l, i) * C(l, k) * B(k, j);
// 			      }
// 			  }
// 		      }
// 		  }
// 		}
// 
// 	// KRATOS_WATCH(temp)
// 	temp *= area * 0.33333333333333333333333;
// 	for (int ii = 0; ii < nodes_number; ii++) {
// 	    int row = ii * (dof + 1);
// 	    int loc_row = ii*dof;
// 	    for (int jj = 0; jj < nodes_number; jj++) {
// 		int column = jj * (dof + 1);
// 		int loc_column = jj*dof;
// 
// 		K(row, column) += temp(loc_row, loc_column);
// 		K(row, column + 1) += temp(loc_row, loc_column + 1);
// 		K(row + 1, column + 1) += temp(loc_row + 1, loc_column + 1);
// 		K(row + 1, column) += temp(loc_row + 1, loc_column);
// 	    }
// 	}
// 
// 	KRATOS_CATCH("")
//     }
    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateAdvectiveTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double, 3 >& N, const double tauone, const double tautwo, const double time, const double area) {
	KRATOS_TRY
	
	array_1d<double, 2 > ms_adv_vel;

	//calculate mean advective velocity and taus
	const array_1d<double, 3 > & adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double, 3 > & adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double, 3 > & adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);

	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	ms_adv_vel[0] = N[0]*(adv_vel0[0] - mesh_vel0[0]) + N[1]*(adv_vel1[0] - mesh_vel1[0]) + N[2]*(adv_vel2[0] - mesh_vel2[0]);
	ms_adv_vel[1] = N[0]*(adv_vel0[1] - mesh_vel0[1]) + N[1]*(adv_vel1[1] - mesh_vel1[1]) + N[2]*(adv_vel2[1] - mesh_vel2[1]);

	//calculate convective term
	int nodes_number = 3;
	int dof = 2;
	int matsize = dof*nodes_number;

	boost::numeric::ublas::bounded_matrix<double, 2, 6 > conv_opr = ZeroMatrix(dof, matsize);
	boost::numeric::ublas::bounded_matrix<double, 6, 2 > shape_func = ZeroMatrix(matsize, dof);

	for (int ii = 0; ii < nodes_number; ii++) {
	    int column = ii*dof;
	    conv_opr(0, column) = DN_DX(ii, 0) * ms_adv_vel[0] + DN_DX(ii, 1) * ms_adv_vel[1];
	    conv_opr(1, column + 1) = conv_opr(0, column);

	    shape_func(column, 0) = N[ii];
	    shape_func(column + 1, 1) = shape_func(column, 0);
	}
	boost::numeric::ublas::bounded_matrix<double, 6, 6 > temp_convterm = ZeroMatrix(matsize, matsize);
	temp_convterm = prod(shape_func, conv_opr);

	for (int ii = 0; ii < nodes_number; ii++) {
	    int row = ii * (dof + 1);
	    int loc_row = ii*dof;
	    for (int jj = 0; jj < nodes_number; jj++) {
		int column = jj * (dof + 1);
		int loc_column = jj*dof;

		K(row, column) += area * density * temp_convterm(loc_row, loc_column);
		K(row + 1, column + 1) += area * density * temp_convterm(loc_row + 1, loc_column + 1);
	    }
	}

	KRATOS_CATCH("")

    }
    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculatePressureTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double, 3 > & N, const double time, const double area) {
	KRATOS_TRY
	int nodes_number = 3;
	int dof = 2;

	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);


	for (int ii = 0; ii < nodes_number; ii++) {
	    int row = ii * (dof + 1);
	    for (int jj = 0; jj < nodes_number; jj++) {
		int column = jj * (dof + 1) + dof;
		//**************************************************************
		//Elemental gradient of pressure term (momentum equation)
		K(row, column) -= area * N(jj) * DN_DX(ii, 0);
		K(row + 1, column) -= area * N(jj) * DN_DX(ii, 1);
		//**************************************************************
		//Elemental divergence terms (continuity equation)
		//Fomulation n1  int( q * rho * Div( u ))
		K(column, row) += area * density * N(jj) * DN_DX(ii, 0);
		K(column, row + 1) += area * density * N(jj) * DN_DX(ii, 1);
	    }
	}
#ifdef COMPRESSIBLE_MODEL_ML__PN1_PN__
	// 	////compressibility term & assemble Lumped Mass Matrix
	double c2_inv = 1e-5;
	double compressibility_coef = 0.33333333333333333333333 * area * c2_inv /time;

	for (int ii = 0; ii < nodes_number; ii++) {
	    int row = ii * (dof + 1) + dof;

		//**************************************************************
		//Elemental compressibility terms (continuity equation)
		// int( q * rho * 1/(rho c^2)(p)
		K(row, row) += compressibility_coef ;

	 }
#endif	
#ifdef COMPRESSIBLE_MODEL_MLPN1_MCPN
	// 	////compressibility term & assemble Lumped Mass Matrix
	double c2_inv = 1e-5;
	double compressibility_coef = 0.33333333333333333333333 * area * c2_inv /time;

	for (int ii = 0; ii < nodes_number; ii++) {
	    int row = ii * (dof + 1) + dof;

		//**************************************************************
		//Elemental compressibility terms (continuity equation)
		// int( q * rho * 1/(rho c^2)(p)
		K(row, row) += compressibility_coef ;

	 }
#endif
// //compressibility term & assemble Consistent Mass Matrix (Attention: Rank zero matrix, it may couse problems!)	
// 	double c2_inv = 1e-6;
// 	double compressibility_coef = area * c2_inv /time;
// // KRATOS_WATCH("time")
// // KRATOS_WATCH(time)
// // KRATOS_WATCH("compressibility_coef LHS")
// // KRATOS_WATCH(compressibility_coef)
// // KRATOS_WATCH("K before")
// // KRATOS_WATCH(K)
// 	for (int ii = 0; ii < nodes_number; ii++) {
// 	    int row = ii * (dof + 1) + dof;
// 
// 	    for (int jj = 0; jj < nodes_number; jj++) {
// 		int column = jj * (dof + 1) + dof;
// 		//**************************************************************
// 		//Elemental compressibility terms (continuity equation)
// 		// int( q * rho * 1/(rho c^2)(p)
// 		K(row, column) += compressibility_coef * N(ii) * N(jj);
// 
// 	    }
// 	}
// KRATOS_WATCH("K after")
// KRATOS_WATCH(K)



	KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateDivStblTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const double tautwo, const double area) {
	KRATOS_TRY
	int nodes_number = 3;
	int dof = 2;
	int matsize = dof*nodes_number;

	boost::numeric::ublas::bounded_matrix<double, 1, 6 > div_opr = ZeroMatrix(1, matsize);
	for (int ii = 0; ii < nodes_number; ii++) {
	    int index = dof*ii;
	    div_opr(0, index) = DN_DX(ii, 0);
	    div_opr(0, index + 1) = DN_DX(ii, 1);
	}


	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);


	boost::numeric::ublas::bounded_matrix<double, 6, 6 > temp_div = ZeroMatrix(matsize, matsize);
	temp_div = tautwo * prod(trans(div_opr), div_opr);

	for (int ii = 0; ii < nodes_number; ii++) {
	    int row = ii * (dof + 1);
	    int loc_row = ii*dof;
	    for (int jj = 0; jj < nodes_number; jj++) {
		int column = jj * (dof + 1);
		int loc_column = jj*dof;

		K(row, column) += area * density * temp_div(loc_row, loc_column);
		K(row, column + 1) += area * density * temp_div(loc_row, loc_column + 1);
		K(row + 1, column) += area * density * temp_div(loc_row + 1, loc_column);
		K(row + 1, column + 1) += area * density * temp_div(loc_row + 1, loc_column + 1);
	    }
	}

	KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateAdvStblAllTerms(MatrixType& K, VectorType& F, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double, 3 > & N, const double tauone, const double time, const double area) {
	KRATOS_TRY
		const array_1d<double, 3 > & adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double, 3 > & adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double, 3 > & adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);

	array_1d<double,2> ms_adv_vel;
	ms_adv_vel[0] = N[0]*(adv_vel0[0] - mesh_vel0[0]) + N[1]*(adv_vel1[0] - mesh_vel1[0]) + N[2]*(adv_vel2[0] - mesh_vel2[0]);
	ms_adv_vel[1] = N[0]*(adv_vel0[1] - mesh_vel0[1]) + N[1]*(adv_vel1[1] - mesh_vel1[1]) + N[2]*(adv_vel2[1] - mesh_vel2[1]);

	//ms_adv_vel[0] = 0.0;
	//ms_adv_vel[1] = 0.0;

	//calculate convective term
	int nodes_number = 3;
	int dof = 2;
	int matsize = dof*nodes_number;

	boost::numeric::ublas::bounded_matrix<double, 2, 6 > conv_opr = ZeroMatrix(dof, matsize);
	boost::numeric::ublas::bounded_matrix<double, 6, 2 > shape_func = ZeroMatrix(matsize, dof);

	for (int ii = 0; ii < nodes_number; ii++) {
	    int column = ii*dof;
	    conv_opr(0, column) = DN_DX(ii, 0) * ms_adv_vel[0] + DN_DX(ii, 1) * ms_adv_vel[1];
	    conv_opr(1, column + 1) = conv_opr(0, column);

	    shape_func(column, 0) = N[ii];
	    shape_func(column + 1, 1) = shape_func(column, 0);
	}
	//build (a.grad V)(ro*a.grad U) stabilization term & assemble
	boost::numeric::ublas::bounded_matrix<double, 6, 6 > adv_stblterm = ZeroMatrix(matsize, matsize);
	adv_stblterm = tauone * prod(trans(conv_opr), conv_opr);


	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	for (int ii = 0; ii < nodes_number; ii++) {
	    int row = ii * (dof + 1);
	    int loc_row = ii*dof;
	    for (int jj = 0; jj < nodes_number; jj++) {
		int column = jj * (dof + 1);
		int loc_column = jj*dof;

		K(row, column) += area * density * adv_stblterm(loc_row, loc_column);
		K(row, column + 1) += area * density * adv_stblterm(loc_row, loc_column + 1);
		K(row + 1, column) += area * density * adv_stblterm(loc_row + 1, loc_column);
		K(row + 1, column + 1) += area * density * adv_stblterm(loc_row + 1, loc_column + 1);
	    }
	}
	//build 1*tau1*(a.grad V)(grad P) & 1*tau1*(grad q)(ro*a.grad U) stabilization terms & assemble
	boost::numeric::ublas::bounded_matrix<double, 6, 3 > grad_stblterm = ZeroMatrix(matsize, nodes_number);
	grad_stblterm = tauone * prod(trans(conv_opr), trans(DN_DX));

	for (int ii = 0; ii < nodes_number; ii++) {
	    int row = ii * (dof + 1);
	    int loc_row = ii*dof;
	    for (int jj = 0; jj < nodes_number; jj++) {
		int column = jj * (dof + 1) + dof;

		K(row, column) += area * grad_stblterm(loc_row, jj);
		K(row + 1, column) += area * grad_stblterm(loc_row + 1, jj);

		K(column, row) += area * density * grad_stblterm(loc_row, jj);
		K(column, row + 1) += area * density * grad_stblterm(loc_row + 1, jj);
	    }
	}


	//build (1.0*a.grad V) (Fbody) stabilization term & assemble
	array_1d<double, 2 > bdf = ZeroVector(2);
	const array_1d<double, 3 > & bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double, 3 > & bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double, 3 > & bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);


	bdf[0] = N[0]*(density * bdf0[0]) + N[1]*(density * bdf1[0]) + N[2]*(density * bdf2[0]);
	bdf[1] = N[0]*(density * bdf0[1]) + N[1]*(density * bdf1[1]) + N[2]*(density * bdf2[1]);


	array_1d<double, 6 > fbd_stblterm = ZeroVector(matsize);
	fbd_stblterm = tauone * prod(trans(conv_opr), bdf);

	for (int ii = 0; ii < nodes_number; ++ii) {
	    int index = ii * (dof + 1);
	    int loc_index = ii*dof;
	    F[index] += area * fbd_stblterm[loc_index];
	    F[index + 1] += area * fbd_stblterm[loc_index + 1];
	}
	
// 	//build (a.grad V)(ro*grad nu.grad U) stabilization term & assemble
// 	boost::numeric::ublas::bounded_matrix<double, 2, 6 > visc_opr = ZeroMatrix(dof, matsize);
// 
// 	for (int ii = 0; ii < nodes_number; ii++) {
// 	    int column = ii*dof;
// 	    visc_opr(0, column) = DN_DX(ii, 0) * ms_adv_vel[0] + DN_DX(ii, 1) * ms_adv_vel[1];
// 	    visc_opr(1, column + 1) = visc_opr(0, column);
// 
// 	}
// 	boost::numeric::ublas::bounded_matrix<double, 6, 6 > visc_stblterm = ZeroMatrix(matsize, matsize);
// 	visc_stblterm = tauone * prod(trans(conv_opr), visc_opr);
// 
// 
// 	double density;
// 	double mu;
// 	calculatedensity(GetGeometry(), density, mu);
// 
// 	for (int ii = 0; ii < nodes_number; ii++) {
// 	    int row = ii * (dof + 1);
// 	    int loc_row = ii*dof;
// 	    for (int jj = 0; jj < nodes_number; jj++) {
// 		int column = jj * (dof + 1);
// 		int loc_column = jj*dof;
// 
// 		K(row, column) += area * density * adv_stblterm(loc_row, loc_column);
// 		K(row, column + 1) += area * density * adv_stblterm(loc_row, loc_column + 1);
// 		K(row + 1, column) += area * density * adv_stblterm(loc_row + 1, loc_column);
// 		K(row + 1, column + 1) += area * density * adv_stblterm(loc_row + 1, loc_column + 1);
// 	    }
// 	}
	KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateAdvMassStblTerms(MatrixType& M, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double, 3 > & N, const double tauone, const double area) {
	KRATOS_TRY
		const array_1d<double, 3 > & adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double, 3 > & adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double, 3 > & adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);

	array_1d<double,2> ms_adv_vel;
	ms_adv_vel[0] = N[0]*(adv_vel0[0] - mesh_vel0[0]) + N[1]*(adv_vel1[0] - mesh_vel1[0]) + N[2]*(adv_vel2[0] - mesh_vel2[0]);
	ms_adv_vel[1] = N[0]*(adv_vel0[1] - mesh_vel0[1]) + N[1]*(adv_vel1[1] - mesh_vel1[1]) + N[2]*(adv_vel2[1] - mesh_vel2[1]);

	//ms_adv_vel[0] = 0.0;
	//ms_adv_vel[1] = 0.0;

	//calculate convective term
	int nodes_number = 3;
	int dof = 2;
	int matsize = dof*nodes_number;

	//calculate density
	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	boost::numeric::ublas::bounded_matrix<double, 2, 6 > conv_opr = ZeroMatrix(dof, matsize);
	boost::numeric::ublas::bounded_matrix<double, 6, 2 > shape_func = ZeroMatrix(matsize, dof);

	for (int ii = 0; ii < nodes_number; ii++) {
	    int column = ii*dof;
	    conv_opr(0, column) = DN_DX(ii, 0) * ms_adv_vel[0] + DN_DX(ii, 1) * ms_adv_vel[1];
	    conv_opr(1, column + 1) = conv_opr(0, column);

	    shape_func(column, 0) = N[ii];
	    shape_func(column + 1, 1) = shape_func(column, 0);
	}


	//tau1*ro*Nacc.(1.0*a.grad V)
	boost::numeric::ublas::bounded_matrix<double, 6, 6 > temp_convterm = ZeroMatrix(matsize, matsize);
	temp_convterm = prod(trans(conv_opr), trans(shape_func));

	double fac = tauone*density;

	for (int ii = 0; ii < nodes_number; ii++) {
	    int row = ii * (dof + 1);
	    int loc_row = ii*dof;
	    for (int jj = 0; jj < nodes_number; jj++) {
		int column = jj * (dof + 1);
		int loc_column = jj*dof;

		M(row, column) += area * fac * temp_convterm(loc_row, loc_column);
		M(row + 1, column + 1) += area * fac * temp_convterm(loc_row + 1, loc_column + 1);
	    }
	}


	KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateGradStblAllTerms(MatrixType& K, VectorType& F, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX,const array_1d<double,3>& N, const double time, const double tauone, const double area) {
	KRATOS_TRY
		int nodes_number = 3;
	int dof = 2;

	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	//build 1*(grad q . grad p) stabilization term & assemble
	boost::numeric::ublas::bounded_matrix<double, 3, 3 > gard_opr = ZeroMatrix(nodes_number, nodes_number);
	gard_opr = 1.0 * tauone * prod(DN_DX, trans(DN_DX));

	for (int ii = 0; ii < nodes_number; ii++) {
	    int row = ii * (dof + 1) + dof;

	    for (int jj = 0; jj < nodes_number; jj++) {
		int column = jj * (dof + 1) + dof;

		K(row, column) += area * gard_opr(ii, jj);

	    }
	}

	//build 1*(grad q) (Fbody ) stabilization term & assemble
	array_1d<double, 2 > bdf = ZeroVector(2);
	const array_1d<double, 3 > bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double, 3 > bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
	const array_1d<double, 3 > bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);


	bdf[0] = N[0]*(density * bdf0[0]) + N[1]*(density * bdf1[0]) + N[2]*(density * bdf2[0]);
	bdf[1] = N[0]*(density * bdf0[1]) + N[1]*(density * bdf1[1]) + N[2]*(density * bdf2[1]);

	array_1d<double, 3 > fbd_stblterm = ZeroVector(nodes_number);
	fbd_stblterm = tauone * prod(DN_DX, bdf);


	for (int ii = 0; ii < nodes_number; ++ii) {
	    int index = ii * (dof + 1) + dof;
	    F[index] += area * fbd_stblterm[ii];
	}

	KRATOS_CATCH("")

    }
    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateGradMassStblTerms(MatrixType& M, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double,3>& N, const double tauone, const double area) {
	KRATOS_TRY
		int nodes_number = 3;
	int dof = 2;

	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	//build 1*tau1*ro Nacc grad q)
	double fac = tauone*density;
	for (int ii = 0; ii < nodes_number; ii++) {
	    int row = ii * (dof + 1);
	    for (int jj = 0; jj < nodes_number; jj++) {
		int column = jj * (dof + 1) + dof;

		//K(row,column) += -1*area * fac* N(ii) * DN_DX(jj,0);
		M(column, row) += area * fac * N[ii] * DN_DX(jj, 0);

		//K(row + 1,column) += -1*area * fac* N(ii) * DN_DX(jj,1);
		M(column, row + 1) += area * fac * N[ii] * DN_DX(jj, 1);
	    }
	}

	KRATOS_CATCH("")

    }
    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::AddBodyForceAndMomentum(VectorType& F, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double, 3 > & N, const double time, const double area, const double tauone, const double tautwo) {
	KRATOS_TRY
	int nodes_number = 3;
	int dof = 2;


	//double lump_mass_fac = area * 0.333333333333333333333333;

	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	//body  & momentum term force
	for (int ii = 0; ii < nodes_number; ii++) {
	    int index = ii * (dof + 1);
// 	    int loc_index = ii * dof;
	    const array_1d<double, 2 > bdf = GetGeometry()[ii].FastGetSolutionStepValue(BODY_FORCE);


	    F[index] += area * N[ii] * density * bdf[0];
	    F[index + 1] += area * N[ii] * density * bdf[1];



	}
#ifdef COMPRESSIBLE_MODEL_ML__PN1_PN__
	////compressibility term & assemble Lumped Mass Matrix
	double c2_inv = 1e-5;
	double compressibility_coef = 0.33333333333333333333333 * area * c2_inv /time;
// KRATOS_WATCH("time")
// KRATOS_WATCH(time)
// KRATOS_WATCH("compressibility_coef RHS")
// KRATOS_WATCH(compressibility_coef)
// KRATOS_WATCH("F before")
// KRATOS_WATCH(F)
//      // int( q * rho * 1/(rho c^2)(p)
//      //already divided by density 
	for (int ii = 0; ii < nodes_number; ++ii) {
	    int index = ii * (dof + 1) + dof;
	    F[index] += compressibility_coef * N(ii) * GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE,1);
	}
#endif
#ifdef COMPRESSIBLE_MODEL_MLPN1_MCPN
	//compressibility term & assemble Ml * 1/(k Dt) Pn+1 - Mc * 1/(k Dt) Pn
	double pressure = 0.0;
	const double pressure0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
	const double pressure1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
	const double pressure2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
	
	double c2_inv = 1e-5;
	double compressibility_coef = 0.08333333333333333333333 * area * c2_inv /time;
// KRATOS_WATCH("time")
// KRATOS_WATCH(time)
// KRATOS_WATCH("compressibility_coef RHS")
// KRATOS_WATCH(compressibility_coef)
// KRATOS_WATCH("F before")
// KRATOS_WATCH(F)
//      // int( q * rho * 1/(rho c^2)(p)
//      //already divided by density 
	array_1d<double, 3 > aux_F = ZeroVector(3);

	aux_F[0] = 2.0 * pressure0 + 1.0 * pressure1 + 1.0 * pressure2;
	aux_F[1] = 1.0 * pressure0 + 2.0 * pressure1 + 1.0 * pressure2;
	aux_F[2] = 1.0 * pressure0 + 1.0 * pressure1 + 2.0 * pressure2;
	for (int ii = 0; ii < nodes_number; ++ii) {
	    int index = ii * (dof + 1) + dof;
	    F[index] += aux_F[ii] * compressibility_coef;
	}
#endif

// 	////compressibility term & assemble Consistent Mass Matrix (Attention: Rank zero matrix, it may couse problems!)
// 	double pressure = 0.0;
// 	const double pressure0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
// 	const double pressure1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
// 	const double pressure2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
// 	
// 	double c2_inv = 1e-6;
// 	double compressibility_coef = area * c2_inv /time;
// // KRATOS_WATCH("time")
// // KRATOS_WATCH(time)
// // KRATOS_WATCH("compressibility_coef RHS")
// // KRATOS_WATCH(compressibility_coef)
// // KRATOS_WATCH("F before")
// // KRATOS_WATCH(F)
// //      // int( q * rho * 1/(rho c^2)(p)
// //      //already divided by density 
// 	pressure = N[0]* pressure0 + N[1]* pressure1 + N[2]* pressure2;
// 	for (int ii = 0; ii < nodes_number; ++ii) {
// 	    int index = ii * (dof + 1) + dof;
// 	    F[index] += compressibility_coef * N(ii) * pressure;
// 	}

	KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateResidual(const MatrixType& K, VectorType& F, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double m_coef, const double area) {
	KRATOS_TRY

	int nodes_number = 3;
	int dof = 2;

	array_1d<double, 9 > UP = ZeroVector(9);
	for (int ii = 0; ii < nodes_number; ++ii) {
	    int index = ii * (dof + 1);
	    UP[index] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY)[0];
	    UP[index + 1] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY)[1];
	    UP[index + 2] = GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE);
	}

	/*Rest Viscous Forces from dampmatrix*/
	double mu;
	double density;
	calculatedensity(GetGeometry(), density, mu);

	boost::numeric::ublas::bounded_matrix<double, 3, 6 > B = ZeroMatrix(3, 6);
//         boost::numeric::ublas::bounded_matrix<double, 3, 3 > C = ZeroMatrix(3, 3);
	array_1d<double, 6 > auxDevStressVector = ZeroVector(6);
	array_1d<double, 3 > grad_sym_vel = ZeroVector(3);
	double app_mu;
	double app_mu_derivative;
	double gamma_dot;
	array_1d<double, 3 > aux_1 = ZeroVector(3);

	//calculating operator B
	CalculateB(B, DN_DX);
	// KRATOS_WATCH(B)
	/*Calculating residual vector */
	F -= prod(K, UP);

	/*Add Bt*sigma_dev (Viscous Forces)*/
	//sigma dev intern

// 	Bingham Fluid:
	CalculateApparentViscosity(app_mu, app_mu_derivative, grad_sym_vel, gamma_dot, B, mu, m_coef);
// 	Newtonian Fluid: Leave Decommented the CalculateApparentviscosity (we need grad_sym_vel) and decomment the following line
//	Remember to modify CalculateViscousTerm and CalculateTau
// 	app_mu = mu;
	
	aux_1 = 2.0 * app_mu * grad_sym_vel;
	aux_1[2] *= 0.5; //considering Voigt notation for the gradient of velocity (alternative to the C matrix of the viscous term.

	mDevStress = 0.5 * aux_1[0] * aux_1[0] + 0.5 * aux_1[1] * aux_1[1] + aux_1[2] * aux_1[2];
// 	mDevStress =  aux_1[0] * aux_1[0] +  aux_1[1] * aux_1[1] + 2 * aux_1[2] * aux_1[2];

	mDevStress = sqrt(mDevStress);
	
	
	auxDevStressVector = prod(trans(B), aux_1);
	/*TO DECOMMENT*/

	for (int ii = 0; ii < nodes_number; ii++) {
	    int row = ii * (dof + 1);
	    int loc_row = ii * dof;

	    F[row] -= auxDevStressVector[loc_row] * area;
	    F[row + 1] -= auxDevStressVector[loc_row + 1] * area;

	}


	KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************

//     void NoNewtonianASGS2D::CalculateResidual(const MatrixType& K, VectorType& F, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double m_coef, const double area) {
// 	KRATOS_TRY
// 
// 	int nodes_number = 3;
// 	int dof = 2;
// 
// 	array_1d<double, 9 > UP = ZeroVector(9);
// 	for (int ii = 0; ii < nodes_number; ++ii) {
// 	    int index = ii * (dof + 1);
// 	    UP[index] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY)[0];
// 	    UP[index + 1] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY)[1];
// 	    UP[index + 2] = GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE);
// 	}
// 
// 	/*Rest Viscous Forces from dampmatrix*/
// 	double mu;
// 	double density;
// 	calculatedensity(GetGeometry(), density, mu);
// 
// 	boost::numeric::ublas::bounded_matrix<double, 3, 6 > B = ZeroMatrix(3, 6);
// //         boost::numeric::ublas::bounded_matrix<double, 3, 3 > C = ZeroMatrix(3, 3);
// 	array_1d<double, 6 > auxDevStressVector = ZeroVector(6);
// 	array_1d<double, 3 > grad_sym_vel = ZeroVector(3);
// 	array_1d<double, 3 > app_mu = ZeroVector(3);
// 	double app_mu_derivative;
// 	double gamma_dot;
// 	array_1d<double, 3 > aux_1 = ZeroVector(3);
// 
// 	//calculating operator B
// 	CalculateB(B, DN_DX);
// 	// KRATOS_WATCH(B)
// 	/*Calculating residual vector */
// 	F -= prod(K, UP);
// 
// 	/*Add Bt*sigma_dev (Viscous Forces)*/
// 	//sigma dev intern
// 
// // 	Bingham Fluid:
// 	CalculateApparentViscosity(app_mu, app_mu_derivative, grad_sym_vel, gamma_dot, B, mu, m_coef);
// // 	Newtonian Fluid: Leave Decommented the CalculateApparentviscosity (we need grad_sym_vel) and decomment the following line
// //	Remember to modify CalculateViscousTerm and CalculateTau
// // 	app_mu = mu;
// 	double average_app_mu = 0.0;
// 	for(unsigned int ii=0.0, ii> nodes_number, ii++){
// 	    average_app_mu += app_mu[ii];
// 	}
// 	average_app_mu *= 0.3333333333333333333333333333;
// 	aux_1 = 2.0 * average_app_mu * grad_sym_vel;	
// 	aux_1[2] *= 0.5; //considering Voigt notation for the gradient of velocity (alternative to the C matrix of the viscous term.
// 
// 	mDevStress = 0.5 * aux_1[0] * aux_1[0] + 0.5 * aux_1[1] * aux_1[1] + aux_1[2] * aux_1[2];
// // 	mDevStress =  aux_1[0] * aux_1[0] +  aux_1[1] * aux_1[1] + 2 * aux_1[2] * aux_1[2];
// 
// 	mDevStress = sqrt(mDevStress);
// 	
// 	
// 	auxDevStressVector = prod(trans(B), aux_1);
// 	/*TO DECOMMENT*/
// 	double fac = area * 0.3333333333333333333333333333;
// 	for (int ii = 0; ii < nodes_number; ii++) {
// 	    int row = ii * (dof + 1);
// 	    int loc_row = ii * dof;
// 
// 	    F[row] -= auxDevStressVector[loc_row] * fac;
// 	    F[row + 1] -= auxDevStressVector[loc_row + 1] * fac;
// 
// 	}
// 
// 
// 	KRATOS_CATCH("")
//     }
    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::ComputeProjections(array_1d<double, 6 > & adv_proj, array_1d<double, 3 > & div_proj, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const double tauone, const double tautwo, const array_1d<double, 3 > & N, const double area, const double time) {
	unsigned int number_of_nodes = GetGeometry().PointsNumber();
	unsigned int dim = 2;

	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	const array_1d<double, 3 > & adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double, 3 > & adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double, 3 > & adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);


	array_1d<double,2> ms_adv_vel;
	ms_adv_vel[0] = N[0]*(adv_vel0[0] - mesh_vel0[0]) + N[1]*(adv_vel1[0] - mesh_vel1[0]) + N[2]*(adv_vel2[0] - mesh_vel2[0]);
	ms_adv_vel[1] = N[0]*(adv_vel0[1] - mesh_vel0[1]) + N[1]*(adv_vel1[1] - mesh_vel1[1]) + N[2]*(adv_vel2[1] - mesh_vel2[1]);


	double const_adv_proj_X = 0.0;
	double const_adv_proj_Y = 0.0;
	double mean_div_proj = 0.0;
	array_1d<double, 3 > mean_new_vel = ZeroVector(3);
	array_1d<double, 3 > mean_old_vel = ZeroVector(3);
	array_1d<double, 3 > mean_bdf = ZeroVector(3);

	for (unsigned int i = 0; i < number_of_nodes; i++) {

	    //int index = i*dim;
	    double pr = GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
	    const array_1d<double, 3 > & vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
	    const array_1d<double, 3 > & old_vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1);
	    const array_1d<double, 3 > & bdf = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);

	    //const array_1d<double,2>& bdf = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
	    // to consider the jump gradp/ro is calculated
	    //pr = pr/density;

	    //adv_proj = PI(ro*dv/dt + ro*a.gradU + gradP - f) considering lumped mass PI() = ()
	    //calculate constant part of RES ->ro*a.gradU + gradP
	    const_adv_proj_X += (pr * DN_DX(i, 0) + density * (ms_adv_vel[0] * DN_DX(i, 0) + ms_adv_vel[1] * DN_DX(i, 1)) * vel[0]);
	    const_adv_proj_Y += (pr * DN_DX(i, 1) + density * (ms_adv_vel[0] * DN_DX(i, 0) + ms_adv_vel[1] * DN_DX(i, 1)) * vel[1]);

	    //div_proj = PI(ro*divU)
	    mean_div_proj += density * (DN_DX(i, 0) * vel[0] + DN_DX(i, 1) * vel[1]);

	    //calcuale mean velocity and body force
	    mean_new_vel += 0.3333333333333333333333333333 * vel;
	    mean_old_vel += 0.3333333333333333333333333333 * old_vel;
	    mean_bdf += 0.3333333333333333333333333333 * bdf;
	}


	for (unsigned int i = 0; i < number_of_nodes; i++) {
	    int index = i*dim;

	    adv_proj[index] = area * N[i]*(density * (mean_new_vel[0] - mean_old_vel[0]) / time + const_adv_proj_X - density * mean_bdf[0]);
	    adv_proj[index + 1] = area * N[i]*(density * (mean_new_vel[1] - mean_old_vel[1]) / time + const_adv_proj_Y - density * mean_bdf[1]);

	    div_proj[i] = area * N[i] * density*mean_div_proj;

	    //update projections
	    array_1d<double, 3 > & advtermproj = GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
	    advtermproj[0] += adv_proj[index];
	    advtermproj[1] += adv_proj[index + 1];

	    double& divtermproj = GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ);
	    divtermproj += div_proj[i];



	    //calculate nodal area

	    GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += 0.333333333333333333 * area;

	}

  
    }

    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::AddProjectionForces(VectorType& F, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const double area, const double tauone, const double tautwo) {
	unsigned int number_of_nodes = GetGeometry().PointsNumber();
	unsigned int dim = 2;
	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	const array_1d<double, 3 >& advproj_0 = GetGeometry()[0].FastGetSolutionStepValue(ADVPROJ);
	const array_1d<double, 3 >& advproj_1 = GetGeometry()[1].FastGetSolutionStepValue(ADVPROJ);
	const array_1d<double, 3 >& advproj_2 = GetGeometry()[2].FastGetSolutionStepValue(ADVPROJ);

	const double div_proj_0 = GetGeometry()[0].FastGetSolutionStepValue(DIVPROJ);
	const double div_proj_1 = GetGeometry()[1].FastGetSolutionStepValue(DIVPROJ);
	const double div_proj_2 = GetGeometry()[2].FastGetSolutionStepValue(DIVPROJ);



	//mean values
	double mean_x_adv = 0.3333333333333333 * (advproj_0[0] + advproj_1[0] + advproj_2[0]);
	double mean_y_adv = 0.3333333333333333 * (advproj_0[1] + advproj_1[1] + advproj_2[1]);

	double mean_div = 0.3333333333333333 * (div_proj_0 + div_proj_1 + div_proj_2);

	for (unsigned int ii = 0; ii < number_of_nodes; ii++) {
	    int index = ii * (dim + 1);
	    const array_1d<double, 3 > & vel = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY);

	    const array_1d<double, 3 > & mesh_vel = GetGeometry()[ii].FastGetSolutionStepValue(MESH_VELOCITY);
	    array_1d<double, 2 > adv_vel = ZeroVector(2);
	    adv_vel[0] = vel[0] - mesh_vel[0];
	    adv_vel[1] = vel[1] - mesh_vel[1];

	    //tauone*ro*(xi,a.gradv)
	    double proj;

	    proj = mean_x_adv * (adv_vel[0] * DN_DX(ii, 0) + adv_vel[1] * DN_DX(ii, 1));
	    F[index] += (tauone * 1.0 * area * proj);

	    proj = mean_y_adv * (adv_vel[0] * DN_DX(ii, 0) + adv_vel[1] * DN_DX(ii, 1));
	    F[index + 1 ] += (tauone * 1.0 * area * proj);


	    //tauone*(xi,gradq)
	    proj = (mean_x_adv * DN_DX(ii, 0) + mean_y_adv * DN_DX(ii, 1));
	    F[index + 2] += (tauone * area * proj);

	    //tautwo*(divv)

	    F[index] += (tautwo * area * mean_div * DN_DX(ii, 0));
	    F[index + 1] += (tautwo * area * mean_div * DN_DX(ii, 1));

	}




    }

    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) {
	KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
	unsigned int dim = 2;
	unsigned int node_size = dim + 1;


	if (rResult.size() != number_of_nodes * node_size)
	    rResult.resize(number_of_nodes * node_size, false);

	for (unsigned int i = 0; i < number_of_nodes; i++) {
	    rResult[i * node_size] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
	    rResult[i * node_size + 1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
	    rResult[i * node_size + 2] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
	}
	KRATOS_CATCH("")

    }

    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {
	KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
	unsigned int dim = 2;
	unsigned int node_size = dim + 1;


	if (ElementalDofList.size() != number_of_nodes * node_size)
	    ElementalDofList.resize(number_of_nodes * node_size);

	for (unsigned int i = 0; i < number_of_nodes; i++) {
	    ElementalDofList[i * node_size] = GetGeometry()[i].pGetDof(VELOCITY_X);
	    ElementalDofList[i * node_size + 1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
	    ElementalDofList[i * node_size + 2] = GetGeometry()[i].pGetDof(PRESSURE);
	}
	KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateB(
	    boost::numeric::ublas::bounded_matrix<double, 3, 6 > & B,
	    const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX) {
	KRATOS_TRY
		unsigned int dim = 2;
	unsigned int node_size = dim + 1;

	for (unsigned int i = 0; i < node_size; i++) {
	    unsigned int index = dim*i;

	    B(0, index) = DN_DX(i, 0);
	    B(0, index + 1) = 0.0;
	    B(1, index) = 0.0;
	    B(1, index + 1) = DN_DX(i, 1);
	    B(2, index) = DN_DX(i, 1);
	    B(2, index + 1) = DN_DX(i, 0);


	    //CalculateBi(Bi,F,DN_DX,i);
	    //MathUtils<double>::WriteMatrix(B,Bi,0,index);
	}
	KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateGradSymVel(array_1d<double, 3 > & grad_sym_vel, double & gamma_dot,
	    const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & B) {
	KRATOS_TRY
	unsigned int dim = 2;
	unsigned int nodes_number = dim + 1;

	array_1d<double, 6 > U;
	for (unsigned int ii = 0; ii < nodes_number; ++ii) {
	    int index = ii * (dim);
	    const array_1d<double,3>& vel =  GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY);
	    U[index] = vel[0];
	    U[index + 1] = vel[1];
	}

	grad_sym_vel = prod(B, U);

// Norm of the gradient of velocity:
//         gamma_dot = grad_sym_vel[0] * grad_sym_vel[0] + grad_sym_vel[1] * grad_sym_vel[1] + 0.5 * grad_sym_vel[2] * grad_sym_vel[2];
// Gamma dot found in literature!!!:
	gamma_dot = 2.0 * grad_sym_vel[0] * grad_sym_vel[0] + 2.0 * grad_sym_vel[1] * grad_sym_vel[1] +  grad_sym_vel[2] * grad_sym_vel[2];
	gamma_dot = sqrt(gamma_dot);

// 	if(gamma_dot < 1e-5){
// 	  gamma_dot=1e-5;
// 	}
	
	KRATOS_CATCH("")
    }



    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateApparentViscosity(double & app_mu, double & app_mu_derivative,
	    array_1d<double, 3 >&  grad_sym_vel, double & gamma_dot,
	    const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & B,
	    const double & mu, const double & m_coef) {
	KRATOS_TRY
	app_mu = 0.0;

// 	double m_coef = 3000;
// KRATOS_WATCH(m_coef)	
	double aux_1;
	CalculateGradSymVel(grad_sym_vel, gamma_dot, B);
	
	
	  // The yield is variable: it decreases where water is present
	  unsigned int nodes_number = 3;
	  double yield = 0.0;
	  double friction_angle_tangent = 0.0; 
	  double water_pressure = 0.0;
	  double gamma_dot_inv;
// 	  double solid_pressure = 0.0;
// 	  double seepage_drag_x = 0.0;
	  
	  

	  
	for (unsigned int ii = 0; ii < nodes_number; ++ii) {
	      if(GetGeometry()[ii].FastGetSolutionStepValue(WATER_PRESSURE) >= 0.0){
		    water_pressure +=  GetGeometry()[ii].FastGetSolutionStepValue(WATER_PRESSURE);		    
	      }
// 	      yield +=  GetGeometry()[ii].FastGetSolutionStepValue(YIELD_STRESS);
	      friction_angle_tangent += GetGeometry()[ii].FastGetSolutionStepValue(INTERNAL_FRICTION_ANGLE);
	      if(GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE) >= 0.0){
		    yield +=  GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE);
// 		    solid_pressure +=  GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE);
	      }
       	      yield +=  GetGeometry()[ii].FastGetSolutionStepValue(YIELD_STRESS);//this is the COHESION of Mohr Coulomb failure criteria!!!
// 	      seepage_drag_x += GetGeometry()[ii].FastGetSolutionStepValue(SEEPAGE_DRAG_X);
	      //CHECK NEGATIVE VALUES....

	  }
	  friction_angle_tangent /= nodes_number;
	  water_pressure /= nodes_number;
	  yield /= nodes_number;
// KRATOS_WATCH(friction_angle_tangent)
// KRATOS_WATCH(yield)
// 	  solid_pressure /= nodes_number;
// 	  seepage_drag_x /= nodes_number;
// 	  pay attention: negative yield stress meaningful
// 	  if(yield < solid_pressure) //otherwise we are enshuring a minimum resistance of the material
// 	    yield = solid_pressure; 
// 	  else
// 	    yield -= 10.0 * seepage_drag_x; //to make more week the downstream toe, where the seepage drag in x in higher.
// 	  
	  if(water_pressure < yield){
	      yield -= water_pressure;
	      yield *= friction_angle_tangent;
	  }
	  else
	      yield = 0.0;
	  
// 	  if(yield < 100.0)
//A sort of COHESION!!!!! not present at a micro level but possibly present at a macro level
// 	    yield += 100.0;
	  
	  
// KRATOS_WATCH(yield)

#ifdef EXPONENCIAL_MODEL
////////EXPONENCIAL MODEL
	if (gamma_dot > 1e-3) {
	    aux_1 = 1.0 - exp(-(m_coef * gamma_dot));
	    app_mu = mu + (yield / gamma_dot) * aux_1;
// 			gamma_dot_inv = 1.0/gamma_dot;
	    if (app_mu < mu) {
		KRATOS_ERROR(std::logic_error, "!!!!!!!!!!!  APPARENT VISCOSITY < VISCOSITY !!!!!!!!", this->Id());
	    }
	} else {
	    app_mu = mu + yield * m_coef ;
// // 			gamma_dot_inv = 0.0;
	}
#else	
// ////////BILINEAR MODEL
      double mu_s = 1e7;
      double gamma_dot_lim = yield/mu_s;
      
      if(gamma_dot >= gamma_dot_lim){
	if (gamma_dot_lim <= 1e-16){
	  app_mu = mu ; //newtonian fluid, no rigid behavior. yield ~ 0;
	}
	else{
	  app_mu = (mu*(gamma_dot-gamma_dot_lim) + yield)/gamma_dot ;
	}
      }
      else{
	app_mu = mu_s ;
      }
#endif

/////// UPPERBOUND IN APPARENT VISCOSITY /////////
//       if(app_mu > 1e4){
// 	app_mu = 1e4;
//       }


// // // // 	if (gamma_dot <= 1e-10) gamma_dot_inv=1e10;
// // // // 	else  
// 	gamma_dot_inv= 1.0 / gamma_dot;
// // // //         app_mu_derivative = yield * gamma_dot_inv*(- gamma_dot_inv + exp(-(m_coef * gamma_dot))*(gamma_dot_inv + m_coef));
//  	app_mu_derivative = - yield * gamma_dot_inv * gamma_dot_inv * (1 - exp(-(m_coef * gamma_dot))*(1 - m_coef * gamma_dot));
// // // //	app_mu_derivative = - yield * gamma_dot_inv * gamma_dot_inv;

	/*

//Calculating nodal yield and nodal app_mu before the calculation of elemental apparent viscosity
	  // The yield is variable: it decreases where water is present
	  unsigned int nodes_number = 3;
	  array_1d<double, 3 > yield = ZeroVector(3);
	  double friction_angle_tangent = 0.0;
	  array_1d<double, 3 > water_pressure = ZeroVector(3);
	  array_1d<double, 3 > nodal_app_mu = ZeroVector(3);

	for (unsigned int ii = 0; ii < nodes_number; ++ii) {
	      if(GetGeometry()[ii].FastGetSolutionStepValue(WATER_PRESSURE) >= 0.0){
		    water_pressure[ii] =  GetGeometry()[ii].FastGetSolutionStepValue(WATER_PRESSURE);    
	      }
	      friction_angle_tangent += GetGeometry()[ii].FastGetSolutionStepValue(INTERNAL_FRICTION_ANGLE);
	      if(GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE) >= 0.0){
		    yield[ii] =  GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE);
	      }
	  }
  
	  friction_angle_tangent /= nodes_number;
// 	  water_pressure /= nodes_number;
// 	  yield /= nodes_number;

	for (unsigned int ii = 0; ii < nodes_number; ++ii) {
	      if(water_pressure[ii] < yield[ii]){
		      yield[ii] -= water_pressure[ii];
	      }
	      else{
		yield[ii] = 0.0;
	      }
	}
	yield *= friction_angle_tangent;



//////////EXPONENCIAL FORM
	if (gamma_dot > 1e-10) {
		aux_1 = 1.0 - exp(-(m_coef * gamma_dot));
		for (unsigned int ii = 0; ii < nodes_number; ++ii) {
		    nodal_app_mu[ii] = mu + (yield[ii] / gamma_dot) * aux_1;
		}
	}
	else {
	      for (unsigned int ii = 0; ii < nodes_number; ++ii) {
		nodal_app_mu[ii] = mu + yield[ii]*m_coef;
	      }
	}
	
	for (unsigned int ii = 0; ii < nodes_number; ++ii) {
	    app_mu += nodal_app_mu[ii];
	}
	app_mu /= nodes_number;
*/	
	
// 	for (unsigned int ii = 0; ii < nodes_number; ++ii) {
// 	  if(GetGeometry()[ii].Y()> (-GetGeometry()[ii].X() + 3.0))
// 	    app_mu = 150.0;
// 	}
// 	if (gamma_dot <= 1e-10) gamma_dot_inv=1e10;
// 	else  gamma_dot_inv= 1.0/gamma_dot;
// 	app_mu_derivative = yield * gamma_dot_inv*(- gamma_dot_inv + exp(-(m_coef * gamma_dot))*(gamma_dot_inv + m_coef));

	KRATOS_CATCH("")
    }

   //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::CalculateNodalApparentViscosity(
	    array_1d<double, 3 >&  grad_sym_vel, double & gamma_dot,
	    const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & B,
	    const double & mu, const double & m_coef) {
	KRATOS_TRY
// 	app_mu = ZeroVector(3);

	double aux_1;
	CalculateGradSymVel(grad_sym_vel, gamma_dot, B);
	
	
	  // The yield is variable: it decreases where water is present
	  unsigned int nodes_number = 3;
	  array_1d<double, 3 > yield = ZeroVector(3); //nodal yield
	  array_1d<double, 3 > water_pressure = ZeroVector(3); //nodal_water_pressure
	  double friction_angle_tangent = 0.0; 
	  
	for (unsigned int ii = 0; ii < nodes_number; ++ii) {
	      if(GetGeometry()[ii].FastGetSolutionStepValue(WATER_PRESSURE) >= 0.0){
		    water_pressure[ii] +=  GetGeometry()[ii].FastGetSolutionStepValue(WATER_PRESSURE);		    
	      }
	      friction_angle_tangent += GetGeometry()[ii].FastGetSolutionStepValue(INTERNAL_FRICTION_ANGLE);
	      if(GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE) >= 0.0){
		    yield[ii] +=  GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE);
	      }
	  }
	  friction_angle_tangent /= nodes_number;

	for (unsigned int ii = 0; ii < nodes_number; ++ii) {
	  if(water_pressure[ii] < yield[ii]){
	      yield[ii] -= water_pressure[ii];
	      yield[ii] *= friction_angle_tangent;
	  }
	  else
	      yield[ii] = 0.0;
	}
#ifdef EXPONENCIAL_MODEL
////////EXPONENCIAL MODEL
      for(unsigned int ii = 0; ii < nodes_number; ++ii) {
	double & app_mu = GetGeometry()[ii].FastGetSolutionStepValue(EFFECTIVE_VISCOSITY);
	if (gamma_dot > 1e-3) {
	    aux_1 = 1.0 - exp(-(m_coef * gamma_dot));
	    app_mu = mu + (yield[ii] / gamma_dot) * aux_1;
// 	    if (app_mu < mu) {
// 		KRATOS_ERROR(std::logic_error, "!!!!!!!!!!!  APPARENT VISCOSITY < VISCOSITY !!!!!!!!", this->Id());
// 	    }
	} else {
	    app_mu = mu + yield[ii] * m_coef ;
// // 			gamma_dot_inv = 0.0;
	}
	 /////// UPPERBOUND IN APPARENT VISCOSITY /////////
// 	if(app_mu > 1e4){
// 	  app_mu = 1e4;
// 	}
      }
#else	
// ////////BILINEAR MODEL
      double mu_s = 1e7;
      double average_yield = 0.0;
      for(unsigned int ii=0, ii<nodes_number; i++){
	average_yield += yield[ii];
      }
      averaged_yield /= nodes_number;
      double gamma_dot_lim = averaged_yield/mu_s;
      
      
      for(unsigned int ii = 0; ii < nodes_number; ++ii) { 
	double & app_mu = GetGeometry()[ii].FastGetSolutionStepValue(EFFECTIVE_VISCOSITY);
	if(gamma_dot >= gamma_dot_lim){
	  if (gamma_dot_lim <= 1e-10){
	    app_mu = mu ; //newtonian fluid, no rigid behavior. yield ~ 0;
	  }
	  else{
	    app_mu = (mu*(gamma_dot-gamma_dot_lim) + yield[ii])/gamma_dot ;
	  }
	}
	else{
	  app_mu = mu_s ;
	}
// 	  /////// UPPERBOUND IN APPARENT VISCOSITY /////////
// 	if(app_mu > 1e5){
// 	  app_mu = 1e5;
// 	}
      }
#endif


	KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
	    array_1d<double, 3 > & Output,
	    const ProcessInfo& rCurrentProcessInfo) {

	array_1d<double, 6 > adv_proj = ZeroVector(6);
	array_1d<double, 3 > div_proj = ZeroVector(3);

	double delta_t = rCurrentProcessInfo[DELTA_TIME];
	
	boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
	array_1d<double, 3 > N;
//         array_1d<double, 2 > ms_adv_vel;

	//getting data for the given geometry
	double Area;
	GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

	double tauone;
	double tautwo;
	CalculateTau(DN_DX,N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);

	ComputeProjections(adv_proj, div_proj, DN_DX, tauone, tautwo, N, Area, delta_t);


    }

    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) {

	double delta_t = rCurrentProcessInfo[DELTA_TIME];
	const double m_coef = rCurrentProcessInfo[M];
	
	boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
	boost::numeric::ublas::bounded_matrix<double, 3, 6 > B = ZeroMatrix(3, 6);
	array_1d<double, 3 > grad_sym_vel = ZeroVector(3);
	array_1d<double, 3 > N;
//         array_1d<double, 2 > ms_adv_vel;

	//getting data for the given geometry
	double Area;
	GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);
	double tauone;
	double tautwo;
	
	double app_mu;
	double mu;
	double density;
      //  double average_app_mu;
	double app_mu_derivative;
	double gamma_dot;

	calculatedensity(GetGeometry(), density, mu);
	CalculateB(B, DN_DX);    
	CalculateApparentViscosity(app_mu, app_mu_derivative, grad_sym_vel, gamma_dot, B, mu, m_coef);
// 	  average_app_mu = (app_mu[0]+app_mu[1]+app_mu[2])*0.333333333333333333333333;
	CalculateTau(DN_DX,N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);
	
	if (rVariable == THAWONE) {
	    for (unsigned int PointNumber = 0;
		    PointNumber < 1;
		    PointNumber++) {
		rValues[PointNumber] = tauone;
	    }
	}
	if (rVariable == THAWTWO) {
	    for (unsigned int PointNumber = 0;
		    PointNumber < 1; PointNumber++) {

		rValues[PointNumber] = tautwo;
	    }
	}
	if (rVariable == IS_WATER_ELEMENT) {
	    for (unsigned int PointNumber = 0;
		    PointNumber < 1; PointNumber++) {
		//	KRATOS_WATCH(this->GetValue(IS_WATER));
		//	KRATOS_WATCH(this->Info());
		rValues[PointNumber] = this->GetValue(IS_WATER_ELEMENT);
	    }
	}
//PROVISIONALbegin---only for debugging
	if (rVariable == EQ_STRAIN_RATE) {//gamma dot
// 	  boost::numeric::ublas::bounded_matrix<double, 3, 6 > B = ZeroMatrix(3, 6);
// 	  array_1d<double, 3 > grad_sym_vel = ZeroVector(3);
// 	  double gamma_dot;
// // 	  double Area;
// // 	  GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);
// 
// 	  CalculateB(B, DN_DX);      
// 
// 	  CalculateGradSymVel(grad_sym_vel, gamma_dot, B);
	    
	    for (unsigned int PointNumber = 0;
		    PointNumber < 1; PointNumber++) {
		rValues[PointNumber] = gamma_dot;

	    }
	}
	if (rVariable == MU) {//app mu
/*	  boost::numeric::ublas::bounded_matrix<double, 3, 6 > B = ZeroMatrix(3, 6);
	  array_1d<double, 3 > grad_sym_vel = ZeroVector(3);
  	  double app_mu = 0.0;
// 	  double gamma_dot = 0.0;
	  double mu;
	  double density;
	//  double average_app_mu;
	  double app_mu_derivative;
	  double gamma_dot;
// 	  double Area;
// 	  GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);
	  calculatedensity(GetGeometry(), density, mu);
	  CalculateB(B, DN_DX);    
	  CalculateApparentViscosity(app_mu, app_mu_derivative, grad_sym_vel, gamma_dot, B, mu, m_coef);
// 	  average_app_mu = (app_mu[0]+app_mu[1]+app_mu[2])*0.333333333333333333333333;*/
	    for (unsigned int PointNumber = 0;
		    PointNumber < 1; PointNumber++) {
		rValues[PointNumber] = app_mu;

	    }
	}
	if (rVariable == TAU) {//SQRT(0.5 DEVSTRESS:DEVSTRESS)
	    
	    for (unsigned int PointNumber = 0;
		    PointNumber < 1; PointNumber++) {
		rValues[PointNumber] = mDevStress;

	    }
	}	
	
//PROVISIONALend---only for debugging
    }

      
    
    //*************************************************************************************
    //*************************************************************************************

    void NoNewtonianASGS2D::calculatedensity(Geometry< Node < 3 > > geom, double& density, double& viscosity) {

	double kk = 0.0;
	density = 0.0;
	viscosity = 0.0;
	for (int ii = 0; ii < 3; ++ii) {
	    kk++;
	    density += geom[ii].FastGetSolutionStepValue(DENSITY);
	    viscosity += geom[ii].FastGetSolutionStepValue(VISCOSITY);
	}

	density /= kk;
	viscosity /= kk;
	//Here we calculate Dynamic viscosity from Kinemeatic viscosity
	viscosity *= density;

    }
    //*************************************************************************************
    //*************************************************************************************

    void NoNewtonianASGS2D::CalculateTau(const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX, const array_1d<double,3>& N, double& tauone, double& tautwo, const double time, const double area, const ProcessInfo& rCurrentProcessInfo) {
	KRATOS_TRY
		//calculate mean advective velocity and taus
		const array_1d<double, 3 > & adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double, 3 > & adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
	const array_1d<double, 3 > & adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, 0);
	const array_1d<double, 3 > & mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);


	array_1d<double,2> ms_adv_vel;
	ms_adv_vel[0] = N[0]*(adv_vel0[0] - mesh_vel0[0]) + N[1]*(adv_vel1[0] - mesh_vel1[0]) + N[2]*(adv_vel2[0] - mesh_vel2[0]);
	ms_adv_vel[1] = N[0]*(adv_vel0[1] - mesh_vel0[1]) + N[1]*(adv_vel1[1] - mesh_vel1[1]) + N[2]*(adv_vel2[1] - mesh_vel2[1]);

	//ms_adv_vel[0] = 0.0;
	//ms_adv_vel[1] = 0.0;


	double advvel_norm = ms_adv_vel[0] * ms_adv_vel[0] + ms_adv_vel[1] * ms_adv_vel[1];
	advvel_norm = sqrt(advvel_norm);

	double ele_length = 2.0 * sqrt(area / 3.00);

	double mu;
	double gamma_dot;
	//const double mu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
	//const double mu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
	//const double mu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
	//mu = 0.333333333333333333333333*(mu0 + mu1 + mu2);

	double density;
	calculatedensity(GetGeometry(), density, mu);

/*provisional*/
	boost::numeric::ublas::bounded_matrix<double, 3, 6 > B = ZeroMatrix(3, 6);
	array_1d<double, 3 > grad_sym_vel = ZeroVector(3); 
	double app_mu_derivative;

	CalculateB(B, msDN_DX);
//         for (int ii = 0; ii < nodes_number; ++ii) {
//             yield += GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE);
//         }

  
	//Bingham
	const double m_coef = rCurrentProcessInfo[M];
	CalculateApparentViscosity(mu, app_mu_derivative, grad_sym_vel, gamma_dot, B, mu, m_coef);
	//Newtonian: comment the CalculateApparentViscosity funcion and nothing more (remember to modify CalculateResidual and CalculateViscousTerm
	//do nothing --> we don't need the calculation of grad_sym_vel in this case!!!
	
	
/*provisional*/	
	const double dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
	tauone = 1.0 / (dyn_st_beta / time + 4.0 * mu / (ele_length * ele_length * density) + 2.0 * advvel_norm  / ele_length);

// 	tautwo = mu / density + 1.0 * ele_length * advvel_norm / 2.0;
	tautwo = 0.0;



// 	//being lagrangian no advective velocity should be present!
// 	if(advvel_norm> 1e-10)
// 	  KRATOS_ERROR(std::logic_error,"LAGRANGIAN TEST: advective velocity is not zero!!!!!!","");

	KRATOS_CATCH("")


    }


    //*************************************************************************************
    //*************************************************************************************

    void NoNewtonianASGS2D::GetFirstDerivativesVector(Vector& values, int Step) {
	const unsigned int number_of_nodes = GetGeometry().size();
	const unsigned int dim = GetGeometry().WorkingSpaceDimension();
	unsigned int MatSize = number_of_nodes * (dim + 1);
	if (values.size() != MatSize) values.resize(MatSize, false);
	for (unsigned int i = 0; i < number_of_nodes; i++) {
	    unsigned int index = i * (dim + 1);
	    values[index] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_X, Step);
	    values[index + 1] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y, Step);
	    values[index + 2] = GetGeometry()[i].GetSolutionStepValue(PRESSURE, Step);

	}
    }
    //************************************************************************************
    //************************************************************************************

    void NoNewtonianASGS2D::GetSecondDerivativesVector(Vector& values, int Step) {
	const unsigned int number_of_nodes = GetGeometry().size();
	const unsigned int dim = GetGeometry().WorkingSpaceDimension();
	unsigned int MatSize = number_of_nodes * (dim + 1);
	if (values.size() != MatSize) values.resize(MatSize, false);
	for (unsigned int i = 0; i < number_of_nodes; i++) {
	    unsigned int index = i * (dim + 1);
	    values[index] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X, Step);
	    values[index + 1] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y, Step);
	    values[index + 2] = 0.0;
	}
    }
  

    
    //************************************************************************************
    //************************************************************************************
      
    
} // Namespace Kratos


