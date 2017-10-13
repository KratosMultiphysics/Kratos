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
//   Last modified by:    $Author: pavel $
//   Date:                $Date: 2009-01-21 14:15:02 $
//   Revision:            $Revision: 1.0 $
//
//

//#define GRADPN_FORM
//#define STOKES

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/v-p-precond_2d.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos {

    //************************************************************************************
    //************************************************************************************

    VP_PRECOND2D::VP_PRECOND2D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    VP_PRECOND2D::VP_PRECOND2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {

    }

    Element::Pointer VP_PRECOND2D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {

        KRATOS_TRY
        return Element::Pointer(new VP_PRECOND2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }

    VP_PRECOND2D::~VP_PRECOND2D() {
    }

    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY

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

        boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N; 

        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

        double tauone;
        double tautwo;
        CalculateTau(N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);
	
	unsigned int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

	if (FractionalStepNumber==2)
	{

		//unsigned int matsize = nodes_number * (dim + 1);	

		//viscous term
		CalculateViscousTerm(rLeftHandSideMatrix, DN_DX, Area, rCurrentProcessInfo);
		
		CalculateAdvectiveTerm(rLeftHandSideMatrix, DN_DX,N, tauone, tautwo, delta_t, Area);

		//calculate pressure term
		CalculatePressureTerm(rLeftHandSideMatrix, DN_DX, N, delta_t, Area, rCurrentProcessInfo);

		//compute projections
		//unsigned int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];
		//add stabilization only in the correction step
		
		//stabilization terms
		//CalculateDivStblTerm(rLeftHandSideMatrix, DN_DX, tautwo, Area);
		CalculateAdvStblAllTerms(rLeftHandSideMatrix, rRightHandSideVector, DN_DX, N, tauone, delta_t, Area, rCurrentProcessInfo);
		CalculateGradStblAllTerms(rLeftHandSideMatrix, rRightHandSideVector, DN_DX,N, delta_t, tauone, Area);

		//THIS IS BECAUSE STATIC SCHEME JUST CALLS CalculateLocalSystem -> for correction step we must do the time integration here 
		

		CalculateMassContribution(rLeftHandSideMatrix, delta_t, Area, rCurrentProcessInfo);

		//double alpha_bossak=-0.2;
		double lump_mass_fac = (1.0/delta_t)*Area * 0.333333333333333333333333;
		double density;
		double mu;
		calculatedensity(GetGeometry(), density, mu);
		array_1d<double, 2 > temp_vec; 			
				
		for (int i=0;i<3;i++)
			{
			temp_vec[0] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X,1);
			temp_vec[1] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y,1);
			int index=i*3;
			rRightHandSideVector[index]+=density*lump_mass_fac*temp_vec[0];
			rRightHandSideVector[index+1]+=density*lump_mass_fac*temp_vec[1];
			}	
		
		//add stablilization terms due to advective term (a)grad(V) * ro*Acce
		CalculateAdvMassStblTerms(rLeftHandSideMatrix, rRightHandSideVector, DN_DX, N, tauone, Area, rCurrentProcessInfo);
		//add stablilization terms due to grad term grad(q) * ro*Acce
		CalculateGradMassStblTerms(rLeftHandSideMatrix, rRightHandSideVector, DN_DX,N, tauone, Area, rCurrentProcessInfo);

		CalculateResidual(rLeftHandSideMatrix, rRightHandSideVector);

		
		
		
		//AddPressureMassTermsPrediction(rLeftHandSideMatrix, rRightHandSideVector, Area, rCurrentProcessInfo );
	
	}



        //add body force and momentum
        AddBodyForceAndMomentum(rRightHandSideVector, DN_DX, N, delta_t, Area, tauone, tautwo);

        //add projections
        //if (rCurrentProcessInfo[OSS_SWITCH] == 1)
        //    AddProjectionForces(rRightHandSideVector, DN_DX, Area, tauone, tautwo);


        KRATOS_CATCH("")
    }
    //***********************************************************************************
    //**************************************************************************************

    void VP_PRECOND2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY

        MatrixType temp = Matrix();
        CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::CalculateMassContribution(MatrixType& K, const double time, const double area, const ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY
        double lump_mass_fac = area * 0.333333333333333333333333;
        double density;
        double mu;
        calculatedensity(GetGeometry(), density, mu);

	double delta_t = rCurrentProcessInfo[DELTA_TIME];

        int nodes_number = 3; 
        int dof = 2;

	int FractionalStepNumber =  rCurrentProcessInfo[FRACTIONAL_STEP];
	double factor=1.0;
	//double alpha_bossak=-0.2;

	if (FractionalStepNumber==2)
		factor=1.0/delta_t;
		
        for (int nd = 0; nd < nodes_number; nd++) {
            int row = nd * (dof + 1);
            for (int jj = 0; jj < dof; jj++)
                K(row + jj, row + jj) += density  * lump_mass_fac * factor;

		//add pressure mass
		double SV2=10.0;
		
		if (FractionalStepNumber==1)
			K(row + dof , row + dof ) += (1.0/SV2)*lump_mass_fac;
		
        }

        KRATOS_CATCH("")

    }
    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY
	unsigned int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];
	if (FractionalStepNumber==1)
	{
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
		CalculateTau(N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);

		CalculateMassContribution(rMassMatrix, delta_t, Area, rCurrentProcessInfo);	
	
		Vector dummy=ZeroVector(9);
		//add stablilization terms due to advective term (a)grad(V) * ro*Acce
		CalculateAdvMassStblTerms(rMassMatrix, dummy, DN_DX, N, tauone, Area, rCurrentProcessInfo);
		//add stablilization terms due to grad term grad(q) * ro*Acce
		//next one stabilizaes the continuity equation
		CalculateGradMassStblTerms(rMassMatrix, dummy, DN_DX,N, tauone, Area, rCurrentProcessInfo);
	}
        KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::CalculateLocalVelocityContribution(MatrixType& rDampingMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY
        unsigned int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

	if (FractionalStepNumber==1)
	{
		int nodes_number = 3;
		int dim = 2;
		unsigned int matsize = nodes_number * (dim + 1);

		if (rDampingMatrix.size1() != matsize)
		    rDampingMatrix.resize(matsize, matsize, false); //false says not to preserve existing storage!!


		noalias(rDampingMatrix) = ZeroMatrix(matsize, matsize);
		double delta_t = rCurrentProcessInfo[DELTA_TIME];

		boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
		array_1d<double, 3 > N; 
	
		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

		//viscous term
		CalculateViscousTerm(rDampingMatrix, DN_DX, Area, rCurrentProcessInfo);

		//Advective term
		double tauone;
		double tautwo;
		CalculateTau(N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);

		CalculateAdvectiveTerm(rDampingMatrix, DN_DX,N, tauone, tautwo, delta_t, Area);

		//calculate pressure term
		CalculatePressureTerm(rDampingMatrix, DN_DX, N, delta_t, Area, rCurrentProcessInfo);

		//compute projections
		//unsigned int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];
		//add stabilization only in the correction step
				
		//stabilization terms
		CalculateDivStblTerm(rDampingMatrix, DN_DX, tautwo, Area);
		CalculateAdvStblAllTerms(rDampingMatrix, rRightHandSideVector, DN_DX, N, tauone, delta_t, Area, rCurrentProcessInfo);
		CalculateGradStblAllTerms(rDampingMatrix, rRightHandSideVector, DN_DX,N, delta_t, tauone, Area);
		//KRATOS_WATCH(rRightHandSideVector);
	
		
		CalculateResidual(rDampingMatrix, rRightHandSideVector);
		//AddPressureMassTermsPrediction(rDampingMatrix, rRightHandSideVector, Area, rCurrentProcessInfo );
	
	}
        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::CalculateViscousTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const double area, const ProcessInfo& CurrentProcessInfo) {
        KRATOS_TRY
        
	double mu;        
        double density;
        calculatedensity(GetGeometry(), density, mu);

	//unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

	//Matrix TEMP=ZeroMatrix(9,9);
        //nu = nu/density;

        int nodes_number = 3;
        int dof = 2;

        for (int ii = 0; ii < nodes_number; ii++) {
            int row = ii * (dof + 1);
            for (int jj = 0; jj < nodes_number; jj++) {
                int column = jj * (dof + 1);
                K(row, column) += mu * area * (DN_DX(ii, 0) * DN_DX(jj, 0) + DN_DX(ii, 1) * DN_DX(jj, 1));
                K(row + 1, column + 1) += mu * area * (DN_DX(ii, 0) * DN_DX(jj, 0) + DN_DX(ii, 1) * DN_DX(jj, 1));
            }
        }
	
	
        KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::CalculateAdvectiveTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double, 3 > & N, const double tauone, const double tautwo, const double time, const double area) {
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
        boost::numeric::ublas::bounded_matrix<double, 6, 6 > temp_convterm = ZeroMatrix(matsize, matsize);
        temp_convterm = prod(shape_func, conv_opr);

        //double fac = tauone/time;
        //fac = 0.0;
        //temp_convterm *= ((1 + fac*density)*density); // For the simplicity of implementation the stabilization term tau1*ro/deltat*(a.gradV U(n+1,t)) is added

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

    void VP_PRECOND2D::CalculatePressureTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double, 3 > & N, const double time, const double area, const ProcessInfo& CurrentProcessInfo) {
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
	        K(row, column) -=  area * N(jj) * DN_DX(ii, 0);
	        K(column, row) +=  area *density * N(jj) * DN_DX(ii, 0);

	        K(row + 1, column) -=  area * N(jj) * DN_DX(ii, 1);
	        K(column, row + 1) +=  area * density * N(jj) * DN_DX(ii, 1);
	    }
	}

	

        KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::CalculateDivStblTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const double tautwo, const double area) {
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

                K(row, column) +=  area * density * temp_div(loc_row, loc_column);
                K(row, column + 1) +=  area * density * temp_div(loc_row, loc_column + 1);
                K(row + 1, column) +=  area * density * temp_div(loc_row + 1, loc_column);
                K(row + 1, column + 1) +=  area * density * temp_div(loc_row + 1, loc_column + 1);
            }
        }

        KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::CalculateAdvStblAllTerms(MatrixType& K, VectorType& F, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double, 3 > & N, const double tauone, const double time, const double area, const ProcessInfo& CurrentProcessInfo) {
        KRATOS_TRY

        //unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
	array_1d<double, 2 > ms_adv_vel;
	
                const array_1d<double, 3 > & adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 0);
        const array_1d<double, 3 > & mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY, 0);
        const array_1d<double, 3 > & mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, 0);
        const array_1d<double, 3 > & mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);

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

                K(row, column) +=  area * density * adv_stblterm(loc_row, loc_column);
                K(row, column + 1) +=  area * density * adv_stblterm(loc_row, loc_column + 1);
		
                K(row + 1, column) +=  area * density * adv_stblterm(loc_row + 1, loc_column);
                K(row + 1, column + 1) +=  area * density * adv_stblterm(loc_row + 1, loc_column + 1);
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

                K(row, column) +=  area * 1.0 * grad_stblterm(loc_row, jj);
                K(row + 1, column) +=  area * 1.0 * grad_stblterm(loc_row + 1, jj);
		
		//if (FractionalStepNumber==2) //THIS TERM IS ESSENTIAL FOR GETTING THE CORRECT SOLUTION IN THE PREDICTION STEP!
		//{					
		K(column, row) +=  area * density * grad_stblterm(loc_row, jj);
		K(column, row + 1) +=  area * density * grad_stblterm(loc_row + 1, jj);
		//}		
		
            }
        }

        
        //build (1.0*a.grad V) (Fbody) stabilization term & assemble
        array_1d<double, 2 > bdf = ZeroVector(2);
        const array_1d<double, 2 > bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
        const array_1d<double, 2 > bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
        const array_1d<double, 2 > bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);


        bdf[0] = N[0]*(density * bdf0[0]) + N[1]*(density * bdf1[0]) + N[2]*(density * bdf2[0]);
        bdf[1] = N[0]*(density * bdf0[1]) + N[1]*(density * bdf1[1]) + N[2]*(density * bdf2[1]);


        array_1d<double, 6 > fbd_stblterm = ZeroVector(matsize);
        fbd_stblterm = tauone *  prod(trans(conv_opr), bdf);

        for (int ii = 0; ii < nodes_number; ++ii) {
            int index = ii * (dof + 1);
            int loc_index = ii*dof;
            F[index] +=  area * fbd_stblterm[loc_index];
            F[index + 1] +=  area * fbd_stblterm[loc_index + 1];
        }
        KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::CalculateAdvMassStblTerms(MatrixType& M, VectorType& F, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double, 3 > & N, const double tauone, const double area, const ProcessInfo& CurrentProcessInfo) {
        KRATOS_TRY
        
	array_1d<double, 2 > ms_adv_vel;

	    const array_1d<double, 3 > & adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 0);
        const array_1d<double, 3 > & mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY, 0);
        const array_1d<double, 3 > & mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, 0);
        const array_1d<double, 3 > & mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);

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

	//Pavel

	double factor=1.0;
	double delta_t=CurrentProcessInfo[DELTA_TIME];
	int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
	//double SV2=10.0;
	if (FractionalStepNumber==1)
		factor=1.0;
	//double alpha_bossak=-0.2;
	if (FractionalStepNumber==2)
		factor=1.0/delta_t;

	boost::numeric::ublas::bounded_matrix<double, 9, 9 > TEMP = ZeroMatrix(9, 9);
	
        for (int ii = 0; ii < nodes_number; ii++) {
            int row = ii * (dof + 1);
            int loc_row = ii*dof;
            for (int jj = 0; jj < nodes_number; jj++) {
                int column = jj * (dof + 1);
                int loc_column = jj*dof;

                M(row, column) += area * fac * temp_convterm(loc_row, loc_column) * factor;
		TEMP(row, column)+=area * fac * temp_convterm(loc_row, loc_column) * factor;

                M(row + 1, column + 1) += area * fac * temp_convterm(loc_row + 1, loc_column + 1)* factor;
                TEMP(row + 1, column + 1) += area * fac * temp_convterm(loc_row + 1, loc_column + 1)* factor;
            }
        }

	//And now the for time integration in the second step: Pavel: -> we have to add also the RHS term by hand (as we are using linear static strategy)
	if (FractionalStepNumber==2)
		{

		//Vector Old_Vels=ZeroVector(9);
		array_1d<double, 9 > Old_Vels;
		for (int i=0;i<3;i++)
		{
			int index=i*3;
			Old_Vels[index]=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X, 1);
			Old_Vels[index+1]=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y, 1);

		}
		F+=prod(TEMP, Old_Vels);
		}


        KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::CalculateGradStblAllTerms(MatrixType& K, VectorType& F, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double,3>& N, const double time, const double tauone, const double area) {
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

                K(row, column) += area *  gard_opr(ii, jj);

            }
        }
       


        //build 1*(grad q) (Fbody ) stabilization term & assemble
        array_1d<double, 2 > bdf = ZeroVector(2);
        const array_1d<double, 2 > bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
        const array_1d<double, 2 > bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
        const array_1d<double, 2 > bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);


        bdf[0] = N[0]*(density * bdf0[0]) + N[1]*(density * bdf1[0]) + N[2]*(density * bdf2[0]);
        bdf[1] = N[0]*(density * bdf0[1]) + N[1]*(density * bdf1[1]) + N[2]*(density * bdf2[1]);

        array_1d<double, 3 > fbd_stblterm = ZeroVector(nodes_number);
        fbd_stblterm = tauone * prod(DN_DX, bdf);


        for (int ii = 0; ii < nodes_number; ++ii) {
            int index = ii * (dof + 1) + dof;
            F[index] +=  area * fbd_stblterm[ii];
        }





        KRATOS_CATCH("")

    }
    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::CalculateGradMassStblTerms(MatrixType& M, VectorType& F, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N,const double tauone,const double area, const ProcessInfo& CurrentProcessInfo) {
        KRATOS_TRY
                int nodes_number = 3;
        int dof = 2;

	unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
	double delta_t = CurrentProcessInfo[DELTA_TIME];


        double density;
        double mu;
        calculatedensity(GetGeometry(), density, mu);

        //build 1*tau1*ro Nacc grad q)
        double fac = tauone*density;
	
	//double SV2=100.0;
	double factor=1.0;
	//if (FractionalStepNumber==1)
	//	factor=1.0/SV2;

	if (FractionalStepNumber==2)
		factor=1.0/(delta_t);

	boost::numeric::ublas::bounded_matrix<double, 9, 9 > TEMP = ZeroMatrix(9,9);

        for (int ii = 0; ii < nodes_number; ii++) {
            int row = ii * (dof + 1);
            for (int jj = 0; jj < nodes_number; jj++) {
                int column = jj * (dof + 1) + dof;

                //K(row,column) += -1*area * fac* N(ii) * DN_DX(jj,0);
                M(column, row) += area * fac * N[ii] * DN_DX(jj, 0)* factor;
                TEMP(column, row) += area * fac * N[ii] * DN_DX(jj, 0)* factor;

                //K(row + 1,column) += -1*area * fac* N(ii) * DN_DX(jj,1);
                M(column, row + 1) += area * fac * N[ii] * DN_DX(jj, 1)* factor;
                TEMP(column, row + 1) += area * fac * N[ii] * DN_DX(jj, 1)* factor;
            }
        }	      
	
	//And now the for time integration in the second step: Pavel: -> we have to add also the RHS term by hand (as we are using linear static strategy)
	if (FractionalStepNumber==2)
		{
	//Vector Old_Vels=ZeroVector(9);
		array_1d<double,9> Old_Vels = ZeroVector(9);
		for (int i=0;i<3;i++)
		{
			int index=i*3;
			Old_Vels[index]=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X, 1);
			Old_Vels[index+1]=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y, 1);

		}
		F+=prod(TEMP, Old_Vels);
		}




        KRATOS_CATCH("")

    }
    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::AddBodyForceAndMomentum(VectorType& F,const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX,  const array_1d<double, 3 > & N, const double time, const double area, const double tauone, const double tautwo) {
        KRATOS_TRY
                int nodes_number = 3;
        int dof = 2;


        //double lump_mass_fac = area * 0.333333333333333333333333;

        double density;
        double mu;
        calculatedensity(GetGeometry(), density, mu);

        //for Arhenious
        int matsize = dof*nodes_number;
        boost::numeric::ublas::bounded_matrix<double, 1, 6 > div_opr = ZeroMatrix(1, matsize);
        for (int ii = 0; ii < nodes_number; ii++) {
            int index = dof*ii;
            div_opr(0, index) = DN_DX(ii, 0);
            div_opr(0, index + 1) = DN_DX(ii, 1);
        }
        const double ar_0 = GetGeometry()[0].FastGetSolutionStepValue(ARRHENIUS);
        const double ar_1 = GetGeometry()[1].FastGetSolutionStepValue(ARRHENIUS);
        const double ar_2 = GetGeometry()[2].FastGetSolutionStepValue(ARRHENIUS);

        double mean_ar = 0.333333333333333333 * (ar_0 + ar_1 + ar_2);


        //body  & momentum term force
        for (int ii = 0; ii < nodes_number; ii++) {
            int index = ii * (dof + 1);
            int loc_index = ii * dof;
            const array_1d<double, 2 > bdf = GetGeometry()[ii].FastGetSolutionStepValue(BODY_FORCE);


            F[index] += area * N[ii] * density * bdf[0];
            F[index + 1] += area * N[ii] * density * bdf[1];


            //arrhenius
            F[index + 2] += (area * N[ii] * mean_ar);
            F[index] += tautwo * area * mean_ar * div_opr(0, loc_index);
            F[index + 1] += tautwo * area * mean_ar * div_opr(0, loc_index + 1);
// KRATOS_WATCH(density);
        }


        KRATOS_CATCH("")
    }
	

    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::CalculateResidual(const MatrixType& K, VectorType& F) {
        KRATOS_TRY

                int nodes_number = 3;
        int dof = 2;


        array_1d<double, 9 > UP = ZeroVector(9);
        for (int ii = 0; ii < nodes_number; ++ii) {
            int index = ii * (dof + 1);
            UP[index] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[0];
            UP[index + 1] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[1];
            UP[index + 2] = GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE, 0);
        }

        F -= prod(K, UP);
	//KRATOS_WATCH(F)
        KRATOS_CATCH("")
    }
   
    

    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) {
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

    void VP_PRECOND2D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {
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

    void VP_PRECOND2D::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
            array_1d<double, 3 > & Output,
            const ProcessInfo& rCurrentProcessInfo) {

        array_1d<double, 6 > adv_proj = ZeroVector(6);
        array_1d<double, 3 > div_proj = ZeroVector(3);
	boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N;
	
        double delta_t = rCurrentProcessInfo[DELTA_TIME];

        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

        double tauone;
        double tautwo;
        CalculateTau(N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);

        //ComputeProjections(adv_proj, div_proj, DN_DX, tauone, tautwo, N, Area, delta_t);
        

    }

    //************************************************************************************
    //************************************************************************************

    void VP_PRECOND2D::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) {

        double delta_t = rCurrentProcessInfo[DELTA_TIME];
	boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N;
	
        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);
        double tauone;
        double tautwo;
        CalculateTau(N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);
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

    }
    //*************************************************************************************
    //*************************************************************************************

    void VP_PRECOND2D::calculatedensity(Geometry< Node < 3 > > geom, double& density, double& viscosity) {

        double kk = 0.0;
	density = 0.0;
	viscosity = 0.0;
        for(int ii=0;ii<3;++ii)
                        {
                                kk++;
                                density +=geom[ii].FastGetSolutionStepValue(DENSITY);
				viscosity +=geom[ii].FastGetSolutionStepValue(VISCOSITY);
                        }

        density/=kk;
	viscosity/=kk	;

      

	//Here we calculate Dynamic viscosity from Kinemeatic viscosity
	viscosity *= density;

    }
    //*************************************************************************************
    //*************************************************************************************

    void VP_PRECOND2D::CalculateTau(const array_1d<double,3>& N, double& tauone, double& tautwo, const double time, const double area, const ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY
	//unsigned int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];
	
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
        //const double mu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
        //const double mu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
        //const double mu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
        //mu = 0.333333333333333333333333*(mu0 + mu1 + mu2);

        double density;
        calculatedensity(GetGeometry(), density, mu);

        const double dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
        tauone = 1.0 / (dyn_st_beta / time + 4.0 * mu / (ele_length * ele_length * density) + 2.0 * advvel_norm  / ele_length) ;

        tautwo = mu / density + 1.0 * ele_length * advvel_norm / 2.0;
	

        KRATOS_CATCH("")

 
    }


    //*************************************************************************************
    //*************************************************************************************

    void VP_PRECOND2D::GetFirstDerivativesVector(Vector& values, int Step) {
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

    void VP_PRECOND2D::GetSecondDerivativesVector(Vector& values, int Step) {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = number_of_nodes * (dim + 1);
        if (values.size() != MatSize) values.resize(MatSize, false);
        for (unsigned int i = 0; i < number_of_nodes; i++) {
            unsigned int index = i * (dim + 1);
            values[index] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X, Step);
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y, Step);
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue(PRESSURE_DT,Step);
        }
    }
    //************************************************************************************
    //************************************************************************************
} // Namespace Kratos


