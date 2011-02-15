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
#include "custom_elements/explicit_asgs_compressible_2d.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{

	//************************************************************************************
	//************************************************************************************
	ExplicitASGSCompressible2D::ExplicitASGSCompressible2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: ASGS2D(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
		
	}

	//************************************************************************************
	//************************************************************************************
	ExplicitASGSCompressible2D::ExplicitASGSCompressible2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: ASGS2D(NewId, pGeometry, pProperties)
	{
			
	}

	Element::Pointer ExplicitASGSCompressible2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		
		KRATOS_TRY
		return Element::Pointer(new ExplicitASGSCompressible2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	ExplicitASGSCompressible2D::~ExplicitASGSCompressible2D()
	{
	}

    //************************************************************************************
    //************************************************************************************

    void ExplicitASGSCompressible2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY

                int nodes_number = 3;
        int dim = 2;
        unsigned int matsize = nodes_number * (dim + 1);

        if (rLeftHandSideMatrix.size1() != matsize)
            rLeftHandSideMatrix.resize(matsize, matsize,false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != matsize)
            rRightHandSideVector.resize(matsize,false); //false says not to preserve existing storage!!


        noalias(rLeftHandSideMatrix) = ZeroMatrix(matsize, matsize);
        noalias(rRightHandSideVector) = ZeroVector(matsize);

        double delta_t = rCurrentProcessInfo[DELTA_TIME];



        //getting data for the given geometry
        double Area;
	     boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2); 
        array_1d<double,3> N = ZeroVector(3);
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

        double tauone;
        double tautwo;
        CalculateTau(N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);


        //add body force and momentum
        AddBodyForceAndMomentum(rRightHandSideVector, DN_DX, N, delta_t, Area, tauone, tautwo);


        KRATOS_CATCH("")
    }
	//************************************************************************************
	//************************************************************************************
	void ExplicitASGSCompressible2D::CalculateLocalVelocityContribution(MatrixType& rDampMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
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
        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
        array_1d<double,3> N = ZeroVector(3);
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);
	

	//viscous term	
	CalculateViscousTerm(rDampMatrix, DN_DX, Area);
// 	CalcualteDCOperatior(rDampMatrix, DN_DX, Area);
	//Advective term
	double tauone;
	double tautwo;
	CalculateTau(N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);



	//CalculateAdvectiveTerm(rDampMatrix, DN_DX, tauone, tautwo, delta_t, Area);

	//calculate pressure term
	CalculatePressureTerm(rDampMatrix, DN_DX, N, delta_t,Area);

	//compute projections
	//stabilization terms
	//KRATOS_WATCH("ñññññññññ calculate stabilizing terms ñññññññññ");
	 double VC2;
	 CalculateSoundVelocity(GetGeometry(), VC2);
     double tau_div = tautwo*VC2;

	CalculateDivStblTerm(rDampMatrix, DN_DX, tau_div, Area);
	//CalculateAdvStblAllTerms(rDampMatrix,rRightHandSideVector, DN_DX, N, tauone,delta_t, Area);
	CalculateGradStblAllTerms(rDampMatrix,rRightHandSideVector,DN_DX,N, delta_t, tauone, Area);

	CalculateDivPdotStblTerms(rDampMatrix,rRightHandSideVector, DN_DX,N,  delta_t, tautwo, Area);


	CalculateResidual(rDampMatrix, rRightHandSideVector);

		KRATOS_CATCH("")
	}
   
	//************************************************************************************
	//************************************************************************************
	void ExplicitASGSCompressible2D::CalculateMassContribution(MatrixType& K,const double time,const double area)
	{
		KRATOS_TRY
	double lump_mass_fac = area * 0.333333333333333333333333;
	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	//update density and calculate sound velocity

	int nodes_number = 3;
	int dof = 2;
	for ( int nd = 0; nd< nodes_number; nd++)
	    {
		int row = nd*(dof + 1);
		for( int jj=0; jj< dof; jj++)
			K(row + jj, row + jj) += density*lump_mass_fac;

		//add pressure mass
			K(row + dof , row + dof ) += lump_mass_fac;
	    }
	
	
		KRATOS_CATCH("")
	
	}
	//************************************************************************************
	//************************************************************************************
	void ExplicitASGSCompressible2D::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			//lumped
			unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int NumberOfNodes = GetGeometry().size();
		unsigned int MatSize = (dimension + 1) * NumberOfNodes;
		if(rMassMatrix.size1() != MatSize)
			rMassMatrix.resize(MatSize,MatSize,false);

		noalias(rMassMatrix) = ZeroMatrix(MatSize,MatSize);
	double delta_t= rCurrentProcessInfo[DELTA_TIME];
		

	//getting data for the given geometry
	double Area;
        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
        array_1d<double,3> N = ZeroVector(3);
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);

	//Calculate tau
	//double tauone;
	//double tautwo;
	//CalculateTau(N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);

	CalculateMassContribution(rMassMatrix,delta_t,Area); 
		//add stablilization terms due to advective term (a)grad(V) * ro*Acce
	//CalculateAdvMassStblTerms(rMassMatrix, DN_DX, N,tauone,Area);
		//add stablilization terms due to grad term grad(q) * ro*Acce
	//CalculateGradMassStblTerms(rMassMatrix, DN_DX, tauone,Area);
		//add compressible stabilization terms
	//CalculateCompressibleStblTerms(rMassMatrix, DN_DX,N, tautwo,Area);
	
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
    void ExplicitASGSCompressible2D::CalculatePressureTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double, 3 > & N, const double time, const double area) {
        KRATOS_TRY
                int nodes_number = 3;
        int dof = 2;

        double density;
        double mu;
        calculatedensity(GetGeometry(), density, mu);

	double VC2;
	CalculateSoundVelocity(GetGeometry(), VC2);
// KRATOS_WATCH(VC2);
        for (int ii = 0; ii < nodes_number; ii++) {
            int row = ii * (dof + 1);
            for (int jj = 0; jj < nodes_number; jj++) {
                int column = jj * (dof + 1) + dof;

                K(row, column) -=  area * N(jj) * DN_DX(ii, 0);
                K(column, row) += VC2 * area * density * N(jj) * DN_DX(ii, 0);

                K(row + 1, column) -=  area * N(jj) * DN_DX(ii, 1);
                K(column, row + 1) += VC2 * area * density * N(jj) * DN_DX(ii, 1);
            }
        }


        KRATOS_CATCH("")
    }


	//*************************************************************************************
	//*************************************************************************************
    void ExplicitASGSCompressible2D::CalculateGradStblAllTerms(MatrixType& K, VectorType& F, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX,const array_1d<double,3>& N, const double time, const double tauone, const double area) {
        KRATOS_TRY
                int nodes_number = 3;
        int dof = 2;


        double density;
        double mu;
        calculatedensity(GetGeometry(), density, mu);
	
	double VC2;
	CalculateSoundVelocity(GetGeometry(), VC2);	

        //build 1*(grad q . grad p) stabilization term & assemble
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > gard_opr = ZeroMatrix(nodes_number, nodes_number);
        gard_opr =  prod(DN_DX, trans(DN_DX));

        for (int ii = 0; ii < nodes_number; ii++) {
            int row = ii * (dof + 1) + dof;

            for (int jj = 0; jj < nodes_number; jj++) {
                int column = jj * (dof + 1) + dof;

                K(row, column) += area * tauone * gard_opr(ii, jj);

            }
        }

        //build 1*(grad q) (Fbody ) stabilization term & assemble
        array_1d<double, 2 > bdf = ZeroVector(2);
        const array_1d<double, 3 > bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
        const array_1d<double, 3 > bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
        const array_1d<double, 3 > bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);

        //add acceleration
        const array_1d<double,3>& acce0 = GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION);
        const array_1d<double,3>& acce1 = GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION);
        const array_1d<double,3>& acce2 = GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION);


        bdf[0] = N[0]*(density * (bdf0[0] - acce0[0])) + N[1]*(density * (bdf1[0] - acce1[0])) + N[2]*(density * (bdf2[0] - acce2[0]));
        bdf[1] = N[0]*(density * (bdf0[1] - acce0[1])) + N[1]*(density * (bdf1[1] - acce1[1])) + N[2]*(density * (bdf2[1] - acce2[1]));


        array_1d<double, 3 > fbd_stblterm = ZeroVector(nodes_number);
        fbd_stblterm = tauone * prod(DN_DX, bdf);


        for (int ii = 0; ii < nodes_number; ++ii) {
            int index = ii * (dof + 1) + dof;
            F[index] +=  area * fbd_stblterm[ii];
        }

        KRATOS_CATCH("")
    }
	//*************************************************************************************
	//*************************************************************************************
	void ExplicitASGSCompressible2D::CalculateDivPdotStblTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const array_1d<double,3>& N, const double time,const double tautwo,const double area)
	  {
	    KRATOS_TRY
	//tau*div(V).P_dot
	  double stbl_fac = tautwo * area;
          int nodes_number = 3;
          int dof = 2;

         double mean_pressure_rate = N(0) * (GetGeometry()[0].FastGetSolutionStepValue(AIR_PRESSURE_DT));
	 for( int ii=1; ii < nodes_number; ++ii)
            mean_pressure_rate +=  N(ii) * (GetGeometry()[ii].FastGetSolutionStepValue(AIR_PRESSURE_DT));

	 mean_pressure_rate *= stbl_fac;

        for (int ii = 0; ii < nodes_number; ++ii) {
            int index = ii * (dof + 1);
            F[index] -=  DN_DX(ii,0) * mean_pressure_rate;
            F[index+1] -=  DN_DX(ii,1) * mean_pressure_rate;
        }
// 	  double VC2;
// 	  CalculateSoundVelocity(GetGeometry(), VC2);
// 	  double stbl_fac = tautwo/VC2 * volume;
// 	int nodes_number = 4;
// 	int dof = 3;
// 	for ( int ii = 0; ii < nodes_number; ii++)
// 	    {
// 		//LHS contribution
// 		int row = ii*(dof+1);
// 		for( int jj=0; jj < nodes_number; jj++)
// 		   {
// 			int column = jj*(dof+1) + dof;
// 
// 			K(row,column) +=  stbl_fac* N(jj) * DN_DX(ii,0);
// 			K(row + 1,column) +=  stbl_fac* N(jj) * DN_DX(ii,1);
// 			K(row + 2,column) +=  stbl_fac* N(jj) * DN_DX(ii,2);
// 
// 		   }
// 	    }

        KRATOS_CATCH("")
	  }
    	//*************************************************************************************
	//*************************************************************************************
	  void ExplicitASGSCompressible2D::GetFirstDerivativesVector(Vector& values, int Step)
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
			values[index + 2] = GetGeometry()[i].GetSolutionStepValue(AIR_PRESSURE,Step);

		}
	}
	//************************************************************************************
	//************************************************************************************
	  void ExplicitASGSCompressible2D::GetSecondDerivativesVector(Vector& values, int Step)
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
			values[index + 2] = GetGeometry()[i].GetSolutionStepValue(AIR_PRESSURE_DT,Step);
		}
	
	}
	//************************************************************************************
	//************************************************************************************
	void ExplicitASGSCompressible2D::CalculateSoundVelocity(Geometry< Node<3> > geom, double& vc2)
	{

	/*	vc2 = 1441.000 * 1441.000;
	  double K1 = 2070000000;
          double K2 = 7.15;

	double mean_rho_w = 0.0;
	double mean_old_rho_w = 0.0;
	double mean_old_pr_w = 0.0;

	for(int ii = 0; ii<3; ii++)
	   {
		mean_rho_w += GetGeometry()[ii].FastGetSolutionStepValue(DENSITY );
		mean_old_rho_w += GetGeometry()[ii].FastGetSolutionStepValue(DENSITY,1 );
		mean_old_pr_w += GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE,1);
	   }
	
	mean_rho_w *= 0.333333333333333333333333333333333333;
	mean_old_rho_w *=0.333333333333333333333333333333333333;
	mean_old_pr_w *=0.333333333333333333333333333333333333;


	double alpha = (mean_old_pr_w * K2 + K1)/mean_old_rho_w;
	
	vc2 = alpha*pow(mean_rho_w/mean_old_rho_w, K2-1.0);*/


	 double mean_vc2=0.0;
	 double mean_vc=0.0;
		for(int ii = 0; ii<3; ii++)
		    {
			mean_vc = GetGeometry()[ii].FastGetSolutionStepValue(AIR_SOUND_VELOCITY ) ;
			mean_vc2 += mean_vc * mean_vc;
		    }

	 vc2 =mean_vc2*0.333333333333333333333333;

vc2 = 10.0;
	}


	//************************************************************************************
	//************************************************************************************
	void ExplicitASGSCompressible2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
	  void ExplicitASGSCompressible2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
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

          void ExplicitASGSCompressible2D::CalculateTau(const array_1d<double,3>& N,double& tauone, double& tautwo, const double time, const double area, const ProcessInfo& rCurrentProcessInfo)
          {

            KRATOS_TRY
        array_1d<double, 2 > ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension

        ms_adv_vel[0] = 0.0;
        ms_adv_vel[1] = 0.0;


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

        double VC2;
	CalculateSoundVelocity(GetGeometry(), VC2);
	double vc = sqrt(VC2);

	double int_time = ele_length/vc;

        const double dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
        //tauone = 1.0 / (dyn_st_beta / time + 4.0 * mu / (ele_length * ele_length * density) + 2.0 * advvel_norm  / ele_length);
       // tauone = 1.0 / (dyn_st_beta / int_time + 4.0 * mu / (ele_length * ele_length * density) + 2.0 * advvel_norm  / ele_length);
         //tauone = time*VC2;
            tauone = time;
            tautwo = time;
	    
	    
            KRATOS_CATCH("")

              }
    
	//*************************************************************************************
	//*************************************************************************************
	void ExplicitASGSCompressible2D::calculatedensity(Geometry< Node<3> > geom, double& density, double& viscosity)
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

  //density = 10.0;

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
        //************************************************************************************
    //************************************************************************************

    void ExplicitASGSCompressible2D::AddBodyForceAndMomentum(VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX, const array_1d<double, 3 > & N, const double time, const double area, const double tauone, const double tautwo) {
        KRATOS_TRY
                int nodes_number = 3;
        int dof = 2;


        //double lump_mass_fac = area * 0.333333333333333333333333;

        double density;
        double mu;
        calculatedensity(GetGeometry(), density, mu);
//KRATOS_WATCH(density);
        //body  & momentum term force
        for (int ii = 0; ii < nodes_number; ii++) {
            int index = ii * (dof + 1);
            int loc_index = ii * dof;
            const array_1d<double, 2 > bdf = GetGeometry()[ii].FastGetSolutionStepValue(BODY_FORCE);


            F[index] += area * N[ii] * density * bdf[0];
            F[index + 1] += area * N[ii] * density * bdf[1];

        }


        KRATOS_CATCH("")
    } 
    //************************************************************************************
    //************************************************************************************

    void ExplicitASGSCompressible2D::CalculateResidual(const MatrixType& K, VectorType& F) 
	{
	    KRATOS_TRY

		    int nodes_number = 3;
	    int dof = 2;


	    array_1d<double, 9 > UP = ZeroVector(9);
	    for (int ii = 0; ii < nodes_number; ++ii) {
		int index = ii * (dof + 1);
		UP[index] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[0];
		UP[index + 1] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[1];
		UP[index + 2] = GetGeometry()[ii].FastGetSolutionStepValue(AIR_PRESSURE, 0);
	    }

	    F -= prod(K, UP);

	    KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************

    void ExplicitASGSCompressible2D::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
            array_1d<double, 3 > & Output,
            const ProcessInfo& rCurrentProcessInfo) 
        {
	  Output = ZeroVector(3);

	double Area;
        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
        array_1d<double,3> N = ZeroVector(3);
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);


	  double lump_mass_fac = Area * 0.333333333333333333333333;
	  double density;
	  double mu;

	  calculatedensity(GetGeometry(), density, mu);

	  //fill velocity and pressure mass

	  Output[0] = density*lump_mass_fac;
	  Output[1] = lump_mass_fac;
//KRATOS_WATCH(Output[0]);

       }
	//************************************************************************************
	//************************************************************************************
        void ExplicitASGSCompressible2D::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
       {

	   Output = 100.0;
	   double Area = 0.0;

	  boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX = ZeroMatrix(3, 2);
	  array_1d<double, 3 > N = ZeroVector(3); //dimension = number of nodes
	  GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

	double calc_t = 0.0;
        for (unsigned int ii = 0; ii < 3; ii++)
        {
            double inv_max_h = DN_DX(ii,0)*DN_DX(ii,0) + DN_DX(ii,1)*DN_DX(ii,1);
            double VC = GetGeometry()[ii].FastGetSolutionStepValue(AIR_SOUND_VELOCITY ) ;
	    inv_max_h = sqrt(inv_max_h);
	    
VC = 10.0;	    
	    calc_t = 1.0/(inv_max_h * VC );


	    if( calc_t < Output)
		  Output = calc_t;
              
        }
//KRATOS_WATCH(Output);

        }
        //************************************************************************************
    //************************************************************************************
//      void ExplicitASGSCompressible2D::CalcualteDCOperatior(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double area)
//       {
//         KRATOS_TRY 
//         
// 	double art_visc = 0.0;
// 	CalculateArtifitialViscosity(art_visc,DN_DX);
// 	
//         double mu;
//         double density;
//         calculatedensity(GetGeometry(), density, mu);
// 
//        double fac = art_visc * density * area;
//         //nu = nu/density;
// 
//         int nodes_number = 3;
//         int dof = 2;
// 
//         for (int ii = 0; ii < nodes_number; ii++) {
//             int row = ii * (dof + 1);
//             for (int jj = 0; jj < nodes_number; jj++) {
//                 int column = jj * (dof + 1);
//                 K(row, column) +=  fac * (DN_DX(ii, 0) * DN_DX(jj, 0) + 0.5 * DN_DX(ii, 1) * DN_DX(jj, 1));
//                 K(row, column + 1) +=  fac * 0.5 * DN_DX(ii, 1) * DN_DX(jj, 0);		
// 		K(row + 1, column ) +=  fac * 0.5 * DN_DX(ii, 0) * DN_DX(jj, 1);	
//                 K(row + 1, column + 1) +=  fac * (DN_DX(ii, 1) * DN_DX(jj, 1) + 0.5*DN_DX(ii, 0) * DN_DX(jj, 0));
//             }
//         }
// 	    KRATOS_CATCH("")
//       }
    //************************************************************************************
    //************************************************************************************
//      void ExplicitASGSCompressible2D::CalculateArtifitialViscosity(double& art_visc ,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX)
//  	{
// 	    KRATOS_TRY    
// 	  
//          int nodes_number = 3;
//          for( int ii=1; ii < nodes_number; ++ii)
// 	 {
// 	    const array_1d<double,3>& vel = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY);
// 	    
// 	 }  
// 	  double div_vel = DN_DX(0,0) + 
// 
// 
// 	  double H=0.0;
//           CalculateCharectristicLength(H);
// 	  
// 
// 	  
// 	  
//      
// 	    KRATOS_CATCH("")
// 	}  
	
    //************************************************************************************
    //************************************************************************************
/*     void ExplicitASGSCompressible2D::CalculateCharectristicLength(double& ch_length, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX )
 	{
	    KRATOS_TRY 
	    
  	  GeometryType::JacobiansType J;
	  GetGeometry().Jacobian(J); 

	  Matrix CC;
	  CC.resize(2,2,false);	  	
	  
	  Matrix DD;
	  DD.resize(2,2,false);	  	  

	  noalias(CC) = prod(J[0],trans(J[0]));	 
	  DD=CC;
	  
	  double det = CC(1,1)*CC(0,0) - CC(0,1)*CC(1,0);
	  
	  if(det == 0.0)
	      	KRATOS_ERROR(std::logic_error,"ZERO DETERMINANT IN ARTIFICIAL VISCOSITY ",det)
	  else 
	     det = 1.0/det;
	  
          double zarf = CC(0,0);	  
	  CC(0,0) = det*CC(1,1);
	  CC(1,1) = det*zarf;
	  CC(0,1) = -CC(0,1)*det;
	  CC(1,0) = -CC(1,0)*det;
	  KRATOS_WATCH( prod(DD,CC));	  

	  nodes_number = 3;
	  array_1d<double,3>& mean_acc =  GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION);
          for (int ii = 0; ii < nodes_number; ii++) {	  
	         const array_1d<double,3>& acce0 = GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION);  
	  }
	    KRATOS_CATCH("")
	}*/	  	    
} // Namespace Kratos


