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
#include "custom_elements/asgs_compressible_2d.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
        namespace ASGSCompressible2Dauxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
        #pragma omp threadprivate(DN_DX)

        array_1d<double,3> N = ZeroVector(3); //dimension = number of nodes
        #pragma omp threadprivate(N)

        array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
        #pragma omp threadprivate(ms_adv_vel)

    }
    using  namespace ASGSCompressible2Dauxiliaries;


	//************************************************************************************
	//************************************************************************************
	ASGSCompressible2D::ASGSCompressible2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: ASGS2D(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
		
	}

	//************************************************************************************
	//************************************************************************************
	ASGSCompressible2D::ASGSCompressible2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: ASGS2D(NewId, pGeometry, pProperties)
	{
			
	}

	Element::Pointer ASGSCompressible2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		
		KRATOS_TRY
		return Element::Pointer(new ASGSCompressible2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	ASGSCompressible2D::~ASGSCompressible2D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void ASGSCompressible2D::CalculateLocalVelocityContribution(MatrixType& rDampMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
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
	
	//Advective term
	double tauone;
	double tautwo;
	CalculateTau(tauone, tautwo, delta_t, Area, rCurrentProcessInfo);

        double VC2;
	CalculateSoundVelocity(GetGeometry(), VC2);

	CalculateAdvectiveTerm(rDampMatrix, DN_DX, tauone, tautwo, delta_t, Area);
		
	//calculate pressure term
	CalculatePressureTerm(rDampMatrix, DN_DX, N, delta_t,Area);

	//compute projections
	
	//stabilization terms
	//KRATOS_WATCH("ñññññññññ calculate stabilizing terms ñññññññññ");
        double tau_div = tautwo*VC2;
	CalculateDivStblTerm(rDampMatrix, DN_DX, tau_div, Area);
	//CalculateAdvStblAllTerms(rDampMatrix,rRightHandSideVector, DN_DX, N, tauone,delta_t, Area);
	CalculateGradStblAllTerms(rDampMatrix,rRightHandSideVector,DN_DX, delta_t, tauone, Area);
	//KRATOS_WATCH(rRightHandSideVector);


        //add body force and momentum
       // AddBodyForceAndMomentum(rRightHandSideVector, N, delta_t, Area, tauone, tautwo);

		KRATOS_CATCH("")
	}
        
	//************************************************************************************
	//************************************************************************************
	void ASGSCompressible2D::CalculateMassContribution(MatrixType& K,const double time,const double area)
	{
		KRATOS_TRY
	double lump_mass_fac = area * 0.333333333333333333333333;
	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	//update density and calculate sound velocity
	//double VC2;
	//CalculateSoundVelocity(GetGeometry(), VC2);

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
	void ASGSCompressible2D::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
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
        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
        array_1d<double,3> N = ZeroVector(3);
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);

	//Calculate tau
	double tauone;
	double tautwo;
	CalculateTau(tauone, tautwo, delta_t, Area, rCurrentProcessInfo);

	CalculateMassContribution(rMassMatrix,delta_t,Area); 
		//add stablilization terms due to advective term (a)grad(V) * ro*Acce
	//CalculateAdvMassStblTerms(rMassMatrix, DN_DX, N,tauone,Area);
		//add stablilization terms due to grad term grad(q) * ro*Acce
	CalculateGradMassStblTerms(rMassMatrix, DN_DX, tauone,Area);
		//add compressible stabilization terms
	CalculateCompressibleStblTerms(rMassMatrix, DN_DX,N, tautwo,Area);
	
		KRATOS_CATCH("")
	}
            //************************************************************************************
    //************************************************************************************

    void ASGSCompressible2D::CalculateGradMassStblTerms(MatrixType& M, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const double tauone, const double area) {
        KRATOS_TRY
                int nodes_number = 3;
        int dof = 2;

     //   double density;
       // double mu;
       // calculatedensity(GetGeometry(), density, mu);

        double VC2;
	CalculateSoundVelocity(GetGeometry(), VC2);

        //build 1*tau1*ro Nacc grad q)
        double fac = area*tauone*VC2;
        for (int ii = 0; ii < nodes_number; ii++) {
            int row = ii * (dof + 1);
            for (int jj = 0; jj < nodes_number; jj++) {
                int column = jj * (dof + 1) + dof;

                //K(row,column) += -1*area * fac* N(ii) * DN_DX(jj,0);
                M(column, row) +=  fac * N[ii] * DN_DX(jj, 0);

                //K(row + 1,column) += -1*area * fac* N(ii) * DN_DX(jj,1);
                M(column, row + 1) +=  fac * N[ii] * DN_DX(jj, 1);
            }
        }

        KRATOS_CATCH("")

    }

	//*************************************************************************************
	//*************************************************************************************
    void ASGSCompressible2D::CalculatePressureTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double, 3 > & N, const double time, const double area) {
        KRATOS_TRY
                int nodes_number = 3;
        int dof = 2;

        double density;
        double mu;
        calculatedensity(GetGeometry(), density, mu);

	double VC2;
	CalculateSoundVelocity(GetGeometry(), VC2);

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

        //************************************************************************************
    //************************************************************************************

    void ASGSCompressible2D::CalculateGradStblAllTerms(MatrixType& K, VectorType& F, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const double time, const double tauone, const double area) {
        KRATOS_TRY
                int nodes_number = 3;
        int dof = 2;


        double VC2;
	CalculateSoundVelocity(GetGeometry(), VC2);

        //build 1*(grad q . grad p) stabilization term & assemble
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > gard_opr = ZeroMatrix(nodes_number, nodes_number);
        gard_opr =  prod(DN_DX, trans(DN_DX));

        for (int ii = 0; ii < nodes_number; ii++) {
            int row = ii * (dof + 1) + dof;

            for (int jj = 0; jj < nodes_number; jj++) {
                int column = jj * (dof + 1) + dof;

                K(row, column) += area * tauone*VC2 * gard_opr(ii, jj);

            }
        }
        KRATOS_CATCH("")
    }
    	//*************************************************************************************
	//*************************************************************************************
	  void ASGSCompressible2D::GetFirstDerivativesVector(Vector& values, int Step)
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
	  void ASGSCompressible2D::GetSecondDerivativesVector(Vector& values, int Step)
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
	void ASGSCompressible2D::CalculateSoundVelocity(Geometry< Node<3> > geom, double& vc2)
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



	}
	//*************************************************************************************
	//*************************************************************************************
	void ASGSCompressible2D::CalculateCompressibleStblTerms(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,array_1d<double,3> N,const double tautwo,const double area)
	{

	  //double VC2;
	  //CalculateSoundVelocity(GetGeometry(), VC2);
	  double stbl_fac = tautwo * area;
	int nodes_number = 3;
	int dof = 2;
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		//LHS contribution
		int row = ii*(dof+1);
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;

			K(row,column) +=  stbl_fac* N(jj) * DN_DX(ii,0);
			K(row + 1,column) +=  stbl_fac* N(jj) * DN_DX(ii,1);

		   }
	    }

	}

	//************************************************************************************
	//************************************************************************************
	void ASGSCompressible2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
	  void ASGSCompressible2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
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

          void ASGSCompressible2D::CalculateTau(double& tauone, double& tautwo, const double time, const double area, const ProcessInfo& rCurrentProcessInfo)
          {

            KRATOS_TRY

            tauone = time;
            tautwo = time;


            KRATOS_CATCH("")

              }
    
	//*************************************************************************************
	//*************************************************************************************
	void ASGSCompressible2D::calculatedensity(Geometry< Node<3> > geom, double& density, double& viscosity)
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
	//************************************************************************************
	//************************************************************************************




} // Namespace Kratos


