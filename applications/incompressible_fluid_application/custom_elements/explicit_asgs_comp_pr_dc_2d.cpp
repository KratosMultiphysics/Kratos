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
#include "custom_elements/explicit_asgs_comp_pr_dc_2d.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{


	//************************************************************************************
	//************************************************************************************
	ExplicitASGSCOMPPRDC2D::ExplicitASGSCOMPPRDC2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: ExplicitASGSCompressible2D(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	ExplicitASGSCOMPPRDC2D::ExplicitASGSCOMPPRDC2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: ExplicitASGSCompressible2D(NewId, pGeometry, pProperties)
	{
		
	}

	Element::Pointer ExplicitASGSCOMPPRDC2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		
		KRATOS_TRY
		return Element::Pointer(new ExplicitASGSCOMPPRDC2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	ExplicitASGSCOMPPRDC2D::~ExplicitASGSCOMPPRDC2D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void ExplicitASGSCOMPPRDC2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
				rResult[i*node_size+2] = GetGeometry()[i].GetDof(WATER_PRESSURE).EquationId();
			}
		KRATOS_CATCH("")
		
	}

	//************************************************************************************
	//************************************************************************************
	  void ExplicitASGSCOMPPRDC2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
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
				ElementalDofList[i*node_size+2] = GetGeometry()[i].pGetDof(WATER_PRESSURE);
			}
		KRATOS_CATCH("");
	}


	//*************************************************************************************
	//************************************************************************************

	  void ExplicitASGSCOMPPRDC2D::GetFirstDerivativesVector(Vector& values, int Step)
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
			values[index + 2] = GetGeometry()[i].GetSolutionStepValue(WATER_PRESSURE,Step);

		}
	}
	//************************************************************************************
	//************************************************************************************
	  void ExplicitASGSCOMPPRDC2D::GetSecondDerivativesVector(Vector& values, int Step)
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
			values[index + 2] = GetGeometry()[i].GetSolutionStepValue(WATER_PRESSURE_DT,Step);
		}
	
	}
	//*************************************************************************************
	//*************************************************************************************
	void ExplicitASGSCOMPPRDC2D::calculatedensity(Geometry< Node<3> > geom, double& density, double& viscosity)
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
	
	const double rho0 = geom[0].FastGetSolutionStepValue(DENSITY_WATER);
	const double rho1 = geom[1].FastGetSolutionStepValue(DENSITY_WATER);
	const double rho2 = geom[2].FastGetSolutionStepValue(DENSITY_WATER);
	 density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );

	const double visc0 = geom[0].FastGetSolutionStepValue(VISCOSITY_WATER);
	const double visc1 = geom[1].FastGetSolutionStepValue(VISCOSITY_WATER);
	const double visc2 = geom[2].FastGetSolutionStepValue(VISCOSITY_WATER);
	 viscosity = 0.3333333333333333333333*(visc0 + visc1 + visc2 );
//density = 1000.0;
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
	//************************************************************************************
	//************************************************************************************
	void ExplicitASGSCOMPPRDC2D::CalculateSoundVelocity(Geometry< Node<3> > geom, double& vc2)
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
			mean_vc = GetGeometry()[ii].FastGetSolutionStepValue(WATER_SOUND_VELOCITY ) ;
			mean_vc2 += mean_vc * mean_vc;
		    }

	 vc2 =mean_vc2*0.333333333333333333333333;



	}
	//*************************************************************************************
	//*************************************************************************************

    void ExplicitASGSCOMPPRDC2D::CalculateResidual(const MatrixType& K, VectorType& F) 
	{
	    KRATOS_TRY

		    int nodes_number = 3;
	    int dof = 2;


	    array_1d<double, 9 > UP = ZeroVector(9);
	    for (int ii = 0; ii < nodes_number; ++ii) {
		int index = ii * (dof + 1);
		UP[index] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[0];
		UP[index + 1] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[1];
		UP[index + 2] = GetGeometry()[ii].FastGetSolutionStepValue(WATER_PRESSURE, 0);
	    }

	    F -= prod(K, UP);

	    KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
        void ExplicitASGSCOMPPRDC2D::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
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
            double VC = GetGeometry()[ii].FastGetSolutionStepValue(WATER_SOUND_VELOCITY ) ;

	    inv_max_h = sqrt(inv_max_h);
	    calc_t = 1.0/(inv_max_h * VC );
   
	    if( calc_t < Output)
		  Output = calc_t;
              
        }


        }
	//*************************************************************************************
	//*************************************************************************************
	void ExplicitASGSCOMPPRDC2D::CalculateDivPdotStblTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const array_1d<double,3>& N, const double time,const double tautwo,const double area)
	  {
	    KRATOS_TRY
	//tau*div(V).P_dot
	  double stbl_fac = tautwo * area;
          int nodes_number = 3;
          int dof = 2;

	  //N1=N2=N3=0.33333333333333333333333333333333333
         double mean_pressure_rate = GetGeometry()[0].FastGetSolutionStepValue(WATER_PRESSURE_DT);
	 for( int ii=1; ii < nodes_number; ++ii)
            mean_pressure_rate += GetGeometry()[ii].FastGetSolutionStepValue(WATER_PRESSURE_DT);

         mean_pressure_rate *= 0.33333333333333333333333333333333333333333333333;
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

} // Namespace Kratos


