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
#include "custom_elements/explicit_asgs_comp_pr_dc_3d.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{


	//************************************************************************************
	//************************************************************************************
	ExplicitASGSCOMPPRDC3D::ExplicitASGSCOMPPRDC3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: ExplicitASGSCompressible3D(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	ExplicitASGSCOMPPRDC3D::ExplicitASGSCOMPPRDC3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: ExplicitASGSCompressible3D(NewId, pGeometry, pProperties)
	{
		
	}

	Element::Pointer ExplicitASGSCOMPPRDC3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		
		KRATOS_TRY
		return Element::Pointer(new ExplicitASGSCOMPPRDC3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	ExplicitASGSCOMPPRDC3D::~ExplicitASGSCOMPPRDC3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void ExplicitASGSCOMPPRDC3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 3;
		unsigned int node_size = dim+1;

		
			if(rResult.size() != number_of_nodes*node_size)
				rResult.resize(number_of_nodes*node_size,false);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				rResult[i*node_size] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
				rResult[i*node_size+1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
				rResult[i*node_size+2] = GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
				rResult[i*node_size+3] = GetGeometry()[i].GetDof(WATER_PRESSURE).EquationId();
			}
		KRATOS_CATCH("")
		
	}

	//************************************************************************************
	//************************************************************************************
	  void ExplicitASGSCOMPPRDC3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 3;
		unsigned int node_size = dim+1;

			if(ElementalDofList.size() != number_of_nodes*node_size)
				ElementalDofList.resize(number_of_nodes*node_size);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				ElementalDofList[i*node_size] = GetGeometry()[i].pGetDof(VELOCITY_X);
				ElementalDofList[i*node_size+1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
				ElementalDofList[i*node_size+2] = GetGeometry()[i].pGetDof(VELOCITY_Z);
				ElementalDofList[i*node_size+3] = GetGeometry()[i].pGetDof(WATER_PRESSURE);
			}
		KRATOS_CATCH("");
	}


	//*************************************************************************************
	//************************************************************************************

	  void ExplicitASGSCOMPPRDC3D::GetFirstDerivativesVector(Vector& values, int Step)
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
			values[index + 2] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Z,Step);
			values[index + 3] = GetGeometry()[i].GetSolutionStepValue(WATER_PRESSURE,Step);

		}
	}
	//************************************************************************************
	//************************************************************************************
	  void ExplicitASGSCOMPPRDC3D::GetSecondDerivativesVector(Vector& values, int Step)
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
			values[index + 2] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Z,Step);
			values[index + 3] = GetGeometry()[i].GetSolutionStepValue(WATER_PRESSURE_DT,Step);
		}
	
	}
	//*************************************************************************************
	//*************************************************************************************
	void ExplicitASGSCOMPPRDC3D::calculatedensity(Geometry< Node<3> > geom, double& density, double& viscosity)
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
	const double rho3 = geom[3].FastGetSolutionStepValue(DENSITY_WATER);
	 density = 0.25*(rho0 + rho1 + rho2 + rho3 );

	const double visc0 = geom[0].FastGetSolutionStepValue(VISCOSITY_WATER);
	const double visc1 = geom[1].FastGetSolutionStepValue(VISCOSITY_WATER);
	const double visc2 = geom[2].FastGetSolutionStepValue(VISCOSITY_WATER);
	const double visc3 = geom[3].FastGetSolutionStepValue(VISCOSITY_WATER);
	 viscosity = 0.25*(visc0 + visc1 + visc2 + visc3);


    
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
	void ExplicitASGSCOMPPRDC3D::CalculateSoundVelocity(Geometry< Node<3> > geom, double& vc2)
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
		for(int ii = 0; ii<4; ii++)
		    {
			mean_vc = GetGeometry()[ii].FastGetSolutionStepValue(WATER_SOUND_VELOCITY ) ;
			mean_vc2 += mean_vc * mean_vc;
		    }

	 vc2 =mean_vc2*0.25;

    // KRATOS_WATCH("THIS IS A WATER ELEMENT");
      //    KRATOS_WATCH(vc2);

	}
	//*************************************************************************************
	//*************************************************************************************


    void ExplicitASGSCOMPPRDC3D::CalculateResidual(const MatrixType& K, VectorType& F) {
        KRATOS_TRY

                int nodes_number = 4;
        int dof = 3;


        array_1d<double, 16 > UP = ZeroVector(16);
        for (int ii = 0; ii < nodes_number; ++ii) {
            int index = ii * (dof + 1);
            UP[index] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[0];
            UP[index + 1] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[1];
            UP[index + 2] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[2];
            UP[index + 3] = GetGeometry()[ii].FastGetSolutionStepValue(WATER_PRESSURE, 0);
        }

        F -= prod(K, UP);

        KRATOS_CATCH("")
    }

    //*************************************************************************************
    //*************************************************************************************
        void ExplicitASGSCOMPPRDC3D::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
       {

	   Output = 100.0;
           double volume;

	  boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX = ZeroMatrix(3, 2);
	  array_1d<double, 4 > N = ZeroVector(4); //dimension = number of nodes
	  GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, volume);

	double calc_t = 0.0;
        for (unsigned int ii = 0; ii < 4; ii++)
        {
            double inv_max_h = DN_DX(ii,0)*DN_DX(ii,0) + DN_DX(ii,1)*DN_DX(ii,1) + DN_DX(ii,2)*DN_DX(ii,2);
            double VC = GetGeometry()[ii].FastGetSolutionStepValue(WATER_SOUND_VELOCITY ) ;
	    inv_max_h = sqrt(inv_max_h);
	    calc_t = 1.0/(inv_max_h * VC );
	    
	    if( calc_t < Output)
		  Output = calc_t;
              
        }


        }

    //*************************************************************************************
    //*************************************************************************************
    void ExplicitASGSCOMPPRDC3D::CalculateTau(const array_1d<double,4>& N,double& tauone, double& tautwo, const double time, const double volume, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
                //calculate mean advective velocity and taus

//         const array_1d<double, 3 > & adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 0);
//         const array_1d<double, 3 > & mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
//         const array_1d<double, 3 > & adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY, 0);
//         const array_1d<double, 3 > & mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
//         const array_1d<double, 3 > & adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, 0);
//         const array_1d<double, 3 > & mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
//         const array_1d<double, 3 > & adv_vel3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY, 0);
//         const array_1d<double, 3 > & mesh_vel3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);
        array_1d<double, 3 > ms_adv_vel = ZeroVector(3); //dimesion coincides with space dimension

//         ms_adv_vel[0] = N[0]*(adv_vel0[0] - mesh_vel0[0]) + N[1]*(adv_vel1[0] - mesh_vel1[0]) + N[2]*(adv_vel2[0] - mesh_vel2[0]) + N[3]*(adv_vel3[0] - mesh_vel3[0]);
//         ms_adv_vel[1] = N[0]*(adv_vel0[1] - mesh_vel0[1]) + N[1]*(adv_vel1[1] - mesh_vel1[1]) + N[2]*(adv_vel2[1] - mesh_vel2[1]) + N[3]*(adv_vel3[1] - mesh_vel3[1]);
//         ms_adv_vel[2] = N[0]*(adv_vel0[2] - mesh_vel0[2]) + N[1]*(adv_vel1[2] - mesh_vel1[2]) + N[2]*(adv_vel2[2] - mesh_vel2[2]) + N[3]*(adv_vel3[2] - mesh_vel3[2]);

        ms_adv_vel[0] = 0.0;
        ms_adv_vel[1] = 0.0;
        ms_adv_vel[2] = 0.0;

        double advvel_norm = ms_adv_vel[0] * ms_adv_vel[0] + ms_adv_vel[1] * ms_adv_vel[1]+ ms_adv_vel[2] * ms_adv_vel[2];
        advvel_norm = sqrt(advvel_norm);

        double ele_length = pow(12*volume,0.333333333333333333333);  
        ele_length = 2.0/3.0 * ele_length * sqrt(3.00);

        double mu;
        //const double mu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
        //const double mu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
        //const double mu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
        //mu = 0.333333333333333333333333*(mu0 + mu1 + mu2);

        double density;
        calculatedensity(GetGeometry(), density, mu);


        const double dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];

        double VC2;
	CalculateSoundVelocity(GetGeometry(), VC2);
	double vc = sqrt(VC2);

	double int_time = ele_length/vc;
	double length2 = ele_length * ele_length;

       // tauone = 1.0 / (dyn_st_beta / time + 4.0 * mu / (length2 * density) + 2.0 * advvel_norm / ele_length);
        tauone = 1.0 / (dyn_st_beta / int_time + 4.0 * mu / (length2 * density) + 2.0 * advvel_norm / ele_length);

// std::cout << Id() <<" advvel_norm: " << advvel_norm << " " << "ele_length: " << ele_length << std::endl;
// std::cout << "mu density time " << mu << ""<< density << ""<< time << std::endl;

        //tautwo = mu / density + 1.0 * ele_length * advvel_norm / 2.0;
        //tautwo = int_time;
        tautwo = ele_length*vc;


// std::cout << "TAUs "<< tauone << " " << tautwo << std::endl;

        KRATOS_CATCH("")


    }
	//*************************************************************************************
	//*************************************************************************************
	void ExplicitASGSCOMPPRDC3D::CalculateDivPdotStblTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,4,3>& DN_DX,const array_1d<double,4>& N, const double time,const double tautwo,const double volume)
	  {
        KRATOS_TRY
	//tau*div(V).P_dot
      //double VC2;
      //CalculateSoundVelocity(GetGeometry(), VC2);
	  double stbl_fac = tautwo * volume;

          int nodes_number = 4;
          int dof = 3;

	  //N1=N2=N3=N4 = 0.25
         double mean_pressure_rate = 0.25 * (GetGeometry()[0].FastGetSolutionStepValue(WATER_PRESSURE_DT));
	 for( int ii=1; ii < nodes_number; ++ii)
            mean_pressure_rate +=  0.25 * (GetGeometry()[ii].FastGetSolutionStepValue(WATER_PRESSURE_DT));

	 mean_pressure_rate *= stbl_fac;

        for (int ii = 0; ii < nodes_number; ++ii) {
            int index = ii * (dof + 1);
            F[index] -=  DN_DX(ii,0) * mean_pressure_rate;
            F[index+1] -=  DN_DX(ii,1) * mean_pressure_rate;
            F[index+2] -=  DN_DX(ii,2) * mean_pressure_rate;
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
    //************************************************************************************
    //************************************************************************************
     void ExplicitASGSCOMPPRDC3D::CalculateArtifitialViscosity(double& Vel_art_visc,double& Pr_art_visc ,const boost::numeric::ublas::bounded_matrix<double,4,3>&DN_DX)
 	{
	    KRATOS_TRY    
	  
	 Vel_art_visc = 0.0; 
	 Pr_art_visc = 0.0;
         int nodes_number = 4;
	 double div_vel = 0.0;
         for( int ii=0; ii < nodes_number; ++ii)
	 {
	    const array_1d<double,3>& vel = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY);
	    div_vel += (vel[0]*DN_DX(ii,0) + vel[1]*DN_DX(ii,1) + vel[2]*DN_DX(ii,2));
	    
	 }  
	 
	 double H=0.0;
	 double norm_grad_p = 0.0;
	 
	double mu;
        double density;
        calculatedensity(GetGeometry(), density, mu);
	 
	 if( div_vel < 0.0)
	    {
             CalculateCharectristicLength(H,DN_DX,norm_grad_p);	
             Vel_art_visc = 1.4* abs(div_vel) * pow(H,2);

	     Pr_art_visc = 1.0*sqrt(norm_grad_p/density) * pow(H,1.5);
	    } 
	   
	   this->GetValue(VEL_ART_VISC)=Vel_art_visc;
	   this->GetValue(PR_ART_VISC)=Pr_art_visc;
	    
	    
	    KRATOS_CATCH("")
	}
    //************************************************************************************
    //************************************************************************************
     void ExplicitASGSCOMPPRDC3D::CalculateCharectristicLength(double& ch_length, const boost::numeric::ublas::bounded_matrix<double,4,3>& DN_DX,double& norm_grad )
 	{
	    KRATOS_TRY 
	    
  	  GeometryType::JacobiansType J;
	  GetGeometry().Jacobian(J); 

	  Matrix CC;
	  CC.resize(3,3,false);	  	
	  
	  Matrix inverse_CC;
	  inverse_CC.resize(3,3,false);	  	  

	  noalias(CC) = prod(J[0],trans(J[0]));	 
	  
	  double det=0.0;
	  
	 MathUtils<double>::InvertMatrix3(CC,inverse_CC ,det);  	  	  

	  int nodes_number = 4;
	  array_1d<double,3> mean_acc =  GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION);
	  double rho = GetGeometry()[0].FastGetSolutionStepValue(DENSITY_WATER);
	  double pr = GetGeometry()[0].FastGetSolutionStepValue(WATER_PRESSURE);
	  
	  array_1d<double,3> grad_rho,grad_pr;
	  grad_rho[0] = DN_DX(0,0)*rho;
	  grad_rho[1] = DN_DX(0,1)*rho;
	  grad_rho[2] = DN_DX(0,2)*rho;	  
	  
	  grad_pr[0] = 0.25*DN_DX(0,0)*pr;
	  grad_pr[1] = 0.25*DN_DX(0,1)*pr;
	  grad_pr[2] = 0.25*DN_DX(0,2)*pr;		  
		  
          for (int ii = 1; ii < nodes_number; ii++) {	  
	        mean_acc += GetGeometry()[ii].FastGetSolutionStepValue(ACCELERATION);  
		rho = GetGeometry()[ii].FastGetSolutionStepValue(DENSITY_WATER);
	        pr = GetGeometry()[ii].FastGetSolutionStepValue(WATER_PRESSURE);
		
		grad_rho[0] += DN_DX(ii,0)*rho;
		grad_rho[1] += DN_DX(ii,1)*rho;	
		grad_rho[2] += DN_DX(ii,2)*rho;	
		
		grad_pr[0] += 0.25*DN_DX(ii,0)*pr;
		grad_pr[1] += 0.25*DN_DX(ii,1)*pr;
		grad_pr[2] += 0.25*DN_DX(ii,2)*pr;			
		
	  }
	  mean_acc *= 0.25;
	  grad_rho *= 0.25;
 
	  
	  double norm_acc =  MathUtils<double>::Norm3(mean_acc);
	  double norm_grad_rho =  MathUtils<double>::Norm3(grad_rho);
	  
	  norm_grad =  MathUtils<double>::Norm3(grad_pr);

	  array_1d<double,3>  n_dir= ZeroVector(3);
	  if(norm_acc != 0.0)
	             n_dir = 0.75/norm_acc*mean_acc;
	  if(norm_grad_rho != 0.0)
	             n_dir +=  0.25/norm_grad_rho*grad_rho;
	  
	  double norm_n_dir =  MathUtils<double>::Norm3(n_dir);	
 	  
	  if(norm_n_dir != 0.0)
	    {
	      n_dir/=norm_n_dir;	  
	      array_1d<double,3>  CC_n;	  
	      CC_n = prod(inverse_CC,n_dir);	  
	      double denom = n_dir[0]*CC_n[0] + n_dir[1]*CC_n[1] + n_dir[2]*CC_n[2];

	      if(denom <= 0.0)
		    KRATOS_ERROR(std::logic_error,"CalculateCharectristicLength zero or negative denominator ",denom)	
	      else
		    ch_length = 2.0/sqrt(denom);
	    }
	   else
	         ch_length = 0.0;
	  
	    KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************

    void ExplicitASGSCOMPPRDC3D::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
            array_1d<double, 3 > & Output,
            const ProcessInfo& rCurrentProcessInfo) 
        {
        KRATOS_TRY

	  Output = ZeroVector(3);

	  //getting data for the given geometry
// 	  double volume = GetGeometry().Volume();
        //getting data for the given geometry
        double volume;
        boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX = ZeroMatrix(4, 3);
        array_1d<double, 4 > N = ZeroVector(4); //dimension = number of nodes

        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, volume);

	  double lump_mass_fac = volume * 0.25;
	  double density;
	  double mu;
	  calculatedensity(GetGeometry(), density, mu);
	//update density and calculate sound velocity
	double VC2;
	CalculateSoundVelocity(GetGeometry(), VC2);

	  //fill velocity and pressure mass
	  Output[0] = density*lump_mass_fac;
	  Output[1] = lump_mass_fac;
 	  //Output[1] = lump_mass_fac/VC2;

        KRATOS_CATCH("")
       }
	//************************************************************************************
	//************************************************************************************


} // Namespace Kratos


