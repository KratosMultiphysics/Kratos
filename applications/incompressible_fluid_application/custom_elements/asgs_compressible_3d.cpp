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
#include "custom_elements/asgs_compressible_3d.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
  /*      namespace ASGSCompressible2Dauxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
        #pragma omp threadprivate(DN_DX)

        array_1d<double,3> N = ZeroVector(3); //dimension = number of nodes
        #pragma omp threadprivate(N)

        array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
        #pragma omp threadprivate(ms_adv_vel)

    }
    using  namespace ASGSCompressible2Dauxiliaries;*/


	//************************************************************************************
	//************************************************************************************
	ASGSCompressible3D::ASGSCompressible3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: ASGS3D(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
		
	}

	//************************************************************************************
	//************************************************************************************
	ASGSCompressible3D::ASGSCompressible3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: ASGS3D(NewId, pGeometry, pProperties)
	{
			
	}

	Element::Pointer ASGSCompressible3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		
		KRATOS_TRY
		return Element::Pointer(new ASGSCompressible3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	ASGSCompressible3D::~ASGSCompressible3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void ASGSCompressible3D::CalculateLocalVelocityContribution(MatrixType& rDampMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
	int nodes_number = 4;
	int dim = 3;
	unsigned int matsize = nodes_number*(dim+1);

	if(rDampMatrix.size1() != matsize)
			rDampMatrix.resize(matsize,matsize,false); //false says not to preserve existing storage!!


	noalias(rDampMatrix) = ZeroMatrix(matsize,matsize); 

	double delta_t= rCurrentProcessInfo[DELTA_TIME];
	
	
	//getting data for the given geometry
	double Volume;
        boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX = ZeroMatrix(4,3);
        array_1d<double,4> N = ZeroVector(4);
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Volume);
	
	if(Volume <= 0.0)
	{
	   // std::cout << "element " << Id() << " has negative volume!!! " << Volume << std::endl;
	    noalias(rRightHandSideVector) = ZeroVector(matsize);
	}
	else
	{
	      //viscous term	
	      CalculateViscousTerm(rDampMatrix, DN_DX, Volume);
	      
	      //Advective term
	      double tauone;
	      double tautwo;
	      CalculateTau(tauone, tautwo, delta_t, Volume, rCurrentProcessInfo);

 	      CalculateAdvectiveTerm(rDampMatrix, DN_DX, tauone, tautwo, delta_t, Volume);


	      //calculate pressure term
	      CalculatePressureTerm(rDampMatrix, DN_DX, N, delta_t,Volume);

	      //compute projections
	      
	      //stabilization terms
	      //KRATOS_WATCH("ñññññññññ calculate stabilizing terms ñññññññññ");
	      CalculateDivStblTerm(rDampMatrix, DN_DX, tautwo, Volume);
			    
	      CalculateAdvStblAllTerms(rDampMatrix,rRightHandSideVector, DN_DX, N, tauone,delta_t, Volume);
		

	      CalculateGradStblAllTerms(rDampMatrix,rRightHandSideVector,DN_DX, delta_t, tauone, Volume);


	      //new stabilization added
	      //CalculateLCSMassContribution(rRightHandSideVector, delta_t,Volume);
		
	      //KRATOS_WATCH(rRightHandSideVector);

// double air_auxL = 0.0;
// double air_auxR = 0.0;
// for(unsigned int i=0; i<16; i++)
// {
//     air_auxR += rRightHandSideVector[i];
//     //if(CalculateStiffnessMatrixFlag==true)	
//       for(unsigned int j=0; j<16; j++)
// 	  air_auxL += rDampMatrix(i,j);
// }

// std::cout << Id() << " " << air_auxR << " " << air_auxL << std::endl;
// 					    double aaa = norm_2(row(DN_DX,0));
// 					    double bbb = norm_2(row(DN_DX,1));
// 					    double ccc = norm_2(row(DN_DX,2));
// 					    double ddd = norm_2(row(DN_DX,3));
// std::cout << "GRAD DN_DX "<< aaa << " " << bbb << " " << ccc << " " << ddd << std::endl;
// std::cout << "TAUs "<< tauone << " " << tautwo << std::endl;
// std::cout << "Volume: "<< Volume <<  std::endl;
// std::cout << "Volume: "<<  std::endl;
	  }

		KRATOS_CATCH("")
	}
        
	//************************************************************************************
	//************************************************************************************
	void ASGSCompressible3D::CalculateMassContribution(MatrixType& K,const double time,const double volume)
	{
		KRATOS_TRY
	double lump_mass_fac = volume * 0.25;
	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	//update density and calculate sound velocity
	double VC2;
	CalculateSoundVelocity(GetGeometry(), VC2);

	int nodes_number = 4;
	int dof = 3;
	for ( int nd = 0; nd< nodes_number; nd++)
	    {
		int row = nd*(dof + 1);
		for( int jj=0; jj< dof; jj++)
			K(row + jj, row + jj) += density/1.0*lump_mass_fac;

		//add pressure mass
			K(row + dof , row + dof ) += lump_mass_fac/VC2;
	    }
	
	
		KRATOS_CATCH("")
	
	}
	//************************************************************************************
	//************************************************************************************
	void ASGSCompressible3D::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
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
	double Volume;
        boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX = ZeroMatrix(4,3);
        array_1d<double,4> N = ZeroVector(4);
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Volume);

	//if(Volume <= 0.0)
	 //   std::cout << "element " << Id() << " has negative volume!!! " << Volume << std::endl;
	if(Volume > 0.0)
	{
	    //Calculate tau
	    double tauone;
	    double tautwo;
	    CalculateTau(tauone, tautwo, delta_t, Volume, rCurrentProcessInfo);

	    CalculateMassContribution(rMassMatrix,delta_t,Volume); 
		    //add stablilization terms due to advective term (a)grad(V) * ro*Acce
	    CalculateAdvMassStblTerms(rMassMatrix, DN_DX, N,tauone,Volume);
		    //add stablilization terms due to grad term grad(q) * ro*Acce
	    CalculateGradMassStblTerms(rMassMatrix, DN_DX, tauone,Volume);
		    //add compressible stabilization terms
	    CalculateCompressibleStblTerms(rMassMatrix, DN_DX,N, tautwo,Volume);
	}
	
	
		KRATOS_CATCH("")
	}

	//*************************************************************************************
	//*************************************************************************************

	  void ASGSCompressible3D::GetFirstDerivativesVector(Vector& values, int Step)
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
			values[index + 3] = GetGeometry()[i].GetSolutionStepValue(AIR_PRESSURE,Step);

		}
	}
	//************************************************************************************
	//************************************************************************************
	  void ASGSCompressible3D::GetSecondDerivativesVector(Vector& values, int Step)
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
			values[index + 3] = GetGeometry()[i].GetSolutionStepValue(AIR_PRESSURE_DT,Step);
		}
	
	}
	//************************************************************************************
	//************************************************************************************
	void ASGSCompressible3D::CalculateSoundVelocity(Geometry< Node<3> > geom, double& vc2)
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
			mean_vc = GetGeometry()[ii].FastGetSolutionStepValue(AIR_SOUND_VELOCITY ) ;
			mean_vc2 += mean_vc * mean_vc;
		    }

	 vc2 =mean_vc2*0.25;

   // KRATOS_WATCH("THIS IS An AIR ELEMENT");
             // KRATOS_WATCH(vc2);

	}
	//*************************************************************************************
	//*************************************************************************************
	void ASGSCompressible3D::CalculateCompressibleStblTerms(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,4,3>& DN_DX,array_1d<double,4> N,const double tautwo,const double volume)
	{

	  double VC2;
	  CalculateSoundVelocity(GetGeometry(), VC2);
	  double stbl_fac = tautwo/VC2 * volume;
	int nodes_number = 4;
	int dof = 3;
	for ( int ii = 0; ii < nodes_number; ii++)
	    {
		//LHS contribution
		int row = ii*(dof+1);
		for( int jj=0; jj < nodes_number; jj++)
		   {
			int column = jj*(dof+1) + dof;

			K(row,column) +=  stbl_fac* N(jj) * DN_DX(ii,0);
			K(row + 1,column) +=  stbl_fac* N(jj) * DN_DX(ii,1);
			K(row + 2,column) +=  stbl_fac* N(jj) * DN_DX(ii,2);

		   }
	    }

	}

	//************************************************************************************
	//************************************************************************************
	void ASGSCompressible3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
				rResult[i*node_size+3] = GetGeometry()[i].GetDof(AIR_PRESSURE).EquationId();
			}
		KRATOS_CATCH("")
			
	}

	//************************************************************************************
	//************************************************************************************
	  void ASGSCompressible3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
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
				ElementalDofList[i*node_size+3] = GetGeometry()[i].pGetDof(AIR_PRESSURE);
			}
		KRATOS_CATCH("");

	}
	//*************************************************************************************
	//*************************************************************************************
	void ASGSCompressible3D::calculatedensity(Geometry< Node<3> > geom, double& density, double& viscosity)
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
	const double rho3 = geom[3].FastGetSolutionStepValue(DENSITY_AIR);
	 density = 0.25*(rho0 + rho1 + rho2 + rho3);

	const double visc0 = geom[0].FastGetSolutionStepValue(VISCOSITY_AIR);
	const double visc1 = geom[1].FastGetSolutionStepValue(VISCOSITY_AIR);
	const double visc2 = geom[2].FastGetSolutionStepValue(VISCOSITY_AIR);
	const double visc3 = geom[3].FastGetSolutionStepValue(VISCOSITY_AIR);
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
	void ASGSCompressible3D::CalculateLCSMassContribution(VectorType& rhs,const double time,const double volume)
	{
		KRATOS_TRY
	double lump_mass_fac = volume * 0.25;
	double C_mass_fac = volume/20.0;
	double density;
	double mu;
	calculatedensity(GetGeometry(), density, mu);

	//update density and calculate sound velocity
	double VC2;
	CalculateSoundVelocity(GetGeometry(), VC2);

	int nodes_number = 4;
	int dof = 3;

        boost::numeric::ublas::bounded_matrix<double,4,4> CLSMass = ZeroMatrix(4,4);
	for ( int nd = 0; nd< nodes_number; nd++)
	    {


		CLSMass(nd , nd ) += lump_mass_fac/VC2;
		CLSMass(nd , nd) += -C_mass_fac/VC2;
// 		int row = nd*(dof + 1);
		for( int jj=0; jj< nodes_number; jj++)
			  CLSMass(nd , jj) += -C_mass_fac/VC2;			  
// 
// 		//add pressure mass
// 			CLSMass(row + dof , row + dof ) += lump_mass_fac/VC2;
	    }
	
	double alpha_fac = -0.1;
	Vector time_drv;
	GetSecondDerivativesVector(time_drv, 1);
        array_1d<double,4> pr_acc = ZeroVector(4);
	pr_acc(0) = alpha_fac*time_drv(3); 
	pr_acc(1) = alpha_fac*time_drv(7); 
	pr_acc(2) = alpha_fac*time_drv(11); 
	pr_acc(3) = alpha_fac*time_drv(15); 

        array_1d<double,4> rhs_ctr = ZeroVector(4);

	noalias(rhs_ctr) = prod(CLSMass,pr_acc);

        for (int ii = 0; ii < nodes_number; ++ii) {
            int index = ii * (dof + 1) + dof;
            rhs[index] +=  rhs_ctr(ii);
        }
	
	
		KRATOS_CATCH("")
	
	}

    //*************************************************************************************
    //*************************************************************************************

    void ASGSCompressible3D::CalculateTau(double& tauone, double& tautwo, const double time, const double volume, const ProcessInfo& rCurrentProcessInfo)
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
	double length2 = ele_length * ele_length;

        tauone = 1.0 / (dyn_st_beta / time + 4.0 * mu / (ele_length * ele_length * density) + 2.0 * advvel_norm / ele_length);
// std::cout << Id() <<" advvel_norm: " << advvel_norm << " " << "ele_length: " << ele_length << std::endl;
// std::cout << "mu density time " << mu << ""<< density << ""<< time << std::endl;

        tautwo = mu / density + 1.0 * ele_length * advvel_norm / 2.0;



        KRATOS_CATCH("")


    }


    //*************************************************************************************
    //*************************************************************************************


} // Namespace Kratos


