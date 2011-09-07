/*
* File:   edgebased_levelset.h
* Author: rrossi
*
* Created on July 31, 2009, 10:51 AM
*/

/*
==============================================================================
KratosPFEMApplication
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
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2009-01-14 16:24:38 $
//   Revision:            $Revision: 1.11 $
//
//


#if !defined(KRATOS_EDGEBASED_FLUID_SOLVER_H_INCLUDED)
#define  KRATOS_EDGEBASED_FLUID_SOLVER_H_INCLUDED

//#define SPLIT_OSS
//#define SYMM_PRESS


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// #include <omp.h>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
//#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
#include "incompressible_fluid_application.h"


namespace Kratos {

    template<unsigned int TDim, class MatrixContainer, class TSparseSpace, class TLinearSolver>
    class FluidSolver {
    public:
	//name for the self defined structure
	typedef EdgesStructureType<TDim> CSR_Tuple;
	typedef vector<CSR_Tuple> EdgesVectorType;

	//name for row start and column index vectors
	typedef vector<unsigned int> IndicesVectorType;
	//defining matrix type for test calculations
	typedef vector< array_1d<double, TDim> > CalcVectorType;
	//defining type for local storage of nodal values
	typedef vector<double> ValuesVectorType;

	//defining types for matrix operations
	typedef typename TSparseSpace::MatrixType TSystemMatrixType;
	typedef typename TSparseSpace::VectorType TSystemVectorType;

        typedef std::size_t SizeType;

	//constructor and destructor

	FluidSolver(MatrixContainer& mr_matrix_container,
		ModelPart& mr_model_part,
		const double viscosity,
		const double density,
		const Vector body_force,
		bool use_mass_correction,
		double edge_detection_angle,
		double stabdt_pressure_factor,
		double stabdt_convection_factor,
		double tau2_factor,
		bool assume_constant_dp
		)
	: mr_matrix_container(mr_matrix_container),
	  mr_model_part(mr_model_part),
	  mstabdt_pressure_factor(stabdt_pressure_factor),
	    mstabdt_convection_factor(stabdt_convection_factor),
	    medge_detection_angle(edge_detection_angle),
	    mtau2_factor(tau2_factor),
	    massume_constant_dp(assume_constant_dp)

	{
	    mViscosity = viscosity;

	    noalias(mBodyForce) = body_force;
	    mRho = density;

	    mdelta_t_avg = 1000.0;

	    max_dt = 1.0;

	    muse_mass_correction = use_mass_correction;
		
	    mWallLawIsActive = false;

//            for (unsigned int i = 0; i < TDim; i++) mBodyForce[i] = 0;
//            mBodyForce[1] = -9.81;
//
//            mRho = 1000.0;




	};

	~FluidSolver() {
	};

	//***********************************
	//function to initialize fluid solver

	void Initialize(
		) {
	    KRATOS_TRY


	  //get number of nodes
	    unsigned int n_nodes = mr_model_part.Nodes().size();
	    unsigned int n_edges = mr_matrix_container.GetNumberEdges();
	    //size data vectors
	    mWork.resize(n_nodes);
	    mvel_n.resize(n_nodes);
	    mvel_n1.resize(n_nodes);
	    mPn.resize(n_nodes);
	    mPn1.resize(n_nodes);
	    mHmin.resize(n_nodes);
	    mHavg.resize(n_nodes);
	    mNodalFlag.resize(n_nodes);

	    mTauPressure.resize(n_nodes);
	    mTauConvection.resize(n_nodes);
	    mTau2.resize(n_nodes);
	    mPi.resize(n_nodes);
	    mXi.resize(n_nodes);
	    mx.resize(n_nodes);

	    mEdgeDimensions.resize(n_edges);

	    //convection variables
	    mBeta.resize(n_nodes);

	    mdiv_error.resize(n_nodes);
	    mr_matrix_container.SetToZero(mdiv_error);


//	    ValuesVectorType external_pressure;
//	    external_pressure.resize(n_nodes);

	    //read velocity and pressure data from Kratos
	    mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes());
	    mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, mr_model_part.Nodes());
	    mr_matrix_container.FillOldScalarFromDatabase(PRESSURE, mPn, mr_model_part.Nodes());
	    mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mvel_n, mr_model_part.Nodes());
	    mr_matrix_container.FillCoordinatesFromDatabase(mx, mr_model_part.Nodes());
	    //set flag for first time step
	    mFirstStep = true;

    
	    //loop to categorize boundary nodes
	    std::vector< unsigned int> tempFixedVelocities;
	    std::vector< array_1d<double,TDim> > tempFixedVelocitiesValues;
	    std::vector< unsigned int> tempPressureOutletList;
	    for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
		    inode != mr_model_part.NodesEnd();
		    inode++) {
		int index = inode->FastGetSolutionStepValue(AUX_INDEX);
		if (inode->IsFixed(VELOCITY_X)) //note that the variables can be either all fixed or no one fixed
		{
		    if (inode->IsFixed(VELOCITY_Y) == false || inode->IsFixed(VELOCITY_Z) == false) {
			std::cout << "error found on the fixity of node " << inode->Id() << std::endl;
			KRATOS_ERROR(std::logic_error, "velocities can be either all fixed or none fixed", "")
		    }


		    tempFixedVelocities.push_back(index);
		    tempFixedVelocitiesValues.push_back(mvel_n1[index]);
		}

		if (inode->IsFixed(PRESSURE)) {
		    tempPressureOutletList.push_back(index);
//		    mPressureOutlet.push_back(external_pressure[index]);
		}
	    }
	    mFixedVelocities.resize(tempFixedVelocities.size(),false);
	    mFixedVelocitiesValues.resize(tempFixedVelocitiesValues.size(),false);
	    mPressureOutletList.resize(tempPressureOutletList.size(),false);
	    
	    #pragma omp parallel for
	    for(unsigned int i=0; i<tempFixedVelocities.size(); i++)
	    {
	      mFixedVelocities[i] = tempFixedVelocities[i];
	      mFixedVelocitiesValues[i] = tempFixedVelocitiesValues[i];
	    }
	    #pragma omp parallel for
	    for(unsigned int i=0; i<tempPressureOutletList.size(); i++)
	    {
	      mPressureOutletList[i] = tempPressureOutletList[i];
	    }	

	    //compute slip normals and fill SlipList
	    CalculateNormals(mr_model_part.Conditions());
	    mr_matrix_container.WriteVectorToDatabase(NORMAL, mSlipNormal, mr_model_part.Nodes());

	    if(TDim == 3)
		DetectEdges3D(mr_model_part.Conditions());


	    //determine number of edges and entries
	    unsigned int n_nonzero_entries = 2 * n_edges + n_nodes;
	    //allocate memory for variables
	    mL.resize(n_nodes, n_nodes, n_nonzero_entries);

	    //loop over all nodes
	    for (unsigned int i_node = 0; i_node < n_nodes; i_node++) {
		//flag for considering diagonal matrix elements
		bool flag = 0;

		//loop over all neighbours
		for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
		    //get global index of neighbouring node j
		    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
		    //define matrix structure row by row (the order does matter!)
		    if ((j_neighbour > i_node) && (flag == 0)) {
			//add diagonal/nodal contribution
			mL.push_back(i_node, i_node, 0.0);
			flag = 1;
		    }
		    //add non-diagonal/edge contribution
		    mL.push_back(i_node, j_neighbour, 0.0);
		}
		//if diagonal element is the last non-zero element of the row
		if (flag == 0)
		    mL.push_back(i_node, i_node, 0.0);
	    }



	    //compute minimum length of the surrounding edges
	    CalculateEdgeLengths(mr_model_part.Nodes());


		      //set the pressure projection to the body force value
		      array_1d<double,3> temp = mRho * mBodyForce;
		      for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
			      inode != mr_model_part.NodesEnd();
			      inode++)
			  inode->FastGetSolutionStepValue(PRESS_PROJ) = temp;

	    KRATOS_CATCH("")
	}
	
	


	//***************************************
	//function to set adequate time step size

	double ComputeTimeStep(const double CFLNumber, const double MaxDt)
	{
	    KRATOS_TRY

	    //save the maximum time step
	    max_dt = MaxDt;

	    //local variable for time step size 
	    double delta_t = 1e10;

	    mdelta_t_avg = 1e10;

	    //getting value of current velocity and of viscosity
	    mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes());

	    //*******************
	    //loop over all nodes
	    double n_nodes = mvel_n1.size();
	    for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
	    {
		const array_1d<double, TDim>& v_i = mvel_n1[i_node];
		const double havg_i = mHavg[i_node];
		const double hmin_i = mHmin[i_node];

		double vel_norm = norm_2(v_i);


		//use CFL condition to compute time step size
		double delta_t_i = CFLNumber * 1.0 / (2.0 * vel_norm /hmin_i + 4.0 * mViscosity / (hmin_i * hmin_i) );
		double delta_t_i_avg = 1.0 / (2.0 * vel_norm /havg_i + 4.0 * mViscosity / (havg_i * havg_i) );

		//considering the most restrictive case of neighbor's velocities with similar direction but opposite sense.
		//loop over all neighbours
		for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
		    //get global index of neighbouring node j
		    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];

		    const array_1d<double, TDim>& v_j = mvel_n1[j_neighbour];

		    double v_diff_norm = 0.0;
		    for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
		    {
			double temp = v_i[l_comp] - v_j[l_comp];
			v_diff_norm += temp*temp;
		    }
		    v_diff_norm = sqrt(v_diff_norm);

		    double delta_t_j = CFLNumber * 1.0 / (2.0 * v_diff_norm /hmin_i + 4.0 * mViscosity / (hmin_i * hmin_i));

		    if (delta_t_j < delta_t_i)
			    delta_t_i = delta_t_j;
		}

		//choose the overall minimum of delta_t_i
		if (delta_t_i < delta_t)
		    delta_t = delta_t_i;

		if(delta_t_i_avg < mdelta_t_avg)
		  mdelta_t_avg = delta_t_i_avg;

	    }
	    //*******************
	    //perform MPI syncronization of the dt (minimum should be kept)

	    return delta_t;

	    KRATOS_CATCH("")
	}

	void UpdateFixedVelocityValues()
	{
	    KRATOS_TRY

	    //read velocity and pressure data from Kratos
	    ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
	    mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes);

	    int fixed_size = mFixedVelocities.size();
	    #pragma omp parallel for firstprivate(fixed_size)
	    for (int i_velocity = 0; i_velocity < fixed_size; i_velocity++)
	    {
		unsigned int i_node = mFixedVelocities[i_velocity];
		array_1d<double, TDim>& u_i_fix = mFixedVelocitiesValues[i_velocity];
		const array_1d<double, TDim>& u_i = mvel_n1[i_node];

		for (unsigned int comp = 0; comp < TDim; comp++)
			u_i_fix[comp] = u_i[comp];
	    }
	    KRATOS_CATCH("");
	}

	//**********************************************************************************
	//function to solve fluid equations - fractional step 1: compute fractional momentum

      void SolveStep1() {
	    KRATOS_TRY

	    //PREREQUISITES

	    //variables for node based data handling
	    ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
	    int n_nodes = rNodes.size();
	    //storage of nodal values in local variables
	    CalcVectorType rhs;
	    rhs.resize(n_nodes);


	    //read velocity and pressure data from Kratos
	    mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes);
	    mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mvel_n, rNodes);

	    mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, rNodes);
	    mr_matrix_container.FillOldScalarFromDatabase(PRESSURE, mPn, rNodes);

	    //read time step size from Kratos
	    ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
	    double delta_t = CurrentProcessInfo[DELTA_TIME];


	    //compute intrinsic time
//            double time_inv = 1.0 / delta_t;
	    double time_inv_avg = 1.0/mdelta_t_avg;

	    double stabdt_pressure_factor  = mstabdt_pressure_factor;
	    double stabdt_convection_factor  = mstabdt_convection_factor;
	    double tau2_factor = mtau2_factor;
            KRATOS_WATCH(stabdt_pressure_factor);
	    #pragma omp parallel for firstprivate(time_inv_avg,stabdt_pressure_factor,stabdt_convection_factor,tau2_factor)
	    for (int i_node = 0; i_node < n_nodes; i_node++) {
		double& h_avg_i = mHavg[i_node];
		array_1d<double, TDim>& a_i = mvel_n1[i_node];
		const double nu_i = mViscosity;

		double vel_norm = norm_2(a_i);

		double tau = 1.0 / (2.0 * vel_norm / h_avg_i + stabdt_pressure_factor*time_inv_avg + (4.0*nu_i) / (h_avg_i * h_avg_i) );
		double tau_conv = 1.0 / (2.0 * vel_norm / h_avg_i + stabdt_convection_factor*time_inv_avg + (4.0*nu_i) / (h_avg_i * h_avg_i) );
		mTauPressure[i_node] = tau;
		mTauConvection[i_node] = tau_conv;

		mTau2[i_node] = (mViscosity + h_avg_i*vel_norm*0.5)*tau2_factor;
	    }

	    //calculating the convective projection
	    #pragma omp parallel for
	    for (int i_node = 0; i_node < n_nodes; i_node++) {
		array_1d<double, TDim>& pi_i = mPi[i_node]; 

		//setting to zero
		for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
		    pi_i[l_comp] = 0.0;

		array_1d<double, TDim> a_i = mvel_n1[i_node];
		const array_1d<double, TDim>& U_i = mvel_n1[i_node];

		for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
		    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
		    array_1d<double, TDim> a_j = mvel_n1[j_neighbour];
		    const array_1d<double, TDim>& U_j = mvel_n1[j_neighbour];

		    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

		    edge_ij.Add_ConvectiveContribution(pi_i, a_i, U_i, a_j, U_j);
		}

		const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];

		for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
		    pi_i[l_comp] *= m_inv;
	    }




	    mr_matrix_container.AssignVectorToVector(mvel_n, mWork); //mWork = mvel_n

	    //first step of Runge Kutta
	    mr_matrix_container.AssignVectorToVector(mvel_n, mvel_n1); //mvel_n1 = mvel_n

	    mr_matrix_container.SetToZero(rhs);
	    CalculateRHS(mvel_n1, mPn, mvel_n1, rhs);
	    mr_matrix_container.Add_Minv_value(mWork, mWork, delta_t / 6.0, mr_matrix_container.GetInvertedMass(), rhs);
	    mr_matrix_container.Add_Minv_value(mvel_n1, mvel_n, 0.5 * delta_t, mr_matrix_container.GetInvertedMass(), rhs);
	    ApplyVelocityBC(mvel_n1);

//            KRATOS_WATCH("end of first stage")
//            double vnorm2 = 0.0;
//                                for( int i = 0; i< rNodes.size(); i++)
//                                    vnorm2 += pow(mvel_n1[i][0],2) + pow(mvel_n1[i][1],2) + pow(mvel_n1[i][2],2);
//                                KRATOS_WATCH(sqrt(vnorm2));

	    //second step
	    mr_matrix_container.SetToZero(rhs);
	    CalculateRHS(mvel_n1, mPn, mvel_n1, rhs);
	    mr_matrix_container.Add_Minv_value(mWork, mWork, delta_t / 3.0, mr_matrix_container.GetInvertedMass(), rhs);
	    mr_matrix_container.Add_Minv_value(mvel_n1, mvel_n, 0.5 * delta_t, mr_matrix_container.GetInvertedMass(), rhs);
	    ApplyVelocityBC(mvel_n1);

//                        KRATOS_WATCH("end of second stage")
//            vnorm2 = 0.0;
//                                for( int i = 0; i< rNodes.size(); i++)
//                                    vnorm2 += pow(mvel_n1[i][0],2) + pow(mvel_n1[i][1],2) + pow(mvel_n1[i][2],2);
//                                KRATOS_WATCH(sqrt(vnorm2));

	    //third step
	    mr_matrix_container.SetToZero(rhs);
	    CalculateRHS(mvel_n1, mPn, mvel_n1, rhs);
	    mr_matrix_container.Add_Minv_value(mWork, mWork, delta_t / 3.0, mr_matrix_container.GetInvertedMass(), rhs);
	    mr_matrix_container.Add_Minv_value(mvel_n1, mvel_n, delta_t, mr_matrix_container.GetInvertedMass(), rhs);
	    ApplyVelocityBC(mvel_n1);

//                        KRATOS_WATCH("end of thir stage")
//            vnorm2 = 0.0;
//                                for( int i = 0; i< rNodes.size(); i++)
//                                    vnorm2 += pow(mvel_n1[i][0],2) + pow(mvel_n1[i][1],2) + pow(mvel_n1[i][2],2);
//                                KRATOS_WATCH(sqrt(vnorm2));

	    //fourth step
	    mr_matrix_container.SetToZero(rhs);
	    CalculateRHS(mvel_n1, mPn, mvel_n1, rhs);
	    mr_matrix_container.Add_Minv_value(mWork, mWork, delta_t / 6.0, mr_matrix_container.GetInvertedMass(), rhs);

	    //compute right-hand side
	    mr_matrix_container.AssignVectorToVector(mWork, mvel_n1);
	    ApplyVelocityBC(mvel_n1);

//            KRATOS_WATCH("end of Step1")
//            vnorm2 = 0.0;
//                                for( int i = 0; i< rNodes.size(); i++)
//                                    vnorm2 += pow(mvel_n1[i][0],2) + pow(mvel_n1[i][1],2) + pow(mvel_n1[i][2],2);
//                                KRATOS_WATCH(sqrt(vnorm2));

	    KRATOS_CATCH("")
	}



	//*********************************************************************
	//function to calculate right-hand side of fractional momentum equation

	void CalculateRHS(
		const CalcVectorType& vel,
		const ValuesVectorType& pressure,
		const CalcVectorType& convective_velocity,
		CalcVectorType& rhs) 
	{
	    KRATOS_TRY

	    int n_nodes = vel.size();

	    //calculating the RHS
	    array_1d<double, TDim> stab_low;
	    array_1d<double, TDim> stab_high;
	    const double nu_i = mViscosity;
	    const double nu_j = mViscosity;
	    double inverse_rho = 1.0 / mRho;
	    #pragma omp parallel for private(stab_low,stab_high)
	    for (int i_node = 0; i_node < n_nodes; i_node++) {
		    array_1d<double, TDim>& rhs_i = rhs[i_node];
		    const array_1d<double, TDim>& f_i = mBodyForce;
		    array_1d<double, TDim> a_i = convective_velocity[i_node];

		    const array_1d<double, TDim>& U_i = vel[i_node];
		    const array_1d<double, TDim>& pi_i = mPi[i_node];
		    const double& p_i = pressure[i_node];

		    double edge_tau = mTauConvection[i_node];

		    //initializing with the external forces (e.g. gravity)
		    double& m_i = mr_matrix_container.GetLumpedMass()[i_node];
		    for (unsigned int comp = 0; comp < TDim; comp++)
			rhs_i[comp] = m_i  * f_i[comp] ;

		    //convective term
		    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
			unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
			    array_1d<double, TDim> a_j = convective_velocity[j_neighbour];
			    const array_1d<double, TDim>& U_j = vel[j_neighbour];
			    const array_1d<double, TDim>& pi_j = mPi[j_neighbour];
			    const double& p_j = pressure[j_neighbour];

			    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

			    edge_ij.Sub_ConvectiveContribution(rhs_i, a_i, U_i, a_j, U_j);
			    edge_ij.Sub_grad_p(rhs_i, p_i*inverse_rho, p_j * inverse_rho);
			    edge_ij.Sub_ViscousContribution(rhs_i, U_i, nu_i, U_j, nu_j);

			    //add stabilization
			    edge_ij.CalculateConvectionStabilization_LOW(stab_low, a_i, U_i, a_j, U_j);
			    edge_ij.CalculateConvectionStabilization_HIGH(stab_high, a_i, pi_i, a_j, pi_j);

//                            double beta = 1.0;
//                             double beta = beta_i;
//                             if(beta_j > beta)
//                                 beta = beta_j;
//                            beta = 1.0; 

//                            edge_ij.Sub_StabContribution(rhs_i, edge_tau*beta, 1.0, stab_low, stab_high);
//                             edge_ij.Sub_StabContribution(rhs_i, edge_tau, (1.0-beta), stab_low, stab_high);
			    edge_ij.Sub_StabContribution(rhs_i, edge_tau, 1.0, stab_low, stab_high);
			    

			    //add tau2 term
//                             boost::numeric::ublas::bounded_matrix<double,TDim,TDim>& LL = edge_ij.LaplacianIJ;
//                             for (unsigned int k_comp = 0; k_comp < TDim; k_comp++)
//                             {
//                                 double aaa = 0.0;
//                                 for (unsigned int m_comp = 0; m_comp < TDim; m_comp++)
//                                     aaa +=  LL(k_comp,m_comp) * (U_j[m_comp] - U_i[m_comp]);
//                                 rhs_i[k_comp] -= tau2_i*aaa;
//                             }


		    }

		

	    }

	    //apply wall resistance
	  if(mWallLawIsActive == true)
		  ComputeWallResistance(vel,rhs);

	    ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
	    mr_matrix_container.WriteVectorToDatabase(VELOCITY, mvel_n1, rNodes);
	    KRATOS_CATCH("")
	}

	//*************************************************************************
	//function to solve fluid equations - fractional step 2: calculate pressure

	void SolveStep2(typename TLinearSolver::Pointer pLinearSolver) {
	    KRATOS_TRY
	    
	    //PREREQUISITES

	    //allocate memory for variables
	    ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
	    int n_nodes = rNodes.size();
	    //unknown and right-hand side vector
	    TSystemVectorType dp, rhs;
	    dp.resize(n_nodes);
	    rhs.resize(n_nodes);
	    array_1d<double, TDim> dU_i, dU_j, work_array;
	    //read time step size from Kratos
	    ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
	    double delta_t = CurrentProcessInfo[DELTA_TIME];

#ifdef _OPENMP
	    double time_inv = 0.0; //1.0/delta_t;

	    //read the pressure projection from the database
#endif
	    mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, mr_model_part.Nodes());
            mr_matrix_container.FillOldScalarFromDatabase(PRESSURE, mPn, mr_model_part.Nodes());
	    mr_matrix_container.FillVectorFromDatabase(PRESS_PROJ, mXi, rNodes);
	    mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes);
	    mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mvel_n, rNodes);

	    #pragma omp parallel for firstprivate(time_inv)
	    for (int i_node = 0; i_node < n_nodes; i_node++) {

		double& rhs_i = rhs[i_node];
		rhs_i = 0.0;
		const double& p_i = mPn1[i_node];
		const double& p_old_i = mPn[i_node];
		const array_1d<double, TDim>& U_i_curr = mvel_n1[i_node];

		array_1d<double, TDim>& xi_i = mXi[i_node];

		double l_ii = 0.0;

		//loop over all neighbours
		for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
		    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
		    const double& p_j = mPn1[j_neighbour];
		    const double& p_old_j = mPn[j_neighbour];
		    const array_1d<double, TDim>& U_j_curr = mvel_n1[j_neighbour];
		    const array_1d<double, TDim>& xi_j = mXi[j_neighbour];

		    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

#ifdef SYMM_PRESS
		    double edge_tau = 0.5 * (mTauPressure[i_node] + mTauPressure[j_neighbour]);
#else
		    double edge_tau = mTauPressure[i_node];
#endif



		    //compute laplacian operator
		    double sum_l_ikjk;
		    edge_ij.CalculateScalarLaplacian(sum_l_ikjk);
//                    sum_l_ikjk *= 2.0;
                    
		    double sum_l_ikjk_onlydt = sum_l_ikjk * (2.0*delta_t);
		    sum_l_ikjk *= (2.0*delta_t + edge_tau);

		    //assemble right-hand side
		    //pressure contribution
		    rhs_i -= sum_l_ikjk * (p_j - p_i);
		    rhs_i += sum_l_ikjk_onlydt * (p_old_j - p_old_i);


		    //calculating the divergence of the fract vel
		    edge_ij.Sub_D_v(rhs_i, U_i_curr*mRho, U_j_curr * mRho);

		    //high order stabilizing term
		    double temp = 0.0;
		    edge_ij.Add_div_v(temp, xi_i, xi_j);
		    rhs_i += edge_tau * temp;

		    //assemble laplacian matrix
		    mL(i_node, j_neighbour) = sum_l_ikjk;
		    l_ii -= sum_l_ikjk;
		}

		mL(i_node, i_node) = l_ii;
	    }

	    if(muse_mass_correction == true)
	    {
                std::cout << "****************************************" << std::endl;
		#pragma omp parallel for
		for (int i_node = 0; i_node < n_nodes; i_node++)
		{
		    double& rhs_i = rhs[i_node];
		    rhs_i -= mdiv_error[i_node];
		}
	    }

// 	    //find the max diagonal term
// 	    double max_diag = 0.0;
// 	    for (int i_node = 0; i_node < n_nodes; i_node++) {
// 		double L_diag = mL(i_node, i_node);
// 		if (fabs(L_diag) > fabs(max_diag)) max_diag = L_diag;
// 	    }




	    //respect pressure boundary conditions by penalization
            double huge = 1e20;
           for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++) {
               unsigned int i_node = mPressureOutletList[i_pressure];
               mL(i_node, i_node) = huge;
               rhs[i_node] = 0.0;
           }
// 	    for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++) {
// 		unsigned int i_node = mPressureOutletList[i_pressure];
// 		mL(i_node, i_node) = max_diag;
// 		rhs[i_node] = 0.0;
// 		for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
// 		{
// 		    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
// 		    mL(i_node, j_neighbour) = 0.0;
// 		}
// 	    }

	    //set starting vector for iterative solvers
	    for (int i_node = 0; i_node < n_nodes; i_node++)
		dp[i_node] = 0.0;

            //compute row scaling factors
            TSystemVectorType scaling_factors(n_nodes);
            double* Lvalues = mL.value_data().begin();
            SizeType* Lrow_indices = mL.index1_data().begin();
            SizeType* Lcol_indices = mL.index2_data().begin();

            for (SizeType k = 0; k < mL.size1(); k++)
            {
                double t = 0.0;
                SizeType col_begin = Lrow_indices[k];
                SizeType col_end = Lrow_indices[k+1];

                for (SizeType j=col_begin; j<col_end; j++)
                    if( Lcol_indices[j] == k)
                    {
                        t = fabs(Lvalues[j]);
                    }
//                        t += Lvalues[j]*Lvalues[j];

//                t = sqrt(t);
                scaling_factors[k] = 1.0/sqrt(t);
            }

            for (SizeType k = 0; k < mL.size1(); k++)
            {
                SizeType col_begin = Lrow_indices[k];
                SizeType col_end = Lrow_indices[k+1];
                double k_factor = scaling_factors[k];

                rhs[k] *= k_factor;

                for (SizeType j=col_begin; j<col_end; j++)
                {
                    Lvalues[j] *= scaling_factors[Lcol_indices[j]] * k_factor;
                }
            }

//            double huge = 1e20;
//           for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++) {
//               unsigned int i_node = mPressureOutletList[i_pressure];
//               mL(i_node, i_node) = 1.0;
//               rhs[i_node] = 0.0;
//           }

//            KRATOS_WATCH(norm_2(rhs));
//            KRATOS_WATCH(norm_frobenius(mL));
	    pLinearSolver->Solve(mL, dp, rhs);

            //apply inverse scaling
            for (unsigned int k = 0; k < dp.size(); k++)
                dp[k] *= scaling_factors[k];

	    KRATOS_WATCH(*pLinearSolver)
//KRATOS_WATCH(norm_2(dp));

	    //update pressure
	    for (int i_node = 0; i_node < n_nodes; i_node++)
		mPn1[i_node] += dp[i_node];

	    //write pressure and density to Kratos
	    mr_matrix_container.WriteScalarToDatabase(PRESSURE, mPn1, rNodes);


	    //compute pressure proj for the next step

	    #pragma omp parallel for firstprivate(time_inv), private(work_array)
	    for (int i_node = 0; i_node < n_nodes; i_node++) {
		array_1d<double, TDim>& xi_i = mXi[i_node];
		for (unsigned int comp = 0; comp < TDim; comp++)
		    xi_i[comp] = 0.0;


		    const double& p_i = mPn1[i_node];

		    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
			//get global index of neighbouring node j
			unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];

			    const double& p_j = mPn1[j_neighbour];

			    //projection of pressure gradients
			    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

			    edge_ij.Add_grad_p(xi_i, p_i, p_j);
		    }

		    const double& m_inv = mr_matrix_container.GetInvertedMass()[i_node];
		    for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
			xi_i[l_comp] *= m_inv;
		

	    }

	    mr_matrix_container.WriteVectorToDatabase(PRESS_PROJ, mXi, rNodes);

	    KRATOS_CATCH("")
	}

	//**********************************************************************************
	//function to solve fluid equations - fractional step 3: correct fractional momentum

	void SolveStep3() {
	    KRATOS_TRY
	    //get number of nodes
	    ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();

	    int n_nodes = rNodes.size();

	    //define work array
	    array_1d<double, TDim> correction;
	    //read time step size from Kratos
	    ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
	    double delta_t = CurrentProcessInfo[DELTA_TIME];

	    double factor = 0.5;
	    if(massume_constant_dp == true)
		factor = 1.0;

	    //compute end of step momentum
	    double rho_inv = 1.0 / mRho;
	    #pragma omp parallel for private(correction) firstprivate(delta_t,rho_inv,factor)
	    for (int i_node = 0; i_node < n_nodes; i_node++) {

		    array_1d<double, TDim>& U_i_curr = mvel_n1[i_node];
		    double delta_p_i = (mPn1[i_node] - mPn[i_node]) * rho_inv*factor;
		    const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];

		    //setting to zero
		    for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
			correction[l_comp] = 0.0;

		    //compute edge contributions dt*M^(-1)Gp
		    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
			unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
			double delta_p_j = (mPn1[j_neighbour] - mPn[j_neighbour]) * rho_inv*factor;

			CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

			edge_ij.Sub_grad_p(correction, delta_p_i, delta_p_j);
		    }
		    //compute prefactor
		    double coefficient = delta_t * m_inv;

		    //correct fractional momentum
		    for (unsigned int comp = 0; comp < TDim; comp++)
			U_i_curr[comp] += coefficient * correction[comp];
		
	    }

	    ApplyVelocityBC(mvel_n1);

	    //write velocity of time step n+1 to Kratos
	    mr_matrix_container.WriteVectorToDatabase(VELOCITY, mvel_n1, rNodes);



	    //calculate the error on the divergence
	    if(muse_mass_correction == true)
	    {
		#pragma omp parallel for private(correction) firstprivate(delta_t,rho_inv)
		for (int i_node = 0; i_node < n_nodes; i_node++)
		{
		    double& div_i_err = mdiv_error[i_node];
		    div_i_err = 0.0;

			const array_1d<double, TDim>& U_i_curr = mvel_n1[i_node];

			//compute edge contributions dt*M^(-1)Gp
			for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
			{
			    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
			    array_1d<double, TDim>& U_j_curr = mvel_n1[j_neighbour];

			    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

			    edge_ij.Add_D_v(div_i_err, U_i_curr*mRho, U_j_curr * mRho);
			}
		    
		}
	    }

//            KRATOS_WATCH("end of step3")
//            double vnorm2 = 0.0;
//                                for( int i = 0; i< rNodes.size(); i++)
//                                    vnorm2 += pow(mvel_n1[i][0],2) + pow(mvel_n1[i][1],2) + pow(mvel_n1[i][2],2);
//                                KRATOS_WATCH(sqrt(vnorm2));

	    KRATOS_CATCH("")
	}


	//************************************

	void ApplyVelocityBC(CalcVectorType& VelArray) {
	    KRATOS_TRY


	    if(mWallLawIsActive == false)
	    {
		//apply conditions on corner edges
		int edge_size = medge_nodes_direction.size();
		#pragma omp parallel for firstprivate(edge_size)
		for (int i = 0; i < edge_size; i++)
		{
		    int i_node = medge_nodes[i];
		    const array_1d<double, TDim>& direction = medge_nodes_direction[i];
			array_1d<double, TDim>& U_i = VelArray[i_node];
			double temp=0.0;
			for (unsigned int comp = 0; comp < TDim; comp++)
			    temp += U_i[comp] * direction[comp];

			for (unsigned int comp = 0; comp < TDim; comp++)
			    U_i[comp] = direction[comp]*temp;
		    
		}

		//apply conditions on corners
		int corner_size = mcorner_nodes.size();
		for (int i = 0; i < corner_size; i++)
		{
		    int i_node = mcorner_nodes[i];

		    array_1d<double, TDim>& U_i = VelArray[i_node];
			for (unsigned int comp = 0; comp < TDim; comp++)
			    U_i[comp] = 0.0;
		}
	    }


	    //slip condition
	    int slip_size = mSlipBoundaryList.size();
	    #pragma omp parallel for firstprivate(slip_size)
	    for (int i_slip = 0; i_slip < slip_size; i_slip++)
	    {
		unsigned int i_node = mSlipBoundaryList[i_slip];
		    array_1d<double, TDim>& U_i = VelArray[i_node];
		    array_1d<double, TDim>& an_i = mSlipNormal[i_node];
		    double projection_length = 0.0;
		    double normalization = 0.0;
		    for (unsigned int comp = 0; comp < TDim; comp++) {
			projection_length += U_i[comp] * an_i[comp];
			normalization += an_i[comp] * an_i[comp];
		    }
		    projection_length /= normalization;
		    //tangential momentum as difference between original and normal momentum
		    for (unsigned int comp = 0; comp < TDim; comp++)
			U_i[comp] -= projection_length * an_i[comp];
		
	    }

	    //fixed condition
	    int fixed_size = mFixedVelocities.size();
	    #pragma omp parallel for firstprivate(fixed_size)
	    for (int i_velocity = 0; i_velocity < fixed_size; i_velocity++)
	    {
		unsigned int i_node = mFixedVelocities[i_velocity];

		    const array_1d<double, TDim>& u_i_fix = mFixedVelocitiesValues[i_velocity];
		    array_1d<double, TDim>& u_i = VelArray[i_node];

		    for (unsigned int comp = 0; comp < TDim; comp++)
			u_i[comp] = u_i_fix[comp];
		
	    }



	    KRATOS_CATCH("")
	}



	//**************************************
	//function to calculate the area normals

	void CalculateNormals(ModelPart::ConditionsContainerType& rConditions) {
	    KRATOS_TRY

	    //calculate area normals face-by-face
	    array_1d<double, 3 > area_normal;
	    //2D case
	    if (TDim == 2)
	    {
		for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
		    CalculateNormal2D(cond_it, area_normal);
	    }//3D case
            else if (TDim == 3)
            {
                //help vectors for cross product
                array_1d<double, 3 > v1;
                array_1d<double, 3 > v2;
                for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
                    CalculateNormal3D(cond_it, area_normal, v1, v2);
            }

            //(re)initialize normals
            unsigned int n_nodes = mNodalFlag.size();
            mSlipNormal.resize(n_nodes);
            std::vector<bool> is_slip(n_nodes);
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                noalias(mSlipNormal[i_node]) = ZeroVector(TDim);
                is_slip[i_node] = false;
            }


            //loop over all faces
            const double node_factor = 1.0 / TDim;
            for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
            {
                //get geometry data of the face
                Geometry<Node < 3 > >& face_geometry = cond_it->GetGeometry();

                //reference for area normal of the face
                array_1d<double, 3 > & face_normal = cond_it->GetValue(NORMAL);

                //slip condition
                if (cond_it->GetValue(IS_STRUCTURE) == true)
                    for (unsigned int if_node = 0; if_node < TDim; if_node++)
                    {
                        unsigned int i_node = static_cast<unsigned int> (face_geometry[if_node].FastGetSolutionStepValue(AUX_INDEX));
                        array_1d<double, TDim>& slip_normal = mSlipNormal[i_node];
                        is_slip[i_node] = true;
                        for (unsigned int comp = 0; comp < TDim; comp++)
                        {
                            slip_normal[comp] += node_factor * face_normal[comp];
                        }
                    }
            }

            //fill the list of slip nodes
	    std::vector< unsigned int> tempmSlipBoundaryList;
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                if (is_slip[i_node] == true)
                    tempmSlipBoundaryList.push_back(i_node);
            }
            mSlipBoundaryList.resize(tempmSlipBoundaryList.size(),false);
	    #pragma omp parallel for
	    for(unsigned int i=0; i<tempmSlipBoundaryList.size(); i++)
	      mSlipBoundaryList[i] = tempmSlipBoundaryList[i];



            KRATOS_CATCH("")
        }











        //*******************************
        //function to free dynamic memory

        void Clear()
        {
            KRATOS_TRY
            mWork.clear();
            mvel_n.clear();
            mvel_n1.clear();
            mPn.clear();
            mPn1.clear();
            mHmin.clear();
            mHavg.clear();
            mSlipNormal.clear();
            mNodalFlag.clear();
            mFixedVelocities.clear();
            mFixedVelocitiesValues.clear();
            mPressureOutletList.clear();
            //	    mPressureOutlet.clear();
            mSlipBoundaryList.clear();
            mL.clear();
            mTauPressure.clear();
            mTauConvection.clear();
            mTau2.clear();

            mBeta.clear();

            mdiv_error.clear();

            KRATOS_CATCH("")
        }

        void ActivateWallResistance(double Ywall)
        {
            mWallLawIsActive = true;
            mY_wall = Ywall;
        }

        void ComputePressureStabilization()
        {
            KRATOS_TRY

            //PREREQUISITES

            //variables for node based data handling
            ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
            int n_nodes = rNodes.size();
            //storage of nodal values in local variables
            CalcVectorType rhs;
            rhs.resize(n_nodes);


            //read velocity and pressure data from Kratos
            mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes);

            //read time step size from Kratos
//            ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
//            double delta_t = CurrentProcessInfo[DELTA_TIME];


            //compute intrinsic time
            //            double time_inv = 1.0 / delta_t;
            double time_inv_avg = 1.0 / mdelta_t_avg;

            double stabdt_pressure_factor = mstabdt_pressure_factor;

            KRATOS_WATCH(stabdt_pressure_factor);
#pragma omp parallel for firstprivate(time_inv_avg,stabdt_pressure_factor)
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                double& h_avg_i = mHavg[i_node];
                array_1d<double, TDim>& a_i = mvel_n1[i_node];
                const double nu_i = mViscosity;

                double vel_norm = norm_2(a_i);

                double tau = 1.0 / (2.0 * vel_norm / h_avg_i + stabdt_pressure_factor * time_inv_avg + (4.0 * nu_i) / (h_avg_i * h_avg_i));
                mTauPressure[i_node] = tau;
            }

            KRATOS_CATCH("");
        }


        	//*********************************************************************
	//function to calculate right-hand side of fractional momentum equation

	void ViscosityCorrectionStep()
	{
	    KRATOS_TRY

	    int n_nodes = mvel_n1.size();
            ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
            mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes);

            ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
            double delta_t = CurrentProcessInfo[DELTA_TIME];

            CalcVectorType rhs;
	    rhs.resize(n_nodes);

	    //calculating the RHS
//	    double inverse_rho = 1.0 / mRho;
	    #pragma omp parallel for 
	    for (int i_node = 0; i_node < n_nodes; i_node++)
            {
		    array_1d<double, TDim>& rhs_i = rhs[i_node];
		    const array_1d<double, TDim>& U_i = mvel_n1[i_node];

		    //initializing with the external forces (e.g. gravity)
//		    double& m_i = mr_matrix_container.GetLumpedMass()[i_node];
		    for (unsigned int comp = 0; comp < TDim; comp++)
			rhs_i[comp] = 0.0 ;

		    //convective term
		    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
			unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
			    const array_1d<double, TDim>& U_j = mvel_n1[j_neighbour];

			    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

			    edge_ij.Sub_ViscousContribution(rhs_i, U_i, mViscosity, U_j, mViscosity);
		    }
	    }

            //correcting the velocity
 	    mr_matrix_container.Add_Minv_value(mvel_n1, mvel_n1, delta_t, mr_matrix_container.GetInvertedMass(), rhs);
	    ApplyVelocityBC(mvel_n1);

	    mr_matrix_container.WriteVectorToDatabase(VELOCITY, mvel_n1, rNodes);
	    KRATOS_CATCH("")
	}

        void ComputeViscousForces()
	{
	    KRATOS_TRY

            if (mr_model_part.NodesBegin()->SolutionStepsDataHas(FORCE) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----FORCE---- variable!!!!!! ERROR", "");

	    int n_nodes = mvel_n1.size();
            ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
            mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes);

            ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
            double delta_t = CurrentProcessInfo[DELTA_TIME];

            CalcVectorType rhs;
	    rhs.resize(n_nodes);

	    //calculating the RHS
//	    double inverse_rho = 1.0 / mRho;
	    #pragma omp parallel for
	    for (int i_node = 0; i_node < n_nodes; i_node++)
            {
		    array_1d<double, TDim>& rhs_i = rhs[i_node];
		    const array_1d<double, TDim>& U_i = mvel_n1[i_node];

		    //initializing with the external forces (e.g. gravity)
//		    double& m_i = mr_matrix_container.GetLumpedMass()[i_node];
		    for (unsigned int comp = 0; comp < TDim; comp++)
			rhs_i[comp] = 0.0 ;

		    //convective term
		    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
			unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
			    const array_1d<double, TDim>& U_j = mvel_n1[j_neighbour];

			    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

			    edge_ij.Sub_ViscousContribution(rhs_i, U_i, mViscosity, U_j, mViscosity);
		    }

                    const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];

                    for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                        rhs_i[l_comp] *= m_inv;
	    }

	    mr_matrix_container.WriteVectorToDatabase(FORCE, rhs, rNodes);
	    KRATOS_CATCH("")
	}

    private:
        MatrixContainer& mr_matrix_container;
        ModelPart& mr_model_part;

        bool muse_mass_correction;

        //parameters controlling the wall law
        bool mWallLawIsActive;
        bool mY_wall;

        //parameters for controlling the usage of the delta time in the stabilization
        double mstabdt_pressure_factor;
        double mstabdt_convection_factor;
        double medge_detection_angle;
        double mtau2_factor;
        bool massume_constant_dp;

        //nodal values
        //velocity vector U at time steps n and n+1
        CalcVectorType mWork, mvel_n, mvel_n1, mx;
        //pressure vector p at time steps n and n+1
        ValuesVectorType mPn, mPn1;

        //minimum length of the edges surrounding edges surrounding each nodal point
        ValuesVectorType mHmin;
        ValuesVectorType mHavg;
        CalcVectorType mEdgeDimensions;

        //area normal
        CalcVectorType mSlipNormal;
        //projection terms
        CalcVectorType mPi, mXi;

        //flag for first time step
        bool mFirstStep;

        //flag to differentiate interior and boundary nodes
        ValuesVectorType mNodalFlag;
        //lists of nodes with different types of boundary conditions
        IndicesVectorType mSlipBoundaryList, mPressureOutletList, mFixedVelocities;
        CalcVectorType mFixedVelocitiesValues;
        //	ValuesVectorType mPressureOutlet;

        //intrinsic time step size
        ValuesVectorType mTauPressure;
        ValuesVectorType mTauConvection;
        ValuesVectorType mTau2;

        ValuesVectorType mdiv_error;

        //variables for resolving pressure equation
        //laplacian matrix
        TSystemMatrixType mL;

        //constant variables
        double mRho;
        double mViscosity;
        array_1d<double, TDim> mBodyForce;


        //variables for convection
        ValuesVectorType mBeta;

        //variables for edge BCs
        IndicesVectorType medge_nodes;
        CalcVectorType medge_nodes_direction;
        IndicesVectorType mcorner_nodes;


        double mdelta_t_avg;
        double max_dt;


        //***********************************************************
        //functions to calculate area normals for boundary conditions

        void CalculateNormal2D(ModelPart::ConditionsContainerType::iterator cond_it, array_1d<double, 3 > & area_normal)
        {
            Geometry<Node < 3 > >& face_geometry = (cond_it)->GetGeometry();

            area_normal[0] = face_geometry[1].Y() - face_geometry[0].Y();
            area_normal[1] = -(face_geometry[1].X() - face_geometry[0].X());
            area_normal[2] = 0.00;

            noalias((cond_it)->GetValue(NORMAL)) = area_normal;
        }

        void CalculateNormal3D(ModelPart::ConditionsContainerType::iterator cond_it, array_1d<double, 3 > & area_normal, array_1d<double, 3 > & v1, array_1d<double, 3 > & v2)
        {
            Geometry<Node < 3 > >& face_geometry = (cond_it)->GetGeometry();

            v1[0] = face_geometry[1].X() - face_geometry[0].X();
            v1[1] = face_geometry[1].Y() - face_geometry[0].Y();
            v1[2] = face_geometry[1].Z() - face_geometry[0].Z();

            v2[0] = face_geometry[2].X() - face_geometry[0].X();
            v2[1] = face_geometry[2].Y() - face_geometry[0].Y();
            v2[2] = face_geometry[2].Z() - face_geometry[0].Z();

            MathUtils<double>::CrossProduct(area_normal, v1, v2);
            area_normal *= -0.5;

            noalias((cond_it)->GetValue(NORMAL)) = area_normal;
        }


        //*********************************************************
        //function to calculate minimum length of surrounding edges

        void CalculateEdgeLengths(ModelPart::NodesContainerType& rNodes)
        {
            KRATOS_TRY

                    //get number of nodes
                    unsigned int n_nodes = rNodes.size();
            //reserve memory for storage of nodal coordinates
            std::vector< array_1d<double, 3 > > position;
            position.resize(n_nodes);

            //get position of all nodes
            for (typename ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
            {
                //get the global index of the node
                unsigned int i_node = static_cast<unsigned int> (node_it->FastGetSolutionStepValue(AUX_INDEX));
                //save its coordinates locally
                noalias(position[i_node]) = node_it->Coordinates();

                //initialize minimum edge length with relatively big values
                mHmin[i_node] = 1e10;
            }

            ValuesVectorType& aaa = mr_matrix_container.GetHmin();
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                mHmin[i_node] = aaa[i_node];
            }

            //take unstructured meshes into account
            if (TDim == 2)
            {
                for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
                {
                    double& h_i = mHavg[i_node];
                    double& m_i = mr_matrix_container.GetLumpedMass()[i_node];
                    // 						double& rho_i = mRho[i_node];

                    h_i = sqrt(2.0 * m_i);
                }
            } else if (TDim == 3)
            {
                for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
                {
                    double& h_i = mHavg[i_node];
                    double& m_i = mr_matrix_container.GetLumpedMass()[i_node];
                    // 						double& rho_i = mRho[i_node];

                    h_i = pow(6.0 * m_i, 1.0 / 3.0);
                }
            }

            //compute edge coordinates
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                array_1d<double, 3 > & pos_i = position[i_node];

                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    array_1d<double, 3 > & pos_j = position[j_neighbour];

                    array_1d<double, TDim>& l_k = mEdgeDimensions[csr_index];
                    for (unsigned int comp = 0; comp < TDim; comp++)
                        l_k[comp] = pos_i[comp] - pos_j[comp];
                }
            }

            KRATOS_CATCH("")
        }







        //**************************************

        void CornerDectectionHelper(Geometry< Node < 3 > >& face_geometry,
                const array_1d<double, 3 > & face_normal,
                const double An,
                const WeakPointerVector<Condition>& neighb,
                const unsigned int i1,
                const unsigned int i2,
                const unsigned int neighb_index,
                std::vector<unsigned int>& edge_nodes,
                CalcVectorType& cornern_list
                )
        {
            double acceptable_angle = 45.0 / 180.0 * 3.1; //angles of less than 45 deg will be accepted
            double acceptable_cos = cos(acceptable_angle);

            if (face_geometry[i1].Id() < face_geometry[i2].Id()) //we do this to add the face ones
            {
                const array_1d<double, 3 > & neighb_normal = neighb[neighb_index].GetValue(NORMAL);
                double neighb_An = norm_2(neighb_normal);

                double cos_normal = 1.0 / (An * neighb_An) * inner_prod(face_normal, neighb_normal);

                //if the angle is too big between the two normals then the edge in the middle is a corner
                if (cos_normal < acceptable_cos)
                {
                    array_1d<double, 3 > edge = face_geometry[i2].Coordinates() - face_geometry[i1].Coordinates();
                    double temp = norm_2(edge);
                    edge /= temp;

                    int index1 = face_geometry[i1].FastGetSolutionStepValue(AUX_INDEX);
                    int index2 = face_geometry[i2].FastGetSolutionStepValue(AUX_INDEX);

                    edge_nodes[index1] += 1;
                    edge_nodes[index2] += 1;

                    double sign1 = inner_prod(cornern_list[index1], edge);
                    if (sign1 >= 0)
                        cornern_list[index1] += edge;
                    else
                        cornern_list[index1] -= edge;

                    double sign2 = inner_prod(cornern_list[index2], edge);
                    if (sign2 >= 0)
                        cornern_list[index2] += edge;
                    else
                        cornern_list[index2] -= edge;


                }



            }


        }

        //function to calculate the area normals

        void DetectEdges3D(ModelPart::ConditionsContainerType& rConditions)
        {
            KRATOS_TRY

            //calculate area normals face-by-face
            array_1d<double, 3 > area_normal;

            //(re)initialize normals
            unsigned int n_nodes = mNodalFlag.size();
            std::vector<unsigned int> temp_edge_nodes(n_nodes);
            CalcVectorType temp_cornern_list(n_nodes);
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                temp_edge_nodes[i_node] = 0.0;
                noalias(temp_cornern_list[i_node]) = ZeroVector(TDim);
            }

            //loop over all faces
            //            const double node_factor = 1.0 / TDim;
            for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
            {
                //get geometry data of the face
                Geometry<Node < 3 > >& face_geometry = cond_it->GetGeometry();

                //reference for area normal of the face
                const array_1d<double, 3 > & face_normal = cond_it->GetValue(NORMAL);
                double An = norm_2(face_normal);

                unsigned int current_id = cond_it->Id();

                //slip condition
                if (cond_it->GetValue(IS_STRUCTURE) == 1.0) //this is a slip face --> now look for its neighbours
                {
                    const WeakPointerVector<Condition>& neighb = cond_it->GetValue(NEIGHBOUR_CONDITIONS);

                    //check for neighbour zero
                    if (neighb[0].Id() != current_id) //check if the neighbour exists
                        CornerDectectionHelper(face_geometry, face_normal, An, neighb, 1, 2, 0, temp_edge_nodes, temp_cornern_list);

                    //check for neighbour one
                    if (neighb[1].Id() != current_id) //check if the neighbour exists
                        CornerDectectionHelper(face_geometry, face_normal, An, neighb, 2, 0, 1, temp_edge_nodes, temp_cornern_list);

                    //check for neighbour two
                    if (neighb[2].Id() != current_id) //check if the neighbour exists
                        CornerDectectionHelper(face_geometry, face_normal, An, neighb, 0, 1, 2, temp_edge_nodes, temp_cornern_list);

                }
            }

            //fill the list of edge_nodes
            std::vector<unsigned int> tempmedge_nodes;
            std::vector< array_1d<double,TDim> > tempmedge_nodes_direction;
            std::vector<unsigned int> tempmcorner_nodes;
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                if (temp_edge_nodes[i_node] == 2) //node is a edge_node
                {
                    tempmedge_nodes.push_back(i_node);
                    array_1d<double, TDim>& node_edge = temp_cornern_list[i_node];

                    node_edge /= norm_2(node_edge);
                    tempmedge_nodes_direction.push_back(node_edge);
                } else if (temp_edge_nodes[i_node] > 2)
                    tempmcorner_nodes.push_back(i_node);
            }
            medge_nodes.resize(tempmedge_nodes.size(),false);
            medge_nodes_direction.resize(tempmedge_nodes_direction.size(),false);
            mcorner_nodes.resize(tempmcorner_nodes.size(),false);
	    #pragma omp parallel for
	    for (unsigned int i = 0; i < tempmedge_nodes.size(); i++)
	    {
	      medge_nodes[i] = tempmedge_nodes[i];
	      medge_nodes_direction[i] = tempmedge_nodes_direction[i];
	    }
	    #pragma omp parallel for
	    for (unsigned int i = 0; i < tempmcorner_nodes.size(); i++)
	    {
	      mcorner_nodes[i] = tempmcorner_nodes[i];
	    }



            for (unsigned int i = 0; i < mcorner_nodes.size(); i++)
            {
                KRATOS_WATCH(mcorner_nodes[i]);

            }


            KRATOS_CATCH("")
        }

        void ComputeWallResistance(
                const CalcVectorType& vel,
                CalcVectorType& rhs
                )
        {
            //parameters:
            double k = 0.41;
            double B = 5.1;
            double density = mRho;
            double mu = mViscosity;
            double toll = 1e-6;
            double ym = mY_wall; //0.0825877; //0.0093823
            double y_plus_incercept = 10.9931899;
            unsigned int itmax = 100;

            if (mu == 0)
                KRATOS_ERROR(std::logic_error, "it is not possible to use the wall law with 0 viscosity", "");

            //slip condition
            int slip_size = mSlipBoundaryList.size();
#pragma omp parallel for firstprivate(slip_size,B,density,mu,toll,ym,y_plus_incercept,itmax)
            for (int i_slip = 0; i_slip < slip_size; i_slip++)
            {
                unsigned int i_node = mSlipBoundaryList[i_slip];

                array_1d<double, TDim>& rhs_i = rhs[i_node];
                const array_1d<double, TDim>& U_i = vel[i_node];
                const array_1d<double, TDim>& an_i = mSlipNormal[i_node];

                //compute the modulus of the velocity
                double mod_vel = 0.0;
                double area = 0.0;
                for (unsigned int comp = 0; comp < TDim; comp++)
                {
                    mod_vel += U_i[comp] * U_i[comp];
                    area += an_i[comp] * an_i[comp];
                }
                mod_vel = sqrt(mod_vel);
                area = sqrt(area);

                //now compute the skin friction
                double mod_uthaw = sqrt(mod_vel * mu / ym);
                const double y_plus = ym * mod_uthaw / mu;

                if (y_plus > y_plus_incercept)
                {
                    //begin cicle to calculate the real u_thaw's module:
                    unsigned int it = 0;
                    double dx = 1e10;
                    //                        KRATOS_WATCH(fabs(dx));
                    while (fabs(dx) > toll * mod_uthaw && it < itmax)
                    {
                        double a = 1.0 / k;
                        double temp = a * log(ym * mod_uthaw / mu) + B;
                        double y = mod_uthaw * (temp) - mod_vel;
                        double y1 = temp + a;
                        dx = y / y1;
                        mod_uthaw -= dx;
                        it = it + 1;
                    }

                    //                         KRATOS_WATCH(toll*mod_uthaw);
                    //                         KRATOS_WATCH(area);
                    //                        KRATOS_WATCH(it);
                    if (it == itmax)
                        std::cout << "attention max number of iterations exceeded in wall law computation" << std::endl;


                }
                //                    else
                //                    {
                //                        for (unsigned int comp = 0; comp < TDim; comp++)
                //                            rhs_i[comp] -= U_i[comp] * area * mu  / (density*ym) ;
                //                    }

                if (mod_vel > 1e-12)
                    for (unsigned int comp = 0; comp < TDim; comp++)
                        rhs_i[comp] -= U_i[comp] * area * mod_uthaw * mod_uthaw * density / (mod_vel);




            }
        }





    };
} //namespace Kratos
//#undef SYMM_PRESS
#endif //KRATOS_EDGEBASED_FLUID_SOLVER_H_INCLUDED defined


