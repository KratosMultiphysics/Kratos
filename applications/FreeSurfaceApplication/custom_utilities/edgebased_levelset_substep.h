// Kratos Multi-Physics
// 
// Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement: 
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 	
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY 
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2009-01-14 16:24:38 $
//   Revision:            $Revision: 1.11 $
//
//
#if !defined(KRATOS_EDGEBASED_LEVELSET_SUBSTEP_FLUID_SOLVER_H_INCLUDED)
#define  KRATOS_EDGEBASED_LEVELSET_SUBSTEP_FLUID_SOLVER_H_INCLUDED

// #define DEBUG_OUTPUT

//#define SPLIT_OSS
// #define SYMM_PRESS
// System includes
#include <string>
#include <iostream>
#include <algorithm>
// #include <omp.h>
// External includes
// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/node.h"
//#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
#include "free_surface_application.h"
#include "custom_utilities/edge_data_c2c.h"

namespace Kratos
{
template<unsigned int TDim, class MatrixContainer, class TSparseSpace, class TLinearSolver>
class EdgeBasedLevelSetSubstep
{
public:
    //name for the self defined structure
    typedef EdgesStructureTypeC2C<TDim> CSR_Tuple;
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
    EdgeBasedLevelSetSubstep (MatrixContainer& mr_matrix_container,
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
        : mr_matrix_container (mr_matrix_container),
          mr_model_part (mr_model_part),
          mstabdt_pressure_factor (stabdt_pressure_factor),
          mstabdt_convection_factor (stabdt_convection_factor),
          medge_detection_angle (edge_detection_angle),
          mtau2_factor (tau2_factor),
          massume_constant_dp (assume_constant_dp)
    {
        for (ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); it++)
            it->FastGetSolutionStepValue (VISCOSITY) = viscosity;

	mMolecularViscosity = viscosity;
// 	    mViscosity = viscosity;
        noalias (mBodyForce) = body_force;
        mRho = density;
        mdelta_t_avg = 1000.0;
        max_dt = 1.0;
        muse_mass_correction = use_mass_correction;
        mshock_coeff = 0.7;
        mWallLawIsActive = false;
        mnumsubsteps=5;
	mmax_dt = 0.0;
        mcorner_coefficient = 30.0; //50.0;
        medge_coefficient = 2.0; //30.0; //10.0;
//            for (unsigned int i = 0; i < TDim; i++) mBodyForce[i] = 0;
//            mBodyForce[1] = -9.81;
//
//            mRho = 1000.0;
		std::cout << "Edge based level set substep solver is created" << std::endl;
    };
    ~EdgeBasedLevelSetSubstep()
    {
    };
	void SetBodyForce( const Vector& body_force)
	{
		noalias(mBodyForce) = body_force;
		KRATOS_WATCH(mBodyForce);
	}

    //***********************************
    //function to initialize fluid solver
    void Initialize (
    )
    {
        KRATOS_TRY
        //get number of nodes
        unsigned int n_nodes = mr_model_part.Nodes().size();
        unsigned int n_edges = mr_matrix_container.GetNumberEdges();
        //size data vectors
        mViscosity.resize (n_nodes);
        mr_matrix_container.SetToZero (mViscosity);
        mWork.resize (n_nodes);
        mr_matrix_container.SetToZero (mWork);
        mvel_n.resize (n_nodes);
        mr_matrix_container.SetToZero (mvel_n);
        mvel_n1.resize (n_nodes);
        mr_matrix_container.SetToZero (mvel_n1);
        mPn.resize (n_nodes);
        mr_matrix_container.SetToZero (mPn);
        mPn1.resize (n_nodes);
        mr_matrix_container.SetToZero (mPn1);
        mHmin.resize (n_nodes);
        mr_matrix_container.SetToZero (mHmin);
        mHavg.resize (n_nodes);
        mr_matrix_container.SetToZero (mHavg);
        mNodalFlag.resize (n_nodes);
        mr_matrix_container.SetToZero (mNodalFlag);
        mdistances.resize (n_nodes);
        mr_matrix_container.SetToZero (mdistances);
        mTauPressure.resize (n_nodes);
        mr_matrix_container.SetToZero (mTauPressure);
        mTauConvection.resize (n_nodes);
        mr_matrix_container.SetToZero (mTauConvection);
        mTau2.resize (n_nodes);
        mr_matrix_container.SetToZero (mTau2);
        mPi.resize (n_nodes);
        mr_matrix_container.SetToZero (mPi);
        mXi.resize (n_nodes);
        mr_matrix_container.SetToZero (mXi);
        mx.resize (n_nodes);
        mr_matrix_container.SetToZero (mx);
        mEdgeDimensions.resize (n_edges);
        mr_matrix_container.SetToZero (mEdgeDimensions);
        //convection variables
        mBeta.resize (n_nodes);
        mr_matrix_container.SetToZero (mBeta);
        mPiConvection.resize (n_nodes);
        mr_matrix_container.SetToZero (mPiConvection);
        mphi_n.resize (n_nodes);
        mr_matrix_container.SetToZero (mphi_n);
        mphi_n1.resize (n_nodes);
        mr_matrix_container.SetToZero (mphi_n1);
        mEps.resize (n_nodes);
        mr_matrix_container.SetToZero (mEps);
//         mD.resize(n_nodes);
// 	mr_matrix_container.SetToZero(mD);
        mA.resize (n_nodes);
        mr_matrix_container.SetToZero (mA);
        mB.resize (n_nodes);
        mr_matrix_container.SetToZero (mB);
        mdiv_error.resize (n_nodes);
        mr_matrix_container.SetToZero (mdiv_error);
        mWallReductionFactor.resize (n_nodes);
        mr_matrix_container.SetToZero (mWallReductionFactor);
        mdiag_stiffness.resize (n_nodes);
        mr_matrix_container.SetToZero (mdiag_stiffness);
        mis_slip.resize (n_nodes);
        mis_visited.resize (n_nodes);
        macc.resize (n_nodes);
        mr_matrix_container.SetToZero (macc);
//	    ValuesVectorType external_pressure;
//	    external_pressure.resize(n_nodes);
        //read velocity and pressure data from Kratos
        mr_matrix_container.FillScalarFromDatabase (VISCOSITY, mViscosity, mr_model_part.Nodes() );
        mr_matrix_container.FillVectorFromDatabase (VELOCITY, mvel_n1, mr_model_part.Nodes() );
        mr_matrix_container.FillScalarFromDatabase (PRESSURE, mPn1, mr_model_part.Nodes() );
        mr_matrix_container.FillOldScalarFromDatabase (PRESSURE, mPn, mr_model_part.Nodes() );
        mr_matrix_container.FillOldVectorFromDatabase (VELOCITY, mvel_n, mr_model_part.Nodes() );
        mr_matrix_container.FillCoordinatesFromDatabase (mx, mr_model_part.Nodes() );
        //set flag for first time step
        mFirstStep = true;
        //loop to categorize boundary nodes
        std::vector< unsigned int> tempFixedVelocities;
        std::vector< array_1d<double,TDim> > tempFixedVelocitiesValues;
        std::vector< unsigned int> tempPressureOutletList;
	std::vector< unsigned int> tempDistanceList;
        for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode != mr_model_part.NodesEnd();
                inode++)
        {
            int index = inode->FastGetSolutionStepValue (AUX_INDEX);
            if (inode->IsFixed (VELOCITY_X) ) //note that the variables can be either all fixed or no one fixed
            {
                if (inode->IsFixed (VELOCITY_Y) == false || inode->IsFixed (VELOCITY_Z) == false)
                {
                    std::cout << "error found on the fixity of node " << inode->Id() << std::endl;
                    KRATOS_THROW_ERROR (std::logic_error, "velocities can be either all fixed or none fixed", "")
                }
                tempFixedVelocities.push_back (index);
                tempFixedVelocitiesValues.push_back (mvel_n1[index]);
            }
            if (inode->IsFixed (DISTANCE) )
                tempDistanceList.push_back (index);
            if (inode->IsFixed (PRESSURE) )
            {
                tempPressureOutletList.push_back (index);
//		    mPressureOutlet.push_back(external_pressure[index]);
            }
        }
        mFixedVelocities.resize (tempFixedVelocities.size(),false);
        mFixedVelocitiesValues.resize (tempFixedVelocitiesValues.size(),false);
        mPressureOutletList.resize (tempPressureOutletList.size(),false);
        mDistanceBoundaryList.resize (tempDistanceList.size(),false);
        mDistanceValuesList.resize (tempDistanceList.size(),false);
        #pragma omp parallel for
        for (int i=0; i<static_cast<int> (tempFixedVelocities.size() ); i++)
        {
            mFixedVelocities[i] = tempFixedVelocities[i];
            mFixedVelocitiesValues[i] = tempFixedVelocitiesValues[i];
        }
        #pragma omp parallel for
        for (int i=0; i<static_cast<int> (tempPressureOutletList.size() ); i++)
        {
            mPressureOutletList[i] = tempPressureOutletList[i];
        }
        for (int i=0; i<static_cast<int> (tempDistanceList.size() ); i++)
        {
            mDistanceBoundaryList[i] = tempDistanceList[i];
        }
        //compute slip normals and fill SlipList
        CalculateNormals (mr_model_part.Conditions() );
        mr_matrix_container.WriteVectorToDatabase (NORMAL, mSlipNormal, mr_model_part.Nodes() );
        if (TDim == 3)
            DetectEdges3D (mr_model_part.Conditions() );
        //print number of nodes corresponding to the different types of boundary conditions
        //				KRATOS_WATCH(mFixedVelocities.size())
        //				KRATOS_WATCH(mPressureOutletList.size())
        //				KRATOS_WATCH(mSlipBoundaryList.size())
        //determine number of edges and entries
        unsigned int n_nonzero_entries = 2 * n_edges + n_nodes;
        //allocate memory for variables
        mL.resize (n_nodes, n_nodes, n_nonzero_entries);
        int number_of_threads= OpenMPUtils::GetNumThreads();
        std::vector<int> row_partition (number_of_threads);
        OpenMPUtils::DivideInPartitions (n_nodes,number_of_threads,row_partition);
        for (int k = 0; k < number_of_threads; k++)
        {
            #pragma omp parallel
            if (OpenMPUtils::ThisThread() == k)
            {
                for (int i_node = static_cast<int> (row_partition[k]); i_node < static_cast<int> (row_partition[k + 1]); i_node++)
                {
                    //loop over all nodes
// 	    for (unsigned int i_node = 0; i_node < n_nodes; i_node++) {
                    //flag for considering diagonal matrix elements
                    bool flag = 0;
                    //loop over all neighbours
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
                    {
                        //get global index of neighbouring node j
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                        //define matrix structure row by row (the order does matter!)
                        if ( (static_cast<int> (j_neighbour) > i_node) && (flag == 0) )
                        {
                            //add diagonal/nodal contribution
                            mL.push_back (i_node, i_node, 0.0);
                            flag = 1;
                        }
                        //add non-diagonal/edge contribution
                        mL.push_back (i_node, j_neighbour, 0.0);
                    }
                    //if diagonal element is the last non-zero element of the row
                    if (flag == 0)
                        mL.push_back (i_node, i_node, 0.0);
                }
            }
        }
        //compute minimum length of the surrounding edges
        CalculateEdgeLengths (mr_model_part.Nodes() );
        //set the pressure projection to the body force value
        array_1d<double,3> temp = mRho * mBodyForce;
        for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode != mr_model_part.NodesEnd();
                inode++)
            inode->FastGetSolutionStepValue (PRESS_PROJ) = temp;
        mr_matrix_container.FillScalarFromDatabase (POROSITY, mEps, mr_model_part.Nodes() );
        //verify that neither h_min nor havg are 0
        for (unsigned int i_node=0; i_node<mHmin.size(); i_node++)
	{
            if (mHmin[i_node] < 1e-20) KRATOS_THROW_ERROR ( std::logic_error,"hmin too small on node ",i_node+1)
                if (mHavg[i_node] < 1e-20) KRATOS_THROW_ERROR ( std::logic_error,"havg too small on node ",i_node+1)
                    if (mHmin[i_node] > 1e20) KRATOS_THROW_ERROR ( std::logic_error,"hmin too big on node ",i_node+1)
                        if (mHavg[i_node] > 1e20) KRATOS_THROW_ERROR ( std::logic_error,"havg too big on node ",i_node+1)
	}
        for (ModelPart::ElementsContainerType::iterator it=mr_model_part.ElementsBegin(); it!=mr_model_part.ElementsEnd(); it++)
	{
	  if (it->Id() < 1)
	  {
                KRATOS_THROW_ERROR (std::logic_error, "Element found with Id 0 or negative","")
	  }
	  double elem_vol = 0.0;
            if (TDim == 2)
	    elem_vol = it->GetGeometry().Area();
            else
	    elem_vol = it->GetGeometry().Volume();
	  if (elem_vol <= 0)
	  {
	      std::cout << "error on element -> " << it->Id() << std::endl;
                KRATOS_THROW_ERROR (std::logic_error, "Area can not be lesser than 0","")
	  }
	}
        KRATOS_CATCH ("")
    }
    void SetShockCapturingCoefficient (double coeff)
    {
        mshock_coeff = coeff;
    }

    void GatherValues()
    {
        KRATOS_TRY

        mr_matrix_container.FillScalarFromDatabase (VISCOSITY, mViscosity, mr_model_part.Nodes() );
        mr_matrix_container.FillScalarFromDatabase (POROSITY, mEps, mr_model_part.Nodes() );
        mr_matrix_container.FillScalarFromDatabase (PRESSURE, mPn1, mr_model_part.Nodes() );
        mr_matrix_container.FillScalarFromDatabase (DISTANCE, mdistances, mr_model_part.Nodes() );
        mr_matrix_container.FillVectorFromDatabase (VELOCITY, mvel_n1, mr_model_part.Nodes() );
        mr_matrix_container.FillVectorFromDatabase(PRESS_PROJ, mXi, mr_model_part.Nodes());

        mr_matrix_container.FillOldVectorFromDatabase (VELOCITY, mvel_n, mr_model_part.Nodes() );
        mr_matrix_container.FillOldScalarFromDatabase (PRESSURE, mPn, mr_model_part.Nodes() );

        KRATOS_CATCH("")
    }

    //***************************************
    //function to set adequate time step size
    double ComputeTimeStep (const double CFLNumber, const double MaxDt)
    {
        KRATOS_TRY
        //save the maximum time step
        max_dt = MaxDt;
        //local variable for time step size
        //getting value of current velocity and of viscosity
//          mr_matrix_container.FillScalarFromDatabase (VISCOSITY, mViscosity, mr_model_part.Nodes() );
//         mr_matrix_container.FillScalarFromDatabase (POROSITY, mEps, mr_model_part.Nodes() );
//         mr_matrix_container.FillScalarFromDatabase (PRESSURE, mPn1, mr_model_part.Nodes() );
//         mr_matrix_container.FillScalarFromDatabase (DISTANCE, mdistances, mr_model_part.Nodes() );
//         mr_matrix_container.FillVectorFromDatabase (VELOCITY, mvel_n1, mr_model_part.Nodes() );
//            mr_matrix_container.FillVectorFromDatabase(PRESS_PROJ, mXi, mr_model_part.Nodes());
//
//         mr_matrix_container.FillOldVectorFromDatabase (VELOCITY, mvel_n, mr_model_part.Nodes() );
//         mr_matrix_container.FillOldScalarFromDatabase (PRESSURE, mPn, mr_model_part.Nodes() );

//         mr_matrix_container.FillScalarFromDatabase(DIAMETER, mD, mr_model_part.Nodes());
//         mr_matrix_container.FillScalarFromDatabase (LIN_DARCY_COEF, mA, mr_model_part.Nodes() );
//         mr_matrix_container.FillScalarFromDatabase (NONLIN_DARCY_COEF, mB, mr_model_part.Nodes() );
//            double delta_t_i = delta_t;
        //*******************
        //loop over all nodes
        int n_nodes = static_cast<int>(mvel_n1.size());

        unsigned int n_proc = OpenMPUtils::GetNumThreads();
        Vector dt_avg_vec(n_proc,1e10);
        Vector dt_vec(n_proc,1e10);
        Vector dt_avg_novisc_vec(n_proc,1e10);

        #pragma omp parallel for firstprivate(n_nodes)
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            unsigned int my_id = OpenMPUtils::ThisThread();
            double& delta_t = dt_vec[my_id];
            double& mdelta_t_avg = dt_avg_vec[my_id];
            double& delta_t_avg_novisc = dt_avg_novisc_vec[my_id];

            const array_1d<double, TDim>& v_i = mvel_n1[i_node];
            const double havg_i = mHavg[i_node];
            const double hmin_i = mHmin[i_node];
            const double eps_i = mEps[i_node];
            //const double d_i = mD[i_node];
            double nu = mViscosity[i_node];
//            const double lindarcy_i = mA[i_node];
//            const double nonlindarcy_i = mB[i_node];
            double vel_norm = norm_2 (v_i);
            //double porosity_coefficient = ComputePorosityCoefficient(nu, vel_norm, eps_i, d_i);
//            double porosity_coefficient = ComputePorosityCoefficient( vel_norm, eps_i, lindarcy_i, nonlindarcy_i);
            vel_norm /= eps_i;
            //use CFL condition to compute time step size
            double delta_t_i = 1.0 / (vel_norm /hmin_i + nu / (hmin_i * hmin_i) /*+ porosity_coefficient*/);
            double delta_t_i_avg = 1.0 / (vel_norm /havg_i + nu / (havg_i * havg_i) /*+ porosity_coefficient*/);
            double delta_t_i_avg_novisc = 1.0 / (2.0 * vel_norm /havg_i );
            //considering the most restrictive case of neighbor's velocities with similar direction but opposite sense.
            //loop over all neighbours
            for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
            {
                //get global index of neighbouring node j
                unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                const array_1d<double, TDim>& v_j = mvel_n1[j_neighbour];
                double v_diff_norm = 0.0;
                for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                {
                    double temp = v_i[l_comp] - v_j[l_comp];
                    v_diff_norm += temp*temp;
                }
                v_diff_norm = sqrt (v_diff_norm);
                v_diff_norm /= eps_i;
                double delta_t_j = 1.0 / (v_diff_norm /havg_i + 4.0 * nu / (havg_i * havg_i) );
//                 double delta_t_j = 1.0 / (2.0 * v_diff_norm /hmin_i + 4.0 * nu / (hmin_i * hmin_i) );
                double delta_t_j_avg_novisc =  1.0 / (2.0 * v_diff_norm /havg_i );
                if (delta_t_j < delta_t_i)
                    delta_t_i = delta_t_j;
                if (delta_t_j_avg_novisc < delta_t_i_avg_novisc)
                    delta_t_i_avg_novisc = delta_t_j_avg_novisc;
//                    if ((v_i_par >= 0.0 && v_j_par <= 0.0) || (v_i_par <= 0.0 && v_j_par >= 0.0))
//                    {
//                        double delta_t_j = CFLNumber * 1.0 / (2.0 * norm_2(v_diff) /hmin_i + 4.0 * mViscosity / (hmin_i * hmin_i));
////                        double delta_t_j = CFLNumber / ((fabs(v_i_par) + fabs(v_j_par)) / mHmin[i_node] + 2.0 * mViscosity / (mHmin[i_node] * mHmin[i_node]));
//                        // 							KRATOS_WATCH(delta_t_j);
//                        // 							KRATOS_WATCH(delta_t_i);
//                        if (delta_t_j < delta_t_i)
//                            delta_t_i = delta_t_j;
//                    }
            }
            //choose the overall minimum of delta_t_i
            if (delta_t_i < delta_t)
                delta_t = delta_t_i;
            if (delta_t_i_avg < mdelta_t_avg)
                mdelta_t_avg = delta_t_i_avg;
            if (delta_t_i_avg_novisc < delta_t_avg_novisc)
                delta_t_avg_novisc = delta_t_i_avg_novisc;
        }

        //finalizing parallel computations
        double delta_t = dt_vec[0];
        mdelta_t_avg = dt_avg_vec[0];
        double delta_t_avg_novisc = dt_avg_novisc_vec[0];
        for(unsigned int i=1; i<dt_vec.size(); i++)
        {
            if(delta_t > dt_vec[i]) delta_t = dt_vec[i];
            if(mdelta_t_avg > dt_vec[i]) mdelta_t_avg = dt_avg_vec[i];
            if(delta_t_avg_novisc > dt_vec[i]) delta_t_avg_novisc = dt_avg_novisc_vec[i];
        }



        //take into account wall law in the estimation
//         int slip_size = mSlipBoundaryList.size();
//         for (int i_slip = 0; i_slip < slip_size; i_slip++)
//         {
//             unsigned int i_node = mSlipBoundaryList[i_slip];
// 	    double nu = mViscosity[i_node];
//
// 	    double delta_t_i = 0.25*mY_wall*mY_wall/nu;
//
// 	    // Reducing wall friction for the large element near wall. Pooyan.
// 		double reducing_factor = 1.00;
//                 double h_min = mHavg[i_node];
// 		if(mY_wall < h_min)
//                      reducing_factor = mY_wall / h_min;
// 		delta_t_i /= reducing_factor;
//
// 	    if (delta_t_i < delta_t)
//                 delta_t = delta_t_i;
// 	}
// 	mdelta_t_avg = delta_t; //this should not be done ... remove it or decide what to do...
        delta_t_avg_novisc *= CFLNumber;
        //
        mnumsubsteps = ceil (delta_t_avg_novisc/delta_t);
//         mnumsubsteps += 1; //this is for security
//            delta_t *= CFLNumber;
        if (mnumsubsteps <= 1)
        {
            mnumsubsteps=1;
            delta_t_avg_novisc = delta_t;
        }
        //std::cout << "mdelta_t_avg ="  << mdelta_t_avg <<std::endl;
        //std::cout << "delta_t ="  << delta_t <<std::endl;
        //std::cout << "mnumsubsteps ="  << mnumsubsteps <<std::endl;
        delta_t = delta_t_avg_novisc;
//            delta_t *= CFLNumber;
        //*******************
        //perform MPI syncronization of the dt (minimum should be kept)
        return delta_t;
        KRATOS_CATCH ("")
    }

    void ApplySmagorinsky (double MolecularViscosity, double Cs)
    {
        if (Cs != 0)
      {
            if (TDim == 3)
                ApplySmagorinsky3D (MolecularViscosity, Cs);
	 else
                KRATOS_THROW_ERROR (std::logic_error,"smagorinsky not yet implemented in 2D","");
      }
    }

    void UpdateFixedVelocityValues()
    {
        KRATOS_TRY
        //read velocity and pressure data from Kratos
	  // ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
//         mr_matrix_container.FillVectorFromDatabase (VELOCITY, mvel_n1, rNodes);
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
        KRATOS_CATCH ("");
    }
    //**********************************************************************************
    //function to solve fluid equations - fractional step 1: compute fractional momentum
    void SolveStep1()
    {
        KRATOS_TRY
        //PREREQUISITES
        //variables for node based data handling
        ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
        int n_nodes = rNodes.size();
        //storage of nodal values in local variables
        CalcVectorType rhs;
        rhs.resize (n_nodes);
        //read velocity and pressure data from Kratos
//         mr_matrix_container.FillVectorFromDatabase (VELOCITY, mvel_n1, rNodes);
//         mr_matrix_container.FillOldVectorFromDatabase (VELOCITY, mvel_n, rNodes);
//         mr_matrix_container.FillScalarFromDatabase (VISCOSITY, mViscosity, rNodes);
//         mr_matrix_container.FillScalarFromDatabase (PRESSURE, mPn1, rNodes);
//         mr_matrix_container.FillOldScalarFromDatabase (PRESSURE, mPn, rNodes);
        mr_matrix_container.FillScalarFromDatabase (DISTANCE, mdistances, mr_model_part.Nodes() );
//         mr_matrix_container.FillScalarFromDatabase(DIAMETER, mD, mr_model_part.Nodes());
        mr_matrix_container.FillScalarFromDatabase (POROSITY, mEps, mr_model_part.Nodes() );
        mr_matrix_container.FillScalarFromDatabase (LIN_DARCY_COEF, mA, mr_model_part.Nodes() );
        mr_matrix_container.FillScalarFromDatabase (NONLIN_DARCY_COEF, mB, mr_model_part.Nodes() );
        //read time step size from Kratos
        ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
        double delta_t = CurrentProcessInfo[DELTA_TIME];
        //compute intrinsic time
         double time_inv_avg = 1.0/mdelta_t_avg;
// 	if(mmax_dt < mdelta_t_avg) mmax_dt = mdelta_t_avg;
// 	double time_inv_avg = 1.0/mmax_dt;
        double stabdt_pressure_factor  = mstabdt_pressure_factor;
        double stabdt_convection_factor  = mstabdt_convection_factor;


        //double tau2_factor = mtau2_factor;
        #pragma omp parallel for firstprivate(time_inv_avg,stabdt_pressure_factor,stabdt_convection_factor)
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            // 					double& h_i = mHavg[i_node];
            double& h_avg_i = mHavg[i_node];
	    double& h_min_i = mHmin[i_node];

            array_1d<double, TDim>& a_i = mvel_n1[i_node];
            const double nu_i = mViscosity[i_node];
            const double eps_i = mEps[i_node];
            //const double d_i = mD[i_node];
            const double lindarcy_i = mA[i_node];
            const double nonlindarcy_i = mB[i_node];
            double vel_norm = norm_2 (a_i);
            //double porosity_coefficient = ComputePorosityCoefficient(nu_i, vel_norm, eps_i, d_i);
            double porosity_coefficient = ComputePorosityCoefficient (vel_norm, eps_i, lindarcy_i, nonlindarcy_i);
            vel_norm /= eps_i;
            double tau = 1.0 / (2.0 * vel_norm / h_min_i + stabdt_pressure_factor*time_inv_avg + (4.0*nu_i) / (h_avg_i * h_avg_i) + porosity_coefficient);
            double tau_conv = 1.0 / (2.0 * vel_norm / h_min_i + stabdt_convection_factor*time_inv_avg );
            mTauPressure[i_node] = tau;
            mTauConvection[i_node] = tau_conv;
//             mTau2[i_node] = (nu_i + h_avg_i*vel_norm*0.5) *tau2_factor;
        }

//         //smoothen the tau press - mTau2 used as temp var
//         #pragma omp parallel for
//         for (int i_node = 0; i_node < n_nodes; i_node++)
//         {
//             double& tau = mTau2[i_node]; //******************
//             tau = mTauPressure[i_node];
// 	    double counter = 1.0;
//              //const double& p_i = pressure[i_node];
//             for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
//             {
//                 unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
// 		tau += mTauPressure[j_neighbour];
// 		counter+=1.0;
//             }
//             tau/=counter;
//         }
//
//         mTauPressure = mTau2;



        //calculating the convective projection
        #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            array_1d<double, TDim>& pi_i = mPi[i_node]; //******************
            //setting to zero
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                pi_i[l_comp] = 0.0;
            array_1d<double, TDim> a_i = mvel_n1[i_node];
            const array_1d<double, TDim>& U_i = mvel_n1[i_node];
            const double& eps_i = mEps[i_node];
            a_i /= eps_i;
            //const double& p_i = pressure[i_node];
            for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
            {
                unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                array_1d<double, TDim> a_j = mvel_n1[j_neighbour];
                const array_1d<double, TDim>& U_j = mvel_n1[j_neighbour];
                const double& eps_j = mEps[j_neighbour];
                a_j /= eps_j;
                CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues() [csr_index];
                edge_ij.Add_ConvectiveContribution (pi_i, a_i, U_i, a_j, U_j);
//                    edge_ij.Add_grad_p(pi_i, p_i, p_j);
            }
//              const double m_inv = mr_matrix_container.GetInvertedMass() [i_node];
//              for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
//                  pi_i[l_comp] *= m_inv;
        }

                int inout_size = mInOutBoundaryList.size();
        //#pragma omp parallel for firstprivate(slip_size)
        for (int i = 0; i < inout_size; i++)
        {
            unsigned int i_node = mInOutBoundaryList[i];
//             double dist = mdistances[i_node];
//             if (dist <= 0.0)
//             {
                const array_1d<double, TDim>& U_i = mvel_n1[i_node];
                const array_1d<double, TDim>& an_i = mInOutNormal[i_node];
                double projection_length = 0.0;
                //double Ain = 0.0;
                for (unsigned int comp = 0; comp < TDim; comp++)
                {
                    projection_length += U_i[comp] * an_i[comp];
                }

		array_1d<double, TDim>& pi_i = mPi[i_node];

		for (unsigned int comp = 0; comp < TDim; comp++)
                    pi_i[comp] += projection_length * U_i[comp] ;
//             }
        }

                #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            array_1d<double, TDim>& pi_i = mPi[i_node];

            const double m_inv = mr_matrix_container.GetInvertedMass() [i_node];
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                pi_i[l_comp] *= m_inv;

        }


//        //completing with boundary integrals
//	//loop over all faces
//	for (ModelPart::ConditionsContainerType::iterator cond_it = mr_model_part.ConditionsBegin(); cond_it != mr_model_part.ConditionsEnd(); cond_it++)
//	{
//		//get geometry data of the face
//		Geometry<Node < 3 > >& face_geometry = cond_it->GetGeometry();
//
//		//reference for area normal of the face
//		array_1d<double, 3 > & face_normal = cond_it->GetValue(NORMAL);
//		double A = norm_2(face_normal);
//
//		unsigned int i_node0 = static_cast<unsigned int> (face_geometry[0].FastGetSolutionStepValue(AUX_INDEX));
//		unsigned int i_node1 = static_cast<unsigned int> (face_geometry[1].FastGetSolutionStepValue(AUX_INDEX));
//		unsigned int i_node2 = static_cast<unsigned int> (face_geometry[2].FastGetSolutionStepValue(AUX_INDEX));
//
//		if(face_geometry[0].IsFixed(VELOCITY_X) && face_geometry[1].IsFixed(VELOCITY_X) && face_geometry[2].IsFixed(VELOCITY_X))
//		{
//
//		//KRATOS_WATCH(cond_it->Id());
//		//	    if (static_cast<bool>(cond_it->GetValue(IS_STRUCTURE)) == false)
//		//{
// 		const array_1d<double,TDim>& v_0 = mvel_n1[i_node0];
//		const array_1d<double,TDim>& v_1 = mvel_n1[i_node1];
//		const array_1d<double,TDim>& v_2 = mvel_n1[i_node2];
//		double An0 = inner_prod(v_0,face_normal) / (A*mEps[i_node0]);
//		double An1 = inner_prod(v_1,face_normal) / (A*mEps[i_node1]);
//		double An2 = inner_prod(v_2,face_normal) / (A*mEps[i_node2]);
//		//KRATOS_WATCH(face_normal);
//		mPi[i_node0] -= ((2.0*An0+An1+An2)*0.5*0.333333333333333333333333333333*0.5)*face_normal;
//		mPi[i_node1] -= ((An0+2.0*An1+An2)*0.5*0.333333333333333333333333333333*0.5)*face_normal;
//		mPi[i_node2] -= ((An0+An1+2.0*An2)*0.5*0.333333333333333333333333333333*0.5)*face_normal;
//		}
//		//}
//	}
//
//
//        //calculating the convective projection
//        #pragma omp parallel for
//        for (int i_node = 0; i_node < n_nodes; i_node++)
//        {
//            array_1d<double, TDim>& pi_i = mPi[i_node]; //******************
//            const double m_inv = mr_matrix_container.GetInvertedMass() [i_node];
//            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
//                pi_i[l_comp] *= m_inv;
//        }
//
//
//         KRATOS_WATCH("step1 before rk loop")
// 	KRATOS_WATCH(mnumsubsteps)
// 	KRATOS_WATCH(mPn)
// 	KRATOS_WATCH(mPn1)
// 	KRATOS_WATCH(mPi)
// 	KRATOS_WATCH(mvel_n1)
// 	KRATOS_WATCH(mvel_n)

#ifdef DEBUG_OUTPUT
        KRATOS_WATCH("before RK of step1 - new")
        double aux_v=0.0;
        for (int i_node = 0; i_node < mvel_n1.size(); i_node++)
            aux_v += inner_prod(mvel_n1[i_node],mvel_n1[i_node]);
        double aux_oldv=0.0;
        for (int i_node = 0; i_node < mvel_n1.size(); i_node++)
            aux_oldv += inner_prod(mvel_n[i_node],mvel_n[i_node]);
        double aux_pi=0.0;
        for (int i_node = 0; i_node < mvel_n1.size(); i_node++)
            aux_pi += inner_prod(mPi[i_node],mPi[i_node]);

        KRATOS_WATCH(inner_prod(mPn,mPn));
        KRATOS_WATCH(aux_v);
        KRATOS_WATCH(aux_oldv);
        KRATOS_WATCH(aux_pi);
        KRATOS_WATCH(inner_prod(mdistances,mdistances));
        KRATOS_WATCH(inner_prod(mViscosity,mViscosity));
#endif

        CalcVectorType auxn = mvel_n;

        double n_substeps = mnumsubsteps+1;

        double reduced_it = 0;

        double energy_initial = 0.0;
        double energy_final = 1.0;

        //compute initial kinetic energy
        #pragma omp parallel for firstprivate(n_nodes) reduction(+:energy_initial)
        for (int i_node = 0; i_node < n_nodes; i_node++)
            if (mdistances[i_node] <= 0.0)
                energy_initial += mr_matrix_container.GetLumpedMass()[i_node] * inner_prod(mvel_n[i_node],mvel_n[i_node]);

        //KRATOS_WATCH(energy_initial)
        // KRATOS_WATCH(n_substeps)

        while(reduced_it++ < 2 )
        {
        double delta_t_substep = delta_t/n_substeps;
        for (unsigned int substep = 0; substep<n_substeps; substep++)
        {
            //std::cout << "substep " << substep+1 << " of " << n_substeps << std::endl;
            mr_matrix_container.AssignVectorToVector (mvel_n, mWork); //mWork = mvel_n
            //first step of Runge Kutta
            mr_matrix_container.AssignVectorToVector (mvel_n, mvel_n1); //mvel_n1 = mvel_n
            mr_matrix_container.SetToZero (rhs);
            CalculateRHS (mvel_n1, mPn, mvel_n1, rhs,mdiag_stiffness);
            Add_Effective_Inverse_Multiply (mWork, mWork, delta_t_substep / 6.0, mr_matrix_container.GetLumpedMass(),mdiag_stiffness,rhs);
            Add_Effective_Inverse_Multiply (mvel_n1, mvel_n, 0.5 * delta_t_substep, mr_matrix_container.GetLumpedMass(),mdiag_stiffness, rhs);
            ApplyVelocityBC (mvel_n1);
            //second step
            mr_matrix_container.SetToZero (rhs);
            CalculateRHS (mvel_n1, mPn, mvel_n1, rhs,mdiag_stiffness);
            Add_Effective_Inverse_Multiply (mWork, mWork, delta_t_substep / 3.0, mr_matrix_container.GetLumpedMass(),mdiag_stiffness, rhs);
            Add_Effective_Inverse_Multiply (mvel_n1, mvel_n, 0.5 * delta_t_substep, mr_matrix_container.GetLumpedMass(),mdiag_stiffness, rhs);
            ApplyVelocityBC (mvel_n1);
            //third step
            mr_matrix_container.SetToZero (rhs);
            CalculateRHS (mvel_n1, mPn, mvel_n1, rhs,mdiag_stiffness);
            Add_Effective_Inverse_Multiply (mWork, mWork, delta_t_substep / 3.0, mr_matrix_container.GetLumpedMass(),mdiag_stiffness, rhs);
            Add_Effective_Inverse_Multiply (mvel_n1, mvel_n, delta_t_substep, mr_matrix_container.GetLumpedMass(),mdiag_stiffness, rhs);
            ApplyVelocityBC (mvel_n1);
            //fourth step
            mr_matrix_container.SetToZero (rhs);
            CalculateRHS (mvel_n1, mPn, mvel_n1, rhs,mdiag_stiffness);
            Add_Effective_Inverse_Multiply (mWork, mWork, delta_t_substep / 6.0, mr_matrix_container.GetLumpedMass(),mdiag_stiffness, rhs);
            //compute right-hand side
            mr_matrix_container.AssignVectorToVector (mWork, mvel_n1);
            ApplyVelocityBC (mvel_n1);
            //prepare for next step
            mr_matrix_container.AssignVectorToVector (mvel_n1, mvel_n);

        }

            energy_final = 0.0;
            //compute initial kinetic energy
            #pragma omp parallel for firstprivate(n_nodes) reduction(+:energy_final)
            for (int i_node = 0; i_node < n_nodes; i_node++)
                if (mdistances[i_node] <= 0.0)
                    energy_final += mr_matrix_container.GetLumpedMass()[i_node] * inner_prod(mvel_n1[i_node],mvel_n1[i_node]);

            //put back the original velocity at step n
            mr_matrix_container.AssignVectorToVector (auxn, mvel_n);

            if(energy_final < 1.5*energy_initial) break;
            else n_substeps*=10;

            if(reduced_it > 1)
            {
                KRATOS_WATCH(energy_initial)
                KRATOS_WATCH(energy_final)
                KRATOS_WATCH(n_substeps)
            }
        }


//         mr_matrix_container.WriteVectorToDatabase (VELOCITY, mvel_n1,  mr_model_part.Nodes() );
//         KRATOS_WATCH("end of step1")
// 	KRATOS_WATCH(mvel_n1)
// 	KRATOS_WATCH(mvel_n)
#ifdef DEBUG_OUTPUT
        KRATOS_WATCH("end of step1 - new")
        aux_v=0.0;
        for (int i_node = 0; i_node < mvel_n1.size(); i_node++)
            aux_v += inner_prod(mvel_n1[i_node],mvel_n1[i_node]);
        double aux_xi=0.0;
        for (int i_node = 0; i_node < mvel_n1.size(); i_node++)
            aux_xi += inner_prod(mXi[i_node],mXi[i_node]);

        KRATOS_WATCH(inner_prod(mPn,mPn));
        KRATOS_WATCH(inner_prod(mdistances,mdistances));
        KRATOS_WATCH(inner_prod(mViscosity,mViscosity));

        KRATOS_WATCH(aux_v);
        KRATOS_WATCH(aux_xi);
#endif
        KRATOS_CATCH ("")
    }
    //*********************************************************************
    //function to calculate right-hand side of fractional momentum equation
    void CalculateRHS (
        const CalcVectorType& vel,
        const ValuesVectorType& pressure,
        const CalcVectorType& convective_velocity,
        CalcVectorType& rhs,
	ValuesVectorType& diag_stiffness)
    {
        KRATOS_TRY
        int n_nodes = vel.size();
        //perform MPI syncronization
        //calculating the RHS
        array_1d<double, TDim> stab_low;
        array_1d<double, TDim> stab_high;
        double inverse_rho = 1.0 / mRho;
        #pragma omp parallel for private(stab_low,stab_high)
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            double dist = mdistances[i_node];
            if (dist <= 0.0) //node is inside domain ---- if outside do nothing
            {
                const double nu_i = mViscosity[i_node];
                const double nu_j = nu_i;
                array_1d<double, TDim>& rhs_i = rhs[i_node];
                const array_1d<double, TDim>& f_i = mBodyForce;
                array_1d<double, TDim> a_i = convective_velocity[i_node];
                const array_1d<double, TDim>& U_i = vel[i_node];
                const array_1d<double, TDim>& pi_i = mPi[i_node];
                const double& p_i = pressure[i_node];
                const double& eps_i = mEps[i_node];
                const double lindarcy_i = mA[i_node];
                const double nonlindarcy_i = mB[i_node];
                double edge_tau = mTauConvection[i_node];
                a_i /= eps_i;
                //initializing with the external forces (e.g. gravity)
                double& m_i = mr_matrix_container.GetLumpedMass() [i_node];
                for (unsigned int comp = 0; comp < TDim; comp++)
                    rhs_i[comp] = m_i * eps_i * f_i[comp] ;
                //applying the effect of the porosity
                double porosity_coefficient = ComputePorosityCoefficient ( norm_2 (U_i), eps_i, lindarcy_i, nonlindarcy_i);
                diag_stiffness[i_node]= m_i * porosity_coefficient;

                //std::cout << i_node << "rhs =" << rhs_i << "after adding body force" << std::endl;
                //convective term
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                    array_1d<double, TDim> a_j = convective_velocity[j_neighbour];
                    const array_1d<double, TDim>& U_j = vel[j_neighbour];
                    const array_1d<double, TDim>& pi_j = mPi[j_neighbour];
                    const double& p_j = pressure[j_neighbour];
                    const double& eps_j = mEps[j_neighbour];
//                             const double& beta_j = mBeta[j_neighbour];
                    a_j /= eps_j;
                    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues() [csr_index];
                    edge_ij.Sub_ConvectiveContribution (rhs_i, a_i, U_i, a_j, U_j);
                    //std::cout << i_node << "rhs =" << rhs_i << "after convective contrib" << std::endl;
                    //take care! we miss including a B.C.  for the external pressure
                    //edge_ij.Add_Gp (rhs_i,p_i*inverse_rho,p_j*inverse_rho);
                    edge_ij.Sub_grad_p(rhs_i, p_i*inverse_rho*eps_i, p_j * inverse_rho*eps_i);
                    //                        edge_ij.Add_grad_p(rhs_i, p_i*inverse_rho, p_j * inverse_rho);
                    //std::cout << i_node << "rhs =" << rhs_i << "after Gp" << std::endl;
                    edge_ij.Sub_ViscousContribution (rhs_i, U_i, nu_i, U_j, nu_j);
// 		    edge_ij.Add_ViscousContribution(rhs_i, U_i, nu_i, U_j, nu_j);
                    //std::cout << i_node << "rhs =" << rhs_i << "after viscous" << std::endl;
                    //add stabilization
                    edge_ij.CalculateConvectionStabilization_LOW (stab_low, a_i, U_i, a_j, U_j);
//                            edge_ij.CalculateConvectionStabilization_LOW(stab_low, a_i, U_i,p_i, a_j, U_j,p_j);
                    edge_ij.CalculateConvectionStabilization_HIGH (stab_high, a_i, pi_i, a_j, pi_j);
//                            double beta = 1.0;
//                             double beta = beta_i;
//                             if(beta_j > beta)
//                                 beta = beta_j;
//                            beta = 1.0;
//                            edge_ij.Sub_StabContribution(rhs_i, edge_tau*beta, 1.0, stab_low, stab_high);
//                             edge_ij.Sub_StabContribution(rhs_i, edge_tau, (1.0-beta), stab_low, stab_high);
                    edge_ij.Sub_StabContribution (rhs_i, edge_tau, 1.0, stab_low, stab_high);
                }
                //                                                std::cout << i_node << "rhs =" << rhs_i << std::endl;
            }
        }

        int inout_size = mInOutBoundaryList.size();
        //#pragma omp parallel for firstprivate(slip_size)
        for (int i = 0; i < inout_size; i++)
        {
            unsigned int i_node = mInOutBoundaryList[i];
//             double dist = mdistances[i_node];
//             if (dist <= 0.0)
//             {
                const array_1d<double, TDim>& U_i = mvel_n1[i_node];
                const array_1d<double, TDim>& an_i = mInOutNormal[i_node];
                double projection_length = 0.0;
		double Ain = 0.0;
                for (unsigned int comp = 0; comp < TDim; comp++)
                {
                    projection_length += U_i[comp] * an_i[comp];
		    Ain += an_i[comp]*an_i[comp];
                }

		array_1d<double, TDim>& rhs_i = rhs[i_node];

		for (unsigned int comp = 0; comp < TDim; comp++)
                    rhs_i[comp] += projection_length * U_i[comp] ;
//             }
        }

        /*        		for (int i = 0; i < mSlipBoundaryList.size(); i++)
                {
        			int i_node = mSlipBoundaryList[i];
        			double dist = mdistances[i_node];
                    if (dist <= 0.0) //node is inside domain ---- if outside do nothing
                    {
						const double& p_i = pressure[i_node];
        				const array_1d<double,3>& Ani = mSlipNormal[i_node];
                        array_1d<double, TDim>& rhs_i = rhs[i_node];
        				array_1d<double, TDim> temp;
        				for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
        					temp[l_comp] = 0.0;
                         for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                        {
                            unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
        					if(mdistances[j_neighbour] <= 0.0 && mis_slip[j_neighbour] == true)
        					{
        						//const double& p_j = pressure[j_neighbour];
        						array_1d<double,3> Anj = mSlipNormal[j_neighbour];
								Anj /= norm_2(Anj);
        						for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
        							temp[l_comp] += p_i*Anj[l_comp];
        					}
        				}
        				//take out part in the direction of Ani
        				double Ai = norm_2(Ani);
        				double aux = 0.0;
        				for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
        					aux += temp[l_comp]*Ani[l_comp];
        				for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
        					temp[l_comp] -= aux *Ani[l_comp] / (Ai*Ai);
        				for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
        					rhs_i[l_comp] -= 0.25*Ai*temp[l_comp];
        			}
        		}*/
//       KRATOS_WATCH("finished**************************************************")		*/

        /*
	//correction to the pressure graient
	                //loop over all faces
        CalcVectorType press_correction(vel.size());
	  mr_matrix_container.SetToZero(press_correction);
        //        mr_matrix_container.SetToZero(slip_area);
	for (ModelPart::ConditionsContainerType::iterator cond_it = mr_model_part.ConditionsBegin(); cond_it != mr_model_part.ConditionsEnd(); cond_it++)
	{
		//get geometry data of the face
		Geometry<Node < 3 > >& face_geometry = cond_it->GetGeometry();

		//reference for area normal of the face
		array_1d<double, 3 > & face_normal = cond_it->GetValue(NORMAL);
		double A = norm_2(face_normal);

		unsigned int i_node0 = static_cast<unsigned int> (face_geometry[0].FastGetSolutionStepValue(AUX_INDEX));
		unsigned int i_node1 = static_cast<unsigned int> (face_geometry[1].FastGetSolutionStepValue(AUX_INDEX));
		unsigned int i_node2 = static_cast<unsigned int> (face_geometry[2].FastGetSolutionStepValue(AUX_INDEX));

	    if (static_cast<bool>(cond_it->GetValue(IS_STRUCTURE)) == true)
	    {

		const double& p_0 = pressure[i_node0];
		const double& p_1 = pressure[i_node1];
        		const double& p_2 = pressure[i_node2];

		//TODO: we should only keep the part orthogonal to the external normal on each node!!!!
		press_correction[i_node0] -= ((2.0*p_0+p_1+p_2)*0.5*0.333333333333333333333333333333*0.5*inverse_rho)*face_normal;
		press_correction[i_node1] -= ((p_0+2.0*p_1+p_2)*0.5*0.333333333333333333333333333333*0.5*inverse_rho)*face_normal;
		press_correction[i_node2] -= ((p_0+p_1+2.0*p_2)*0.5*0.333333333333333333333333333333*0.5*inverse_rho)*face_normal;

	    }
	    else
	    {

		const array_1d<double,TDim>& v_0 = vel[i_node0];
		const array_1d<double,TDim>& v_1 = vel[i_node1];
        		const array_1d<double,TDim>& v_2 = vel[i_node2];
		double An0 = inner_prod(v_0,face_normal) / (A*A);
		double An1 = inner_prod(v_1,face_normal) / (A*A);
		double An2 = inner_prod(v_2,face_normal) / (A*A);

		rhs[i_node0] -= ((2.0*An0+An1+An2)*0.5*0.333333333333333333333333333333*0.5)*face_normal;
		rhs[i_node1] -= ((An0+2.0*An1+An2)*0.5*0.333333333333333333333333333333*0.5)*face_normal;
		rhs[i_node2] -= ((An0+An1+2.0*An2)*0.5*0.333333333333333333333333333333*0.5)*face_normal;

	    }
	}

	//slip condition
	int slip_size = mSlipBoundaryList.size();
	#pragma omp parallel for firstprivate(slip_size)
	for (int i_slip = 0; i_slip < slip_size; i_slip++)
	{
	    unsigned int i_node = mSlipBoundaryList[i_slip];
	    double dist = mdistances[i_node];
	    if (dist <= 0.0 && mis_slip[i_node] == true)
	    {
		array_1d<double, TDim>& rhs_i = rhs[i_node];
        // 		array_1d<double, TDim>& an_i = mSlipNormal[i_node];
        // 		double normalization = 0.0;
        // 		for (unsigned int comp = 0; comp < TDim; comp++)
        // 		{
        // 		    normalization += an_i[comp] * an_i[comp];
        // 		}
        // 		normalization = sqrt(normalization);
		array_1d<double,TDim>& press_corr_i = press_correction[i_node];
		for (unsigned int comp = 0; comp < TDim; comp++)
 		    rhs_i[comp] +=  press_corr_i[comp];

		//we should remove here the normal component!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    }
	}


        */


        //apply wall resistance
        if (mWallLawIsActive == true)
            ComputeWallResistance (vel,diag_stiffness);
//         ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
//         mr_matrix_container.WriteVectorToDatabase (VELOCITY, mvel_n1, rNodes);
        KRATOS_CATCH ("")
    }
    //*************************************************************************
    //function to solve fluid equations - fractional step 2: calculate pressure
    int SolveStep2 (typename TLinearSolver::Pointer pLinearSolver)
    {
        KRATOS_TRY
//         typedef Node < 3 > PointType;
//         typedef PointerVector<PointType > PointVector;
//         typedef PointVector::iterator PointIterator;
        #pragma omp parallel for
        for ( int i_node = 0; i_node < static_cast<int>(mr_model_part.Nodes().size()); i_node++)
            mis_visited[i_node] = 0;

        int layer_counter = -1;
        boost::numeric::ublas::vector<int> layers(mr_model_part.Nodes().size());
        boost::numeric::ublas::vector<int> layer_limits(3);

        //Re-generate a container with LAYER 0 and LAYER 1 after convection of the free surface
        layer_limits[0] = 0;
		#pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(mr_model_part.Nodes().size()); i_node++)
        {
            if(mdistances[i_node] < 0.0)
            {
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    //get global index of neighbouring node j
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    if(mdistances[j_neighbour] >= 0.0 && mis_visited[i_node] == 0 )
                    {
						#pragma omp critical 
                        layers[++layer_counter] = i_node;
						
                        mis_visited[i_node] = 1;
						break;
                        }
                    }
                }
            else
                mPn1[i_node] = 0.0;
        }
        layer_limits[1] = layer_counter;

        for(unsigned int i=0; i<static_cast<unsigned int>(layer_limits[1]); i++)
        {
            unsigned int i_node = layers[i];

            for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
            {
                unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                if( mdistances[j_neighbour] >= 0.0 && mis_visited[j_neighbour] == 0)
                {
                    layers[layer_counter++] = j_neighbour;
                    mis_visited[j_neighbour] = 2;
                }
            }
        }
        layer_limits[2] = layer_counter;

        int return_value = 0;

        //on the first layer outside the pressure is set to a value such that on the free surface the pressure is approx 0
		#pragma omp parallel for
        for( int iii=static_cast<int>(layer_limits[1]); iii<static_cast<int>(layer_limits[2]); iii++)
        {
            unsigned int i_node = layers[iii];
            array_1d<double, TDim> grad_d;
            for (unsigned int comp = 0; comp < TDim; comp++)
                grad_d[comp] = 0.0;
            double dist_i = mdistances[i_node];
            for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
            {
                //get global index of neighbouring node j
                unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                const double& dist_j  = mdistances[j_neighbour];
                //projection of pressure gradients
                CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues() [csr_index];
                edge_ij.Add_grad_p (grad_d, dist_i, dist_j);
            }
            const double& m_inv = mr_matrix_container.GetInvertedMass() [i_node];
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                grad_d[l_comp] *= m_inv;
            double norm_grad = norm_2 (grad_d);
            if (norm_grad < 2.0)
            {
                if(dist_i < 0.01*mHavg[i_node] )
                    dist_i = 0.0;
                else if(dist_i > 2.0*mHavg[i_node] )
                {
                    KRATOS_WATCH("distance is much larger than expected!!")
                    dist_i = 2.0*mHavg[i_node];
                }

                if(norm_grad > 0.001)
                {
                grad_d /= norm_grad; //this is the direction of the gradient of the distances
                grad_d *= dist_i; //this is the vector with the distance of node_i from the closest point on the free surface
                }
                else
                {
                    KRATOS_WATCH("norm grad is very small!!!!")
                    grad_d *= 0.0;
                }


                const array_1d<double, TDim>& press_grad = mXi[i_node]; //iii->FastGetSolutionStepValue (PRESS_PROJ);
                double pestimate = inner_prod (press_grad,grad_d);
                mPn1[i_node] = pestimate;
// 			    KRATOS_WATCH("peastimate step2")
// 			    KRATOS_WATCH(iii->Id())
// 			    KRATOS_WATCH(grad_d)
// 			    KRATOS_WATCH(press_grad)
// 			    KRATOS_WATCH(pestimate)
            }
            else
            {
                std::cout << "attention gradient of distance much greater than 1 on node:" << i_node  <<std::endl;
                return_value = -1;
//                 return -1;
                double avg_number = 0.0;
                double pavg = 0.0;
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    if( mis_visited[j_neighbour] == 1)
                    {
                        pavg += mPn1[j_neighbour];
                        avg_number += 1.0;
                    }
                }
                if (avg_number == 0)
                    KRATOS_THROW_ERROR (std::logic_error,"can not happen that the extrapolation node has no neighbours","");
                mPn1[i_node] = pavg/avg_number;
            }
        }
        //if a node is very close to the free surface (relatively to the element size) fix the pressure on it
//            for(ModelPart::NodesContainerType::iterator iii = mr_model_part.NodesBegin(); iii!=mr_model_part.NodesEnd(); iii++)
//            {
//                unsigned int i_node = iii->FastGetSolutionStepValue(AUX_INDEX);
//
//                double dist = mdistances[i_node];
//                if(dist > 0.0 && dist < 0.01*mHavg[i_node])
//                    iii->FastGetSolutionStepValue(PRESSURE) = 0.0;
//
//            }
        //PREREQUISITES
        //allocate memory for variables
        ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
        int n_nodes = rNodes.size();
        //unknown and right-hand side vector
        TSystemVectorType dp, rhs;
        dp.resize (n_nodes);
        rhs.resize (n_nodes);
        array_1d<double, TDim> dU_i, dU_j, work_array;
        //read time step size from Kratos
        ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
        double delta_t = CurrentProcessInfo[DELTA_TIME];
#ifdef _OPENMP
//        double time_inv = 0.0; //1.0/delta_t;
        //read the pressure projection from the database
#endif
//         mr_matrix_container.FillOldScalarFromDatabase (PRESSURE, mPn, mr_model_part.Nodes() );
//         mr_matrix_container.FillScalarFromDatabase (PRESSURE, mPn1, mr_model_part.Nodes() );
//         mr_matrix_container.FillVectorFromDatabase (PRESS_PROJ, mXi, rNodes);
//         mr_matrix_container.FillVectorFromDatabase (VELOCITY, mvel_n1, rNodes);
        //for (int i_node = 0; i_node < n_nodes; i_node++)
        //    std::cout << mvel_n1[i_node] << std::endl;
        //loop over all nodes
//            double rho_inv = 1.0 / mRho;
        #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            double& rhs_i = rhs[i_node];
            rhs_i = 0.0;
            const double& p_i = mPn1[i_node];
            const double& p_old_i = mPn[i_node];
            const array_1d<double, TDim>& U_i_curr = mvel_n1[i_node];
//                const double& eps_i = mEps[i_node];
            array_1d<double, TDim>& xi_i = mXi[i_node];
            double l_ii = 0.0;
//                double div_i = 0.0;
            //loop over all neighbours
            for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
            {
                unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                const double& p_j = mPn1[j_neighbour];
                const double& p_old_j = mPn[j_neighbour];
                const array_1d<double, TDim>& U_j_curr = mvel_n1[j_neighbour];
                const array_1d<double, TDim>& xi_j = mXi[j_neighbour];
//                    const double& eps_j = mEps[j_neighbour];
                CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues() [csr_index];
#ifdef SYMM_PRESS
                double edge_tau = 0.25* (mTauPressure[i_node] + mTauPressure[j_neighbour]);
#else
                double edge_tau = 0.5*mTauPressure[i_node];
#endif
                //     						double edge_tau = CalculateEdgeTau(time_inv,h_i,a_i,h_j,a_j);
                //
                if (edge_tau < delta_t) edge_tau=delta_t;
                //compute laplacian operator
                double sum_l_ikjk;
                edge_ij.CalculateScalarLaplacian (sum_l_ikjk);
//                    double sum_l_ikjk_onlystab = sum_l_ikjk * (edge_tau);
                double sum_l_ikjk_onlydt = sum_l_ikjk * (delta_t);
                sum_l_ikjk *= (delta_t + edge_tau);
                //assemble right-hand side
                //pressure contribution
//                    rhs_i -= sum_l_ikjk_onlystab * (p_j - p_i);
                rhs_i -= sum_l_ikjk * (p_j - p_i);
                rhs_i += sum_l_ikjk_onlydt * (p_old_j - p_old_i);
                //calculating the divergence of the fract vel
//                    edge_ij.Sub_D_v(div_i, U_i_curr*mRho*eps_i, U_j_curr * mRho*eps_j);
                edge_ij.Sub_D_v (rhs_i, U_i_curr*mRho, U_j_curr * mRho);
                //       						edge_ij.Sub_D_v(rhs_i,a_i*rho_i,a_j*rho_i);
                //high order stabilizing term
                double temp = 0.0;
                // 						edge_ij.Add_div_v(temp,mTauPressure[i_node]*xi_i,mTauPressure[j_neighbour]*xi_j);
                edge_ij.Add_div_v (temp, xi_i, xi_j);
                rhs_i += edge_tau * temp;
                //assemble laplacian matrix
                mL (i_node, j_neighbour) = sum_l_ikjk;
                l_ii -= sum_l_ikjk;
            }
//                //area correction to prevent mass loss
//                rhs_i -= mdiv_error[i_node];
//                rhs_i += div_i * eps_i;
            mL (i_node, i_node) = l_ii;
        }
        if (muse_mass_correction == true)
        {
            #pragma omp parallel for
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                double& rhs_i = rhs[i_node];
                rhs_i -= mdiv_error[i_node];
            }
        }
        //find the max diagonal term
        double max_diag = 0.0;
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            double L_diag = mL (i_node, i_node);
            if (fabs (L_diag) > fabs (max_diag) ) max_diag = L_diag;
        }
        max_diag *= 1e10;
//         if (max_diag < 1e20) max_diag=1e20;
        //respect pressure boundary conditions by penalization
//            double huge = max_diag * 1e6;
//            for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++) {
//                unsigned int i_node = mPressureOutletList[i_pressure];
//                mL(i_node, i_node) = huge;
//                rhs[i_node] = 0.0;
//            }
        for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++)
        {
            unsigned int i_node = mPressureOutletList[i_pressure];
            mL (i_node, i_node) = max_diag;
            rhs[i_node] = 0.0;
            for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
            {
                unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                mL (i_node, j_neighbour) = 0.0;
            }
        }
        //modification for level_set
        //				mr_matrix_container.FillScalarFromDatabase(DISTANCE, mdistances, mr_model_part.Nodes());
        //				for (unsigned int i_dist = 0; i_dist < mdistances.size(); i_dist++)
        //				{
        //					if(mdistances[i_dist] >= 0)
        //					{
        //						mL(i_dist, i_dist) = huge;
        //						rhs[i_dist] = 0.0;
        //					}
        //				}
        #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            if (mdistances[i_node] >= 0)
            {
                mL (i_node, i_node) = max_diag;
                rhs[i_node] = 0.0;
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                    mL (i_node, j_neighbour) = 0.0;
                }
            }
            else
            {
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                    if (mdistances[j_neighbour] >= 0)
                        mL (i_node, j_neighbour) = 0.0;
                }
            }
        }
//	    for (int i_node = 0; i_node < n_nodes; i_node++)
//	    {
//	      if(  fabs(mL(i_node, i_node)) < 1e-20)
//	      {
//		mL(i_node, i_node)=max_diag;
//		rhs[i_node] = 0.0;
//		KRATOS_WATCH("arghhhhhhhhhhhhhhhhhhhhhhhhhhhhhh");
//	      }
//	    }
//compute row scaling factors
        TSystemVectorType scaling_factors (n_nodes);
        double* Lvalues = mL.value_data().begin();
        SizeType* Lrow_indices = mL.index1_data().begin();
        SizeType* Lcol_indices = mL.index2_data().begin();
        #pragma omp parallel for
        for (int k = 0; k < static_cast< int> (mL.size1() ); k++)
        {
            double t = 0.0;
            SizeType col_begin = Lrow_indices[k];
            SizeType col_end = Lrow_indices[k+1];
            for (SizeType j=col_begin; j<col_end; j++)
                if ( static_cast<int> (Lcol_indices[j]) == k)
                {
                    t = fabs (Lvalues[j]);
                    break;
                }
//                        t += Lvalues[j]*Lvalues[j];
//                t = sqrt(t);
            scaling_factors[k] = 1.0/sqrt (t);
        }
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int> (mL.size1() ); k++)
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
        //set starting vector for iterative solvers
        #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; i_node++)
            dp[i_node] = 0.0;
        //KRATOS_WATCH(rhs);
        //solve linear equation system L dp = rhs
        pLinearSolver->Solve (mL, dp, rhs);
        //KRATOS_WATCH(*pLinearSolver)
        //update pressure
        #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; i_node++)
            mPn1[i_node] += dp[i_node]*scaling_factors[i_node];
        //				for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++)
        //				{
        //					unsigned int i_node = mPressureOutletList[i_pressure];
        //					mPn1[i_node] = mPressureOutlet[i_pressure];
        //				}
        //write pressure and density to Kratos
        mr_matrix_container.WriteScalarToDatabase (PRESSURE, mPn1, rNodes);
        //compute pressure proj for the next step
        #pragma omp parallel for  private(work_array)
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            array_1d<double, TDim>& xi_i = mXi[i_node];
            for (unsigned int comp = 0; comp < TDim; comp++)
                xi_i[comp] = 0.0;
            double dist = mdistances[i_node];
            if (dist <= 0.0) //node is inside domain ---- if outside do nothing
            {
                const double& p_i = mPn1[i_node];
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
                {
                    //get global index of neighbouring node j
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                    const double& p_j = mPn1[j_neighbour];
                    //projection of pressure gradients
                    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues() [csr_index];
                    edge_ij.Add_grad_p (xi_i, p_i, p_j);
                }
                const double& m_inv = mr_matrix_container.GetInvertedMass() [i_node];
                for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                    xi_i[l_comp] *= m_inv;
            }
        }
        mr_matrix_container.WriteVectorToDatabase (PRESS_PROJ, mXi, rNodes);

// 	        KRATOS_WATCH("end of step2")
// 	KRATOS_WATCH(mPn)
// 	KRATOS_WATCH(mPn1)
// 	KRATOS_WATCH(mXi)
#ifdef DEBUG_OUTPUT

        KRATOS_WATCH("end of step2 - new")
        double aux_v=0.0;
        for (int i_node = 0; i_node < mvel_n1.size(); i_node++)
            aux_v += inner_prod(mvel_n1[i_node],mvel_n1[i_node]);
        double aux_xi=0.0;
        for (int i_node = 0; i_node < mvel_n1.size(); i_node++)
            aux_xi += inner_prod(mXi[i_node],mXi[i_node]);

        KRATOS_WATCH(inner_prod(mPn1,mPn1));
        KRATOS_WATCH(aux_v);
        KRATOS_WATCH(aux_xi);
#endif

        return return_value;
        KRATOS_CATCH ("")
    }
    //**********************************************************************************
    //function to solve fluid equations - fractional step 3: correct fractional momentum
    void SolveStep3()
    {
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
        if (massume_constant_dp == true)
            factor = 1.0;
        //compute end of step momentum
        double rho_inv = 1.0 / mRho;
        #pragma omp parallel for private(correction) firstprivate(delta_t,rho_inv,factor)
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            double dist = mdistances[i_node];
            if (dist < 0.0) //node is inside domain ---- if outside do nothing
            {
                array_1d<double, TDim>& U_i_curr = mvel_n1[i_node];
                double delta_p_i = (mPn1[i_node] - mPn[i_node]) * rho_inv*factor;
//                const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];
                //setting to zero
                for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                    correction[l_comp] = 0.0;
                //compute edge contributions dt*M^(-1)Gp
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                    double delta_p_j = (mPn1[j_neighbour] - mPn[j_neighbour]) * rho_inv*factor;
                    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues() [csr_index];
                    // 							edge_ij.Sub_grad_p(correction,delta_p_i,delta_p_j);
                    edge_ij.Sub_grad_p(correction, delta_p_i, delta_p_j);
                    //                        edge_ij.Add_grad_p(correction, delta_p_i, delta_p_j);
                    //edge_ij.Add_Gp (correction,delta_p_i,delta_p_j);
                    // 							edge_ij.Sub_Gp(correction,delta_p_i,delta_p_j);
                }
                //compute prefactor
//                 double coefficient = delta_t * m_inv;
                const double m = mr_matrix_container.GetLumpedMass() [i_node];
                const double&  d = mdiag_stiffness[i_node];

                //correct fractional momentum
                for (unsigned int comp = 0; comp < TDim; comp++)
		{
                    U_i_curr[comp] += delta_t / (m + delta_t*d) * correction[comp];
                }
        }
        }

//         //imit acceleration
// 	#pragma omp parallel for
//         for(int i_node = 0; i_node < n_nodes; i_node++)
// 	{
// 	      array_1d<double,TDim>& acc = macc[i_node];
// 	      array_1d<double,TDim>& v1  = mvel_n1[i_node];
// 	      array_1d<double,TDim>& v   = mvel_n[i_node];
//
// 	      for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
// 		  acc[l_comp] = (v1[l_comp] - v[l_comp])/delta_t;
//
// 	      //limit accelerations to a maximum=100m/s/2
// 	      const double max_acc = 200;
// 	      double acc_norm = norm_2(acc);
// 	      if(acc_norm > max_acc)
// 	      {
// 		std::cout << "########################### acc norm " << acc_norm <<std::endl;
//
// 		  acc *= max_acc/acc_norm;
// 		  v1 = v;
// 		  v1 += delta_t*acc;
// 	      }
// 	}


        ApplyVelocityBC (mvel_n1);

        //save acceleration
        #pragma omp parallel for
        for(int i_node = 0; i_node < n_nodes; i_node++)
        {
            array_1d<double,TDim>& acc = macc[i_node];
            array_1d<double,TDim>& v1  = mvel_n1[i_node];
            array_1d<double,TDim>& v   = mvel_n[i_node];

            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                acc[l_comp] = (v1[l_comp] - v[l_comp])/delta_t;

        }

        //write velocity of time step n+1 to Kratos
        //calculate the error on the divergence
        if (muse_mass_correction == true)
        {
            #pragma omp parallel for private(correction) firstprivate(delta_t,rho_inv)
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                const double dist = mdistances[i_node];
                double& div_i_err = mdiv_error[i_node];
                div_i_err = 0.0;
                if (dist < 0.0) //node is inside domain ---- if outside do nothing
                {
                    const array_1d<double, TDim>& U_i_curr = mvel_n1[i_node];
                    //compute edge contributions dt*M^(-1)Gp
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                        array_1d<double, TDim>& U_j_curr = mvel_n1[j_neighbour];
                        CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues() [csr_index];
                        edge_ij.Add_D_v (div_i_err, U_i_curr*mRho, U_j_curr * mRho);
                    }
                }
            }
        }


#ifdef DEBUG_OUTPUT

        KRATOS_WATCH("end of step 3")
        double aux=0.0;
        for (int i_node = 0; i_node < n_nodes; i_node++)
            aux += inner_prod(mvel_n1[i_node],mvel_n1[i_node]);

        KRATOS_WATCH(inner_prod(mPn1,mPn1));
        KRATOS_WATCH(aux);
#endif


        mr_matrix_container.WriteVectorToDatabase (VELOCITY, mvel_n1, rNodes);


        KRATOS_CATCH ("")
    }
    void ApplyDistanceBC()
    {
        KRATOS_TRY
        //slip condition
        int size = mDistanceBoundaryList.size();
        #pragma omp parallel for firstprivate(size)
        for (int i_dist = 0; i_dist < size; i_dist++)
        {
            unsigned int i_node = mDistanceBoundaryList[i_dist];
            double& dist = mdistances[i_node];
	    dist = mDistanceValuesList[i_dist];
        }
        //fix the distance if velocity goes inwards
//         int slip_size = mSlipBoundaryList.size();
//         #pragma omp parallel for firstprivate(slip_size)
//         for (int i_slip = 0; i_slip < slip_size; i_slip++)
//             {
//             unsigned int i_node = mSlipBoundaryList[i_slip];
//             double dist = mphi_n[i_node];
// //             if(dist > 0.0)
// //             {
//             array_1d<double, TDim>& U_i = mvel_n1[i_node];
//             array_1d<double, TDim>& an_i = mSlipNormal[i_node];
//             double projection_length = 0.0;
//             double normalization = 0.0;
//             for (unsigned int comp = 0; comp < TDim; comp++)
//             {
//                 projection_length += U_i[comp] * an_i[comp];
//              }
//             if(projection_length > 0.0)
//                 dist = mphi_n[i_node];
// //              }
//         }
        KRATOS_CATCH ("")
    }
    //************************************
    void ApplyVelocityBC (CalcVectorType& VelArray)
    {
        KRATOS_TRY
//         if(mWallLawIsActive == false)
//         {
// 			std::cout << "applying corners condition" << std::endl;
//             apply conditions on corner edges
//             int edge_size = medge_nodes_direction.size();
//             #pragma omp parallel for firstprivate(edge_size)
//             for (int i = 0; i < edge_size; i++)
//             {
//                 int i_node = medge_nodes[i];
//                 const array_1d<double, TDim>& direction = medge_nodes_direction[i];
//                 double dist = mdistances[i_node];
//
//                 if(dist <= 0.0)
//                 {
//                      array_1d<double, TDim>& U_i = VelArray[i_node];
// // 		     for (unsigned int comp = 0; comp < TDim; comp++)
// //                         U_i[comp] = 0.0;
//
//                     double temp=0.0;
//                     for (unsigned int comp = 0; comp < TDim; comp++)
//                         temp += U_i[comp] * direction[comp];
//
//                     for (unsigned int comp = 0; comp < TDim; comp++)
//                         U_i[comp] = direction[comp]*temp;
//         }
//             }
// // //
//             //apply conditions on corners
//             int corner_size = mcorner_nodes.size();
//             for (int i = 0; i < corner_size; i++)
//             {
//                 int i_node = mcorner_nodes[i];
//
//                 array_1d<double, TDim>& U_i = VelArray[i_node];
//                 for (unsigned int comp = 0; comp < TDim; comp++)
//                     U_i[comp] = 0.0;
//             }



//             //apply conditions on corners
             int corner_size = mcorner_nodes.size();
             for (int i = 0; i < corner_size; i++)
             {
                 int i_node = mcorner_nodes[i];

                 array_1d<double, TDim>& U_i = VelArray[i_node];


// 		 if(mdistances[i_node] <= 0.0)
// 		  {
			array_1d<double, TDim> aux;
			for (unsigned int comp = 0; comp < TDim; comp++)
			    aux[comp] = 0.0;

			double counter = 0.0;
		      for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
		      {
			  //get global index of neighbouring node j
			  unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
			  const double& dist_j  = mdistances[j_neighbour];
			  array_1d<double, TDim>& vj = VelArray[j_neighbour];

			  if(dist_j <=  0 && mis_slip[j_neighbour] == false)
			  {
			    counter += 1.0;
			    for (unsigned int comp = 0; comp < TDim; comp++)
				aux[comp] += vj[comp];
			    }

		      }

            if(counter != 0.0)
			  for (unsigned int comp = 0; comp < TDim; comp++)
				U_i[comp] = aux[comp]/counter;
// 		  }
             }


//         }
        //slip condition
        int slip_size = mSlipBoundaryList.size();
        #pragma omp parallel for firstprivate(slip_size)
        for (int i_slip = 0; i_slip < slip_size; i_slip++)
        {
            unsigned int i_node = mSlipBoundaryList[i_slip];
            double dist = mdistances[i_node];
            if (dist <= 0.0)
            {
                array_1d<double, TDim>& U_i = VelArray[i_node];
                array_1d<double, TDim>& an_i = mSlipNormal[i_node];
                double projection_length = 0.0;
                double normalization = 0.0;
                for (unsigned int comp = 0; comp < TDim; comp++)
                {
                    projection_length += U_i[comp] * an_i[comp];
                    normalization += an_i[comp] * an_i[comp];
                }
                projection_length /= normalization;
                //tangential momentum as difference between original and normal momentum
                for (unsigned int comp = 0; comp < TDim; comp++)
                    U_i[comp] -= projection_length * an_i[comp];
            }
        }

//                 //loop over all faces
//         ValuesVectorType vel_correction(VelArray.size());
// //        CalcVectorType slip_area(VelArray.size());
//         int iterations = 10;
//         for(unsigned int i=0;i<iterations; i++)
//         {
//             mr_matrix_container.SetToZero(vel_correction);
//     //        mr_matrix_container.SetToZero(slip_area);
//             for (ModelPart::ConditionsContainerType::iterator cond_it = mr_model_part.ConditionsBegin(); cond_it != mr_model_part.ConditionsEnd(); cond_it++)
//             {
//                 if (static_cast<bool>(cond_it->GetValue(IS_STRUCTURE)) == true)
//                 {
//                     //get geometry data of the face
//                     Geometry<Node < 3 > >& face_geometry = cond_it->GetGeometry();
//
//                     //reference for area normal of the face
//                     array_1d<double, 3 > & face_normal = cond_it->GetValue(NORMAL);
//                     double n_area = norm_2(face_normal) / static_cast<double>(TDim);
//
//                     unsigned int i_node0 = static_cast<unsigned int> (face_geometry[0].FastGetSolutionStepValue(AUX_INDEX));
//                     unsigned int i_node1 = static_cast<unsigned int> (face_geometry[1].FastGetSolutionStepValue(AUX_INDEX));
//                     unsigned int i_node2 = static_cast<unsigned int> (face_geometry[2].FastGetSolutionStepValue(AUX_INDEX));
//
//                     const array_1d<double, TDim>& U_0 = VelArray[i_node0];
//                     const array_1d<double, TDim>& U_1 = VelArray[i_node1];
//                     const array_1d<double, TDim>& U_2 = VelArray[i_node2];
//
// 		    double vn0=0.0;
// 		    double vn1=0.0;
// 		    double vn2=0.0;
//                     if(mdistances[i_node0] <= 0 && face_geometry[0].IsFixed(VELOCITY_X) == false) vn0 = inner_prod(U_0,face_normal);
//                     if(mdistances[i_node1] <= 0 && face_geometry[1].IsFixed(VELOCITY_X) == false) vn1 = inner_prod(U_1,face_normal);
//                     if(mdistances[i_node2] <= 0 && face_geometry[2].IsFixed(VELOCITY_X) == false) vn2 = inner_prod(U_2,face_normal);
//
//                     double edge01 = 0.5*(vn0+vn1)*0.333333333333333333333333333333*0.5;
//                     double edge02 = 0.5*(vn0+vn2)*0.333333333333333333333333333333*0.5;
//                     double edge12 = 0.5*(vn2+vn2)*0.333333333333333333333333333333*0.5;
//
//                     vel_correction[i_node0] += edge01 + edge02;
//                     vel_correction[i_node1] += edge01  + edge12;
//                     vel_correction[i_node2] += edge02 + edge12;
//
// /*		     double tmp = 0.333333333333333333333333333333333*0.333333333333333333333333333333333*(vn0+vn1+vn2);
//                      vel_correction[i_node0] += tmp;
//                     vel_correction[i_node1] += tmp;
//                     vel_correction[i_node2] += tmp;       */
//                 }
//             }
//
//            //slip condition
//             int slip_size = mSlipBoundaryList.size();
//             #pragma omp parallel for firstprivate(slip_size)
//             for (int i_slip = 0; i_slip < slip_size; i_slip++)
//             {
//                 unsigned int i_node = mSlipBoundaryList[i_slip];
//                 double dist = mdistances[i_node];
//                 if (dist <= 0.0)
//                 {
//                     array_1d<double, TDim>& U_i = VelArray[i_node];
//                     array_1d<double, TDim>& an_i = mSlipNormal[i_node];
//                     double normalization = 0.0;
//                     for (unsigned int comp = 0; comp < TDim; comp++)
//                     {
//                         normalization += an_i[comp] * an_i[comp];
//                     }
//                      //tangential momentum as difference between original and normal momentum
//                     double coeff = vel_correction[i_node] / normalization;
//                     for (unsigned int comp = 0; comp < TDim; comp++)
//                         U_i[comp] +=  coeff * an_i[comp];
//                 }
//             }
//         }


        //fixed condition
        int fixed_size = mFixedVelocities.size();
        #pragma omp parallel for firstprivate(fixed_size)
        for (int i_velocity = 0; i_velocity < fixed_size; i_velocity++)
        {
            unsigned int i_node = mFixedVelocities[i_velocity];
            double dist = mdistances[i_node];
            if (dist <= 0.0)
            {
                const array_1d<double, TDim>& u_i_fix = mFixedVelocitiesValues[i_velocity];
                array_1d<double, TDim>& u_i = VelArray[i_node];
                for (unsigned int comp = 0; comp < TDim; comp++)
                    u_i[comp] = u_i_fix[comp];
            }
        }
        KRATOS_CATCH ("")
    }
    //********************************
    //function to compute coefficients
    void ExtrapolateValues (unsigned int extrapolation_layers)
    {
        KRATOS_TRY
        //ensure that corner nodes are wet if all of the nodes around them have a negative distance
//         typedef Node < 3 > PointType;
//         typedef PointerVector<PointType > PointVector;
//         typedef PointVector::iterator PointIterator;
        mr_matrix_container.FillScalarFromDatabase (DISTANCE, mdistances,mr_model_part.Nodes() );

        #pragma omp parallel for
        for ( int i_node = 0; i_node < static_cast<int>(mr_model_part.Nodes().size()); i_node++)
            mis_visited[i_node] = 0.0;

        boost::numeric::ublas::vector<int> layers(mr_model_part.Nodes().size(),-1);
// 	std::vector<int> layer_color(mr_model_part.Nodes().size(),-1000);
        boost::numeric::ublas::vector<int> layer_limits(extrapolation_layers+1);


        layer_limits[0] = 0;
        int layer_counter = -1;
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>( mr_model_part.Nodes().size()); i_node++)
        {
            if(mdistances[i_node] < 0.0)
        {
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
            {
                    //get global index of neighbouring node j
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    if(mdistances[j_neighbour] >= 0.0 && mis_visited[i_node] == 0)
                {
                        #pragma omp critical
                        layers[++layer_counter] = i_node;

                        mis_visited[i_node] = 1;
                        break;
                        }
                    }
                }
            else
            {
                mvel_n1[i_node] = ZeroVector (TDim);
                mvel_n[i_node] = ZeroVector (TDim);
                mPn[i_node] = 0.0;
                mPn1[i_node] = 0.0;
                mXi[i_node] = ZeroVector (TDim);
            }

        }
        layer_limits[1] = layer_counter;

        //fill the following layers by neighbour relationships
        //each layer fills the following
        for (unsigned int il = 0; il < extrapolation_layers - 1; il++)
        {
            //parallelization not trivial
            for(unsigned int iii = static_cast<unsigned int>(layer_limits[il]); iii<static_cast<unsigned int>(layer_limits[il+1]); iii++)
            {
                unsigned int i_node = layers[iii];
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    if(mdistances[j_neighbour] >= 0.0 && mis_visited[j_neighbour] == 0)
                    {
                        layers[layer_counter++] = j_neighbour;
                        mis_visited[j_neighbour] = il+2;
                    }
                }
            }
            layer_limits[il+2] = layer_counter;
        }

        array_1d<double, TDim > aux, aux_proj;

        //ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
        //double delta_t = CurrentProcessInfo[DELTA_TIME];

        //fill the pressure projection on the first layer inside the fluid
        //by extrapolating from the pressure projection on the layer -1 (the first layer completely inside the domain)
        #pragma omp parallel for
        for(int i=layer_limits[0]; i<layer_limits[1]; i++)
        {
            unsigned int i_node = layers[i];
            noalias (aux_proj) = ZeroVector (TDim);
            double avg_number = 0.0;

            for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
            {
                unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                if( mis_visited[j_neighbour] == 0)
                {
                    const array_1d<double, TDim > & inside_press_grad = mXi[j_neighbour];
                    noalias (aux_proj) += inside_press_grad;
                    avg_number += 1.0;
                }
            }
            if (avg_number != 0.0) //this case means that it has some neighbours that are completely internal
            {
                aux_proj /= avg_number;
                noalias (mXi[i_node] ) = aux_proj;
            }
            else //case in which there is not a layer of nodes completely internal
            {
                array_1d<double,TDim>& xi = mXi[i_node];
                noalias ( xi ) = mRho*mBodyForce;
                noalias ( xi ) -= mRho*macc[i_node];
            }
        }

        //perform extrapolation layer by layer by making an average
        //of the neighbours of lower order
        /*        KRATOS_WATCH(extrapolation_layers)
        		for (unsigned int il = 0; il < extrapolation_layers; il++)
        			std::cout << layer_limits[il] << " ";
        		std::cout << std::endl;
        	std::cout << std::endl;

        	for (unsigned int il = 0; il < extrapolation_layers; il++)
        	{
        	  std::cout << "level = " << il << " nneighb = " << layer_limits[il+1] - layer_limits[il] << " -- ";
        	    for(unsigned int iii = layer_limits[il]; iii<layer_limits[il+1]; iii++)
        	      std::cout << layers[iii] << " ";

        	    std::cout << std::endl;
        	}
        	std::cout << std::endl;

                std::cout << " printing is visited " << std::endl;
        	  for (unsigned int i_node = 0; i_node < mr_model_part.Nodes().size(); i_node++)
        	    std::cout << mis_visited[i_node] << std::endl;
        	 std::cout << std::endl;*/

        for (int il = 1; il < static_cast<int>(extrapolation_layers); il++)
        {
            //parallelization of this loop not trivial
            for(int iii = layer_limits[il]; iii<layer_limits[il+1]; iii++)
            {

                unsigned int i_node = layers[iii];
                noalias (aux) = ZeroVector (TDim);
                noalias (aux_proj) = ZeroVector (TDim);
                double avg_number = 0.0;
                double pavg = 0.0;

                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    if (mis_visited[j_neighbour] < (il + 1) && mis_visited[j_neighbour] != 0)
                    {
                        const array_1d<double, TDim >& direction_vec = mEdgeDimensions[csr_index];
// 			noalias (direction_vec) -= coords_bottom;
                        const array_1d<double, TDim >& press_grad = mXi[j_neighbour]; //i->FastGetSolutionStepValue (PRESS_PROJ);
                        double temp = inner_prod (direction_vec, press_grad);
                        double pestimate = mPn[j_neighbour] + temp;
                        pavg += pestimate;
                        noalias (aux_proj) += press_grad;
                        noalias (aux) += mvel_n1[j_neighbour]; //i->FastGetSolutionStepValue (VELOCITY);
                        avg_number += 1.0;
                    }
                }
                if (avg_number != 0.0)
                {
                    aux /= avg_number;
                    pavg /= avg_number;
                    aux_proj /= avg_number;

// 		    KRATOS_WATCH(avg_number);
// 		    KRATOS_WATCH(aux);
// 		    KRATOS_WATCH(pavg);
// 		    KRATOS_WATCH(aux_proj);
                }
                else
                {
                    KRATOS_THROW_ERROR (std::runtime_error, "error in extrapolation:: no neighbours find on a extrapolation layer -- impossible", "");
                    //                                                    KRATOS_THROW_ERROR(std:logic_error,"error in extrapolation:: no neighbours find on a extrapolation layer -- impossible","");
                }
                mvel_n1[i_node] = aux;
                mvel_n[i_node] = aux;
                mPn[i_node] = pavg;
//  		mPn1[i_node] = pavg;
                mXi[i_node] = aux_proj;
            }
        }

        //mark nodes on which we will have to solve for convection
        //mark all of internal nodes
        #pragma omp parallel for
        for ( int i_node = 0; i_node < static_cast<int>(mr_model_part.Nodes().size()); i_node++)
        {
            if (mdistances[i_node] <= 0.0)
                mis_visited[i_node] = 1.0;
            else
                mis_visited[i_node] = 0.0;
        }

        //now mark all of the nodes up to the extrapolation layers - 1
        for (unsigned int il = 0; il < extrapolation_layers-1; il++)
        {
            #pragma omp parallel for
            for( int iii = static_cast<int>(layer_limits[il]); iii<static_cast<int>(layer_limits[il+1]); iii++)
            {
                unsigned int i_node = layers[iii];
                mis_visited[i_node] = 1.0;
            }
        }
        ApplyVelocityBC (mvel_n1);
// 	mr_matrix_container.WriteVectorToDatabase (VELOCITY, mvel_n1,  mr_model_part.Nodes() );
// 		KRATOS_WATCH("end of Extrapolate Values ")
// 	KRATOS_WATCH(mvel_n1)
// 	KRATOS_WATCH(mPn)
// 	KRATOS_WATCH(mPn1)
// 	KRATOS_WATCH(mXi)
// 	KRATOS_WATCH(mdistances)
#ifdef DEBUG_OUTPUT

        KRATOS_WATCH("end of extrapolate values - new")
        double aux_v=0.0;
        for (int i_node = 0; i_node < mvel_n1.size(); i_node++)
            aux_v += inner_prod(mvel_n1[i_node],mvel_n1[i_node]);
        double aux_xi=0.0;
        for (int i_node = 0; i_node < mvel_n1.size(); i_node++)
            aux_xi += inner_prod(mXi[i_node],mXi[i_node]);

        KRATOS_WATCH(inner_prod(mPn1,mPn1));
        KRATOS_WATCH(aux_v);
        KRATOS_WATCH(aux_xi);
#endif
        KRATOS_CATCH ("")
    }
    void ChangeSignToDistance()
    {
        KRATOS_TRY
        for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode != mr_model_part.NodesEnd();
                inode++)
        {
            double dist = inode->FastGetSolutionStepValue (DISTANCE);
            inode->FastGetSolutionStepValue (DISTANCE) = -dist;
        }
        KRATOS_CATCH ("")
    }
    void MarkNodesByDistance (double min, double max)
    {
        KRATOS_TRY
        #pragma omp parallel for
        for ( int i_node = 0; i_node < static_cast<int>(mr_model_part.Nodes().size()); i_node++)
        {
            double& dist = mdistances[i_node];
            if ( dist > min && dist < max )
                mis_visited[i_node] = 1.0;
            else
                mis_visited[i_node] = 0.0;
        }

        KRATOS_CATCH ("")
    }
    void SaveScalarVariableToOldStep (Variable<double>& rVar)
    {
        KRATOS_TRY
        for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode != mr_model_part.NodesEnd();
                inode++)
        {
            inode->FastGetSolutionStepValue (rVar, 1) = inode->FastGetSolutionStepValue (rVar);
        }
        KRATOS_CATCH ("")
    }
    void MarkExternalAndMixedNodes()
    {
        KRATOS_TRY

        #pragma omp parallel for
        for ( int i_node = 0; i_node < static_cast<int>(mr_model_part.Nodes().size()); i_node++)
            mis_visited[i_node] = 0;

        for (unsigned int i_node = 0; i_node < mr_model_part.Nodes().size(); i_node++)
        {
            if(mdistances[i_node] > 0.0)
        {
                mis_visited[i_node] = 1;
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
            {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    mis_visited[j_neighbour] = 1;
                }
            }
        }

        KRATOS_CATCH ("")
    }
    void MarkInternalAndMixedNodes()
    {
        KRATOS_TRY
        #pragma omp parallel for
        for ( int i_node = 0; i_node < static_cast<int>(mr_model_part.Nodes().size()); i_node++)
            mis_visited[i_node] = 0;

        for (unsigned int i_node = 0; i_node < mr_model_part.Nodes().size(); i_node++)
        {
            if(mdistances[i_node] <= 0.0)
        {
                mis_visited[i_node] = 1;
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
            {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    mis_visited[j_neighbour] = 1;
                }
            }
        }

        KRATOS_CATCH ("")
    }
    void MarkInternalNodes()
    {
        KRATOS_TRY

        #pragma omp parallel for
        for ( int i_node = 0; i_node < static_cast<int>(mr_model_part.Nodes().size()); i_node++)
        {
            if(mdistances[i_node] <= 0.0)
                mis_visited[i_node] = 1;
            else
                mis_visited[i_node] = 0;
        }

        KRATOS_CATCH ("")
    }
    //**************************************
    //function to calculate the area normals
    void CalculateNormals (ModelPart::ConditionsContainerType& rConditions)
    {
        KRATOS_TRY
        //calculate area normals face-by-face
        array_1d<double, 3 > area_normal;
        //2D case
        if (TDim == 2)
        {
            for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
                CalculateNormal2D (cond_it, area_normal);
        }//3D case
        else if (TDim == 3)
        {
            //help vectors for cross product
            array_1d<double, 3 > v1;
            array_1d<double, 3 > v2;
            for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
                CalculateNormal3D (cond_it, area_normal, v1, v2);
        }
// area_normal *= -1; //CHAPUZA: REMOVE!!!s
        //(re)initialize normals
        unsigned int n_nodes = mNodalFlag.size();
        mInOutNormal.resize (n_nodes);
        mSlipNormal.resize (n_nodes);
        for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
        {
            noalias (mSlipNormal[i_node]) = ZeroVector (TDim);
            mis_slip[i_node] = false;
            noalias (mInOutNormal[i_node]) = ZeroVector (TDim);
        }
        //loop over all faces
        const double node_factor = 1.0 / TDim;
        for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
        {
            //get geometry data of the face
            Geometry<Node < 3 > >& face_geometry = cond_it->GetGeometry();
            //reference for area normal of the face
            array_1d<double, 3 > & face_normal = cond_it->GetValue (NORMAL);
            //slip condition
            if (static_cast<bool> (cond_it->GetValue (IS_STRUCTURE) ) == true)
                for (unsigned int if_node = 0; if_node < TDim; if_node++)
                {
                    unsigned int i_node = static_cast<unsigned int> (face_geometry[if_node].FastGetSolutionStepValue (AUX_INDEX) );
                    array_1d<double, TDim>& slip_normal = mSlipNormal[i_node];
                    mis_slip[i_node] = true;
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
            if (mis_slip[i_node] == true)
                tempmSlipBoundaryList.push_back (i_node);
            mis_slip[i_node] = false;
        }
        mSlipBoundaryList.resize (tempmSlipBoundaryList.size(),false);
        #pragma omp parallel for
        for (int i=0; i<static_cast<int> (tempmSlipBoundaryList.size() ); i++)
            mSlipBoundaryList[i] = tempmSlipBoundaryList[i];

        //check that all of the normals are not zero
        for (int i=0; i<static_cast<int> (mSlipBoundaryList.size() ); i++)
        {
            unsigned int i_node = mSlipBoundaryList[i];
            double tmp = norm_2(mSlipNormal[i_node]);
            if(tmp < 1e-20)
                KRATOS_THROW_ERROR(std::logic_error,"found a slip node with zero normal on node with id",i_node+1)
            }


        //loop over all faces to fill inlet outlet
        for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
        {
            //get geometry data of the face
            Geometry<Node < 3 > >& face_geometry = cond_it->GetGeometry();
            //reference for area normal of the face
            array_1d<double, 3 > & face_normal = cond_it->GetValue (NORMAL);
	    bool is_inlet_or_outlet = false;
            if (cond_it->GetValue (IS_STRUCTURE) == 0) is_inlet_or_outlet = true;
	    else
	    {
	      for (unsigned int if_node = 0; if_node < TDim; if_node++)
                    if (face_geometry[if_node].IsFixed (VELOCITY_X)  )
		    is_inlet_or_outlet = true;
	    }
            //slip condition
            if (is_inlet_or_outlet) //the opposite of the loop before
                for (unsigned int if_node = 0; if_node < TDim; if_node++)
                {
                    unsigned int i_node = static_cast<unsigned int> (face_geometry[if_node].FastGetSolutionStepValue (AUX_INDEX) );
                    array_1d<double, TDim>& inout_normal = mInOutNormal[i_node];
                    mis_slip[i_node] = true; //reutilize it!
                    for (unsigned int comp = 0; comp < TDim; comp++)
                    {
                        inout_normal[comp] += node_factor * face_normal[comp];
                    }
                }
        }
        
//        KRATOS_WATCH( mInOutNormal[7-1] );
//        KRATOS_THROW_ERROR(std::logic_error,"remove line 2318 " ,"");
       
        std::vector< unsigned int> tempmInOutBoundaryList;
        for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
        {
            if (mis_slip[i_node] == true)
                tempmInOutBoundaryList.push_back (i_node);
        }
        mInOutBoundaryList.resize (tempmInOutBoundaryList.size(),false);
        #pragma omp parallel for
        for (int i=0; i<static_cast<int> (tempmInOutBoundaryList.size() ); i++)
            mInOutBoundaryList[i] = tempmInOutBoundaryList[i];
        //store for future use the list of slip nodes
        #pragma omp parallel for
        for (int i=0; i<static_cast<int> (mis_slip.size() ); i++)
            mis_slip[ i ] = false;
        #pragma omp parallel for
        for (int i=0; i<static_cast<int> (mSlipBoundaryList.size() ); i++)
            mis_slip[ mSlipBoundaryList[i] ] = true;
        KRATOS_CATCH ("")
    }
    //*******************************
    //function to free dynamic memory
    void Clear()
    {
        KRATOS_TRY
        mViscosity.clear();
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
        mPiConvection.clear();
        mphi_n.clear();
        mphi_n1.clear();
        mEps.clear();
//         mD.clear();
        mA.clear();
        mB.clear();
        mdiv_error.clear();
	mWallReductionFactor.clear();
	mdiag_stiffness.clear();
        mis_slip.clear();
        mis_visited.clear();
        macc.clear();

        KRATOS_CATCH ("")
    }
    void ConvectDistance()
    {
        KRATOS_TRY
        //variables for node based data handling
        ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
        int n_nodes = rNodes.size();
        //storage of nodal values in local variables
        ValuesVectorType rhs, WorkConvection;
        rhs.resize (n_nodes);
        WorkConvection.resize (n_nodes);
        ValuesVectorType active_nodes;
        active_nodes.resize (n_nodes);
//         mr_matrix_container.FillScalarFromDatabase (POROSITY, mEps, mr_model_part.Nodes() );
        //read variables from Kratos
//         mr_matrix_container.FillVectorFromDatabase (VELOCITY, mvel_n1, mr_model_part.Nodes() );
//         mr_matrix_container.FillOldVectorFromDatabase (VELOCITY, mvel_n, mr_model_part.Nodes() );
        mr_matrix_container.FillScalarFromDatabase (DISTANCE, mphi_n1, mr_model_part.Nodes() );
        mr_matrix_container.FillOldScalarFromDatabase (DISTANCE, mphi_n, mr_model_part.Nodes() );
	//get the "fresh" values to be fixed_size
        for (unsigned int i=0; i< mDistanceValuesList.size(); i++)
	{
	  mDistanceValuesList[ i ] = mphi_n1[ mDistanceBoundaryList[i] ];
	}
        //mr_matrix_container.AssignVectorToVector(mphi_n1, mphi_n); //mWork = mphi_n
        //            //chapuza
        //            //set the distance to zero when it tries to go out of the pressure boundary
        //            int pressure_size = mPressureOutletList.size();
        //            #pragma omp parallel for firstprivate(pressure_size)
        //            for (int iii = 0; iii < pressure_size; iii++)
        //            {
        //                unsigned int i_node = mPressureOutletList[iii];
        //                mphi_n1[i_node] = fabs(mphi_n1[i_node]);
        //                mphi_n[i_node] = fabs(mphi_n[i_node]);
        //            }
        //create and fill a vector of nodes for which we want to convect the velocity
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            active_nodes[i_node] = mis_visited[i_node];
        }
//         ComputeConvectiveProjection(mPiConvection,mphi_n1,mEps,mvel_n1);
// 	ComputeLimitor(mPiConvection,mphi_n1,mBeta,mvel_n1,mEdgeDimensions);
        //            mr_matrix_container.WriteScalarToDatabase(TEMPERATURE, active_nodes, rNodes);
        //read time step size from Kratos
        ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
        double delta_t = CurrentProcessInfo[DELTA_TIME];
        double n_substeps = mnumsubsteps;
//            del
        double delta_t_substep = delta_t/n_substeps;
        for (unsigned int substep = 0; substep<n_substeps; substep++)
        {
            mr_matrix_container.AssignVectorToVector (mphi_n, WorkConvection); //mWork = mphi_n
            //first step of Runge Kutta
            // 				mr_matrix_container.AssignVectorToVector(mphi_n,mphi_n1); //mphi_n1 = mphi_n
            mr_matrix_container.SetToZero (rhs);
            ComputeConvectiveProjection (mPiConvection,mphi_n1,mEps,mvel_n1);
            ComputeLimitor (mPiConvection,mphi_n1,mBeta,mvel_n1,mEdgeDimensions);
            CalculateRHS_convection (mphi_n1, mvel_n1, rhs, active_nodes);
            mr_matrix_container.Add_Minv_value (WorkConvection, WorkConvection, delta_t_substep / 6.0, mr_matrix_container.GetInvertedMass(), rhs);
            mr_matrix_container.Add_Minv_value (mphi_n1, mphi_n, 0.5 * delta_t_substep, mr_matrix_container.GetInvertedMass(), rhs);
            ApplyDistanceBC();
            //second step
            mr_matrix_container.SetToZero (rhs);
            ComputeConvectiveProjection (mPiConvection,mphi_n1,mEps,mvel_n1);
            ComputeLimitor (mPiConvection,mphi_n1,mBeta,mvel_n1,mEdgeDimensions);
            CalculateRHS_convection (mphi_n1, mvel_n1, rhs, active_nodes);
            mr_matrix_container.Add_Minv_value (WorkConvection, WorkConvection, delta_t_substep / 3.0, mr_matrix_container.GetInvertedMass(), rhs);
            mr_matrix_container.Add_Minv_value (mphi_n1, mphi_n, 0.5 * delta_t_substep, mr_matrix_container.GetInvertedMass(), rhs);
	    ApplyDistanceBC();
            //third step
            mr_matrix_container.SetToZero (rhs);
            ComputeConvectiveProjection (mPiConvection,mphi_n1,mEps,mvel_n1);
            ComputeLimitor (mPiConvection,mphi_n1,mBeta,mvel_n1,mEdgeDimensions);
            CalculateRHS_convection (mphi_n1, mvel_n1, rhs, active_nodes);
            mr_matrix_container.Add_Minv_value (WorkConvection, WorkConvection, delta_t_substep / 3.0, mr_matrix_container.GetInvertedMass(), rhs);
            mr_matrix_container.Add_Minv_value (mphi_n1, mphi_n, delta_t_substep, mr_matrix_container.GetInvertedMass(), rhs);
	    ApplyDistanceBC();
            //fourth step
            mr_matrix_container.SetToZero (rhs);
            ComputeConvectiveProjection (mPiConvection,mphi_n1,mEps,mvel_n1);
            ComputeLimitor (mPiConvection,mphi_n1,mBeta,mvel_n1,mEdgeDimensions);
            CalculateRHS_convection (mphi_n1, mvel_n1, rhs, active_nodes);
            mr_matrix_container.Add_Minv_value (WorkConvection, WorkConvection, delta_t_substep / 6.0, mr_matrix_container.GetInvertedMass(), rhs);
	    ApplyDistanceBC();
            //compute right-hand side
            mr_matrix_container.AssignVectorToVector (WorkConvection, mphi_n1);
            mr_matrix_container.AssignVectorToVector (mphi_n1, mphi_n);
        }
        //            // make sure that boundary nodes that are very close to the free surface get wet
        //            int slip_size = mSlipBoundaryList.size();
        //            #pragma omp parallel for firstprivate(slip_size)
        //            for (int i_slip = 0; i_slip < slip_size; i_slip++) {
        //                unsigned int i_node = mSlipBoundaryList[i_slip];
        //                const double& h_i = mHmin[i_node];
        //                double& dist_i = mphi_n1[i_node];
        //
        //                if(dist_i > 0.0 && dist_i < 0.5*h_i)
        //                {
        //                    //loop to all the edges surrounding node I
        //                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
        //                    {
        //                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
        //                        if(mphi_n1[j_neighbour] <= 0.0)
        //                            dist_i = -0.01 * h_i;
        //                    }
        //                }
        //
        //            }
        //            int fixed_size = mFixedVelocities.size();
        //            #pragma omp parallel for firstprivate(fixed_size)
        //            for (int i_velocity = 0; i_velocity < fixed_size; i_velocity++) {
        //                unsigned int i_node = mFixedVelocities[i_velocity];
        //                const double& h_i = mHmin[i_node];
        //                double& dist_i = mphi_n1[i_node];
        //
        //                if(dist_i > 0.0 && dist_i < 0.5*h_i)
        //                {
        //                    //loop to all the edges surrounding node I
        //                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
        //                    {
        //                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
        //                        if(mphi_n1[j_neighbour] <= 0.0)
        //                            dist_i = -0.01 * h_i;
        //                    }
        //                }
        //            }
        //wetten corner nodes if needed
        int corner_size = mcorner_nodes.size();
        for (int i = 0; i < corner_size; i++)
        {
            int i_node = mcorner_nodes[i];
            bool to_be_wettened = true;
            double min_dist = 0.0;
            for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
            {
                unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                double neighb_dist = mphi_n1[j_neighbour];
                if (min_dist > neighb_dist)
                    min_dist = neighb_dist;
                if (neighb_dist >= 0.0)
                {
                    to_be_wettened=false;
                }
            }
            if (to_be_wettened==true)
                mphi_n1[i_node] = min_dist;
        }
        mr_matrix_container.WriteScalarToDatabase (DISTANCE, mphi_n1, mr_model_part.Nodes() );
        KRATOS_CATCH ("")
    }
    void ReduceTimeStep (ModelPart& rModelPart, double NewTime)
    {
        KRATOS_TRY
        /*
            double current_time = rModelPart.GetProcessInfo()[TIME];
            double current_delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
            double old_time = current_time - current_delta_time;
            double new_reduced_time = NewTtime;
            double new_delta_time = new_reduced_time - old_time;
            rModelPart.GetProcessInfo()[TIME] = new_reduced_time;
            rModelPart.GetProcessInfo()[DELTA_TIME] = new_delta_time;
            //now copy the database from the old step on the top of the current step
            int step_data_size = ThisModelPart.GetNodalSolutionStepDataSize();
            double* current_data = (pnode)->SolutionStepData().Data(0);
            double* old_data     = (pnode)->SolutionStepData().Data(1);
            for (int j = 0; j < step_data_size; j++)
                current_data[j] = old_data[j];
         */
        rModelPart.OverwriteSolutionStepData (1, 0);
        rModelPart.GetProcessInfo().SetCurrentTime (NewTime);
        KRATOS_CATCH ("error in reducing the time step")
    }
    bool CheckDistanceConvection()
    {
        int n_large_distance_gradient = 0;
        array_1d<double, TDim> grad_d;
        ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
        int n_nodes = rNodes.size();
        //calculate gradient of distance on the nodes and count occurrences of large gradients (that indicate a failure)
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            double dist = mdistances[i_node];
            if (dist <= 0.0)
            {
                for (unsigned int comp = 0; comp < TDim; comp++)
                    grad_d[comp] = 0.0;
                double dist_i = mdistances[i_node];
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
                {
                    //get global index of neighbouring node j
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                    const double& dist_j = mdistances[j_neighbour];
                    //projection of pressure gradients
                    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues() [csr_index];
                    edge_ij.Add_grad_p (grad_d, dist_i, dist_j);
                }
                const double& m_inv = mr_matrix_container.GetInvertedMass() [i_node];
                for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                    grad_d[l_comp] *= m_inv;
                double norm_grad = norm_2 (grad_d);
                if (norm_grad > 1.5) //large gradient found
                    n_large_distance_gradient += 1;
            }
        }
        if (n_large_distance_gradient != 0)
        {
            bool success = false;
            return success;
        }
        else
        {
            bool success = true;
            return success;
        }
    }
    void ActivateWallResistance (double Ywall)
    {
        mWallLawIsActive = true;
        mY_wall = Ywall;
	double max_angle_overall = 0.0;
	//compute wall reduction factor
	         //slip condition
        int slip_size = mSlipBoundaryList.size();
        #pragma omp parallel for firstprivate(slip_size)
        for (int i_slip = 0; i_slip < slip_size; i_slip++)
        {
	    unsigned int i_node = mSlipBoundaryList[i_slip];
            /*            const array_1d<double, TDim>& an_i = mSlipNormal[i_node];
            double AI = norm_2(an_i);
                        array_1d<double,TDim> nI = an_i/AI;
                        double min_dot_prod = 1.0;
                        for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
            {
                            unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                            const array_1d<double, TDim>& an_j = mSlipNormal[j_neighbour];
                double AJ = norm_2(an_j);
                            if(AJ > 1e-20) //...a slip node!
                            {
                                double tmp = 0.0;
                                for (unsigned int comp = 0; comp < TDim; comp++)
                                    tmp += nI[comp] * an_j[comp];
                                tmp /= AJ;
                    tmp = fabs(tmp);
                                if(tmp < min_dot_prod) min_dot_prod = tmp;
                            }
                        }
                        double max_angle = acos(min_dot_prod);
            // 	    max_angle *= 2.0;
            // 	    if(max_angle > 3.1415926*0.5) max_angle = 3.1415926*0.5;
                        if(max_angle > max_angle_overall) max_angle_overall = max_angle;*/
            mWallReductionFactor[i_node] = 1.0; //sin(max_angle) + 0.1; // pow(sin(max_angle),6) * 10.0 /** 100.0*/ ;
	}
	std::cout << "max angle between normals found in the model = " << max_angle_overall << std::endl;
// 	mr_matrix_container.WriteScalarToDatabase(YOUNG_MODULUS, mWallReductionFactor, mr_model_part.Nodes());
         //slip condition
//         int slip_size = mSlipBoundaryList.size();
//         #pragma omp parallel for firstprivate(slip_size)
//         for (int i_slip = 0; i_slip < slip_size; i_slip++)
//         {
// 	    unsigned int i_node = mSlipBoundaryList[i_slip];
// 	    double h = mHavg[i_node];
// 	    if(mY_wall < h)
// 	      mWallReductionFactor[i_node] = mY_wall/h;
// 	}
//
        int edge_size = medge_nodes.size();
        #pragma omp parallel for firstprivate(edge_size)
        for (int i = 0; i < edge_size; i++)
        {
            int i_node = medge_nodes[i];
            mWallReductionFactor[i_node] = medge_coefficient; //10.0;
    }
//
// 	//apply conditions on corners
        int corner_size = mcorner_nodes.size();
        for (int i = 0; i < corner_size; i++)
    {
            int i_node = mcorner_nodes[i];
            mWallReductionFactor[i_node] = mcorner_coefficient; //50.0;
        }
    }
    void ActivateClassicalWallResistance (double Ywall)
    {
        mWallLawIsActive = true;
        mY_wall = Ywall;
        for (unsigned int i = 0; i < mWallReductionFactor.size(); i++)
	    mWallReductionFactor[i] = 1.0 ;
    }
    double ComputeVolumeVariation()
    {
        ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
        double dt = CurrentProcessInfo[DELTA_TIME];
        //slip condition
        int inout_size = mInOutBoundaryList.size();
        double vol_var = 0.0;
        //#pragma omp parallel for firstprivate(slip_size)
        for (int i = 0; i < inout_size; i++)
        {
            unsigned int i_node = mInOutBoundaryList[i];
            double dist = mdistances[i_node];
            if (dist <= 0.0)
            {
                const array_1d<double, TDim>& U_i = mvel_n1[i_node];
                const array_1d<double, TDim>& an_i = mInOutNormal[i_node];
                double projection_length = 0.0;
                for (unsigned int comp = 0; comp < TDim; comp++)
                {
                    projection_length += U_i[comp] * an_i[comp];
                }
                vol_var += projection_length;
            }
        }
        return -vol_var * dt;
    }
    double ComputeWetVolume()
    {
        KRATOS_TRY
        mr_matrix_container.FillScalarFromDatabase (DISTANCE, mdistances, mr_model_part.Nodes() );
        //slip condition
        double wet_volume = 0.0;
        //#pragma omp parallel for firstprivate(slip_size)
        for (int i = 0; i < static_cast<int> (mdistances.size() ); i++)
        {
            double dist = mdistances[i];
            const double m = mr_matrix_container.GetLumpedMass() [i];
	    double porosity = mEps[i];
            if (dist <= 0.0)
            {
                wet_volume += m/porosity;
            }
        }
        return wet_volume;
        KRATOS_CATCH ("");
    }
    double ComputeTotalVolume()
    {
        KRATOS_TRY
        mr_matrix_container.FillScalarFromDatabase (DISTANCE, mdistances, mr_model_part.Nodes() );
        //slip condition
        double volume = 0.0;
        //#pragma omp parallel for firstprivate(slip_size)
        for (int i = 0; i < static_cast<int> (mdistances.size() ); i++)
        {
            const double m = mr_matrix_container.GetLumpedMass() [i];
            double porosity = mEps[i];
            volume += m/porosity;
        }
        return volume;
        KRATOS_CATCH ("");
    }
    void DiscreteVolumeCorrection (double expected_volume, double measured_volume)
    {
        double volume_error = expected_volume - measured_volume;
        if (measured_volume < expected_volume)
        {
            double layer_volume = 0.0;
            std::vector<unsigned int> first_outside;
            int n_nodes = mdistances.size();
            // find list of the first nodes outside of the fluid and compute their volume
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                double dist = mdistances[i_node];
                if (dist > 0.0) //node is outside domain
                {
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                        if (mdistances[j_neighbour] <= 0.0)
                        {
                            const double nodal_mass = 1.0 / mr_matrix_container.GetInvertedMass() [i_node];
                            if (nodal_mass < volume_error - layer_volume)
                            {
                                first_outside.push_back (i_node);
                                layer_volume += nodal_mass;
				                break;
                            }
                            //const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];
                            //layer_volume += 1.0/m_inv;
                        }
                    }
                }
            }
//				std::cout << ", layer_volume: " << layer_volume  << std::endl;
//              if (measured_volume + layer_volume <= expected_volume)
            {
                // mark the nodes in the outside layer with a small negative distance
                for (unsigned int i=0; i<first_outside.size(); i++)
                {
                    unsigned int i_node = first_outside[i];
                    mdistances[i_node] = -mHavg[i_node];
                }
            }
        }
        mr_matrix_container.WriteScalarToDatabase (DISTANCE, mdistances, mr_model_part.Nodes() );
        //if (measured_volume < expected_volume)
        //         {
        //             double layer_volume = 0.0;
        //             std::vector<unsigned int> first_outside;
        //             int n_nodes = mdistances.size();
        //             //find list of the first nodes outside of the fluid and compute their volume
        //             for (int i_node = 0; i_node < n_nodes; i_node++)
        //             {
        //                 double dist = mdistances[i_node];
        //                 if (dist > 0.0) //node is outside domain
        //                 {
        //                     for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
        //                     {
        //                         unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
        //                         if(mdistances[j_neighbour] <= 0.0)
        //                         {
        //                             first_outside.push_back(i_node);
        //                             const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];
        //                             layer_volume += 1.0/m_inv;
        //                         }
        //                     }
        //                 }
        //             }
        //             if (measured_volume + layer_volume <= expected_volume)
        //             {
        //                 //mark the nodes in the outside layer with a small negative distance
        //                 for(unsigned int i=0; i<first_outside.size(); i++)
        //                 {
        //                     unsigned int i_node = first_outside[i];
        //                     mdistances[i_node] = -mHavg[i_node];
        //                 }
        //             }
        //         }
        //         mr_matrix_container.WriteScalarToDatabase(DISTANCE, mdistances, mr_model_part.Nodes());
    }
    void SetWallReductionCoefficients (double corner_coefficient, double edge_coefficient)
        {
        mcorner_coefficient = corner_coefficient;
        medge_coefficient = edge_coefficient;
    }
    void ContinuousVolumeCorrection (double expected_volume, double measured_volume)
    {
			double volume_error = expected_volume - measured_volume;
			if (volume_error == 0.0)
				return ;
           if (measured_volume < expected_volume)
            {
             double layer_volume = 0.0;
                std::vector<unsigned int> first_outside;
                int n_nodes = mdistances.size();
                // find list of the first nodes outside of the fluid and compute their volume
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{
					double dist = mdistances[i_node];
					bool is_bubble = true;
					bool is_first_outside = false;
					if (dist > 0.0) //node is outside domain
					{
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
						{
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                        if (mdistances[j_neighbour] <= 0.0)
							{
								is_first_outside = true;
							}
							else
								is_bubble = false;
						}
					}
                if (is_first_outside && !is_bubble)
					{
                    const double nodal_mass = 1.0 / mr_matrix_container.GetInvertedMass() [i_node];
                    first_outside.push_back (i_node);
						layer_volume += nodal_mass;
//                     if(nodal_mass > volume_error - layer_volume)
//                     {
//                         extra_volume += nodal_mass;
//                     }
						}
					}
//				std::cout << ", layer_volume: " << layer_volume  << std::endl;
				if (layer_volume == 0.00)
					return;
				double ratio = volume_error / layer_volume;
            if (ratio > 1.0) ratio = 1.0;
//             KRATOS_WATCH (ratio);
            if (ratio < 0.1)  // NO correction for less than 10% error
                return;
	    double average_layer_h = 0.0;
            for (unsigned int i=0; i<first_outside.size(); i++)
            {
                unsigned int i_node = first_outside[i];
                    average_layer_h += mHavg[i_node];
            }
            average_layer_h /= static_cast<double> (first_outside.size() );
	    for (int i_node = 0; i_node < n_nodes; i_node++)
	        mdistances[i_node] -=  average_layer_h* ratio;
//             if((ratio < 1.00))
//             {
//                 // mark the nodes in the outside layer with a small negative distance
//                 for(unsigned int i=0; i<first_outside.size(); i++)
//                 {
//                     unsigned int i_node = first_outside[i];
//                     mdistances[i_node] -= mHavg[i_node] * ratio;
//                 }
//             }
//             else
//             {
//                 // mark the nodes in the outside layer with a small negative distance
//                 for(unsigned int i=0; i<first_outside.size(); i++)
//                 {
//                     unsigned int i_node = first_outside[i];
//                     mdistances[i_node] = -mHavg[i_node];
//                 }
//             }
                }
        mr_matrix_container.WriteScalarToDatabase (DISTANCE, mdistances, mr_model_part.Nodes() );

			return;
        }
//        void FindBubbles()
//         {
// 			int n_nodes = mdistances.size();
//         ValuesVectorType last_air (n_nodes);
//         mr_matrix_container.SetToZero (last_air);
//         mr_matrix_container.FillScalarFromDatabase (LAST_AIR, last_air, mr_model_part.Nodes() );
// 			const int max_bubble_nodes = 12;
// 			const int min_bubble_nodes = 2;
//         #pragma omp parallel for
//         for ( int i_node = 0; i_node < static_cast<int>(mr_model_part.Nodes().size()); i_node++)
//             mis_visited[i_node] = 0;
// 
// 			// loop over the nodes to find a outside node.
// 			for (int i_node = 0; i_node < n_nodes; i_node++)
// 			{
// 				double dist = mdistances[i_node];
//             if ( (mis_visited[i_node] == 0) && (dist > 0.0) ) // node is outside the domain and has not visited yet
// 				{
//                 std::vector<int> outside_nodes (n_nodes,0);
// 					outside_nodes[0] = i_node;
//                 mis_visited[i_node] = 1;
// 					int n_outside = 1;
//                 for (int i = 0 ; i < n_outside ; i++) // loop over founded outside nodes. NOTE: n_outside is increasing inside the loop
// 					{
// 						int this_node = outside_nodes[i];
// 						// loop over neighbours of this node
//                     for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [this_node]; csr_index != mr_matrix_container.GetRowStartIndex() [this_node + 1]; csr_index++)
// 						{
//                         unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
//                         if ( (mis_visited[j_neighbour] == 0) && (mdistances[j_neighbour] >= 0.0) ) // the neighbour node is outside the fluid and not visited yet
// 							{
// 									outside_nodes[n_outside] = j_neighbour;
// 									n_outside++;
// 							}
//                         mis_visited[j_neighbour] = 1;
// 						}
// 					}
// 					//KRATOS_WATCH(i_node);
// 					//KRATOS_WATCH(n_outside);
// 					//KRATOS_WATCH(is_first_outside);
//                 if ( (n_outside <= max_bubble_nodes) && (n_outside >= min_bubble_nodes) )
// 					{
// 						//KRATOS_WATCH(i_node);
// 						//KRATOS_WATCH(n_outside);
//                     for (int i = 0 ; i < n_outside ; i++)
// 							last_air[outside_nodes[i]] = 1.00;
// 					}
// 				}
// 			}
//         mr_matrix_container.WriteScalarToDatabase (LAST_AIR, last_air, mr_model_part.Nodes() );
//         }
// 
//        void FindColdShots()
//         {
// 			int n_nodes = mdistances.size();
// 			ValuesVectorType cold_shots(n_nodes);
// 
// 			mr_matrix_container.SetToZero(cold_shots);
// 
// 			mr_matrix_container.FillScalarFromDatabase(LAST_AIR, cold_shots, mr_model_part.Nodes());
// 
// 			std::vector<bool> is_first_outside(n_nodes, 0);
// 
//             std::vector<unsigned int> first_outside;
// 
// 			// find list of the first nodes outside of the fluid
// 			for (int i_node = 0; i_node < n_nodes; i_node++)
// 			{
// 				double dist = mdistances[i_node];
// 				if (dist > 0.0) //node is outside domain
// 				{
// 					for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
// 					{
// 						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
// 						if(mdistances[j_neighbour] <= 0.0)
// 						{
// 							is_first_outside[i_node] = true;
// 							first_outside.push_back(i_node);
// 							break;
// 						}
// 					}
// 				}
// 			}
// 
// 
// 			std::vector<bool> is_cold_shot(is_first_outside);
// 
// 			// Now we check if all the neighbours of the first_outside nodes are first outside or inside and mark it as a possible cold shot
// 			for(unsigned int i=0; i<first_outside.size(); i++)
// 			{
// 				unsigned int i_node = first_outside[i];
// 				for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
// 				{
// 					unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
// 					if(!is_first_outside[j_neighbour])
// 					{
// 						is_cold_shot[i_node] = false;
// 						break;
// 					}
// 				}
// 			}
// 
// 
// 			//Now we have the possible cold shots and is time to check the gradient of convection
// 			for(unsigned int i=0; i<first_outside.size(); i++)
// 			{
// 				unsigned int i_node = first_outside[i];
// 				if(is_cold_shot[i_node])
// 				{
// 					for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
// 					{
// 						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
// 						if(mdistances[j_neighbour] <= 0.0)
// 						{
// 							
// 						}
// 					}
// 				}
// 			}
// 
// 
// 				
// 			// Adding the founded cold shots to the previous ones.
// 			for(int i_node = 0; i_node < n_nodes; i_node++)
// 				if(is_cold_shot[i_node])
// 					cold_shots[i_node]=1.00;
// 				
//             mr_matrix_container.WriteScalarToDatabase(LAST_AIR, cold_shots, mr_model_part.Nodes());
//         }

    void CalculatePorousResistanceLaw(unsigned int res_law)
    {
        //variables for node based data handling
// 	    ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
// 	    mr_matrix_container.FillScalarFromDatabase(DIAMETER, mD, mr_model_part.Nodes());
// 	    mr_matrix_container.FillScalarFromDatabase(POROSITY, mEps, mr_model_part.Nodes());
//  	    mr_matrix_container.FillScalarFromDatabase(LIN_DARCY_COEF, mA, mr_model_part.Nodes());
//  	    mr_matrix_container.FillScalarFromDatabase(NONLIN_DARCY_COEF, mB, mr_model_part.Nodes());
// 	  const double nu_i = mViscosity;
        if (res_law == 1)
        {
// 	  KRATOS_WATCH("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Calculating Ergun Darcy coefficients ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            /* if the chosen resistance law is ERGUN calculate Ergun A and B*/
            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                    inode != mr_model_part.NodesEnd();
                    inode++)
            {
                const double eps = inode->FastGetSolutionStepValue (POROSITY);
// 		KRATOS_WATCH("POROSITY	")
// 		KRATOS_WATCH(eps)
                const double d = inode->FastGetSolutionStepValue (DIAMETER);
// 		KRATOS_WATCH("DIAMETER	")
// 		KRATOS_WATCH(d)
//
// 		KRATOS_WATCH("VISCOSITY	")
// 		KRATOS_WATCH(mViscosity)
                double& a = inode-> FastGetSolutionStepValue (LIN_DARCY_COEF);
                double& b = inode-> FastGetSolutionStepValue (NONLIN_DARCY_COEF);
                if (eps < 1.0)
                {
                    double k_inv = 150.0 * (1.0 - eps) * (1.0 - eps) / (eps * eps * eps * d * d);
                    a = mViscosity * k_inv;
                    b = (1.75 / eps) * sqrt (k_inv / (150.0 * eps) );
// 		KRATOS_WATCH("PERMEABILITY	")
// 		KRATOS_WATCH(k_inv)
// 		KRATOS_WATCH("LIN DARCY COEFFICIENT	")
// 		KRATOS_WATCH(a)
// 		KRATOS_WATCH("NONLIN DARCY COEFFICIENT	")
// 		KRATOS_WATCH(b)
                }
                else
                {
                    a = 0;
                    b = 0;
                }
            }
        }
        else
        {
            /* whether it is a Custom Resistance law or NO resistance law is present ---> set to zero A and B for non porous nodes*/
            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                    inode != mr_model_part.NodesEnd();
                    inode++)
            {
                const double eps = inode->FastGetSolutionStepValue (POROSITY); /*reading from kratos database*/
                double& a = inode-> FastGetSolutionStepValue (LIN_DARCY_COEF); /*changing kratos database*/
                double& b = inode-> FastGetSolutionStepValue (NONLIN_DARCY_COEF);	/*changing kratos database*/
                if (eps == 1.0)
                {
                    a = 0;
                    b = 0;
                }
            }
        }
        mr_matrix_container.FillScalarFromDatabase (LIN_DARCY_COEF, mA, mr_model_part.Nodes() ); /*filling edgebased database reading from kratos database*/
        mr_matrix_container.FillScalarFromDatabase (NONLIN_DARCY_COEF, mB, mr_model_part.Nodes() ); /*filling edgebased database reading from kratos database*/
    }
private:
    double mMolecularViscosity;
    double mcorner_coefficient;
    double medge_coefficient;
    double mmax_dt;
    MatrixContainer& mr_matrix_container;
    ModelPart& mr_model_part;
    int mnumsubsteps;
    bool muse_mass_correction;
    //parameters controlling the wall law
    bool mWallLawIsActive;
    double mY_wall;
    //parameters for controlling the usage of the delta time in the stabilization
    double mstabdt_pressure_factor;
    double mstabdt_convection_factor;
    double medge_detection_angle;
    double mtau2_factor;
    bool massume_constant_dp;
    //nodal values
    ValuesVectorType mViscosity;
    //velocity vector U at time steps n and n+1
    CalcVectorType mWork, mvel_n, mvel_n1, mx, macc;
    //pressure vector p at time steps n and n+1
    ValuesVectorType mPn, mPn1;
    //coefficients
    ValuesVectorType mdistances;
    //minimum length of the edges surrounding edges surrounding each nodal point
    ValuesVectorType mHmin;
    ValuesVectorType mHavg;
    CalcVectorType mEdgeDimensions;
    //area normal
    CalcVectorType mSlipNormal;
    CalcVectorType mInOutNormal;
    //projection terms
    CalcVectorType mPi, mXi;
    //flag for first time step
    bool mFirstStep;
    //flag to differentiate interior and boundary nodes
    ValuesVectorType mNodalFlag;
    ValuesVectorType mWallReductionFactor;
    //lists of nodes with different types of boundary conditions
    IndicesVectorType mSlipBoundaryList, mPressureOutletList, mFixedVelocities, mInOutBoundaryList,mDistanceBoundaryList;
    ValuesVectorType mDistanceValuesList;
    CalcVectorType mFixedVelocitiesValues;
    //	ValuesVectorType mPressureOutlet;
    //intrinsic time step size
    ValuesVectorType mTauPressure;
    ValuesVectorType mTauConvection;
    ValuesVectorType mTau2;
    ValuesVectorType mdiv_error;
    boost::numeric::ublas::vector<bool> mis_slip;
    boost::numeric::ublas::vector<int> mis_visited;
    //variables for resolving pressure equation
    //laplacian matrix
    TSystemMatrixType mL;
    //constant variables
    double mRho;
    array_1d<double, TDim> mBodyForce;
    //variables for convection
    ValuesVectorType mphi_n;
    ValuesVectorType mphi_n1;
    CalcVectorType mPiConvection;
    ValuesVectorType mBeta;
    //variables for edge BCs
    IndicesVectorType medge_nodes;
    CalcVectorType medge_nodes_direction;
    IndicesVectorType mcorner_nodes;
    ValuesVectorType mEps;
    ValuesVectorType mdiag_stiffness;
//     ValuesVectorType mD;
    ValuesVectorType mA;
    ValuesVectorType mB;
    double mdelta_t_avg;
    double max_dt;
    double mshock_coeff;
    //***********************************************************
    //functions to calculate area normals for boundary conditions
    void CalculateNormal2D (ModelPart::ConditionsContainerType::iterator cond_it, array_1d<double, 3 > & area_normal)
    {
        Geometry<Node < 3 > >& face_geometry = (cond_it)->GetGeometry();
        area_normal[0] = face_geometry[1].Y() - face_geometry[0].Y();
        area_normal[1] = - (face_geometry[1].X() - face_geometry[0].X() );
        area_normal[2] = 0.00;
        noalias ( (cond_it)->GetValue (NORMAL) ) = area_normal;
    }
    void CalculateNormal3D (ModelPart::ConditionsContainerType::iterator cond_it, array_1d<double, 3 > & area_normal, array_1d<double, 3 > & v1, array_1d<double, 3 > & v2)
    {
        Geometry<Node < 3 > >& face_geometry = (cond_it)->GetGeometry();
        v1[0] = face_geometry[1].X() - face_geometry[0].X();
        v1[1] = face_geometry[1].Y() - face_geometry[0].Y();
        v1[2] = face_geometry[1].Z() - face_geometry[0].Z();
        v2[0] = face_geometry[2].X() - face_geometry[0].X();
        v2[1] = face_geometry[2].Y() - face_geometry[0].Y();
        v2[2] = face_geometry[2].Z() - face_geometry[0].Z();
        MathUtils<double>::CrossProduct (area_normal, v1, v2);
        area_normal *= -0.5;
        noalias ( (cond_it)->GetValue (NORMAL) ) = area_normal;
    }
    //*********************************************************
    //function to calculate minimum length of surrounding edges
    void CalculateEdgeLengths (ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY
        //get number of nodes
        unsigned int n_nodes = rNodes.size();
        //reserve memory for storage of nodal coordinates
        std::vector< array_1d<double, TDim > > position;
        position.resize (n_nodes);
        //get position of all nodes
        for (typename ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
        {
            //get the global index of the node
            unsigned int i_node = static_cast<unsigned int> (node_it->FastGetSolutionStepValue (AUX_INDEX) );
            //save its coordinates locally
            noalias (position[i_node]) = node_it->Coordinates();
            //initialize minimum edge length with relatively big values
//             mHmin[i_node] = 1e10;
        }
        ValuesVectorType& aaa = mr_matrix_container.GetHmin();
        for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
        {
            mHmin[i_node] = aaa[i_node];
            if (aaa[i_node] == 0.0)
                KRATOS_THROW_ERROR (std::logic_error,"found a 0 hmin on node",i_node);
        }
        //take unstructured meshes into account
        if (TDim == 2)
        {
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                double& h_i = mHavg[i_node];
                double& m_i = mr_matrix_container.GetLumpedMass() [i_node];
                // 						double& rho_i = mRho[i_node];
                h_i = sqrt (2.0 * m_i);
            }
        }
        else if (TDim == 3)
        {
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                double& h_i = mHavg[i_node];
                double& m_i = mr_matrix_container.GetLumpedMass() [i_node];
                // 						double& rho_i = mRho[i_node];
                h_i = pow (6.0 * m_i, 1.0 / 3.0);
            }
        }
        //compute edge coordinates
        for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
        {
            array_1d<double, TDim > & pos_i = position[i_node];
            for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
            {
                unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                array_1d<double, TDim > & pos_j = position[j_neighbour];
                array_1d<double, TDim>& l_k = mEdgeDimensions[csr_index];
                for (unsigned int comp = 0; comp < TDim; comp++)
                    l_k[comp] = pos_i[comp] - pos_j[comp];
            }
        }
        KRATOS_CATCH ("")
    }
    //*********************************************************************
    //function to calculate right-hand side of fractional momentum equation
    void CalculateRHS_convection (
        const ValuesVectorType& mphi,
        const CalcVectorType& convective_velocity,
        ValuesVectorType& rhs,
        ValuesVectorType& active_nodes
    )
    {
        KRATOS_TRY
        int n_nodes = mphi.size();
        //            //calculating the convective projection
        //#pragma omp parallel for
        //            for (int i_node = 0; i_node < n_nodes; i_node++)
        //            {
        //
        //                double& pi_i = mPiConvection[i_node];
        //                const double& phi_i = mphi[i_node];
        //
        //                //set to zero the projection
        //                pi_i = 0;
        //                if (active_nodes[i_node] != 0.0)
        //                {
        //
        //                    const array_1d<double, TDim>& a_i = convective_velocity[i_node];
        //
        //                    //loop to all the edges surrounding node I
        //                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
        //                    {
        //                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
        //
        //                        if (active_nodes[j_neighbour] != 0.0)
        //                        {
        //                            const array_1d<double, TDim>& a_j = convective_velocity[j_neighbour];
        //                            const double& phi_j = mphi[j_neighbour];
        //
        //                            CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
        //
        //                            edge_ij.Add_ConvectiveContribution(pi_i, a_i, phi_i, a_j, phi_j);
        //                        }
        //                    }
        //
        //                    //apply inverted mass matrix
        //                    const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];
        //                    pi_i *= m_inv;
        //                }
        //                // KRATOS_WATCH(pi_i);
        //                //                                num = fabs(num);
        //                //                                if(num > norm_vI*0.0001)
        //                //                                    mBeta[i_node] = 1.0 - num/denom;
        //                //                                else
        //                //                                    mBeta[i_node] = 1.0;
        //
        //            }
        //perform MPI syncronization
        //calculating the RHS
        double stab_low;
        double stab_high;
        array_1d<double, TDim> a_i;
        array_1d<double, TDim> a_j;
        #pragma omp parallel for private(stab_low,stab_high,a_i,a_j)
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            double& rhs_i = rhs[i_node];
            const double& h_i = mHavg[i_node];
            const double& phi_i = mphi[i_node];
            noalias (a_i) = convective_velocity[i_node];
            a_i /= mEps[i_node];
            const array_1d<double, TDim>& proj_i = mPiConvection[i_node];
            //                 const double& pi_i = mPiConvection[i_node];
            double pi_i = proj_i[0] * a_i[0];
            for (unsigned int l_comp = 1; l_comp < TDim; l_comp++)
                pi_i += proj_i[l_comp] * a_i[l_comp];
            //                double beta = mBeta[i_node];
            rhs_i = 0.0;
            if (active_nodes[i_node] != 0.0)
            {
                const double& beta = mBeta[i_node];
                double norm_a = a_i[0] * a_i[0];
                for (unsigned int l_comp = 1; l_comp < TDim; l_comp++)
                    norm_a += a_i[l_comp] * a_i[l_comp];
                norm_a = sqrt (norm_a);
                //loop to all the edges surrounding node I
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                    if (active_nodes[j_neighbour] != 0.0)
                    {
                        //double& rhs_j = rhs[j_neighbour];
                        const double& phi_j = mphi[j_neighbour];
                        noalias (a_j) = convective_velocity[j_neighbour];
                        a_j /= mEps[j_neighbour];
                        //                             const double& pi_j = mPiConvection[j_neighbour];
                        const array_1d<double, TDim>& proj_j = mPiConvection[j_neighbour];
                        double pi_j = proj_j[0] * a_i[0];
                        for (unsigned int l_comp = 1; l_comp < TDim; l_comp++)
                            pi_j += proj_j[l_comp] * a_i[l_comp];
                        CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues() [csr_index];
                        //convection operator
                        edge_ij.Sub_ConvectiveContribution (rhs_i, a_i, phi_i, a_j, phi_j); //esto funciona
                        //			    edge_ij.Sub_D_v(rhs_i, a_i*phi_i, a_i*phi_j);
                        //calculate stabilization part
                        edge_ij.CalculateConvectionStabilization_LOW (stab_low, a_i, phi_i, a_j, phi_j);
                        double edge_tau = mTauConvection[i_node];
                        edge_ij.CalculateConvectionStabilization_HIGH (stab_high, a_i, pi_i, a_j, pi_j);
                        edge_ij.Sub_StabContribution (rhs_i, edge_tau, 1.0, stab_low, stab_high);
                        double coeff = 0.5 * mshock_coeff; //=0.7*0.5;
                        double laplacian_ij = 0.0;
                        edge_ij.CalculateScalarLaplacian (laplacian_ij);
                        double capturing = laplacian_ij * (phi_j - phi_i);
                        //    			    rhs_i-= coeff*capturing*beta*norm_a*h_i;
                        double aaa = 0.0;
                        for (unsigned int k_comp = 0; k_comp < TDim; k_comp++)
                            for (unsigned int m_comp = 0; m_comp < TDim; m_comp++)
                                aaa += a_i[k_comp] * a_i[m_comp] * edge_ij.LaplacianIJ (k_comp, m_comp);
                        if (norm_a > 1e-10)
                        {
                            aaa /= (norm_a * norm_a);
                            double capturing2 = aaa * (phi_j - phi_i);
                            if (fabs (capturing) > fabs (capturing2) )
                                rhs_i -= coeff * (capturing - capturing2) * beta * norm_a * h_i;
                        }
                    }
                }
            }
            // KRATOS_WATCH(rhs_i);
        }

//                 int inout_size = mInOutBoundaryList.size();
//         //#pragma omp parallel for firstprivate(slip_size)
//         for (int i = 0; i < inout_size; i++)
//         {
//             unsigned int i_node = mInOutBoundaryList[i];
//             double dist = mdistances[i_node];
//             if (dist <= 0.0)
//             {
//                 const array_1d<double, TDim>& U_i = mvel_n1[i_node];
//                 const array_1d<double, TDim>& an_i = mInOutNormal[i_node];
//                 double projection_length = 0.0;
// 		double Ain = 0.0;
//                 for (unsigned int comp = 0; comp < TDim; comp++)
//                 {
//                     projection_length += U_i[comp] * an_i[comp];
// 		    Ain += an_i[comp]*an_i[comp];
//                 }
//
// 		double& rhs_i = rhs[i_node];
//
//                rhs_i += projection_length * mphi[i_node];
//             }
//         }
//         int inout_size = mInOutBoundaryList.size();
//         double vol_var = 0.0;
//         //#pragma omp parallel for firstprivate(slip_size)
//         for (int i = 0; i < inout_size; i++)
//         {
//             unsigned int i_node = mInOutBoundaryList[i];
//             double dist = mdistances[i_node];
// //             if (dist <= 0.0)
// //             {
//                 const array_1d<double, TDim>& U_i = mvel_n1[i_node];
//                 const array_1d<double, TDim>& an_i = mInOutNormal[i_node];
// 		 double A = norm_2(an_i);
//
// 		double projection_length = 0.0;
//                 for (unsigned int comp = 0; comp < TDim; comp++)
//                 {
//                     projection_length += U_i[comp] * an_i[comp];
//                 }
//
//                 double& rhs_i = rhs[i_node];
// //                 if(projection_length > 0) //outlet
// // 		  rhs_i += A;
// // 		else
// 		  rhs_i -= A;
//
// // 	    }
// 	}
        KRATOS_CATCH ("")
    }
    //**************************************
    void CornerDectectionHelper (Geometry< Node < 3 > >& face_geometry,
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
        double acceptable_cos = cos (acceptable_angle);
        if (face_geometry[i1].Id() < face_geometry[i2].Id() ) //we do this to add the face ones
        {
            const array_1d<double, 3 > & neighb_normal = neighb[neighb_index].GetValue (NORMAL);
            double neighb_An = norm_2 (neighb_normal);
            double cos_normal = 1.0 / (An * neighb_An) * inner_prod (face_normal, neighb_normal);
            //if the angle is too big between the two normals then the edge in the middle is a corner
            if (cos_normal < acceptable_cos)
            {
                array_1d<double, 3 > edge = face_geometry[i2].Coordinates() - face_geometry[i1].Coordinates();
                double temp = norm_2 (edge);
                edge /= temp;
                int index1 = face_geometry[i1].FastGetSolutionStepValue (AUX_INDEX);
                int index2 = face_geometry[i2].FastGetSolutionStepValue (AUX_INDEX);
                edge_nodes[index1] += 1;
                edge_nodes[index2] += 1;
//                double sign1 = inner_prod (cornern_list[index1], edge);
                double sign1 = 0.0;
                for(unsigned int i = 0 ; i < edge.size() ; i++)
                    {sign1 += cornern_list[index1][i]*edge[i];}

                if (sign1 >= 0)
                {    for(unsigned int i = 0 ; i < edge.size() ; i++)
                        cornern_list[index1][i] += edge[i];
                }
                else
                {    for(unsigned int i = 0 ; i < edge.size() ; i++)
                        cornern_list[index1][i] -= edge[i];
                }

                double sign2 = inner_prod(cornern_list[index2], edge);
                if (sign2 >= 0)
                {    for(unsigned int i = 0 ; i < edge.size() ; i++)
                        cornern_list[index2][i] += edge[i];
                }
                else
                 {   for(unsigned int i = 0 ; i < edge.size() ; i++)
                        cornern_list[index2][i] -= edge[i];
                  }

            }
        }
    }
    //function to calculate the area normals
    void DetectEdges3D (ModelPart::ConditionsContainerType& rConditions)
    {
        KRATOS_TRY
        //calculate area normals face-by-face
        array_1d<double, 3 > area_normal;
        //(re)initialize normals
        unsigned int n_nodes = mNodalFlag.size();
        std::vector<unsigned int> temp_edge_nodes (n_nodes);
        CalcVectorType temp_cornern_list (n_nodes);
        for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
        {
            temp_edge_nodes[i_node] = 0.0;
            noalias (temp_cornern_list[i_node]) = ZeroVector (TDim);
        }
        //loop over all faces
        //            const double node_factor = 1.0 / TDim;
        for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
        {
            //get geometry data of the face
            Geometry<Node < 3 > >& face_geometry = cond_it->GetGeometry();
            //reference for area normal of the face
            const array_1d<double, 3 > & face_normal = cond_it->GetValue (NORMAL);
            double An = norm_2 (face_normal);
            unsigned int current_id = cond_it->Id();
            //slip condition
            if (cond_it->GetValue (IS_STRUCTURE) == 1.0) //this is a slip face --> now look for its neighbours
            {
                const WeakPointerVector<Condition>& neighb = cond_it->GetValue (NEIGHBOUR_CONDITIONS);
                //check for neighbour zero
                if (neighb[0].Id() != current_id) //check if the neighbour exists
                    CornerDectectionHelper (face_geometry, face_normal, An, neighb, 1, 2, 0, temp_edge_nodes, temp_cornern_list);
                //check for neighbour one
                if (neighb[1].Id() != current_id) //check if the neighbour exists
                    CornerDectectionHelper (face_geometry, face_normal, An, neighb, 2, 0, 1, temp_edge_nodes, temp_cornern_list);
                //check for neighbour two
                if (neighb[2].Id() != current_id) //check if the neighbour exists
                    CornerDectectionHelper (face_geometry, face_normal, An, neighb, 0, 1, 2, temp_edge_nodes, temp_cornern_list);
            }
        }
        //            ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
        //            mr_matrix_container.WriteVectorToDatabase(ACCELERATION, temp_cornern_list, rNodes);
        //fill the list of edge_nodes
        std::vector<unsigned int> tempmedge_nodes;
        std::vector< array_1d<double,TDim> > tempmedge_nodes_direction;
        std::vector<unsigned int> tempmcorner_nodes;
        for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
        {
            if (temp_edge_nodes[i_node] == 2) //node is a edge_node
            {
                tempmedge_nodes.push_back (i_node);
                array_1d<double, TDim>& node_edge = temp_cornern_list[i_node];
                node_edge /= norm_2 (node_edge);
                tempmedge_nodes_direction.push_back (node_edge);
            }
            else if (temp_edge_nodes[i_node] > 2)
                tempmcorner_nodes.push_back (i_node);
        }
        medge_nodes.resize (tempmedge_nodes.size(),false);
        medge_nodes_direction.resize (tempmedge_nodes_direction.size(),false);
        mcorner_nodes.resize (tempmcorner_nodes.size(),false);
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int> (tempmedge_nodes.size() ); i++)
        {
            medge_nodes[i] = tempmedge_nodes[i];
            medge_nodes_direction[i] = tempmedge_nodes_direction[i];
        }
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int> (tempmcorner_nodes.size() ); i++)
        {
            mcorner_nodes[i] = tempmcorner_nodes[i];
        }
        for (unsigned int i = 0; i < mcorner_nodes.size(); i++)
        {
            KRATOS_WATCH (mcorner_nodes[i]);
        }
        KRATOS_CATCH ("")
    }
//     double ComputePorosityCoefficient(const double& viscosity, const double& vel_norm, const double& eps, const double& d)
//     {
//         //             const double d = 0.01; //to be changed
//         double linear;
//         double non_linear;
//         if (eps < 1.0)
//         {
//             double k_inv = 150.0 * (1.0 - eps)*(1.0 - eps) / (eps * eps * eps * d * d);
//             linear = eps * viscosity * k_inv;
//             non_linear = (1.75 * vel_norm) * sqrt(k_inv / (150.0 * eps));
//             //             double linear = viscosity * k_inv;
//             //             double non_linear = (1.75 * vel_norm / eps) * sqrt(k_inv / (150.0 * eps));
//         }
//         else
//         {
//             linear = 0.0;
//             non_linear = 0.0;
//         }
//         return linear + non_linear;
//     }
    double ComputePorosityCoefficient (const double& vel_norm, const double& eps, const double& a, const double& b)
    {
        double linear;
        double non_linear;
        linear = eps * a;
        non_linear =  eps * b * vel_norm;
        return linear + non_linear;
    }
    void LaplacianSmooth (ValuesVectorType& to_be_smoothed, ValuesVectorType& aux)
    {
        ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
        int n_nodes = rNodes.size();
        #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            double dist = mdistances[i_node];
            double correction = 0.0;
            const double& origin_i = to_be_smoothed[i_node];
            if (dist <= 0.0) //node is inside domain ---- if outside do nothing
            {
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                    const double& origin_j = to_be_smoothed[j_neighbour];
                    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues() [csr_index];
                    double l_ikjk;
                    edge_ij.CalculateScalarLaplacian (l_ikjk);
                    correction += l_ikjk * (origin_j - origin_i);
                }
            }
            aux[i_node] = origin_i - correction;
        }
        #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; i_node++)
            to_be_smoothed[i_node] = aux[i_node];
    }

    void ComputeWallResistance (
        const CalcVectorType& vel,
	ValuesVectorType& diag_stiffness
//         CalcVectorType& rhs
    )
    {
        //parameters:
      //        double k = 0.41;
      //        double B = 5.1;
//         double density = mRho;
//        double toll = 1e-6;
        double ym = mY_wall; //0.0825877; //0.0093823
//        double y_plus_incercept = 10.9931899;
//        unsigned int itmax = 100;
        if (mViscosity[0] == 0)
            KRATOS_THROW_ERROR (std::logic_error, "it is not possible to use the wall law with 0 viscosity", "");
        /*        //slip condition
        int slip_size = mSlipBoundaryList.size();
        #pragma omp parallel for firstprivate(slip_size,B,toll,ym,y_plus_incercept,itmax)
        for (int i_slip = 0; i_slip < slip_size; i_slip++)
        {
            unsigned int i_node = mSlipBoundaryList[i_slip];
            double dist = mdistances[i_node];
            if (dist <= 0.0)
            {
                double nu = mViscosity[i_node];
                //array_1d<double, TDim>& rhs_i = rhs[i_node];
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
                double mod_uthaw = sqrt(mod_vel * nu / ym);
                double y_plus = ym * mod_uthaw / nu;
                if (y_plus > y_plus_incercept)
                {
                    //begin cicle to calculate the real u_thaw's module:
                    unsigned int it = 0;
                    double dx = 1e10;
                    //                        KRATOS_WATCH(fabs(dx));
                    while ( (fabs(dx) > toll * mod_uthaw) && (it < itmax) )
                    {
                        double a = 1.0 / k;
                        double temp = a * log(ym * mod_uthaw / nu) + B;
                        double y = mod_uthaw * (temp) - mod_vel;
                        double y1 = temp + a;
                        dx = y / y1;
                        mod_uthaw -= dx;
                        it = it + 1;
                    }
                    if (it == itmax)
                        std::cout << "attention max number of iterations exceeded in wall law computation" << std::endl;
                }
                        double tau = mod_uthaw * mod_uthaw  ;
                        tau *= mWallReductionFactor[i_node];
		if (mod_vel > 1e-9)
                            diag_stiffness[i_node] = tau * area / mod_vel;*/

        /*		        int slip_size = mSlipBoundaryList.size();
        #pragma omp parallel for firstprivate(slip_size,B,toll,ym,y_plus_incercept,itmax)
        for (int i_slip = 0; i_slip < slip_size; i_slip++)
        {
            unsigned int i_node = mSlipBoundaryList[i_slip];
            double dist = mdistances[i_node];
            if (dist <= 0.0)
            {
                double nu = mViscosity[i_node];
                //array_1d<double, TDim>& rhs_i = rhs[i_node];
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
                mod_vel = sqrt (mod_vel);
                area = sqrt (area);
                diag_stiffness[i_node] = area * mod_vel /pow(1.0/k*log(100) + B,2) * mWallReductionFactor[ i_node ];
            }
            else
                diag_stiffness[i_node] = 0.0;
        }*/

		//slip condition
        int slip_size = mSlipBoundaryList.size();
#pragma omp parallel for firstprivate(slip_size,ym)
        for (int i_slip = 0; i_slip < slip_size; i_slip++)
        {
            unsigned int i_node = mSlipBoundaryList[i_slip];
            double dist = mdistances[i_node];
            if (dist <= 0.0)
            {
                double nu = mMolecularViscosity; //mViscosity[i_node];
                //array_1d<double, TDim>& rhs_i = rhs[i_node];
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
                mod_vel = sqrt (mod_vel);
                area = sqrt (area);
				//the 0.1 is such that the dissipation is as for the linear case for a velocity of 10m/s
                diag_stiffness[i_node] = area * nu * mod_vel/ (ym ) * mWallReductionFactor[ i_node ] ;
        }
            else
	      {
                diag_stiffness[i_node] = 0.0 ;
	      }
    }


//         //apply higher resistance normally to the edges
//         int edge_size = medge_nodes_direction.size();
// 	#pragma omp parallel for firstprivate(edge_size)
// 	for (int i = 0; i < edge_size; i++)
// 	{
// 	    int i_node = medge_nodes[i];
// 	    double dist = mdistances[i_node];
//
// 	    if(dist <= 0.0)
// 	    {
// 		  double nu = mViscosity[i_node];
// 		  const array_1d<double, TDim>& an_i = mSlipNormal[i_node];
//
// 		  //compute the modulus of the velocity
// 		  double area = 0.0;
// 		  for (unsigned int comp = 0; comp < TDim; comp++)
// 		  {
// 		      area += an_i[comp] * an_i[comp];
// 		  }
// 		  area = sqrt (area);
//
// 		  diag_stiffness[i_node] += area * nu  / (ym ) ;
//
// 	    }
// 	}
//
// 	int corner_size = mcorner_nodes.size();
//         for (int i = 0; i < corner_size; i++)
//         {
//             int i_node = mcorner_nodes[i];
// 	    double nu = mViscosity[i_node];
//             mWallReductionFactor[i_node] = mcorner_coefficient; //50.0;
// 	    const double m = mr_matrix_container.GetLumpedMass()[i_node];
// 	    diag_stiffness[i_node] += 100.0*m * nu / (ym ) ;
//         }



    }

    void ApplySmagorinsky3D (double MolecularViscosity, double Cs)
    {
      KRATOS_TRY
      ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
      //calculating the RHS
        array_1d<double, TDim> grad_vx;
        array_1d<double, TDim> grad_vy;
	array_1d<double, TDim> grad_vz;
	int n_nodes = rNodes.size();
        mr_matrix_container.FillVectorFromDatabase (VELOCITY, mvel_n1, rNodes);
        array_1d<double, TDim> stab_high;
        #pragma omp parallel for private(grad_vx,grad_vy,grad_vz)
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
	  //set to zero the gradients
	  for (unsigned int comp = 0; comp < TDim; comp++)
	  {
                    grad_vx[comp] = 0.0 ;
		    grad_vy[comp] = 0.0 ;
		    grad_vz[comp] = 0.0 ;
	  }
	  //compute node by node the gradients
	  const array_1d<double, TDim>& U_i = mvel_n1[i_node];
	  const double h = mHmin[i_node];
            const double m_inv = mr_matrix_container.GetInvertedMass() [i_node];
            for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
                {
                unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                    const array_1d<double, TDim>& U_j = mvel_n1[j_neighbour];
                CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues() [csr_index];
                edge_ij.Add_grad_p (grad_vx, U_i[0], U_j[0]);
                edge_ij.Add_grad_p (grad_vy, U_i[1], U_j[1]);
                edge_ij.Add_grad_p (grad_vz, U_i[2], U_j[2]);
		}
	   //finalize computation of the gradients
	  //set to zero the gradients
	  for (unsigned int comp = 0; comp < TDim; comp++)
	  {
                    grad_vx[comp] *= m_inv ;
		    grad_vy[comp] *= m_inv ;
		    grad_vz[comp] *= m_inv ;
	  }
	  //symmetrize and multiply by 2
            grad_vx[0] *= 2.0;
	    grad_vy[1] *= 2.0;
			if(TDim > 2)
		        grad_vz[2] *= 2.0;
            grad_vx[1] += grad_vy[0];
			if(TDim > 2)
				grad_vx[2] += grad_vz[0];
			if(TDim > 2)
	            grad_vy[2] += grad_vz[1];
            grad_vy[0] += grad_vx[1];
            grad_vz[0] += grad_vx[2];
            grad_vz[1] += grad_vy[2];


	  //compute smagorinsky term
	  double aux = 0.0;
	  for (unsigned int comp = 0; comp < TDim; comp++)
	  {
                    aux += grad_vx[comp] * grad_vx[comp] ;
		    aux += grad_vy[comp] * grad_vy[comp] ;
		    aux += grad_vz[comp] * grad_vz[comp] ;
            }
            aux *= 0.5;
            if (aux < 0.0 ) aux=0.0;
            double turbulent_viscosity = Cs*h*h*sqrt (aux) /**MolecularViscosity*/;
// 	  KRATOS_WATCH(aux);
// 	  KRATOS_WATCH(turbulent_viscosity);
	  mViscosity[i_node] = turbulent_viscosity + MolecularViscosity;
	}
        mr_matrix_container.WriteScalarToDatabase (VISCOSITY, mViscosity, rNodes);
        KRATOS_CATCH ("");
    }

    void Add_Effective_Inverse_Multiply (
        CalcVectorType& destination,
        const CalcVectorType& origin1,
        const double value,
        const ValuesVectorType& mass,
	const ValuesVectorType& diag_stiffness,
        const CalcVectorType& origin
    )
    {
        KRATOS_TRY
        int loop_size = destination.size();
        #pragma omp parallel for
        for (int i_node = 0; i_node < loop_size; i_node++)
        {
            array_1d<double, TDim>& dest = destination[i_node];
            const double m = mass[i_node];
	    const double d = diag_stiffness[i_node];
            const array_1d<double, TDim>& origin_vec1 = origin1[i_node];
            const array_1d<double, TDim>& origin_value = origin[i_node];

            for (unsigned int comp = 0; comp < TDim; comp++)
                dest[comp] = value / (m + value*d) * ( m/value * origin_vec1[comp] +  origin_value[comp] );
        }
        KRATOS_CATCH ("")
    }

    void ComputeConvectiveProjection (
        CalcVectorType& mPiConvection,
        const ValuesVectorType& mphi_n1,
        const ValuesVectorType& mEps,
        const CalcVectorType& mvel_n1
    )
    {
        int n_nodes = mPiConvection.size();
        //calculating the convective projection
        array_1d<double, TDim> a_i;
        array_1d<double, TDim> a_j;
        #pragma omp parallel for  private(a_i,a_j)
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            array_1d<double, TDim>& pi_i = mPiConvection[i_node];
            // 		    setting to zero the projection
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                pi_i[l_comp] = 0.0;
            /*		    if (active_nodes[i_node] != 0.0)
                                {*/
            const double& phi_i = mphi_n1[i_node];
            noalias (a_i) = mvel_n1[i_node];
            a_i /= mEps[i_node];
            // 			  loop to all the edges surrounding node I
            for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
            {
                unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                noalias (a_j) = mvel_n1[j_neighbour];
                a_j /= mEps[j_neighbour];
                const double& phi_j = mphi_n1[j_neighbour];
                CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues() [csr_index];
                edge_ij.Add_grad_p (pi_i, phi_i, phi_j);
// 		if(i_node == 3255)
// 		{
// 		    KRATOS_WATCH(j_neighbour)
// 		    KRATOS_WATCH(pi_i)
// 		    KRATOS_WATCH(mEps[i_node])
// 		    KRATOS_WATCH(mEps[j_neighbour])
// 		    KRATOS_WATCH(phi_i)
// 		    KRATOS_WATCH(phi_j)
// 		    KRATOS_WATCH(a_i)
// 		    KRATOS_WATCH(a_j)
// 		    KRATOS_WATCH(mr_matrix_container.GetInvertedMass()[i_node])
// 		    KRATOS_WATCH(edge_ij.Ni_DNj)
//
// 		}
            }
            // 			  apply inverted mass matrix
            const double m_inv = mr_matrix_container.GetInvertedMass() [i_node];
            for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                pi_i[l_comp] *= m_inv;
// 	    std::cout << i_node << " " << pi_i << " " << mvel_n1[i_node] << " " << phi_i <<std::endl;
//             for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
//                 if(std::isnan(pi_i[l_comp]))
// 		  KRATOS_WATCH(m_inv);
            // 		    }
        }
    }

    void ComputeLimitor (
        CalcVectorType& mPiConvection,
        const ValuesVectorType& mphi_n1,
        ValuesVectorType& mBeta,
        const CalcVectorType& mvel_n1,
        const CalcVectorType& mEdgeDimensions
    )
    {
        int n_nodes = mPiConvection.size();
        #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            const array_1d<double, TDim>& pi_i = mPiConvection[i_node];
            const double& p_i = mphi_n1[i_node];
            double& beta_i = mBeta[i_node];
            beta_i = 0.0;
            double n = 0.0;
            for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex() [i_node]; csr_index != mr_matrix_container.GetRowStartIndex() [i_node + 1]; csr_index++)
            {
                unsigned int j_neighbour = mr_matrix_container.GetColumnIndex() [csr_index];
                const double& p_j = mphi_n1[j_neighbour];
                const array_1d<double, TDim>& l_k = mEdgeDimensions[csr_index];
                const array_1d<double, TDim>& pi_j = mPiConvection[j_neighbour];
                //                                 double proj = 0.0;
                //                                 for (unsigned int comp = 0; comp < TDim; comp++)
                //                                      proj += 0.5*l_k[comp]*(pi_i[comp]+pi_j[comp]);
                //                                 double beta = fabs((p_i - p_j - proj)/(fabs(p_i-p_j)+fabs(proj)+1e-4));
                double proj = 0.0;
                for (unsigned int comp = 0; comp < TDim; comp++)
                    proj += 0.5 * l_k[comp]* (pi_i[comp] + pi_j[comp]);
                // 							proj += dir[comp]*pi_i[comp];
                double numerator = fabs (fabs (p_j - p_i) - fabs (proj) );
                double denom = fabs (fabs (p_j - p_i) + 1e-6);
                beta_i += numerator / denom;
                n += 1.0;
            }
            beta_i /= n;
            if (beta_i > 1.0)
                beta_i = 1.0;
        }
    }
};
} //namespace Kratos
#undef SYMM_PRESS
#endif //KRATOS_EDGEBASED_LEVELSET_SUBSTEP_FLUID_SOLVER_H_INCLUDED defined
