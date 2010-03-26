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


#if !defined(KRATOS_EDGEBASED_LEVELSET_FLUID_SOLVER_H_INCLUDED)
#define  KRATOS_EDGEBASED_LEVELSET_FLUID_SOLVER_H_INCLUDED

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
    class EdgeBasedLevelSet {
    public:
        //name for the self defined structure
        typedef EdgesStructureType<TDim> CSR_Tuple;
        typedef std::vector<CSR_Tuple> EdgesVectorType;

        //name for row start and column index vectors
        typedef std::vector<unsigned int> IndicesVectorType;
        //defining matrix type for test calculations
        typedef std::vector< array_1d<double, TDim> > CalcVectorType;
        //defining type for local storage of nodal values
        typedef std::vector<double> ValuesVectorType;

        //defining types for matrix operations
        typedef typename TSparseSpace::MatrixType TSystemMatrixType;
        typedef typename TSparseSpace::VectorType TSystemVectorType;

        //constructor and destructor

        EdgeBasedLevelSet(MatrixContainer& mr_matrix_container,
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
            medge_detection_angle(edge_detection_angle),
          mstabdt_pressure_factor(stabdt_pressure_factor),
            mstabdt_convection_factor(stabdt_convection_factor),
            mtau2_factor(tau2_factor),
            massume_constant_dp(assume_constant_dp)

        {
            //to be changed
            mViscosity = viscosity;

            noalias(mBodyForce) = body_force;
            mRho = density;

            mdelta_t_avg = 1000.0;

            max_dt = 1.0;

            muse_mass_correction = use_mass_correction;

//            for (unsigned int i = 0; i < TDim; i++) mBodyForce[i] = 0;
//            mBodyForce[1] = -9.81;
//
//            mRho = 1000.0;




        };

        ~EdgeBasedLevelSet() {
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
            mdistances.resize(n_nodes);

            mTauPressure.resize(n_nodes);
            mTauConvection.resize(n_nodes);
            mTau2.resize(n_nodes);
            mPi.resize(n_nodes);
            mXi.resize(n_nodes);
            mx.resize(n_nodes);

            mEdgeDimensions.resize(n_edges);

            //convection variables
            mBeta.resize(n_nodes);
            mPiConvection.resize(n_nodes);
            mphi_n.resize(n_nodes);
            mphi_n1.resize(n_nodes);

            mEps.resize(n_nodes);
	   mD.resize(n_nodes);

            mdiv_error.resize(n_nodes);
            mr_matrix_container.SetToZero(mdiv_error);


            ValuesVectorType external_pressure;
            external_pressure.resize(n_nodes);

            //read velocity and pressure data from Kratos
            mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, mr_model_part.Nodes());
            mr_matrix_container.FillOldScalarFromDatabase(PRESSURE, mPn, mr_model_part.Nodes());
            mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mvel_n, mr_model_part.Nodes());
            mr_matrix_container.FillCoordinatesFromDatabase(mx, mr_model_part.Nodes());
            //set flag for first time step
            mFirstStep = true;

            //loop to categorize boundary nodes
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


                    mFixedVelocities.push_back(index);
                    mFixedVelocitiesValues.push_back(mvel_n1[index]);
                }

                if (inode->IsFixed(PRESSURE)) {
                    mPressureOutletList.push_back(index);
                    mPressureOutlet.push_back(external_pressure[index]);
                }
            }

            //compute slip normals and fill SlipList
            CalculateNormals(mr_model_part.Conditions());
            mr_matrix_container.WriteVectorToDatabase(NORMAL, mSlipNormal, mr_model_part.Nodes());

            if(TDim == 3)
                DetectEdges3D(mr_model_part.Conditions());

            //print number of nodes corresponding to the different types of boundary conditions
            //				KRATOS_WATCH(mFixedVelocities.size())
            //				KRATOS_WATCH(mPressureOutletList.size())
            //				KRATOS_WATCH(mSlipBoundaryList.size())


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


            //            //set the pressure projection to the body force value
            //            array_1d<double,3> temp = mRho * mBodyForce;
            //            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
            //                    inode != mr_model_part.NodesEnd();
            //                    inode++)
            //                inode->FastGetSolutionStepValue(PRESS_PROJ) = temp;

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
//            mr_matrix_container.FillVectorFromDatabase(PRESS_PROJ, mXi, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(POROSITY, mEps, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(DIAMETER, mD, mr_model_part.Nodes());


//            double delta_t_i = delta_t;

            //*******************
            //loop over all nodes
            double n_nodes = mvel_n1.size();
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                const array_1d<double, TDim>& v_i = mvel_n1[i_node];
                const double havg_i = mHavg[i_node];
                const double hmin_i = mHmin[i_node];
                const double eps_i = mEps[i_node];
	       const double d_i = mD[i_node];

                double vel_norm = norm_2(v_i);

                double porosity_coefficient = ComputePorosityCoefficient(mViscosity, vel_norm, eps_i, d_i);
                vel_norm /= eps_i;

                //use CFL condition to compute time step size
                double delta_t_i = CFLNumber * 1.0 / (2.0 * vel_norm /hmin_i + 4.0 * mViscosity / (hmin_i * hmin_i) + porosity_coefficient);
                double delta_t_i_avg = 1.0 / (2.0 * vel_norm /havg_i + 4.0 * mViscosity / (havg_i * havg_i) + porosity_coefficient);

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
                    v_diff_norm /= eps_i;
                    double delta_t_j = CFLNumber * 1.0 / (2.0 * v_diff_norm /hmin_i + 4.0 * mViscosity / (hmin_i * hmin_i));

                    if (delta_t_j < delta_t_i)
                            delta_t_i = delta_t_j;
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

            mr_matrix_container.FillScalarFromDatabase(DISTANCE, mdistances, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(DIAMETER, mD, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(POROSITY, mEps, mr_model_part.Nodes());

            //read time step size from Kratos
            ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
            double delta_t = CurrentProcessInfo[DELTA_TIME];
           //read the prescribed values of velocity
//            int fixed_size = mD.size();
// //            #pragma omp parallel for firstprivate(fixed_size)
//            for (int i_d = 0; i_d < fixed_size; i_d++)
//            {
//                    KRATOS_WATCH(mD[i_d]);
//            }


//            //read the prescribed values of velocity
//            int fixed_size = mFixedVelocities.size();
////            #pragma omp parallel for firstprivate(fixed_size)
//            for (int i_velocity = 0; i_velocity < fixed_size; i_velocity++)
//            {
//                unsigned int i_node = mFixedVelocities[i_velocity];
//                     array_1d<double, TDim>& u_i_fix = mFixedVelocitiesValues[i_velocity];
//                    const array_1d<double, TDim>& u_i = mvel_n1[i_node];
//                    KRATOS_WATCH(mvel_n1[i_node-1]);
//                    KRATOS_WATCH(mvel_n1[i_node]);
//                    KRATOS_WATCH(mvel_n1[i_node+1]);
//                    for (unsigned int comp = 0; comp < TDim; comp++)
//                        u_i_fix[comp] = u_i[comp];
//            }
//            KRATOS_WATCH("AAAAAAAAAAAAAAAAA")
//            mFixedVelocities.resize(0);
//            mFixedVelocitiesValues.resize(0);
//            mFixedVelocities.clear();
//            mFixedVelocitiesValues.clear();
//
//            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
//                    inode != mr_model_part.NodesEnd();
//                    inode++) {
//                int index = inode->FastGetSolutionStepValue(AUX_INDEX);
//                if (inode->IsFixed(VELOCITY_X)) //note that the variables can be either all fixed or no one fixed
//                {
//                    KRATOS_WATCH(index);
//                    KRATOS_WATCH(inode->FastGetSolutionStepValue(VELOCITY));
//                    mFixedVelocities.push_back(index);
//                    mFixedVelocitiesValues.push_back(mvel_n1[index]);
//
//
//                }
//            }
//
//            int fixed_size = mFixedVelocities.size();
////            #pragma omp parallel for firstprivate(fixed_size)
//            for (int i_velocity = 0; i_velocity < fixed_size; i_velocity++)
//            {
////                                    mFixedVelocitiesValues[i_velocity][0] = 10.0;
////                    mFixedVelocitiesValues[i_velocity][1] = 10.0;
////                    mFixedVelocitiesValues[i_velocity][2] = 10.0;
//                KRATOS_WATCH( mFixedVelocitiesValues[i_velocity] );
//            }



            //compute intrinsic time
//            double time_inv = 1.0 / delta_t;
            double time_inv_avg = 1.0/mdelta_t_avg;

            const double stabdt_pressure_factor  = mstabdt_pressure_factor;
            const double stabdt_convection_factor  = mstabdt_convection_factor;
	    const double tau2_factor = mtau2_factor;

//            const double max_dt_inv = 1.0 / max_dt;

//            double max_dt_inv_coeff = 0.1 / max_dt;
//            const double max_dt_coeff = 10.0*max_dt;

            #pragma omp parallel for firstprivate(time_inv_avg,stabdt_pressure_factor,stabdt_convection_factor,tau2_factor)
            for (int i_node = 0; i_node < n_nodes; i_node++) {
                // 					double& h_i = mHavg[i_node];
                double& h_avg_i = mHavg[i_node];
                array_1d<double, TDim>& a_i = mvel_n1[i_node];
                const double nu_i = mViscosity;
                const double eps_i = mEps[i_node];
                const double d_i = mD[i_node];

                double vel_norm = norm_2(a_i);

                double porosity_coefficient = ComputePorosityCoefficient(mViscosity, vel_norm, eps_i, d_i);
                vel_norm /= eps_i;

//                double tau = 1.0 / (2.0 * vel_norm / h_avg_i + time_inv_avg + (4.0*nu_i) / (h_avg_i * h_avg_i) + porosity_coefficient);
//                double denom = (2.0 * vel_norm / h_avg_i + (4.0*nu_i) / (h_avg_i * h_avg_i) + porosity_coefficient);
//                double tau = 0.0;
//                if(denom > max_dt_inv_coeff)
//                    tau = max_dt_coeff;
//                else
//                    tau = 1.0/denom;

//                double tau = 1.0 / (2.0 * vel_norm / h_avg_i + max_dt_inv + (4.0*nu_i) / (h_avg_i * h_avg_i) + porosity_coefficient);
                 double tau = 1.0 / (2.0 * vel_norm / h_avg_i + stabdt_pressure_factor*time_inv_avg + (4.0*nu_i) / (h_avg_i * h_avg_i) + porosity_coefficient);
//                 double tau = 1.0 / (2.0 * vel_norm / h_avg_i + 0.01*time_inv_avg + (4.0*nu_i) / (h_avg_i * h_avg_i) + porosity_coefficient);
                double tau_conv = 1.0 / (2.0 * vel_norm / h_avg_i + stabdt_convection_factor*time_inv_avg + (4.0*nu_i) / (h_avg_i * h_avg_i) + porosity_coefficient);
                mTauPressure[i_node] = tau;
                mTauConvection[i_node] = tau_conv;

                 mTau2[i_node] = (mViscosity + h_avg_i*vel_norm*0.5)*mtau2_factor;

//                mTauPressure[i_node] = 1.0 / (2.0 * vel_norm / mHavg[i_node] + (4.0*nu_i) / (mHavg[i_node] * mHavg[i_node]));
//                mTauConvection[i_node] = 1.0 / (2.0 * vel_norm / h_i + time_inv + (4.0*nu_i) / (h_i * h_i));

////                mTauPressure[i_node] = 1.0 / (2.0 * vel_norm / h_i + 0.01 * time_inv + 4.0 * nu_i / (h_i * h_i));
////                //                 mTauPressure[i_node] = delta_t;
////                mTauConvection[i_node] = 1.0 / (2.0 * vel_norm / h_i + 0.01 * time_inv + 4.0 * nu_i / (h_i * h_i));

                //					if (mTauPressure[i_node] < delta_t)
                //						mTauPressure[i_node] =  delta_t;
                //					else if(mTauPressure[i_node] > 100.0*delta_t)
                //						mTauPressure[i_node] = 100.0*delta_t;
            }

////            //the tau is set to 1/dt on the corner nodes
////            //apply conditions on corners
////            int corner_size = mcorner_nodes.size();
////            for (int i = 0; i < corner_size; i++)
////            {
////                int i_node = mcorner_nodes[i];
////                mTauPressure[i_node] = mdelta_t_avg;
////                mTauConvection[i_node] = mdelta_t_avg;
////            }

            //laplacian smoothing on the taus
            //note here that we use TauConvection as a temporary vector
//            LaplacianSmooth(mTauConvection, mTauPressure);
//            LaplacianSmooth(mTauPressure, mTauConvection);
//            mr_matrix_container.AssignVectorToVector(mTauPressure, mTauConvection);



            //calculating the convective projection
            #pragma omp parallel for
            for (int i_node = 0; i_node < n_nodes; i_node++) {
                array_1d<double, TDim>& pi_i = mPi[i_node]; //******************

                //setting to zero
                for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                    pi_i[l_comp] = 0.0;

                array_1d<double, TDim> a_i = mvel_n1[i_node];
                const array_1d<double, TDim>& U_i = mvel_n1[i_node];
//                 const double& p_i = mPn1[i_node];
		const double& eps_i = mEps[i_node];

		a_i /= eps_i;

                //const double& p_i = pressure[i_node];

                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    array_1d<double, TDim> a_j = mvel_n1[j_neighbour];
                    const array_1d<double, TDim>& U_j = mvel_n1[j_neighbour];
                    const double& eps_j = mEps[j_neighbour];

		    a_j /= eps_j;

                    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

                    edge_ij.Add_ConvectiveContribution(pi_i, a_i, U_i, a_j, U_j);

//                    edge_ij.Add_grad_p(pi_i, p_i, p_j);
                }

                const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];

                for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                    pi_i[l_comp] *= m_inv;
            }

            //calculating limitor
            #pragma omp parallel for
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {

                    const array_1d<double, TDim>& pi_i = mXi[i_node];
                    const double& p_i = mPn[i_node];
                    double& beta_i = mBeta[i_node];
                    beta_i = 0.0;

                        for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                        {
                            unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];

                                const double& p_j = mPn[j_neighbour];
                                const array_1d<double, TDim>& l_k = mEdgeDimensions[csr_index];
                                const array_1d<double, TDim>& pi_j = mXi[j_neighbour];

                                double proj = 0.0;
                                for (unsigned int comp = 0; comp < TDim; comp++)
                                     proj += 0.5*l_k[comp]*(pi_i[comp]+pi_j[comp]);
                                double beta = fabs((p_i - p_j - proj)/(fabs(p_i-p_j)+fabs(proj)+1e-4));

                                if(beta_i < beta)
                                    beta_i = beta;
                        }

            }


            mr_matrix_container.AssignVectorToVector(mvel_n, mWork); //mWork = mvel_n

            //first step of Runge Kutta
            mr_matrix_container.AssignVectorToVector(mvel_n, mvel_n1); //mvel_n1 = mvel_n

            mr_matrix_container.SetToZero(rhs);
            CalculateRHS(mvel_n1, mPn, mvel_n1, rhs);
            mr_matrix_container.Add_Minv_value(mWork, mWork, delta_t / 6.0, mr_matrix_container.GetInvertedMass(), rhs);
            mr_matrix_container.Add_Minv_value(mvel_n1, mvel_n, 0.5 * delta_t, mr_matrix_container.GetInvertedMass(), rhs);
            ApplyVelocityBC(mvel_n1);

            //second step
            mr_matrix_container.SetToZero(rhs);
            CalculateRHS(mvel_n1, mPn, mvel_n1, rhs);
            mr_matrix_container.Add_Minv_value(mWork, mWork, delta_t / 3.0, mr_matrix_container.GetInvertedMass(), rhs);
            mr_matrix_container.Add_Minv_value(mvel_n1, mvel_n, 0.5 * delta_t, mr_matrix_container.GetInvertedMass(), rhs);
            ApplyVelocityBC(mvel_n1);

            //third step
            mr_matrix_container.SetToZero(rhs);
            CalculateRHS(mvel_n1, mPn, mvel_n1, rhs);
            mr_matrix_container.Add_Minv_value(mWork, mWork, delta_t / 3.0, mr_matrix_container.GetInvertedMass(), rhs);
            mr_matrix_container.Add_Minv_value(mvel_n1, mvel_n, delta_t, mr_matrix_container.GetInvertedMass(), rhs);
            ApplyVelocityBC(mvel_n1);

            //fourth step
            mr_matrix_container.SetToZero(rhs);
            CalculateRHS(mvel_n1, mPn, mvel_n1, rhs);
            mr_matrix_container.Add_Minv_value(mWork, mWork, delta_t / 6.0, mr_matrix_container.GetInvertedMass(), rhs);

            //compute right-hand side
            mr_matrix_container.AssignVectorToVector(mWork, mvel_n1);
            ApplyVelocityBC(mvel_n1);




            KRATOS_CATCH("")
        }



        //*********************************************************************
        //function to calculate right-hand side of fractional momentum equation

        void CalculateRHS(
                const CalcVectorType& vel,
                const ValuesVectorType& pressure,
                const CalcVectorType& convective_velocity,
                CalcVectorType& rhs) {
            KRATOS_TRY

                    int n_nodes = vel.size();






            //perform MPI syncronization

            //calculating the RHS
            array_1d<double, TDim> stab_low;
            array_1d<double, TDim> stab_high;
            const double nu_i = mViscosity;
            const double nu_j = mViscosity;
            double inverse_rho = 1.0 / mRho;
#pragma omp parallel for private(stab_low,stab_high)
            for (int i_node = 0; i_node < n_nodes; i_node++) {
                double dist = mdistances[i_node];
                if (dist <= 0.0) //node is inside domain ---- if outside do nothing
                {
                    array_1d<double, TDim>& rhs_i = rhs[i_node];
                    const array_1d<double, TDim>& f_i = mBodyForce;
                    array_1d<double, TDim> a_i = convective_velocity[i_node];
                    const double& beta_i = mBeta[i_node];
                    const array_1d<double, TDim>& U_i = vel[i_node];
                    const array_1d<double, TDim>& pi_i = mPi[i_node];
                    const double& p_i = pressure[i_node];
                    const double& eps_i = mEps[i_node];
                    const double& d_i = mD[i_node];
// KRATOS_WATCH("before");
// KRATOS_WATCH(d_i);
                    const double& tau2_i = mTau2[i_node];

                    double edge_tau = mTauConvection[i_node];
                    a_i /= eps_i;

                    //double& h_i = mHmin[i_node];

                    //initializing with the external forces (e.g. gravity)
                    double& m_i = mr_matrix_container.GetLumpedMass()[i_node];
                    for (unsigned int comp = 0; comp < TDim; comp++)
                        rhs_i[comp] = m_i * eps_i * f_i[comp] ;

                    //applying the effect of the porosity
                    double porosity_coefficient = ComputePorosityCoefficient(mViscosity,norm_2(U_i),eps_i, d_i);
// KRATOS_WATCH("after");
// KRATOS_WATCH(d_i);
                    for (unsigned int comp = 0; comp < TDim; comp++)
                        rhs_i[comp] -= m_i * porosity_coefficient * U_i[comp];

                    //std::cout << i_node << "rhs =" << rhs_i << "after adding body force" << std::endl;
                    //convective term
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                            array_1d<double, TDim> a_j = convective_velocity[j_neighbour];
                            const array_1d<double, TDim>& U_j = vel[j_neighbour];
                            const array_1d<double, TDim>& pi_j = mPi[j_neighbour];
                            const double& p_j = pressure[j_neighbour];
                            const double& eps_j = mEps[j_neighbour];
                            const double& beta_j = mBeta[j_neighbour];
                            a_j /= eps_j;

                            CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

                            edge_ij.Sub_ConvectiveContribution(rhs_i, a_i, U_i, a_j, U_j);
                            //std::cout << i_node << "rhs =" << rhs_i << "after convective contrib" << std::endl;

                            //take care! we miss including a B.C.  for the external pressure
                            //							edge_ij.Add_Gp(rhs_i,p_i*inverse_rho,p_j*inverse_rho);
                            edge_ij.Sub_grad_p(rhs_i, p_i*inverse_rho*eps_i, p_j * inverse_rho*eps_i);
                            //                        edge_ij.Add_grad_p(rhs_i, p_i*inverse_rho, p_j * inverse_rho);
                            //std::cout << i_node << "rhs =" << rhs_i << "after Gp" << std::endl;

                            edge_ij.Sub_ViscousContribution(rhs_i, U_i, nu_i, U_j, nu_j);
                            //std::cout << i_node << "rhs =" << rhs_i << "after viscous" << std::endl;

                            //add stabilization
                            edge_ij.CalculateConvectionStabilization_LOW(stab_low, a_i, U_i, a_j, U_j);
//                            edge_ij.CalculateConvectionStabilization_LOW(stab_low, a_i, U_i,p_i, a_j, U_j,p_j);
                             edge_ij.CalculateConvectionStabilization_HIGH(stab_high, a_i, pi_i, a_j, pi_j);

//                            double beta = 1.0;
                            double beta = beta_i;
                            if(beta_j > beta)
                                beta = beta_j;
//                            beta = 1.0; 

//                            edge_ij.Sub_StabContribution(rhs_i, edge_tau*beta, 1.0, stab_low, stab_high);
                            edge_ij.Sub_StabContribution(rhs_i, edge_tau, (1.0-beta), stab_low, stab_high);
                            //std::cout << i_node << "rhs =" << rhs_i << "after stab" << std::endl;

                            //add tau2 term
                            boost::numeric::ublas::bounded_matrix<double,TDim,TDim>& LL = edge_ij.LaplacianIJ;
                            for (unsigned int k_comp = 0; k_comp < TDim; k_comp++)
                            {
                                double aaa = 0.0;
                                for (unsigned int m_comp = 0; m_comp < TDim; m_comp++)
                                    aaa +=  LL(k_comp,m_comp) * (U_i[m_comp] - U_i[m_comp]);
                                rhs_i[k_comp] -= tau2_i*aaa;
                            }


                    }

                    //                                                std::cout << i_node << "rhs =" << rhs_i << std::endl;
                }

            }

            //apply wall resistance
           if(mWallLawIsActive == true)
                    ComputeWallResistance(vel,rhs);

            //boundary integrals --> finishing the calculation of the pressure gradient
            //				int loop_size1 = mPressureOutletList.size();
            //				#pragma omp parallel for
            //				for (int i_pressure = 0; i_pressure < loop_size1; i_pressure++)
            //				{
            //					unsigned int i_node = mPressureOutletList[i_pressure];
            //					array_1d<double, TDim>& rhs_i = rhs[i_node];
            //					const double& p_ext_i = mPressureOutlet[i_pressure];
            //					const array_1d<double, TDim>& an_i = mPressureNormal[i_node];
            //
            // 					for (unsigned int comp = 0; comp < TDim; comp++)
            // 						rhs_i[comp] -= an_i[comp] *  p_ext_i;
            //				}

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
            mr_matrix_container.FillVectorFromDatabase(PRESS_PROJ, mXi, rNodes);
            mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes);
            //for (int i_node = 0; i_node < n_nodes; i_node++)
            //    std::cout << mvel_n1[i_node] << std::endl;

            //loop over all nodes
//            double rho_inv = 1.0 / mRho;
#pragma omp parallel for firstprivate(time_inv)
            for (int i_node = 0; i_node < n_nodes; i_node++) {

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
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    const double& p_j = mPn1[j_neighbour];
                    const double& p_old_j = mPn[j_neighbour];
                    const array_1d<double, TDim>& U_j_curr = mvel_n1[j_neighbour];
                    const array_1d<double, TDim>& xi_j = mXi[j_neighbour];
//                    const double& eps_j = mEps[j_neighbour];

                    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

#ifdef SYMM_PRESS
                    double edge_tau = 0.5 * (mTauPressure[i_node] + mTauPressure[j_neighbour]);
#else
                    double edge_tau = mTauPressure[i_node];
#endif
                    //     						double edge_tau = CalculateEdgeTau(time_inv,h_i,a_i,h_j,a_j);
                    //



                    //compute laplacian operator
                    double sum_l_ikjk;
                    edge_ij.CalculateScalarLaplacian(sum_l_ikjk);
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
                    edge_ij.Sub_D_v(rhs_i, U_i_curr*mRho, U_j_curr * mRho);
                    //       						edge_ij.Sub_D_v(rhs_i,a_i*rho_i,a_j*rho_i);

                    //high order stabilizing term
                    double temp = 0.0;
                    // 						edge_ij.Add_div_v(temp,mTauPressure[i_node]*xi_i,mTauPressure[j_neighbour]*xi_j);
                    edge_ij.Add_div_v(temp, xi_i, xi_j);
                    rhs_i += edge_tau * temp;

                    //assemble laplacian matrix
                    mL(i_node, j_neighbour) = sum_l_ikjk;
                    l_ii -= sum_l_ikjk;
                }

//                //area correction to prevent mass loss
//                rhs_i -= mdiv_error[i_node];

//                rhs_i += div_i * eps_i;
                mL(i_node, i_node) = l_ii;
            }

            if(muse_mass_correction == true)
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
            for (int i_node = 0; i_node < n_nodes; i_node++) {
                double L_diag = mL(i_node, i_node);
                if (fabs(L_diag) > fabs(max_diag)) max_diag = L_diag;
            }




            //respect pressure boundary conditions by penalization
//            double huge = max_diag * 1e6;
//            for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++) {
//                unsigned int i_node = mPressureOutletList[i_pressure];
//                mL(i_node, i_node) = huge;
//                rhs[i_node] = 0.0;
//            }
            for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++) {
                unsigned int i_node = mPressureOutletList[i_pressure];
                mL(i_node, i_node) = max_diag;
                rhs[i_node] = 0.0;
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    mL(i_node, j_neighbour) = 0.0;
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

            for (int i_node = 0; i_node < n_nodes; i_node++) {
                if (mdistances[i_node] >= 0) {
                    mL(i_node, i_node) = max_diag;
                    rhs[i_node] = 0.0;
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        mL(i_node, j_neighbour) = 0.0;
                    }
                } else {
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        if (mdistances[j_neighbour] >= 0)
                            mL(i_node, j_neighbour) = 0.0;
                    }
                }
            }





            //set starting vector for iterative solvers
            for (int i_node = 0; i_node < n_nodes; i_node++)
                dp[i_node] = 0.0;
            //KRATOS_WATCH(rhs);
            //solve linear equation system L dp = rhs
            pLinearSolver->Solve(mL, dp, rhs);
            KRATOS_WATCH(*pLinearSolver)


                    //update pressure
            for (int i_node = 0; i_node < n_nodes; i_node++)
                mPn1[i_node] += dp[i_node];
            //				for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++)
            //				{
            //					unsigned int i_node = mPressureOutletList[i_pressure];
            //					mPn1[i_node] = mPressureOutlet[i_pressure];
            //				}

            //write pressure and density to Kratos
            mr_matrix_container.WriteScalarToDatabase(PRESSURE, mPn1, rNodes);


            //compute pressure proj for the next step

            #pragma omp parallel for firstprivate(time_inv), private(work_array)
            for (int i_node = 0; i_node < n_nodes; i_node++) {
                array_1d<double, TDim>& xi_i = mXi[i_node];
                for (unsigned int comp = 0; comp < TDim; comp++)
                    xi_i[comp] = 0.0;

                double dist = mdistances[i_node];
                if (dist <= 0.0) //node is inside domain ---- if outside do nothing
                {

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
                double dist = mdistances[i_node];
                if (dist < 0.0) //node is inside domain ---- if outside do nothing
                {
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

                        // 							edge_ij.Sub_grad_p(correction,delta_p_i,delta_p_j);
                        edge_ij.Sub_grad_p(correction, delta_p_i, delta_p_j);
                        //                        edge_ij.Add_grad_p(correction, delta_p_i, delta_p_j);
                        //                                                        edge_ij.Add_Gp(correction,delta_p_i,delta_p_j);
                        // 							edge_ij.Sub_Gp(correction,delta_p_i,delta_p_j);
                    }
                    //compute prefactor
                    double coefficient = delta_t * m_inv;

                    //correct fractional momentum
                    for (unsigned int comp = 0; comp < TDim; comp++)
                        U_i_curr[comp] += coefficient * correction[comp];
                }
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
                    const double dist = mdistances[i_node];
                    double& div_i_err = mdiv_error[i_node];
                    div_i_err = 0.0;
                    if (dist < 0.0) //node is inside domain ---- if outside do nothing
                    {
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
            }

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
                    double dist = mdistances[i_node];

                    if(dist <= 0.0)
                    {
                         array_1d<double, TDim>& U_i = VelArray[i_node];
                         double temp=0.0;
                         for (unsigned int comp = 0; comp < TDim; comp++)
                            temp += U_i[comp] * direction[comp];

                         for (unsigned int comp = 0; comp < TDim; comp++)
                             U_i[comp] = direction[comp]*temp;
                    }
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
                double dist = mdistances[i_node];
                if(dist <= 0.0)
                {
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
            }

            //fixed condition
            int fixed_size = mFixedVelocities.size();
#pragma omp parallel for firstprivate(fixed_size)
            for (int i_velocity = 0; i_velocity < fixed_size; i_velocity++)
            {
                unsigned int i_node = mFixedVelocities[i_velocity];
                double dist = mdistances[i_node];
                if(dist <= 0.0)
                {
                    const array_1d<double, TDim>& u_i_fix = mFixedVelocitiesValues[i_velocity];
                    array_1d<double, TDim>& u_i = VelArray[i_node];

                    for (unsigned int comp = 0; comp < TDim; comp++)
                        u_i[comp] = u_i_fix[comp];
                }
            }



            KRATOS_CATCH("")
        }



        //********************************
        //function to compute coefficients

        void ExtrapolateValues(unsigned int extrapolation_layers)
        {
            KRATOS_TRY

                    typedef Node < 3 > PointType;
            typedef PointerVector<PointType > PointVector;
            typedef PointVector::iterator PointIterator;

            mr_matrix_container.FillScalarFromDatabase(DISTANCE, mdistances,mr_model_part.Nodes());
//            mr_matrix_container.FillVectorFromDatabase(PRESS_PROJ, mXi,mr_model_part.Nodes());
//
//            //loop on all the slip nodes and Set the pressure projection to -BodyForce if it has neighbours with distance greater than 0
//            int slip_size = mSlipBoundaryList.size();
//            #pragma omp parallel for firstprivate(slip_size)
//            for (int i_slip = 0; i_slip < slip_size; i_slip++)
//            {
//                unsigned int i_node = mSlipBoundaryList[i_slip];
//                double dist = mdistances[i_node];
//
//
//                if(dist <= 0.0)
//                {
//                    int nout = 0;
//                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
//                    {
//                        //get global index of neighbouring node j
//                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
//                        const double& dist_j  = mdistances[j_neighbour];
//
//                        if(dist_j > 0)
//                            nout++;
//                    }
//
//                    if(nout > 0) mXi[i_node] += mRho*mBodyForce;
//                }
//            }
//
//            mr_matrix_container.WriteVectorToDatabase(PRESS_PROJ, mXi,mr_model_part.Nodes());


            //reset is visited flag
            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                    inode != mr_model_part.NodesEnd();
                    inode++) {
                inode->GetValue(IS_VISITED) = 0;
            }


            //generate a container with the layers to be extrapolated
            std::vector< PointVector > layers(extrapolation_layers);

            //detect the nodes inside the fluid surface
            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                    inode != mr_model_part.NodesEnd();
                    inode++) {
                if (inode->FastGetSolutionStepValue(DISTANCE) <= 0.0) //candidates are only the ones inside the fluid domain
                {
                    WeakPointerVector< Node < 3 > >& neighb_nodes = inode->GetValue(NEIGHBOUR_NODES);
                    for (WeakPointerVector< Node < 3 > >::iterator i = neighb_nodes.begin(); i != neighb_nodes.end(); i++) {
                        if (i->FastGetSolutionStepValue(DISTANCE) > 0) //add the node as free surface if one of its neighb is outside
                        {
                            if (inode->GetValue(IS_VISITED) == 0) {
                                layers[0].push_back(*(inode.base()));
                                inode->GetValue(IS_VISITED) = 1;
                            }
                        }
                    }
                } else {
                    //set everything to zero
                    noalias(inode->FastGetSolutionStepValue(VELOCITY)) = ZeroVector(3);
                    inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
                    noalias(inode->FastGetSolutionStepValue(VELOCITY, 1)) = ZeroVector(3);
                    inode->FastGetSolutionStepValue(PRESSURE, 1) = 0.0;

                    noalias(inode->FastGetSolutionStepValue(PRESS_PROJ)) = ZeroVector(3);
                    noalias(inode->FastGetSolutionStepValue(PRESS_PROJ, 1)) = ZeroVector(3);

                }
            }



            //fill the following layers by neighbour relationships
            //each layer fills the following
            for (unsigned int il = 0; il < extrapolation_layers - 1; il++) {
                for (PointIterator iii = (layers[il]).begin(); iii != (layers[il]).end(); iii++) {
                    WeakPointerVector< Node < 3 > >& neighb_nodes = iii->GetValue(NEIGHBOUR_NODES);
                    for (WeakPointerVector< Node < 3 > >::iterator jjj = neighb_nodes.begin(); jjj != neighb_nodes.end(); jjj++) //destination = origin1 + value * Minv*origin
                    {

                        if (jjj->FastGetSolutionStepValue(DISTANCE) > 0 &&
                                jjj->GetValue(IS_VISITED) == 0.0) {
                            layers[il + 1].push_back(Node < 3 > ::Pointer(*(jjj.base())));
                            jjj->GetValue(IS_VISITED) = double(il + 2.0);
                        }
                    }
                }
            }


            array_1d<double, 3 > aux, aux_proj;

            //TESTING!!!
            //fill the pressure projection on the first layer inside the fluid
            //by extrapolating from the pressure projection on the layer -1 (the first layer completely inside the domain)
            for (PointIterator iii = (layers[0]).begin(); iii != (layers[0]).end(); iii++)
            {
                    noalias(aux_proj) = ZeroVector(3);
                    double avg_number = 0.0;

                    WeakPointerVector< Node < 3 > >& neighb_nodes = iii->GetValue(NEIGHBOUR_NODES);
                    for (WeakPointerVector< Node < 3 > >::iterator i = neighb_nodes.begin(); i != neighb_nodes.end(); i++)
                    {
                        if (i->GetValue(IS_VISITED) == 0) //the node will be considered for extrapolation only if completely inside
                        {
                            const array_1d<double, 3 > & inside_press_grad = i->FastGetSolutionStepValue(PRESS_PROJ);
                            noalias(aux_proj) += inside_press_grad;
                            avg_number += 1.0;
                        }
                    }

                    if (avg_number != 0.0) //this case means that it has some neighbours that are completely internal
                    {
                        aux_proj /= avg_number;
                        noalias(iii->FastGetSolutionStepValue(PRESS_PROJ)) = aux_proj;
                    }
                    else //case in which there is not a layer of nodes completely internal
                    {
                        noalias(iii->FastGetSolutionStepValue(PRESS_PROJ)) = mRho*mBodyForce;
                    }
            }

            //perform extrapolation layer by layer by making an average
            //of the neighbours of lower order
            for (unsigned int il = 1; il < extrapolation_layers; il++)
            {
                //                                    std::cout << "layer " << il << std::endl;
                for (PointIterator iii = layers[il].begin(); iii != layers[il].end(); iii++) {
                    //                                            std::cout << iii->Id() << " " << std::endl;

                    const array_1d<double, 3 > & coords_top = iii->Coordinates();

                    //extrapolate the average velocity
                    noalias(aux) = ZeroVector(3);
                    noalias(aux_proj) = ZeroVector(3);
                    double avg_number = 0.0;

                    double pavg = 0.0;

                    WeakPointerVector< Node < 3 > >& neighb_nodes = iii->GetValue(NEIGHBOUR_NODES);
                    for (WeakPointerVector< Node < 3 > >::iterator i = neighb_nodes.begin(); i != neighb_nodes.end(); i++)
                    {
                        if (i->GetValue(IS_VISITED) < (il + 1) && i->GetValue(IS_VISITED) != 0) {
                            const array_1d<double, 3 > & coords_bottom = i->Coordinates();
                            array_1d<double, 3 > direction_vec = coords_top;
                            noalias(direction_vec) -= coords_bottom;
                            const array_1d<double, 3 > & press_grad = i->FastGetSolutionStepValue(PRESS_PROJ);
                            double temp = inner_prod(direction_vec, press_grad);
                            double pestimate = i->FastGetSolutionStepValue(PRESSURE,1) + temp;
                            pavg += pestimate;
                            noalias(aux_proj) += press_grad;

                            noalias(aux) += i->FastGetSolutionStepValue(VELOCITY);
                            avg_number += 1.0;
                        }
                    }



                    if (avg_number != 0.0) {
                        aux /= avg_number;
                        pavg /= avg_number;
                        aux_proj /= avg_number;
                    } else {
                        KRATOS_ERROR(std::runtime_error, "error in extrapolation:: no neighbours find on a extrapolation layer -- impossible", "");
                        //                                                    KRATOS_ERROR(std:logic_error,"error in extrapolation:: no neighbours find on a extrapolation layer -- impossible","");
                    }

                    noalias(iii->FastGetSolutionStepValue(VELOCITY)) = aux;
                    noalias(iii->FastGetSolutionStepValue(VELOCITY, 1)) = aux;

                    iii->FastGetSolutionStepValue(PRESSURE, 1) = pavg;

                    noalias(iii->FastGetSolutionStepValue(PRESS_PROJ)) = aux_proj;
                    noalias(iii->FastGetSolutionStepValue(PRESS_PROJ, 1)) = aux_proj;


                }
            }

            //on the first layer outside the pressure is set to a value such that on the free surface the pressure is approx 0
            for (PointIterator iii = layers[1].begin(); iii != layers[1].end(); iii++)
            {
                    //get the node
                    unsigned int i_node = iii->FastGetSolutionStepValue(AUX_INDEX);

                    array_1d<double, TDim> grad_d;
                    for (unsigned int comp = 0; comp < TDim; comp++)
                        grad_d[comp] = 0.0;

                    double dist_i = mdistances[i_node];

                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        //get global index of neighbouring node j
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];

                            const double& dist_j  = mdistances[j_neighbour];

                            //projection of pressure gradients
                            CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

                            edge_ij.Add_grad_p(grad_d, dist_i, dist_j);
                    }

                    const double& m_inv = mr_matrix_container.GetInvertedMass()[i_node];
                    for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                        grad_d[l_comp] *= m_inv;

                    double norm_grad = norm_2(grad_d);

                        if(norm_grad < 100.0)
                        {
                            grad_d /= norm_grad; //this is the direction of the gradient of the distances

                            grad_d *= dist_i; //this is the vector with the distance of node_i from the closest point on the free surface

                            const array_1d<double, TDim> press_grad = iii->FastGetSolutionStepValue(PRESS_PROJ);
                            double pestimate = inner_prod(press_grad,grad_d);

                            iii->FastGetSolutionStepValue(PRESSURE) = pestimate;
                        }
                        else
                        {
                            std::cout << "attention gradient of distance much greater than 1 on node:" << i_node  <<std::endl;
                            double avg_number = 0.0;

                            double pavg = 0.0;

                            WeakPointerVector< Node < 3 > >& neighb_nodes = iii->GetValue(NEIGHBOUR_NODES);
                            for (WeakPointerVector< Node < 3 > >::iterator i = neighb_nodes.begin(); i != neighb_nodes.end(); i++)
                            {
                                if (i->GetValue(IS_VISITED) == 1) {
                                    pavg += i->FastGetSolutionStepValue(PRESSURE);
                                    avg_number += 1.0;
                                }
                            }

                            if(avg_number == 0)
                                KRATOS_ERROR(std::logic_error,"can not happen that the extrapolation node has no neighbours","");

                            iii->FastGetSolutionStepValue(PRESSURE) = pavg/avg_number;

                        }

            }


            //set the pressure to zero on the outer layers (>2)
            for (unsigned int il = 2; il < extrapolation_layers; il++)
            {
                for (PointIterator iii = layers[il].begin(); iii != layers[il].end(); iii++)

                {
                    iii->FastGetSolutionStepValue(PRESSURE) = 0.0;
                }
            }





            //mark nodes on which we will have to solve for convection
            //mark all of internal nodes
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.NodesBegin();
            for (unsigned int i_node = 0; i_node < mr_model_part.Nodes().size(); i_node++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin+i_node;
                if(it->FastGetSolutionStepValue(DISTANCE) <= 0.0)
                    it->GetValue(IS_VISITED) = 1.0;
                else
                    it->GetValue(IS_VISITED) = 0.0;
            }

            //now mark all of the nodes up to the extrapolation layers - 1
            for (unsigned int il = 0; il < extrapolation_layers-1; il++)
                for (PointIterator iii = layers[il].begin(); iii != layers[il].end(); iii++)
                    iii->GetValue(IS_VISITED) = 1.0;

            mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes());
            ApplyVelocityBC(mvel_n1);
            mr_matrix_container.WriteVectorToDatabase(VELOCITY, mvel_n1,  mr_model_part.Nodes());



            KRATOS_CATCH("")
        }

        void ChangeSignToDistance() {
            KRATOS_TRY

            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                    inode != mr_model_part.NodesEnd();
                    inode++) {
                double dist = inode->FastGetSolutionStepValue(DISTANCE);
                inode->FastGetSolutionStepValue(DISTANCE) = -dist;
            }

            KRATOS_CATCH("")
        }

        void MarkNodesByDistance(double min, double max) {
            KRATOS_TRY

            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                    inode != mr_model_part.NodesEnd();
                    inode++) {
                double dist = inode->FastGetSolutionStepValue(DISTANCE);
                if (dist > min && dist < max)
                    inode->GetValue(IS_VISITED) = 1;
                else
                    inode->GetValue(IS_VISITED) = 0;
            }

            KRATOS_CATCH("")
        }

        void SaveScalarVariableToOldStep(Variable<double>& rVar) {
            KRATOS_TRY

            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                    inode != mr_model_part.NodesEnd();
                    inode++) {
                inode->FastGetSolutionStepValue(rVar, 1) = inode->FastGetSolutionStepValue(rVar);
            }

            KRATOS_CATCH("")
        }

        void MarkExternalAndMixedNodes() {
            KRATOS_TRY

            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                    inode != mr_model_part.NodesEnd();
                    inode++) {
                inode->GetValue(IS_VISITED) = 0;
            }

            //detect the nodes inside the fluid surface
            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                    inode != mr_model_part.NodesEnd();
                    inode++) {
                if (inode->FastGetSolutionStepValue(DISTANCE) > 0.0) //candidates are only the ones inside the fluid domain
                {
                    inode->GetValue(IS_VISITED) = 1;
                    WeakPointerVector< Node < 3 > >& neighb_nodes = inode->GetValue(NEIGHBOUR_NODES);
                    for (WeakPointerVector< Node < 3 > >::iterator i = neighb_nodes.begin(); i != neighb_nodes.end(); i++) {
                        i->GetValue(IS_VISITED) = 1;
                    }
                }
            }
            KRATOS_CATCH("")
        }

        void MarkInternalAndMixedNodes() {
            KRATOS_TRY

            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                    inode != mr_model_part.NodesEnd();
                    inode++) {
                inode->GetValue(IS_VISITED) = 0;
            }

            //detect the nodes inside the fluid surface
            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                    inode != mr_model_part.NodesEnd();
                    inode++) {
                if (inode->FastGetSolutionStepValue(DISTANCE) <= 0.0) //candidates are only the ones inside the fluid domain
                {
                    inode->GetValue(IS_VISITED) = 1;
                    WeakPointerVector< Node < 3 > >& neighb_nodes = inode->GetValue(NEIGHBOUR_NODES);
                    for (WeakPointerVector< Node < 3 > >::iterator i = neighb_nodes.begin(); i != neighb_nodes.end(); i++) {
                        i->GetValue(IS_VISITED) = 1;
                    }
                }
            }
            KRATOS_CATCH("")
        }

        void MarkInternalNodes() {
            KRATOS_TRY

            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                    inode != mr_model_part.NodesEnd();
                    inode++) {
                inode->GetValue(IS_VISITED) = 0;
            }

            //detect the nodes inside the fluid surface
            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                    inode != mr_model_part.NodesEnd();
                    inode++) {
                if (inode->FastGetSolutionStepValue(DISTANCE) <= 0.0) //candidates are only the ones inside the fluid domain
                {
                    inode->GetValue(IS_VISITED) = 1;
                }
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
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                if (is_slip[i_node] == true)
                    mSlipBoundaryList.push_back(i_node);
            }



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
            mPressureOutlet.clear();
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
	   mD.clear();
            mdiv_error.clear();

            KRATOS_CATCH("")
        }

        void ConvectDistance()
        {
            KRATOS_TRY

            //variables for node based data handling
            ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
            int n_nodes = rNodes.size();

            //storage of nodal values in local variables
            ValuesVectorType rhs, WorkConvection;
            rhs.resize(n_nodes);
            WorkConvection.resize(n_nodes);

            ValuesVectorType active_nodes;
            active_nodes.resize(n_nodes);

            mr_matrix_container.FillScalarFromDatabase(POROSITY, mEps, mr_model_part.Nodes());


            //read variables from Kratos
            mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes());
            mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mvel_n, mr_model_part.Nodes());

            mr_matrix_container.FillScalarFromDatabase(DISTANCE, mphi_n1, mr_model_part.Nodes());
            mr_matrix_container.FillOldScalarFromDatabase(DISTANCE, mphi_n, mr_model_part.Nodes());
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

            //calculating the convective projection
            array_1d<double, TDim> a_i;
            array_1d<double, TDim> a_j;
            #pragma omp parallel for private(a_i,a_j)
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                double& pi_i = mPiConvection[i_node];
                const double& phi_i = mphi_n1[i_node];
                //set to zero the projection
                pi_i = 0.0;
                if (active_nodes[i_node] != 0.0)
                {
                    a_i = mvel_n1[i_node];
                    a_i /= mEps[i_node];

                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];

                        if (active_nodes[j_neighbour] != 0.0)
                        {
                            noalias(a_j) = mvel_n1[j_neighbour];
                            a_j /= mEps[j_neighbour];

                            const double& phi_j = mphi_n1[j_neighbour];
                            CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
                            edge_ij.Add_ConvectiveContribution(pi_i, a_i, phi_i, a_j, phi_j);
                        }
                    }
                    //apply inverted mass matrix
                    const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];
                    pi_i *= m_inv;
                }
            }


            //create and fill a vector of nodes for which we want to convect the velocity
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                ModelPart::NodesContainerType::iterator it_begin = mr_model_part.NodesBegin();
                active_nodes[i_node] = (it_begin + i_node)->GetValue(IS_VISITED);
            }
//            mr_matrix_container.WriteScalarToDatabase(TEMPERATURE, active_nodes, rNodes);

            //read time step size from Kratos
            ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
            double delta_t = CurrentProcessInfo[DELTA_TIME];

            mr_matrix_container.AssignVectorToVector(mphi_n, WorkConvection); //mWork = mphi_n

            //first step of Runge Kutta
            // 				mr_matrix_container.AssignVectorToVector(mphi_n,mphi_n1); //mphi_n1 = mphi_n
            mr_matrix_container.SetToZero(rhs);
            CalculateRHS_convection(mphi_n1, mvel_n1, rhs, active_nodes);
            mr_matrix_container.Add_Minv_value(WorkConvection, WorkConvection, delta_t / 6.0, mr_matrix_container.GetInvertedMass(), rhs);
            mr_matrix_container.Add_Minv_value(mphi_n1, mphi_n, 0.5 * delta_t, mr_matrix_container.GetInvertedMass(), rhs);

            //second step
            mr_matrix_container.SetToZero(rhs);
            CalculateRHS_convection(mphi_n1, mvel_n1, rhs, active_nodes);
            mr_matrix_container.Add_Minv_value(WorkConvection, WorkConvection, delta_t / 3.0, mr_matrix_container.GetInvertedMass(), rhs);
            mr_matrix_container.Add_Minv_value(mphi_n1, mphi_n, 0.5 * delta_t, mr_matrix_container.GetInvertedMass(), rhs);

            //third step
            mr_matrix_container.SetToZero(rhs);
            CalculateRHS_convection(mphi_n1, mvel_n1, rhs, active_nodes);
            mr_matrix_container.Add_Minv_value(WorkConvection, WorkConvection, delta_t / 3.0, mr_matrix_container.GetInvertedMass(), rhs);
            mr_matrix_container.Add_Minv_value(mphi_n1, mphi_n, delta_t, mr_matrix_container.GetInvertedMass(), rhs);

            //fourth step
            mr_matrix_container.SetToZero(rhs);
            CalculateRHS_convection(mphi_n1, mvel_n1, rhs, active_nodes);
            mr_matrix_container.Add_Minv_value(WorkConvection, WorkConvection, delta_t / 6.0, mr_matrix_container.GetInvertedMass(), rhs);

            //compute right-hand side
            mr_matrix_container.AssignVectorToVector(WorkConvection, mphi_n1);


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



            mr_matrix_container.WriteScalarToDatabase(DISTANCE, mphi_n1, mr_model_part.Nodes());



            KRATOS_CATCH("")
        }

void ReduceTimeStep(ModelPart& rModelPart, double NewTime)
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
    rModelPart.OverwriteSolutionStepData(1,0);
    rModelPart.GetProcessInfo().SetCurrentTime(NewTime);

    KRATOS_CATCH("error in reducing the time step")

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

                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        //get global index of neighbouring node j
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];

                            const double& dist_j  = mdistances[j_neighbour];

                            //projection of pressure gradients
                            CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

                            edge_ij.Add_grad_p(grad_d, dist_i, dist_j);
                    }

                    const double& m_inv = mr_matrix_container.GetInvertedMass()[i_node];
                    for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                        grad_d[l_comp] *= m_inv;

                    double norm_grad = norm_2(grad_d);

                    if(norm_grad > 1.5) //large gradient found
                        n_large_distance_gradient += 1;
                }
            }

            if(n_large_distance_gradient != 0)
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


void ActivateWallResistance(double Ywall)
        {
            mWallLawIsActive = true;
            mY_wall = Ywall;
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
        //coefficients
        ValuesVectorType mdistances;
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
        ValuesVectorType mPressureOutlet;

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
        ValuesVectorType mphi_n;
        ValuesVectorType mphi_n1;
        ValuesVectorType mPiConvection;
        ValuesVectorType mBeta;

        //variables for edge BCs
        IndicesVectorType medge_nodes;
        CalcVectorType medge_nodes_direction;
        IndicesVectorType mcorner_nodes;

        ValuesVectorType mEps;
        ValuesVectorType mD;

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





        //*********************************************************************
        //function to calculate right-hand side of fractional momentum equation

        void CalculateRHS_convection(
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
                const double& phi_i = mphi[i_node];
                noalias(a_i) = convective_velocity[i_node];
                a_i /= mEps[i_node];


                const double& pi_i = mPiConvection[i_node];
//                double beta = mBeta[i_node];
                rhs_i = 0.0;

                if (active_nodes[i_node] != 0.0)
                {

                    //loop to all the edges surrounding node I
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];

                        if (active_nodes[j_neighbour] != 0.0)
                        {

                            //double& rhs_j = rhs[j_neighbour];
                            const double& phi_j = mphi[j_neighbour];
                            noalias(a_j) = convective_velocity[j_neighbour];
                            a_j /= mEps[j_neighbour];

                            const double& pi_j = mPiConvection[j_neighbour];

                            CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

                            //convection operator
//                            edge_ij.Sub_ConvectiveContribution(rhs_i, a_i, phi_i, a_j, phi_j); //esto funciona
                            edge_ij.Sub_D_v(rhs_i, a_i*phi_i, a_i*phi_j);


                            //calculate stabilization part
                            edge_ij.CalculateConvectionStabilization_LOW(stab_low, a_i, phi_i, a_j, phi_j);

                            double edge_tau = mTauConvection[i_node];

                            edge_ij.CalculateConvectionStabilization_HIGH(stab_high, a_i, pi_i, a_j, pi_j);

                            edge_ij.Sub_StabContribution(rhs_i, edge_tau, 1.0, stab_low, stab_high);
                        }
                    }
                }

                // KRATOS_WATCH(rhs_i);

            }


            KRATOS_CATCH("")
        }


        //**************************************
        void CornerDectectionHelper(Geometry< Node<3> >& face_geometry,
                                    const array_1d<double,3>& face_normal,
                                    const double An,
                                    const WeakPointerVector<Condition>& neighb,
                                    const unsigned int i1,
                                    const unsigned int i2,
                                    const unsigned int neighb_index,
                                    std::vector<unsigned int>& edge_nodes,
                                    CalcVectorType& cornern_list
                                    )
        {
            double acceptable_angle = 45.0/180.0*3.1; //angles of less than 45 deg will be accepted
            double acceptable_cos = cos(acceptable_angle);

            if(face_geometry[i1].Id() < face_geometry[i2].Id()) //we do this to add the face ones
                        {
                            const array_1d<double, 3 > & neighb_normal = neighb[neighb_index].GetValue(NORMAL);
                            double neighb_An = norm_2(neighb_normal);

                            double cos_normal = 1.0/(An*neighb_An) * inner_prod(face_normal,neighb_normal);

                            //if the angle is too big between the two normals then the edge in the middle is a corner
                            if(cos_normal < acceptable_cos)
                            {
                                array_1d<double, 3 > edge = face_geometry[i2].Coordinates() -  face_geometry[i1].Coordinates();
                                double temp = norm_2(edge);
                                edge/=temp;

                                int index1 = face_geometry[i1].FastGetSolutionStepValue(AUX_INDEX);
                                int index2 = face_geometry[i2].FastGetSolutionStepValue(AUX_INDEX);

                                edge_nodes[index1] += 1;
                                edge_nodes[index2] += 1;

                                cornern_list[index1] += edge;
                                cornern_list[index2] += edge;


                            }



                        }


        }

        //function to calculate the area normals
        void DetectEdges3D(ModelPart::ConditionsContainerType& rConditions) {
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
                    if(neighb[0].Id() != current_id) //check if the neighbour exists
                        CornerDectectionHelper(face_geometry,face_normal,An,neighb,1,2, 0,temp_edge_nodes,temp_cornern_list );

                    //check for neighbour one
                    if(neighb[0].Id() != current_id) //check if the neighbour exists
                        CornerDectectionHelper(face_geometry,face_normal,An,neighb,2,0, 1,temp_edge_nodes,temp_cornern_list );

                    //check for neighbour two
                    if(neighb[0].Id() != current_id) //check if the neighbour exists
                        CornerDectectionHelper(face_geometry,face_normal,An,neighb,0,1, 2,temp_edge_nodes,temp_cornern_list );

                }
            }

//            ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
//            mr_matrix_container.WriteVectorToDatabase(ACCELERATION, temp_cornern_list, rNodes);


            //fill the list of edge_nodes
            medge_nodes.resize(0);
            medge_nodes_direction.resize(0);
            mcorner_nodes.resize(0);
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                if (temp_edge_nodes[i_node] == 2) //node is a edge_node
                {
                    medge_nodes.push_back(i_node);
                    array_1d<double,TDim>& node_edge = temp_cornern_list[i_node];
                    node_edge /= norm_2(node_edge);
                    medge_nodes_direction.push_back(node_edge);
                }
                else if (temp_edge_nodes[i_node] > 2)
                    mcorner_nodes.push_back(i_node);
            }



            for (unsigned int i = 0; i < mcorner_nodes.size(); i++)
            {
                KRATOS_WATCH(mcorner_nodes[i]);

            }


            KRATOS_CATCH("")
        }

        double ComputePorosityCoefficient(const double& viscosity, const double& vel_norm, const double& eps, const double& d)
        {
//             const double d = 0.01; //to be changed

            double k_inv = 150.0 * (1.0 - eps)*(1.0 - eps) / (eps * eps * eps * d * d);
            double linear = viscosity * k_inv;
            double non_linear = (1.75 * vel_norm / eps) * sqrt(k_inv / (150.0 * eps));
            return linear + non_linear;
        }

        void LaplacianSmooth(ValuesVectorType& target, const ValuesVectorType& origin)
        {
            ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
            int n_nodes = rNodes.size();
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                double dist = mdistances[i_node];
                double& correction = target[i_node];
                correction = 0.0;

                if (dist <= 0.0) //node is inside domain ---- if outside do nothing
                {
                    const double& origin_i = origin[i_node];
                    double& target_i = target[i_node];
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        const double& origin_j = origin[j_neighbour];

                        CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

                        double l_ikjk;
                        edge_ij.CalculateScalarLaplacian(l_ikjk);

                        correction += l_ikjk * (origin_j - origin_i);
                    }

                    correction += origin_i;
                }
            }
        }


void ComputeWallResistance(
                const CalcVectorType& vel,
                CalcVectorType& rhs
                )
        {
            //parameters:
            const double k = 0.41;
            const double B = 5.1;
            const double density = mRho;
            const double mu = mViscosity;
            const double toll = 1e-6;
            const double ym = mY_wall; //0.0825877; //0.0093823
            const double y_plus_incercept = 10.9931899;
            const unsigned int itmax = 100;

            if (mu == 0)
                KRATOS_ERROR(std::logic_error, "it is not possible to use the wall law with 0 viscosity", "");

            //slip condition
            int slip_size = mSlipBoundaryList.size();
#pragma omp parallel for firstprivate(slip_size,B,density,mu,toll,ym,y_plus_incercept,itmax)
            for (int i_slip = 0; i_slip < slip_size; i_slip++)
            {
                unsigned int i_node = mSlipBoundaryList[i_slip];
                double dist = mdistances[i_node];
                if (dist <= 0.0)
                {
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
        }




    };
} //namespace Kratos

#endif //KRATOS_EDGEBASED_LEVELSET_FLUID_SOLVER_H_INCLUDED defined


