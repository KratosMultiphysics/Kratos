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
#define KRATOS_EDGEBASED_LEVELSET_SUBSTEP_FLUID_SOLVER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

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
#include "utilities/reduction_utilities.h"

namespace Kratos
{
    template <unsigned int TDim, class MatrixContainer, class TSparseSpace, class TLinearSolver>
    class EdgeBasedLevelSetSubstep
    {
    public:
        // name for the self defined structure
        typedef EdgesStructureTypeC2C<TDim> CSR_Tuple;
        typedef vector<CSR_Tuple> EdgesVectorType;
        // name for row start and column index vectors
        typedef vector<unsigned int> IndicesVectorType;
        // defining matrix type for test calculations
        typedef vector<array_1d<double, TDim>> CalcVectorType;
        // defining type for local storage of nodal values
        typedef vector<double> ValuesVectorType;
        // defining types for matrix operations
        typedef typename TSparseSpace::MatrixType TSystemMatrixType;
        typedef typename TSparseSpace::VectorType TSystemVectorType;
        typedef std::size_t SizeType;
        // constructor and destructor
        EdgeBasedLevelSetSubstep(MatrixContainer &mr_matrix_container,
                                 ModelPart &mr_model_part,
                                 const double viscosity,
                                 const double density,
                                 bool use_mass_correction,
                                 double stabdt_pressure_factor,
                                 double stabdt_convection_factor,
                                 double tau2_factor,
                                 bool assume_constant_dp)
            : mr_matrix_container(mr_matrix_container),
              mr_model_part(mr_model_part),
              mstabdt_pressure_factor(stabdt_pressure_factor),
              mstabdt_convection_factor(stabdt_convection_factor),
              mtau2_factor(tau2_factor),
              massume_constant_dp(assume_constant_dp)
        {
            for (ModelPart::NodesContainerType::iterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); it++)
                it->FastGetSolutionStepValue(VISCOSITY) = viscosity;

            mRho = density;
            mdelta_t_avg = 1000.0;
            max_dt = 1.0;
            muse_mass_correction = use_mass_correction;
            mshock_coeff = 0.7;
            mWallLawIsActive = false;
            mnumsubsteps = 5;
            mmax_dt = 0.0;
            mcorner_coefficient = 30.0;
            medge_coefficient = 2.0;
            std::cout << "Edge based level set substep solver is created" << std::endl;
        };
        ~EdgeBasedLevelSetSubstep(){};

        //***********************************
        // function to initialize fluid solver
        void Initialize()
        {
            KRATOS_TRY
            // get number of nodes
            unsigned int n_nodes = mr_model_part.Nodes().size();
            unsigned int n_edges = mr_matrix_container.GetNumberEdges();
            // size data vectors
            mBodyForce.resize(n_nodes);
            mr_matrix_container.SetToZero(mBodyForce);
            mViscosity.resize(n_nodes);
            mr_matrix_container.SetToZero(mViscosity);
            mWork.resize(n_nodes);
            mr_matrix_container.SetToZero(mWork);
            mvel_n.resize(n_nodes);
            mr_matrix_container.SetToZero(mvel_n);
            mvel_n1.resize(n_nodes);
            mr_matrix_container.SetToZero(mvel_n1);
            mPn.resize(n_nodes);
            mr_matrix_container.SetToZero(mPn);
            mPn1.resize(n_nodes);
            mr_matrix_container.SetToZero(mPn1);
            mHmin.resize(n_nodes);
            mr_matrix_container.SetToZero(mHmin);
            mHavg.resize(n_nodes);
            mr_matrix_container.SetToZero(mHavg);
            mNodalFlag.resize(n_nodes);
            mr_matrix_container.SetToZero(mNodalFlag);
            mdistances.resize(n_nodes);
            mr_matrix_container.SetToZero(mdistances);
            mTauPressure.resize(n_nodes);
            mr_matrix_container.SetToZero(mTauPressure);
            mTauConvection.resize(n_nodes);
            mr_matrix_container.SetToZero(mTauConvection);
            mTau2.resize(n_nodes);
            mr_matrix_container.SetToZero(mTau2);
            mPi.resize(n_nodes);
            mr_matrix_container.SetToZero(mPi);
            mXi.resize(n_nodes);
            mr_matrix_container.SetToZero(mXi);
            mx.resize(n_nodes);
            mr_matrix_container.SetToZero(mx);
            mEdgeDimensions.resize(n_edges);
            mr_matrix_container.SetToZero(mEdgeDimensions);
            // convection variables
            mBeta.resize(n_nodes);
            mr_matrix_container.SetToZero(mBeta);
            mPiConvection.resize(n_nodes);
            mr_matrix_container.SetToZero(mPiConvection);
            mphi_n.resize(n_nodes);
            mr_matrix_container.SetToZero(mphi_n);
            mphi_n1.resize(n_nodes);
            mr_matrix_container.SetToZero(mphi_n1);
            mEps.resize(n_nodes);
            mr_matrix_container.SetToZero(mEps);
            mA.resize(n_nodes);
            mr_matrix_container.SetToZero(mA);
            mB.resize(n_nodes);
            mr_matrix_container.SetToZero(mB);
            mdiv_error.resize(n_nodes);
            mr_matrix_container.SetToZero(mdiv_error);
            mWallReductionFactor.resize(n_nodes);
            mr_matrix_container.SetToZero(mWallReductionFactor);
            mdiag_stiffness.resize(n_nodes);
            mr_matrix_container.SetToZero(mdiag_stiffness);
            mis_slip.resize(n_nodes);
            mis_visited.resize(n_nodes);
            macc.resize(n_nodes);
            mr_matrix_container.SetToZero(macc);

            // read velocity and pressure data from Kratos
            mr_matrix_container.FillVectorFromDatabase(BODY_FORCE, mBodyForce, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(VISCOSITY, mViscosity, mr_model_part.Nodes());
            mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, mr_model_part.Nodes());
            mr_matrix_container.FillOldScalarFromDatabase(PRESSURE, mPn, mr_model_part.Nodes());
            mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mvel_n, mr_model_part.Nodes());
            mr_matrix_container.FillCoordinatesFromDatabase(mx, mr_model_part.Nodes());
            // set flag for first time step
            mFirstStep = true;
            // loop to categorize boundary nodes
            std::vector<unsigned int> tempFixedVelocities;
            std::vector<array_1d<double, TDim>> tempFixedVelocitiesValues;
            std::vector<unsigned int> tempPressureOutletList;
            std::vector<unsigned int> tempDistanceList;
            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                 inode != mr_model_part.NodesEnd();
                 inode++)
            {
                int index = inode->FastGetSolutionStepValue(AUX_INDEX);

                if (inode->Is(INLET))
                {
                    tempFixedVelocities.push_back(index);
                    tempFixedVelocitiesValues.push_back(mvel_n1[index]);
                }

                if (inode->IsFixed(DISTANCE))
                {
                    tempDistanceList.push_back(index);
                }

                if (inode->Is(OUTLET) || inode->IsFixed(PRESSURE))
                {
                    tempPressureOutletList.push_back(index);
                }
            }
            mFixedVelocities.resize(tempFixedVelocities.size(), false);
            mFixedVelocitiesValues.resize(tempFixedVelocitiesValues.size(), false);
            mPressureOutletList.resize(tempPressureOutletList.size(), false);
            mDistanceBoundaryList.resize(tempDistanceList.size(), false);
            mDistanceValuesList.resize(tempDistanceList.size(), false);

            IndexPartition<unsigned int>(tempFixedVelocities.size()).for_each([&](unsigned int i){
                mFixedVelocities[i] = tempFixedVelocities[i];
                mFixedVelocitiesValues[i] = tempFixedVelocitiesValues[i];
            });

            IndexPartition<unsigned int>(tempPressureOutletList.size()).for_each([&](unsigned int i){
                mPressureOutletList[i] = tempPressureOutletList[i];
            });

            for (int i = 0; i < static_cast<int>(tempDistanceList.size()); i++)
            {
                mDistanceBoundaryList[i] = tempDistanceList[i];
            }
            // compute slip normals and fill SlipList
            CalculateNormals(mr_model_part.Conditions());
            mr_matrix_container.WriteVectorToDatabase(NORMAL, mSlipNormal, mr_model_part.Nodes());
            if constexpr (TDim == 3)
                DetectEdges3D(mr_model_part.Conditions());

            // determine number of edges and entries
            unsigned int n_nonzero_entries = 2 * n_edges + n_nodes;
            // allocate memory for variables
            mL.resize(n_nodes, n_nodes, n_nonzero_entries);
            int number_of_threads = ParallelUtilities::GetNumThreads();
            std::vector<int> row_partition(number_of_threads);
            OpenMPUtils::DivideInPartitions(n_nodes, number_of_threads, row_partition);
            for (int k = 0; k < number_of_threads; k++)
            {
                #pragma omp parallel
                if (OpenMPUtils::ThisThread() == k)
                {
                    for (int i_node = static_cast<int>(row_partition[k]); i_node < static_cast<int>(row_partition[k + 1]); i_node++)
                    {
                        // loop over all nodes
                        // flag for considering diagonal matrix elements
                        bool flag = 0;
                        // loop over all neighbours
                        for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                        {
                            // get global index of neighbouring node j
                            unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                            // define matrix structure row by row (the order does matter!)
                            if ((static_cast<int>(j_neighbour) > i_node) && (flag == 0))
                            {
                                // add diagonal/nodal contribution
                                mL.push_back(i_node, i_node, 0.0);
                                flag = 1;
                            }
                            // add non-diagonal/edge contribution
                            mL.push_back(i_node, j_neighbour, 0.0);
                        }
                        // if diagonal element is the last non-zero element of the row
                        if (flag == 0)
                            mL.push_back(i_node, i_node, 0.0);
                    }
                }
            }
            // compute minimum length of the surrounding edges
            CalculateEdgeLengths(mr_model_part.Nodes());

            array_1d<double, 3> temp_body_force;

            // set the pressure projection to the body force value
            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                 inode != mr_model_part.NodesEnd();
                 inode++)
            {
                temp_body_force = mRho * inode->FastGetSolutionStepValue(BODY_FORCE);
                inode->FastGetSolutionStepValue(PRESS_PROJ) = temp_body_force;
            }

            mr_matrix_container.FillScalarFromDatabase(POROSITY, mEps, mr_model_part.Nodes());
            // verify that neither h_min nor havg are 0
            for (unsigned int i_node = 0; i_node < mHmin.size(); i_node++)
            {
                KRATOS_ERROR_IF(mHmin[i_node] < 1.0e-15) << "hmin too small on node " << i_node + 1 << std::endl;
                KRATOS_ERROR_IF(mHmin[i_node] < 1.0e-15) << "hmin too small on node " << i_node + 1 << std::endl;
                KRATOS_ERROR_IF(mHmin[i_node] > 1.0e15) << "hmin too big on node " << i_node + 1 << std::endl;
                KRATOS_ERROR_IF(mHmin[i_node] > 1.0e15) << "havg too big on node " << i_node + 1 << std::endl;
            }
            for (ModelPart::ElementsContainerType::iterator it = mr_model_part.ElementsBegin(); it != mr_model_part.ElementsEnd(); it++)
            {
                KRATOS_ERROR_IF(it->Id() < 1) << "Element found with Id 0 or negative" << std::endl;

                double elem_vol = 0.0;
                if constexpr (TDim == 2)
                    elem_vol = it->GetGeometry().Area();
                else
                    elem_vol = it->GetGeometry().Volume();

                KRATOS_ERROR_IF(elem_vol <= 0) << "error on element -> " << it->Id() << ". Area can not be lesser than 0" << std::endl;
            }
            KRATOS_CATCH("")
        }

        void SetShockCapturingCoefficient(double coeff)
        {
            mshock_coeff = coeff;
        }

        void GatherValues()
        {
            KRATOS_TRY

            mr_matrix_container.FillVectorFromDatabase(BODY_FORCE, mBodyForce, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(VISCOSITY, mViscosity, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(POROSITY, mEps, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(DISTANCE, mdistances, mr_model_part.Nodes());
            mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes());
            mr_matrix_container.FillVectorFromDatabase(PRESS_PROJ, mXi, mr_model_part.Nodes());

            mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mvel_n, mr_model_part.Nodes());
            mr_matrix_container.FillOldScalarFromDatabase(PRESSURE, mPn, mr_model_part.Nodes());

            KRATOS_CATCH("")
        }

        //***************************************
        // function to set adequate time step size
        double ComputeTimeStep(const double CFLNumber, const double MaxDt)
        {
            KRATOS_TRY
            // save the maximum time step
            max_dt = MaxDt;

            // loop over all nodes
            int n_nodes = static_cast<int>(mvel_n1.size());

            unsigned int n_proc = ParallelUtilities::GetNumThreads();
            Vector dt_avg_vec(n_proc, 1e10);
            Vector dt_vec(n_proc, 1e10);
            Vector dt_avg_novisc_vec(n_proc, 1e10);

            #pragma omp parallel for firstprivate(n_nodes)
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                unsigned int my_id = OpenMPUtils::ThisThread();
                double &delta_t = dt_vec[my_id];
                double &mdelta_t_avg = dt_avg_vec[my_id];
                double &delta_t_avg_novisc = dt_avg_novisc_vec[my_id];

                const array_1d<double, TDim> &v_i = mvel_n1[i_node];
                const double havg_i = mHavg[i_node];
                const double hmin_i = mHmin[i_node];
                const double eps_i = mEps[i_node];
                double nu = mViscosity[i_node];
                double vel_norm = norm_2(v_i);

                vel_norm /= eps_i;

                // use CFL condition to compute time step size
                double delta_t_i = 1.0 / (vel_norm / hmin_i + nu / (hmin_i * hmin_i));
                const double delta_t_i_avg = 1.0 / (vel_norm / havg_i + nu / (havg_i * havg_i));
                double delta_t_i_avg_novisc = 1.0 / (2.0 * vel_norm / havg_i);

                // considering the most restrictive case of neighbor's velocities with similar direction but opposite sense.
                // loop over all neighbours
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    // get global index of neighbouring node j
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    const array_1d<double, TDim> &v_j = mvel_n1[j_neighbour];
                    double v_diff_norm = 0.0;
                    for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                    {
                        double temp = v_i[l_comp] - v_j[l_comp];
                        v_diff_norm += temp * temp;
                    }
                    v_diff_norm = sqrt(v_diff_norm);
                    v_diff_norm /= eps_i;
                    double delta_t_j = 1.0 / (v_diff_norm / havg_i + 4.0 * nu / (havg_i * havg_i));
                    double delta_t_j_avg_novisc = 1.0 / (2.0 * v_diff_norm / havg_i);
                    if (delta_t_j < delta_t_i)
                        delta_t_i = delta_t_j;
                    if (delta_t_j_avg_novisc < delta_t_i_avg_novisc)
                        delta_t_i_avg_novisc = delta_t_j_avg_novisc;
                }
                // choose the overall minimum of delta_t_i
                if (delta_t_i < delta_t)
                    delta_t = delta_t_i;
                if (delta_t_i_avg < mdelta_t_avg)
                    mdelta_t_avg = delta_t_i_avg;
                if (delta_t_i_avg_novisc < delta_t_avg_novisc)
                    delta_t_avg_novisc = delta_t_i_avg_novisc;
            }

            // finalizing parallel computations
            double delta_t = dt_vec[0];
            mdelta_t_avg = dt_avg_vec[0];
            double delta_t_avg_novisc = dt_avg_novisc_vec[0];
            for (unsigned int i = 1; i < dt_vec.size(); i++)
            {
                if (delta_t > dt_vec[i])
                    delta_t = dt_vec[i];
                if (mdelta_t_avg > dt_vec[i])
                    mdelta_t_avg = dt_avg_vec[i];
                if (delta_t_avg_novisc > dt_vec[i])
                    delta_t_avg_novisc = dt_avg_novisc_vec[i];
            }

            delta_t_avg_novisc *= CFLNumber;

            mnumsubsteps = ceil(delta_t_avg_novisc / delta_t);

            if (mnumsubsteps <= 1)
            {
                mnumsubsteps = 1;
                delta_t_avg_novisc = delta_t;
            }

            delta_t = delta_t_avg_novisc;

            return delta_t;
            KRATOS_CATCH("")
        }

        void ApplySmagorinsky(double MolecularViscosity, double Cs)
        {
            if (Cs != 0)
            {
                if constexpr (TDim == 3)
                    ApplySmagorinsky3D(MolecularViscosity, Cs);
                else
                    KRATOS_THROW_ERROR(std::logic_error, "smagorinsky not yet implemented in 2D", "");
            }
        }

        void UpdateFixedVelocityValues()
        {
            KRATOS_TRY
            // read velocity and pressure data from Kratos
            int fixed_size = mFixedVelocities.size();

            #pragma omp parallel for firstprivate(fixed_size)
            for (int i_velocity = 0; i_velocity < fixed_size; i_velocity++)
            {
                unsigned int i_node = mFixedVelocities[i_velocity];
                array_1d<double, TDim> &u_i_fix = mFixedVelocitiesValues[i_velocity];
                const array_1d<double, TDim> &u_i = mvel_n1[i_node];
                for (unsigned int comp = 0; comp < TDim; comp++)
                    u_i_fix[comp] = u_i[comp];
            }
            KRATOS_CATCH("");
        }

        //**********************************************************************************
        // function to solve fluid equations - fractional step 1: compute fractional momentum
        void SolveStep1()
        {
            KRATOS_TRY
            // PREREQUISITES
            // variables for node based data handling
            ModelPart::NodesContainerType &rNodes = mr_model_part.Nodes();
            int n_nodes = rNodes.size();
            // storage of nodal values in local variables
            CalcVectorType rhs;
            rhs.resize(n_nodes);

            mr_matrix_container.FillScalarFromDatabase(DISTANCE, mdistances, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(POROSITY, mEps, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(LIN_DARCY_COEF, mA, mr_model_part.Nodes());
            mr_matrix_container.FillScalarFromDatabase(NONLIN_DARCY_COEF, mB, mr_model_part.Nodes());

            // read time step size from Kratos
            ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
            double delta_t = CurrentProcessInfo[DELTA_TIME];

            // compute intrinsic time
            double time_inv_avg = 1.0 / mdelta_t_avg;

            double stabdt_pressure_factor = mstabdt_pressure_factor;
            double stabdt_convection_factor = mstabdt_convection_factor;

            #pragma omp parallel for firstprivate(time_inv_avg, stabdt_pressure_factor, stabdt_convection_factor)
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                double &h_avg_i = mHavg[i_node];
                double &h_min_i = mHmin[i_node];

                array_1d<double, TDim> &a_i = mvel_n1[i_node];
                const double nu_i = mViscosity[i_node];
                const double eps_i = mEps[i_node];
                const double lindarcy_i = mA[i_node];
                const double nonlindarcy_i = mB[i_node];
                double vel_norm = norm_2(a_i);
                double porosity_coefficient = ComputePorosityCoefficient(vel_norm, eps_i, lindarcy_i, nonlindarcy_i);
                vel_norm /= eps_i;
                double tau = 1.0 / (2.0 * vel_norm / h_min_i + stabdt_pressure_factor * time_inv_avg + (4.0 * nu_i) / (h_avg_i * h_avg_i) + porosity_coefficient);
                double tau_conv = 1.0 / (2.0 * vel_norm / h_min_i + stabdt_convection_factor * time_inv_avg);
                mTauPressure[i_node] = tau;
                mTauConvection[i_node] = tau_conv;
            }

            // calculating the convective projection
            IndexPartition<unsigned int>(n_nodes).for_each([&](unsigned int i_node){
                array_1d<double, TDim> &pi_i = mPi[i_node];
                // setting to zero
                for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                    pi_i[l_comp] = 0.0;
                array_1d<double, TDim> a_i = mvel_n1[i_node];
                const array_1d<double, TDim> &U_i = mvel_n1[i_node];
                const double &eps_i = mEps[i_node];
                a_i /= eps_i;

                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    array_1d<double, TDim> a_j = mvel_n1[j_neighbour];
                    const array_1d<double, TDim> &U_j = mvel_n1[j_neighbour];
                    const double &eps_j = mEps[j_neighbour];
                    a_j /= eps_j;
                    CSR_Tuple &edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
                    edge_ij.Add_ConvectiveContribution(pi_i, a_i, U_i, a_j, U_j);
                }
            });

            int inout_size = mInOutBoundaryList.size();

            for (int i = 0; i < inout_size; i++)
            {
                unsigned int i_node = mInOutBoundaryList[i];
                const array_1d<double, TDim> &U_i = mvel_n1[i_node];
                const array_1d<double, TDim> &an_i = mInOutNormal[i_node];
                double projection_length = 0.0;
                for (unsigned int comp = 0; comp < TDim; comp++)
                {
                    projection_length += U_i[comp] * an_i[comp];
                }

                array_1d<double, TDim> &pi_i = mPi[i_node];

                for (unsigned int comp = 0; comp < TDim; comp++)
                    pi_i[comp] += projection_length * U_i[comp];
            }

            IndexPartition<unsigned int>(n_nodes).for_each([&](unsigned int i_node){
                array_1d<double, TDim> &pi_i = mPi[i_node];

                const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];
                for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                    pi_i[l_comp] *= m_inv;
            });

            CalcVectorType auxn = mvel_n;

            double n_substeps = mnumsubsteps + 1;
            double reduced_it = 0;

            std::vector<int> initial_data_vector(n_nodes);

            auto initial_assembled_vector = block_for_each<AccumReduction<int>>(initial_data_vector, [&](int& i_node){
                double temp;
                if (mdistances[i_node] <= 0.0)
                    temp = mr_matrix_container.GetLumpedMass()[i_node] * inner_prod(mvel_n[i_node], mvel_n[i_node]);
                else
                    temp = 0.0;
                return temp;
            });

            double energy_initial = accumulate(initial_assembled_vector.begin(), initial_assembled_vector.end(), 0);

            while (reduced_it++ < 2)
            {
                double delta_t_substep = delta_t / n_substeps;
                for (unsigned int substep = 0; substep < n_substeps; substep++)
                {
                    // std::cout << "substep " << substep+1 << " of " << n_substeps << std::endl;
                    mr_matrix_container.AssignVectorToVector(mvel_n, mWork); // mWork = mvel_n
                    // first step of Runge Kutta
                    mr_matrix_container.AssignVectorToVector(mvel_n, mvel_n1); // mvel_n1 = mvel_n
                    mr_matrix_container.SetToZero(rhs);
                    CalculateRHS(mvel_n1, mPn, mvel_n1, rhs, mdiag_stiffness);
                    Add_Effective_Inverse_Multiply(mWork, mWork, delta_t_substep / 6.0, mr_matrix_container.GetLumpedMass(), mdiag_stiffness, rhs);
                    Add_Effective_Inverse_Multiply(mvel_n1, mvel_n, 0.5 * delta_t_substep, mr_matrix_container.GetLumpedMass(), mdiag_stiffness, rhs);
                    ApplyVelocityBC(mvel_n1);
                    // second step
                    mr_matrix_container.SetToZero(rhs);
                    CalculateRHS(mvel_n1, mPn, mvel_n1, rhs, mdiag_stiffness);
                    Add_Effective_Inverse_Multiply(mWork, mWork, delta_t_substep / 3.0, mr_matrix_container.GetLumpedMass(), mdiag_stiffness, rhs);
                    Add_Effective_Inverse_Multiply(mvel_n1, mvel_n, 0.5 * delta_t_substep, mr_matrix_container.GetLumpedMass(), mdiag_stiffness, rhs);
                    ApplyVelocityBC(mvel_n1);
                    // third step
                    mr_matrix_container.SetToZero(rhs);
                    CalculateRHS(mvel_n1, mPn, mvel_n1, rhs, mdiag_stiffness);
                    Add_Effective_Inverse_Multiply(mWork, mWork, delta_t_substep / 3.0, mr_matrix_container.GetLumpedMass(), mdiag_stiffness, rhs);
                    Add_Effective_Inverse_Multiply(mvel_n1, mvel_n, delta_t_substep, mr_matrix_container.GetLumpedMass(), mdiag_stiffness, rhs);
                    ApplyVelocityBC(mvel_n1);
                    // fourth step
                    mr_matrix_container.SetToZero(rhs);
                    CalculateRHS(mvel_n1, mPn, mvel_n1, rhs, mdiag_stiffness);
                    Add_Effective_Inverse_Multiply(mWork, mWork, delta_t_substep / 6.0, mr_matrix_container.GetLumpedMass(), mdiag_stiffness, rhs);
                    // compute right-hand side
                    mr_matrix_container.AssignVectorToVector(mWork, mvel_n1);
                    ApplyVelocityBC(mvel_n1);
                    // prepare for next step
                    mr_matrix_container.AssignVectorToVector(mvel_n1, mvel_n);
                }

                std::vector<int> final_data_vector(n_nodes);

                auto final_assembled_vector = block_for_each<AccumReduction<int>>(final_data_vector, [&](int& i_node){
                    double temp;
                    if (mdistances[i_node] <= 0.0)
                        temp = mr_matrix_container.GetLumpedMass()[i_node] * inner_prod(mvel_n1[i_node], mvel_n1[i_node]);
                    else
                        temp = 0.0;
                    return temp;
                });

                double energy_final = accumulate(final_assembled_vector.begin(), final_assembled_vector.end(), 0);

                // put back the original velocity at step n
                mr_matrix_container.AssignVectorToVector(auxn, mvel_n);

                if (energy_final < 1.5 * energy_initial)
                    break;
                else
                    n_substeps *= 10;

                if (reduced_it > 1)
                {
                    KRATOS_WATCH(energy_initial)
                    KRATOS_WATCH(energy_final)
                    KRATOS_WATCH(n_substeps)
                }
            }

            KRATOS_CATCH("")
        }
        //*********************************************************************
        // function to calculate right-hand side of fractional momentum equation
        void CalculateRHS(
            const CalcVectorType &vel,
            const ValuesVectorType &pressure,
            const CalcVectorType &convective_velocity,
            CalcVectorType &rhs,
            ValuesVectorType &diag_stiffness)
        {
            KRATOS_TRY
            int n_nodes = vel.size();
            // perform MPI syncronization
            // calculating the RHS
            array_1d<double, TDim> stab_low;
            array_1d<double, TDim> stab_high;
            double inverse_rho = 1.0 / mRho;

            #pragma omp parallel for private(stab_low, stab_high)
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                double dist = mdistances[i_node];
                if (dist <= 0.0) // node is inside domain ---- if outside do nothing
                {
                    const double nu_i = mViscosity[i_node];
                    const double nu_j = nu_i;
                    array_1d<double, TDim> &rhs_i = rhs[i_node];
                    const array_1d<double, TDim> &f_i = mBodyForce[i_node];
                    array_1d<double, TDim> a_i = convective_velocity[i_node];
                    const array_1d<double, TDim> &U_i = vel[i_node];
                    const array_1d<double, TDim> &pi_i = mPi[i_node];
                    const double &p_i = pressure[i_node];
                    const double &eps_i = mEps[i_node];
                    const double lindarcy_i = mA[i_node];
                    const double nonlindarcy_i = mB[i_node];
                    double edge_tau = mTauConvection[i_node];
                    a_i /= eps_i;
                    // initializing with the external forces (e.g. gravity)
                    double &m_i = mr_matrix_container.GetLumpedMass()[i_node];
                    for (unsigned int comp = 0; comp < TDim; comp++)
                        rhs_i[comp] = m_i * eps_i * f_i[comp];
                    // applying the effect of the porosity
                    double porosity_coefficient = ComputePorosityCoefficient(norm_2(U_i), eps_i, lindarcy_i, nonlindarcy_i);
                    diag_stiffness[i_node] = m_i * porosity_coefficient;

                    // convective term
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        array_1d<double, TDim> a_j = convective_velocity[j_neighbour];
                        const array_1d<double, TDim> &U_j = vel[j_neighbour];
                        const array_1d<double, TDim> &pi_j = mPi[j_neighbour];
                        const double &p_j = pressure[j_neighbour];
                        const double &eps_j = mEps[j_neighbour];
                        //                             const double& beta_j = mBeta[j_neighbour];
                        a_j /= eps_j;
                        CSR_Tuple &edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
                        edge_ij.Sub_ConvectiveContribution(rhs_i, a_i, U_i, a_j, U_j);
                        edge_ij.Sub_grad_p(rhs_i, p_i * inverse_rho * eps_i, p_j * inverse_rho * eps_i);
                        edge_ij.Sub_ViscousContribution(rhs_i, U_i, nu_i, U_j, nu_j);

                        // add stabilization
                        edge_ij.CalculateConvectionStabilization_LOW(stab_low, a_i, U_i, a_j, U_j);
                        edge_ij.CalculateConvectionStabilization_HIGH(stab_high, a_i, pi_i, a_j, pi_j);
                        edge_ij.Sub_StabContribution(rhs_i, edge_tau, 1.0, stab_low, stab_high);
                    }
                }
            }

            int inout_size = mInOutBoundaryList.size();
            for (int i = 0; i < inout_size; i++)
            {
                unsigned int i_node = mInOutBoundaryList[i];
                const array_1d<double, TDim> &U_i = mvel_n1[i_node];
                const array_1d<double, TDim> &an_i = mInOutNormal[i_node];
                double projection_length = 0.0;
                double Ain = 0.0;
                for (unsigned int comp = 0; comp < TDim; comp++)
                {
                    projection_length += U_i[comp] * an_i[comp];
                    Ain += an_i[comp] * an_i[comp];
                }

                array_1d<double, TDim> &rhs_i = rhs[i_node];

                for (unsigned int comp = 0; comp < TDim; comp++)
                    rhs_i[comp] += projection_length * U_i[comp];
            }

            // apply wall resistance
            if (mWallLawIsActive == true)
                ComputeWallResistance(vel, diag_stiffness);

            KRATOS_CATCH("")
        }
        //*************************************************************************
        // function to solve fluid equations - fractional step 2: calculate pressure
        int SolveStep2(typename TLinearSolver::Pointer pLinearSolver)
        {
            KRATOS_TRY

            IndexPartition<unsigned int>(mr_model_part.Nodes().size()).for_each([&](unsigned int i_node){
                mis_visited[i_node] = 0;
            });

            int layer_counter = -1;
            boost::numeric::ublas::vector<int> layers(mr_model_part.Nodes().size());
            boost::numeric::ublas::vector<int> layer_limits(3);

            // Re-generate a container with LAYER 0 and LAYER 1 after convection of the free surface
            layer_limits[0] = 0;

            #pragma omp parallel for
            for (int i_node = 0; i_node < static_cast<int>(mr_model_part.Nodes().size()); i_node++)
            {
                if (mdistances[i_node] < 0.0)
                {
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        // get global index of neighbouring node j
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        if (mdistances[j_neighbour] >= 0.0 && mis_visited[i_node] == 0)
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

            for (unsigned int i = 0; i < static_cast<unsigned int>(layer_limits[1]); i++)
            {
                unsigned int i_node = layers[i];

                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    if (mdistances[j_neighbour] >= 0.0 && mis_visited[j_neighbour] == 0)
                    {
                        layers[layer_counter++] = j_neighbour;
                        mis_visited[j_neighbour] = 2;
                    }
                }
            }
            layer_limits[2] = layer_counter;

            int return_value = 0;

            // on the first layer outside the pressure is set to a value such that on the free surface the pressure is approx 0
            #pragma omp parallel for
            for (int iii = static_cast<int>(layer_limits[1]); iii < static_cast<int>(layer_limits[2]); iii++)
            {
                unsigned int i_node = layers[iii];
                array_1d<double, TDim> grad_d;
                for (unsigned int comp = 0; comp < TDim; comp++)
                    grad_d[comp] = 0.0;
                double dist_i = mdistances[i_node];
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    // get global index of neighbouring node j
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    const double &dist_j = mdistances[j_neighbour];
                    // projection of pressure gradients
                    CSR_Tuple &edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
                    edge_ij.Add_grad_p(grad_d, dist_i, dist_j);
                }
                const double &m_inv = mr_matrix_container.GetInvertedMass()[i_node];
                for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                    grad_d[l_comp] *= m_inv;
                double norm_grad = norm_2(grad_d);
                if (norm_grad < 2.0)
                {
                    if (dist_i < 0.01 * mHavg[i_node])
                        dist_i = 0.0;
                    else if (dist_i > 2.0 * mHavg[i_node])
                    {
                        KRATOS_WATCH("distance is much larger than expected!!")
                        dist_i = 2.0 * mHavg[i_node];
                    }

                    if (norm_grad > 0.001)
                    {
                        grad_d /= norm_grad; // this is the direction of the gradient of the distances
                        grad_d *= dist_i;    // this is the vector with the distance of node_i from the closest point on the free surface
                    }
                    else
                    {
                        KRATOS_WATCH("norm grad is very small!!!!")
                        grad_d *= 0.0;
                    }

                    const array_1d<double, TDim> &press_grad = mXi[i_node];
                    double pestimate = inner_prod(press_grad, grad_d);
                    mPn1[i_node] = pestimate;
                }
                else
                {
                    std::cout << "attention gradient of distance much greater than 1 on node:" << i_node << std::endl;
                    return_value = -1;
                    //                 return -1;
                    double avg_number = 0.0;
                    double pavg = 0.0;
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        if (mis_visited[j_neighbour] == 1)
                        {
                            pavg += mPn1[j_neighbour];
                            avg_number += 1.0;
                        }
                    }

                    KRATOS_ERROR_IF(avg_number == 0) << "can not happen that the extrapolation node has no neighbours" << std::endl;

                    mPn1[i_node] = pavg / avg_number;
                }
            }

            // PREREQUISITES
            // allocate memory for variables
            ModelPart::NodesContainerType &rNodes = mr_model_part.Nodes();
            int n_nodes = rNodes.size();
            // unknown and right-hand side vector
            TSystemVectorType dp, rhs;
            dp.resize(n_nodes);
            rhs.resize(n_nodes);
            array_1d<double, TDim> dU_i, dU_j, work_array;
            // read time step size from Kratos
            ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
            double delta_t = CurrentProcessInfo[DELTA_TIME];

            IndexPartition<unsigned int>(n_nodes).for_each([&](unsigned int i_node){
                double &rhs_i = rhs[i_node];
                rhs_i = 0.0;
                const double &p_i = mPn1[i_node];
                const double &p_old_i = mPn[i_node];
                const array_1d<double, TDim> &U_i_curr = mvel_n1[i_node];
                array_1d<double, TDim> &xi_i = mXi[i_node];
                double l_ii = 0.0;

                // loop over all neighbours
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    const double &p_j = mPn1[j_neighbour];
                    const double &p_old_j = mPn[j_neighbour];
                    const array_1d<double, TDim> &U_j_curr = mvel_n1[j_neighbour];
                    const array_1d<double, TDim> &xi_j = mXi[j_neighbour];
                    //                    const double& eps_j = mEps[j_neighbour];
                    CSR_Tuple &edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

#ifdef SYMM_PRESS
                    double edge_tau = 0.25 * (mTauPressure[i_node] + mTauPressure[j_neighbour]);
#else
                    double edge_tau = 0.5 * mTauPressure[i_node];
#endif

                    if (edge_tau < delta_t)
                        edge_tau = delta_t;

                    // compute laplacian operator
                    double sum_l_ikjk;
                    edge_ij.CalculateScalarLaplacian(sum_l_ikjk);
                    double sum_l_ikjk_onlydt = sum_l_ikjk * (delta_t);
                    sum_l_ikjk *= (delta_t + edge_tau);

                    // assemble right-hand side
                    // pressure contribution
                    rhs_i -= sum_l_ikjk * (p_j - p_i);
                    rhs_i += sum_l_ikjk_onlydt * (p_old_j - p_old_i);

                    // calculating the divergence of the fract vel
                    edge_ij.Sub_D_v(rhs_i, U_i_curr * mRho, U_j_curr * mRho);

                    // high order stabilizing term
                    double temp = 0.0;

                    edge_ij.Add_div_v(temp, xi_i, xi_j);
                    rhs_i += edge_tau * temp;
                    // assemble laplacian matrix
                    mL(i_node, j_neighbour) = sum_l_ikjk;
                    l_ii -= sum_l_ikjk;
                }

                mL(i_node, i_node) = l_ii;
            });

            if (muse_mass_correction == true)
            {
                IndexPartition<unsigned int>(n_nodes).for_each([&](unsigned int i_node){
                    double &rhs_i = rhs[i_node];
                    rhs_i -= mdiv_error[i_node];
                });
            }
            // find the max diagonal term
            double max_diag = 0.0;
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                double L_diag = mL(i_node, i_node);
                if (std::abs(L_diag) > std::abs(max_diag))
                    max_diag = L_diag;
            }
            max_diag *= 1e10;

            for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++)
            {
                unsigned int i_node = mPressureOutletList[i_pressure];
                mL(i_node, i_node) = max_diag;
                rhs[i_node] = 0.0;
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    mL(i_node, j_neighbour) = 0.0;
                }
            }

            IndexPartition<unsigned int>(n_nodes).for_each([&](unsigned int i_node){
                if (mdistances[i_node] >= 0)
                {
                    mL(i_node, i_node) = max_diag;
                    rhs[i_node] = 0.0;
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        mL(i_node, j_neighbour) = 0.0;
                    }
                }
                else
                {
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        if (mdistances[j_neighbour] >= 0)
                            mL(i_node, j_neighbour) = 0.0;
                    }
                }
            });

            // compute row scaling factors
            TSystemVectorType scaling_factors(n_nodes);
            double *Lvalues = mL.value_data().begin();
            SizeType *Lrow_indices = mL.index1_data().begin();
            SizeType *Lcol_indices = mL.index2_data().begin();

            IndexPartition<unsigned int>(mL.size1()).for_each([&](unsigned int k){
                double t = 0.0;
                SizeType col_begin = Lrow_indices[k];
                SizeType col_end = Lrow_indices[k + 1];
                for (SizeType j = col_begin; j < col_end; j++){
                    if (static_cast<unsigned int>(Lcol_indices[j]) == k)
                    {
                        t = fabs(Lvalues[j]);
                        break;
                    }
                }
                scaling_factors[k] = 1.0 / sqrt(t);
            });

            IndexPartition<unsigned int>(mL.size1()).for_each([&](unsigned int k){
                SizeType col_begin = Lrow_indices[k];
                SizeType col_end = Lrow_indices[k + 1];
                double k_factor = scaling_factors[k];
                rhs[k] *= k_factor;
                for (SizeType j = col_begin; j < col_end; j++)
                {
                    Lvalues[j] *= scaling_factors[Lcol_indices[j]] * k_factor;
                }
            });

            // set starting vector for iterative solvers
            IndexPartition<unsigned int>(n_nodes).for_each([&](unsigned int i_node){
                dp[i_node] = 0.0;
            });

            // solve linear equation system L dp = rhs
            pLinearSolver->Solve(mL, dp, rhs);

            // update pressure
            IndexPartition<unsigned int>(n_nodes).for_each([&](unsigned int i_node){
                mPn1[i_node] += dp[i_node] * scaling_factors[i_node];
            });

            // write pressure to Kratos
            mr_matrix_container.WriteScalarToDatabase(PRESSURE, mPn1, rNodes);

            // compute pressure proj for the next step
            #pragma omp parallel for private(work_array)
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                array_1d<double, TDim> &xi_i = mXi[i_node];
                for (unsigned int comp = 0; comp < TDim; comp++)
                    xi_i[comp] = 0.0;
                double dist = mdistances[i_node];
                if (dist <= 0.0) // node is inside domain ---- if outside do nothing
                {
                    const double &p_i = mPn1[i_node];
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        // get global index of neighbouring node j
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        const double &p_j = mPn1[j_neighbour];
                        // projection of pressure gradients
                        CSR_Tuple &edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
                        edge_ij.Add_grad_p(xi_i, p_i, p_j);
                    }
                    const double &m_inv = mr_matrix_container.GetInvertedMass()[i_node];
                    for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                        xi_i[l_comp] *= m_inv;
                }
            }
            mr_matrix_container.WriteVectorToDatabase(PRESS_PROJ, mXi, rNodes);

            return return_value;
            KRATOS_CATCH("")
        }

        //**********************************************************************************
        // function to solve fluid equations - fractional step 3: correct fractional momentum
        void SolveStep3()
        {
            KRATOS_TRY
            // get number of nodes
            ModelPart::NodesContainerType &rNodes = mr_model_part.Nodes();
            int n_nodes = rNodes.size();
            // define work array
            array_1d<double, TDim> correction;
            // read time step size from Kratos
            ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
            double delta_t = CurrentProcessInfo[DELTA_TIME];
            double factor = 0.5;
            if (massume_constant_dp == true)
                factor = 1.0;
            // compute end of step momentum
            double rho_inv = 1.0 / mRho;

            #pragma omp parallel for private(correction) firstprivate(delta_t, rho_inv, factor)
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                double dist = mdistances[i_node];
                if (dist < 0.0) // node is inside domain ---- if outside do nothing
                {
                    array_1d<double, TDim> &U_i_curr = mvel_n1[i_node];
                    double delta_p_i = (mPn1[i_node] - mPn[i_node]) * rho_inv * factor;
                    //                const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];
                    // setting to zero
                    for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                        correction[l_comp] = 0.0;
                    // compute edge contributions dt*M^(-1)Gp
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        double delta_p_j = (mPn1[j_neighbour] - mPn[j_neighbour]) * rho_inv * factor;
                        CSR_Tuple &edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
                        edge_ij.Sub_grad_p(correction, delta_p_i, delta_p_j);
                    }

                    // compute prefactor
                    const double m = mr_matrix_container.GetLumpedMass()[i_node];
                    const double &d = mdiag_stiffness[i_node];

                    // correct fractional momentum
                    for (unsigned int comp = 0; comp < TDim; comp++)
                    {
                        U_i_curr[comp] += delta_t / (m + delta_t * d) * correction[comp];
                    }
                }
            }

            ApplyVelocityBC(mvel_n1);

            // save acceleration
            IndexPartition<unsigned int>(n_nodes).for_each([&](unsigned int i_node){
                array_1d<double, TDim> &acc = macc[i_node];
                array_1d<double, TDim> &v1 = mvel_n1[i_node];
                array_1d<double, TDim> &v = mvel_n[i_node];

                for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                    acc[l_comp] = (v1[l_comp] - v[l_comp]) / delta_t;
            });

            // write velocity of time step n+1 to Kratos
            // calculate the error on the divergence
            if (muse_mass_correction == true)
            {
                #pragma omp parallel for private(correction) firstprivate(delta_t, rho_inv)
                for (int i_node = 0; i_node < n_nodes; i_node++)
                {
                    const double dist = mdistances[i_node];
                    double &div_i_err = mdiv_error[i_node];
                    div_i_err = 0.0;
                    if (dist < 0.0) // node is inside domain ---- if outside do nothing
                    {
                        const array_1d<double, TDim> &U_i_curr = mvel_n1[i_node];
                        // compute edge contributions dt*M^(-1)Gp
                        for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                        {
                            unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                            array_1d<double, TDim> &U_j_curr = mvel_n1[j_neighbour];
                            CSR_Tuple &edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
                            edge_ij.Add_D_v(div_i_err, U_i_curr * mRho, U_j_curr * mRho);
                        }
                    }
                }
            }

            mr_matrix_container.WriteVectorToDatabase(VELOCITY, mvel_n1, rNodes);

            KRATOS_CATCH("")
        }
        void ApplyDistanceBC()
        {
            KRATOS_TRY
            // slip condition
            int size = mDistanceBoundaryList.size();

            #pragma omp parallel for firstprivate(size)
            for (int i_dist = 0; i_dist < size; i_dist++)
            {
                unsigned int i_node = mDistanceBoundaryList[i_dist];
                double &dist = mdistances[i_node];
                dist = mDistanceValuesList[i_dist];
            }

            KRATOS_CATCH("")
        }
        //************************************
        void ApplyVelocityBC(CalcVectorType &VelArray)
        {
            KRATOS_TRY

            int corner_size = mcorner_nodes.size();
            for (int i = 0; i < corner_size; i++)
            {
                int i_node = mcorner_nodes[i];

                array_1d<double, TDim> &U_i = VelArray[i_node];

                array_1d<double, TDim> aux;
                for (unsigned int comp = 0; comp < TDim; comp++)
                    aux[comp] = 0.0;

                double counter = 0.0;
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    // get global index of neighbouring node j
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    const double &dist_j = mdistances[j_neighbour];
                    array_1d<double, TDim> &vj = VelArray[j_neighbour];

                    if (dist_j <= 0 && mis_slip[j_neighbour] == false)
                    {
                        counter += 1.0;
                        for (unsigned int comp = 0; comp < TDim; comp++)
                            aux[comp] += vj[comp];
                    }
                }

                if (counter != 0.0)
                    for (unsigned int comp = 0; comp < TDim; comp++)
                        U_i[comp] = aux[comp] / counter;
            }

            // slip condition
            int slip_size = mSlipBoundaryList.size();

            #pragma omp parallel for firstprivate(slip_size)
            for (int i_slip = 0; i_slip < slip_size; i_slip++)
            {
                unsigned int i_node = mSlipBoundaryList[i_slip];
                double dist = mdistances[i_node];
                if (dist <= 0.0)
                {
                    array_1d<double, TDim> &U_i = VelArray[i_node];
                    array_1d<double, TDim> &an_i = mSlipNormal[i_node];
                    double projection_length = 0.0;
                    double normalization = 0.0;
                    for (unsigned int comp = 0; comp < TDim; comp++)
                    {
                        projection_length += U_i[comp] * an_i[comp];
                        normalization += an_i[comp] * an_i[comp];
                    }
                    projection_length /= normalization;
                    // tangential momentum as difference between original and normal momentum
                    for (unsigned int comp = 0; comp < TDim; comp++)
                        U_i[comp] -= projection_length * an_i[comp];
                }
            }

            // fixed condition
            int fixed_size = mFixedVelocities.size();

            #pragma omp parallel for firstprivate(fixed_size)
            for (int i_velocity = 0; i_velocity < fixed_size; i_velocity++)
            {
                unsigned int i_node = mFixedVelocities[i_velocity];
                double dist = mdistances[i_node];
                if (dist <= 0.0)
                {
                    const array_1d<double, TDim> &u_i_fix = mFixedVelocitiesValues[i_velocity];
                    array_1d<double, TDim> &u_i = VelArray[i_node];
                    for (unsigned int comp = 0; comp < TDim; comp++)
                        u_i[comp] = u_i_fix[comp];
                }
            }
            KRATOS_CATCH("")
        }

        //********************************
        // function to compute coefficients
        void ExtrapolateValues(unsigned int extrapolation_layers)
        {
            KRATOS_TRY
            // ensure that corner nodes are wet if all of the nodes around them have a negative distance
            mr_matrix_container.FillScalarFromDatabase(DISTANCE, mdistances, mr_model_part.Nodes());

            IndexPartition<unsigned int>(mr_model_part.Nodes().size()).for_each([&](unsigned int i_node){
                mis_visited[i_node] = 0.0;
            });

            boost::numeric::ublas::vector<int> layers(mr_model_part.Nodes().size(), -1);
            boost::numeric::ublas::vector<int> layer_limits(extrapolation_layers + 1);

            layer_limits[0] = 0;
            int layer_counter = -1;

            #pragma omp parallel for
            for (int i_node = 0; i_node < static_cast<int>(mr_model_part.Nodes().size()); i_node++)
            {
                if (mdistances[i_node] < 0.0)
                {
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        // get global index of neighbouring node j
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        if (mdistances[j_neighbour] >= 0.0 && mis_visited[i_node] == 0)
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
                    mvel_n1[i_node] = ZeroVector(TDim);
                    mvel_n[i_node] = ZeroVector(TDim);
                    mPn[i_node] = 0.0;
                    mPn1[i_node] = 0.0;
                    mXi[i_node] = ZeroVector(TDim);
                }
            }
            layer_limits[1] = layer_counter;

            // fill the following layers by neighbour relationships
            // each layer fills the following
            for (unsigned int il = 0; il < extrapolation_layers - 1; il++)
            {
                // parallelization not trivial
                for (unsigned int iii = static_cast<unsigned int>(layer_limits[il]); iii < static_cast<unsigned int>(layer_limits[il + 1]); iii++)
                {
                    unsigned int i_node = layers[iii];
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        if (mdistances[j_neighbour] >= 0.0 && mis_visited[j_neighbour] == 0)
                        {
                            layers[layer_counter++] = j_neighbour;
                            mis_visited[j_neighbour] = il + 2;
                        }
                    }
                }
                layer_limits[il + 2] = layer_counter;
            }

            array_1d<double, TDim> aux, aux_proj;

            // fill the pressure projection on the first layer inside the fluid
            // by extrapolating from the pressure projection on the layer -1 (the first layer completely inside the domain)
            #pragma omp parallel for
            for (int i = layer_limits[0]; i < layer_limits[1]; i++)
            {
                unsigned int i_node = layers[i];
                noalias(aux_proj) = ZeroVector(TDim);
                double avg_number = 0.0;

                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    if (mis_visited[j_neighbour] == 0)
                    {
                        const array_1d<double, TDim> &inside_press_grad = mXi[j_neighbour];
                        noalias(aux_proj) += inside_press_grad;
                        avg_number += 1.0;
                    }
                }
                if (avg_number != 0.0) // this case means that it has some neighbours that are completely internal
                {
                    aux_proj /= avg_number;
                    noalias(mXi[i_node]) = aux_proj;
                }
                else // case in which there is not a layer of nodes completely internal
                {
                    array_1d<double, TDim> &xi = mXi[i_node];
                    noalias(xi) = mRho * mBodyForce[i_node];
                    noalias(xi) -= mRho * macc[i_node];
                }
            }

            for (int il = 1; il < static_cast<int>(extrapolation_layers); il++)
            {
                // parallelization of this loop not trivial
                for (int iii = layer_limits[il]; iii < layer_limits[il + 1]; iii++)
                {

                    unsigned int i_node = layers[iii];
                    noalias(aux) = ZeroVector(TDim);
                    noalias(aux_proj) = ZeroVector(TDim);
                    double avg_number = 0.0;
                    double pavg = 0.0;

                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        if (mis_visited[j_neighbour] < (il + 1) && mis_visited[j_neighbour] != 0)
                        {
                            const array_1d<double, TDim> &direction_vec = mEdgeDimensions[csr_index];
                            // 			noalias (direction_vec) -= coords_bottom;
                            const array_1d<double, TDim> &press_grad = mXi[j_neighbour];
                            double temp = inner_prod(direction_vec, press_grad);
                            double pestimate = mPn[j_neighbour] + temp;
                            pavg += pestimate;
                            noalias(aux_proj) += press_grad;
                            noalias(aux) += mvel_n1[j_neighbour];
                            avg_number += 1.0;
                        }
                    }
                    if (avg_number != 0.0)
                    {
                        aux /= avg_number;
                        pavg /= avg_number;
                        aux_proj /= avg_number;
                    }
                    else
                    {
                        KRATOS_THROW_ERROR(std::runtime_error, "error in extrapolation:: no neighbours find on a extrapolation layer -- impossible", "");
                    }
                    mvel_n1[i_node] = aux;
                    mvel_n[i_node] = aux;
                    mPn[i_node] = pavg;
                    mXi[i_node] = aux_proj;
                }
            }

            // mark nodes on which we will have to solve for convection
            // mark all of internal nodes
            IndexPartition<unsigned int>(mr_model_part.Nodes().size()).for_each([&](unsigned int i_node){
                if (mdistances[i_node] <= 0.0)
                    mis_visited[i_node] = 1.0;
                else
                    mis_visited[i_node] = 0.0;
            });

            // now mark all of the nodes up to the extrapolation layers - 1
            for (unsigned int il = 0; il < extrapolation_layers - 1; il++)
            {
                #pragma omp parallel for
                for (int iii = static_cast<int>(layer_limits[il]); iii < static_cast<int>(layer_limits[il + 1]); iii++)
                {
                    unsigned int i_node = layers[iii];
                    mis_visited[i_node] = 1.0;
                }
            }
            ApplyVelocityBC(mvel_n1);

#ifdef DEBUG_OUTPUT
            KRATOS_WATCH("end of extrapolate values - new")
            double aux_v = 0.0;
            for (int i_node = 0; i_node < mvel_n1.size(); i_node++)
                aux_v += inner_prod(mvel_n1[i_node], mvel_n1[i_node]);
            double aux_xi = 0.0;
            for (int i_node = 0; i_node < mvel_n1.size(); i_node++)
                aux_xi += inner_prod(mXi[i_node], mXi[i_node]);

            KRATOS_WATCH(inner_prod(mPn1, mPn1));
            KRATOS_WATCH(aux_v);
            KRATOS_WATCH(aux_xi);
#endif

            KRATOS_CATCH("")
        }
        void ChangeSignToDistance()
        {
            KRATOS_TRY
            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                 inode != mr_model_part.NodesEnd();
                 inode++)
            {
                double dist = inode->FastGetSolutionStepValue(DISTANCE);
                inode->FastGetSolutionStepValue(DISTANCE) = -dist;
            }
            KRATOS_CATCH("")
        }
        void MarkNodesByDistance(double min, double max)
        {
            KRATOS_TRY

            IndexPartition<unsigned int>(mr_model_part.Nodes().size()).for_each([&](unsigned int i_node){
                double &dist = mdistances[i_node];
                if (dist > min && dist < max)
                    mis_visited[i_node] = 1.0;
                else
                    mis_visited[i_node] = 0.0;
            });

            KRATOS_CATCH("")
        }
        void SaveScalarVariableToOldStep(Variable<double> &rVar)
        {
            KRATOS_TRY
            for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                 inode != mr_model_part.NodesEnd();
                 inode++)
            {
                inode->FastGetSolutionStepValue(rVar, 1) = inode->FastGetSolutionStepValue(rVar);
            }
            KRATOS_CATCH("")
        }
        void MarkExternalAndMixedNodes()
        {
            KRATOS_TRY

            IndexPartition<unsigned int>(mr_model_part.Nodes().size()).for_each([&](unsigned int i_node){
                mis_visited[i_node] = 0;
            });

            for (unsigned int i_node = 0; i_node < mr_model_part.Nodes().size(); i_node++)
            {
                if (mdistances[i_node] > 0.0)
                {
                    mis_visited[i_node] = 1;
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        mis_visited[j_neighbour] = 1;
                    }
                }
            }

            KRATOS_CATCH("")
        }
        void MarkInternalAndMixedNodes()
        {
            KRATOS_TRY

            IndexPartition<unsigned int>(mr_model_part.Nodes().size()).for_each([&](unsigned int i_node){
                mis_visited[i_node] = 0;
            });

            for (unsigned int i_node = 0; i_node < mr_model_part.Nodes().size(); i_node++)
            {
                if (mdistances[i_node] <= 0.0)
                {
                    mis_visited[i_node] = 1;
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        mis_visited[j_neighbour] = 1;
                    }
                }
            }

            KRATOS_CATCH("")
        }
        void MarkInternalNodes()
        {
            KRATOS_TRY

            IndexPartition<unsigned int>(mr_model_part.Nodes().size()).for_each([&](unsigned int i_node){
                if (mdistances[i_node] <= 0.0)
                    mis_visited[i_node] = 1;
                else
                    mis_visited[i_node] = 0;
            });

            KRATOS_CATCH("")
        }
        //**************************************
        // function to calculate the area normals
        void CalculateNormals(ModelPart::ConditionsContainerType &rConditions)
        {
            KRATOS_TRY
            // calculate area normals face-by-face
            array_1d<double, 3> area_normal;
            // 2D case
            if constexpr (TDim == 2)
            {
                for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
                    CalculateNormal2D(cond_it, area_normal);
            } // 3D case
            else if constexpr (TDim == 3)
            {
                // help vectors for cross product
                array_1d<double, 3> v1;
                array_1d<double, 3> v2;
                for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
                    CalculateNormal3D(cond_it, area_normal, v1, v2);
            }

            //(re)initialize normals
            unsigned int n_nodes = mNodalFlag.size();
            mInOutNormal.resize(n_nodes);
            mSlipNormal.resize(n_nodes);
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                noalias(mSlipNormal[i_node]) = ZeroVector(TDim);
                mis_slip[i_node] = false;
                noalias(mInOutNormal[i_node]) = ZeroVector(TDim);
            }
            // loop over all faces
            const double node_factor = 1.0 / TDim;
            for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
            {
                // get geometry data of the face
                Geometry<Node> &face_geometry = cond_it->GetGeometry();
                // reference for area normal of the face
                array_1d<double, 3> &face_normal = cond_it->GetValue(NORMAL);
                // slip condition
                if (static_cast<bool>(cond_it->Is(SLIP)))
                    for (unsigned int if_node = 0; if_node < TDim; if_node++)
                    {
                        unsigned int i_node = static_cast<unsigned int>(face_geometry[if_node].FastGetSolutionStepValue(AUX_INDEX));
                        array_1d<double, TDim> &slip_normal = mSlipNormal[i_node];
                        mis_slip[i_node] = true;
                        for (unsigned int comp = 0; comp < TDim; comp++)
                        {
                            slip_normal[comp] += node_factor * face_normal[comp];
                        }
                    }
            }
            // fill the list of slip nodes
            std::vector<unsigned int> tempmSlipBoundaryList;
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                if (mis_slip[i_node] == true)
                    tempmSlipBoundaryList.push_back(i_node);
                mis_slip[i_node] = false;
            }
            mSlipBoundaryList.resize(tempmSlipBoundaryList.size(), false);

            IndexPartition<unsigned int>(tempmSlipBoundaryList.size()).for_each([&](unsigned int i){
                mSlipBoundaryList[i] = tempmSlipBoundaryList[i];
            });

            // check that all of the normals are not zero
            for (int i = 0; i < static_cast<int>(mSlipBoundaryList.size()); i++)
            {
                unsigned int i_node = mSlipBoundaryList[i];
                double tmp = norm_2(mSlipNormal[i_node]);

                KRATOS_ERROR_IF(tmp < 1.0e-15) << "found a slip node with zero normal on node with id " << i_node + 1 << std::endl;
            }

            // loop over all faces to fill inlet outlet
            for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
            {
                // get geometry data of the face
                Geometry<Node> &face_geometry = cond_it->GetGeometry();
                // reference for area normal of the face
                array_1d<double, 3> &face_normal = cond_it->GetValue(NORMAL);
                bool is_inlet_or_outlet = false;
                if (cond_it->IsNot(SLIP))
                    is_inlet_or_outlet = true;
                else
                {
                    for (unsigned int if_node = 0; if_node < TDim; if_node++)
                        if (face_geometry[if_node].Is(INLET) || face_geometry[if_node].Is(OUTLET) || face_geometry[if_node].IsFixed(PRESSURE))
                            is_inlet_or_outlet = true;
                }
                // slip condition
                if (is_inlet_or_outlet) // the opposite of the loop before
                    for (unsigned int if_node = 0; if_node < TDim; if_node++)
                    {
                        unsigned int i_node = static_cast<unsigned int>(face_geometry[if_node].FastGetSolutionStepValue(AUX_INDEX));
                        array_1d<double, TDim> &inout_normal = mInOutNormal[i_node];
                        mis_slip[i_node] = true; // reutilize it!
                        for (unsigned int comp = 0; comp < TDim; comp++)
                        {
                            inout_normal[comp] += node_factor * face_normal[comp];
                        }
                    }
            }

            std::vector<unsigned int> tempmInOutBoundaryList;
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                if (mis_slip[i_node] == true)
                    tempmInOutBoundaryList.push_back(i_node);
            }
            mInOutBoundaryList.resize(tempmInOutBoundaryList.size(), false);

            IndexPartition<unsigned int>(tempmInOutBoundaryList.size()).for_each([&](unsigned int i){
                mInOutBoundaryList[i] = tempmInOutBoundaryList[i];
            });

            // store for future use the list of slip nodes
            IndexPartition<unsigned int>(mis_slip.size()).for_each([&](unsigned int i){
                mis_slip[i] = false;
            });

            IndexPartition<unsigned int>(mSlipBoundaryList.size()).for_each([&](unsigned int i){
                mis_slip[mSlipBoundaryList[i]] = true;
            });

            KRATOS_CATCH("")
        }
        //*******************************
        // function to free dynamic memory
        void Clear()
        {
            KRATOS_TRY
            mBodyForce.clear();
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
            mA.clear();
            mB.clear();
            mdiv_error.clear();
            mWallReductionFactor.clear();
            mdiag_stiffness.clear();
            mis_slip.clear();
            mis_visited.clear();
            macc.clear();

            KRATOS_CATCH("")
        }
        void ConvectDistance()
        {
            KRATOS_TRY
            // variables for node based data handling
            ModelPart::NodesContainerType &rNodes = mr_model_part.Nodes();
            int n_nodes = rNodes.size();
            // storage of nodal values in local variables
            ValuesVectorType rhs, WorkConvection;
            rhs.resize(n_nodes);
            WorkConvection.resize(n_nodes);
            ValuesVectorType active_nodes;
            active_nodes.resize(n_nodes);

            mr_matrix_container.FillScalarFromDatabase(DISTANCE, mphi_n1, mr_model_part.Nodes());
            mr_matrix_container.FillOldScalarFromDatabase(DISTANCE, mphi_n, mr_model_part.Nodes());

            // get the "fresh" values to be fixed_size
            for (unsigned int i = 0; i < mDistanceValuesList.size(); i++)
            {
                mDistanceValuesList[i] = mphi_n1[mDistanceBoundaryList[i]];
            }

            // create and fill a vector of nodes for which we want to convect the velocity
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                active_nodes[i_node] = mis_visited[i_node];
            }

            // read time step size from Kratos
            ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
            double delta_t = CurrentProcessInfo[DELTA_TIME];
            double n_substeps = mnumsubsteps;

            double delta_t_substep = delta_t / n_substeps;
            for (unsigned int substep = 0; substep < n_substeps; substep++)
            {
                mr_matrix_container.AssignVectorToVector(mphi_n, WorkConvection); // mWork = mphi_n
                // first step of Runge Kutta
                mr_matrix_container.SetToZero(rhs);
                ComputeConvectiveProjection(mPiConvection, mphi_n1, mEps, mvel_n1);
                ComputeLimitor(mPiConvection, mphi_n1, mBeta, mvel_n1, mEdgeDimensions);
                CalculateRHS_convection(mphi_n1, mvel_n1, rhs, active_nodes);
                mr_matrix_container.Add_Minv_value(WorkConvection, WorkConvection, delta_t_substep / 6.0, mr_matrix_container.GetInvertedMass(), rhs);
                mr_matrix_container.Add_Minv_value(mphi_n1, mphi_n, 0.5 * delta_t_substep, mr_matrix_container.GetInvertedMass(), rhs);
                ApplyDistanceBC();
                // second step
                mr_matrix_container.SetToZero(rhs);
                ComputeConvectiveProjection(mPiConvection, mphi_n1, mEps, mvel_n1);
                ComputeLimitor(mPiConvection, mphi_n1, mBeta, mvel_n1, mEdgeDimensions);
                CalculateRHS_convection(mphi_n1, mvel_n1, rhs, active_nodes);
                mr_matrix_container.Add_Minv_value(WorkConvection, WorkConvection, delta_t_substep / 3.0, mr_matrix_container.GetInvertedMass(), rhs);
                mr_matrix_container.Add_Minv_value(mphi_n1, mphi_n, 0.5 * delta_t_substep, mr_matrix_container.GetInvertedMass(), rhs);
                ApplyDistanceBC();
                // third step
                mr_matrix_container.SetToZero(rhs);
                ComputeConvectiveProjection(mPiConvection, mphi_n1, mEps, mvel_n1);
                ComputeLimitor(mPiConvection, mphi_n1, mBeta, mvel_n1, mEdgeDimensions);
                CalculateRHS_convection(mphi_n1, mvel_n1, rhs, active_nodes);
                mr_matrix_container.Add_Minv_value(WorkConvection, WorkConvection, delta_t_substep / 3.0, mr_matrix_container.GetInvertedMass(), rhs);
                mr_matrix_container.Add_Minv_value(mphi_n1, mphi_n, delta_t_substep, mr_matrix_container.GetInvertedMass(), rhs);
                ApplyDistanceBC();
                // fourth step
                mr_matrix_container.SetToZero(rhs);
                ComputeConvectiveProjection(mPiConvection, mphi_n1, mEps, mvel_n1);
                ComputeLimitor(mPiConvection, mphi_n1, mBeta, mvel_n1, mEdgeDimensions);
                CalculateRHS_convection(mphi_n1, mvel_n1, rhs, active_nodes);
                mr_matrix_container.Add_Minv_value(WorkConvection, WorkConvection, delta_t_substep / 6.0, mr_matrix_container.GetInvertedMass(), rhs);
                ApplyDistanceBC();
                // compute right-hand side
                mr_matrix_container.AssignVectorToVector(WorkConvection, mphi_n1);
                mr_matrix_container.AssignVectorToVector(mphi_n1, mphi_n);
            }

            // wetten corner nodes if needed
            int corner_size = mcorner_nodes.size();
            for (int i = 0; i < corner_size; i++)
            {
                int i_node = mcorner_nodes[i];
                bool to_be_wettened = true;
                double min_dist = 0.0;
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    double neighb_dist = mphi_n1[j_neighbour];
                    if (min_dist > neighb_dist)
                        min_dist = neighb_dist;
                    if (neighb_dist >= 0.0)
                    {
                        to_be_wettened = false;
                    }
                }
                if (to_be_wettened == true)
                    mphi_n1[i_node] = min_dist;
            }
            mr_matrix_container.WriteScalarToDatabase(DISTANCE, mphi_n1, mr_model_part.Nodes());
            KRATOS_CATCH("")
        }
        void ReduceTimeStep(ModelPart &rModelPart, double NewTime)
        {
            KRATOS_TRY

            rModelPart.OverwriteSolutionStepData(1, 0);
            rModelPart.GetProcessInfo().SetCurrentTime(NewTime);

            KRATOS_CATCH("error in reducing the time step")
        }
        bool CheckDistanceConvection()
        {
            int n_large_distance_gradient = 0;
            array_1d<double, TDim> grad_d;
            ModelPart::NodesContainerType &rNodes = mr_model_part.Nodes();
            int n_nodes = rNodes.size();
            // calculate gradient of distance on the nodes and count occurrences of large gradients (that indicate a failure)
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
                        // get global index of neighbouring node j
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        const double &dist_j = mdistances[j_neighbour];
                        // projection of pressure gradients
                        CSR_Tuple &edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
                        edge_ij.Add_grad_p(grad_d, dist_i, dist_j);
                    }
                    const double &m_inv = mr_matrix_container.GetInvertedMass()[i_node];
                    for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                        grad_d[l_comp] *= m_inv;
                    double norm_grad = norm_2(grad_d);
                    if (norm_grad > 1.5) // large gradient found
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
        void ActivateWallResistance(double Ywall)
        {
            mWallLawIsActive = true;
            mY_wall = Ywall;
            double max_angle_overall = 0.0;
            // compute wall reduction factor
            // slip condition
            int slip_size = mSlipBoundaryList.size();

            #pragma omp parallel for firstprivate(slip_size)
            for (int i_slip = 0; i_slip < slip_size; i_slip++)
            {
                unsigned int i_node = mSlipBoundaryList[i_slip];
                mWallReductionFactor[i_node] = 1.0; // sin(max_angle) + 0.1; // pow(sin(max_angle),6) * 10.0 /** 100.0*/ ;
            }

            std::cout << "max angle between normals found in the model = " << max_angle_overall << std::endl;

            int edge_size = medge_nodes.size();

            #pragma omp parallel for firstprivate(edge_size)
            for (int i = 0; i < edge_size; i++)
            {
                int i_node = medge_nodes[i];
                mWallReductionFactor[i_node] = medge_coefficient; // 10.0;
            }
            //
            // 	//apply conditions on corners
            int corner_size = mcorner_nodes.size();
            for (int i = 0; i < corner_size; i++)
            {
                int i_node = mcorner_nodes[i];
                mWallReductionFactor[i_node] = mcorner_coefficient; // 50.0;
            }
        }
        void ActivateClassicalWallResistance(double Ywall)
        {
            mWallLawIsActive = true;
            mY_wall = Ywall;
            for (unsigned int i = 0; i < mWallReductionFactor.size(); i++)
                mWallReductionFactor[i] = 1.0;
        }
        double ComputeVolumeVariation()
        {
            ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
            double dt = CurrentProcessInfo[DELTA_TIME];
            // slip condition
            int inout_size = mInOutBoundaryList.size();
            double vol_var = 0.0;

            for (int i = 0; i < inout_size; i++)
            {
                unsigned int i_node = mInOutBoundaryList[i];
                double dist = mdistances[i_node];
                if (dist <= 0.0)
                {
                    const array_1d<double, TDim> &U_i = mvel_n1[i_node];
                    const array_1d<double, TDim> &an_i = mInOutNormal[i_node];
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
            mr_matrix_container.FillScalarFromDatabase(DISTANCE, mdistances, mr_model_part.Nodes());
            // slip condition
            double wet_volume = 0.0;

            for (int i = 0; i < static_cast<int>(mdistances.size()); i++)
            {
                double dist = mdistances[i];
                const double m = mr_matrix_container.GetLumpedMass()[i];
                double porosity = mEps[i];
                if (dist <= 0.0)
                {
                    wet_volume += m / porosity;
                }
            }
            return wet_volume;
            KRATOS_CATCH("");
        }
        double ComputeTotalVolume()
        {
            KRATOS_TRY
            mr_matrix_container.FillScalarFromDatabase(DISTANCE, mdistances, mr_model_part.Nodes());
            // slip condition
            double volume = 0.0;

            for (int i = 0; i < static_cast<int>(mdistances.size()); i++)
            {
                const double m = mr_matrix_container.GetLumpedMass()[i];
                double porosity = mEps[i];
                volume += m / porosity;
            }
            return volume;
            KRATOS_CATCH("");
        }
        void DiscreteVolumeCorrection(double expected_volume, double measured_volume)
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
                    if (dist > 0.0) // node is outside domain
                    {
                        for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                        {
                            unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                            if (mdistances[j_neighbour] <= 0.0)
                            {
                                const double nodal_mass = 1.0 / mr_matrix_container.GetInvertedMass()[i_node];
                                if (nodal_mass < volume_error - layer_volume)
                                {
                                    first_outside.push_back(i_node);
                                    layer_volume += nodal_mass;
                                    break;
                                }
                                // const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];
                                // layer_volume += 1.0/m_inv;
                            }
                        }
                    }
                    // mark the nodes in the outside layer with a small negative distance
                    for (unsigned int i = 0; i < first_outside.size(); i++)
                    {
                        unsigned int i_node = first_outside[i];
                        mdistances[i_node] = -mHavg[i_node];
                    }
                }
            }
            mr_matrix_container.WriteScalarToDatabase(DISTANCE, mdistances, mr_model_part.Nodes());
        }
        void SetWallReductionCoefficients(double corner_coefficient, double edge_coefficient)
        {
            mcorner_coefficient = corner_coefficient;
            medge_coefficient = edge_coefficient;
        }
        void ContinuousVolumeCorrection(double expected_volume, double measured_volume)
        {
            double volume_error = expected_volume - measured_volume;
            if (volume_error == 0.0)
                return;
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
                    if (dist > 0.0) // node is outside domain
                    {
                        for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                        {
                            unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
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
                        const double nodal_mass = 1.0 / mr_matrix_container.GetInvertedMass()[i_node];
                        first_outside.push_back(i_node);
                        layer_volume += nodal_mass;
                    }
                }
                if (layer_volume == 0.00)
                    return;
                double ratio = volume_error / layer_volume;
                if (ratio > 1.0)
                    ratio = 1.0;
                if (ratio < 0.1) // NO correction for less than 10% error
                    return;
                double average_layer_h = 0.0;
                for (unsigned int i = 0; i < first_outside.size(); i++)
                {
                    unsigned int i_node = first_outside[i];
                    average_layer_h += mHavg[i_node];
                }
                average_layer_h /= static_cast<double>(first_outside.size());
                for (int i_node = 0; i_node < n_nodes; i_node++)
                    mdistances[i_node] -= average_layer_h * ratio;
            }
            mr_matrix_container.WriteScalarToDatabase(DISTANCE, mdistances, mr_model_part.Nodes());

            return;
        }

        void CalculatePorousResistanceLaw(unsigned int res_law)
        {
            if (res_law == 1)
            {
                /* if the chosen resistance law is ERGUN calculate Ergun A and B*/
                for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                     inode != mr_model_part.NodesEnd();
                     inode++)
                {
                    const double eps = inode->FastGetSolutionStepValue(POROSITY);
                    const double d = inode->FastGetSolutionStepValue(DIAMETER);
                    double &a = inode->FastGetSolutionStepValue(LIN_DARCY_COEF);
                    double &b = inode->FastGetSolutionStepValue(NONLIN_DARCY_COEF);
                    if (eps < 1.0)
                    {
                        double k_inv = 150.0 * (1.0 - eps) * (1.0 - eps) / (eps * eps * eps * d * d);
                        a = mViscosity * k_inv;
                        b = (1.75 / eps) * sqrt(k_inv / (150.0 * eps));
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
                    const double eps = inode->FastGetSolutionStepValue(POROSITY);   /*reading from kratos database*/
                    double &a = inode->FastGetSolutionStepValue(LIN_DARCY_COEF);    /*changing kratos database*/
                    double &b = inode->FastGetSolutionStepValue(NONLIN_DARCY_COEF); /*changing kratos database*/
                    if (eps == 1.0)
                    {
                        a = 0;
                        b = 0;
                    }
                }
            }
            mr_matrix_container.FillScalarFromDatabase(LIN_DARCY_COEF, mA, mr_model_part.Nodes());    /*filling edgebased database reading from kratos database*/
            mr_matrix_container.FillScalarFromDatabase(NONLIN_DARCY_COEF, mB, mr_model_part.Nodes()); /*filling edgebased database reading from kratos database*/
        }

    private:
        double mcorner_coefficient;
        double medge_coefficient;
        double mmax_dt;
        MatrixContainer &mr_matrix_container;
        ModelPart &mr_model_part;
        int mnumsubsteps;
        bool muse_mass_correction;
        // parameters controlling the wall law
        bool mWallLawIsActive;
        double mY_wall;
        // parameters for controlling the usage of the delta time in the stabilization
        double mstabdt_pressure_factor;
        double mstabdt_convection_factor;
        double mtau2_factor;
        bool massume_constant_dp;
        // nodal values
        CalcVectorType mBodyForce;
        ValuesVectorType mViscosity;
        // velocity vector U at time steps n and n+1
        CalcVectorType mWork, mvel_n, mvel_n1, mx, macc;
        // pressure vector p at time steps n and n+1
        ValuesVectorType mPn, mPn1;
        // coefficients
        ValuesVectorType mdistances;
        // minimum length of the edges surrounding edges surrounding each nodal point
        ValuesVectorType mHmin;
        ValuesVectorType mHavg;
        CalcVectorType mEdgeDimensions;
        // area normal
        CalcVectorType mSlipNormal;
        CalcVectorType mInOutNormal;
        // projection terms
        CalcVectorType mPi, mXi;
        // flag for first time step
        bool mFirstStep;
        // flag to differentiate interior and boundary nodes
        ValuesVectorType mNodalFlag;
        ValuesVectorType mWallReductionFactor;
        // lists of nodes with different types of boundary conditions
        IndicesVectorType mSlipBoundaryList, mPressureOutletList, mFixedVelocities, mInOutBoundaryList, mDistanceBoundaryList;
        ValuesVectorType mDistanceValuesList;
        CalcVectorType mFixedVelocitiesValues;
        //	ValuesVectorType mPressureOutlet;
        // intrinsic time step size
        ValuesVectorType mTauPressure;
        ValuesVectorType mTauConvection;
        ValuesVectorType mTau2;
        ValuesVectorType mdiv_error;
        boost::numeric::ublas::vector<bool> mis_slip;
        boost::numeric::ublas::vector<int> mis_visited;
        // variables for resolving pressure equation
        // laplacian matrix
        TSystemMatrixType mL;
        // constant variables
        double mRho;
        // variables for convection
        ValuesVectorType mphi_n;
        ValuesVectorType mphi_n1;
        CalcVectorType mPiConvection;
        ValuesVectorType mBeta;
        // variables for edge BCs
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
        // functions to calculate area normals for boundary conditions
        void CalculateNormal2D(ModelPart::ConditionsContainerType::iterator cond_it, array_1d<double, 3> &area_normal)
        {
            Geometry<Node> &face_geometry = (cond_it)->GetGeometry();
            area_normal[0] = face_geometry[1].Y() - face_geometry[0].Y();
            area_normal[1] = -(face_geometry[1].X() - face_geometry[0].X());
            area_normal[2] = 0.00;
            noalias((cond_it)->GetValue(NORMAL)) = area_normal;
        }
        void CalculateNormal3D(ModelPart::ConditionsContainerType::iterator cond_it, array_1d<double, 3> &area_normal, array_1d<double, 3> &v1, array_1d<double, 3> &v2)
        {
            Geometry<Node> &face_geometry = (cond_it)->GetGeometry();

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
        // function to calculate minimum length of surrounding edges
        void CalculateEdgeLengths(ModelPart::NodesContainerType &rNodes)
        {
            KRATOS_TRY
            // get number of nodes
            unsigned int n_nodes = rNodes.size();
            // reserve memory for storage of nodal coordinates
            std::vector<array_1d<double, TDim>> position;
            position.resize(n_nodes);
            // get position of all nodes
            for (typename ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
            {
                // get the global index of the node
                unsigned int i_node = static_cast<unsigned int>(node_it->FastGetSolutionStepValue(AUX_INDEX));
                // save its coordinates locally
                noalias(position[i_node]) = node_it->Coordinates();
            }
            ValuesVectorType &aaa = mr_matrix_container.GetHmin();
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                mHmin[i_node] = aaa[i_node];

                KRATOS_ERROR_IF(aaa[i_node] == 0.0) << "found a 0 hmin on node " << i_node << std::endl;
            }
            // take unstructured meshes into account
            if constexpr (TDim == 2)
            {
                for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
                {
                    double &h_i = mHavg[i_node];
                    double &m_i = mr_matrix_container.GetLumpedMass()[i_node];
                    h_i = sqrt(2.0 * m_i);
                }
            }
            else if constexpr (TDim == 3)
            {
                for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
                {
                    double &h_i = mHavg[i_node];
                    double &m_i = mr_matrix_container.GetLumpedMass()[i_node];
                    h_i = pow(6.0 * m_i, 1.0 / 3.0);
                }
            }
            // compute edge coordinates
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                array_1d<double, TDim> &pos_i = position[i_node];
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    array_1d<double, TDim> &pos_j = position[j_neighbour];
                    array_1d<double, TDim> &l_k = mEdgeDimensions[csr_index];
                    for (unsigned int comp = 0; comp < TDim; comp++)
                        l_k[comp] = pos_i[comp] - pos_j[comp];
                }
            }
            KRATOS_CATCH("")
        }
        //*********************************************************************
        // function to calculate right-hand side of fractional momentum equation
        void CalculateRHS_convection(
            const ValuesVectorType &mphi,
            const CalcVectorType &convective_velocity,
            ValuesVectorType &rhs,
            ValuesVectorType &active_nodes)
        {
            KRATOS_TRY
            int n_nodes = mphi.size();

            // calculating the RHS
            double stab_low;
            double stab_high;
            array_1d<double, TDim> a_i;
            array_1d<double, TDim> a_j;

            #pragma omp parallel for private(stab_low, stab_high, a_i, a_j)
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                double &rhs_i = rhs[i_node];
                const double &h_i = mHavg[i_node];
                const double &phi_i = mphi[i_node];
                noalias(a_i) = convective_velocity[i_node];
                a_i /= mEps[i_node];
                const array_1d<double, TDim> &proj_i = mPiConvection[i_node];
                double pi_i = proj_i[0] * a_i[0];
                for (unsigned int l_comp = 1; l_comp < TDim; l_comp++)
                    pi_i += proj_i[l_comp] * a_i[l_comp];
                rhs_i = 0.0;
                if (active_nodes[i_node] != 0.0)
                {
                    const double &beta = mBeta[i_node];
                    double norm_a = a_i[0] * a_i[0];
                    for (unsigned int l_comp = 1; l_comp < TDim; l_comp++)
                        norm_a += a_i[l_comp] * a_i[l_comp];
                    norm_a = sqrt(norm_a);
                    // loop to all the edges surrounding node I
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        if (active_nodes[j_neighbour] != 0.0)
                        {
                            // double& rhs_j = rhs[j_neighbour];
                            const double &phi_j = mphi[j_neighbour];
                            noalias(a_j) = convective_velocity[j_neighbour];
                            a_j /= mEps[j_neighbour];
                            const array_1d<double, TDim> &proj_j = mPiConvection[j_neighbour];
                            double pi_j = proj_j[0] * a_i[0];
                            for (unsigned int l_comp = 1; l_comp < TDim; l_comp++)
                                pi_j += proj_j[l_comp] * a_i[l_comp];
                            CSR_Tuple &edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
                            // convection operator
                            edge_ij.Sub_ConvectiveContribution(rhs_i, a_i, phi_i, a_j, phi_j);
                            // calculate stabilization part
                            edge_ij.CalculateConvectionStabilization_LOW(stab_low, a_i, phi_i, a_j, phi_j);
                            double edge_tau = mTauConvection[i_node];
                            edge_ij.CalculateConvectionStabilization_HIGH(stab_high, a_i, pi_i, a_j, pi_j);
                            edge_ij.Sub_StabContribution(rhs_i, edge_tau, 1.0, stab_low, stab_high);
                            double coeff = 0.5 * mshock_coeff; //=0.7*0.5;
                            double laplacian_ij = 0.0;
                            edge_ij.CalculateScalarLaplacian(laplacian_ij);
                            double capturing = laplacian_ij * (phi_j - phi_i);
                            double aaa = 0.0;
                            for (unsigned int k_comp = 0; k_comp < TDim; k_comp++)
                                for (unsigned int m_comp = 0; m_comp < TDim; m_comp++)
                                    aaa += a_i[k_comp] * a_i[m_comp] * edge_ij.LaplacianIJ(k_comp, m_comp);
                            if (norm_a > 1e-10)
                            {
                                aaa /= (norm_a * norm_a);
                                double capturing2 = aaa * (phi_j - phi_i);
                                if (fabs(capturing) > fabs(capturing2))
                                    rhs_i -= coeff * (capturing - capturing2) * beta * norm_a * h_i;
                            }
                        }
                    }
                }
                // KRATOS_WATCH(rhs_i);
            }

            KRATOS_CATCH("")
        }
        //**************************************
        void CornerDectectionHelper(Geometry<Node> &face_geometry,
                                    const array_1d<double, 3> &face_normal,
                                    const double An,
                                    const GlobalPointersVector<Condition> &neighb,
                                    const unsigned int i1,
                                    const unsigned int i2,
                                    const unsigned int neighb_index,
                                    std::vector<unsigned int> &edge_nodes,
                                    CalcVectorType &cornern_list)
        {
            double acceptable_angle = 45.0 / 180.0 * 3.1; // angles of less than 45 deg will be accepted
            double acceptable_cos = cos(acceptable_angle);
            if (face_geometry[i1].Id() < face_geometry[i2].Id()) // we do this to add the face ones
            {
                const array_1d<double, 3> &neighb_normal = neighb[neighb_index].GetValue(NORMAL);
                double neighb_An = norm_2(neighb_normal);
                double cos_normal = 1.0 / (An * neighb_An) * inner_prod(face_normal, neighb_normal);
                // if the angle is too big between the two normals then the edge in the middle is a corner
                if (cos_normal < acceptable_cos)
                {
                    array_1d<double, 3> edge = face_geometry[i2].Coordinates() - face_geometry[i1].Coordinates();
                    double temp = norm_2(edge);
                    edge /= temp;
                    int index1 = face_geometry[i1].FastGetSolutionStepValue(AUX_INDEX);
                    int index2 = face_geometry[i2].FastGetSolutionStepValue(AUX_INDEX);
                    edge_nodes[index1] += 1;
                    edge_nodes[index2] += 1;
                    //                double sign1 = inner_prod (cornern_list[index1], edge);
                    double sign1 = 0.0;
                    for (unsigned int i = 0; i < edge.size(); i++)
                    {
                        sign1 += cornern_list[index1][i] * edge[i];
                    }

                    if (sign1 >= 0)
                    {
                        for (unsigned int i = 0; i < edge.size(); i++)
                            cornern_list[index1][i] += edge[i];
                    }
                    else
                    {
                        for (unsigned int i = 0; i < edge.size(); i++)
                            cornern_list[index1][i] -= edge[i];
                    }

                    double sign2 = inner_prod(cornern_list[index2], edge);
                    if (sign2 >= 0)
                    {
                        for (unsigned int i = 0; i < edge.size(); i++)
                            cornern_list[index2][i] += edge[i];
                    }
                    else
                    {
                        for (unsigned int i = 0; i < edge.size(); i++)
                            cornern_list[index2][i] -= edge[i];
                    }
                }
            }
        }
        // function to calculate the area normals
        void DetectEdges3D(ModelPart::ConditionsContainerType &rConditions)
        {
            KRATOS_TRY
            // calculate area normals face-by-face
            array_1d<double, 3> area_normal;
            //(re)initialize normals
            unsigned int n_nodes = mNodalFlag.size();
            std::vector<unsigned int> temp_edge_nodes(n_nodes);
            CalcVectorType temp_cornern_list(n_nodes);
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                temp_edge_nodes[i_node] = 0.0;
                noalias(temp_cornern_list[i_node]) = ZeroVector(TDim);
            }
            // loop over all faces
            for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
            {
                // get geometry data of the face
                Geometry<Node> &face_geometry = cond_it->GetGeometry();
                // reference for area normal of the face
                const array_1d<double, 3> &face_normal = cond_it->GetValue(NORMAL);
                double An = norm_2(face_normal);
                unsigned int current_id = cond_it->Id();
                // slip condition
                if (cond_it->Is(SLIP)) // this is a slip face --> now look for its neighbours
                {
                    const GlobalPointersVector<Condition> &neighb = cond_it->GetValue(NEIGHBOUR_CONDITIONS);
                    // check for neighbour zero
                    if (neighb[0].Id() != current_id) // check if the neighbour exists
                        CornerDectectionHelper(face_geometry, face_normal, An, neighb, 1, 2, 0, temp_edge_nodes, temp_cornern_list);
                    // check for neighbour one
                    if (neighb[1].Id() != current_id) // check if the neighbour exists
                        CornerDectectionHelper(face_geometry, face_normal, An, neighb, 2, 0, 1, temp_edge_nodes, temp_cornern_list);
                    // check for neighbour two
                    if (neighb[2].Id() != current_id) // check if the neighbour exists
                        CornerDectectionHelper(face_geometry, face_normal, An, neighb, 0, 1, 2, temp_edge_nodes, temp_cornern_list);
                }
            }

            // fill the list of edge_nodes
            std::vector<unsigned int> tempmedge_nodes;
            std::vector<array_1d<double, TDim>> tempmedge_nodes_direction;
            std::vector<unsigned int> tempmcorner_nodes;
            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                if (temp_edge_nodes[i_node] == 2) // node is a edge_node
                {
                    tempmedge_nodes.push_back(i_node);
                    array_1d<double, TDim> &node_edge = temp_cornern_list[i_node];
                    node_edge /= norm_2(node_edge);
                    tempmedge_nodes_direction.push_back(node_edge);
                }
                else if (temp_edge_nodes[i_node] > 2)
                    tempmcorner_nodes.push_back(i_node);
            }
            medge_nodes.resize(tempmedge_nodes.size(), false);
            medge_nodes_direction.resize(tempmedge_nodes_direction.size(), false);
            mcorner_nodes.resize(tempmcorner_nodes.size(), false);

            IndexPartition<unsigned int>(tempmedge_nodes.size()).for_each([&](unsigned int i){
                medge_nodes[i] = tempmedge_nodes[i];
                medge_nodes_direction[i] = tempmedge_nodes_direction[i];
            });

            IndexPartition<unsigned int>(tempmcorner_nodes.size()).for_each([&](unsigned int i){
                mcorner_nodes[i] = tempmcorner_nodes[i];
            });

            for (unsigned int i = 0; i < mcorner_nodes.size(); i++)
            {
                KRATOS_WATCH(mcorner_nodes[i]);
            }
            KRATOS_CATCH("")
        }

        double ComputePorosityCoefficient(const double &vel_norm, const double &eps, const double &a, const double &b)
        {
            double linear;
            double non_linear;
            linear = eps * a;
            non_linear = eps * b * vel_norm;
            return linear + non_linear;
        }
        void LaplacianSmooth(ValuesVectorType &to_be_smoothed, ValuesVectorType &aux)
        {
            ModelPart::NodesContainerType &rNodes = mr_model_part.Nodes();
            int n_nodes = rNodes.size();

            IndexPartition<unsigned int>(n_nodes).for_each([&](unsigned int i_node){
                double dist = mdistances[i_node];
                double correction = 0.0;
                const double &origin_i = to_be_smoothed[i_node];
                if (dist <= 0.0) // node is inside domain ---- if outside do nothing
                {
                    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                    {
                        unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                        const double &origin_j = to_be_smoothed[j_neighbour];
                        CSR_Tuple &edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
                        double l_ikjk;
                        edge_ij.CalculateScalarLaplacian(l_ikjk);
                        correction += l_ikjk * (origin_j - origin_i);
                    }
                }
                aux[i_node] = origin_i - correction;
            });

            IndexPartition<unsigned int>(n_nodes).for_each([&](unsigned int i_node){
                to_be_smoothed[i_node] = aux[i_node];
            });
        }

        void ComputeWallResistance(
            const CalcVectorType &vel,
            ValuesVectorType &diag_stiffness
        )
        {

            double ym = mY_wall;


            // slip condition
            int slip_size = mSlipBoundaryList.size();

            #pragma omp parallel for firstprivate(slip_size, ym)
            for (int i_slip = 0; i_slip < slip_size; i_slip++)
            {
                unsigned int i_node = mSlipBoundaryList[i_slip];
                double dist = mdistances[i_node];
                if (dist <= 0.0)
                {
                    KRATOS_ERROR_IF(mViscosity[i_node] == 0) << "it is not possible to use the wall law with 0 viscosity" << std::endl;

                    double nu = mViscosity[i_node];

                    const array_1d<double, TDim> &U_i = vel[i_node];
                    const array_1d<double, TDim> &an_i = mSlipNormal[i_node];

                    // compute the modulus of the velocity
                    double mod_vel = 0.0;
                    double area = 0.0;
                    for (unsigned int comp = 0; comp < TDim; comp++)
                    {
                        mod_vel += U_i[comp] * U_i[comp];
                        area += an_i[comp] * an_i[comp];
                    }
                    mod_vel = sqrt(mod_vel);
                    area = sqrt(area);

                    // the 0.1 is such that the dissipation is as for the linear case for a velocity of 10m/s
                    diag_stiffness[i_node] = area * nu * mod_vel / (ym)*mWallReductionFactor[i_node];
                }
                else
                {
                    diag_stiffness[i_node] = 0.0;
                }
            }
        }

        void ApplySmagorinsky3D(double MolecularViscosity, double Cs)
        {
            KRATOS_TRY
            ModelPart::NodesContainerType &rNodes = mr_model_part.Nodes();
            // calculating the RHS
            array_1d<double, TDim> grad_vx;
            array_1d<double, TDim> grad_vy;
            array_1d<double, TDim> grad_vz;
            int n_nodes = rNodes.size();
            mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes);
            array_1d<double, TDim> stab_high;

            #pragma omp parallel for private(grad_vx, grad_vy, grad_vz)
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                // set to zero the gradients
                for (unsigned int comp = 0; comp < TDim; comp++)
                {
                    grad_vx[comp] = 0.0;
                    grad_vy[comp] = 0.0;
                    grad_vz[comp] = 0.0;
                }
                // compute node by node the gradients
                const array_1d<double, TDim> &U_i = mvel_n1[i_node];
                const double h = mHmin[i_node];
                const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    const array_1d<double, TDim> &U_j = mvel_n1[j_neighbour];
                    CSR_Tuple &edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
                    edge_ij.Add_grad_p(grad_vx, U_i[0], U_j[0]);
                    edge_ij.Add_grad_p(grad_vy, U_i[1], U_j[1]);
                    edge_ij.Add_grad_p(grad_vz, U_i[2], U_j[2]);
                }
                // finalize computation of the gradients
                // set to zero the gradients
                for (unsigned int comp = 0; comp < TDim; comp++)
                {
                    grad_vx[comp] *= m_inv;
                    grad_vy[comp] *= m_inv;
                    grad_vz[comp] *= m_inv;
                }
                // symmetrize and multiply by 2
                grad_vx[0] *= 2.0;
                grad_vy[1] *= 2.0;
                if constexpr (TDim > 2)
                    grad_vz[2] *= 2.0;
                grad_vx[1] += grad_vy[0];
                if constexpr (TDim > 2)
                    grad_vx[2] += grad_vz[0];
                if constexpr (TDim > 2)
                    grad_vy[2] += grad_vz[1];
                grad_vy[0] += grad_vx[1];
                grad_vz[0] += grad_vx[2];
                grad_vz[1] += grad_vy[2];

                // compute smagorinsky term
                double aux = 0.0;
                for (unsigned int comp = 0; comp < TDim; comp++)
                {
                    aux += grad_vx[comp] * grad_vx[comp];
                    aux += grad_vy[comp] * grad_vy[comp];
                    aux += grad_vz[comp] * grad_vz[comp];
                }
                aux *= 0.5;
                if (aux < 0.0)
                    aux = 0.0;
                double turbulent_viscosity = Cs * h * h * sqrt(aux);
                mViscosity[i_node] = turbulent_viscosity + MolecularViscosity;
            }
            mr_matrix_container.WriteScalarToDatabase(VISCOSITY, mViscosity, rNodes);
            KRATOS_CATCH("");
        }

        void Add_Effective_Inverse_Multiply(
            CalcVectorType &destination,
            const CalcVectorType &origin1,
            const double value,
            const ValuesVectorType &mass,
            const ValuesVectorType &diag_stiffness,
            const CalcVectorType &origin)
        {
            KRATOS_TRY
            int loop_size = destination.size();


            IndexPartition<unsigned int>(loop_size).for_each([&](unsigned int i_node){
                array_1d<double, TDim> &dest = destination[i_node];
                const double m = mass[i_node];
                const double d = diag_stiffness[i_node];
                const array_1d<double, TDim> &origin_vec1 = origin1[i_node];
                const array_1d<double, TDim> &origin_value = origin[i_node];

                for (unsigned int comp = 0; comp < TDim; comp++)
                    dest[comp] = value / (m + value * d) * (m / value * origin_vec1[comp] + origin_value[comp]);
            });

            KRATOS_CATCH("")
        }

        void ComputeConvectiveProjection(
            CalcVectorType &mPiConvection,
            const ValuesVectorType &mphi_n1,
            const ValuesVectorType &mEps,
            const CalcVectorType &mvel_n1)
        {
            int n_nodes = mPiConvection.size();
            // calculating the convective projection
            array_1d<double, TDim> a_i;
            array_1d<double, TDim> a_j;

            #pragma omp parallel for private(a_i, a_j)
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                array_1d<double, TDim> &pi_i = mPiConvection[i_node];
                // 		    setting to zero the projection
                for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                    pi_i[l_comp] = 0.0;

                const double &phi_i = mphi_n1[i_node];
                noalias(a_i) = mvel_n1[i_node];
                a_i /= mEps[i_node];
                // 			  loop to all the edges surrounding node I
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    noalias(a_j) = mvel_n1[j_neighbour];
                    a_j /= mEps[j_neighbour];
                    const double &phi_j = mphi_n1[j_neighbour];
                    CSR_Tuple &edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
                    edge_ij.Add_grad_p(pi_i, phi_i, phi_j);
                }
                // 			  apply inverted mass matrix
                const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];
                for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
                    pi_i[l_comp] *= m_inv;
            }
        }

        void ComputeLimitor(
            CalcVectorType &mPiConvection,
            const ValuesVectorType &mphi_n1,
            ValuesVectorType &mBeta,
            const CalcVectorType &mvel_n1,
            const CalcVectorType &mEdgeDimensions)
        {
            int n_nodes = mPiConvection.size();

            IndexPartition<unsigned int>(n_nodes).for_each([&](unsigned int i_node){
                const array_1d<double, TDim> &pi_i = mPiConvection[i_node];
                const double &p_i = mphi_n1[i_node];
                double &beta_i = mBeta[i_node];
                beta_i = 0.0;
                double n = 0.0;
                for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
                {
                    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
                    const double &p_j = mphi_n1[j_neighbour];
                    const array_1d<double, TDim> &l_k = mEdgeDimensions[csr_index];
                    const array_1d<double, TDim> &pi_j = mPiConvection[j_neighbour];
                    double proj = 0.0;
                    for (unsigned int comp = 0; comp < TDim; comp++)
                        proj += 0.5 * l_k[comp] * (pi_i[comp] + pi_j[comp]);
                    // 							proj += dir[comp]*pi_i[comp];
                    double numerator = fabs(fabs(p_j - p_i) - fabs(proj));
                    double denom = fabs(fabs(p_j - p_i) + 1e-6);
                    beta_i += numerator / denom;
                    n += 1.0;
                }
                beta_i /= n;
                if (beta_i > 1.0)
                    beta_i = 1.0;
            });
        }
    };
} // namespace Kratos
#undef SYMM_PRESS
#endif // KRATOS_EDGEBASED_LEVELSET_SUBSTEP_FLUID_SOLVER_H_INCLUDED defined
