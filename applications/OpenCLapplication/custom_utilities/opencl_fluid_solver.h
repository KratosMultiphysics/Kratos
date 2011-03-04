/*
==============================================================================
KratosOpenCLApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Farshid Mossaiby
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
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
//   Last Modified by:    $Author: mossaiby $
//   Date:                $Date: 2011-02-24 15:23:33 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_OPENCL_EDGEBASED_LEVELSET_FLUID_SOLVER_H_INCLUDED)
#define KRATOS_OPENCL_EDGEBASED_LEVELSET_FLUID_SOLVER_H_INCLUDED

// TODO: What are these?
//#define SPLIT_OSS
//#define SYMM_PRESS


// System includes
#include <string>
#include <iostream>
#include <algorithm>


// External includes


// Project includes
#include "includes/model_part.h"

// TODO: Remove unneeded ones
//#include "includes/node.h"
//#include "geometries/geometry.h"
//#include "utilities/geometry_utilities.h"
//#include "incompressible_fluid_application.h"


namespace Kratos
{
	//
	// OpenCLFluidSolver3D
	//
	// OpenCL based 3D fluid solver

	class OpenCLFluidSolver3D
	{
		public:

			//
			// Used types

			typedef cl_uint *IndicesVectorType;
			typedef cl_double3 *CalcVectorType;
			typedef cl_double *ValuesVectorType;
			//typedef std::vector <unsigned int> IndicesSTDVectorType;  // TODO: Is this OK?
			//typedef std::vector <cl_double3> CalcSTDVectorType;  // TODO: Is this OK?
			//typedef std::vector <double> ValuesSTDVectorType;  // TODO: Is this OK?

			//
			// OpenCLFluidSolver3D
			//
			// Constructor

			OpenCLFluidSolver3D(
				OpenCLMatrixContainer &matrix_container,
				ModelPart &model_part,
				const double viscosity,
				const double density,
				const Vector body_force,
				bool use_mass_correction,
				double edge_detection_angle,
				double stabdt_pressure_factor,
				double stabdt_convection_factor,
				double tau2_factor,
				bool assume_constant_dp
			):
				mrDeviceGroup(matrix_container.GetDeviceGroup()),
				mr_matrix_container(matrix_container),
				mr_model_part(model_part),
				mstabdt_pressure_factor(stabdt_pressure_factor),
				mstabdt_convection_factor(stabdt_convection_factor),
				medge_detection_angle(edge_detection_angle),
				mtau2_factor(tau2_factor),
				massume_constant_dp(assume_constant_dp)
			{
				mViscosity = viscosity;

				noalias(mBodyForce) = body_force;
				mRho = density;

				mdelta_t_avg = 1000.00;

				max_dt = 1.00;

				muse_mass_correction = use_mass_correction;

				mWallLawIsActive = false;
			}

			//
			// ~OpenCLFluidSolver3D
			//
			// Destructor

			~OpenCLFluidSolver3D()
			{
				// Nothing to do!
			}

			//
			// Initialize
			//
			// Initializes the solver

			void Initialize()
			{
				KRATOS_TRY

				// Get no. of nodes
				n_nodes = mr_model_part.Nodes().size();
				n_edges = mr_matrix_container.GetNumberEdges();

				// TODO: Check if these are all needed

				// Size data vectors
				AllocateArray(&mWork, n_nodes);
				AllocateArray(&mvel_n, n_nodes);
				AllocateArray(&mvel_n1, n_nodes);
				AllocateArray(&mPn, n_nodes);
				AllocateArray(&mPn1, n_nodes);
				AllocateArray(&mHmin, n_nodes);
				AllocateArray(&mHavg, n_nodes);
				AllocateArray(&mNodalFlag, n_nodes);

				AllocateArray(&mTauPressure, n_nodes);
				AllocateArray(&mTauConvection, n_nodes);
				AllocateArray(&mTau2, n_nodes);
				AllocateArray(&mPi, n_nodes);
				AllocateArray(&mXi, n_nodes);
				AllocateArray(&mx, n_nodes);

				AllocateArray(&mEdgeDimensions, n_edges);

				// Convection variables
				AllocateArray(&mBeta, n_nodes);
				AllocateArray(&mdiv_error, n_nodes);

				// RHS
				AllocateArray(&mrhs, n_nodes);

				// Allocating buffers on OpenCL device
				mbWork = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);
				mbvel_n = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);
				mbvel_n1 = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);

				mbPn = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);
				mbPn1 = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);
				mbHmin = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);
				mbHavg = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);

				mbNodalFlag = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);

				mbTauPressure = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);
				mbTauConvection = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);
				mbTau2 = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);

				mbPi = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);
				mbXi = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);
				mbx = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);

				mbEdgeDimensions = mrDeviceGroup.CreateBuffer(n_edges * sizeof(cl_double3), CL_MEM_READ_WRITE);

				mbBeta = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);
				mbdiv_error = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);

				// RHS
				mbrhs = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);

				mr_matrix_container.SetToZero(mbdiv_error);

				// Read velocity and pressure data from Kratos

				// TODO: It seems that we do not need this as Solve() does this at first step
				/*
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes());
				mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, mr_model_part.Nodes());
				mr_matrix_container.FillOldScalarFromDatabase(PRESSURE, mPn, mr_model_part.Nodes());
				mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mvel_n, mr_model_part.Nodes());
				*/
				mr_matrix_container.FillCoordinatesFromDatabase(mx, mr_model_part.Nodes(), mbx);

				// Set flag for first time step
				mFirstStep = true;

				// Loop to categorize boundary nodes
				for (
					ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
					inode != mr_model_part.NodesEnd();
					inode++)
				{
					int index = inode -> FastGetSolutionStepValue(AUX_INDEX);

					// Note that the variables can be either all fixed or no one fixed
					if (inode -> IsFixed(VELOCITY_X))
					{
						if (inode -> IsFixed(VELOCITY_Y) == false || inode -> IsFixed(VELOCITY_Z) == false)
						{
							std::cout << "Error found on the fixity of node " << inode -> Id() << std::endl;
							KRATOS_ERROR(std::logic_error, "Velocities can be either all fixed or none fixed", "")
						}

						// TODO: Is this OK?
						//mFixedVelocities.push_back(index);
						//mFixedVelocitiesValues.push_back(mvel_n1[index]);
						mFixedVelocitiesList[mFixedVelocitiesListLength] = index;
						mFixedVelocitiesValuesList[mFixedVelocitiesListLength] = mvel_n1[index];
						mFixedVelocitiesListLength++;
					}

					if (inode -> IsFixed(PRESSURE))
					{
						// TODO: Is this OK?
						//mPressureOutletList.push_back(index);
						mPressureOutletList[mPressureOutletListLength] = index;
						mPressureOutletListLength++;
					}
				}

				// Compute slip normals and fill SlipList
				CalculateNormals(mr_model_part.Conditions());
				mr_matrix_container.WriteVectorToDatabase(NORMAL, mSlipNormal, mr_model_part.Nodes(), mbSlipNormal);

				DetectEdges3D(mr_model_part.Conditions());

				// Determine number of edges and entries
				unsigned int n_nonzero_entries = 2 * n_edges + n_nodes;

				// Allocate memory for variables
				//mL.resize(n_nodes, n_nodes, n_nonzero_entries);  // TODO: Fix this!

				// Loop over all nodes
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					// Flag for considering diagonal matrix elements
					bool flag = 0;

					// Loop over all neighbours
					for (
						unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node];
						csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1];
						csr_index++)
					{
						// Get global index of neighbouring node j
						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];

						// Define matrix structure row by row (the order does matter!)
						if ((j_neighbour > i_node) && (flag == 0))
						{
							// Add diagonal/nodal contribution
							//mL.push_back(i_node, i_node, 0.00);  // TODO: Fix this!
							flag = 1;
						}

						// Add non-diagonal/edge contribution
						//mL.push_back(i_node, j_neighbour, 0.00);  // TODO: Fix this!
					}

					// If diagonal element is the last non-zero element of the row
					if (flag == 0)
						//mL.push_back(i_node, i_node, 0.0);  // TODO: Fix this!
						;
				}

				// Compute minimum length of the surrounding edges
				CalculateEdgeLengths(mr_model_part.Nodes());

				// Set the pressure projection to the body force value
				array_1d <double, 3> temp = mRho * mBodyForce;
				for (
					ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
					inode != mr_model_part.NodesEnd();
					inode++)
				{
					inode->FastGetSolutionStepValue(PRESS_PROJ) = temp;
				}

				KRATOS_CATCH("")
			}

			//
			// ComputeTimeStep
			//
			// Function to set adequate time step size

			double ComputeTimeStep(const double CFLNumber, const double MaxDt)
			{
				KRATOS_TRY

				// Save the maximum time step
				max_dt = MaxDt;

				// Local variable for time step size
				double delta_t = 1e10;

				mdelta_t_avg = 1e10;

				// Getting value of current velocity and of viscosity
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes(), mbvel_n1);

				// Loop over all nodes
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					const double havg_i = mHavg[i_node];
					const double hmin_i = mHmin[i_node];

					double vel_norm = Norm2_3(mvel_n1[i_node]);


					//use CFL condition to compute time step size
					double delta_t_i = CFLNumber / (2.00 * vel_norm / hmin_i + 4.00 * mViscosity / (hmin_i * hmin_i));
					double delta_t_i_avg = 1.00 / (2.00 * vel_norm / havg_i + 4.00 * mViscosity / (havg_i * havg_i));

					// Considering the most restrictive case of neighbor's velocities with similar direction but opposite sense.

					// Loop over all neighbours
					for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
					{
						// Get global index of neighbouring node j
						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];

						double v_diff_norm = 0.00;
						for (unsigned int l_comp = 0; l_comp < 3; l_comp++)
						{
							double temp = KRATOS_OCL_COMP(mvel_n1[i_node], l_comp) - KRATOS_OCL_COMP(mvel_n1[j_neighbour], l_comp);
							v_diff_norm += temp * temp;
						}

						v_diff_norm = sqrt(v_diff_norm);

						double delta_t_j = CFLNumber / (2.00 * v_diff_norm / hmin_i + 4.00 * mViscosity / (hmin_i * hmin_i));

						if (delta_t_j < delta_t_i)
							delta_t_i = delta_t_j;
					}

					// Choose the overall minimum of delta_t_i
					if (delta_t_i < delta_t)
						delta_t = delta_t_i;

					if (delta_t_i_avg < mdelta_t_avg)
					  mdelta_t_avg = delta_t_i_avg;

				}

				return delta_t;

				KRATOS_CATCH("")
			}

			//
			// UpdateFixedVelocityValues
			//
			//

            void UpdateFixedVelocityValues()
            {
                KRATOS_TRY

                // Read velocity and pressure data from Kratos
                ModelPart::NodesContainerType &rNodes = mr_model_part.Nodes();
                mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes, mbvel_n1);

                unsigned int fixed_size = mFixedVelocitiesListLength; // TODO: Is this OK? //mFixedVelocities.size();

                #pragma omp parallel for firstprivate(fixed_size)
                for (unsigned int i_velocity = 0; i_velocity < fixed_size; i_velocity++)
                {
                    unsigned int i_node = mFixedVelocitiesList[i_velocity];
                    mFixedVelocitiesValuesList[i_velocity] = mvel_n1[i_node];
                }

                // TODO: Should we update anything on GPU?

                KRATOS_CATCH("")
            }

            //
            // SolveStep1
            //
            // Function to solve fluid equations - fractional step 1: compute fractional momentum

            void SolveStep1()
            {
				KRATOS_TRY

				// Prerequisites

				// Variables for node based data handling
				ModelPart::NodesContainerType &rNodes = mr_model_part.Nodes();

				// Read velocity and pressure data from Kratos
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes, mbvel_n1);
				mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mvel_n, rNodes, mbvel_n);

				mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, rNodes, mbPn1);
				mr_matrix_container.FillOldScalarFromDatabase(PRESSURE, mPn, rNodes, mbPn);

				// Read time step size from Kratos
				ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
				double delta_t = CurrentProcessInfo[DELTA_TIME];

				// Compute intrinsic time
				double time_inv_avg = 1.00 / mdelta_t_avg;

				double stabdt_pressure_factor = mstabdt_pressure_factor;
				double stabdt_convection_factor = mstabdt_convection_factor;
				double tau2_factor = mtau2_factor;

				// TODO: Call kernels here!

				KRATOS_CATCH("")
            }

            //
            // CalculateNormals
            //
            // Function to calculate the area normals

			void CalculateNormals(ModelPart::ConditionsContainerType &rConditions)
			{
				KRATOS_TRY

				// Calculate area normals face-by-face
				array_1d <double, 3> area_normal;

				// Helper vectors for cross product
				array_1d <double, 3> v1, v2;

				for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
				{
					CalculateNormal3D(cond_it, area_normal, v1, v2);
				}

				// (Re)initialize normals
				//unsigned int n_nodes = mNodalFlag.size();  // TODO: Is this OK?
				//mSlipNormal.resize(n_nodes);  // TODO: Fix this!
				std::vector<bool> is_slip(n_nodes);

				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					KRATOS_OCL_ZERO_VECTOR(mSlipNormal[i_node]);  // TODO: Fix this!
					is_slip[i_node] = false;
				}

				// Loop over all faces
				const double node_factor = 1.00 / 3.00;
				for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
				{
					// Get geometry data of the face
					Geometry <Node <3> > &face_geometry = cond_it -> GetGeometry();

					// Reference for area normal of the face
					array_1d <double, 3> &face_normal = cond_it -> GetValue(NORMAL);

					// Slip condition
					if (cond_it -> GetValue(IS_STRUCTURE) == true)
					{
						for (unsigned int if_node = 0; if_node < 3; if_node++)
						{
							unsigned int i_node = static_cast <unsigned int> (face_geometry[if_node].FastGetSolutionStepValue(AUX_INDEX));
							is_slip[i_node] = true;
							for (unsigned int comp = 0; comp < 3; comp++)
							{
								KRATOS_OCL_COMP(mSlipNormal[i_node], comp) += node_factor * face_normal[comp];
							}
						}
					}

				}

				// Fill the list of slip nodes
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					if (is_slip[i_node] == true)
					{
						// TODO: Is this OK?
						//mSlipBoundaryList.push_back(i_node);
						mSlipBoundaryList[mSlipBoundaryListLength] = i_node;
						mSlipBoundaryListLength++;
					}
				}

				KRATOS_CATCH("")
			}

			//
			// Clear
			//
			// Function to free dynamic memory

			void Clear()
			{
				KRATOS_TRY

				// TODO: Fix these! Copy from Initialize

				/*mWork.clear();
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

				mdiv_error.clear();*/

				KRATOS_CATCH("")
			}

			//
			// ActivateWallResistance
			//

			void ActivateWallResistance(double Ywall)
			{
				mWallLawIsActive = true;
				mY_wall = Ywall;
			}

		private:

			// OpenCL stuff
			OpenCL::DeviceGroup &mrDeviceGroup;
			cl_uint mbWork, mbvel_n, mbvel_n1, mbPn, mbPn1, mbHmin, mbHavg, mbNodalFlag, mbTauPressure, mbTauConvection, mbTau2, mbPi, mbXi, mbx, mbEdgeDimensions, mbBeta, mbdiv_error, mbSlipNormal, mbrhs;  // TODO: Fix this!

			// Matrix container
			OpenCLMatrixContainer &mr_matrix_container;

			// Associated model part
			ModelPart &mr_model_part;

			// No. of nodes
			unsigned int n_nodes;

			// No. of edges
			unsigned int n_edges;

			bool muse_mass_correction;

			// Parameters controlling the wall law
			bool mWallLawIsActive;
			bool mY_wall;

			// Parameters for controlling the usage of the delta time in the stabilization
			double mstabdt_pressure_factor;
			double mstabdt_convection_factor;
			double medge_detection_angle;
			double mtau2_factor;
			bool massume_constant_dp;

			// Nodal values

			// Velocity vector U at time steps n and n+1
			CalcVectorType mWork, mvel_n, mvel_n1, mx;

			// Pressure vector p at time steps n and n+1
			ValuesVectorType mPn, mPn1;

			// Minimum length of the edges surrounding edges surrounding each nodal point
			ValuesVectorType mHmin;
			ValuesVectorType mHavg;
			CalcVectorType mEdgeDimensions;

			// Area normal
			CalcVectorType mSlipNormal;

			// Projection terms
			CalcVectorType mPi, mXi;

			// Flag for first time step
			bool mFirstStep;

			// Flag to differentiate interior and boundary nodes
			ValuesVectorType mNodalFlag;

			// Lists of nodes with different types of boundary conditions
			IndicesVectorType mSlipBoundaryList, mPressureOutletList, mFixedVelocitiesList;  // TODO: These are all used on GPU, try to avoid using std::vector
			CalcVectorType mFixedVelocitiesValuesList;

			unsigned int mSlipBoundaryListLength, mPressureOutletListLength, mFixedVelocitiesListLength;

			// Intrinsic time step size
			ValuesVectorType mTauPressure;
			ValuesVectorType mTauConvection;
			ValuesVectorType mTau2;

			ValuesVectorType mdiv_error;

			// Storage of nodal values in local variables
			CalcVectorType mrhs;

			// Variables for resolving pressure equation

			// Laplacian matrix
			//TSystemMatrixType mL;  // TODO: Fix this! Should we use ViennaCL types?

			// Constant variables
			double mRho;
			double mViscosity;
			array_1d <double, 3> mBodyForce;

			// Variables for convection
			ValuesVectorType mBeta;

			// Variables for edge BCs
			IndicesVectorType medge_nodesList;  // TODO: Fix this!
			CalcVectorType medge_nodes_directionList;  // TODO: Is this OK?
			IndicesVectorType mcorner_nodesList;  // TODO: Fix this!

			unsigned int medge_nodesListLength, medge_nodes_directionListLength, mcorner_nodesListLength;

			double mdelta_t_avg;
			double max_dt;

			//
			// CalculateNormal3D
			//
			// Function to calculate area normals for boundary conditions

			void CalculateNormal3D(ModelPart::ConditionsContainerType::iterator cond_it, array_1d <double, 3> &area_normal, array_1d <double, 3> &v1, array_1d <double, 3> &v2)
			{
				Geometry <Node <3> > &face_geometry = (cond_it) -> GetGeometry();

				v1[0] = face_geometry[1].X() - face_geometry[0].X();
				v1[1] = face_geometry[1].Y() - face_geometry[0].Y();
				v1[2] = face_geometry[1].Z() - face_geometry[0].Z();

				v2[0] = face_geometry[2].X() - face_geometry[0].X();
				v2[1] = face_geometry[2].Y() - face_geometry[0].Y();
				v2[2] = face_geometry[2].Z() - face_geometry[0].Z();

				MathUtils<double>::CrossProduct(area_normal, v1, v2);
				area_normal *= -0.5;

				noalias((cond_it) -> GetValue(NORMAL)) = area_normal;
			}

			//
			// CalculateEdgeLengths
			//
			// Function to calculate minimum length of surrounding edges

			void CalculateEdgeLengths(ModelPart::NodesContainerType &rNodes)
			{
				KRATOS_TRY

				// Get number of nodes
				unsigned int n_nodes = rNodes.size();

				// Reserve memory for storage of nodal coordinates
				std::vector <array_1d <double, 3> > position;  // TODO: Is this OK?
				position.resize(n_nodes);

				// Get position of all nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of the node
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Save its coordinates locally
					//noalias(position[i_node]) = node_it -> Coordinates();
					position[i_node] = node_it -> Coordinates();  // TODO: Is this OK?

					// Initialize minimum edge length with relatively big values
					mHmin[i_node] = 1e10;
				}

				ValuesVectorType TempHmin = mr_matrix_container.GetHmin();
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					mHmin[i_node] = TempHmin[i_node];
				}

				// Take unstructured meshes into account
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					double &h_i = mHavg[i_node];
					double &m_i = mr_matrix_container.GetLumpedMass()[i_node];

					h_i = pow(6.0 * m_i, 1.0 / 3.0);
				}

				// Compute edge coordinates
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					array_1d <double, 3> &pos_i = position[i_node];

					for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
					{
						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
						array_1d <double, 3> &pos_j = position[j_neighbour];

						KRATOS_OCL_COMP_0(mEdgeDimensions[csr_index]) = pos_i[0] - pos_j[0];
						KRATOS_OCL_COMP_1(mEdgeDimensions[csr_index]) = pos_i[1] - pos_j[1];
						KRATOS_OCL_COMP_2(mEdgeDimensions[csr_index]) = pos_i[2] - pos_j[2];
						KRATOS_OCL_COMP_3(mEdgeDimensions[csr_index]) = 0.00;
					}
				}

				KRATOS_CATCH("")
			}

			//
			// CornerDectectionHelper
			//

			void CornerDectectionHelper(Geometry <Node<3> > &face_geometry,
				const array_1d <double, 3> &face_normal,
				const double An,
				const WeakPointerVector <Condition> &neighb,
				const unsigned int i1,
				const unsigned int i2,
				const unsigned int neighb_index,
				std::vector <unsigned int> &edge_nodes,
				std::vector <array_1d <double, 3> > &cornern_list)
			{
				double acceptable_angle = 45.00 / 180.00 * 3.1;  // Angles of less than 45 deg will be accepted
				double acceptable_cos = cos(acceptable_angle);

				if(face_geometry[i1].Id() < face_geometry[i2].Id()) // We do this to add the face ones
				{
					const array_1d <double, 3> &neighb_normal = neighb[neighb_index].GetValue(NORMAL);
					double neighb_An = norm_2(neighb_normal);

					double cos_normal = 1.00 / (An * neighb_An) * inner_prod(face_normal, neighb_normal);

					// If the angle is too big between the two normals then the edge in the middle is a corner
					if (cos_normal < acceptable_cos)
					{
						array_1d <double, 3> edge = face_geometry[i2].Coordinates() - face_geometry[i1].Coordinates();
						double temp = norm_2(edge);
						edge /= temp;

						int index1 = face_geometry[i1].FastGetSolutionStepValue(AUX_INDEX);
						int index2 = face_geometry[i2].FastGetSolutionStepValue(AUX_INDEX);

						edge_nodes[index1]++;
						edge_nodes[index2]++;

						double sign1 = inner_prod(cornern_list[index1], edge);

						if (sign1 >= 0)
							cornern_list[index1] += edge;
						else
							cornern_list[index1] -= edge;

						double sign2 = inner_prod(cornern_list[index2],edge);

						if (sign2 >= 0)
						  cornern_list[index2] += edge;
						else
						  cornern_list[index2] -= edge;
					}
				}
			}

			//
			// DetectEdges3D
			//
			// Function to calculate the area normals


			void DetectEdges3D(ModelPart::ConditionsContainerType& rConditions)
			{
				KRATOS_TRY

				// Calculate area normals face-by-face
				array_1d <double, 3> area_normal;

				// (Re)initialize normals
				//unsigned int n_nodes = mNodalFlag.size();  // TODO: Is this OK?
				std::vector <unsigned int> temp_edge_nodes(n_nodes);
				std::vector <array_1d <double, 3> > temp_cornern_list(n_nodes);

				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					temp_edge_nodes[i_node] = 0.00;
					noalias(temp_cornern_list[i_node]) = ZeroVector(3);
				}

				// Loop over all faces
				for (ModelPart::ConditionsContainerType::iterator cond_it = rConditions.begin(); cond_it != rConditions.end(); cond_it++)
				{
					// Get geometry data of the face
					Geometry <Node <3> > &face_geometry = cond_it -> GetGeometry();

					// Reference for area normal of the face
					const array_1d <double, 3> &face_normal = cond_it -> GetValue(NORMAL);
					double An = norm_2(face_normal);

					unsigned int current_id = cond_it -> Id();

					// Slip condition
					if (cond_it -> GetValue(IS_STRUCTURE) == 1.00) // This is a slip face, now look for its neighbours
					{
						const WeakPointerVector <Condition> &neighb = cond_it -> GetValue(NEIGHBOUR_CONDITIONS);

						// Check for neighbour zero
						if (neighb[0].Id() != current_id) // Check if the neighbour exists
						CornerDectectionHelper(face_geometry, face_normal, An, neighb, 1, 2, 0, temp_edge_nodes, temp_cornern_list);

						// Check for neighbour one
						if (neighb[1].Id() != current_id) // Check if the neighbour exists
						CornerDectectionHelper(face_geometry, face_normal, An, neighb, 2, 0, 1, temp_edge_nodes, temp_cornern_list);

						// Check for neighbour two
						if (neighb[2].Id() != current_id) // Check if the neighbour exists
						CornerDectectionHelper(face_geometry, face_normal, An, neighb, 0, 1, 2, temp_edge_nodes, temp_cornern_list);
					}
				}

				// Fill the list of edge_nodes
				// TODO: These lists should be allocated somewhere, Initialize?
				//medge_nodes.resize(0);
				//medge_nodes_direction.resize(0);
				//mcorner_nodes.resize(0);

				medge_nodesListLength = 0;
				medge_nodes_directionListLength = 0;
				mcorner_nodesListLength = 0;

				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					if (temp_edge_nodes[i_node] == 2) // Node is a edge_node
					{
						// TODO: Is this OK?
						//medge_nodes.push_back(i_node);
						medge_nodesList[medge_nodesListLength] = i_node;
						medge_nodesListLength++;

						array_1d <double, 3> &node_edge = temp_cornern_list[i_node];

						node_edge /= norm_2(node_edge);

						// TODO: Is this OK?
						//medge_nodes_direction.push_back(node_edge);
						KRATOS_OCL_COMP_0(medge_nodes_directionList[medge_nodes_directionListLength]) = node_edge[0];
						KRATOS_OCL_COMP_1(medge_nodes_directionList[medge_nodes_directionListLength]) = node_edge[1];
						KRATOS_OCL_COMP_2(medge_nodes_directionList[medge_nodes_directionListLength]) = node_edge[2];
						KRATOS_OCL_COMP_3(medge_nodes_directionList[medge_nodes_directionListLength]) = 0.00;
						medge_nodes_directionListLength++;
					}
					else if (temp_edge_nodes[i_node] > 2)
					{
						// TODO: Is this OK?
						// mcorner_nodes.push_back(i_node);
						mcorner_nodesList[mcorner_nodesListLength] = i_node;
						mcorner_nodesListLength++;
					}
				}

				// TODO: Fix this when needed
				//for (unsigned int i = 0; i < mcorner_nodesListLength; i++)
				//{
				//	KRATOS_WATCH(mcorner_nodes[i]);
				//}

				KRATOS_CATCH("")
			}
	};

}  // Namespace Kratos

#endif  // KRATOS_OPENCL_EDGEBASED_LEVELSET_FLUID_SOLVER_H_INCLUDED
