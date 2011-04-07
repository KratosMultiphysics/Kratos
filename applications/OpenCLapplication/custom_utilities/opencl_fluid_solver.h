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


// TODO: Just for test, delete
#define VIENNACL_HAVE_UBLAS 1

// Project includes
#include "includes/model_part.h"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/gmres.hpp"

// TODO: Just for test, delete
#include "includes/matrix_market_interface.h"

// TODO: Remove unneeded ones
//#include "includes/node.h"
//#include "geometries/geometry.h"
//#include "utilities/geometry_utilities.h"
//#include "incompressible_fluid_application.h"


namespace Kratos
{
	//
	// ViennaCLHelper
	//
	// ViennaCL helper class

	class ViennaCLHelper
	{
		public:

			//
			// ViennaCLHelper
			//
			// Constructor

			ViennaCLHelper(OpenCL::DeviceGroup &_DeviceGroup)
			{
				// Set our own context, device and command queue
				viennacl::ocl::setup_context(0, _DeviceGroup.Contexts[0], _DeviceGroup.DeviceIDs[0], _DeviceGroup.CommandQueues[0]);
			}
	};

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

			// TODO: Check for effect of alignment in GPU types while optimizing the performance

			typedef CompressedMatrix HostMatrixType;
			typedef viennacl::compressed_matrix <double> DeviceMatrixType;

			typedef Vector HostVectorType;
			typedef viennacl::vector <double> DeviceVectorType;

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
				mViennaCLHelper(matrix_container.GetDeviceGroup()),
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

				// Loading OpenCL program and defining appropriate macros
				// TODO: Remove if not needed (also remove in .cl file)
				mpOpenCLFluidSolver = mrDeviceGroup.BuildProgramFromFile(
					"opencl_fluid_solver.cl",
					"-cl-fast-relaxed-math"

#ifdef SPLIT_OSS

					" -DSPLIT_OSS"

#endif

#ifdef SYMM_PRESS

					" -DSYMM_PRESS"

#endif

				);

				// Register kernels
				mkAddVectorInplace = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "AddVectorInplace");
				mkSubVectorInplace = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "SubVectorInplace");

				mkSolveStep1_1 = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "SolveStep1_1");
				mkSolveStep1_2 = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "SolveStep1_2");

				mkSolveStep2_1 = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "SolveStep2_1");
				mkSolveStep2_2 = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "SolveStep2_2");
				mkSolveStep2_3 = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "SolveStep2_3");

				mkSolveStep3_1 = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "SolveStep3_1");
				mkSolveStep3_2 = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "SolveStep3_2");

				mkCalculateRHS = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "CalculateRHS");

				mkComputeWallResistance = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "ComputeWallResistance");

				mkApplyVelocityBC_1 = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "ApplyVelocityBC_1");
				mkApplyVelocityBC_2 = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "ApplyVelocityBC_2");
				mkApplyVelocityBC_3 = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "ApplyVelocityBC_3");
				mkApplyVelocityBC_4 = mrDeviceGroup.RegisterKernel(mpOpenCLFluidSolver, "ApplyVelocityBC_4");
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

				// Get no. of nodes and edges
				n_nodes = mr_model_part.Nodes().size();
				n_edges = mr_matrix_container.GetNumberEdges();

				// TODO: Check if these are all needed

				// Size data vectors
				AllocateArray(&mWork, n_nodes);
				AllocateArray(&mvel_n, n_nodes);
				AllocateArray(&mvel_n1, n_nodes);
				AllocateArray(&mPn, n_nodes);
				AllocateArray(&mPn1, n_nodes);
				//AllocateArray(&mHmin, n_nodes);
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
				AllocateArray(&mrhs, n_nodes);  // TODO: It seems that this is not needed

				// Allocating lists
				// TODO: Maximum size is used; if this seems a problem with OpenCL buffers, try fixing this
				// TODO: Decide where to copy these lists to GPU
				AllocateArray(&mSlipNormal, n_nodes);

				AllocateArray(&mSlipBoundaryList, n_nodes);
				AllocateArray(&mPressureOutletList, n_nodes);
				AllocateArray(&mFixedVelocitiesList, n_nodes);

				AllocateArray(&mFixedVelocitiesValuesList, n_nodes);


				AllocateArray(&medge_nodesList, n_nodes);
				AllocateArray(&medge_nodes_directionList, n_nodes);
				AllocateArray(&mcorner_nodesList, n_nodes);


				// Allocating buffers on OpenCL device
				mbWork = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);
				mbvel_n = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);
				mbvel_n1 = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);

				mbPn = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);
				mbPn1 = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);
				mbHmin = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);
				mbHavg = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);

				mbNodalFlag = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);  // TODO: Not used properly

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

				// Lists
				mbSlipNormal = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);

				mbSlipBoundaryList = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_uint), CL_MEM_READ_WRITE);
				mbPressureOutletList = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_uint), CL_MEM_READ_WRITE);
				mbFixedVelocitiesList = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_uint), CL_MEM_READ_WRITE);

				mbFixedVelocitiesValuesList = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);


				mbedge_nodesList = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_uint), CL_MEM_READ_WRITE);
				mbedge_nodes_directionList = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);
				mbcorner_nodesList = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_uint), CL_MEM_READ_WRITE);


				// Lists' lengths
				mSlipBoundaryListLength = 0;
				mPressureOutletListLength = 0;
				mFixedVelocitiesListLength = 0;

				medge_nodesListLength = 0;
				medge_nodes_directionListLength = 0;
				mcorner_nodesListLength = 0;


				mr_matrix_container.SetToZero(mbdiv_error);

				// Read velocity and pressure data from Kratos

				// TODO: Can these be host only?
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes(), mbvel_n1);
				mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mvel_n, mr_model_part.Nodes(), mbvel_n);

				mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, mr_model_part.Nodes(), mbPn1);
				mr_matrix_container.FillOldScalarFromDatabase(PRESSURE, mPn, mr_model_part.Nodes(), mbPn);

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
				mr_matrix_container.WriteVectorToDatabase(NORMAL, mSlipNormal, mr_model_part.Nodes(), mbSlipNormal);  // TODO: Check what type of update is needed

				DetectEdges3D(mr_model_part.Conditions());

				// Determine number of edges and entries
				unsigned int n_nonzero_entries = 2 * n_edges + n_nodes;

				// Allocate memory for variables
				rhs_GPU.resize(n_nodes);
				dp_GPU.resize(n_nodes);

				mL.resize(n_nodes, n_nodes, n_nonzero_entries);

				// TODO: May we assume that the graph of mL is the same as the graph of mr_matrix_container? What does the flag do?

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
							mL.push_back(i_node, i_node, 1.00); // TODO: Just for test, was: 0.00);
							flag = 1;
						}

						// Add non-diagonal/edge contribution
						mL.push_back(i_node, j_neighbour, 0.00);
					}

					// If diagonal element is the last non-zero element of the row
					if (flag == 0)
					{
						mL.push_back(i_node, i_node, 1.00); // TODO: Just for test, was: 0.00);
					}
				}

				// Copy mL to GPU
				copy(mL, mL_GPU);

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

				// TODO: Is here OK to copy to GPU?
				mrDeviceGroup.CopyBuffer(mbFixedVelocitiesList, OpenCL::HostToDevice, OpenCL::VoidPList(1, mFixedVelocitiesList));
				mrDeviceGroup.CopyBuffer(mbFixedVelocitiesValuesList, OpenCL::HostToDevice, OpenCL::VoidPList(1, mFixedVelocitiesValuesList));
				mrDeviceGroup.CopyBuffer(mbPressureOutletList, OpenCL::HostToDevice, OpenCL::VoidPList(1, mPressureOutletList));

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
					const double hmin_i = mr_matrix_container.GetHmin()[i_node];

					double vel_norm = Norm2_3(mvel_n1[i_node]);


					// Use CFL condition to compute time step size
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

                unsigned int fixed_size = mFixedVelocitiesListLength; // TODO: Why firstprivate? //mFixedVelocities.size();

                #pragma omp parallel for firstprivate(fixed_size)
                for (unsigned int i_velocity = 0; i_velocity < fixed_size; i_velocity++)
                {
                    unsigned int i_node = mFixedVelocitiesList[i_velocity];
                    mFixedVelocitiesValuesList[i_velocity] = mvel_n1[i_node];
                }

				// TODO: Should we update anything on GPU?
				mrDeviceGroup.CopyBuffer(mbFixedVelocitiesList, OpenCL::HostToDevice, OpenCL::VoidPList(1, mFixedVelocitiesList));

                KRATOS_CATCH("")
            }

            //
            // CalculateRHS
            //
            // Function to calculate right-hand side of fractional momentum equation

			void CalculateRHS(cl_uint vel_buffer, cl_uint pressure_buffer, cl_uint convective_velocity_buffer, cl_uint rhs_buffer)
			{
				KRATOS_TRY

				// TODO: Is this OK?
				//int n_nodes = vel.size();

				double inverse_rho = 1.0 / mRho;

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS, 0, mbPi);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS, 1, vel_buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS, 2, mr_matrix_container.GetRowStartIndexBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS, 3, mr_matrix_container.GetColumnIndexBuffer());
				mrDeviceGroup.SetImageAsKernelArg(mkCalculateRHS, 4, mr_matrix_container.GetEdgeValuesBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS, 5, mr_matrix_container.GetLumpedMassBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS, 6, convective_velocity_buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS, 7, pressure_buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS, 8, rhs_buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS, 9, mbTauConvection);
				mrDeviceGroup.SetKernelArg(mkCalculateRHS, 10, mBodyForce[0]);
				mrDeviceGroup.SetKernelArg(mkCalculateRHS, 11, mBodyForce[1]);
				mrDeviceGroup.SetKernelArg(mkCalculateRHS, 12, mBodyForce[2]);
				mrDeviceGroup.SetKernelArg(mkCalculateRHS, 13, inverse_rho);
				mrDeviceGroup.SetKernelArg(mkCalculateRHS, 14, mViscosity);
				mrDeviceGroup.SetKernelArg(mkCalculateRHS, 15, n_nodes);
				mrDeviceGroup.SetLocalMemAsKernelArg(mkCalculateRHS, 16, (mrDeviceGroup.WorkGroupSizes[mkCalculateRHS][0] + 1) * sizeof(cl_uint));

				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkCalculateRHS, n_nodes);


				// Apply wall resistance
				if (mWallLawIsActive == true)
				{
					ComputeWallResistance(vel_buffer, rhs_buffer);
				}

				ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
				mr_matrix_container.WriteVectorToDatabase(VELOCITY, mvel_n1, rNodes, mbvel_n1);

				KRATOS_CATCH("")
			}

			//
			// ApplyVelocityBC
			//

			void ApplyVelocityBC(cl_uint VelArray_buffer)
			{
				KRATOS_TRY

				if (mWallLawIsActive == false)
				{
					// Apply conditions on corner edges

					// Calling ApplyVelocityBC_1 OpenCL kernel

					// Setting arguments
					mrDeviceGroup.SetBufferAsKernelArg(mkApplyVelocityBC_1, 0, VelArray_buffer);
					mrDeviceGroup.SetBufferAsKernelArg(mkApplyVelocityBC_1, 1, mbedge_nodes_directionList);
					mrDeviceGroup.SetBufferAsKernelArg(mkApplyVelocityBC_1, 2, mbedge_nodesList);
					mrDeviceGroup.SetKernelArg(mkApplyVelocityBC_1, 3, medge_nodes_directionListLength);


					// Execute OpenCL kernel
					mrDeviceGroup.ExecuteKernel(mkApplyVelocityBC_1, medge_nodes_directionListLength);


					// Apply conditions on corners

					// Calling ApplyVelocityBC_2 OpenCL kernel

					// Setting arguments
					mrDeviceGroup.SetBufferAsKernelArg(mkApplyVelocityBC_2, 0, VelArray_buffer);
					mrDeviceGroup.SetBufferAsKernelArg(mkApplyVelocityBC_2, 1, mbcorner_nodesList);
					mrDeviceGroup.SetKernelArg(mkApplyVelocityBC_2, 2, mcorner_nodesListLength);


					// Execute OpenCL kernel
					mrDeviceGroup.ExecuteKernel(mkApplyVelocityBC_2, mcorner_nodesListLength);
				}

				// Slip condition

				// Calling ApplyVelocityBC_3 OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkApplyVelocityBC_3, 0, VelArray_buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkApplyVelocityBC_3, 1, mbSlipNormal);
				mrDeviceGroup.SetBufferAsKernelArg(mkApplyVelocityBC_3, 2, mbSlipBoundaryList);
				mrDeviceGroup.SetKernelArg(mkApplyVelocityBC_3, 3, mSlipBoundaryListLength);


				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkApplyVelocityBC_3, mSlipBoundaryListLength);


				// Fixed condition

				// Calling ApplyVelocityBC_4 OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkApplyVelocityBC_4, 0, VelArray_buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkApplyVelocityBC_4, 1, mbFixedVelocitiesValuesList);
				mrDeviceGroup.SetBufferAsKernelArg(mkApplyVelocityBC_4, 2, mbFixedVelocitiesList);
				mrDeviceGroup.SetKernelArg(mkApplyVelocityBC_4, 3, mFixedVelocitiesListLength);


				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkApplyVelocityBC_4, mFixedVelocitiesListLength);

				KRATOS_CATCH("")
			}

            //
            // SolveStep1
            //
            // Function to solve fluid equations - fractional step 1: compute fractional momentum

            void SolveStep1()
            {
				KRATOS_TRY

				// TODO: Test only, delete!
				//for (unsigned int i = 0; i < n_nodes; i++)
				//{
				//	std::cout << mr_matrix_container.GetInvertedMass()[i] << "  ";
				//}
				//std::cout << std::endl;

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

				// Calling Solve1_1 OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep1_1, 0, mbHavg);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep1_1, 1, mbvel_n1);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep1_1, 2, mbTauPressure);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep1_1, 3, mbTauConvection);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep1_1, 4, mbTau2);
				mrDeviceGroup.SetKernelArg(mkSolveStep1_1, 5, mViscosity);
				mrDeviceGroup.SetKernelArg(mkSolveStep1_1, 6, time_inv_avg);
				mrDeviceGroup.SetKernelArg(mkSolveStep1_1, 7, mstabdt_pressure_factor);
				mrDeviceGroup.SetKernelArg(mkSolveStep1_1, 8, mstabdt_convection_factor);
				mrDeviceGroup.SetKernelArg(mkSolveStep1_1, 9, mtau2_factor);
				mrDeviceGroup.SetKernelArg(mkSolveStep1_1, 10, n_nodes);

				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkSolveStep1_1, n_nodes);

				// Calling Solve1_2 OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep1_2, 0, mbPi);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep1_2, 1, mbvel_n1);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep1_2, 2, mr_matrix_container.GetRowStartIndexBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep1_2, 3, mr_matrix_container.GetColumnIndexBuffer());
				mrDeviceGroup.SetImageAsKernelArg(mkSolveStep1_2, 4, mr_matrix_container.GetEdgeValuesBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep1_2, 5, mr_matrix_container.GetInvertedMassBuffer());
				mrDeviceGroup.SetKernelArg(mkSolveStep1_2, 6, n_nodes);
				mrDeviceGroup.SetLocalMemAsKernelArg(mkSolveStep1_2, 7, (mrDeviceGroup.WorkGroupSizes[mkSolveStep1_2][0] + 1) * sizeof(cl_uint));

				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkSolveStep1_2, n_nodes);


				mr_matrix_container.AssignVectorToVector(mbvel_n, mbWork);  // mWork = mvel_n

				// First step of Runge Kutta
				mr_matrix_container.AssignVectorToVector(mbvel_n, mbvel_n1);  // mvel_n1 = mvel_n

				mr_matrix_container.SetToZero(mbrhs);
				CalculateRHS(mbvel_n1, mbPn, mbvel_n1, mbrhs);
				mr_matrix_container.Add_Minv_value3(mbWork, mbWork, delta_t / 6.00, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);
				mr_matrix_container.Add_Minv_value3(mbvel_n1, mbvel_n, 0.5 * delta_t, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);
				ApplyVelocityBC(mbvel_n1);

				// TODO: Debugging only, delete this!
				mrDeviceGroup.CopyBuffer(mbvel_n1, OpenCL::DeviceToHost, OpenCL::VoidPList(1, mvel_n1));
				KRATOS_WATCH("end of first stage");

				double vnorm2 = 0.00;
				for (unsigned int i = 0; i < n_nodes; i++)
				{
					vnorm2 += pow(KRATOS_OCL_COMP_0(mvel_n1[i]), 2) + pow(KRATOS_OCL_COMP_1(mvel_n1[i]), 2) + pow(KRATOS_OCL_COMP_2(mvel_n1[i]), 2);
				}

				KRATOS_WATCH(sqrt(vnorm2));

				// Second step
				mr_matrix_container.SetToZero(mbrhs);
				CalculateRHS(mbvel_n1, mbPn, mbvel_n1, mbrhs);
				mr_matrix_container.Add_Minv_value3(mbWork, mbWork, delta_t / 3.00, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);
				mr_matrix_container.Add_Minv_value3(mbvel_n1, mbvel_n, 0.5 * delta_t, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);
				ApplyVelocityBC(mbvel_n1);

				// TODO: Debugging only, delete this!
				mrDeviceGroup.CopyBuffer(mbvel_n1, OpenCL::DeviceToHost, OpenCL::VoidPList(1, mvel_n1));
				KRATOS_WATCH("end of second stage");

                                vnorm2 = 0.00;
				for (unsigned int i = 0; i < n_nodes; i++)
				{
					vnorm2 += pow(KRATOS_OCL_COMP_0(mvel_n1[i]), 2) + pow(KRATOS_OCL_COMP_1(mvel_n1[i]), 2) + pow(KRATOS_OCL_COMP_2(mvel_n1[i]), 2);
				}

				KRATOS_WATCH(sqrt(vnorm2));

				// Third step
				mr_matrix_container.SetToZero(mbrhs);
				CalculateRHS(mbvel_n1, mbPn, mbvel_n1, mbrhs);
				mr_matrix_container.Add_Minv_value3(mbWork, mbWork, delta_t / 3.00, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);
				mr_matrix_container.Add_Minv_value3(mbvel_n1, mbvel_n, delta_t, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);
				ApplyVelocityBC(mbvel_n1);

				// TODO: Debugging only, delete this!
				mrDeviceGroup.CopyBuffer(mbvel_n1, OpenCL::DeviceToHost, OpenCL::VoidPList(1, mvel_n1));
				KRATOS_WATCH("end of third stage");

                                vnorm2 = 0.00;
				for (unsigned int i = 0; i < n_nodes; i++)
				{
					vnorm2 += pow(KRATOS_OCL_COMP_0(mvel_n1[i]), 2) + pow(KRATOS_OCL_COMP_1(mvel_n1[i]), 2) + pow(KRATOS_OCL_COMP_2(mvel_n1[i]), 2);
				}

				KRATOS_WATCH(sqrt(vnorm2));

				// Fourth step
				mr_matrix_container.SetToZero(mbrhs);
				CalculateRHS(mbvel_n1, mbPn, mbvel_n1, mbrhs);
				mr_matrix_container.Add_Minv_value3(mbWork, mbWork, delta_t / 6.00, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);

				// Compute right-hand side
				mr_matrix_container.AssignVectorToVector(mbWork, mbvel_n1);
				ApplyVelocityBC(mbvel_n1);

				// TODO: Debugging only, delete this!
				mrDeviceGroup.CopyBuffer(mbvel_n1, OpenCL::DeviceToHost, OpenCL::VoidPList(1, mvel_n1));
				KRATOS_WATCH("end of step1");

                                 vnorm2 = 0.00;
				for (unsigned int i = 0; i < n_nodes; i++)
				{
					vnorm2 += pow(KRATOS_OCL_COMP_0(mvel_n1[i]), 2) + pow(KRATOS_OCL_COMP_1(mvel_n1[i]), 2) + pow(KRATOS_OCL_COMP_2(mvel_n1[i]), 2);
				}

				KRATOS_WATCH(sqrt(vnorm2));

				KRATOS_CATCH("")
            }

            //
            // SolveStep2
            //
            // Function to solve fluid equations - fractional step 2: calculate pressure

			void SolveStep2()  // TODO: Fix this! Should we get the linear solver as an argument? //typename TLinearSolver::Pointer pLinearSolver
			{

				KRATOS_TRY

				// Prerequisites

				// Allocate memory for variables
				ModelPart::NodesContainerType &rNodes = mr_model_part.Nodes();

				// TODO: Do the resize inside the Initialize
				//dp.resize(n_nodes);
				//rhs.resize(n_nodes);

				// Read time step size from Kratos
				ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
				double delta_t = CurrentProcessInfo[DELTA_TIME];

				// Read the pressure projection from the database
				mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, rNodes, mbPn1);  // TODO: Is this OK? //mr_model_part.Nodes()
				mr_matrix_container.FillVectorFromDatabase(PRESS_PROJ, mXi, rNodes, mbXi);
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes, mbvel_n1);

				// Calling Solve2_1 OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep2_1, 0, mbvel_n1);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep2_1, 1, mbXi);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep2_1, 2, mbTauPressure);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep2_1, 3, mbPn);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep2_1, 4, mbPn1);
				mrDeviceGroup.SetKernelArg(mkSolveStep2_1, 5, mL_GPU.handle1());
				mrDeviceGroup.SetKernelArg(mkSolveStep2_1, 6, mL_GPU.handle2());
				mrDeviceGroup.SetImageAsKernelArg(mkSolveStep2_1, 7, mr_matrix_container.GetEdgeValuesBuffer());
				mrDeviceGroup.SetKernelArg(mkSolveStep2_1, 8, rhs_GPU.handle());
				mrDeviceGroup.SetKernelArg(mkSolveStep2_1, 9, mL_GPU.handle());
				mrDeviceGroup.SetKernelArg(mkSolveStep2_1, 10, mRho);
				mrDeviceGroup.SetKernelArg(mkSolveStep2_1, 11, delta_t);
				mrDeviceGroup.SetKernelArg(mkSolveStep2_1, 12, n_nodes);
				mrDeviceGroup.SetLocalMemAsKernelArg(mkSolveStep2_1, 13, (mrDeviceGroup.WorkGroupSizes[mkSolveStep2_1][0] + 1) * sizeof(cl_uint));

				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkSolveStep2_1, n_nodes);


				if (muse_mass_correction == true)
				{
					// Calling SubVectorInplace OpenCL kernel

					// Setting arguments
					mrDeviceGroup.SetKernelArg(mkSubVectorInplace, 0, rhs_GPU.handle());
					mrDeviceGroup.SetBufferAsKernelArg(mkSubVectorInplace, 1, mbdiv_error);
					mrDeviceGroup.SetKernelArg(mkSubVectorInplace, 2, n_nodes);

					// Execute OpenCL kernel
					mrDeviceGroup.ExecuteKernel(mkSubVectorInplace, n_nodes);
				}

				// Calling SubVectorInplace OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep2_2, 0, mbPressureOutletList);
				mrDeviceGroup.SetKernelArg(mkSolveStep2_2, 1, mL_GPU.handle1());
				mrDeviceGroup.SetKernelArg(mkSolveStep2_2, 2, mL_GPU.handle2());
				mrDeviceGroup.SetKernelArg(mkSolveStep2_2, 3, mL_GPU.handle());
				mrDeviceGroup.SetKernelArg(mkSolveStep2_2, 4, rhs_GPU.handle());
				mrDeviceGroup.SetKernelArg(mkSolveStep2_2, 5, mPressureOutletListLength);


				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkSolveStep2_2, mPressureOutletListLength);


				// TODO: Is this a good thing to do? Can we start from last dp instead?
				// TODO: Maybe we do not need this, because of the way ViennaCL solver is used

				// Set starting vector for iterative solvers
				//dp_GPU.clear();

				// TODO: For debugging ONLY! Delete it!
				HostVectorType rhs(n_nodes);
				viennacl::copy(rhs_GPU, rhs);
				KRATOS_WATCH(norm_2(rhs));

				copy(mL_GPU, mL);
				KRATOS_WATCH(matrix_norm_frobenius <HostMatrixType> :: apply(mL));

				WriteMatrixMarketMatrix("mL.mm", mL, false);
				WriteMatrixMarketVector("rhs.mm", rhs);

				// Calling the ViennaCL solver
				dp_GPU = viennacl::linalg::solve(mL_GPU, rhs_GPU, viennacl::linalg::gmres_tag());  // TODO: Is this OK to hard-code BiCGStab?

				// TODO: For debugging ONLY! Delete it!
				HostVectorType dp(n_nodes);
				HostVectorType dp2(n_nodes);
				viennacl::copy(dp_GPU, dp);
				KRATOS_WATCH(norm_2(dp));

				dp2 = viennacl::linalg::solve(mL, rhs, viennacl::linalg::gmres_tag());  // TODO: Just for test, delete
				WriteMatrixMarketVector("dp2.mm", dp2);
				KRATOS_WATCH(norm_2(dp2));

				// Update pressure

				// Calling AddVectorInplace OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkAddVectorInplace, 0, mbPn1);
				mrDeviceGroup.SetKernelArg(mkAddVectorInplace, 1, dp_GPU.handle());
				mrDeviceGroup.SetKernelArg(mkAddVectorInplace, 2, n_nodes);

				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkAddVectorInplace, n_nodes);

				// Write pressure and density to Kratos
				mr_matrix_container.WriteScalarToDatabase(PRESSURE, mPn1, rNodes, mbPn1);

				// Compute pressure proj for the next step

				// Calling SolveStep2_3 OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep2_3, 0, mbXi);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep2_3, 1, mbPn1);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep2_3, 2, mr_matrix_container.GetRowStartIndexBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep2_3, 3, mr_matrix_container.GetColumnIndexBuffer());
				mrDeviceGroup.SetImageAsKernelArg(mkSolveStep2_3, 4, mr_matrix_container.GetEdgeValuesBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep2_3, 5, mr_matrix_container.GetInvertedMassBuffer());
				mrDeviceGroup.SetKernelArg(mkSolveStep2_3, 6, n_nodes);
				mrDeviceGroup.SetLocalMemAsKernelArg(mkSolveStep2_3, 7, (mrDeviceGroup.WorkGroupSizes[mkSolveStep2_3][0] + 1) * sizeof(cl_uint));

				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkSolveStep2_3, n_nodes);


				mr_matrix_container.WriteVectorToDatabase(PRESS_PROJ, mXi, rNodes, mbXi);

				KRATOS_CATCH("")
			}

			//
			// SolveStep3
			//
			// Function to solve fluid equations - fractional step 3: correct fractional momentum

			void SolveStep3()
			{
				KRATOS_TRY

				ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();

				// Read time step size from Kratos
				ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
				double delta_t = CurrentProcessInfo[DELTA_TIME];

				double factor = 0.5;

				if (massume_constant_dp == true)
				{
					factor = 1.00;
				}

				// Compute end of step momentum
				double rho_inv = 1.00 / mRho;

				// Calling SolveStep3_1 OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep3_1, 0, mbvel_n1);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep3_1, 1, mbPn);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep3_1, 2, mbPn1);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep3_1, 3, mr_matrix_container.GetRowStartIndexBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep3_1, 4, mr_matrix_container.GetColumnIndexBuffer());
				mrDeviceGroup.SetImageAsKernelArg(mkSolveStep3_1, 5, mr_matrix_container.GetEdgeValuesBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep3_1, 6, mr_matrix_container.GetInvertedMassBuffer());
				mrDeviceGroup.SetKernelArg(mkSolveStep3_1, 7, rho_inv);
				mrDeviceGroup.SetKernelArg(mkSolveStep3_1, 8, factor);
				mrDeviceGroup.SetKernelArg(mkSolveStep3_1, 9, delta_t);
				mrDeviceGroup.SetKernelArg(mkSolveStep3_1, 10, n_nodes);
				mrDeviceGroup.SetLocalMemAsKernelArg(mkSolveStep3_1, 11, (mrDeviceGroup.WorkGroupSizes[mkSolveStep3_1][0] + 1) * sizeof(cl_uint));

				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkSolveStep3_1, n_nodes);


				ApplyVelocityBC(mbvel_n1);

				// Write velocity of time step n+1 to Kratos
				mr_matrix_container.WriteVectorToDatabase(VELOCITY, mvel_n1, rNodes, mbvel_n1);


				// Calculate the error on the divergence
				if (muse_mass_correction == true)
				{
					// Calling SolveStep3_2 OpenCL kernel

					// Setting arguments
					mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep3_2, 0, mbvel_n1);
					mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep3_2, 1, mbdiv_error);
					mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep3_2, 2, mr_matrix_container.GetRowStartIndexBuffer());
					mrDeviceGroup.SetBufferAsKernelArg(mkSolveStep3_2, 3, mr_matrix_container.GetColumnIndexBuffer());
					mrDeviceGroup.SetImageAsKernelArg(mkSolveStep3_2, 4, mr_matrix_container.GetEdgeValuesBuffer());
					mrDeviceGroup.SetKernelArg(mkSolveStep3_2, 5, mRho);
					mrDeviceGroup.SetKernelArg(mkSolveStep3_2, 6, n_nodes);
					mrDeviceGroup.SetLocalMemAsKernelArg(mkSolveStep3_2, 7, (mrDeviceGroup.WorkGroupSizes[mkSolveStep3_2][0] + 1) * sizeof(cl_uint));

					// Execute OpenCL kernel
					mrDeviceGroup.ExecuteKernel(mkSolveStep3_2, n_nodes);
				}

				// TODO: Debugging only, delete this!
				mrDeviceGroup.CopyBuffer(mbvel_n1, OpenCL::DeviceToHost, OpenCL::VoidPList(1, mvel_n1));
				KRATOS_WATCH("end of step3");

				double vnorm2 = 0.00;
				for (unsigned int i = 0; i < n_nodes; i++)
				{
					vnorm2 += pow(KRATOS_OCL_COMP_0(mvel_n1[i]), 2) + pow(KRATOS_OCL_COMP_1(mvel_n1[i]), 2) + pow(KRATOS_OCL_COMP_2(mvel_n1[i]), 2);
				}

				KRATOS_WATCH(sqrt(vnorm2));


				KRATOS_CATCH("")
			}

			//
			//
			//

			void ComputeWallResistance(cl_uint vel_buffer, cl_uint rhs_buffer)
			{
				// Parameters:
				double density = mRho;
				double mu = mViscosity;
				double ym = mY_wall;

				unsigned int slip_size = mSlipBoundaryListLength;

				if (mu == 0.00)
				{
					KRATOS_ERROR(std::logic_error, "It is not possible to use the wall law with zero viscosity", "");
				}

				// Slip condition

				// Calling ComputeWallResistance OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkComputeWallResistance, 0, vel_buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkComputeWallResistance, 1, rhs_buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkComputeWallResistance, 2, mbSlipNormal);  // TODO: This is not allocated yet
				mrDeviceGroup.SetBufferAsKernelArg(mkComputeWallResistance, 3, mbSlipBoundaryList);  // TODO: This is not allocated yet
				mrDeviceGroup.SetKernelArg(mkComputeWallResistance, 4, density);
				mrDeviceGroup.SetKernelArg(mkComputeWallResistance, 5, mu);
				mrDeviceGroup.SetKernelArg(mkComputeWallResistance, 6, ym);
				mrDeviceGroup.SetKernelArg(mkComputeWallResistance, 7, slip_size);


				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkSolveStep3_1, slip_size);
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

				// TODO: Is here OK to copy to GPU?
				mrDeviceGroup.CopyBuffer(mbSlipNormal, OpenCL::HostToDevice, OpenCL::VoidPList(1, mSlipNormal));  // TODO: Is this needed?
				mrDeviceGroup.CopyBuffer(mbSlipBoundaryList, OpenCL::HostToDevice, OpenCL::VoidPList(1, mSlipBoundaryList));

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
				// TODO: mL, dp, rhs and their GPU counterparts

				/*mWork.clear();
				mvel_n.clear();
				mvel_n1.clear();
				mPn.clear();
				mPn1.clear();
				//mHmin.clear();
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
			cl_uint mbWork, mbvel_n, mbvel_n1, mbPn, mbPn1, mbHmin, mbHavg, mbNodalFlag, mbTauPressure, mbTauConvection, mbTau2, mbPi, mbXi, mbx, mbEdgeDimensions, mbBeta, mbdiv_error, mbSlipNormal, mbSlipBoundaryList, mbPressureOutletList, mbedge_nodes_directionList, mbedge_nodesList, mbcorner_nodesList, mbFixedVelocitiesList, mbFixedVelocitiesValuesList, mbrhs;

			cl_uint mpOpenCLFluidSolver, mkAddVectorInplace, mkSubVectorInplace, mkSolveStep1_1, mkSolveStep1_2, mkSolveStep2_1, mkSolveStep2_2, mkSolveStep2_3, mkSolveStep3_1, mkSolveStep3_2, mkCalculateRHS, mkComputeWallResistance, mkApplyVelocityBC_1, mkApplyVelocityBC_2, mkApplyVelocityBC_3, mkApplyVelocityBC_4;

			// Matrix container
			OpenCLMatrixContainer &mr_matrix_container;

			// Associated model part
			ModelPart &mr_model_part;

			// ViennaCL helper
			ViennaCLHelper mViennaCLHelper;

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
			//ValuesVectorType mHmin;
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
			IndicesVectorType mSlipBoundaryList, mPressureOutletList, mFixedVelocitiesList;
			CalcVectorType mFixedVelocitiesValuesList;

			unsigned int mSlipBoundaryListLength, mPressureOutletListLength, mFixedVelocitiesListLength;

			// Intrinsic time step size
			ValuesVectorType mTauPressure;
			ValuesVectorType mTauConvection;
			ValuesVectorType mTau2;

			ValuesVectorType mdiv_error;

			// Storage of nodal values in local variables
			CalcVectorType mrhs;  // TODO: It seems that this is not needed

			// Variables for resolving pressure equation

			// Laplacian matrix
			HostMatrixType mL;
			DeviceMatrixType mL_GPU;

			// Vectors on GPU
			DeviceVectorType dp_GPU, rhs_GPU;

			// Constant variables
			double mRho;
			double mViscosity;
			array_1d <double, 3> mBodyForce;

			// Variables for convection
			ValuesVectorType mBeta;

			// Variables for edge BCs
			IndicesVectorType medge_nodesList;
			CalcVectorType medge_nodes_directionList;
			IndicesVectorType mcorner_nodesList;

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
					noalias(position[i_node]) = node_it -> Coordinates();

					// Initialize minimum edge length with relatively big values
					//mHmin[i_node] = 1e10;
				}

				//ValuesVectorType TempHmin = mr_matrix_container.GetHmin();
				//for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				//{
				//	mHmin[i_node] = TempHmin[i_node];
				//}

				// Take unstructured meshes into account
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					double &h_i = mHavg[i_node];
					double &m_i = mr_matrix_container.GetLumpedMass()[i_node];

					h_i = pow(6.0 * m_i, 1.0 / 3.0);
				}

                                //transfer mHavg to gpu
                                mrDeviceGroup.CopyBuffer(mbHavg, OpenCL::HostToDevice, OpenCL::VoidPList(1, mHavg));

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
					temp_edge_nodes[i_node] = 0;
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
				for (unsigned int i = 0; i < mcorner_nodesListLength; i++)
				{
					KRATOS_WATCH(mcorner_nodesList[i]);
				}

				// TODO: Is here OK to copy to GPU?
				mrDeviceGroup.CopyBuffer(mbedge_nodesList, OpenCL::HostToDevice, OpenCL::VoidPList(1, medge_nodesList));
				mrDeviceGroup.CopyBuffer(mbedge_nodes_directionList, OpenCL::HostToDevice, OpenCL::VoidPList(1, medge_nodes_directionList));
				mrDeviceGroup.CopyBuffer(mbcorner_nodesList, OpenCL::HostToDevice, OpenCL::VoidPList(1, mcorner_nodesList));

				KRATOS_CATCH("")
			}
	};

}  // Namespace Kratos

#endif  // KRATOS_OPENCL_EDGEBASED_LEVELSET_FLUID_SOLVER_H_INCLUDED
