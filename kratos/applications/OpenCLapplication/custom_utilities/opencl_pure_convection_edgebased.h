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
//   Date:                $Date: 2010-09-30 18:13:33 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_OPENCL_PURE_CONVECTION_EDGEBASED_SOLVER_H_INCLUDED)
#define KRATOS_OPENCL_PURE_CONVECTION_EDGEBASED_SOLVER_H_INCLUDED

// TODO: What are these?
#define SPLIT_OSS
//#define SYMM_PRESS


// System includes
#include <string>
#include <iostream>
#include <algorithm>


// External includes


// Project includes
#include "includes/model_part.h"


int64_t timeNanos()
{
	struct timespec tp;

	clock_gettime(CLOCK_MONOTONIC, &tp);
	return (unsigned long long) tp.tv_sec * (1000ULL * 1000ULL * 1000ULL) + (unsigned long long) tp.tv_nsec;
}

namespace Kratos
{

	//
	// OpenCLPureConvectionEdgeBased3D
	//
	// OpenCL based pure convection edge based solver

	class OpenCLPureConvectionEdgeBased3D
	{
		public:

			//
			// Used types

			typedef cl_double3 *CalcVectorType;
			typedef cl_double *ValuesVectorType;

			//
			// OpenCLPureConvectionEdgeBased3D
			//
			// Constructor

			OpenCLPureConvectionEdgeBased3D(OpenCLMatrixContainer &matrix_container, ModelPart &model_part): mr_matrix_container(matrix_container), mrDeviceGroup(mr_matrix_container.GetDeviceGroup()), mr_model_part(model_part)
			{
				// Loading OpenCL program
				// TODO: Add optimization flags here
				mpOpenCLPureConvectionEdgeBased = mrDeviceGroup.BuildProgramFromFile("opencl_pure_convection_edgebased.cl", "-cl-fast-relaxed-math");

				// Register kernels
				mkSolve1 = mrDeviceGroup.RegisterKernel(mpOpenCLPureConvectionEdgeBased, "Solve1");
				mkCalculateRHS1 = mrDeviceGroup.RegisterKernel(mpOpenCLPureConvectionEdgeBased, "CalculateRHS1");
				mkCalculateRHS2 = mrDeviceGroup.RegisterKernel(mpOpenCLPureConvectionEdgeBased, "CalculateRHS2");
				mkCalculateRHS3 = mrDeviceGroup.RegisterKernel(mpOpenCLPureConvectionEdgeBased, "CalculateRHS3");
				mkCalculateAdvectiveVelocity = mrDeviceGroup.RegisterKernel(mpOpenCLPureConvectionEdgeBased, "CalculateAdvectiveVelocity");
			}

			//
			// ~OpenCLPureConvectionEdgeBased3D
			//
			// Destructor

			~OpenCLPureConvectionEdgeBased3D()
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

				// TODO: Allocate a single chunk of memory for variables
				// TODO: Use Page-locked memory for faster data transfer to GPU
				// TODO: Order variables, such that variables copied to GPU together are together here too
				// TODO: Account for device address alignment

				// Size data vectors
				AllocateArray(&mWork, n_nodes);

				AllocateArray(&mPi, n_nodes);

				AllocateArray(&mUn, n_nodes);
				AllocateArray(&mUn1, n_nodes);

				AllocateArray(&mphi_n, n_nodes);
				AllocateArray(&mphi_n1, n_nodes);

				AllocateArray(&mA, n_nodes);

				AllocateArray(&mTau, n_nodes);
				AllocateArray(&mBeta, n_nodes);

				AllocateArray(&mx, n_nodes);

				AllocateArray(&mrhs, n_nodes);

				// Allocating buffers on OpenCL device
				mbWork = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);

				mbPi = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);

				mbUn = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);
				mbUn1 = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);

				mbphi_n = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);
				mbphi_n1 = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);

				mbA = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);

				mbTau = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);
				mbBeta = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);

				mbx = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_WRITE);

				mbrhs = mrDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_WRITE);

				// Read variables from database

				// TODO: It seems that we do not need this as Solve() does this at first step
				/*
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mUn1, mr_model_part.Nodes(), mbUn1);
				mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mUn, mr_model_part.Nodes(), mbUn);

				mr_matrix_container.FillScalarFromDatabase(DISTANCE, mphi_n1, mr_model_part.Nodes(), mbphi_n1);
				mr_matrix_container.FillOldScalarFromDatabase(DISTANCE, mphi_n, mr_model_part.Nodes(), mbphi_n);
				*/
				mr_matrix_container.FillCoordinatesFromDatabase(mx, mr_model_part.Nodes(), mbx);

				// TODO: Unused, remove
				// Set flag for first time step
				// mFirstStep = true;

				KRATOS_CATCH("")
			}

			//
			// ComputeTimeStep
			//
			// Function to set adequate time step size

			void ComputeTimeStep(double CFLNumber)
			{
				KRATOS_TRY

				// Local variable for time step size
				double delta_t = 1e10;

				// Getting value of current velocity
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mUn1, mr_model_part.Nodes(), mbUn1, true);

				// Loop over all nodes
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					// Use CFL condition to compute time step size
					double delta_t_i = CFLNumber * mr_matrix_container.GetHmin()[i_node] / Norm2_3(mUn1[i_node]);

					// Choose the overall minimum of delta_t_i
					if (delta_t_i < delta_t)
					{
						delta_t = delta_t_i;
					}
				}

				// Write time step size to Kratos
				ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
				CurrentProcessInfo[DELTA_TIME] = delta_t;

				KRATOS_CATCH("")
			}

			//
			// Solve
			//
			// Function to solve fluid equations - fractional step 1: compute fractional momentum

			void Solve()
			{
				int64_t t0 = timeNanos();

				// Read variables from Kratos
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mUn1, mr_model_part.Nodes(), mbUn1);
				mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mUn, mr_model_part.Nodes(), mbUn);

				mr_matrix_container.FillScalarFromDatabase(DISTANCE, mphi_n1, mr_model_part.Nodes(), mbphi_n1);
				mr_matrix_container.FillOldScalarFromDatabase(DISTANCE, mphi_n, mr_model_part.Nodes(), mbphi_n);

				std::cout << "Reading data from Kratos took " << timeNanos() - t0 << " ns." << std::endl;

				// Read time step size from Kratos
				ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
				double delta_t = CurrentProcessInfo[DELTA_TIME];

				// TODO: This should take place on GPU

				// Compute advective velocity - area average of the current velocity
				double coefficient = 1.00;
				CalculateAdvectiveVelocity(mbUn, mbUn1, mbA, coefficient);

				// Compute intrinsic time
				double time_inv = 1.00 / delta_t;

				// Calling Solve1 OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkSolve1, 0, mr_matrix_container.GetHminBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkSolve1, 1, mbA);
				mrDeviceGroup.SetBufferAsKernelArg(mkSolve1, 2, mbTau);
				mrDeviceGroup.SetKernelArg(mkSolve1, 3, time_inv);
				mrDeviceGroup.SetKernelArg(mkSolve1, 4, n_nodes);

				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkSolve1, n_nodes);

				mr_matrix_container.AssignVectorToVector(mbphi_n, mbWork); // mbWork = mbphi_n

				// First step of Runge Kutta
				mr_matrix_container.SetToZero(mbrhs);
				CalculateRHS(mbphi_n1, mbA, mbrhs);

				mr_matrix_container.Add_Minv_value1(mbWork, mbWork, delta_t / 6.0, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);
				mr_matrix_container.Add_Minv_value1(mbphi_n1, mbphi_n, 0.5 * delta_t, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);

				// Second step
				mr_matrix_container.SetToZero(mbrhs);
				CalculateRHS(mbphi_n1, mbA, mbrhs);

				mr_matrix_container.Add_Minv_value1(mbWork, mbWork, delta_t / 3.0, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);
				mr_matrix_container.Add_Minv_value1(mbphi_n1, mbphi_n, 0.5 * delta_t, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);

				// Third step
				CalculateAdvectiveVelocity(mbUn, mbUn1, mbA, coefficient);

				mr_matrix_container.SetToZero(mbrhs);
				CalculateRHS(mbphi_n1, mbA, mbrhs);

				mr_matrix_container.Add_Minv_value1(mbWork, mbWork, delta_t / 3.0, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);
				mr_matrix_container.Add_Minv_value1(mbphi_n1, mbphi_n, delta_t, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);

				// Fourth step
				CalculateAdvectiveVelocity(mbUn, mbUn1, mbA, coefficient);

				mr_matrix_container.SetToZero(mbrhs);

				int64_t t1 = timeNanos();

				CalculateRHS(mbphi_n1, mbA, mbrhs);

				std::cout << "CalculateRHS took " << timeNanos() - t1 << " ns." << std::endl;

				mr_matrix_container.Add_Minv_value1(mbWork, mbWork, delta_t / 6.0, mr_matrix_container.GetInvertedMassBuffer(), mbrhs);

				// Compute right-hand side
				mr_matrix_container.AssignVectorToVector(mbWork, mbphi_n1);

				int64_t t2 = timeNanos();

				mr_matrix_container.WriteScalarToDatabase(DISTANCE, mphi_n1, mr_model_part.Nodes(), mbphi_n1);

				std::cout << "WriteScalarToDatabase took " << timeNanos() - t2 << " ns." << std::endl;
			}

			//
			// CalculateRHS
			//
			// Function to calculate right-hand side of fractional momentum equation

			void CalculateRHS(cl_uint phi_buffer, cl_uint convective_velocity_buffer, cl_uint rhs_buffer)
			{
				// Calling CalculateRHS1 OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS1, 0, mbPi);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS1, 1, phi_buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS1, 2, mr_matrix_container.GetRowStartIndexBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS1, 3, mr_matrix_container.GetColumnIndexBuffer());
				//mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS1, 4, mr_matrix_container.GetEdgeValuesBuffer());
				mrDeviceGroup.SetImageAsKernelArg(mkCalculateRHS1, 4, mr_matrix_container.GetEdgeValuesBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS1, 5, mr_matrix_container.GetInvertedMassBuffer());
				mrDeviceGroup.SetKernelArg(mkCalculateRHS1, 6, n_nodes);
				mrDeviceGroup.SetLocalMemAsKernelArg(mkCalculateRHS1, 7, (mrDeviceGroup.WorkGroupSizes[mkCalculateRHS1][0] + 1) * sizeof(cl_uint));

				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkCalculateRHS1, n_nodes);

				// Calling CalculateRHS2 OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS2, 0, mbPi);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS2, 1, phi_buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS2, 2, mr_matrix_container.GetRowStartIndexBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS2, 3, mr_matrix_container.GetColumnIndexBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS2, 4, mbx);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS2, 5, mbBeta);
				mrDeviceGroup.SetKernelArg(mkCalculateRHS2, 6, n_nodes);
				mrDeviceGroup.SetLocalMemAsKernelArg(mkCalculateRHS2, 7, (mrDeviceGroup.WorkGroupSizes[mkCalculateRHS2][0] + 1) * sizeof(cl_uint));

				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkCalculateRHS2, n_nodes);

				// Calling CalculateRHS3 OpenCL kernel

				// Setting arguments
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS3, 0, mbPi);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS3, 1, phi_buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS3, 2, mr_matrix_container.GetRowStartIndexBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS3, 3, mr_matrix_container.GetColumnIndexBuffer());
				//mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS3, 4, mr_matrix_container.GetEdgeValuesBuffer());
				mrDeviceGroup.SetImageAsKernelArg(mkCalculateRHS3, 4, mr_matrix_container.GetEdgeValuesBuffer());
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS3, 5, convective_velocity_buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS3, 6, mbBeta);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS3, 7, mbrhs);
				mrDeviceGroup.SetBufferAsKernelArg(mkCalculateRHS3, 8, mbTau);
				mrDeviceGroup.SetKernelArg(mkCalculateRHS3, 9, n_nodes);
				mrDeviceGroup.SetLocalMemAsKernelArg(mkCalculateRHS3, 10, (mrDeviceGroup.WorkGroupSizes[mkCalculateRHS3][0] + 1) * sizeof(cl_uint));

				// Execute OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mkCalculateRHS3, n_nodes);
			}

			//
			// CalculateAdvectiveVelocity
			//
			// CalculateAdvectiveVelocity

			void CalculateAdvectiveVelocity(cl_uint Un_buffer, cl_uint Un1_buffer, cl_uint A_buffer, double coefficient)
			{
				// Special cases
				if (coefficient == 0.00)
				{
					// A = Un
					mrDeviceGroup.CopyBufferToBuffer(Un_buffer, A_buffer);
				}
				else if (coefficient == 1.00)
				{
					// A = Un1
					mrDeviceGroup.CopyBufferToBuffer(Un1_buffer, A_buffer);
				}
				else
				{
					// Setting arguments
					mrDeviceGroup.SetBufferAsKernelArg(mkCalculateAdvectiveVelocity, 0, Un_buffer);
					mrDeviceGroup.SetBufferAsKernelArg(mkCalculateAdvectiveVelocity, 1, Un1_buffer);
					mrDeviceGroup.SetBufferAsKernelArg(mkCalculateAdvectiveVelocity, 2, A_buffer);
					mrDeviceGroup.SetKernelArg(mkCalculateAdvectiveVelocity, 3, coefficient);
					mrDeviceGroup.SetKernelArg(mkCalculateAdvectiveVelocity, 4, n_nodes);

					// Execute OpenCL kernel
					mrDeviceGroup.ExecuteKernel(mkCalculateAdvectiveVelocity, n_nodes);
				}
			}

			//
			// Clear
			//
			// Frees allocated memory

			void Clear()
			{
				KRATOS_TRY

				// Delete OpenCL buffers and Images

				mrDeviceGroup.DeleteBuffer(mbWork);

				mrDeviceGroup.DeleteBuffer(mbPi);

				mrDeviceGroup.DeleteBuffer(mbUn);
				mrDeviceGroup.DeleteBuffer(mbUn1);

				mrDeviceGroup.DeleteBuffer(mbphi_n);
				mrDeviceGroup.DeleteBuffer(mbphi_n1);

				mrDeviceGroup.DeleteBuffer(mbA);

				mrDeviceGroup.DeleteBuffer(mbTau);
				mrDeviceGroup.DeleteBuffer(mbBeta);

				mrDeviceGroup.DeleteBuffer(mbx);

				mrDeviceGroup.DeleteBuffer(mbrhs);

				// Free memory

				FreeArray(&mWork);

				FreeArray(&mPi);

				FreeArray(&mUn);
				FreeArray(&mUn1);

				FreeArray(&mphi_n);
				FreeArray(&mphi_n1);

				FreeArray(&mA);

				FreeArray(&mTau);
				FreeArray(&mBeta);

				FreeArray(&mx);

				KRATOS_CATCH("")
			}

		private:

			// Matrix container
			OpenCLMatrixContainer &mr_matrix_container;

			// OpenCL stuff
			OpenCL::DeviceGroup &mrDeviceGroup;

			// OpenCL buffers
			cl_uint mbWork, mbPi, mbUn, mbUn1, mbphi_n, mbphi_n1, mbA, mbTau, mbBeta, mbx, mbrhs;

			// OpenCL program and kernels
			cl_uint mpOpenCLPureConvectionEdgeBased, mkSolve1, mkCalculateRHS1, mkCalculateRHS2, mkCalculateRHS3, mkCalculateAdvectiveVelocity;

			// Associated model part
			ModelPart &mr_model_part;

			// No. of nodes
			unsigned int n_nodes;

			// TODO: Unused, remove
			// bool msmooth_convective_velocity;
			// bool minclude_shock_capturing;

			// Nodal values

			// Velocity vector U at time steps n and n + 1
			CalcVectorType mUn1, mUn;
			CalcVectorType mPi;

			// Pressure vector p at time steps n and n + 1
			CalcVectorType mx;

			ValuesVectorType mBeta;
			ValuesVectorType mWork;

			// Variable to be convected
			ValuesVectorType mphi_n, mphi_n1;

			// Advective velocity vector
 			CalcVectorType mA;

			// TODO: Unused, remove
			// Flag for first time step
			//bool mFirstStep;

			// Intrinsic time step size
			ValuesVectorType mTau;

			// RHS, used in Solve()
			ValuesVectorType mrhs;

	};

} // Namespace Kratos

#endif // KRATOS_OPENCL_PURE_CONVECTION_EDGEBASED_SOLVER_H_INCLUDED defined
