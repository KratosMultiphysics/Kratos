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

#define SPLIT_OSS
//#define SYMM_PRESS


// System includes
#include <string>
#include <iostream>
#include <algorithm>


// External includes


// Project includes
//#include "includes/define.h"
#include "includes/model_part.h"
//#include "includes/node.h"
//#include "geometries/geometry.h"
//#include "utilities/geometry_utilities.h"
//#include "incompressible_fluid_application.h"
//#include "opencl_interface.h"

namespace Kratos
{

	//
	// OpenCLPureConvectionEdgeBased
	//
	// OpenCL based pure convection edge based solver

	class OpenCLPureConvectionEdgeBased
	{
		public:

			//
			// Used types

			// typedef unsigned int *IndicesVectorType;
			typedef double *CalcVectorType;
			typedef double *ValuesVectorType;

			//
			// OpenCLPureConvectionEdgeBased
			//
			// Constructor

			OpenCLPureConvectionEdgeBased(OpenCLMatrixContainer opencl_matrix_container, ModelPart &model_part): opencl_matrix_container(opencl_matrix_container), mrDeviceGroup(opencl_matrix_container.GetDeviceGroup()), mr_model_part(model_part)
			{
				// Loading program
				// TODO: Add optimization flags here
				mOpenCLPureConvectionEdgeBasedProgram = mrDeviceGroup.BuildProgramFromFile("opencl_pure_convection_edgebased.cl");

				// Register kernels
				mCalculateAdvectiveVelocityKernel = mrDeviceGroup.RegisterKernel(mOpenCLPureConvectionEdgeBasedProgram, "CalculateAdvectiveVelocity");
			}

			//
			// ~OpenCLPureConvectionEdgeBased
			//
			// Destructor

			~OpenCLPureConvectionEdgeBased()
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

					// Allocate a single chunk of memory for variables
					// TODO: Use Page-locked memory for faster data transfer to GPU
					// TODO: Order variables, such that variables copied to GPU together are together here too
					// TODO: Account for device address alignment

#define KRATOS_ASSIGN_AND_ADVANCE_POINTER(P, Size)	P = Temp; Temp += Size

					AllocateArray(&Mem, 21 * n_nodes); // 6 * n_nodes + 5 * (3 * n_nodes)
					double *Temp = Mem;

					// Size data vectors
					KRATOS_ASSIGN_AND_ADVANCE_POINTER(mWork, n_nodes);

					KRATOS_ASSIGN_AND_ADVANCE_POINTER(mPi, n_nodes * 3);

					KRATOS_ASSIGN_AND_ADVANCE_POINTER(mUn, n_nodes * 3);
					KRATOS_ASSIGN_AND_ADVANCE_POINTER(mUn1, n_nodes * 3);

					KRATOS_ASSIGN_AND_ADVANCE_POINTER(mphi_n, n_nodes);
					KRATOS_ASSIGN_AND_ADVANCE_POINTER(mphi_n1, n_nodes);

					KRATOS_ASSIGN_AND_ADVANCE_POINTER(mA, n_nodes * 3);

					KRATOS_ASSIGN_AND_ADVANCE_POINTER(mHmin, n_nodes);
					KRATOS_ASSIGN_AND_ADVANCE_POINTER(mTau, n_nodes);
					KRATOS_ASSIGN_AND_ADVANCE_POINTER(mBeta, n_nodes);

					KRATOS_ASSIGN_AND_ADVANCE_POINTER(mx, n_nodes * 3);

					// Read variables from database

					// TODO: It seems that we do not need this as Solve() does this at first step
					/*
					opencl_matrix_container.FillVectorFromDatabase(VELOCITY, mUn1, mr_model_part.Nodes());
					opencl_matrix_container.FillOldVectorFromDatabase(VELOCITY, mUn, mr_model_part.Nodes());

					opencl_matrix_container.FillScalarFromDatabase(DISTANCE, mphi_n1, mr_model_part.Nodes());
					opencl_matrix_container.FillOldScalarFromDatabase(DISTANCE, mphi_n, mr_model_part.Nodes());
					*/

					opencl_matrix_container.FillCoordinatesFromDatabase(mx, mr_model_part.Nodes());

					// Set flag for first time step
					mFirstStep = true;

					mHmin = opencl_matrix_container.GetHmin();

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
				opencl_matrix_container.FillVectorFromDatabase(VELOCITY, mUn1, mr_model_part.Nodes());

				// Loop over all nodes
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					// Use CFL condition to compute time step size
					double delta_t_i = CFLNumber * mHmin[i_node] / Norm2_3(mUn1, i_node);

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
				// TODO: Correct this

				//ValuesVectorType rhs;
				//rhs.resize(n_nodes);

				// Read variables from Kratos
				opencl_matrix_container.FillVectorFromDatabase(VELOCITY, mUn1, mr_model_part.Nodes());
				opencl_matrix_container.FillOldVectorFromDatabase(VELOCITY, mUn, mr_model_part.Nodes());

				opencl_matrix_container.FillScalarFromDatabase(DISTANCE, mphi_n1, mr_model_part.Nodes());
				opencl_matrix_container.FillOldScalarFromDatabase(DISTANCE, mphi_n, mr_model_part.Nodes());

				// Read time step size from Kratos
				ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
				double delta_t = CurrentProcessInfo[DELTA_TIME];


			}

			//
			// CalculateAdvectiveVelocity
			//
			// CalculateAdvectiveVelocity

			void CalculateAdvectiveVelocity(const CalcVectorType mUn, const CalcVectorType mUn1, CalcVectorType mA, double coefficient)
			{
				// Setting arguments
				// TODO: Set arguments
				// TODO: Add a flag so that we know if we should really set them; we may need to change argument list

				// Call OpenCL kernel
				mrDeviceGroup.ExecuteKernel(mCalculateAdvectiveVelocityKernel, n_nodes);
			}

			//
			// Clear
			//
			// Frees allocated memory

			void Clear()
			{
				KRATOS_TRY

				// We only need to free this, as we allocated one chunk of memory only
				FreeArray(&Mem);

				KRATOS_CATCH("")
			}

		private:

			// Matrix container
			OpenCLMatrixContainer opencl_matrix_container;

			// OpenCL stuff
			OpenCL::DeviceGroup &mrDeviceGroup;

			cl_uint mOpenCLPureConvectionEdgeBasedProgram, mCalculateAdvectiveVelocityKernel;

			// Associated model part
			ModelPart &mr_model_part;

			// Pointer to memory used for variables
			double *Mem;

			// No. of nodes
			unsigned int n_nodes;

			bool msmooth_convective_velocity;
			bool minclude_shock_capturing;

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

			// Minimum length of the edges surrounding each nodal point
			ValuesVectorType mHmin;

			// Flag for first time step
			bool mFirstStep;

			// Intrinsic time step size
			ValuesVectorType mTau;

	};

} //namespace Kratos

#endif // KRATOS_OPENCL_PURE_CONVECTION_EDGEBASED_SOLVER_H_INCLUDED defined


