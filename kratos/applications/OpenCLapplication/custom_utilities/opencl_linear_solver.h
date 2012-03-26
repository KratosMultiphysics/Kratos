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
//   Date:                $Date: 2012-03-23 00:46:43 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_OPENCL_LINEAR_SOLVER_H_INCLUDED)
#define KRATOS_OPENCL_LINEAR_SOLVER_H_INCLUDED


// System includes
#include <cmath>
#include <sstream>


// External includes


// Project includes
#include "opencl_interface.h"


// InnerProd and SpMV kernel parameters range
#define KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE_BITS_MIN 0
#define KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE_BITS_MAX 10

#define KRATOS_OCL_SPMV_CSR_ROWS_PER_WORKGROUP_BITS_MIN 0
#define KRATOS_OCL_SPMV_CSR_ROWS_PER_WORKGROUP_BITS_MAX 10

#define KRATOS_OCL_SPMV_CSR_WORKGROUP_SIZE_BITS_MIN 0
#define KRATOS_OCL_SPMV_CSR_WORKGROUP_SIZE_BITS_MAX 10

// Kernel optimization iteration count
#define KRATOS_OCL_OPTIMIZATION_ITERATION_COUNT_1 10
#define KRATOS_OCL_OPTIMIZATION_ITERATION_COUNT_2 10


namespace Kratos
{

namespace OpenCL
{
	//
	// Timer
	//
	// Gives current system time in ns

	int64_t Timer()
	{
		struct timespec tp;

		clock_gettime(CLOCK_MONOTONIC, &tp);

		return (unsigned long long) tp.tv_sec * (1000ULL * 1000ULL * 1000ULL) + (unsigned long long) tp.tv_nsec;
	}

	//
	// LinearSolverOptimizationParameters
	//
	// A class to hold optimized parameters for a linear solver

	class LinearSolverOptimizationParameters
	{
	public:

		//
		// LinearSolverOptimizationParameters
		//
		// Constructor

		LinearSolverOptimizationParameters(DeviceGroup &DeviceGroup, cl_uint Size): mrDeviceGroup(DeviceGroup), mSize(Size)
		{
			// Nothing to do!
		}

		//
		// ~LinearSolverOptimizationParameters
		//
		// Destructor

		~LinearSolverOptimizationParameters()
		{
			// Nothing to do!
		}

		//
		// OptimizeInnerProd
		//
		// Selects the best InnerProd kernel for given vectors
		// Note: Use only once, before using the class with a LinearSolver

		void OptimizeInnerProd(cl_uint X_Values_Buffer, cl_uint Y_Values_Buffer, cl_uint Z_Values_Buffer, cl_uint WavefrontSize)
		{
			cl_uint mpOpenCLLinearSolver, mkInnerProd;
			int64_t T0, T1, BestTime;

			BestTime = 0x7FFFFFFFFFFFFFFF;

			for (unsigned int i = KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE_BITS_MIN; i < KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE_BITS_MAX; i++)
			{
				std::stringstream OptionsBuilder;

				OptionsBuilder <<
					"-cl-fast-relaxed-math" << " " <<
					"-DKRATOS_OCL_NEED_INNER_PROD" << " " <<
					"-DKRATOS_OCL_INNER_PROD_WORKGROUP_SIZE_BITS=" << i;

				mpOpenCLLinearSolver = mrDeviceGroup.BuildProgramFromFile("opencl_linear_solver.cl", OptionsBuilder.str().c_str());

				// Register kernel
				mkInnerProd = mrDeviceGroup.RegisterKernel(mpOpenCLLinearSolver, "InnerProd");

				//X_Values, Y_Values, Z_Values, N, __local ValueType *Buffer)
				mrDeviceGroup.SetBufferAsKernelArg(mkInnerProd, 0, X_Values_Buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkInnerProd, 1, Y_Values_Buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkInnerProd, 2, Z_Values_Buffer);
				mrDeviceGroup.SetKernelArg(mkInnerProd, 3, mSize);
				mrDeviceGroup.SetLocalMemAsKernelArg(mkInnerProd, 4, (1 << i) * sizeof(cl_double));

				// Warm-up
				for (int k = 0; k < KRATOS_OCL_OPTIMIZATION_ITERATION_COUNT_1; k++)
				{
					mrDeviceGroup.ExecuteKernel(mkInnerProd, mSize);
				}

				mrDeviceGroup.Synchronize();

				T1 = 0;

				// Start timing
				for (int k = 0; k < KRATOS_OCL_OPTIMIZATION_ITERATION_COUNT_2; k++)
				{
					T0 = Timer();

					mrDeviceGroup.ExecuteKernel(mkInnerProd, mSize);
					mrDeviceGroup.Synchronize();

					T1 += Timer() - T0;
				}

				if (T1 < BestTime)
				{
					BestTime = T1;

					mOptimizedInnerProdKernel = mkInnerProd;
					mOptimizedInnerProd2Kernel = mrDeviceGroup.RegisterKernel(mpOpenCLLinearSolver, "InnerProd2");  // With same parameters

					mOptimizedInnerProdKernelLaunchSize = mSize;
					mOptimizedInnerProdKernelBufferSize1 = ((mSize + (1 << i) - 1) / (1 << i)) * sizeof(cl_double);
					mOptimizedInnerProdKernelBufferSize2 = (1 << i) * sizeof(cl_double);
				}
			}
		}

		//
		// OptimizeSpMV
		//
		// Selects the best SpMV kernel for given system of equations
		// Note: Use only once, before using the class with a LinearSolver; fills X with A * B

		void OptimizeSpMV(cl_uint A_RowIndices_Buffer, cl_uint A_Column_Indices_Buffer, cl_uint A_Values_Buffer, cl_uint B_Values_Buffer, cl_uint X_Values_Buffer, cl_uint WavefrontSize)
		{
			cl_uint mpOpenCLLinearSolver, mkSpMVCSR;
			int64_t T0, T1, BestTime;

			BestTime = 0x7FFFFFFFFFFFFFFF;

			for (unsigned int i = KRATOS_OCL_SPMV_CSR_ROWS_PER_WORKGROUP_BITS_MIN; i < KRATOS_OCL_SPMV_CSR_ROWS_PER_WORKGROUP_BITS_MAX; i++)
			{
				for (unsigned int j = KRATOS_OCL_SPMV_CSR_WORKGROUP_SIZE_BITS_MIN; j < KRATOS_OCL_SPMV_CSR_WORKGROUP_SIZE_BITS_MAX; j++)
				{
					if (j >= i)
					{
						std::stringstream OptionsBuilder;

						OptionsBuilder <<
							"-cl-fast-relaxed-math" << " " <<
							"-DKRATOS_OCL_NEED_SPMV_CSR" << " " <<
							"-DKRATOS_OCL_SPMV_CSR_ROWS_PER_WORKGROUP_BITS=" << i << " " <<
							"-DKRATOS_OCL_SPMV_CSR_WORKGROUP_SIZE_BITS=" << j;

						// Check if a local mem barrier is needed
						if (1U << (j - i) > WavefrontSize)
						{
							OptionsBuilder <<
								" " <<
								"-DKRATOS_OCL_SPMV_CSR_USE_LOCAL_MEM_BARRIER";
						}

						mpOpenCLLinearSolver = mrDeviceGroup.BuildProgramFromFile("opencl_linear_solver.cl", OptionsBuilder.str().c_str());

						// Register kernel
						mkSpMVCSR = mrDeviceGroup.RegisterKernel(mpOpenCLLinearSolver, "SpMV_CSR");

						mrDeviceGroup.SetBufferAsKernelArg(mkSpMVCSR, 0, A_RowIndices_Buffer);
						mrDeviceGroup.SetBufferAsKernelArg(mkSpMVCSR, 1, A_Column_Indices_Buffer);
						mrDeviceGroup.SetBufferAsKernelArg(mkSpMVCSR, 2, A_Values_Buffer);
						mrDeviceGroup.SetBufferAsKernelArg(mkSpMVCSR, 3, B_Values_Buffer);
						mrDeviceGroup.SetBufferAsKernelArg(mkSpMVCSR, 4, X_Values_Buffer);
						mrDeviceGroup.SetKernelArg(mkSpMVCSR, 5, mSize);
						mrDeviceGroup.SetLocalMemAsKernelArg(mkSpMVCSR, 6, ((1 << i) + 1) * sizeof(cl_uint));  // Note: Change this when IndexType is changed in opencl_common.cl
						mrDeviceGroup.SetLocalMemAsKernelArg(mkSpMVCSR, 7, (1 << j) * sizeof(cl_double));

						// Warm-up
						for (int k = 0; k < KRATOS_OCL_OPTIMIZATION_ITERATION_COUNT_1; k++)
						{
							mrDeviceGroup.ExecuteKernel(mkSpMVCSR, mSize * (1 << (j - i)) + 1);
						}

						mrDeviceGroup.Synchronize();

						T1 = 0;

						// Start timing
						for (int k = 0; k < KRATOS_OCL_OPTIMIZATION_ITERATION_COUNT_2; k++)
						{
							T0 = Timer();

							mrDeviceGroup.ExecuteKernel(mkSpMVCSR, mSize * (1 << (j - i)) + 1);
							mrDeviceGroup.Synchronize();

							T1 += Timer() - T0;
						}

						if (T1 < BestTime)
						{
							BestTime = T1;

							mOptimizedSpMVKernel = mkSpMVCSR;

							mOptimizedSpMVKernelLaunchSize = mSize * (1 << (j - i)) + 1;
							mOptimizedSpMVKernelBufferSize1 = ((1 << i) + 1) * sizeof(cl_uint);  // Note: Change this when IndexType is changed in opencl_common.cl
							mOptimizedSpMVKernelBufferSize2 = (1 << j) * sizeof(cl_double);
						}
					}
				}
			}
		}

		//
		// GetOptimizedSpMVKernel
		//
		// Returns optimized SpMV kernel

		cl_uint &GetOptimizedSpMVKernel()
		{
			return mOptimizedSpMVKernel;
		}

		//
		// GetOptimizedSpMVKernelLaunchSize
		//
		// Returns optimized SpMV kernel launch size

		cl_uint &GetOptimizedSpMVKernelLaunchSize()
		{
			return mOptimizedSpMVKernelLaunchSize;
		}

		//
		// GetOptimizedSpMVKernelBufferSize1
		//
		// Returns optimized SpMV kernel first buffer size

		cl_uint &GetOptimizedSpMVKernelBufferSize1()
		{
			return mOptimizedSpMVKernelBufferSize1;
		}

		//
		// GetOptimizedSpMVKernelBufferSize2
		//
		// Returns optimized SpMV kernel second buffer size

		cl_uint &GetOptimizedSpMVKernelBufferSize2()
		{
			return mOptimizedSpMVKernelBufferSize2;
		}

		//
		// GetOptimizedInnerProdKernel
		//
		// Returns optimized inner product kernel

		cl_uint &GetOptimizedInnerProdKernel()
		{
			return mOptimizedInnerProdKernel;
		}

		//
		// GetOptimizedInnerProdKernelLaunchSize
		//
		// Returns optimized inner product kernel launch size

		cl_uint &GetOptimizedInnerProdKernelLaunchSize()
		{
			return mOptimizedInnerProdKernelLaunchSize;
		}

		//
		// GetOptimizedInnerProd2Kernel
		//
		// Returns optimized inner product 2 kernel

		cl_uint &GetOptimizedInnerProd2Kernel()
		{
			return mOptimizedInnerProd2Kernel;
		}

		//
		// GetOptimizedInnerProdKernelBufferSize1
		//
		// Returns optimized inner product kernel first buffer size

		cl_uint &GetOptimizedInnerProdKernelBufferSize1()
		{
			return mOptimizedInnerProdKernelBufferSize1;
		}

		//
		// GetOptimizedInnerProdKernelBufferSize2
		//
		// Returns optimized inner product kernel second buffer size

		cl_uint &GetOptimizedInnerProdKernelBufferSize2()
		{
			return mOptimizedInnerProdKernelBufferSize2;
		}

	private:

		DeviceGroup &mrDeviceGroup;
		cl_uint mSize;
		cl_uint mOptimizedSpMVKernel, mOptimizedInnerProdKernel, mOptimizedInnerProd2Kernel;
		cl_uint mOptimizedSpMVKernelLaunchSize, mOptimizedSpMVKernelBufferSize1, mOptimizedSpMVKernelBufferSize2, mOptimizedInnerProdKernelLaunchSize, mOptimizedInnerProdKernelBufferSize1, mOptimizedInnerProdKernelBufferSize2;
	};

	//
	// LinearSolver
	//
	// A class to solve linear systems of equations on OpenCL devices

	class LinearSolver
	{
	public:

		//
		// LinearSolver
		//
		// Constructor

		LinearSolver(DeviceGroup &DeviceGroup, LinearSolverOptimizationParameters &OptimizationParameters, cl_uint Size, unsigned int MaxIterations, double Tolerance):
			mrDeviceGroup(DeviceGroup),
			mOptimizationParameters(OptimizationParameters),
			mSize(Size),
			mMaxIterations(MaxIterations),
			mTolerance(Tolerance),
			mIterationNo(0),
			mAchievedTolerance(0.00)
        {
			// General routines
			mpOpenCLLinearSolverGeneral = mrDeviceGroup.BuildProgramFromFile("opencl_linear_solver.cl", "-cl-fast-relaxed-math -DKRATOS_OCL_GENERAL_KERNELS_ONLY");
			mkUpdateVectorWithBackup32 = mrDeviceGroup.RegisterKernel(mpOpenCLLinearSolverGeneral, "UpdateVectorWithBackup32");
			mkZeroVector2Negate = mrDeviceGroup.RegisterKernel(mpOpenCLLinearSolverGeneral, "ZeroVector3Negate");

			// Temporary vectors needed on GPU and CPU
			mbr = mrDeviceGroup.CreateBuffer(mSize * sizeof(cl_double), CL_MEM_READ_WRITE);
			mbAr = mrDeviceGroup.CreateBuffer(mSize * sizeof(cl_double), CL_MEM_READ_WRITE);
			mbx_old = mrDeviceGroup.CreateBuffer(mSize * sizeof(cl_double), CL_MEM_READ_WRITE);
			mbr_old = mrDeviceGroup.CreateBuffer(mSize * sizeof(cl_double), CL_MEM_READ_WRITE);

			mbReductionBuffer1 = mrDeviceGroup.CreateBuffer(mOptimizationParameters.GetOptimizedInnerProdKernelBufferSize1() * sizeof(cl_double), CL_MEM_READ_WRITE);
			mbReductionBuffer2 = mrDeviceGroup.CreateBuffer(mOptimizationParameters.GetOptimizedInnerProdKernelBufferSize1() * sizeof(cl_double), CL_MEM_READ_WRITE);  // Yes, both are of the same size!

			mReductionBuffer1 = new double[mOptimizationParameters.GetOptimizedInnerProdKernelBufferSize1()];
			mReductionBuffer2 = new double[mOptimizationParameters.GetOptimizedInnerProdKernelBufferSize1()];  // Yes, both are of the same size!
        }

        //
        // ~LinearSolver
        //
        // Destructor

        ~LinearSolver()
        {
			mrDeviceGroup.DeleteBuffer(mbr);
			mrDeviceGroup.DeleteBuffer(mbAr);
			mrDeviceGroup.DeleteBuffer(mbx_old);
			mrDeviceGroup.DeleteBuffer(mbr_old);

			delete [] mReductionBuffer1;
			delete [] mReductionBuffer2;
        }

        //
        // Solve
        //
        // Solves the linear system

        bool Solve(cl_uint A_RowIndices_Buffer, cl_uint A_Column_Indices_Buffer, cl_uint A_Values_Buffer, cl_uint B_Values_Buffer, cl_uint X_Values_Buffer)
		{
			// Initialization

			// x = x_old = r_old = 0.00, r = -b
			mrDeviceGroup.SetBufferAsKernelArg(mkZeroVector2Negate, 0, X_Values_Buffer);
			mrDeviceGroup.SetBufferAsKernelArg(mkZeroVector2Negate, 1, mbx_old);
			mrDeviceGroup.SetBufferAsKernelArg(mkZeroVector2Negate, 2, mbr_old);
			mrDeviceGroup.SetBufferAsKernelArg(mkZeroVector2Negate, 3, mbr);
			mrDeviceGroup.SetBufferAsKernelArg(mkZeroVector2Negate, 4, B_Values_Buffer);
			mrDeviceGroup.SetKernelArg(mkZeroVector2Negate, 5, mSize);

			mrDeviceGroup.ExecuteKernel(mkZeroVector2Negate, mSize);


			mIterationNo = 0;
			double Rho = 1.00, Rho_old, Gamma, Gamma_old, rr, rr_old, rAr;

			while (true)
			{
				// Ar = A * r
				mrDeviceGroup.SetBufferAsKernelArg(mOptimizationParameters.GetOptimizedSpMVKernel(), 0, A_RowIndices_Buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mOptimizationParameters.GetOptimizedSpMVKernel(), 1, A_Column_Indices_Buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mOptimizationParameters.GetOptimizedSpMVKernel(), 2, A_Values_Buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mOptimizationParameters.GetOptimizedSpMVKernel(), 3, mbr);
				mrDeviceGroup.SetBufferAsKernelArg(mOptimizationParameters.GetOptimizedSpMVKernel(), 4, mbAr);
				mrDeviceGroup.SetKernelArg(mOptimizationParameters.GetOptimizedSpMVKernel(), 5, mSize);
				mrDeviceGroup.SetLocalMemAsKernelArg(mOptimizationParameters.GetOptimizedSpMVKernel(), 6, mOptimizationParameters.GetOptimizedSpMVKernelBufferSize1());
				mrDeviceGroup.SetLocalMemAsKernelArg(mOptimizationParameters.GetOptimizedSpMVKernel(), 7, mOptimizationParameters.GetOptimizedSpMVKernelBufferSize2());

				mrDeviceGroup.ExecuteKernel(mOptimizationParameters.GetOptimizedSpMVKernel(), mOptimizationParameters.GetOptimizedSpMVKernelLaunchSize());


				// rr = r.r, rAr = r.Ar

				// Phase 1 on GPU
				mrDeviceGroup.SetBufferAsKernelArg(mOptimizationParameters.GetOptimizedInnerProd2Kernel(), 0, mbr);
				mrDeviceGroup.SetBufferAsKernelArg(mOptimizationParameters.GetOptimizedInnerProd2Kernel(), 1, mbAr);
				mrDeviceGroup.SetBufferAsKernelArg(mOptimizationParameters.GetOptimizedInnerProd2Kernel(), 2, mbReductionBuffer1);
				mrDeviceGroup.SetBufferAsKernelArg(mOptimizationParameters.GetOptimizedInnerProd2Kernel(), 3, mbReductionBuffer2);
				mrDeviceGroup.SetKernelArg(mOptimizationParameters.GetOptimizedInnerProd2Kernel(), 4, mSize);
				mrDeviceGroup.SetLocalMemAsKernelArg(mOptimizationParameters.GetOptimizedInnerProd2Kernel(), 5, mOptimizationParameters.GetOptimizedInnerProdKernelBufferSize2());
				mrDeviceGroup.SetLocalMemAsKernelArg(mOptimizationParameters.GetOptimizedInnerProd2Kernel(), 6, mOptimizationParameters.GetOptimizedInnerProdKernelBufferSize2());  // Yes, both are of the same size!

				mrDeviceGroup.ExecuteKernel(mOptimizationParameters.GetOptimizedInnerProd2Kernel(), mOptimizationParameters.GetOptimizedInnerProdKernelLaunchSize());


				// Phase 2 on CPU
				mrDeviceGroup.CopyBuffer(mbReductionBuffer1, DeviceToHost, VoidPList(1, mReductionBuffer1));
				mrDeviceGroup.CopyBuffer(mbReductionBuffer2, DeviceToHost, VoidPList(1, mReductionBuffer2));

				rr = 0.00;
				rAr = 0.00;

				#pragma omp parallel for reduction(+:rr)
				for (unsigned int i = 0; i < mOptimizationParameters.GetOptimizedInnerProdKernelBufferSize1(); i++)
				{
					rr += mReductionBuffer1[i];
				}

				#pragma omp parallel for reduction(+:rAr)
				for (unsigned int i = 0; i < mOptimizationParameters.GetOptimizedInnerProdKernelBufferSize1(); i++)
				{
					rAr += mReductionBuffer2[i];
				}

				Gamma = rr / rAr;
				mAchievedTolerance = sqrt(rr);

				// Convergence check
				if (rr < mTolerance * mTolerance)
				{
					return true;
				}
				else
				{
					if (mIterationNo > mMaxIterations)
					{
						return false;
					}
				}

				// If not first iteration, Rho = 1 / (1 - (rr * Gamma) / (rr_old * Gamma_old * Rho_old))
				if (mIterationNo++)
				{
					Rho = 1.00 / (1.00 - (rr * Gamma) / (rr_old * Gamma_old * Rho_old));
				}

				// Update vectors

				// x = x_old * (1 - Rho) + (x - r * Gamma) * Rho
				// r = r_old * (1 - Rho) + (r - Ar * Gamma) * Rho

				mrDeviceGroup.SetBufferAsKernelArg(mkUpdateVectorWithBackup32, 0, X_Values_Buffer);
				mrDeviceGroup.SetBufferAsKernelArg(mkUpdateVectorWithBackup32, 1, mbx_old);
				mrDeviceGroup.SetBufferAsKernelArg(mkUpdateVectorWithBackup32, 2, mbr);
				mrDeviceGroup.SetKernelArg(mkUpdateVectorWithBackup32, 3, Rho);
				mrDeviceGroup.SetKernelArg(mkUpdateVectorWithBackup32, 4, 1.00 - Rho);
				mrDeviceGroup.SetKernelArg(mkUpdateVectorWithBackup32, 5, Rho * Gamma);
				mrDeviceGroup.SetBufferAsKernelArg(mkUpdateVectorWithBackup32, 6, mbr);
				mrDeviceGroup.SetBufferAsKernelArg(mkUpdateVectorWithBackup32, 7, mbr_old);
				mrDeviceGroup.SetBufferAsKernelArg(mkUpdateVectorWithBackup32, 8, mbAr);
				mrDeviceGroup.SetKernelArg(mkUpdateVectorWithBackup32, 9, Rho);
				mrDeviceGroup.SetKernelArg(mkUpdateVectorWithBackup32, 10, 1.00 - Rho);
				mrDeviceGroup.SetKernelArg(mkUpdateVectorWithBackup32, 11, Rho * Gamma);
				mrDeviceGroup.SetKernelArg(mkUpdateVectorWithBackup32, 12, mSize);

				mrDeviceGroup.ExecuteKernel(mkUpdateVectorWithBackup32, mSize);


				// Prepare for next iteration
				rr_old = rr;
				Rho_old = Rho;
				Gamma_old = Gamma;
			}
		}

	private:

		DeviceGroup &mrDeviceGroup;
		LinearSolverOptimizationParameters mOptimizationParameters;
		cl_uint mpOpenCLLinearSolverGeneral;
		cl_uint mkUpdateVectorWithBackup32, mkZeroVector2Negate;
		cl_uint mSize;
		cl_uint mbr, mbAr, mbx_old, mbr_old, mbReductionBuffer1, mbReductionBuffer2;
		unsigned int mMaxIterations;
		double mTolerance;
		unsigned int mIterationNo;
		double mAchievedTolerance;
		double *mReductionBuffer1, *mReductionBuffer2;
	};

}

}


#endif  // KRATOS_OPENCL_LINEAR_SOLVER_H_INCLUDED
