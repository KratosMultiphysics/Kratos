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


// External includes


// Project includes
#include "opencl_interface.h"

namespace Kratos
{

namespace OpenCL
{

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
		LinearSolverOptimizationParameters(DeviceGroup &DeviceGroup): mrDeviceGroup(DeviceGroup)
		{
			//
		}

	private:

		DeviceGroup &mrDeviceGroup;
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
		LinearSolver(DeviceGroup &DeviceGroup, cl_int A_RowIndices_Buffer, cl_int A_Column_Indices_Buffer, cl_int A_Values_Buffer, cl_int B_Values_Buffer, unsigned int MaxIterations, double Tolerance):
			mrDeviceGroup(DeviceGroup),
			mOptimized(false),
			mMaxIterations(MaxIterations),
			mIterationNo(0),
			mTolerance(Tolerance)
        {

			// General routines
			mpOpenCLLinearSolverGeneral = mrDeviceGroup.BuildProgramFromFile("opencl_linear_solver.cl", "-cl-fast-relaxed-math -DKRATOS_OCL_GENERAL_KERNELS_ONLY");
			mkUpdateVectorWithBackup32 = mrDeviceGroup.RegisterKernel(mpOpenCLLinearSolverGeneral, "UpdateVectorWithBackup32");
        }

	private:

		DeviceGroup &mrDeviceGroup;
		bool mOptimized;
		cl_uint mpOpenCLLinearSolverGeneral;
		cl_uint mkUpdateVectorWithBackup32;
		unsigned int mMaxIterations;
		unsigned int mIterationNo;
		double mTolerance;
	};

}

}


#endif
