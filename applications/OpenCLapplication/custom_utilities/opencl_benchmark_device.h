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
//   Date:                $Date: 2011-07-04 14:45:31 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_OPENCL_BENCHMARK_DEVICE_H_INCLUDED)
#define KRATOS_OPENCL_BENCHMARK_DEVICE_H_INCLUDED

// System includes

#include "iostream"


// External includes


// Project includes

#include "opencl_interface.h"


//
// Benchmarking constants

#define KRATOS_OCL_BENCHMARK_CHUNKS		10000
#define KRATOS_OCL_BENCHMARK_CHUNKSIZE	512
#define KRATOS_OCL_BENCHMARK_ELEMENTS	(KRATOS_OCL_BENCHMARK_CHUNKS * KRATOS_OCL_BENCHMARK_CHUNKSIZE)
#define KRATOS_OCL_BENCHMARK_TRIES		10

namespace Kratos
{

namespace OpenCL
{

	//
	// Timer
	//
	// Returns system timer in nano-seconds

	int64_t Timer()
	{
		struct timespec tp;

		clock_gettime(CLOCK_MONOTONIC, &tp);

		return (unsigned long long) tp.tv_sec * (1000ULL * 1000ULL * 1000ULL) + (unsigned long long) tp.tv_nsec;
	}

	//
	// BenchmarkDevice
	//
	// Benchmarks an OpenCL device

	void BenchmarkDevice(DeviceGroup &Device)
	{
		cl_uint Buffer1 = Device.CreateBuffer(KRATOS_OCL_BENCHMARK_ELEMENTS * sizeof(cl_float4), CL_MEM_READ_ONLY);
		cl_uint Buffer2 = Device.CreateBuffer(KRATOS_OCL_BENCHMARK_ELEMENTS * sizeof(cl_float4), CL_MEM_READ_ONLY);
		cl_uint Buffer3 = Device.CreateBuffer(KRATOS_OCL_BENCHMARK_ELEMENTS * sizeof(cl_float4), CL_MEM_WRITE_ONLY);

		cl_uint Program = Device.BuildProgramFromFile("opencl_benchmark_device.cl", "-cl-fast-relaxed-math -DKRATOS_OCL_BENCHMARK_CHUNKSIZE=" KRATOS_OCL_STRINGIZE(KRATOS_OCL_BENCHMARK_CHUNKSIZE));

		cl_uint InitializeBuffersKernel = Device.RegisterKernel(Program, "InitializeBuffers");
		cl_uint BenchmarkDeviceKernel = Device.RegisterKernel(Program, "BenchmarkDevice");

		// Set kernel arguments

		Device.SetBufferAsKernelArg(InitializeBuffersKernel, 0, Buffer1);
		Device.SetBufferAsKernelArg(InitializeBuffersKernel, 1, Buffer2);
		Device.SetBufferAsKernelArg(InitializeBuffersKernel, 2, Buffer3);
		Device.SetKernelArg(InitializeBuffersKernel, 3, KRATOS_OCL_BENCHMARK_ELEMENTS);

		Device.SetBufferAsKernelArg(BenchmarkDeviceKernel, 0, Buffer1);
		Device.SetBufferAsKernelArg(BenchmarkDeviceKernel, 1, Buffer2);
		Device.SetBufferAsKernelArg(BenchmarkDeviceKernel, 2, Buffer3);
		Device.SetKernelArg(BenchmarkDeviceKernel, 3, KRATOS_OCL_BENCHMARK_CHUNKS);

		int64_t t0, t1, t2;

		t0 = 0;

		for (unsigned int i = 0; i < KRATOS_OCL_BENCHMARK_TRIES; i++)
		{
			// Initialize buffers

			Device.ExecuteKernel(InitializeBuffersKernel, KRATOS_OCL_BENCHMARK_ELEMENTS);

			// Perform benchmark

			t1 = Timer();

			Device.ExecuteKernel(BenchmarkDeviceKernel, KRATOS_OCL_BENCHMARK_CHUNKS);

			t2 = Timer();

			if ((i == 0) || (t2 - t1 < t0))
			{
				t0 = t2 - t1;
			}
		}

		std::cout <<
			"OpenCL benchmark device" << std::endl <<
			std::endl <<
			"Memory bandwidth measured:" <<
			std::endl <<
			"Add:	" << KRATOS_OCL_BENCHMARK_ELEMENTS * sizeof(cl_float4) / static_cast <double> (t0) << " GB/s" << std::endl;
	}
}

}

#endif  // KRATOS_OPENCL_BENCHMARK_DEVICE_H_INCLUDED
