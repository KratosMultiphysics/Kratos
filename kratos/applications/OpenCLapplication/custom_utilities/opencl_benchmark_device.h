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

#define KRATOS_OCL_BENCHMARK_CHUNKS		1536
#define KRATOS_OCL_BENCHMARK_CHUNKSIZE	4096
#define KRATOS_OCL_BENCHMARK_ELEMENTS	(KRATOS_OCL_BENCHMARK_CHUNKS * KRATOS_OCL_BENCHMARK_CHUNKSIZE)
#define KRATOS_OCL_BENCHMARK_TRIES		50

namespace Kratos
{

namespace OpenCL
{
	//
	// BenchmarkDevice
	//
	// Benchmarks an OpenCL device

	void BenchmarkDevice(DeviceGroup &Device)
	{
		std::cout <<
			std::endl <<
			"OpenCL benchmark device" << std::endl <<
			std::endl;

		cl_uint Buffer1 = Device.CreateBuffer(KRATOS_OCL_BENCHMARK_ELEMENTS * sizeof(cl_float4), CL_MEM_READ_ONLY);
		cl_uint Buffer2 = Device.CreateBuffer(KRATOS_OCL_BENCHMARK_ELEMENTS * sizeof(cl_float4), CL_MEM_READ_ONLY);
		cl_uint Buffer3 = Device.CreateBuffer(KRATOS_OCL_BENCHMARK_ELEMENTS * sizeof(cl_float4), CL_MEM_WRITE_ONLY);

		cl_uint Program = Device.BuildProgramFromFile("opencl_benchmark_device.cl", "-cl-fast-relaxed-math -DKRATOS_OCL_BENCHMARK_CHUNKSIZE=" KRATOS_OCL_STRINGIZE(KRATOS_OCL_BENCHMARK_CHUNKSIZE));

		cl_uint InitializeBuffersKernel = Device.RegisterKernel(Program, "InitializeBuffers");
		cl_uint BenchmarkDeviceKernel1 = Device.RegisterKernel(Program, "BenchmarkDevice1");
		cl_uint BenchmarkDeviceKernel2 = Device.RegisterKernel(Program, "BenchmarkDevice2");

		std::cout <<
			"Workgroup size of the benchmark kernel #1 is: " << Device.WorkGroupSizes[BenchmarkDeviceKernel1][0] << std::endl<<
			std::endl <<
			"Workgroup size of the benchmark kernel #2 is: " << Device.WorkGroupSizes[BenchmarkDeviceKernel1][0] << std::endl<<
			std::endl;

		// Set kernel arguments

		Device.SetBufferAsKernelArg(InitializeBuffersKernel, 0, Buffer1);
		Device.SetBufferAsKernelArg(InitializeBuffersKernel, 1, Buffer2);
		Device.SetBufferAsKernelArg(InitializeBuffersKernel, 2, Buffer3);
		Device.SetKernelArg(InitializeBuffersKernel, 3, KRATOS_OCL_BENCHMARK_ELEMENTS);

		Device.SetBufferAsKernelArg(BenchmarkDeviceKernel1, 0, Buffer1);
		Device.SetBufferAsKernelArg(BenchmarkDeviceKernel1, 1, Buffer2);
		Device.SetBufferAsKernelArg(BenchmarkDeviceKernel1, 2, Buffer3);
		Device.SetKernelArg(BenchmarkDeviceKernel1, 3, KRATOS_OCL_BENCHMARK_CHUNKS);

		Device.SetBufferAsKernelArg(BenchmarkDeviceKernel2, 0, Buffer1);
		Device.SetBufferAsKernelArg(BenchmarkDeviceKernel2, 1, Buffer2);
		Device.SetBufferAsKernelArg(BenchmarkDeviceKernel2, 2, Buffer3);
		Device.SetKernelArg(BenchmarkDeviceKernel2, 3, KRATOS_OCL_BENCHMARK_ELEMENTS);

		unsigned int FloatElemets = sizeof(cl_float4) / sizeof(cl_float) * KRATOS_OCL_BENCHMARK_ELEMENTS;
		float *HostBuffer = new float[FloatElemets];

		int64_t t1, t2, T1, T2, T3;

		T1 = 0;
		T2 = 0;
		T3 = 0;

		for (unsigned int i = 0; i < KRATOS_OCL_BENCHMARK_TRIES; i++)
		{
			// Initialize buffers

			Device.ExecuteKernel(InitializeBuffersKernel, KRATOS_OCL_BENCHMARK_ELEMENTS);

			// Perform benchmark

			t1 = Timer();

			Device.ExecuteKernel(BenchmarkDeviceKernel1, KRATOS_OCL_BENCHMARK_CHUNKS * Device.WorkGroupSizes[BenchmarkDeviceKernel1][0]);
			Device.Synchronize();

			t2 = Timer();

			if ((i == 0) || (t2 - t1 < T1))
			{
				T1 = t2 - t1;
			}
		}

		Device.CopyBuffer(Buffer3, OpenCL::DeviceToHost, VoidPList(1, HostBuffer));

		std::cout <<
			"Test #1:" << std::endl <<
			std::endl;

		for (unsigned int i = 0; i < FloatElemets; i++)
		{
			if (HostBuffer[i] != KRATOS_OCL_BENCHMARK_ELEMENTS)
			{
				std::cout <<
					"Solution invalid at element " << i << ", value: " << HostBuffer[i] << std::endl <<
					std::endl;

				break;
			}
		}

		std::cout <<
			"Memory bandwidth measured:" <<
			std::endl <<
			"Add:	" << 3 * KRATOS_OCL_BENCHMARK_ELEMENTS * sizeof(cl_float4) / static_cast <double> (T1) << " GB/s" << std::endl <<
			3 * KRATOS_OCL_BENCHMARK_ELEMENTS * sizeof(cl_float4) << " bytes of data processed in " << T1 << "ns." << std::endl <<
			std::endl;

		for (unsigned int i = 0; i < KRATOS_OCL_BENCHMARK_TRIES; i++)
		{
			// Initialize buffers

			Device.ExecuteKernel(InitializeBuffersKernel, KRATOS_OCL_BENCHMARK_ELEMENTS);

			// Perform benchmark

			t1 = Timer();

			Device.ExecuteKernel(BenchmarkDeviceKernel2, KRATOS_OCL_BENCHMARK_ELEMENTS);
			Device.Synchronize();

			t2 = Timer();

			if ((i == 0) || (t2 - t1 < T2))
			{
				T2 = t2 - t1;
			}
		}

		Device.CopyBuffer(Buffer3, OpenCL::DeviceToHost, VoidPList(1, HostBuffer));

		std::cout <<
			"Test #2:" << std::endl <<
			std::endl;

		for (unsigned int i = 0; i < FloatElemets; i++)
		{
			if (HostBuffer[i] != KRATOS_OCL_BENCHMARK_ELEMENTS)
			{
				std::cout <<
					"Solution invalid at element " << i << ", value: " << HostBuffer[i] << std::endl <<
					std::endl;

				break;
			}
		}

		std::cout <<
			"Memory bandwidth measured:" <<
			std::endl <<
			"Add:	" << 3 * KRATOS_OCL_BENCHMARK_ELEMENTS * sizeof(cl_float4) / static_cast <double> (T2) << " GB/s" << std::endl <<
			3 * KRATOS_OCL_BENCHMARK_ELEMENTS * sizeof(cl_float4) << " bytes of data processed in " << T2 << "ns." << std::endl <<
			std::endl;

		delete[] HostBuffer;

		std::cout <<
			"Test #3:" << std::endl <<
			std::endl;

		for (unsigned int i = 0; i < KRATOS_OCL_BENCHMARK_TRIES; i++)
		{
			t1 = Timer();

			Device.CopyBufferToBuffer(Buffer1, Buffer3);
			Device.Synchronize();

			t2 = Timer();

			if ((i == 0) || (t2 - t1 < T3))
			{
				T3 = t2 - t1;
			}
		}

		std::cout <<
			"Memory bandwidth measured:" <<
			std::endl <<
			"Copy:	" << 2 * KRATOS_OCL_BENCHMARK_ELEMENTS * sizeof(cl_float4) / static_cast <double> (T3) << " GB/s" << std::endl <<
			2 * KRATOS_OCL_BENCHMARK_ELEMENTS * sizeof(cl_float4) << " bytes of data processed in " << T3 << "ns." << std::endl <<
			std::endl;
	}
}

}

#endif  // KRATOS_OPENCL_BENCHMARK_DEVICE_H_INCLUDED
