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
//   Date:                $Date: 2011-07-04 15:01:01 $
//   Revision:            $Revision: 1.00 $
//
//

//
// opencl_benchmark_device.cl
//
// OpenCL kernels for benchmarking an OpenCL device

#include "opencl_common.cl"


__kernel void InitializeBuffers(__global float4 *Buffer1, __global float4 *Buffer2, __global float4 *Buffer3, const IndexType Elements)
{
	// Get work item index
	const size_t i = get_global_id(0);

	// Check if we are in the range
	if (i < Elements)
	{
		Buffer1[i] = i;
		Buffer2[i] = Elements - i;
		Buffer3[i] = 0.00;
	}
}

__kernel void BenchmarkDevice(__global float4 *Buffer1, __global float4 *Buffer2, __global float4 *Buffer3, const IndexType Chunks)
{
	// Get work item index
	const size_t gid = get_group_id(0);
	const size_t lid = get_local_id(0);

	const size_t stride = get_local_size(0);

	// Check if we are in the range
	if (get_global_id(0) < Chunks)
	{
		uint start = gid * KRATOS_OCL_BENCHMARK_CHUNKSIZE;
		uint end = start + KRATOS_OCL_BENCHMARK_CHUNKSIZE;

		for (uint j = start + lid; j < end; j += stride)
		{
			Buffer3[j] = Buffer1[j] + Buffer2[j];
		}
	}
}
