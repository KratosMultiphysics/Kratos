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
//   Date:                $Date: 2010-09-30 18:26:51 $
//   Revision:            $Revision: 1.00 $
//
//

//
// opencl_edge_data.cl
//
// OpenCL kernels and functions used in opencl_edge_data.h


#include "opencl_common.cl"


// Common kernels

//
// SetToZero
//
// Zeros a vector

__kernel void SetToZero(__global ValueType *Vector, const IndexType n)
{
	// Get work item index
	const size_t id = get_global_id(0);

	// Check if we are in the range
	if (id < n)
	{
		Vector[id] = 0.00;
	}
}

//
// Add_Minv_value1
//
// DestinationVector = Origin1Vector + Value * MinvVector * OriginVector
// double version

__kernel void Add_Minv_value1(__global ValueType *DestinationVector, __global ValueType *Origin1Vector, const ValueType Value, __global ValueType *MinvVector, __global ValueType *OriginVector, const IndexType n)
{
	// Get work item index
	const size_t id = get_global_id(0);

	// Check if we are in the range
	if (id < n)
	{
		DestinationVector[id] = Origin1Vector[id] + Value * MinvVector[id] * OriginVector[id];
	}
}

//
// Add_Minv_value3
//
// DestinationVector = Origin1Vector + Value * MinvVector * OriginVector
// double3 version

__kernel void Add_Minv_value3(__global VectorType *DestinationVector, __global VectorType *Origin1Vector, const ValueType Value, __global VectorType *MinvVector, __global VectorType *OriginVector, const IndexType n)
{
	// Get work item index
	const size_t id = get_global_id(0);

	// Check if we are in the range
	if (id < n)
	{
		DestinationVector[id] = Origin1Vector[id] + Value * MinvVector[id] * OriginVector[id];
	}
}
