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
//   Date:                $Date: 2012-03-23 01:26:13 $
//   Revision:            $Revision: 1.0 $
//
//


#include "opencl_common.cl"


#define KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE (1 << KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE_BITS)

#define KRATOS_OCL_SPMV_CSR_WORKGROUP_SIZE (1 << KRATOS_OCL_SPMV_CSR_WORKGROUP_SIZE_BITS)
#define KRATOS_OCL_SPMV_CSR_ROWS_PER_WORKGROUP (1 << KRATOS_OCL_SPMV_CSR_ROWS_PER_WORKGROUP_BITS)

#define KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE_BITS (KRATOS_OCL_SPMV_CSR_WORKGROUP_SIZE_BITS - KRATOS_OCL_SPMV_CSR_ROWS_PER_WORKGROUP_BITS)
#define KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE (1 << KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE_BITS)


#ifdef KRATOS_OCL_NEED_INNER_PROD

//
// Norm2Squared
//
// Calculates the sum of squares of a given vector
// Note: Final part must be done on host side

__kernel void __attribute__((reqd_work_group_size(KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE, 1, 1))) Norm2Squared(__global ValueType const *X_Values, __global ValueType *Y_Values, IndexType N, __local ValueType *Buffer)
{
	IndexType gid = get_global_id(0);

	// Serial part
	ValueType Accumulator = 0.00;

	while (gid < N)
	{
		Accumulator += X_Values[gid] * X_Values[gid];
		gid += get_global_size(0);
	}

	// Parallel part
	IndexType lid = get_local_id(0);
	Buffer[lid] = Accumulator;

	barrier(CLK_LOCAL_MEM_FENCE);

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 512

	if (lid < 512)
	{
		Buffer[lid] += Buffer[lid + 512];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 256

	if (lid < 256)
	{
		Buffer[lid] += Buffer[lid + 256];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 128

	if (lid < 128)
	{
		Buffer[lid] += Buffer[lid + 128];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 64

	if (lid < 64)
	{
		Buffer[lid] += Buffer[lid + 64];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 32

	if (lid < 32)
	{
		Buffer[lid] += Buffer[lid + 32];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 16

	if (lid < 16)
	{
		Buffer[lid] += Buffer[lid + 16];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 8

	if (lid < 8)
	{
		Buffer[lid] += Buffer[lid + 8];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 4

	if (lid < 4)
	{
		Buffer[lid] += Buffer[lid + 4];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 2

	if (lid < 2)
	{
		Buffer[lid] += Buffer[lid + 2];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 1

	if (lid < 1)
	{
		Buffer[lid] += Buffer[lid + 1];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

	// Store final result
	if (lid == 0)
	{
		Y_Values[get_group_id(0)] = Buffer[0];
	}
}

//
// InnerProd
//
// Calculates the inner product of two given vectors
// Note: Final part must be done on host side

__kernel void __attribute__((reqd_work_group_size(KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE, 1, 1))) InnerProd(__global ValueType const *X_Values, __global ValueType const *Y_Values, __global ValueType *Z_Values, IndexType N, __local ValueType *Buffer)
{
	IndexType gid = get_global_id(0);

	// Serial part
	ValueType Accumulator = 0.00;

	while (gid < N)
	{
		Accumulator += X_Values[gid] * Y_Values[gid];
		gid += get_global_size(0);
	}

	// Parallel part
	IndexType lid = get_local_id(0);
	Buffer[lid] = Accumulator;

	barrier(CLK_LOCAL_MEM_FENCE);

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 512

	if (lid < 512)
	{
		Buffer[lid] += Buffer[lid + 512];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 256

	if (lid < 256)
	{
		Buffer[lid] += Buffer[lid + 256];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 128

	if (lid < 128)
	{
		Buffer[lid] += Buffer[lid + 128];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 64

	if (lid < 64)
	{
		Buffer[lid] += Buffer[lid + 64];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 32

	if (lid < 32)
	{
		Buffer[lid] += Buffer[lid + 32];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 16

	if (lid < 16)
	{
		Buffer[lid] += Buffer[lid + 16];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 8

	if (lid < 8)
	{
		Buffer[lid] += Buffer[lid + 8];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 4

	if (lid < 4)
	{
		Buffer[lid] += Buffer[lid + 4];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 2

	if (lid < 2)
	{
		Buffer[lid] += Buffer[lid + 2];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 1

	if (lid < 1)
	{
		Buffer[lid] += Buffer[lid + 1];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

	// Store final result
	if (lid == 0)
	{
		Z_Values[get_group_id(0)] = Buffer[0];
	}
}

//
// InnerProd2
//
// Calculates the inner product of two given vectors and the first vector with itself
// Note: Final part must be done on host side

__kernel void __attribute__((reqd_work_group_size(KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE, 1, 1))) InnerProd2(__global ValueType const *X_Values, __global ValueType const *Y_Values, __global ValueType *Z_Values1, __global ValueType *Z_Values2, IndexType N, __local ValueType *Buffer1, __local ValueType *Buffer2)
{
	IndexType gid = get_global_id(0);

	// Serial part
	ValueType Accumulator1 = 0.00;
	ValueType Accumulator2 = 0.00;

	while (gid < N)
	{
		Accumulator1 += X_Values[gid] * X_Values[gid];
		Accumulator2 += X_Values[gid] * Y_Values[gid];
		gid += get_global_size(0);
	}

	// Parallel part
	IndexType lid = get_local_id(0);
	Buffer1[lid] = Accumulator1;
	Buffer2[lid] = Accumulator2;

	barrier(CLK_LOCAL_MEM_FENCE);

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 512

	if (lid < 512)
	{
		Buffer1[lid] += Buffer1[lid + 512];
		Buffer2[lid] += Buffer2[lid + 512];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 256

	if (lid < 256)
	{
		Buffer1[lid] += Buffer1[lid + 256];
		Buffer2[lid] += Buffer2[lid + 256];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 128

	if (lid < 128)
	{
		Buffer1[lid] += Buffer1[lid + 128];
		Buffer2[lid] += Buffer2[lid + 128];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 64

	if (lid < 64)
	{
		Buffer1[lid] += Buffer1[lid + 64];
		Buffer2[lid] += Buffer2[lid + 64];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 32

	if (lid < 32)
	{
		Buffer1[lid] += Buffer1[lid + 32];
		Buffer2[lid] += Buffer2[lid + 32];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 16

	if (lid < 16)
	{
		Buffer1[lid] += Buffer1[lid + 16];
		Buffer2[lid] += Buffer2[lid + 16];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 8

	if (lid < 8)
	{
		Buffer1[lid] += Buffer1[lid + 8];
		Buffer2[lid] += Buffer2[lid + 8];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 4

	if (lid < 4)
	{
		Buffer1[lid] += Buffer1[lid + 4];
		Buffer2[lid] += Buffer2[lid + 4];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 2

	if (lid < 2)
	{
		Buffer1[lid] += Buffer1[lid + 2];
		Buffer2[lid] += Buffer2[lid + 2];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if KRATOS_OCL_INNER_PROD_WORKGROUP_SIZE > 1

	if (lid < 1)
	{
		Buffer1[lid] += Buffer1[lid + 1];
		Buffer2[lid] += Buffer2[lid + 1];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

	// Store final result
	if (lid == 0)
	{
		size_t pid = get_group_id(0);

		Z_Values1[pid] = Buffer1[0];
		Z_Values2[pid] = Buffer2[0];
	}
}


#endif  // KRATOS_OCL_NEED_INNER_PROD

#ifdef KRATOS_OCL_NEED_SPMV_CSR


//
// SpMV_CSR
//
// Performs sparse matrix vector product in CSR format

__kernel void __attribute__((reqd_work_group_size(KRATOS_OCL_SPMV_CSR_WORKGROUP_SIZE, 1, 1))) SpMV_CSR(__global IndexType const *A_RowIndices, __global IndexType const *A_ColumnIndices, __global ValueType const *A_Values, __global ValueType const *X_Values, __global ValueType *Y_Values, IndexType N, __local IndexType *Bounds, __local ValueType *Buffer)
{
	const IndexType gid = get_group_id(0);
	const IndexType tid = get_local_id(0);

	const IndexType lgid = tid >> KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE_BITS;
	const IndexType ltid = tid & (KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE - 1);

	const IndexType Row = (gid << KRATOS_OCL_SPMV_CSR_ROWS_PER_WORKGROUP_BITS) + lgid;

	if (Row < N)
	{
		// TODO: There is a BUG here, fix it!

/*		// Read bounds
		if (tid == KRATOS_OCL_SPMV_CSR_WORKGROUP_SIZE - 1)
		{
			Bounds[KRATOS_OCL_SPMV_CSR_ROWS_PER_WORKGROUP] = A_RowIndices[Row + 1];
		}

		if (ltid == 0)
		{
			Bounds[lgid] = A_RowIndices[Row];
		}
*/
		// Zero reduction buffer
		Buffer[tid] = 0.00;
	}

	// __local memory barrier
	barrier(CLK_LOCAL_MEM_FENCE);

	if (Row < N)
	{
		// TODO: It seems that this is making Bank Conflict!

		// Read bounds from __local memory
		const IndexType Start = A_RowIndices[Row];//Bounds[lgid];
		const IndexType End = A_RowIndices[Row + 1];//Bounds[lgid + 1];

		// TODO: Use images!

		// Actual multiplication
		for (IndexType i = Start + ltid; i < End; i += KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE)
		{
			Buffer[tid] += A_Values[i] * X_Values[A_ColumnIndices[i]];
		}
	}

	// __local memory barrier
	barrier(CLK_LOCAL_MEM_FENCE);

	// Reduction of results

#if KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE > 512

	if (Row < N)
	{
		if (ltid < 512)
		{
			Buffer[tid] += Buffer[tid + 512];
		}
	}

#ifdef KRATOS_OCL_SPMV_CSR_USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE > 256

	if (Row < N)
	{
		if (ltid < 256)
		{
			Buffer[tid] += Buffer[tid + 256];
		}
	}

#ifdef KRATOS_OCL_SPMV_CSR_USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE > 128

	if (Row < N)
	{
		if (ltid < 128)
		{
			Buffer[tid] += Buffer[tid + 128];
		}
	}

#ifdef KRATOS_OCL_SPMV_CSR_USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE > 64

	if (Row < N)
	{
		if (ltid < 64)
		{
			Buffer[tid] += Buffer[tid + 64];
		}
	}

#ifdef KRATOS_OCL_SPMV_CSR_USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE > 32

	if (Row < N)
	{
		if (ltid < 32)
		{
			Buffer[tid] += Buffer[tid + 32];
		}
	}

#ifdef KRATOS_OCL_SPMV_CSR_USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE > 16

	if (Row < N)
	{
		if (ltid < 16)
		{
			Buffer[tid] += Buffer[tid + 16];
		}
	}

#ifdef KRATOS_OCL_SPMV_CSR_USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE > 8

	if (Row < N)
	{
		if (ltid < 8)
		{
			Buffer[tid] += Buffer[tid + 8];
		}
	}

#ifdef KRATOS_OCL_SPMV_CSR_USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE > 4

	if (Row < N)
	{
		if (ltid < 4)
		{
			Buffer[tid] += Buffer[tid + 4];
		}
	}

#ifdef KRATOS_OCL_SPMV_CSR_USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE > 2

	if (Row < N)
	{
		if (ltid < 2)
		{
			Buffer[tid] += Buffer[tid + 2];
		}
	}

#ifdef KRATOS_OCL_SPMV_CSR_USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if KRATOS_OCL_SPMV_CSR_LOCAL_WORKGROUP_SIZE > 1

	if (Row < N)
	{
		if (ltid < 1)
		{
			Buffer[tid] += Buffer[tid + 1];
		}
	}

#ifdef KRATOS_OCL_SPMV_CSR_USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

	if (Row < N)
	{
		// Store final result
		if (ltid == 0)
		{
			Y_Values[Row] = Buffer[tid];
		}
	}
}


#endif  // KRATOS_OCL_NEED_SPMV_CSR

#ifdef KRATOS_OCL_NEED_GENERIC_KERNELS


//
// ZeroVector3Negate
//
// Zeros three vectors and updates the fourth with negative of another
// Note: x = 0; y = 0; z = 0; t = -u

__kernel void ZeroVector3Negate(__global ValueType *X_Values, __global ValueType *Y_Values, __global ValueType *Z_Values, __global ValueType *T_Values, __global const ValueType *U_Values, IndexType N)
{
	// Get work item index
	const size_t gid = get_global_id(0);

	// Check if we are in the range
	if (gid < N)
	{
		X_Values[gid] = 0.00;
		Y_Values[gid] = 0.00;
		Z_Values[gid] = 0.00;
		T_Values[gid] = -U_Values[gid];
	}
}

//
// ZeroVectorCopy2
//
// Zeros a vector and updates the two others with another
// Note: x = 0; y = u; z = u

__kernel void ZeroVectorCopy2(__global ValueType *X_Values, __global ValueType *Y_Values, __global ValueType *Z_Values, __global const ValueType *U_Values, IndexType N)
{
	// Get work item index
	const size_t gid = get_global_id(0);

	// Check if we are in the range
	if (gid < N)
	{
		X_Values[gid] = 0.00;
		Y_Values[gid] = U_Values[gid];
		Z_Values[gid] = U_Values[gid];
	}
}

//
// ZeroVector2Copy2
//
// Zeros two vectors and updates the two others with another
// Note: x = 0; y = 0; z = u; t = u

__kernel void ZeroVector2Copy2(__global ValueType *X_Values, __global ValueType *Y_Values, __global ValueType *Z_Values, __global ValueType *T_Values, __global const ValueType *U_Values, IndexType N)
{
	// Get work item index
	const size_t gid = get_global_id(0);

	// Check if we are in the range
	if (gid < N)
	{
		X_Values[gid] = 0.00;
		Y_Values[gid] = 0.00;
		Z_Values[gid] = U_Values[gid];
		T_Values[gid] = U_Values[gid];
	}
}

//
// UpdateVector
//
// Updates a vector
// Note: p = r + Beta * p

__kernel void UpdateVector(__global ValueType *P_Values, __global const ValueType *R_Values, ValueType Beta, IndexType N)
{
	// Get work item index
	const size_t gid = get_global_id(0);

	// Check if we are in the range
	if (gid < N)
	{
		P_Values[gid] = R_Values[gid] + Beta * P_Values[gid];
	}
}

//
// UpdateVector2
//
// Updates two vectors
// Note: x = x + Alpha * p; r = r - Alpha * q

__kernel void UpdateVector2(__global ValueType *X_Values, __global ValueType *R_Values, __global const ValueType *P_Values, __global const ValueType *Q_Values, ValueType Alpha, IndexType N)
{
	// Get work item index
	const size_t gid = get_global_id(0);

	// Check if we are in the range
	if (gid < N)
	{
		X_Values[gid] += Alpha * P_Values[gid];
		R_Values[gid] -= Alpha * Q_Values[gid];
	}
}

//
// UpdateVector3WithBackup2
//
// Updates two vectors with 3 vectors after backing them up in others
// Note: t1 = x1; x1 = a1 * x1 + b1 * y1 + c1 * z1; y1 = t1; t2 = x2; x2 = a2 * x2 + b2 * y2 + c2 * z2; y2 = t2

__kernel void UpdateVector3WithBackup2(__global ValueType *X_Values1, __global ValueType *Y_Values1, __global const ValueType *Z_Values1, ValueType A1, ValueType B1, ValueType C1, __global ValueType *X_Values2, __global ValueType *Y_Values2, __global const ValueType *Z_Values2, ValueType A2, ValueType B2, ValueType C2, IndexType N)
{
	// Get work item index
	const size_t gid = get_global_id(0);

	// Check if we are in the range
	if (gid < N)
	{
		ValueType T;

		T = X_Values1[gid];
		X_Values1[gid] = A1 * X_Values1[gid] + B1 * Y_Values1[gid] + C1 * Z_Values1[gid];
		Y_Values1[gid] = T;

		T = X_Values2[gid];
		X_Values2[gid] = A2 * X_Values2[gid] + B2 * Y_Values2[gid] + C2 * Z_Values2[gid];
		Y_Values2[gid] = T;
	}
}

//
// UpdateVector4
//
// Updates four vectors
// Note: p = r + Beta * p; q = s + Beta * q; x = x + Alpha * p; r = r - Alpha * q;

__kernel void UpdateVector4(__global ValueType *X_Values, __global ValueType *P_Values, __global ValueType *Q_Values, __global ValueType *R_Values, __global const ValueType *S_Values, ValueType Alpha, ValueType Beta, IndexType N)
{
	// Get work item index
	const size_t gid = get_global_id(0);

	// Check if we are in the range
	if (gid < N)
	{
		ValueType TP = R_Values[gid] + Beta * P_Values[gid];
		ValueType TQ = S_Values[gid] + Beta * Q_Values[gid];

		X_Values[gid] += Alpha * TP;
		R_Values[gid] -= Alpha * TQ;

		P_Values[gid] = TP;
		Q_Values[gid] = TQ;
	}
}

//
// MinimalKernel
//
// A minimal kernel to find out the wavefront size

__kernel void MinimalKernel()
{
	// Nothing to do!
}

#endif  // KRATOS_OCL_NEED_GENERIC_KERNELS
