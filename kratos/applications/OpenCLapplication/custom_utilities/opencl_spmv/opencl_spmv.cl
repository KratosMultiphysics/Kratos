#define ROWS_PER_WORKGROUP_BITS 4
#define ROWS_PER_WORKGROUP (1 << ROWS_PER_WORKGROUP_BITS)

#define WORKGROUP_SIZE_BITS 9
#define WORKGROUP_SIZE (1 << WORKGROUP_SIZE_BITS)

#define LOCAL_WORKGROUP_SIZE_BITS (WORKGROUP_SIZE_BITS - ROWS_PER_WORKGROUP_BITS)
#define LOCAL_WORKGROUP_SIZE (1 << LOCAL_WORKGROUP_SIZE_BITS)

// This should not be needed if LOCAL_WORKGROUP_SIZE matches Warp / Wavefront size of GPU (NVIDIA: 32; AMD: 64); Defining reduces performance!
//#define USE_LOCAL_MEM_BARRIER

// See below comment!
//#include "opencl_common.cl"
#include "opencl_enable_fp64.cl"

// TODO: This is just for using Kratos' CompressedMatrix, otherwise no need to use 64 bit indices! [ulong]
typedef ulong IndexType;
typedef double ValueType;

__kernel void CSR_Matrix_Vector_Multiply(__global IndexType *A_RowIndices, __global IndexType *A_ColumnIndices, __global ValueType *A_Values, __global ValueType *X_Values, __global ValueType *Y_Values, IndexType N, __local IndexType *Bounds, __local ValueType *Buffer)
{
	const IndexType gid = get_group_id(0);
	const IndexType tid = get_local_id(0);

	const IndexType lgid = tid >> LOCAL_WORKGROUP_SIZE_BITS;
	const IndexType ltid = tid & (LOCAL_WORKGROUP_SIZE - 1);

	const IndexType Row = (gid << ROWS_PER_WORKGROUP_BITS) + lgid;
	const IndexType stride = LOCAL_WORKGROUP_SIZE;

	if (Row < N)
	{

		// Read bounds
		if (tid == WORKGROUP_SIZE - 1)
		{
			Bounds[ROWS_PER_WORKGROUP] = A_RowIndices[Row + 1];
		}

		if (ltid == 0)
		{
			Bounds[lgid] = A_RowIndices[Row];
		}

		// Zero reduction buffer
		Buffer[tid] = 0.00;

		// __local memory barrier
		barrier(CLK_LOCAL_MEM_FENCE);

		// Read bounds from __local memory
		const IndexType Start = Bounds[lgid];
		const IndexType End = Bounds[lgid + 1];

		// TODO: Use images!

		// Actual multiplication
		for (IndexType i = Start + ltid; i < End; i += LOCAL_WORKGROUP_SIZE)
		{
			Buffer[tid] += A_Values[i] * X_Values[A_ColumnIndices[i]];
		}

		// Reduction of results

#if LOCAL_WORKGROUP_SIZE > 32

		if (ltid < 32)
		{
			Buffer[tid] += Buffer[tid + 32];
		}

#ifdef USE_LOCAL_MEM_BARRIER

		barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if LOCAL_WORKGROUP_SIZE > 16

		if (ltid < 16)
		{
			Buffer[tid] += Buffer[tid + 16];
		}

#ifdef USE_LOCAL_MEM_BARRIER

		barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if LOCAL_WORKGROUP_SIZE > 8

		if (ltid < 8)
		{
			Buffer[tid] += Buffer[tid + 8];
		}

#ifdef USE_LOCAL_MEM_BARRIER

		barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if LOCAL_WORKGROUP_SIZE > 4

		if (ltid < 4)
		{
			Buffer[tid] += Buffer[tid + 4];
		}

#ifdef USE_LOCAL_MEM_BARRIER

		barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if LOCAL_WORKGROUP_SIZE > 2

		if (ltid < 2)
		{
			Buffer[tid] += Buffer[tid + 2];
		}

#ifdef USE_LOCAL_MEM_BARRIER

		barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if LOCAL_WORKGROUP_SIZE > 1

		if (ltid < 1)
		{
			Buffer[tid] += Buffer[tid + 1];
		}

#ifdef USE_LOCAL_MEM_BARRIER

		barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

		// Store final result
		if (ltid == 0)
		{
			Y_Values[Row] = Buffer[tid];
		}
	}
}
