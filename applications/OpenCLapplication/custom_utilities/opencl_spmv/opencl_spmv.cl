#define WORKGROUP_SIZE 256

#define ROWS_PER_WORKGROUP_BITS 4
#define ROWS_PER_WORKGROUP (1 << ROWS_PER_WORKGROUP_BITS)

#define LOCAL_WORKGROUP_SIZE (WORKGROUP_SIZE >> ROWS_PER_WORKGROUP_BITS)

#include "opencl_common.cl"

__kernel void CSR_Matrix_Vector_Multiply(__global IndexType *A_RowIndices, __global IndexType *A_ColumnIndices, __global ValueType *A_Values, __global ValueType *X_Values, __global ValueType *Y_Values, __local IndexType *Bounds, __local ValueType *Buffer)
{
	const size_t gid = get_group_id(0);  // OK
	const size_t tid = get_local_id(0);  // OK

	const size_t lgid = tid / (WORKGROUP_SIZE >> ROWS_PER_WORKGROUP_BITS);  // OK
	const size_t ltid = tid % (WORKGROUP_SIZE >> ROWS_PER_WORKGROUP_BITS);  // OK

	const size_t stride = WORKGROUP_SIZE >> ROWS_PER_WORKGROUP_BITS; // OK

	// Read bounds

	if (tid == 0)
	{
		Bounds[0] = A_RowIndices[gid << ROWS_PER_WORKGROUP_BITS];
	}

	if (tid < ROWS_PER_WORKGROUP)
	{
		Bounds[tid + 1] = A_RowIndices[(gid << ROWS_PER_WORKGROUP_BITS) + tid + 1];
		Buffer[tid] = 0.00;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	const size_t Start = Bounds[lgid];
	const size_t End = Bounds[lgid + 1];

	//

	for (IndexType i = Start + ltid; i < End; i += stride)
	{
		//
		Buffer[lgid] += A_Values[i] * X_Values[A_ColumnIndices[i]];
	}

	//

#if LOCAL_WORKGROUP_SIZE > 32

	if (ltid < 32)
	{
		Buffer[tid] += Buffer[tid + 32];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if LOCAL_WORKGROUP_SIZE > 16

	if (ltid < 16)
	{
		Buffer[tid] += Buffer[tid + 16];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if LOCAL_WORKGROUP_SIZE > 8

	if (ltid < 8)
	{
		Buffer[tid] += Buffer[tid + 8];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if LOCAL_WORKGROUP_SIZE > 4

	if (ltid < 4)
	{
		Buffer[tid] += Buffer[tid + 4];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if LOCAL_WORKGROUP_SIZE > 2

	if (ltid < 2)
	{
		Buffer[tid] += Buffer[tid + 32];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#if LOCAL_WORKGROUP_SIZE > 1

	if (ltid < 1)
	{
		Buffer[tid] += Buffer[tid + 1];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

	//

	if (ltid == 0)
	{
		Y_Values[(gid << ROWS_PER_WORKGROUP_BITS) + lgid] = Buffer[tid];
	}
}
