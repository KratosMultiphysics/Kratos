//#define ROWS_PER_WORKGROUP_BITS 1
//#define ROWS_PER_WORKGROUP (1 << ROWS_PER_WORKGROUP_BITS)

#define ROWS_PER_WORKGROUP 16

// TODO: Fix this!
#define LOCAL_WORKGROUP_SIZE (512 / 16)

//#include "opencl_common.cl"
#include "opencl_enable_fp64.cl"

typedef ulong IndexType;
typedef double ValueType;

//#define KRATOS_OCL_DEBUG


#ifdef KRATOS_OCL_DEBUG

#pragma OPENCL EXTENSION cl_amd_printf: enable

#define KRATOS_OCL_VIEW_INT(x)		printf(#x ": %d\n", x)
#define KRATOS_OCL_VIEW_DOUBLE(x)	printf(#x ": %f\n", x)

#else

#define KRATOS_OCL_VIEW_INT(x)
#define KRATOS_OCL_VIEW_DOUBLE(x)

#endif

__kernel void CSR_Matrix_Vector_Multiply(__global IndexType *A_RowIndices, __global IndexType *A_ColumnIndices, __global ValueType *A_Values, __global ValueType *X_Values, __global ValueType *Y_Values, IndexType N, __local IndexType *Bounds, __local ValueType *Buffer)
{
	const IndexType gid = get_group_id(0);  // OK
	const IndexType tid = get_local_id(0);  // OK

	const IndexType workgroup_size = get_local_size(0);  // OK
	//const IndexType local_workgroup_size = workgroup_size / ROWS_PER_WORKGROUP;

	const IndexType lgid = tid / (workgroup_size / ROWS_PER_WORKGROUP);  // OK
	const IndexType ltid = tid % (workgroup_size / ROWS_PER_WORKGROUP);  // OK

	const IndexType Row = gid * ROWS_PER_WORKGROUP + lgid;
	const IndexType stride = workgroup_size / ROWS_PER_WORKGROUP; // OK

	if (Row < N)
	{

		// Read bounds

		if (tid == workgroup_size - 1)
		{
			Bounds[ROWS_PER_WORKGROUP] = A_RowIndices[Row + 1];
		}

		if (ltid == 0)
		{
			Bounds[lgid] = A_RowIndices[Row];
		}

		Buffer[tid] = 0.00;

		barrier(CLK_LOCAL_MEM_FENCE);

		const IndexType Start = Bounds[lgid];
		const IndexType End = Bounds[lgid + 1];

		//Safe but non-optimized way
		//const IndexType Start = A_RowIndices[gid * ROWS_PER_WORKGROUP + lgid];
		//const IndexType End = A_RowIndices[gid * ROWS_PER_WORKGROUP + lgid + 1];

		//

		for (IndexType i = Start + ltid; i < End; i += stride)
		{
			//
			Buffer[tid] += A_Values[i] * X_Values[A_ColumnIndices[i]];
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
			Buffer[tid] += Buffer[tid + 2];
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
			Y_Values[Row] = Buffer[tid];
		}
	}
}
