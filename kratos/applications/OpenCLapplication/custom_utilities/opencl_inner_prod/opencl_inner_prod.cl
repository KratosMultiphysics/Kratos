#define WORKGROUP_SIZE_BITS 8
#define WORKGROUP_SIZE (1 << WORKGROUP_SIZE_BITS)

// This should not be needed if WORKGROUP_SIZE matches Warp / Wavefront size of GPU (NVIDIA: 32; AMD: 64); Defining reduces performance!
#define USE_LOCAL_MEM_BARRIER

// See below comment!
//#include "opencl_common.cl"
#include "opencl_enable_fp64.cl"

// TODO: This is just for using Kratos' CompressedMatrix, otherwise no need to use 64 bit indices! [ulong]
typedef ulong IndexType;
typedef double ValueType;

#define INNER_PROD_ELEM(N) (X_Values[N] * Y_Values[N])

__kernel void Vector_Vector_Multiply(__global ValueType const *X_Values, __global ValueType const *Y_Values, __global ValueType *Z_Values, IndexType N, __local ValueType *Buffer)
{
	IndexType gid = get_global_id(0);

	// Serial part
	ValueType Accumulator = 0.00;
	while (gid < N)
	{
		Accumulator += INNER_PROD_ELEM(gid);
		gid += get_global_size(0);
	}

	// Parallel part
	IndexType lid = get_local_id(0);
	Buffer[lid] = Accumulator;

	barrier(CLK_LOCAL_MEM_FENCE);

#if WORKGROUP_SIZE > 512

	if (lid < 512)
	{
		Buffer[lid] += Buffer[lid + 512];
	}

#ifdef USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if WORKGROUP_SIZE > 256

	if (lid < 256)
	{
		Buffer[lid] += Buffer[lid + 256];
	}

#ifdef USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if WORKGROUP_SIZE > 128

	if (lid < 128)
	{
		Buffer[lid] += Buffer[lid + 128];
	}

#ifdef USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if WORKGROUP_SIZE > 64

	if (lid < 64)
	{
		Buffer[lid] += Buffer[lid + 64];
	}

#ifdef USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if WORKGROUP_SIZE > 32

	if (lid < 32)
	{
		Buffer[lid] += Buffer[lid + 32];
	}

#ifdef USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if WORKGROUP_SIZE > 16

	if (lid < 16)
	{
		Buffer[lid] += Buffer[lid + 16];
	}

#ifdef USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if WORKGROUP_SIZE > 8

	if (lid < 8)
	{
		Buffer[lid] += Buffer[lid + 8];
	}

#ifdef USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if WORKGROUP_SIZE > 4

	if (lid < 4)
	{
		Buffer[lid] += Buffer[lid + 4];
	}

#ifdef USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if WORKGROUP_SIZE > 2

	if (lid < 2)
	{
		Buffer[lid] += Buffer[lid + 2];
	}

#ifdef USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

#if WORKGROUP_SIZE > 1

	if (lid < 1)
	{
		Buffer[lid] += Buffer[lid + 1];
	}

#ifdef USE_LOCAL_MEM_BARRIER

	barrier(CLK_LOCAL_MEM_FENCE);

#endif

#endif

	// Store final result
	if (lid == 0)
	{
		Z_Values[get_group_id(0)] = Buffer[0];
	}
}
