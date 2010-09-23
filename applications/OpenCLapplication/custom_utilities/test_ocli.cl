#pragma OPENCL EXTENSION cl_amd_fp64: enable

__kernel void Test(__global double *input, __global double *output, double offset)
{
	size_t id = get_global_id(0);
	output[id] = input[id] * input[id] + offset;
}


__kernel void Test1(__global double *x, __global double *y)
{
	int z;

	z += 1;
}

__kernel void Test2(__global double *x, __global double *y)
{
	int z;

	z += 1;
}
