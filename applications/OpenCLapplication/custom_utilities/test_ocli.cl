#pragma OPENCL EXTENSION cl_amd_fp64: enable

/*
__kernel void Test(__global double *input, __global double *output)
{
	size_t id = get_global_id(0);
	output[id] = input[id] * input[id];
}
*/

__kernel void Test(__global double *x, __global double *y)
{
	y[0] = 2 * x[0];
}
