#include "opencl_interface.h"

cl_int Test()
{
	return CL_INVALID_EVENT;
}

int main(int argc, char *argv[])
{
	//KRATOS_OCL_CHECK(CL_INVALID_EVENT);
	//KRATOS_OCL_CHECKED_EXPRESSION(Test());

	std::cout << "Initializing OpenCL device group..." << std::endl;

	Kratos::OpenCL::DeviceGroup OCLDeviceGroup(CL_DEVICE_TYPE_ALL);

	std::cout << "Found " << OCLDeviceGroup.DeviceNo << " device(s)." << std::endl;
	for (int i = 0; i < OCLDeviceGroup.DeviceNo; i++)
	{
		std::cout << "  Device " << i << ": ";
		switch (OCLDeviceGroup.DeviceTypes[i])
		{
			case CL_DEVICE_TYPE_CPU:

				std::cout << "CPU" << std::endl;
				break;

			case CL_DEVICE_TYPE_GPU:

				std::cout << "GPU" << std::endl;
				break;

			case CL_DEVICE_TYPE_ACCELERATOR:

				std::cout << "Accelerator" << std::endl;
				break;

			default:

				std::cout << "Unknown" << std::endl;
		}
	}

	std::cout << "Loading kernel(s) from file test_ocli.cl..." << std::endl;
	OCLDeviceGroup.BuildProgramFromFile("test_ocli.cl");

	std::cout << "Registering kernel Test() and setting arguments..." << std::endl;
	OCLDeviceGroup.RegisterKernel("Test");

	const int DATA_SIZE = 10;

	double inputData[DATA_SIZE]={1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	double results[DATA_SIZE]={0};

	cl_int Err;
	cl_mem input, output;

	// create buffers for the input and ouput
	input = clCreateBuffer(OCLDeviceGroup.Contexts[0], CL_MEM_READ_ONLY, sizeof(double) * DATA_SIZE, NULL, &Err);
	KRATOS_OCL_CHECK(Err);

	output = clCreateBuffer(OCLDeviceGroup.Contexts[0], CL_MEM_WRITE_ONLY, sizeof(double) * DATA_SIZE, NULL, &Err);
	KRATOS_OCL_CHECK(Err);

	// load data into the input buffer
	Err = clEnqueueWriteBuffer(OCLDeviceGroup.CommandQueues[0], input, CL_TRUE, 0, sizeof(double) * DATA_SIZE, inputData, 0, NULL, NULL);
	KRATOS_OCL_CHECK(Err);

	OCLDeviceGroup.SetKernelArg(0, 0, input);
	OCLDeviceGroup.SetKernelArg(0, 1, output);

    std::cout << "Executing kernel..." << std::endl;
	OCLDeviceGroup.ExecuteKernel(0, DATA_SIZE);

   	// copy the results from out of the output buffer
	Err = clEnqueueReadBuffer(OCLDeviceGroup.CommandQueues[0], output, CL_TRUE, 0, sizeof(double) *DATA_SIZE, results, 0, NULL, NULL);
	KRATOS_OCL_CHECK(Err);

    std::cout << "Results:" << std::endl;
    for (int i = 0; i < DATA_SIZE; i++)
    {
        std::cout << results[i] << std::endl;
    }

	return 0;
}
