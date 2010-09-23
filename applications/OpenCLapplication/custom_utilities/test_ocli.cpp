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

	Kratos::OpenCL::DeviceGroup OCLDeviceGroup(CL_DEVICE_TYPE_ALL, true);

	std::cout << "Found " << OCLDeviceGroup.DeviceNo << " device(s)." << std::endl;
	for (cl_uint i = 0; i < OCLDeviceGroup.DeviceNo; i++)
	{
		std::cout << "  Device " << i << ": " << Kratos::OpenCL::DeviceTypeString(OCLDeviceGroup.DeviceTypes[i]) << std::endl;
	}

	std::cout << "Loading kernel(s) from file test_ocli.cl..." << std::endl;
	OCLDeviceGroup.BuildProgramFromFile("test_ocli.cl");

	std::cout << "Registering kernel Test()..." << std::endl;
	OCLDeviceGroup.RegisterKernel("Test");

	const int DataSize = 10;

	double Data[DataSize] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	double Results[DataSize] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

	std::cout << "Creating the buffers..." << std::endl;
	OCLDeviceGroup.CreateBuffer(DataSize * sizeof(double), CL_MEM_READ_ONLY);
	OCLDeviceGroup.CreateBuffer(DataSize * sizeof(double), CL_MEM_WRITE_ONLY);

	std::cout << "Copying the data to the device(s)..." << std::endl;
	OCLDeviceGroup.CopyBuffer(0, 0, Kratos::OpenCL::HostToDevice, Data);

	std::cout << "Setting kernel arguments..." << std::endl;
	OCLDeviceGroup.SetBufferAsKernelArg(0, 0, 0);
	OCLDeviceGroup.SetBufferAsKernelArg(0, 1, 1);
	OCLDeviceGroup.SetKernelArg(0, 2, 1000.00);

	std::cout << "Executing kernel..." << std::endl;
	OCLDeviceGroup.ExecuteKernel(0, DataSize);

	std::cout << "Copying the data from the device(s)..." << std::endl;
	OCLDeviceGroup.CopyBuffer(0, 1, Kratos::OpenCL::DeviceToHost, Results);

	std::cout << "Results:" << std::endl;
	for (int i = 0; i < DataSize; i++)
	{
		std::cout << Results[i] << "  ";
	}

	std::cout << std::endl;

	return 0;
}
