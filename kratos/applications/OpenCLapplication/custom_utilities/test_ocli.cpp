#include "opencl_interface.h"
#include <cmath>

const cl_uint DataSize = 3600000;
const double Eps = 1e-10;

double Input[DataSize], Output[DataSize];

int main(int argc, char *argv[])
{
	std::cout << "Initializing OpenCL device group..." << std::endl;
	Kratos::OpenCL::DeviceGroup OCLDeviceGroup(CL_DEVICE_TYPE_ALL, false);  // Try to find all available devices

	std::cout << "Found " << OCLDeviceGroup.DeviceNo << " device(s)." << std::endl;
	for (cl_uint i = 0; i < OCLDeviceGroup.DeviceNo; i++)
	{
		std::cout << "  Device " << i << ": " << Kratos::OpenCL::DeviceTypeString(OCLDeviceGroup.DeviceTypes[i]) << std::endl;
	}

	std::cout << "Loading kernel(s) from file test_ocli.cl..." << std::endl;
	OCLDeviceGroup.BuildProgramFromFile("test_ocli.cl");

	std::cout << "Registering kernel Test()..." << std::endl;
	OCLDeviceGroup.RegisterKernel("Test");

	for (cl_uint i = 0; i < DataSize; i++)
	{
		Input[i] = i;
		Output[i] = -1.00;
	}

	std::cout << "Creating the buffers..." << std::endl;
	OCLDeviceGroup.CreateBuffer(DataSize / OCLDeviceGroup.DeviceNo * sizeof(double), CL_MEM_READ_ONLY);
	OCLDeviceGroup.CreateBuffer(DataSize / OCLDeviceGroup.DeviceNo * sizeof(double), CL_MEM_WRITE_ONLY);

	std::cout << "Copying the data to the device(s)..." << std::endl;
	for (cl_uint i = 0; i < OCLDeviceGroup.DeviceNo; i++)
	{
		OCLDeviceGroup.CopyBuffer(i, 0, Kratos::OpenCL::HostToDevice, Input + i * DataSize / OCLDeviceGroup.DeviceNo);
	}

	std::cout << "Setting kernel arguments..." << std::endl;
	OCLDeviceGroup.SetBufferAsKernelArg(0, 0, 0);
	OCLDeviceGroup.SetBufferAsKernelArg(0, 1, 1);
	OCLDeviceGroup.SetKernelArg(0, 2, 10.00);

	std::cout << "Executing kernel..." << std::endl;
	OCLDeviceGroup.ExecuteKernel(0, DataSize / OCLDeviceGroup.DeviceNo);

	std::cout << "Copying the data from the device(s)..." << std::endl;
	for (cl_uint i = 0; i < OCLDeviceGroup.DeviceNo; i++)
	{
		OCLDeviceGroup.CopyBuffer(i, 1, Kratos::OpenCL::DeviceToHost, Output + i * DataSize / OCLDeviceGroup.DeviceNo);
	}

	std::cout << "Verification of results..." << std::endl;
	for (cl_uint i = 0; i < DataSize; i++)
	{
		if (Output[i] - (2.00 * sin(i) * cos(i) + 10.00) > Eps)
		{
			std::cout << "Error at position " << i << std::endl;
		}
	}

	std::cout << "Verification finished." << std::endl;

	return 0;
}
